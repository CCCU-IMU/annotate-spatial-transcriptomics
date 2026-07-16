#!/usr/bin/env python3
"""Small dependency-free JSON-schema subset plus artifact/membership validation."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def path_at(root: Path, value: str) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else root / path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def active_registry_rows(rows: list[dict[str, str]], id_field: str) -> list[dict[str, str]]:
    """Return append-only registry rows not named by a successor's supersedes field."""
    superseded = {
        value.strip()
        for row in rows
        for value in re.split(r"[|,;\s]+", row.get("supersedes", ""))
        if value.strip()
    }
    return [row for row in rows if row.get(id_field, "").strip() not in superseded]


def load_json_object(path: Path) -> tuple[dict[str, Any], list[str]]:
    if not path.is_file():
        return {}, [f"missing JSON artifact: {path}"]
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        return {}, [f"unreadable JSON artifact {path}: {exc}"]
    if not isinstance(value, dict) or not value:
        return {}, [f"JSON artifact must be a nonempty object: {path}"]
    return value, []


def _resolve(schema_root: dict, schema: dict) -> dict:
    ref = schema.get("$ref")
    if not ref:
        return schema
    if not ref.startswith("#/"):
        raise ValueError(f"unsupported external schema ref: {ref}")
    value: Any = schema_root
    for token in ref[2:].split("/"):
        value = value[token.replace("~1", "/").replace("~0", "~")]
    return value


def validate_schema(instance: Any, schema: dict, schema_root: dict | None = None, location: str = "$") -> list[str]:
    root = schema_root or schema
    schema = _resolve(root, schema)
    errors: list[str] = []
    if "const" in schema and instance != schema["const"]:
        errors.append(f"{location}: expected constant {schema['const']!r}")
    if "enum" in schema and instance not in schema["enum"]:
        errors.append(f"{location}: value {instance!r} is not in {schema['enum']!r}")
    expected = schema.get("type")
    type_ok = {
        "object": isinstance(instance, dict),
        "array": isinstance(instance, list),
        "string": isinstance(instance, str),
        "number": isinstance(instance, (int, float)) and not isinstance(instance, bool),
        "integer": isinstance(instance, int) and not isinstance(instance, bool),
        "boolean": isinstance(instance, bool),
    }.get(expected, True)
    if not type_ok:
        return errors + [f"{location}: expected {expected}, found {type(instance).__name__}"]
    if isinstance(instance, dict):
        for key in schema.get("required", []):
            if key not in instance:
                errors.append(f"{location}: missing required property {key}")
        if len(instance) < int(schema.get("minProperties", 0)):
            errors.append(f"{location}: object has too few properties")
        for key, child in schema.get("properties", {}).items():
            if key in instance:
                errors.extend(validate_schema(instance[key], child, root, f"{location}.{key}"))
    elif isinstance(instance, list):
        if len(instance) < int(schema.get("minItems", 0)):
            errors.append(f"{location}: array has too few items")
        if schema.get("uniqueItems") and len({json.dumps(value, sort_keys=True) for value in instance}) != len(instance):
            errors.append(f"{location}: array items are not unique")
        if "items" in schema:
            for index, value in enumerate(instance):
                errors.extend(validate_schema(value, schema["items"], root, f"{location}[{index}]"))
    elif isinstance(instance, str):
        if len(instance) < int(schema.get("minLength", 0)):
            errors.append(f"{location}: string is too short")
        if schema.get("pattern") and not re.fullmatch(schema["pattern"], instance):
            errors.append(f"{location}: string does not match required pattern")
    elif isinstance(instance, (int, float)) and "minimum" in schema and instance < schema["minimum"]:
        errors.append(f"{location}: value is below minimum")
    return errors


def validate_json_against_schema(document_path: Path, schema_path: Path) -> tuple[dict[str, Any], list[str]]:
    document, errors = load_json_object(document_path)
    schema, schema_errors = load_json_object(schema_path)
    errors.extend(schema_errors)
    if not errors:
        errors.extend(validate_schema(document, schema))
    return document, errors


def validate_artifact_ref(root: Path, artifact: dict[str, Any], label: str, require_nonempty: bool = True) -> tuple[Path | None, list[str]]:
    errors: list[str] = []
    path_value = str(artifact.get("path", ""))
    expected = str(artifact.get("sha256", ""))
    if not path_value or not re.fullmatch(r"[0-9a-f]{64}", expected):
        return None, [f"{label}: missing path or canonical SHA256"]
    path = path_at(root, path_value)
    if not path.is_file():
        return None, [f"{label}: artifact does not exist: {path_value}"]
    if require_nonempty and path.stat().st_size == 0:
        errors.append(f"{label}: artifact is empty")
    if sha256(path) != expected:
        errors.append(f"{label}: artifact SHA256 is stale")
    return path, errors


def validate_evidence_artifact(root: Path, artifact: dict[str, Any], label: str) -> tuple[Path | None, list[str]]:
    """Reject placeholder artifacts even when their path and hash are valid."""
    path, errors = validate_artifact_ref(root, artifact, label)
    if path is None or errors:
        return path, errors
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".json"):
        value, json_errors = load_json_object(path)
        errors.extend(json_errors)
        if value and set(value) <= {"status", "created_at", "completed_at"}:
            errors.append(f"{label}: JSON contains status metadata but no biological evidence")
    elif any(suffixes.endswith(value) for value in (".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        rows = read_tsv(path)
        if not rows:
            errors.append(f"{label}: tabular artifact has no evidence rows")
        elif len(rows[0]) < 2:
            errors.append(f"{label}: tabular artifact has no evidence columns")
    return path, errors


def membership_ids(root: Path, membership: dict[str, Any], label: str) -> tuple[set[str], list[str]]:
    path, errors = validate_artifact_ref(root, membership, label)
    if path is None or errors:
        return set(), errors
    rows = read_tsv(path)
    if not rows or "cell_id" not in rows[0]:
        return set(), errors + [f"{label}: membership must be a nonempty TSV with cell_id"]
    ids = [row.get("cell_id", "") for row in rows]
    if "" in ids or len(ids) != len(set(ids)):
        errors.append(f"{label}: membership cell_id values must be unique and nonempty")
    expected_n = membership.get("n_observations")
    if not isinstance(expected_n, int) or expected_n != len(ids):
        errors.append(f"{label}: n_observations differs from membership")
    return set(ids), errors


def write_result(path: Path | None, result: dict[str, Any]) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
