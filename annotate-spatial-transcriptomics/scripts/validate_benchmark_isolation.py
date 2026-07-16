#!/usr/bin/env python3
"""Fail closed when a held-out benchmark leaks into reference or Atlas configuration."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(root: Path, config_path: Path, frozen_annotation: Path | None = None) -> dict:
    errors: list[str] = []
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return {"status": "FAIL", "errors": [f"benchmark config unreadable: {exc}"]}
    required = {
        "benchmark_mode": True, "author_labels_locked_until": "annotation_frozen",
        "unblinding_requires_frozen_hash": True,
    }
    for key, expected in required.items():
        if config.get(key) != expected:
            errors.append(f"benchmark config {key} must equal {expected!r}")
    held_out = str(config.get("held_out_dataset_id", "")).strip()
    forbidden = {str(value).strip() for value in config.get("forbidden_reference_ids", []) if str(value).strip()}
    if not held_out or held_out not in forbidden:
        errors.append("held-out dataset must be included in forbidden_reference_ids")
    registry_paths = [
        root / "state/reference_registry.tsv", root / "state/route_attempt_registry.tsv",
        root / "config/atlas_config.json", root / "config/reference_config.json",
    ]
    hits: list[str] = []
    for path in registry_paths:
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        for reference_id in forbidden:
            if reference_id and reference_id in text:
                hits.append(f"{reference_id}@{path.relative_to(root)}")
    if hits:
        errors.append("held-out benchmark leaked into reference/route/Atlas configuration: " + ", ".join(hits))
    frozen_hash = ""
    if frozen_annotation:
        if not frozen_annotation.is_file():
            errors.append("frozen annotation artifact is missing")
        else:
            frozen_hash = sha256(frozen_annotation)
    result = {
        "status": "PASS" if not errors else "FAIL", "benchmark_mode": bool(config.get("benchmark_mode")),
        "held_out_dataset_id": held_out, "forbidden_reference_ids": sorted(forbidden),
        "author_labels_locked": frozen_annotation is None, "frozen_annotation_sha256": frozen_hash,
        "errors": errors,
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--frozen-annotation", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = validate(args.project_root.resolve(), args.config.resolve(), args.frozen_annotation)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
