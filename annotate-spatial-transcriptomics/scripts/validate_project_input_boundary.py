#!/usr/bin/env python3
"""Validate project-local expression ancestry without rereading expression matrices."""

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
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(root: Path) -> dict[str, object]:
    root = root.resolve()
    errors: list[str] = []
    project_path = root / "config/project.json"
    if not project_path.is_file():
        return {"status": "BLOCKED", "errors": ["missing config/project.json"]}
    project = json.loads(project_path.read_text(encoding="utf-8"))
    project_id = str(project.get("project_id", "")).strip()
    sample = str(project.get("sample_id", "")).strip()
    if not project_id:
        errors.append("project_id is missing")
    snapshots = {row.get("snapshot_id", ""): row for row in read_tsv(root / "state/input_snapshot_registry.tsv")}
    artifacts = read_tsv(root / "state/derived_expression_registry.tsv")
    by_id = {row.get("artifact_id", ""): row for row in artifacts}
    for row in artifacts:
        aid = row.get("artifact_id", "")
        if row.get("project_id") != project_id or row.get("sample_id") != sample:
            errors.append(f"{aid}: project/sample boundary mismatch")
        path = Path(row.get("path", ""))
        if not path.is_file() or sha256(path) != row.get("sha256"):
            errors.append(f"{aid}: artifact is missing or stale")
        if not row.get("raw_counts_sha256") or not row.get("analysis_set_sha256"):
            errors.append(f"{aid}: raw-counts or analysis-set hash is missing")
        parent_snapshot = row.get("parent_snapshot_id", "")
        parent_artifact = row.get("parent_artifact_id", "")
        if not parent_snapshot and not parent_artifact:
            errors.append(f"{aid}: ancestry is missing")
        if parent_snapshot and parent_snapshot not in snapshots:
            errors.append(f"{aid}: parent snapshot is not registered")
        if parent_artifact and parent_artifact not in by_id:
            errors.append(f"{aid}: parent derived artifact is not registered")
        cross_project = row.get("source_project_id", "") != project_id
        external = row.get("external_reference", "").lower() == "true"
        reference_purpose = row.get("purpose", "").lower() in {"atlas", "reference", "external_reference"}
        if cross_project and not (external and reference_purpose):
            errors.append(f"{aid}: cross-project expression entered an annotation evidence channel")
        if external and not reference_purpose:
            errors.append(f"{aid}: external-reference flag is invalid for purpose {row.get('purpose','')}")
    return {"status": "PASS" if not errors else "BLOCKED", "project_id": project_id,
            "sample_id": sample, "snapshots": len(snapshots), "derived_artifacts": len(artifacts), "errors": errors}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = validate(args.project_root)
    out = args.out or args.project_root / "provenance/project_input_boundary_validation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
