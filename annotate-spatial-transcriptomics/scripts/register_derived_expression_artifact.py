#!/usr/bin/env python3
"""Register a project-scoped derived expression artifact with fail-closed ancestry."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from registry_io import locked_tsv_update


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--artifact-id", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--path", required=True, type=Path)
    parser.add_argument("--artifact-kind", required=True)
    parser.add_argument("--purpose", required=True)
    parser.add_argument("--parent-snapshot-id", default="")
    parser.add_argument("--parent-artifact-id", default="")
    parser.add_argument("--raw-counts-sha256", required=True)
    parser.add_argument("--analysis-set-sha256", required=True)
    parser.add_argument("--source-project-id", default="")
    parser.add_argument("--external-reference", action="store_true")
    parser.add_argument("--sha256-file", type=Path, help="scheduler-produced sha256sum record for a large derived object")
    args = parser.parse_args()

    root = args.project_root.resolve()
    project_path = root / "config/project.json"
    if not project_path.is_file():
        raise SystemExit("missing config/project.json")
    project = json.loads(project_path.read_text(encoding="utf-8"))
    project_id = str(project.get("project_id", "")).strip()
    if not project_id:
        raise SystemExit("project_id is missing from config/project.json")
    if args.sample != project.get("sample_id"):
        raise SystemExit("sample does not match the annotation project")
    source_project = args.source_project_id or project_id
    cross_project = source_project != project_id
    reference_purpose = args.purpose.lower() in {"atlas", "reference", "external_reference"}
    if cross_project and not (args.external_reference and reference_purpose):
        raise SystemExit("cross-project derived expression is forbidden outside the explicit external-reference channel")
    if args.external_reference and not reference_purpose:
        raise SystemExit("external-reference artifacts may only use atlas/reference purposes")
    if not args.parent_snapshot_id and not args.parent_artifact_id:
        raise SystemExit("a parent snapshot or parent derived artifact is required")

    path = args.path.resolve()
    if not path.is_file():
        raise SystemExit(f"artifact does not exist: {path}")
    if args.sha256_file:
        parts = args.sha256_file.read_text(encoding="utf-8").strip().split()
        if len(parts) < 2 or len(parts[0]) != 64 or any(ch not in "0123456789abcdefABCDEF" for ch in parts[0]):
            raise SystemExit("invalid sha256sum record")
        recorded = Path(parts[-1].lstrip("*"))
        if recorded.resolve() != path:
            raise SystemExit("sha256sum record path does not match --path")
        artifact_sha256 = parts[0].lower()
    else:
        artifact_sha256 = sha256(path)
    row = {
        "artifact_id": args.artifact_id, "project_id": project_id, "sample_id": args.sample,
        "path": str(path), "artifact_kind": args.artifact_kind, "purpose": args.purpose,
        "parent_snapshot_id": args.parent_snapshot_id, "parent_artifact_id": args.parent_artifact_id,
        "raw_counts_sha256": args.raw_counts_sha256, "analysis_set_sha256": args.analysis_set_sha256,
        "source_project_id": source_project, "external_reference": str(args.external_reference).lower(),
        "size_bytes": str(path.stat().st_size), "sha256": artifact_sha256, "status": "frozen",
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    registry = root / "state/derived_expression_registry.tsv"

    def mutate(rows, fields):
        missing = set(row) - set(fields)
        if missing:
            raise SystemExit(f"registry missing required columns: {sorted(missing)}")
        prior = [item for item in rows if item.get("artifact_id") == args.artifact_id]
        if prior:
            if prior[0] != row:
                stable = {k: v for k, v in row.items() if k != "created_at"}
                old = {k: v for k, v in prior[0].items() if k != "created_at"}
                if old != stable:
                    raise SystemExit("artifact ID already exists with different provenance")
            return rows
        rows.append(row)
        return rows

    locked_tsv_update(registry, mutate)
    print(json.dumps({"status": "REGISTERED", "artifact_id": args.artifact_id, "sha256": row["sha256"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
