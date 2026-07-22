#!/usr/bin/env python3
"""Freeze the full-object versus analysis-set boundary as a hash-bound policy."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt", newline="", encoding="utf-8") if path.suffix == ".gz" else path.open(newline="", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--membership", required=True, type=Path,
                        help="TSV/TSV.GZ with unique cell_id and analysis_scope")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve(); membership = args.membership.resolve()
    errors: list[str] = []
    with open_text(membership) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = set(rows[0]) if rows else set()
    if not rows or not {"cell_id", "analysis_scope"}.issubset(fields):
        errors.append("membership must contain cell_id and analysis_scope")
    ids = [row.get("cell_id", "") for row in rows]
    if "" in ids or len(ids) != len(set(ids)):
        errors.append("membership cell_id values must be unique and nonempty")
    allowed = {"analysis_set", "excluded_initial_qc"}
    scopes = [row.get("analysis_scope", "") for row in rows]
    if any(value not in allowed for value in scopes):
        errors.append("analysis_scope must be analysis_set or excluded_initial_qc")
    analysis_n = sum(value == "analysis_set" for value in scopes)
    excluded_n = sum(value == "excluded_initial_qc" for value in scopes)
    if analysis_n == 0:
        errors.append("analysis set is empty")
    payload = {
        "status": "PASS" if not errors else "FAIL",
        "schema_version": "2.0.1",
        "membership_path": str(membership),
        "membership_sha256": sha256(membership),
        "full_object_n": len(rows),
        "analysis_set_n": analysis_n,
        "excluded_initial_qc_n": excluded_n,
        "analysis_scope_values": sorted(allowed),
        "policy": "immutable_full_object_partition",
        "errors": errors,
    }
    out = args.out or root / "provenance/analysis_scope_policy.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
