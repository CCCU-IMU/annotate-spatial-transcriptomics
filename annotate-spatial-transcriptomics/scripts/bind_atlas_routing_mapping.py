#!/usr/bin/env python3
"""Hash-bind an exact target plus held-out-query union for Atlas routing."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    return gzip.open(path, "rt", newline="", encoding="utf-8") if path.suffix == ".gz" else path.open(newline="", encoding="utf-8")


def ids(path: Path, column: str) -> set[str]:
    with open_text(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    values = [row.get(column, "") for row in rows]
    if "" in values or len(values) != len(set(values)):
        raise SystemExit(f"{path}: {column} must be unique and nonempty")
    return set(values)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--calibration-manifest", required=True, type=Path)
    parser.add_argument("--heldout-mapping", required=True, type=Path)
    parser.add_argument("--combined-mapping", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cell-id-col", default="cell_id")
    args = parser.parse_args()
    source = args.calibration_manifest.resolve()
    payload = json.loads(source.read_text(encoding="utf-8"))
    target_record = payload.get("artifacts", {}).get("query_mapping", {})
    target = Path(target_record.get("path", "")).resolve()
    if not target.is_file() or sha256(target) != target_record.get("sha256"):
        raise SystemExit("calibration target mapping is missing or stale")
    target_ids = ids(target, args.cell_id_col)
    heldout_ids = ids(args.heldout_mapping, args.cell_id_col)
    combined_ids = ids(args.combined_mapping, args.cell_id_col)
    if target_ids & heldout_ids:
        raise SystemExit("target and held-out mappings overlap")
    if combined_ids != target_ids | heldout_ids:
        raise SystemExit("combined mapping is not the exact target plus held-out union")
    payload["artifacts"]["query_mapping"] = {
        "path": str(args.combined_mapping.resolve()), "sha256": sha256(args.combined_mapping)
    }
    payload["routing_mapping_binding"] = {
        "status": "PASS",
        "semantics": "exact disjoint union of calibrated target and current-query held-out mapping rows",
        "source_calibration_manifest": str(source),
        "source_calibration_manifest_sha256": sha256(source),
        "target_mapping": {"path": str(target), "sha256": sha256(target), "n": len(target_ids)},
        "heldout_mapping": {"path": str(args.heldout_mapping.resolve()), "sha256": sha256(args.heldout_mapping), "n": len(heldout_ids)},
        "combined_mapping": {"path": str(args.combined_mapping.resolve()), "sha256": sha256(args.combined_mapping), "n": len(combined_ids)},
        "overlap_n": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
