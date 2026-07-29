#!/usr/bin/env python3
"""Plan cohort memory without forking copies of a large Seurat carrier."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def membership_n(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "cell_id" not in (reader.fieldnames or []):
            raise ValueError("cohort membership lacks cell_id")
        seen: set[str] = set()
        for row in reader:
            cell = str(row.get("cell_id", ""))
            if not cell or cell in seen:
                raise ValueError("cohort membership has blank or duplicate cell_id")
            seen.add(cell)
    return len(seen)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rds", required=True, type=Path)
    ap.add_argument("--rds-sha256")
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--resolution-workers", type=int, default=1)
    ap.add_argument("--grid-size", type=int, default=5)
    ap.add_argument("--allocated-memory-gb", type=float)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    n = membership_n(args.membership)
    rds_sha256 = args.rds_sha256 or sha256(args.rds)
    if len(rds_sha256) != 64 or any(
        value not in "0123456789abcdef" for value in rds_sha256.lower()
    ):
        raise SystemExit("--rds-sha256 is not a SHA256 digest")
    rds_gb = args.rds.stat().st_size / 1024 ** 3
    requested_workers = max(1, min(args.resolution_workers, args.grid_size))
    effective_workers = 1 if n >= 20_000 else requested_workers
    # Conservative empirical envelope from the A1/C3 incidents. The first
    # worker owns the carrier, SCT intermediates and graph; each additional
    # resolution worker can duplicate a substantial fraction of all three.
    single_worker_gb = 20.0 + 1.75 * rds_gb + 0.00012 * n
    copied_worker_gb = 1.25 * rds_gb + 0.00006 * n
    estimated_peak_gb = 1.20 * (
        single_worker_gb + max(0, effective_workers - 1) * copied_worker_gb
    )
    recommended_memory_gb = max(
        32, int(math.ceil(estimated_peak_gb / 16.0) * 16)
    )
    recommended_cpu = 64 if n >= 20_000 else 32 if n >= 5_000 else 16
    allocated = args.allocated_memory_gb
    if allocated is None and os.environ.get("AIP_MEMORY_GB"):
        try:
            allocated = float(os.environ["AIP_MEMORY_GB"])
        except ValueError:
            allocated = None
    errors: list[str] = []
    if allocated is not None and estimated_peak_gb > 0.90 * allocated:
        errors.append(
            "predicted cohort peak exceeds 90% of the declared memory allocation"
        )
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema_version": "2.2",
        "status": "PASS" if not errors else "BLOCKED",
        "artifact_role": "cohort_resource_plan",
        "source_rds": {
            "path": str(args.rds.resolve()), "sha256": rds_sha256,
            "n_bytes": args.rds.stat().st_size,
        },
        "membership": {
            "path": str(args.membership.resolve()),
            "sha256": sha256(args.membership), "n_observations": n,
        },
        "requested_resolution_workers": requested_workers,
        "effective_resolution_workers": effective_workers,
        "large_carrier_resolution_workers_forced_to_one": bool(
            n >= 20_000 and requested_workers > 1
        ),
        "estimated_peak_memory_gb": round(estimated_peak_gb, 2),
        "declared_allocation_memory_gb": allocated,
        "recommended_memory_gb": recommended_memory_gb,
        "recommended_cpu": recommended_cpu,
        "parallelism_policy": (
            "parallelize independent cohort jobs; keep one loaded Seurat carrier "
            "serial across resolution evaluation"
        ),
        "errors": errors,
    }
    path = args.out / "cohort_resource_plan.json"
    path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
