#!/usr/bin/env python3
"""Add v1.4 run-ownership fields without rewriting biological decisions."""

from __future__ import annotations

import argparse
import csv
import shutil
from datetime import datetime, timezone
from pathlib import Path


NEW_FIELDS = ["run_id","work_key","execution_fingerprint","sample_id","stage","script","parameters_path","environment","owner_assignment_id","attempt","scheduler_job_id","status","output_root","supersedes_run_id","started_at","finished_at"]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    path = args.project_root / "state/run_registry.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    if all(field in fields for field in NEW_FIELDS):
        print("already_v1.4")
        return 0
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = args.project_root / "state/history" / f"run_registry.pre_v1.4.{stamp}.tsv"
    backup.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, backup)
    migrated = []
    for index, row in enumerate(rows, 1):
        run_id = row.get("run_id", "")
        migrated.append({**row, "work_key": row.get("work_key") or f"legacy:{row.get('stage','unknown')}:{run_id}", "execution_fingerprint": row.get("execution_fingerprint") or "legacy_unknown", "owner_assignment_id": row.get("owner_assignment_id") or "legacy_main_agent", "attempt": row.get("attempt") or "1", "supersedes_run_id": row.get("supersedes_run_id") or ""})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=NEW_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(migrated)
    print(f"migrated {len(migrated)} runs; backup={backup}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
