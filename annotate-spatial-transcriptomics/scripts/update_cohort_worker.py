#!/usr/bin/env python3
"""Master-only atomic assignment/state update for one sample worker."""

from __future__ import annotations

import argparse
import csv
import fcntl
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path


ALLOWED = ["PLANNED", "WORKER_ASSIGNED", "ANALYSIS_RUNNING", "READY_FOR_MASTER_QUALITY_GATE", "MASTER_QUALITY_APPROVED", "SAMPLE_FROZEN", "CROSS_SAMPLE_AUDIT", "COHORT_CONFIRMATION_PENDING", "RELEASE_RUNNING", "RELEASED", "FAILED_PRESERVED"]
ACTIVE = {"WORKER_ASSIGNED", "ANALYSIS_RUNNING", "RELEASE_RUNNING"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort-root", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--worker-id", required=True)
    parser.add_argument("--status", choices=ALLOWED, required=True)
    parser.add_argument("--takeover", action="store_true")
    args = parser.parse_args()
    path = args.cohort_root / "control/worker_registry.tsv"
    lock_path = path.with_suffix(".tsv.lock")
    with lock_path.open("w") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            rows = list(reader)
        matches = [row for row in rows if row["sample_id"] == args.sample]
        if len(matches) != 1:
            raise SystemExit("sample assignment row is missing or duplicated")
        row = matches[0]
        previous = row.get("worker_id", "")
        if previous and previous != args.worker_id and row.get("status") in ACTIVE and not args.takeover:
            raise SystemExit("sample already has an active worker; use --takeover after auditing live jobs")
        if previous and previous != args.worker_id:
            row["supersedes_worker_id"] = previous
            row["worker_generation"] = str(int(row.get("worker_generation", 0) or 0) + 1)
        row["worker_id"] = args.worker_id
        row["status"] = args.status
        now = datetime.now(timezone.utc).isoformat()
        row["assigned_at"] = row.get("assigned_at") or now
        row["updated_at"] = now
        fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
        try:
            with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader(); writer.writerows(rows); handle.flush(); os.fsync(handle.fileno())
            Path(name).replace(path)
        finally:
            Path(name).unlink(missing_ok=True)
    print(args.sample, args.worker_id, args.status)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
