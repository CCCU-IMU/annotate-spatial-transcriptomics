#!/usr/bin/env python3
"""Append or terminally update one machine-auditable route attempt."""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from multiroute_lib import ROUTE_CLASSES, valid_attempt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--record", required=True, type=Path, help="JSON object containing route_attempt_registry fields")
    args = parser.parse_args()
    registry = args.project_root / "state/route_attempt_registry.tsv"
    record = json.loads(args.record.read_text())
    missing = [x for x in ["route_attempt_id", "sample_id", "decision_id", "route_class", "applicability", "status"] if not record.get(x)]
    if missing:
        raise SystemExit(f"missing route fields: {missing}")
    if record["route_class"] not in ROUTE_CLASSES:
        raise SystemExit(f"invalid route_class: {record['route_class']}")
    now = datetime.now(timezone.utc).isoformat()
    record.setdefault("created_at", now)
    if record.get("status") in {"validated", "not_applicable_reviewed"}:
        record.setdefault("closed_at", now)
        ok, errors = valid_attempt(args.project_root, {k: str(v) for k, v in record.items()})
        if not ok:
            raise SystemExit("terminal route attempt is invalid: " + "; ".join(errors))
    lock_path = registry.with_suffix(registry.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("w") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        if not registry.exists():
            raise SystemExit("route_attempt_registry.tsv is missing; initialize or migrate the project")
        with registry.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            rows = list(reader)
        old = [x for x in rows if x.get("route_attempt_id") == record["route_attempt_id"]]
        if old and old[0].get("status") in {"validated", "not_applicable_reviewed"}:
            raise SystemExit("terminal route_attempt_id is immutable; create a superseding ID")
        rows = [x for x in rows if x.get("route_attempt_id") != record["route_attempt_id"]]
        rows.append({field: record.get(field, "") for field in fields})
        fd, tmp_name = tempfile.mkstemp(prefix=registry.name + ".", suffix=".tmp", dir=registry.parent)
        try:
            with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader(); writer.writerows(rows)
                handle.flush(); os.fsync(handle.fileno())
            Path(tmp_name).replace(registry)
        finally:
            Path(tmp_name).unlink(missing_ok=True)
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
    print(record["route_attempt_id"], record["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
