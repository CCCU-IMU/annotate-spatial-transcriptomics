#!/usr/bin/env python3
"""Non-destructively add tiered RCTD evidence fields to a v1.1 project."""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path


RCTD_FIELDS = [
    "rctd_extreme_n",
    "rctd_high_n",
    "rctd_medium_low_n",
    "rctd_fine_return_n",
    "rctd_broad_return_n",
    "independent_fine_evidence",
    "fallback_route_attempt_id",
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    registry = root / "state/route_attempt_registry.tsv"
    if not registry.exists():
        raise SystemExit("route_attempt_registry.tsv is missing; initialize or run the v1-to-v1.1 migration first")

    with registry.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        old_fields = reader.fieldnames or []
        rows = list(reader)
    new_fields = old_fields[:]
    insert_at = new_fields.index("fine_anchor_eligible") if "fine_anchor_eligible" in new_fields else len(new_fields)
    for field in RCTD_FIELDS:
        if field not in new_fields:
            new_fields.insert(insert_at, field)
            insert_at += 1

    backup = registry.with_name("route_attempt_registry.pre_v1_2.tsv")
    if not backup.exists():
        backup.write_bytes(registry.read_bytes())
    fd, tmp_name = tempfile.mkstemp(prefix=registry.name + ".", suffix=".tmp", dir=registry.parent)
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=new_fields, delimiter="\t")
            writer.writeheader()
            writer.writerows([{field: row.get(field, "") for field in new_fields} for row in rows])
            handle.flush(); os.fsync(handle.fileno())
        Path(tmp_name).replace(registry)
    finally:
        Path(tmp_name).unlink(missing_ok=True)

    config_path = root / "config/project.json"
    config = json.loads(config_path.read_text())
    config["framework_version"] = "1.2.0-dev"
    config["tiered_rctd_completion_required"] = True
    config_path.write_text(json.dumps(config, ensure_ascii=False, indent=2) + "\n")
    result = {
        "status": "MIGRATED_REQUIRES_RCTD_ROUTE_REVIEW",
        "framework_version": "1.2.0-dev",
        "route_attempt_rows": len(rows),
        "new_fields": RCTD_FIELDS,
        "backup": str(backup),
        "required_next_actions": [
            "populate tier fields for active applicable RCTD/reference-assisted routes",
            "route every medium/low observation to calibrated atlas/internal-anchor fallback",
            "supersede legacy terminal routes rather than editing their biological conclusions in place",
            "rerun next-iteration and completion audits",
        ],
        "migrated_at": datetime.now(timezone.utc).isoformat(),
    }
    provenance = root / "provenance"
    provenance.mkdir(exist_ok=True)
    (provenance / "migration_v1_1_to_v1_2.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
