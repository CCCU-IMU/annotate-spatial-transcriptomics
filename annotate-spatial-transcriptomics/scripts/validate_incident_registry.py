#!/usr/bin/env python3
"""Validate the standard incident registry and block completion on open incidents."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


REQUIRED = {
    "incident_id", "scheduler_job_id", "failure_class", "failure_stage", "symptom",
    "root_cause", "failure_boundary", "accepted_prior_artifacts", "repair_action",
    "repair_verification", "state_mutated", "biological_labels_changed",
    "skill_prevention_candidate", "regression_test_candidate", "status", "evidence_paths",
}
OPEN_MARKERS = ("open", "pending", "running", "in_progress", "not_repaired", "unresolved")


def validate(path: Path) -> dict:
    errors = []
    if not path.is_file():
        return {"status": "FAIL", "errors": [f"missing incident registry: {path}"], "rows": 0, "open_incidents": []}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = sorted(REQUIRED - set(fields))
        if missing:
            errors.append("missing fields: " + ",".join(missing))
        rows = list(reader)
    ids = [row.get("incident_id", "") for row in rows]
    if any(not value for value in ids):
        errors.append("blank incident_id")
    if len(ids) != len(set(ids)):
        errors.append("duplicate incident_id")
    open_ids = []
    for row in rows:
        status = row.get("status", "").strip().lower()
        if not status or any(marker in status for marker in OPEN_MARKERS):
            open_ids.append(row.get("incident_id", ""))
        if not row.get("failure_boundary", "").strip() or not row.get("accepted_prior_artifacts", "").strip():
            errors.append(f"{row.get('incident_id','?')} lacks boundary/reusable-artifact record")
    if open_ids:
        errors.append(f"{len(open_ids)} incident(s) remain open")
    return {"status": "PASS" if not errors else "FAIL", "errors": errors, "rows": len(rows), "open_incidents": open_ids}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("registry", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = validate(args.registry)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
