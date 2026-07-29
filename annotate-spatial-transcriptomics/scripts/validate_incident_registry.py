#!/usr/bin/env python3
"""Validate the standard incident registry and block completion on open incidents."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


REQUIRED = {
    "incident_id", "scheduler_job_id", "failure_class", "failure_stage", "symptom",
    "root_cause", "failure_boundary", "accepted_prior_artifacts", "repair_action",
    "repair_verification", "state_mutated", "biological_labels_changed",
    "skill_prevention_candidate", "regression_test_candidate", "status", "evidence_paths",
}
OPEN_MARKERS = ("open", "pending", "running", "in_progress", "not_repaired", "unresolved")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def semantic_registry_hash(rows: list[dict[str, str]], fields: list[str]) -> str:
    payload = "\n".join(
        "\t".join(str(row.get(field, "")) for field in fields)
        for row in sorted(rows, key=lambda row: str(row.get("incident_id", "")))
    )
    return hashlib.sha256((payload + "\n").encode("utf-8")).hexdigest()


def validate(path: Path) -> dict:
    errors = []
    if not path.is_file():
        return {
            "status": "FAIL", "errors": [f"missing incident registry: {path}"],
            "rows": 0, "open_incidents": [], "registry": None,
        }
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
    return {
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "rows": len(rows),
        "open_incidents": open_ids,
        "registry": {
            "path": str(path.resolve()),
            "sha256": sha256(path),
            "semantic_sha256": semantic_registry_hash(rows, fields),
            "row_count": len(rows),
        },
        "validated_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("registry", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument(
        "--verify-existing", type=Path,
        help="fail when a previously saved validation no longer binds this registry",
    )
    args = parser.parse_args()
    result = validate(args.registry)
    if args.verify_existing:
        existing = json.loads(args.verify_existing.read_text(encoding="utf-8"))
        observed = existing.get("registry") or {}
        current = result.get("registry") or {}
        if (
            observed.get("sha256") != current.get("sha256")
            or observed.get("semantic_sha256") != current.get("semantic_sha256")
            or observed.get("row_count") != current.get("row_count")
        ):
            result["status"] = "FAIL"
            result["errors"].append("existing incident-registry validation is stale")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
