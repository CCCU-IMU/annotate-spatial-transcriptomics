#!/usr/bin/env python3
"""Fail-closed audit of master/worker cohort ownership and release state."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from scheduler_job_name import validate as validate_scheduler_job_name


ACTIVE = {"WORKER_ASSIGNED", "ANALYSIS_RUNNING", "RELEASE_RUNNING"}


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cohort_root", type=Path)
    args = parser.parse_args()
    root = args.cohort_root.resolve()
    errors: list[str] = []
    samples = read(root / "control/sample_manifest.tsv")
    workers = read(root / "control/worker_registry.tsv")
    gates = read(root / "control/sample_gate_registry.tsv")
    runs = read(root / "control/cohort_run_index.tsv")
    sample_ids = [row.get("sample_id", "") for row in samples]
    if not sample_ids or len(sample_ids) != len(set(sample_ids)) or "" in sample_ids:
        errors.append("sample manifest IDs are empty or duplicated")
    roots = [row.get("sample_root", "") for row in samples]
    if len(roots) != len(set(roots)) or "" in roots:
        errors.append("sample roots are empty or shared")
    if sorted(row.get("sample_id", "") for row in workers) != sorted(sample_ids):
        errors.append("worker registry does not contain exactly one row per sample")
    if sorted(row.get("sample_id", "") for row in gates) != sorted(sample_ids):
        errors.append("sample gate registry does not contain exactly one row per sample")
    active_workers = [row.get("worker_id", "") for row in workers if row.get("status") in ACTIVE]
    if "" in active_workers or len(active_workers) != len(set(active_workers)):
        errors.append("active worker IDs are empty or assigned to multiple samples")
    active_work_keys = [row.get("work_key", "") for row in runs if row.get("status") in {"submitted", "running"}]
    if "" in active_work_keys or len(active_work_keys) != len(set(active_work_keys)):
        errors.append("active cohort work keys are empty or duplicated")
    for line, row in enumerate(runs, 2):
        if row.get("status") not in {"submitted", "running"}:
            continue
        try:
            validate_scheduler_job_name(row.get("scheduler_job_name", ""), 48)
        except ValueError as exc:
            errors.append(f"cohort run index line {line}: invalid scheduler_job_name: {exc}")
    unknown_run_samples = sorted({row.get("sample_id", "") for row in runs}.difference(sample_ids))
    if unknown_run_samples:
        errors.append(f"run index contains unknown samples: {unknown_run_samples}")
    for row in workers:
        registered = Path(row.get("sample_root", "")).resolve()
        sample = next((item for item in samples if item.get("sample_id") == row.get("sample_id")), None)
        if sample and registered != Path(sample["sample_root"]).resolve():
            errors.append(f"worker root drift for {row.get('sample_id')}")
    confirmation_path = root / "state/cohort_confirmation.json"
    if confirmation_path.exists():
        confirmation = json.loads(confirmation_path.read_text(encoding="utf-8"))
        request_path = Path(confirmation.get("request_path", ""))
        if not request_path.is_file():
            errors.append("confirmed cohort request is missing")
        elif hashlib.sha256(request_path.read_bytes()).hexdigest() != confirmation.get("request_sha256"):
            errors.append("cohort confirmation was invalidated by a changed request")
    config = json.loads((root / "config/cohort.json").read_text(encoding="utf-8"))
    if config.get("status") == "RELEASED":
        if not confirmation_path.exists():
            errors.append("released cohort lacks confirmation")
        for row in gates:
            if row.get("status") != "RELEASED":
                errors.append(f"sample {row.get('sample_id')} is not RELEASED")
            audit_path = Path(row.get("release_audit_path", ""))
            if not audit_path.is_absolute():
                audit_path = root / audit_path
            if not audit_path.is_file():
                errors.append(f"sample {row.get('sample_id')} release audit is missing")
                continue
            observed = hashlib.sha256(audit_path.read_bytes()).hexdigest()
            if observed != row.get("release_audit_sha256"):
                errors.append(f"sample {row.get('sample_id')} release audit hash is stale")
            else:
                try:
                    if json.loads(audit_path.read_text(encoding="utf-8")).get("status") != "PASS":
                        errors.append(f"sample {row.get('sample_id')} release audit is not PASS")
                except json.JSONDecodeError:
                    errors.append(f"sample {row.get('sample_id')} release audit is unreadable")
    result = {"status": "PASS" if not errors else "FAIL", "n_samples": len(samples), "n_active_workers": len(active_workers), "errors": errors}
    out = root / "provenance/cohort_state_validation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
