#!/usr/bin/env python3
"""Initialize a master-owned cohort control board with one logical worker per sample."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path


FIELDS = {
    "worker_registry.tsv": ["sample_id", "worker_id", "worker_generation", "wave", "status", "sample_root", "assigned_at", "updated_at", "supersedes_worker_id"],
    "sample_gate_registry.tsv": ["sample_id", "status", "packet_path", "packet_sha256", "completion_gate_path", "completion_gate_sha256", "release_audit_path", "release_audit_sha256", "master_decision", "updated_at"],
    "cohort_run_index.tsv": ["work_key", "sample_id", "stage", "owner_worker_id", "attempt", "scheduler_job_name", "scheduler_job_id", "status", "output_root", "supersedes_work_key", "updated_at"],
    "cohort_event_registry.tsv": ["event_id", "sample_id", "worker_id", "phase", "status", "summary", "artifact", "timestamp"],
}


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path, help="TSV with sample_id and input_root; sample_root is optional")
    parser.add_argument("--cohort-root", required=True, type=Path)
    parser.add_argument("--max-active-workers", type=int, default=1)
    parser.add_argument("--cohort-id", required=True)
    args = parser.parse_args()
    if args.max_active_workers < 1:
        raise SystemExit("--max-active-workers must be positive")
    source = read_manifest(args.manifest)
    if not source or "sample_id" not in source[0] or "input_root" not in source[0]:
        raise SystemExit("manifest requires sample_id and input_root")
    sample_ids = [row["sample_id"].strip() for row in source]
    if any(not item for item in sample_ids) or len(sample_ids) != len(set(sample_ids)):
        raise SystemExit("sample IDs must be nonempty and unique")

    root = args.cohort_root.resolve()
    for name in ["config", "control", "provenance/worker_packets", "state", "shared", "samples"]:
        (root / name).mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat()
    samples = []
    workers = []
    gates = []
    for index, row in enumerate(source):
        sample_root = Path(row.get("sample_root") or root / "samples" / row["sample_id"]).resolve()
        if root not in sample_root.parents and sample_root != root:
            # External roots are allowed only when explicitly present in the input manifest.
            if not row.get("sample_root"):
                raise SystemExit("derived sample root escaped cohort root")
        sample_root.mkdir(parents=True, exist_ok=True)
        wave = index // args.max_active_workers + 1
        samples.append({**row, "sample_root": str(sample_root), "wave": wave})
        workers.append({"sample_id": row["sample_id"], "worker_id": "", "worker_generation": 0, "wave": wave, "status": "PLANNED", "sample_root": str(sample_root), "assigned_at": "", "updated_at": now, "supersedes_worker_id": ""})
        gates.append({"sample_id": row["sample_id"], "status": "PLANNED", "packet_path": "", "packet_sha256": "", "completion_gate_path": "", "completion_gate_sha256": "", "release_audit_path": "", "release_audit_sha256": "", "master_decision": "", "updated_at": now})
    if len({row["sample_root"] for row in samples}) != len(samples):
        raise SystemExit("sample roots must be unique")
    sample_fields = list(source[0]) + [name for name in ["sample_root", "wave"] if name not in source[0]]
    write_tsv(root / "control/sample_manifest.tsv", sample_fields, samples)
    write_tsv(root / "control/worker_registry.tsv", FIELDS["worker_registry.tsv"], workers)
    write_tsv(root / "control/sample_gate_registry.tsv", FIELDS["sample_gate_registry.tsv"], gates)
    write_tsv(root / "control/cohort_run_index.tsv", FIELDS["cohort_run_index.tsv"], [])
    write_tsv(root / "control/cohort_event_registry.tsv", FIELDS["cohort_event_registry.tsv"], [])
    config = {"cohort_id": args.cohort_id, "status": "PLANNED", "master_agent_only_user_facing": True, "one_logical_worker_per_sample": True, "n_samples": len(samples), "max_active_workers": args.max_active_workers, "n_waves": max(row["wave"] for row in samples), "created_at": now}
    (root / "config/cohort.json").write_text(json.dumps(config, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(config, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
