#!/usr/bin/env python3
"""Stop an excessive whole-query P41 workload before jobs are submitted."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from controller_thresholds import load_controller_thresholds
from evidence_schema_lib import sha256
from lineage_controller_lib import read_tsv


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort-manifest", action="append", required=True, type=Path)
    ap.add_argument("--threshold-registry", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    trigger = load_controller_thresholds(args.threshold_registry)[
        "local_split_trigger_policy"
    ]["combined_analysis_split_fraction_review_trigger"]
    all_ids: set[str] = set()
    pending_ids: set[str] = set()
    records: list[dict[str, object]] = []
    for manifest_path in args.cohort_manifest:
        document = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            document.get("stage") != "cluster_cohort_recluster"
            or document.get("status") not in {"PASS", "LOCAL_SPLIT_REQUIRED"}
        ):
            raise SystemExit("local-split workload received a nonterminal cohort")
        base_record = document.get("base_candidate_membership", {})
        pending_record = document.get("pending_local_split_membership", {})
        base_path = Path(str(base_record.get("path", "")))
        pending_path = Path(str(pending_record.get("path", "")))
        if (
            not base_path.is_file() or base_record.get("sha256") != sha256(base_path)
            or not pending_path.is_file()
            or pending_record.get("sha256") != sha256(pending_path)
        ):
            raise SystemExit("local-split workload cohort membership is stale")
        base = {str(row.get("cell_id", "")) for row in read_tsv(base_path)}
        pending = {str(row.get("cell_id", "")) for row in read_tsv(pending_path)}
        cohort_ids = base | pending
        if "" in cohort_ids or base & pending or all_ids & cohort_ids:
            raise SystemExit("cohort memberships are empty, overlapping or duplicated")
        expected = int(document.get("n_observations", -1))
        if len(cohort_ids) != expected:
            raise SystemExit("cohort base/pending memberships do not cover the cohort")
        all_ids.update(cohort_ids)
        pending_ids.update(pending)
        records.append({
            "cohort_id": document.get("cohort_id", ""),
            "n_observations": len(cohort_ids),
            "n_pending_local_split": len(pending),
        })
    fraction = len(pending_ids) / len(all_ids) if all_ids else 1.0
    status = "REVIEW_REQUIRED" if fraction > trigger else "PASS"
    args.out.mkdir(parents=True, exist_ok=True)
    result = {
        "schema_version": "2.2", "status": status,
        "stage": "pre_local_split_workload_audit",
        "n_analysis_observations": len(all_ids),
        "n_pending_local_split": len(pending_ids),
        "pending_local_split_fraction": fraction,
        "review_trigger": trigger,
        "threshold_registry": {
            "path": str(args.threshold_registry.resolve()),
            "sha256": sha256(args.threshold_registry),
        },
        "required_action": (
            "revisit_second_round_resolution_or_split_trigger_before_P41_submission"
            if status == "REVIEW_REQUIRED" else "submit_bounded_local_splits"
        ),
        "cohorts": records,
        "cohort_manifests": [
            {"path": str(path.resolve()), "sha256": sha256(path)}
            for path in args.cohort_manifest
        ],
    }
    output = args.out / "local_split_workload_audit.json"
    output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 2 if status == "REVIEW_REQUIRED" else 0


if __name__ == "__main__":
    raise SystemExit(main())
