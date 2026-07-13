#!/usr/bin/env python3
"""Create a versioned annotation workspace and empty state registries."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path


REGISTRIES = {
    "input_snapshot_registry.tsv": ["snapshot_id", "sample_id", "path", "kind", "size_bytes", "sha256", "status", "created_at"],
    "clustering_decision_ledger.tsv": ["decision_version", "sample_id", "method", "run_id", "parameters", "n_clusters", "quantitative_rank", "marker_review", "spatial_review", "decision", "rationale", "created_at"],
    "cluster_decision_ledger.tsv": ["decision_version", "decision_id", "sample_id", "source_run_id", "source_cluster", "n_observations", "spatial_object_count", "count_interpretation", "broad_label", "fine_label", "state", "confidence", "evidence_status", "validation_status", "validation_artifact", "validation_feature_scope", "route", "route_run_id", "iteration", "target_pool", "fine_anchor_eligible", "next_action", "closure_rationale", "supersedes", "closed", "created_at"],
    "pool_registry.tsv": ["pool_id", "sample_id", "parent_pool_id", "membership_path", "membership_sha256", "n_observations", "purpose", "status", "decision_version", "created_at", "closed_at"],
    "pool_snapshot_registry.tsv": ["snapshot_id", "pool_id", "sample_id", "generation", "membership_path", "membership_sha256", "n_query", "n_anchors", "parent_decision_ids", "supersedes_snapshot_id", "status", "created_at", "closed_at"],
    "route_attempt_registry.tsv": ["route_attempt_id", "sample_id", "decision_id", "pool_snapshot_id", "route_class", "failure_mode", "applicability", "applicability_rationale", "n_query", "n_anchors", "query_only_graph", "depth_matched_validation", "observed_density_spatial_prior", "selected_resolution", "parameters_artifact", "validation_artifact", "outcome_artifact", "n_defined_fine", "n_defined_broad_only", "n_rerouted", "n_interface_retained", "n_qc_retained", "rctd_extreme_n", "rctd_high_n", "rctd_medium_low_n", "rctd_fine_return_n", "rctd_broad_return_n", "independent_fine_evidence", "fallback_route_attempt_id", "fine_anchor_eligible", "status", "supersedes", "created_at", "closed_at"],
    "branch_control_board.tsv": ["branch_id", "sample_id", "parent_decision_id", "pool_snapshot_id", "generation", "run_id", "selected_resolution", "n_query", "membership_sha256", "current_state", "recluster_policy", "terminal", "next_action", "authoritative_artifact", "updated_at"],
    "workflow_event_registry.tsv": ["event_id", "sample_id", "timestamp", "phase", "branch_id", "action", "input_scope", "parameters", "scheduler_job_id", "status", "decision_summary_zh", "artifact", "supersedes_event_id"],
    "annotation_view_registry.tsv": ["view_id", "sample_id", "view", "membership_path", "membership_sha256", "n_observations", "policy", "marker_deg_eligible", "status", "artifact", "created_at"],
    "run_registry.tsv": ["run_id", "sample_id", "stage", "script", "parameters_path", "environment", "scheduler_job_id", "status", "output_root", "started_at", "finished_at"],
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--input-root", required=True, type=Path)
    ap.add_argument("--project-root", required=True, type=Path)
    ap.add_argument("--modality", choices=["spatial", "single-cell"], required=True)
    ap.add_argument("--observation-unit", choices=["cell", "nucleus", "spot", "cellbin"], default="cell")
    args = ap.parse_args()
    root = args.project_root.resolve()
    for d in ["config", "state", "inputs", "runs", "tables", "figures", "spatial_nodes", "spatial_genes", "report", "provenance", "logs"]:
        (root / d).mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat()
    config = {
        "framework_version": "1.2.0-dev", "sample_id": args.sample,
        "input_root": str(args.input_root.resolve()), "project_root": str(root),
        "modality": args.modality, "observation_unit": args.observation_unit,
        "decision_version": "v001", "status": "initialized", "created_at_utc": now,
        "required_dotplot_levels": ["broad", "subtype"],
        "required_dotplot_panels": ["canonical", "data_specific"],
        "required_annotation_views": ["strict", "inclusive", "display"],
        "multi_route_completion_required": True,
    }
    (root / "config/project.json").write_text(json.dumps(config, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for filename, fields in REGISTRIES.items():
        path = root / "state" / filename
        if not path.exists():
            with path.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle, delimiter="\t").writerow(fields)
    state_md = root / "state/annotation_state.md"
    if not state_md.exists():
        state_md.write_text(
            f"# {args.sample} annotation state\n\n- framework: 1.2.0-dev\n- modality: {args.modality}\n"
            f"- observation unit: {args.observation_unit}\n- status: initialized\n- created: {now}\n\n"
            "## Current gate\n\nInput discovery and immutable snapshot are pending.\n",
            encoding="utf-8",
        )
    print(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
