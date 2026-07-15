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
    "route_attempt_registry.tsv": ["route_attempt_id", "sample_id", "decision_id", "pool_snapshot_id", "parent_pool_snapshot_id", "source_state", "route_class", "failure_mode", "applicability", "applicability_rationale", "n_query", "n_anchors", "query_membership_artifact", "query_membership_sha256", "query_only_graph", "depth_matched_validation", "calibration_origin", "calibration_manifest", "observed_density_spatial_prior", "candidate_resolutions", "selected_resolution", "parameters_artifact", "validation_artifact", "outcome_artifact", "residual_qc_membership_artifact", "residual_qc_membership_sha256", "prerequisite_route_attempt_id", "prerequisite_outcome_sha256", "prerequisite_residual_qc_sha256", "n_defined_fine", "n_defined_broad_only", "n_rerouted", "n_interface_retained", "n_qc_retained", "rctd_extreme_n", "rctd_high_n", "rctd_medium_low_n", "rctd_fine_return_n", "rctd_broad_return_n", "independent_fine_evidence", "successor_state", "reroute_membership_artifact", "reroute_membership_sha256", "fallback_route_attempt_id", "fallback_expected_membership_sha256", "fine_anchor_eligible", "status", "supersedes", "created_at", "closed_at"],
    "branch_control_board.tsv": ["branch_id", "sample_id", "parent_decision_id", "pool_snapshot_id", "generation", "run_id", "selected_resolution", "n_query", "membership_sha256", "current_state", "recluster_policy", "terminal", "next_action", "authoritative_artifact", "updated_at"],
    "workflow_event_registry.tsv": ["event_id", "sample_id", "timestamp", "phase", "branch_id", "action", "input_scope", "parameters", "scheduler_job_id", "status", "decision_summary_zh", "artifact", "supersedes_event_id"],
    "annotation_view_registry.tsv": ["view_id", "sample_id", "view", "membership_path", "membership_sha256", "n_observations", "policy", "marker_deg_eligible", "status", "artifact", "created_at"],
    "annotation_support_registry.tsv": ["support_id", "sample_id", "label_level", "broad_label", "fine_label", "n_observations", "confidence", "positive_marker_evidence", "anti_marker_evidence", "resolution_evidence", "spatial_evidence", "literature_context", "route_summary", "source_decision_ids", "validation_artifacts", "status", "supersedes", "created_at"],
    "run_registry.tsv": ["run_id", "work_key", "execution_fingerprint", "sample_id", "stage", "script", "parameters_path", "environment", "owner_assignment_id", "attempt", "scheduler_job_name", "scheduler_job_id", "status", "output_root", "supersedes_run_id", "started_at", "finished_at"],
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--input-root", required=True, type=Path)
    ap.add_argument("--project-root", required=True, type=Path)
    ap.add_argument("--modality", choices=["spatial", "single-cell"], required=True)
    ap.add_argument("--observation-unit", choices=["cell", "nucleus", "spot", "cellbin"], default="cell")
    ap.add_argument("--strategy-preset", choices=["none", "sheep_ovary_same_batch_rfirst"], default="none")
    args = ap.parse_args()
    root = args.project_root.resolve()
    for d in ["config", "state", "inputs", "runs", "tables", "figures", "spatial_nodes", "spatial_genes", "report", "provenance", "provenance/incidents", "logs"]:
        (root / d).mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat()
    config = {
        "framework_version": "1.5.1", "sample_id": args.sample,
        "input_root": str(args.input_root.resolve()), "project_root": str(root),
        "modality": args.modality, "observation_unit": args.observation_unit,
        "decision_version": "v001", "status": "initialized", "created_at_utc": now,
        "strategy_preset_requested": "" if args.strategy_preset == "none" else args.strategy_preset,
        "required_dotplot_levels": ["broad"],
        "subtype_assets_required_if_final_fine_present": True,
        "required_dotplot_panels": ["canonical", "data_specific"],
        "required_annotation_views": ["final"],
        "final_broad_minimum_confidence": "moderate",
        "final_fine_minimum_confidence": "high",
        "multi_route_completion_required": True,
        "post_completion_master_quality_approval_required": True,
        "preconfirmation_lightweight_review_required": True,
    }
    (root / "config/project.json").write_text(json.dumps(config, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for filename, fields in REGISTRIES.items():
        path = root / "state" / filename
        if not path.exists():
            with path.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle, delimiter="\t").writerow(fields)
    incident = root / "provenance/incidents/incident_registry.tsv"
    if not incident.exists():
        with incident.open("w", newline="", encoding="utf-8") as handle:
            csv.writer(handle, delimiter="\t").writerow([
                "incident_id", "detected_at_utc", "run_id", "scheduler_job_id", "scheduler_terminal",
                "failure_class", "failure_stage", "symptom", "root_cause", "failure_boundary",
                "accepted_prior_artifacts", "repair_run_id", "repair_job_id", "repair_action",
                "repair_verification", "state_mutated", "biological_labels_changed",
                "skill_prevention_candidate", "regression_test_candidate", "status", "evidence_paths",
            ])
    state_md = root / "state/annotation_state.md"
    if not state_md.exists():
        state_md.write_text(
            f"# {args.sample} annotation state\n\n- framework: 1.5.1\n- modality: {args.modality}\n"
            f"- observation unit: {args.observation_unit}\n- status: initialized\n- created: {now}\n\n"
            "## Current gate\n\nInput discovery and immutable snapshot are pending.\n",
            encoding="utf-8",
        )
    print(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
