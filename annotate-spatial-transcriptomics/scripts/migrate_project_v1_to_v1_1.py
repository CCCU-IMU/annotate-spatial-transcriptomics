#!/usr/bin/env python3
"""Non-destructively add v1.1 multi-route registries to a v1 project."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


REGISTRIES = {
    "pool_snapshot_registry.tsv": ["snapshot_id","pool_id","sample_id","generation","membership_path","membership_sha256","n_query","n_anchors","parent_decision_ids","supersedes_snapshot_id","status","created_at","closed_at"],
    "route_attempt_registry.tsv": ["route_attempt_id","sample_id","decision_id","pool_snapshot_id","route_class","failure_mode","applicability","applicability_rationale","n_query","n_anchors","query_only_graph","depth_matched_validation","observed_density_spatial_prior","selected_resolution","parameters_artifact","validation_artifact","outcome_artifact","n_defined_fine","n_defined_broad_only","n_rerouted","n_interface_retained","n_qc_retained","rctd_extreme_n","rctd_high_n","rctd_medium_low_n","rctd_fine_return_n","rctd_broad_return_n","independent_fine_evidence","fallback_route_attempt_id","fine_anchor_eligible","status","supersedes","created_at","closed_at"],
    "branch_control_board.tsv": ["branch_id","sample_id","parent_decision_id","pool_snapshot_id","generation","run_id","selected_resolution","n_query","membership_sha256","current_state","recluster_policy","terminal","next_action","authoritative_artifact","updated_at"],
    "workflow_event_registry.tsv": ["event_id","sample_id","timestamp","phase","branch_id","action","input_scope","parameters","scheduler_job_id","status","decision_summary_zh","artifact","supersedes_event_id"],
    "annotation_view_registry.tsv": ["view_id","sample_id","view","membership_path","membership_sha256","n_observations","policy","marker_deg_eligible","status","artifact","created_at"],
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""): h.update(block)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("project_root", type=Path); args = parser.parse_args()
    root = args.project_root.resolve(); state = root / "state"; provenance = root / "provenance"; provenance.mkdir(exist_ok=True)
    config_path = root / "config/project.json"; config = json.loads(config_path.read_text()); backup = root / "config/project.pre_v1_1.json"
    if not backup.exists(): backup.write_bytes(config_path.read_bytes())
    config["framework_version"] = "1.1.0-dev"; config["multi_route_completion_required"] = True; config["required_annotation_views"] = ["strict","inclusive","display"]
    config_path.write_text(json.dumps(config, ensure_ascii=False, indent=2) + "\n")
    created=[]
    for filename, fields in REGISTRIES.items():
        path=state/filename
        if not path.exists():
            with path.open("w",newline="",encoding="utf-8") as handle:csv.writer(handle,delimiter="\t").writerow(fields)
            created.append(filename)
    cell = state/"cell_ledger.tsv.gz" if (state/"cell_ledger.tsv.gz").exists() else state/"cell_ledger.tsv"
    result={"status":"MIGRATED_REQUIRES_REVIEW","source_config_sha256":sha(backup),"created_registries":created,"cell_ledger":str(cell) if cell.exists() else None,"required_next_actions":["audit_multiroute_state","populate route attempts from validated artifacts only","reopen gaps with new decision version","migrate to v1.5 and build one final annotation","rerun completion and release audit"]}
    (provenance/"migration_v1_to_v1_1.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2));return 0


if __name__=="__main__":raise SystemExit(main())
