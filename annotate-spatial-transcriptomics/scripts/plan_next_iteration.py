#!/usr/bin/env python3
"""Build the next-action queue for the direct cohort annotation architecture."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

from validate_direct_lineage_workflow import audit as audit_direct
from validate_profile_role import load_profile
from validate_lineage_signal_coverage import audit as audit_lineage_signals
from validate_project_input_boundary import validate as validate_input_boundary
from validate_broad_class_completeness import validate as validate_broad_completeness


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def active_decisions(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    def decision_id(row: dict[str, str]) -> str:
        return row.get("decision_id") or f"{row.get('decision_version','')}:{row.get('source_run_id','')}:{row.get('source_cluster','')}"

    superseded: set[str] = set()
    for row in rows:
        superseded.update(item for item in re.split(r"[;,\s]+", row.get("supersedes", "").strip()) if item)
    return [row for row in rows if decision_id(row) not in superseded]


PRIORITY_BY_ENTITY = {
    "analysis_set": 0, "cell": 0, "broad_label": 1, "cohort": 1,
    "direct_return": 1, "targeted_cohort": 1, "terminal_residual_qc": 2,
    "route": 2, "annotation_view": 4,
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--context", required=True, type=Path)
    parser.add_argument("--biological-profile", type=Path)
    parser.add_argument("--profile", type=Path, help="deprecated biological-profile alias")
    parser.add_argument("--strict-exit-code", action="store_true", help="return nonzero for expected iteration work")
    args = parser.parse_args()
    root = args.project_root.resolve()
    context = json.loads(args.context.read_text(encoding="utf-8"))
    profile_path = args.biological_profile or args.profile
    profile = load_profile(profile_path, "biological_evidence") if profile_path else {}
    project_path = root / "config/project.json"
    project = json.loads(project_path.read_text(encoding="utf-8")) if project_path.is_file() else {}

    if project.get("routing_model", "direct_cross_lineage_recluster_cohorts") not in {"direct_cross_lineage_recluster_cohorts", "direct_cross_branch_recluster_cohorts", "direct_cross_lineage_recluster_cohorts_global_atlas"}:
        raise SystemExit("legacy pool project detected; migrate it before using the current planner")

    cluster_rows = active_decisions(read_tsv(root / "state/cluster_decision_ledger.tsv"))
    result = audit_direct(root)
    rows: list[dict[str, object]] = []

    if project.get("project_input_boundary_validation_required", False) is True:
        boundary = validate_input_boundary(root)
        for index, message in enumerate(boundary.get("errors", []), start=1):
            rows.append({
                "priority": 0, "source_run_id": "project", "source_cluster": "__PROJECT_INPUT_BOUNDARY__",
                "n_observations": 0, "current_state": "provenance_gap", "broad_label": "", "fine_label": "",
                "required_route": "register_or_repair_project_local_expression_ancestry", "reason": message,
                "target_scope": "all_query_expression_artifacts", "blocked_until": "project input boundary validation PASS",
                "gap_code": "PROJECT_INPUT_BOUNDARY_REQUIRED", "entity_type": "project", "entity_id": f"input_boundary_{index}",
            })

    if project.get("continuous_open_world_lineage_scan_required", False) is True:
        signal_audit = audit_lineage_signals(root)
        for index, message in enumerate(signal_audit.get("errors", []), start=1):
            rows.append({
                "priority": 0, "source_run_id": "project", "source_cluster": "__LINEAGE_SIGNAL_COVERAGE__",
                "n_observations": 0, "current_state": "lineage_signal_gap", "broad_label": "", "fine_label": "",
                "required_route": "continuous_open_world_lineage_scan", "reason": message,
                "target_scope": "whole_tissue_and_every_recluster_boundary",
                "blocked_until": "full catalog-by-cluster scan and explicit signal resolution",
                "gap_code": "LINEAGE_SIGNAL_COVERAGE_REQUIRED", "entity_type": "lineage_signal", "entity_id": f"signal_gap_{index}",
            })

    if profile.get("final_validation"):
        full_feature = root / "provenance/full_feature_validation.json"
        passed = False
        if full_feature.is_file():
            try:
                passed = json.loads(full_feature.read_text(encoding="utf-8")).get("status") == "PASS"
            except json.JSONDecodeError:
                passed = False
        if not passed:
            rows.append({
                "priority": 0, "source_run_id": "project", "source_cluster": "__FULL_FEATURE_VALIDATION__",
                "n_observations": 0, "current_state": "project_gate", "broad_label": "", "fine_label": "",
                "required_route": "full_feature_marker_validation",
                "reason": "final marker/anti-marker and context-specific evidence require a validated full-feature layer",
                "target_scope": "whole_tissue_and_all_recluster_cohorts",
                "blocked_until": "full-feature audit PASS and evidence writeback",
                "gap_code": "FULL_FEATURE_VALIDATION_REQUIRED", "entity_type": "project", "entity_id": "full_feature_validation",
            })

    atlas_rows = [row for row in read_tsv(root / "state/route_attempt_registry.tsv") if row.get("route_class") == "global_atlas_broad_audit" and row.get("status") in {"completed", "validated_done", "closed"}]
    if atlas_rows and project.get("broad_class_completeness_review_required", False) is True:
        catalog = Path(__file__).resolve().parents[1] / "references/profiles/sheep_ovary_candidate_lineage_catalog.json"
        completeness = validate_broad_completeness(root, catalog)
        for index, message in enumerate(completeness.get("errors", []), start=1):
            rows.append({
                "priority": 1, "source_run_id": "project", "source_cluster": "__BROAD_CLASS_COMPLETENESS__",
                "n_observations": 0, "current_state": "post_atlas_review_gap", "broad_label": "", "fine_label": "",
                "required_route": "broad_class_completeness_review", "reason": message,
                "target_scope": "present_and_zero_census_broad_lineages", "blocked_until": "query-derived broad completeness validation PASS",
                "gap_code": "BROAD_CLASS_COMPLETENESS_REQUIRED", "entity_type": "project", "entity_id": f"broad_completeness_{index}",
            })

    seen: set[tuple[str, str]] = set()
    for workflow_gap in result.get("gaps", []):
        if not workflow_gap.get("blocking", True):
            continue
        code = workflow_gap.get("code", "UNSTRUCTURED_GAP")
        entity_type = workflow_gap.get("entity_type", "project")
        target = workflow_gap.get("entity_id", "project")
        route = workflow_gap.get("required_action", "repair_annotation_workflow_gap")
        message = workflow_gap.get("detail", code)
        priority = PRIORITY_BY_ENTITY.get(entity_type, 2)
        key = (code, target)
        if key in seen:
            continue
        seen.add(key)
        rows.append({
            "priority": priority, "source_run_id": "project", "source_cluster": target,
            "n_observations": 0, "current_state": "annotation_workflow_gap", "broad_label": target if route.startswith("broad_class") else "",
            "fine_label": "", "required_route": route, "reason": message,
            "target_scope": target, "blocked_until": route,
            "gap_code": code, "entity_type": entity_type, "entity_id": target,
        })

    for row in cluster_rows:
        if row.get("state") in {"defined_fine", "defined_broad_only"} and row.get("validation_feature_scope") != "full_feature" and profile.get("final_validation"):
            key = ("full_feature_marker_validation", row.get("source_cluster", ""))
            if key not in seen:
                seen.add(key)
                rows.append({
                    "priority": 1, "source_run_id": row.get("source_run_id", ""), "source_cluster": row.get("source_cluster", ""),
                    "n_observations": int(float(row.get("n_observations", 0) or 0)), "current_state": row.get("state", ""),
                    "broad_label": row.get("broad_label", ""), "fine_label": row.get("fine_label", ""),
                    "required_route": "full_feature_marker_validation", "reason": "biological decision lacks full-feature validation",
                    "target_scope": row.get("recluster_cohort_id", ""), "blocked_until": "full-feature positive/anti-marker evidence",
                    "gap_code": "DECISION_FULL_FEATURE_EVIDENCE_REQUIRED", "entity_type": "cluster_decision", "entity_id": row.get("decision_id", row.get("source_cluster", "")),
                })

    rows.sort(key=lambda row: (int(row["priority"]), -int(row["n_observations"])))
    queue_path = root / "state/next_action_queue.tsv"
    queue_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "priority", "source_run_id", "source_cluster", "n_observations", "current_state", "broad_label",
        "fine_label", "required_route", "reason", "target_scope", "blocked_until",
        "gap_code", "entity_type", "entity_id",
    ]
    with queue_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    plan = {
        "status": "ITERATION_REQUIRED" if rows else "READY_FOR_COMPLETION_AUDIT",
        "workflow_model": result.get("workflow_model"),
        "queued_actions": len(rows),
        "queued_observations": sum(int(row["n_observations"]) for row in rows),
        "active_decisions": len(cluster_rows),
        "workflow_audit_status": result.get("status"),
        "context": context,
    }
    provenance = root / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    (provenance / "direct_lineage_workflow_audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (provenance / "iteration_plan.json").write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(plan, ensure_ascii=False, indent=2))
    return 2 if rows and args.strict_exit_code else 0


if __name__ == "__main__":
    raise SystemExit(main())
