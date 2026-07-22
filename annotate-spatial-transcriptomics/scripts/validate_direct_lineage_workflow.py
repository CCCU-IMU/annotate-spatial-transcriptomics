#!/usr/bin/env python3
"""Validate the active broad-cohort → direct-return → terminal-QC workflow."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from evidence_schema_lib import active_registry_rows, path_at, read_tsv, sha256, write_result
from validate_cohort_outcome import validate as validate_cohort_outcome
from validate_direct_return_evidence import validate as validate_direct_return_evidence
from validate_prelabel_broad_evidence import validate as validate_prelabel_broad_evidence
from workflow_contract_lib import active_workflow_contract


FORMAL_GRID = [0.1, 0.2, 0.3, 0.4, 0.6]
TERMINAL = {"validated_done", "not_applicable_reviewed"}
CANONICAL_CONFIDENCE = {"low", "moderate", "high"}
ATLAS_CONFIDENCE = {"high", "moderate_only", "low_reject"}
RETAINED = {"qc_holdout", "qc_reject", "technical_state", "unknown_candidate", "excluded_initial_qc"}


def truth(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "pass", "passed"}


def confidence(value: str) -> str:
    return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")


def resolutions(value: str) -> list[float]:
    try:
        return [float(item) for item in re.split(r"[,;\s]+", value.strip()) if item]
    except ValueError:
        return []


def gap(code: str, entity_type: str, entity_id: str, required_action: str, detail: str, blocking: bool = True) -> dict:
    return {
        "code": code,
        "entity_type": entity_type,
        "entity_id": entity_id,
        "required_action": required_action,
        "blocking": blocking,
        "detail": detail,
    }


def audit(root: Path) -> dict:
    gaps: list[dict] = []

    def add(code: str, entity_type: str, entity_id: str, action: str, detail: str, blocking: bool = True) -> None:
        gaps.append(gap(code, entity_type, entity_id, action, detail, blocking))

    contract = active_workflow_contract(root)
    project_path = root / "config/project.json"
    try:
        project = json.loads(project_path.read_text(encoding="utf-8")) if project_path.is_file() else {}
    except json.JSONDecodeError:
        project = {}
    global_atlas_required = project.get("routing_model") == "direct_cross_lineage_recluster_cohorts_global_atlas"
    cell_path = root / "state/cell_ledger.tsv.gz"
    if not cell_path.exists():
        cell_path = root / "state/cell_ledger.tsv"
    cells = read_tsv(cell_path)
    cohorts = active_registry_rows(read_tsv(root / "state/recluster_cohort_registry.tsv"), "cohort_id")
    returns = active_registry_rows(read_tsv(root / "state/direct_return_registry.tsv"), "return_id")
    routes = active_registry_rows(read_tsv(root / "state/route_attempt_registry.tsv"), "route_attempt_id")
    views = read_tsv(root / "state/annotation_view_registry.tsv")
    if not cells:
        add("CELL_LEDGER_MISSING", "analysis_set", "analysis_set", "write_or_repair_cell_ledger", "cell ledger is missing or empty")
    required_fields = {
        "cell_id", "analysis_scope", "initial_broad_label", "assignment_mode", "final_state",
        "final_broad_label", "final_fine_label", "final_confidence", "fine_anchor_eligible",
        "recluster_cohort_id", "cross_lineage_target",
    }
    if cells and not required_fields.issubset(cells[0]):
        add("CELL_LEDGER_SCHEMA_INVALID", "analysis_set", "analysis_set", "repair_cell_ledger_schema", "missing fields: " + ", ".join(sorted(required_fields - set(cells[0]))))
    analysis = [row for row in cells if row.get("analysis_scope") == "analysis_set"]
    initial_broad = sorted({row.get("initial_broad_label", "").strip() for row in analysis if row.get("initial_broad_label", "").strip()})
    if project.get("prelabel_evidence_freeze_required"):
        prelabel = validate_prelabel_broad_evidence(root, root / "state/cluster_decision_ledger.tsv")
        write_result(root / "provenance/prelabel_broad_evidence_validation.json", prelabel)
        for message in prelabel.get("errors", []):
            add("PRELABEL_EVIDENCE_INVALID", "whole_tissue", "prelabel_broad_evidence", "freeze_label_blind_candidate_matrix_before_initial_labels", message)
    by_broad: dict[str, list[dict[str, str]]] = {}
    for cohort in cohorts:
        cohort_id = cohort.get("cohort_id", "")
        cohort_type = cohort.get("cohort_type", "")
        if cohort_type == "broad_class_recluster":
            by_broad.setdefault(cohort.get("source_broad_label", ""), []).append(cohort)
            if cohort.get("question_mode") != "broad_purity_audit":
                add("BROAD_COHORT_QUESTION_MODE_INVALID", "cohort", cohort_id, "set_broad_purity_audit", "broad-class cohort must not be forced to split")
        elif cohort_type in {"targeted_recluster", "oocyte_targeted_recluster"}:
            if cohort.get("question_mode") != "targeted_mixture":
                add("TARGETED_COHORT_QUESTION_MODE_INVALID", "cohort", cohort_id, "set_targeted_mixture", "targeted cohort must state competing hypotheses")
        else:
            add("COHORT_TYPE_INVALID", "cohort", cohort_id, "repair_cohort_registry", f"unknown cohort type {cohort_type!r}")
        if "qc" in cohort_type.lower() or cohort.get("source_broad_label", "").strip().lower() in {"qc", "qc_holdout", "qc_reject"}:
            add("TERMINAL_QC_RECLUSTER_FORBIDDEN", "cohort", cohort_id, "remove_qc_reclustering_and_freeze_terminal_qc", "terminal residual QC must never be a reclustering cohort")
        if cohort.get("status") not in TERMINAL:
            add("COHORT_NONTERMINAL", "cohort", cohort_id, "finish_or_repair_cohort", "cohort is not terminal")
            continue
        if cohort.get("status") == "not_applicable_reviewed":
            skip_value = cohort.get("underpowered_skip_artifact", "")
            skip_path = path_at(root, skip_value) if skip_value else Path()
            if len(cohort.get("applicability_rationale", "").strip()) < 20 or not skip_path.is_file() or not cohort.get("underpowered_skip_artifact_sha256") or sha256(skip_path) != cohort.get("underpowered_skip_artifact_sha256"):
                add("UNDERPOWERED_SKIP_EVIDENCE_INVALID", "cohort", cohort_id, "record_hash_bound_underpowered_skip", "reviewed skip lacks substantive, current evidence")
            continue
        candidate_grid = resolutions(cohort.get("candidate_resolutions", ""))
        expected_cohort_grid = [float(value) for value in contract.get("cohort_resolution_grid", [])]
        if expected_cohort_grid and candidate_grid != expected_cohort_grid:
            add("FORMAL_GRID_INCOMPLETE" if contract["same_batch_stereopy_cellbin_rfirst"] else "ACTIVE_GRID_INCOMPLETE", "cohort", cohort_id, "run_complete_active_contract_grid", f"cohort grid must equal the active workflow grid {expected_cohort_grid}")
        elif not candidate_grid:
            add("CANDIDATE_GRID_MISSING", "cohort", cohort_id, "run_workflow_profile_grid", "cohort lacks a candidate resolution grid")
        try:
            selected = float(cohort.get("selected_resolution", ""))
        except ValueError:
            selected = None
        if selected is None or selected not in candidate_grid:
            add("SELECTED_RESOLUTION_INVALID", "cohort", cohort_id, "select_integrated_evidence_optimum", "selected resolution is missing or outside the complete grid")
        terminal_outcome = cohort.get("terminal_outcome", "")
        if terminal_outcome not in {"homogeneous_parent_confirmed", "subclusters_adjudicated"}:
            add("COHORT_TERMINAL_OUTCOME_MISSING", "cohort", cohort_id, "adjudicate_cohort_terminal_outcome", "cohort lacks a formal biological endpoint")
        if terminal_outcome == "homogeneous_parent_confirmed" and cohort.get("question_mode") != "broad_purity_audit":
            add("HOMOGENEOUS_ENDPOINT_INVALID", "cohort", cohort_id, "repair_terminal_outcome", "homogeneous parent is only a broad-purity endpoint")
        outcome_value = cohort.get("outcome_artifact", "")
        outcome_path = path_at(root, outcome_value) if outcome_value else Path()
        if not outcome_path.is_file():
            add("COHORT_OUTCOME_MISSING", "cohort", cohort_id, "create_evidence_complete_cohort_outcome", "cohort outcome artifact is missing")
        else:
            validation = validate_cohort_outcome(root, outcome_path, root / "state/recluster_cohort_registry.tsv")
            for message in validation.get("errors", []):
                add("COHORT_OUTCOME_SCHEMA_INVALID", "cohort", cohort_id, "repair_cohort_outcome_evidence", message)

    for label in initial_broad:
        terminal = [row for row in by_broad.get(label, []) if row.get("status") in TERMINAL]
        if len(terminal) != 1:
            add("INITIAL_BROAD_COHORT_MISSING", "broad_label", label, "run_one_broad_class_recluster_or_formal_underpowered_skip", "each initial broad class requires exactly one terminal cohort")

    if contract["same_batch_stereopy_cellbin_rfirst"]:
        whole_grid_path = root / "provenance/whole_tissue_resolution_grid_validation.json"
        try:
            whole_grid = json.loads(whole_grid_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            whole_grid = {}
        embedded = whole_grid.get("active_workflow_contract", {}) if isinstance(whole_grid.get("active_workflow_contract"), dict) else {}
        if (
            whole_grid.get("status") != "PASS"
            or whole_grid.get("scope") != "whole_tissue"
            or whole_grid.get("candidate_resolutions") != FORMAL_GRID
            or embedded.get("workflow_profile_sha256") != contract.get("workflow_profile_sha256")
            or embedded.get("strategy_preset_id") != contract.get("strategy_preset_id")
        ):
            add("WHOLE_TISSUE_GRID_AUDIT_MISSING", "whole_tissue", "whole_tissue", "validate_complete_active_workflow_grid", "verified same-batch StereoPy cellbin project lacks a current whole-tissue 0.1,0.2,0.3,0.4,0.6 grid audit")

    for direct in returns:
        return_id = direct.get("return_id", "")
        if direct.get("status") != "validated_done":
            add("DIRECT_RETURN_NONTERMINAL", "direct_return", return_id, "finish_or_repair_direct_return", "direct return is nonterminal")
        conf = confidence(direct.get("confidence", ""))
        fine = direct.get("target_fine_label", "").strip()
        if conf not in {"moderate", "high"}:
            add("CONFIDENCE_NONCANONICAL", "direct_return", return_id, "migrate_confidence_to_canonical_enum", f"new direct return confidence is {conf!r}")
        if fine and (conf != "high" or not truth(direct.get("fine_anchor_eligible", ""))):
            add("FINE_RETURN_EVIDENCE_INSUFFICIENT", "direct_return", return_id, "return_broad_only_or_supply_high_fine_evidence", "fine return must be high and fine-anchor eligible")
        evidence_value = direct.get("evidence_artifact", "")
        evidence_path = path_at(root, evidence_value) if evidence_value else Path()
        if not evidence_path.is_file():
            add("DIRECT_RETURN_EVIDENCE_MISSING", "direct_return", return_id, "create_direct_return_evidence", "evidence artifact is missing")
        else:
            validation = validate_direct_return_evidence(root, evidence_path, root / "state/direct_return_registry.tsv")
            for message in validation.get("errors", []):
                add("DIRECT_RETURN_EVIDENCE_SCHEMA_INVALID", "direct_return", return_id, "repair_direct_return_evidence", message)
        if direct.get("assignment_mode") == "rctd_assisted":
            if fine and (direct.get("rctd_tier") != "high" or not truth(direct.get("independent_evidence", ""))):
                add("RCTD_FINE_GATE_FAILED", "direct_return", return_id, "return_broad_or_qc", "RCTD fine requires canonical high plus independent evidence")
            if not fine and direct.get("rctd_tier") not in {"high", "moderate"}:
                add("RCTD_BROAD_GATE_FAILED", "direct_return", return_id, "route_to_terminal_qc", "RCTD broad return is below canonical moderate")

    for row_index, cell in enumerate(analysis, 2):
        broad = cell.get("final_broad_label", "").strip()
        fine = cell.get("final_fine_label", "").strip()
        conf = confidence(cell.get("final_confidence", ""))
        if conf and conf not in CANONICAL_CONFIDENCE:
            add("CONFIDENCE_NONCANONICAL", "cell", cell.get("cell_id", str(row_index)), "migrate_confidence_to_canonical_enum", f"cell confidence is {conf!r}")
        if broad and conf not in {"moderate", "high"}:
            add("FINAL_BROAD_BELOW_MODERATE", "cell", cell.get("cell_id", str(row_index)), "return_to_review_or_qc", "final broad label lacks moderate-or-high confidence")
        if fine and (not broad or conf != "high" or not truth(cell.get("fine_anchor_eligible", ""))):
            add("FINAL_FINE_GATE_FAILED", "cell", cell.get("cell_id", str(row_index)), "collapse_to_broad_or_supply_high_fine_evidence", "final fine label lacks high fine-anchor evidence")
        if cell.get("assignment_mode") == "atlas_broad_rescue" and (fine or truth(cell.get("fine_anchor_eligible", ""))):
            add("ATLAS_FINE_LEAK", "cell", cell.get("cell_id", str(row_index)), "remove_atlas_fine_label", "Atlas rescue is broad-only")
        if not broad and cell.get("final_state") not in RETAINED:
            add("FINAL_CELL_UNCOVERED", "cell", cell.get("cell_id", str(row_index)), "write_unique_final_broad_or_retained_state", "cell has no final biological or retained state")

    for route in routes:
        if route.get("route_class") in {"qc_anchor_recluster", "qc_holdout_recluster"}:
            add("TERMINAL_QC_RECLUSTER_FORBIDDEN", "route", route.get("route_attempt_id", ""), "remove_qc_reclustering_and_freeze_terminal_qc", "legacy QC reclustering route is active")
    global_atlas = [row for row in routes if row.get("route_class") == "global_atlas_broad_audit" and row.get("status") in TERMINAL]
    residual_atlas = [row for row in routes if row.get("route_class") == "residual_qc_atlas_review" and row.get("status") in TERMINAL]
    if global_atlas and residual_atlas:
        add("ATLAS_CONTROLLER_DUPLICATED", "route", "atlas", "retain_one_atlas_controller", "global and legacy residual-QC Atlas controllers cannot both be active")
    atlas = global_atlas if global_atlas_required else (global_atlas or residual_atlas)
    if len(atlas) != 1:
        action = "run_or_record_global_all_cell_atlas_audit" if global_atlas_required else "run_or_record_atlas_consensus"
        add("TERMINAL_ATLAS_AUDIT_MISSING", "terminal_residual_qc", "terminal_residual_qc", action, "exactly one terminal Atlas audit is required")
    for route in atlas:
        route_id = route.get("route_attempt_id", "")
        if route.get("applicability") == "not_applicable":
            if len(route.get("applicability_rationale", "").strip()) < 20:
                add("ATLAS_NOT_APPLICABLE_RATIONALE_MISSING", "route", route_id, "record_reference_applicability_audit", "not-applicable Atlas audit lacks a substantive rationale")
            continue
        is_global = route.get("route_class") == "global_atlas_broad_audit"
        expected_source = "analysis_set_after_terminal_qc_freeze" if is_global else "residual_qc_holdout_after_all_cohorts"
        if route.get("source_state") != expected_source:
            add("ATLAS_QUERY_SCOPE_INVALID", "route", route_id, "freeze_required_atlas_query_scope", f"Atlas source must be {expected_source}")
        tiers = {value for value in re.split(r"[|,;\s]+", route.get("atlas_confidence_enum", "")) if value}
        if tiers != ATLAS_CONFIDENCE:
            add("ATLAS_CONFIDENCE_ENUM_INVALID", "route", route_id, "use_high_moderate_only_low_reject", f"Atlas confidence enum is {sorted(tiers)}")
        if route.get("calibration_origin") != "query_like_heldout_current_query_anchors":
            add("ATLAS_CALIBRATION_INVALID", "route", route_id, "calibrate_on_disjoint_query_like_anchors", "Atlas calibration origin is invalid")
        if not truth(route.get("depth_matched_validation", "")) or (not is_global and not truth(route.get("observed_density_spatial_prior", ""))):
            add("ATLAS_CALIBRATION_OR_LEGACY_CONSENSUS_INCOMPLETE", "route", route_id, "complete_query_anchor_calibration_or_legacy_consensus", "Atlas lacks required query-anchor calibration or a legacy residual-QC route lacks its declared spatial prior")
        required_artifacts = ["calibration_manifest", "validation_artifact", "outcome_artifact"]
        if is_global:
            required_artifacts.extend(["parameters_artifact", "concordance_artifact", "cluster_concordance_artifact", "discrepancy_review_artifact"])
        for key in required_artifacts:
            value = route.get(key, "")
            if not value or not path_at(root, value).is_file():
                add("ATLAS_EVIDENCE_ARTIFACT_MISSING", "route", route_id, "repair_atlas_evidence_bundle", f"missing {key}")
        if is_global:
            validation_value = route.get("validation_artifact", "")
            try:
                validation = json.loads(path_at(root, validation_value).read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                validation = {}
            if validation.get("status") != "PASS":
                add("GLOBAL_ATLAS_REVIEW_OPEN", "route", route_id, "close_every_material_disagreement_once", "global Atlas concordance or discrepancy review has not passed")
        if truth(route.get("fine_anchor_eligible", "")):
            add("ATLAS_FINE_ANCHOR_FORBIDDEN", "route", route_id, "set_fine_anchor_false", "Atlas route cannot seed fine discovery")

    final_view = any(row.get("view") == "final" and row.get("status") == "validated" for row in views)
    if not final_view:
        add("FINAL_VIEW_MISSING", "annotation_view", "final", "build_and_validate_single_final_annotation", "single final annotation view is not validated")
    blocking = [item for item in gaps if item["blocking"]]
    return {
        "status": "PASS" if not blocking else "BLOCKED",
        "controller_state": "CONTINUE" if not blocking else "ITERATION_REQUIRED",
        "workflow_model": "broad_class_recluster_then_global_all_cell_atlas_concordance" if global_atlas_required else "broad_class_recluster_cohort_then_direct_return_then_terminal_residual_qc",
        "active_workflow_contract": contract,
        "initial_broad_labels": initial_broad,
        "recluster_cohort_n": len(cohorts),
        "direct_return_n": len(returns),
        "atlas_phase_n": len(atlas),
        "atlas_controller": "global_all_cell" if global_atlas else "legacy_residual_qc",
        "persistent_biological_pools": False,
        "errors": [item["detail"] for item in blocking],
        "gaps": gaps,
        "invalid_attempts": [item for item in gaps if item["code"].endswith("INVALID")],
        "missing_views": [] if final_view else ["final"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--strict-exit-code", action="store_true", help="return nonzero for expected workflow gaps (CI/completion use)")
    args = parser.parse_args()
    root = args.project_root.resolve()
    result = audit(root)
    write_result(args.out or root / "provenance/direct_lineage_workflow_audit.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if args.strict_exit_code and result["status"] != "PASS":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
