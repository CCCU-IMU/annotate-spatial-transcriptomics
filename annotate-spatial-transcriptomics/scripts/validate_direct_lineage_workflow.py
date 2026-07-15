#!/usr/bin/env python3
"""Validate the default pool-free broad→subtype→QC-Atlas workflow."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path


GRID = [0.1, 0.2, 0.3, 0.4, 0.6]
TERMINAL = {"validated_done", "not_applicable_reviewed"}
MODERATE = {"moderate", "moderate_only", "moderate_high", "medium", "medium_high", "high", "high_confidence", "very_high"}
HIGH = {"high", "high_confidence", "very_high"}
RETAINED = {"qc_holdout", "qc_reject", "technical_state", "excluded_initial_qc"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def path_at(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def truth(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "pass", "passed"}


def confidence(value: str) -> str:
    return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")


def resolutions(value: str) -> list[float]:
    try:
        return [float(item) for item in re.split(r"[,;\s]+", value.strip()) if item]
    except ValueError:
        return []


def validate_membership(root: Path, row: dict[str, str], prefix: str = "membership") -> list[str]:
    errors: list[str] = []
    value = row.get(f"{prefix}_path", "") or row.get(f"{prefix}_artifact", "")
    expected = row.get(f"{prefix}_sha256", "")
    path = path_at(root, value) if value else Path()
    if not value or not path.is_file():
        return [f"missing {prefix} artifact"]
    if not expected or sha256(path) != expected:
        errors.append(f"stale {prefix} SHA256")
    members = read_tsv(path)
    ids = [row.get("cell_id", "") for row in members]
    if not members or "cell_id" not in members[0] or "" in ids or len(ids) != len(set(ids)):
        errors.append(f"{prefix} is not a unique nonempty cell_id table")
    try:
        expected_n = int(float(row.get("n_observations", row.get("n_query", 0)) or 0))
        if expected_n != len(ids):
            errors.append(f"{prefix} count differs from registry")
    except ValueError:
        errors.append(f"invalid {prefix} count")
    return errors


def audit(root: Path) -> dict:
    errors: list[str] = []
    project_path = root / "config/project.json"
    context_path = root / "config/biological_context.json"
    project = json.loads(project_path.read_text(encoding="utf-8")) if project_path.is_file() else {}
    context = json.loads(context_path.read_text(encoding="utf-8")) if context_path.is_file() else {}
    species = str(context.get("species", "")).lower()
    tissue = str(context.get("tissue", "")).lower()
    sheep_ovary = (
        project.get("strategy_preset_requested") == "sheep_ovary_same_batch_rfirst"
        or (any(token in species for token in ("sheep", "ovis", "ovine", "羊")) and any(token in tissue for token in ("ovary", "ovarian", "卵巢")))
    )
    cell_path = root / "state/cell_ledger.tsv.gz"
    if not cell_path.exists():
        cell_path = root / "state/cell_ledger.tsv"
    cells = read_tsv(cell_path)
    cohorts = read_tsv(root / "state/recluster_cohort_registry.tsv")
    returns = read_tsv(root / "state/direct_return_registry.tsv")
    routes = read_tsv(root / "state/route_attempt_registry.tsv")
    views = read_tsv(root / "state/annotation_view_registry.tsv")
    if not cells:
        errors.append("cell ledger is missing or empty")
    required_cell_fields = {
        "cell_id", "analysis_scope", "initial_broad_label", "assignment_mode", "final_state",
        "final_broad_label", "final_fine_label", "fine_anchor_eligible", "recluster_cohort_id", "cross_lineage_target",
    }
    if cells and not required_cell_fields.issubset(cells[0]):
        errors.append("cell ledger lacks direct-workflow fields: " + ", ".join(sorted(required_cell_fields - set(cells[0]))))
    analysis = [row for row in cells if row.get("analysis_scope") == "analysis_set"]
    initial_broad = sorted({row.get("initial_broad_label", "").strip() for row in analysis if row.get("initial_broad_label", "").strip()})
    by_broad: dict[str, list[dict[str, str]]] = {}
    for row in cohorts:
        if row.get("cohort_type") == "broad_class_recluster":
            by_broad.setdefault(row.get("source_broad_label", ""), []).append(row)
        if row.get("cohort_type") not in {"broad_class_recluster", "targeted_recluster", "oocyte_targeted_recluster"}:
            errors.append(f"unknown recluster cohort type: {row.get('cohort_id','')}")
        if row.get("status") not in TERMINAL:
            errors.append(f"recluster cohort is nonterminal: {row.get('cohort_id','')}")
        if row.get("status") == "not_applicable_reviewed":
            if len(row.get("applicability_rationale", "").strip()) < 20:
                errors.append(f"small/nonreclustered broad class lacks rationale: {row.get('cohort_id','')}")
        else:
            errors.extend(f"cohort {row.get('cohort_id','')}: {item}" for item in validate_membership(root, row))
            candidate_grid = resolutions(row.get("candidate_resolutions", ""))
            if sheep_ovary and candidate_grid != GRID:
                errors.append(f"cohort {row.get('cohort_id','')} does not use the formal sheep-ovary resolution grid")
            if not sheep_ovary and not candidate_grid:
                errors.append(f"cohort {row.get('cohort_id','')} lacks a candidate resolution grid")
            try:
                if float(row.get("selected_resolution", "")) not in candidate_grid:
                    errors.append(f"cohort {row.get('cohort_id','')} selected resolution outside its candidate grid")
            except ValueError:
                errors.append(f"cohort {row.get('cohort_id','')} lacks a selected resolution")
            outcome = path_at(root, row.get("outcome_artifact", "")) if row.get("outcome_artifact") else Path()
            if not outcome.is_file():
                errors.append(f"cohort {row.get('cohort_id','')} lacks a validated outcome artifact")
    for label in initial_broad:
        terminal = [row for row in by_broad.get(label, []) if row.get("status") in TERMINAL]
        if not terminal:
            errors.append(f"initial broad class lacks a terminal recluster/skip audit: {label}")
    for row in returns:
        rid = row.get("return_id", "")
        if row.get("status") != "validated_done":
            errors.append(f"direct return is nonterminal: {rid}")
        errors.extend(f"direct return {rid}: {item}" for item in validate_membership(root, row))
        if not row.get("target_broad_label", "").strip():
            errors.append(f"direct return lacks target broad label: {rid}")
        conf = confidence(row.get("confidence", ""))
        fine = row.get("target_fine_label", "").strip()
        mode = row.get("assignment_mode", "")
        if conf not in MODERATE:
            errors.append(f"direct broad return is below moderate confidence: {rid}")
        if fine and conf not in HIGH:
            errors.append(f"direct fine return is below high confidence: {rid}")
        evidence = path_at(root, row.get("evidence_artifact", "")) if row.get("evidence_artifact") else Path()
        if not evidence.is_file():
            errors.append(f"direct return lacks evidence artifact: {rid}")
        if mode == "rctd_assisted":
            if fine and (row.get("rctd_tier") != "extreme" or not truth(row.get("independent_evidence", ""))):
                errors.append(f"RCTD-assisted fine return lacks extreme tier plus independent evidence: {rid}")
            if not fine and row.get("rctd_tier") not in {"extreme", "high"}:
                errors.append(f"RCTD broad return is below high tier: {rid}")
            if truth(row.get("fine_anchor_eligible", "")):
                errors.append(f"RCTD-assisted return cannot become a fine anchor: {rid}")
    for index, row in enumerate(analysis, 2):
        broad = row.get("final_broad_label", "").strip()
        fine = row.get("final_fine_label", "").strip()
        broad_conf = confidence(row.get("final_broad_confidence") or row.get("final_confidence", ""))
        fine_conf = confidence(row.get("final_fine_confidence") or row.get("final_confidence", ""))
        if broad and broad_conf not in MODERATE:
            errors.append(f"cell row {index} broad label is below moderate confidence")
        if fine and (not broad or fine_conf not in HIGH or not truth(row.get("fine_anchor_eligible", ""))):
            errors.append(f"cell row {index} fine label lacks a high-confidence eligible parent")
        if row.get("assignment_mode") == "atlas_broad_rescue":
            if fine or truth(row.get("fine_anchor_eligible", "")):
                errors.append(f"cell row {index} Atlas rescue leaked into a fine label/anchor")
        if not broad and row.get("final_state") not in RETAINED:
            errors.append(f"cell row {index} has neither a biological broad label nor an explicit retained state")
    atlas = [row for row in routes if row.get("route_class") == "residual_qc_atlas_review"]
    active_atlas = [row for row in atlas if row.get("status") in TERMINAL]
    if not active_atlas:
        errors.append("residual QC Atlas phase lacks an applicable or zero-query terminal audit")
    for row in active_atlas:
        aid = row.get("route_attempt_id", "")
        if row.get("applicability") == "not_applicable":
            if row.get("status") != "not_applicable_reviewed" or len(row.get("applicability_rationale", "").strip()) < 20:
                errors.append(f"zero-QC Atlas audit lacks rationale: {aid}")
            continue
        if row.get("source_state") != "residual_qc_holdout_after_all_cohorts":
            errors.append(f"Atlas source is not final residual QC: {aid}")
        reference_id = row.get("reference_id", "")
        if sheep_ovary and reference_id != "GSE233801" and not reference_id.startswith("matched_count_level:"):
            errors.append(f"Atlas reference is not the declared sheep priority/matched reference: {aid}")
        if not sheep_ovary and not reference_id:
            errors.append(f"Atlas route lacks a declared reference: {aid}")
        if row.get("calibration_origin") != "query_like_heldout_current_query_anchors":
            errors.append(f"Atlas calibration origin is invalid: {aid}")
        if not truth(row.get("depth_matched_validation", "")) or not truth(row.get("observed_density_spatial_prior", "")):
            errors.append(f"Atlas lacks depth/spatial validation: {aid}")
        errors.extend(f"Atlas {aid}: {item}" for item in validate_membership(root, row, "query_membership"))
        try:
            n_query = int(float(row.get("n_query", 0) or 0))
            n_broad = int(float(row.get("n_defined_broad_only", 0) or 0))
            n_qc = int(float(row.get("n_qc_retained", 0) or 0))
            if n_broad + n_qc != n_query or int(float(row.get("n_defined_fine", 0) or 0)) != 0:
                errors.append(f"Atlas broad-return/QC outcomes do not partition the query: {aid}")
        except ValueError:
            errors.append(f"Atlas route counts are invalid: {aid}")
        if truth(row.get("fine_anchor_eligible", "")):
            errors.append(f"Atlas route is fine-anchor eligible: {aid}")
        for key in ("calibration_manifest", "validation_artifact", "outcome_artifact"):
            value = row.get(key, "")
            if not value or not path_at(root, value).is_file():
                errors.append(f"Atlas route lacks {key}: {aid}")
    if not any(row.get("view") == "final" and row.get("status") == "validated" for row in views):
        errors.append("single final annotation view is not validated")
    return {
        "status": "PASS" if not errors else "BLOCKED",
        "workflow_model": "broad_then_subtype_direct_cross_lineage_then_residual_qc_atlas",
        "profile_mode": "sheep_ovary" if sheep_ovary else "generic_context",
        "initial_broad_labels": initial_broad,
        "recluster_cohort_n": len(cohorts),
        "direct_return_n": len(returns),
        "atlas_phase_n": len(active_atlas),
        "persistent_biological_pools": False,
        "errors": errors,
        "gaps": [{"reason": item} for item in errors],
        "invalid_attempts": [],
        "missing_views": [] if any(row.get("view") == "final" and row.get("status") == "validated" for row in views) else ["final"],
    }


def main() -> int:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = audit(args.project_root.resolve())
    out = args.out or args.project_root / "provenance/direct_lineage_workflow_audit.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
