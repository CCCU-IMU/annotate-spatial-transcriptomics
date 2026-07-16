#!/usr/bin/env python3
"""Audit exact, disjoint cell membership closure across the direct-lineage workflow."""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from evidence_schema_lib import active_registry_rows, membership_ids, path_at, read_tsv, sha256, write_result
from validate_cohort_outcome import validate as validate_cohort_outcome
from validate_direct_return_evidence import validate as validate_direct_return_evidence


RESIDUAL_QC_STATES = {"qc_holdout", "qc_reject"}
RETAINED_STATES = RESIDUAL_QC_STATES | {"technical_state", "unknown_candidate"}
RETURN_OUTCOMES = {"parent_return", "fine_return", "cross_lineage_return"}


def registry_membership(root: Path, row: dict[str, str], prefix: str = "membership") -> tuple[set[str], list[str]]:
    try:
        n = int(float(row.get(f"{prefix}_n_observations", row.get("n_observations", row.get("n_query", 0))) or 0))
    except ValueError:
        n = -1
    ref = {
        "path": row.get(f"{prefix}_path", "") or row.get(f"{prefix}_artifact", ""),
        "sha256": row.get(f"{prefix}_sha256", ""),
        "n_observations": n,
    }
    return membership_ids(root, ref, prefix)


def audit(root: Path) -> dict:
    errors: list[str] = []
    cell_path = root / "state/cell_ledger.tsv.gz"
    if not cell_path.exists():
        cell_path = root / "state/cell_ledger.tsv"
    cells = read_tsv(cell_path)
    cohorts = active_registry_rows(read_tsv(root / "state/recluster_cohort_registry.tsv"), "cohort_id")
    returns = active_registry_rows(read_tsv(root / "state/direct_return_registry.tsv"), "return_id")
    routes = active_registry_rows(read_tsv(root / "state/route_attempt_registry.tsv"), "route_attempt_id")
    analysis = [row for row in cells if row.get("analysis_scope") == "analysis_set"]
    cell_ids = [row.get("cell_id", "") for row in analysis]
    if not analysis or "" in cell_ids or len(cell_ids) != len(set(cell_ids)):
        errors.append("analysis set must be a nonempty unique cell_id collection")
    cell_by_id = {row["cell_id"]: row for row in analysis if row.get("cell_id")}
    analysis_ids = set(cell_by_id)

    broad_members = defaultdict(set)
    for row in analysis:
        label = row.get("initial_broad_label", "").strip()
        if label:
            broad_members[label].add(row["cell_id"])
    cohort_by_id = {row.get("cohort_id", ""): row for row in cohorts if row.get("cohort_id")}
    outcome_by_cohort: dict[str, dict] = {}
    outcome_memberships: dict[tuple[str, str], set[str]] = {}
    targeted_sources: dict[str, set[str]] = {}

    for label, expected in broad_members.items():
        matches = [row for row in cohorts if row.get("cohort_type") == "broad_class_recluster" and row.get("source_broad_label") == label]
        if len(matches) != 1:
            errors.append(f"initial broad label {label!r} must have exactly one broad cohort")
            continue
        cohort = matches[0]
        if cohort.get("status") == "not_applicable_reviewed":
            if not cohort.get("underpowered_skip_artifact", ""):
                errors.append(f"underpowered broad cohort {cohort.get('cohort_id')} lacks formal skip artifact")
            continue
        actual, member_errors = registry_membership(root, cohort)
        errors.extend(f"cohort {cohort.get('cohort_id')}: {message}" for message in member_errors)
        if actual != expected:
            errors.append(f"initial broad membership differs from broad cohort query: {label}")

    for cohort in cohorts:
        if cohort.get("status") == "not_applicable_reviewed":
            continue
        cohort_id = cohort.get("cohort_id", "")
        outcome_value = cohort.get("outcome_artifact", "")
        outcome_path = path_at(root, outcome_value) if outcome_value else Path()
        if not outcome_path.is_file():
            errors.append(f"cohort {cohort_id} lacks outcome artifact")
            continue
        validation = validate_cohort_outcome(root, outcome_path, root / "state/recluster_cohort_registry.tsv")
        errors.extend(f"cohort {cohort_id}: {message}" for message in validation.get("errors", []))
        try:
            outcome = json.loads(outcome_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        outcome_by_cohort[cohort_id] = outcome
        for subcluster in outcome.get("subcluster_outcomes", []):
            members, member_errors = membership_ids(root, subcluster.get("membership", {}), f"{cohort_id}/{subcluster.get('subcluster_id')}")
            errors.extend(member_errors)
            key = (cohort_id, str(subcluster.get("subcluster_id", "")))
            outcome_memberships[key] = members
            if subcluster.get("outcome") == "targeted_successor":
                target = str(subcluster.get("target_cohort_id", ""))
                if target in targeted_sources:
                    errors.append(f"targeted cohort {target} has multiple source subclusters")
                targeted_sources[target] = members

    for cohort in cohorts:
        if cohort.get("cohort_type") not in {"targeted_recluster", "oocyte_targeted_recluster"} or cohort.get("status") == "not_applicable_reviewed":
            continue
        cohort_id = cohort.get("cohort_id", "")
        actual, member_errors = registry_membership(root, cohort)
        errors.extend(f"targeted cohort {cohort_id}: {message}" for message in member_errors)
        expected = targeted_sources.get(cohort_id)
        if expected is None:
            errors.append(f"targeted cohort {cohort_id} is not linked to one explicit source subcluster")
        elif actual != expected:
            errors.append(f"targeted cohort {cohort_id} membership differs from its source subcluster")

    seen_direct: set[str] = set()
    direct_by_source: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    direct_members: dict[str, set[str]] = {}
    for direct in returns:
        return_id = direct.get("return_id", "")
        members, member_errors = registry_membership(root, direct)
        errors.extend(f"direct return {return_id}: {message}" for message in member_errors)
        overlap = seen_direct & members
        if overlap:
            errors.append(f"direct returns overlap for {len(overlap)} observations")
        seen_direct |= members
        direct_members[return_id] = members
        source_key = (direct.get("source_cohort_id", ""), direct.get("source_cluster", ""))
        direct_by_source[source_key].append(direct)
        evidence_value = direct.get("evidence_artifact", "")
        evidence_path = path_at(root, evidence_value) if evidence_value else Path()
        if evidence_path.is_file():
            validation = validate_direct_return_evidence(root, evidence_path, root / "state/direct_return_registry.tsv")
            errors.extend(f"direct return {return_id}: {message}" for message in validation.get("errors", []))
        else:
            errors.append(f"direct return {return_id} lacks evidence artifact")
        for cell_id in members:
            cell = cell_by_id.get(cell_id)
            if cell is None:
                errors.append(f"direct return {return_id} contains a cell outside the analysis set")
                continue
            if cell.get("final_broad_label", "") != direct.get("target_broad_label", ""):
                errors.append(f"direct return {return_id} broad writeback disagrees with cell ledger")
                break
            if cell.get("final_fine_label", "") != direct.get("target_fine_label", ""):
                errors.append(f"direct return {return_id} fine writeback disagrees with cell ledger")
                break

    for cohort_id, outcome in outcome_by_cohort.items():
        for subcluster in outcome.get("subcluster_outcomes", []):
            if subcluster.get("outcome") not in RETURN_OUTCOMES:
                continue
            key = (cohort_id, str(subcluster.get("subcluster_id", "")))
            matches = direct_by_source.get(key, [])
            if len(matches) != 1:
                errors.append(f"return outcome {cohort_id}/{key[1]} must have exactly one direct-return record")
                continue
            if direct_members.get(matches[0].get("return_id", ""), set()) != outcome_memberships.get(key, set()):
                errors.append(f"direct-return membership differs from source outcome: {cohort_id}/{key[1]}")

    # Reconstruct the pre-Atlas terminal residual QC membership from the
    # current ledger. Accepted Atlas returns have a final broad label now, but
    # remain members of the frozen query through their assignment provenance.
    terminal_qc = {
        row["cell_id"] for row in analysis
        if row.get("assignment_mode", "") == "atlas_broad_rescue"
        or (
            not row.get("final_broad_label", "").strip()
            and (row.get("final_state", "") in RESIDUAL_QC_STATES or row.get("state", "") in RESIDUAL_QC_STATES)
        )
    }
    qc_successor_members = set().union(*[
        members for key, members in outcome_memberships.items()
        if any(
            str(item.get("subcluster_id", "")) == key[1] and item.get("outcome") == "qc_successor"
            for item in outcome_by_cohort.get(key[0], {}).get("subcluster_outcomes", [])
        )
    ]) if outcome_memberships else set()
    if not qc_successor_members <= terminal_qc:
        errors.append("one or more QC-successor observations are missing from terminal residual QC provenance")
    global_atlas_rows = [row for row in routes if row.get("route_class") == "global_atlas_broad_audit" and row.get("status") in {"validated_done", "not_applicable_reviewed"}]
    residual_atlas_rows = [row for row in routes if row.get("route_class") == "residual_qc_atlas_review" and row.get("status") in {"validated_done", "not_applicable_reviewed"}]
    project_path = root / "config/project.json"
    try:
        project = json.loads(project_path.read_text(encoding="utf-8")) if project_path.is_file() else {}
    except json.JSONDecodeError:
        project = {}
    global_required = project.get("routing_model") == "direct_cross_lineage_recluster_cohorts_global_atlas"
    if global_required and len(global_atlas_rows) != 1:
        errors.append("v1.7 workflow requires exactly one terminal global all-cell Atlas audit")
    if global_atlas_rows and residual_atlas_rows:
        errors.append("global and legacy residual-QC Atlas controllers cannot both be active")
    atlas_rows = global_atlas_rows or residual_atlas_rows
    if not global_required and len(atlas_rows) != 1:
        errors.append("terminal workflow must have exactly one Atlas audit")
    elif atlas_rows:
        atlas = atlas_rows[0]
        if atlas.get("status") == "not_applicable_reviewed":
            pass
        elif atlas.get("route_class") == "global_atlas_broad_audit":
            query, query_errors = registry_membership(root, atlas, "query_membership")
            errors.extend(f"global Atlas: {message}" for message in query_errors)
            if query != analysis_ids:
                errors.append("global Atlas query is not cell-for-cell identical to the analysis set")
            frozen_qc, qc_errors = registry_membership(root, atlas, "qc_membership")
            errors.extend(f"global Atlas QC: {message}" for message in qc_errors)
            if frozen_qc != terminal_qc:
                errors.append("global Atlas frozen QC membership differs from terminal residual QC")
            accepted, accepted_errors = registry_membership(root, atlas, "accepted_membership")
            rejected, rejected_errors = registry_membership(root, atlas, "rejected_membership")
            errors.extend(f"global Atlas accepted: {message}" for message in accepted_errors)
            errors.extend(f"global Atlas rejected: {message}" for message in rejected_errors)
            if accepted & rejected or accepted | rejected != frozen_qc:
                errors.append("global Atlas QC accepted and rejected memberships do not partition frozen QC")
            ledger_accepted = {row["cell_id"] for row in analysis if row.get("assignment_mode") == "atlas_broad_rescue"}
            ledger_rejected = frozen_qc - ledger_accepted
            if accepted != ledger_accepted or rejected != ledger_rejected:
                errors.append("global Atlas QC outcomes disagree with cell-ledger writeback")
        else:
            query, query_errors = registry_membership(root, atlas, "query_membership")
            errors.extend(f"Atlas: {message}" for message in query_errors)
            if query != terminal_qc:
                errors.append("Atlas query is not cell-for-cell identical to terminal residual QC")
            accepted, accepted_errors = registry_membership(root, atlas, "accepted_membership")
            rejected, rejected_errors = registry_membership(root, atlas, "rejected_membership")
            errors.extend(f"Atlas accepted: {message}" for message in accepted_errors)
            errors.extend(f"Atlas rejected: {message}" for message in rejected_errors)
            if accepted & rejected or accepted | rejected != query:
                errors.append("Atlas accepted and rejected memberships do not exactly partition the query")
            ledger_accepted = {row["cell_id"] for row in analysis if row.get("assignment_mode") == "atlas_broad_rescue"}
            ledger_rejected = {row["cell_id"] for row in analysis if row["cell_id"] in query and row["cell_id"] not in ledger_accepted}
            if accepted != ledger_accepted or rejected != ledger_rejected:
                errors.append("Atlas accepted/rejected results disagree with cell-ledger writeback")

    final_covered = {
        row["cell_id"] for row in analysis
        if bool(row.get("final_broad_label", "").strip()) ^ bool(row.get("final_state", "") in RETAINED_STATES)
    }
    if final_covered != analysis_ids:
        errors.append("final annotation does not uniquely cover every analysis-set observation")
    result = {
        "status": "PASS" if not errors else "FAIL",
        "cell_ledger": str(cell_path),
        "cell_ledger_sha256": sha256(cell_path) if cell_path.is_file() else None,
        "n_analysis_set": len(analysis_ids),
        "n_initial_broad_labels": len(broad_members),
        "n_cohorts": len(cohorts),
        "n_direct_returns": len(returns),
        "n_terminal_residual_qc": len(terminal_qc),
        "atlas_controller": "global_all_cell" if global_atlas_rows else "legacy_residual_qc",
        "errors": errors,
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    result = audit(root)
    write_result(args.out or root / "provenance/annotation_membership_partition_audit.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
