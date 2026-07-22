#!/usr/bin/env python3
"""Validate cohort biological evidence content and exact subcluster partition."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from evidence_schema_lib import (
    membership_ids, path_at, read_tsv, sha256, validate_evidence_artifact,
    validate_json_against_schema, write_result,
)
from lineage_decision_lib import (
    absolute_supported_families,
    cluster_rows as decision_cluster_rows,
    eligible as decision_eligible,
    purity_passes,
    rank as rank_decisions,
    split_families,
    validate_rows as validate_decision_rows,
)


def _project_config(root: Path) -> dict:
    try:
        value = json.loads((root / "config/project.json").read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, json.JSONDecodeError, TypeError):
        return {}


def _validate_subset_purity(root: Path, row: dict, members: set[str]) -> list[str]:
    errors: list[str] = []
    ref = row.get("observation_purity_evidence")
    if not isinstance(ref, dict):
        return [f"subcluster {row['subcluster_id']} supported subset lacks observation-level purity evidence"]
    path, artifact_errors = validate_evidence_artifact(root, ref, f"subcluster {row['subcluster_id']} subset purity evidence")
    errors.extend(artifact_errors)
    if not path or artifact_errors:
        return errors
    try:
        evidence = json.loads(path.read_text(encoding="utf-8"))
        required = {
            "status", "candidate_id", "membership_sha256", "n_observations",
            "positive_family_count", "contradiction_count", "lineage_supported_fraction",
            "strongest_competing_fraction",
        }
        if not isinstance(evidence, dict) or not required <= set(evidence):
            return errors + [f"subcluster {row['subcluster_id']} subset purity evidence lacks numerical content"]
        if evidence["status"] != "PASS" or evidence["candidate_id"] != row.get("target_candidate_id"):
            errors.append(f"subcluster {row['subcluster_id']} subset purity evidence did not pass for the returned candidate")
        if evidence["membership_sha256"] != row["membership"]["sha256"] or int(evidence["n_observations"]) != len(members):
            errors.append(f"subcluster {row['subcluster_id']} subset purity evidence is stale for its membership")
        if int(evidence["positive_family_count"]) < 2 or int(evidence["contradiction_count"]) != 0:
            errors.append(f"subcluster {row['subcluster_id']} subset return lacks two families or has a contradiction")
        if float(evidence["lineage_supported_fraction"]) <= float(evidence["strongest_competing_fraction"]):
            errors.append(f"subcluster {row['subcluster_id']} subset return is not purer than its strongest competitor")
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError):
        errors.append(f"subcluster {row['subcluster_id']} subset purity evidence is unreadable")
    return errors


def validate(root: Path, outcome_path: Path, registry_path: Path | None = None) -> dict:
    schema = Path(__file__).resolve().parents[1] / "schemas/cohort_outcome.schema.json"
    outcome, errors = validate_json_against_schema(outcome_path, schema)
    if errors:
        return {"status": "FAIL", "artifact": str(outcome_path), "errors": errors}
    query, query_errors = membership_ids(root, outcome["query_membership"], "query membership")
    errors.extend(query_errors)
    project = _project_config(root)
    try:
        framework = tuple(int(part) for part in str(project.get("framework_version", "0")).split(".")[:3])
    except ValueError:
        framework = (0,)
    evidence_v2 = outcome.get("schema_version") == "2.0"
    if (project.get("cohort_outcome_schema_version") == "2.0" or framework >= (2, 0, 2)) and not evidence_v2:
        errors.append("framework 2.0.2 requires cohort outcome schema 2.0 with per-subcluster numerical decisions")
    grid = [float(value) for value in outcome["candidate_grid"]]
    resolution_rows = outcome["resolutions"]
    observed = [float(row["resolution"]) for row in resolution_rows]
    if sorted(observed) != sorted(grid) or len(observed) != len(grid):
        errors.append("resolution evidence does not cover the complete candidate grid exactly once")
    if float(outcome["selected_resolution"]) not in grid:
        errors.append("selected resolution is outside the candidate grid")
    selected_assignments: dict[str, set[str]] = {}
    for row in resolution_rows:
        members, member_errors = membership_ids(root, row["membership"], f"resolution {row['resolution']} membership")
        errors.extend(member_errors)
        if members != query:
            errors.append(f"resolution {row['resolution']} membership differs from cohort query")
        membership_path = path_at(root, row["membership"]["path"])
        assignments = read_tsv(membership_path)
        if assignments and "cluster" not in assignments[0]:
            errors.append(f"resolution {row['resolution']} membership lacks cluster assignments")
        elif assignments:
            observed_clusters = len({item.get("cluster", "") for item in assignments if item.get("cluster", "")})
            if observed_clusters != int(row["cluster_count"]):
                errors.append(f"resolution {row['resolution']} cluster_count differs from membership")
            if float(row["resolution"]) == float(outcome["selected_resolution"]):
                for item in assignments:
                    selected_assignments.setdefault(str(item.get("cluster", "")), set()).add(str(item.get("cell_id", "")))
        evidence_path, evidence_errors = validate_evidence_artifact(root, row["evidence_index"], f"resolution {row['resolution']} evidence index")
        errors.extend(evidence_errors)
        if evidence_path and evidence_path.suffix.lower() == ".json":
            try:
                content = json.loads(evidence_path.read_text(encoding="utf-8"))
                if not isinstance(content, (dict, list)) or not content:
                    errors.append(f"resolution {row['resolution']} evidence index has no evidence content")
            except json.JSONDecodeError:
                errors.append(f"resolution {row['resolution']} evidence index is invalid JSON")
    _, full_feature_errors = validate_evidence_artifact(root, outcome["marker_evidence"]["full_feature_scope"]["artifact"], "cohort full-feature evidence")
    errors.extend(full_feature_errors)
    for field in ("spatial_morphology", "source_qc_composition", "marker_evidence"):
        value = outcome.get(field, {})
        if isinstance(value, dict) and set(value) <= {"status", "created_at", "completed_at"}:
            errors.append(f"{field} contains status metadata but no biological evidence")
    selected = float(outcome["selected_resolution"])
    ordered_grid = sorted(grid)
    adjacent = outcome["adjacent_stability"]
    expected_pairs = {(left, right) for left, right in zip(ordered_grid, ordered_grid[1:])}
    observed_pairs = {(float(row.get("left")), float(row.get("right"))) for row in adjacent}
    if observed_pairs != expected_pairs:
        errors.append("adjacent stability must cover every neighboring resolution pair exactly once")
    for row in adjacent:
        if str(row.get("metric", "")).lower() not in {"ari", "adjusted_rand_index", "membership_migration"}:
            errors.append("adjacent stability metric must be ARI or an explicit membership migration metric")
        try:
            value = float(row.get("value"))
            if not -1.0 <= value <= 1.0:
                errors.append("adjacent stability value is outside the valid range")
        except (TypeError, ValueError):
            errors.append("adjacent stability value is not numeric")
    alternatives = {float(row.get("resolution")) for row in outcome["rejected_alternatives"] if "resolution" in row}
    if set(grid) - {selected} != alternatives:
        errors.append("rejected alternatives must explain every non-selected candidate resolution")
    if any(len(str(row.get("reason", "")).strip()) < 10 for row in outcome["rejected_alternatives"]):
        errors.append("every rejected resolution requires a substantive biological rationale")
    families = outcome["marker_evidence"].get("positive_marker_families", [])
    if any(not str(value).strip() for value in families):
        errors.append("positive marker families cannot contain empty entries")
    decision_rows: list[dict[str, str]] = []
    family_rows: list[dict[str, str]] = []
    family_thresholds: dict = {}
    if evidence_v2:
        manifest_ref = outcome.get("broad_family_evidence_manifest")
        validation_ref = outcome.get("broad_family_evidence_validation")
        if not isinstance(manifest_ref, dict) or not isinstance(validation_ref, dict):
            errors.append("cohort outcome schema 2.0 lacks bound broad-family evidence")
        else:
            manifest_path, manifest_errors = validate_evidence_artifact(root, manifest_ref, "cohort broad-family manifest")
            validation_path, validation_errors = validate_evidence_artifact(root, validation_ref, "cohort broad-family validation")
            errors.extend(manifest_errors + validation_errors)
            if manifest_path and validation_path and not manifest_errors and not validation_errors:
                try:
                    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                    validation = json.loads(validation_path.read_text(encoding="utf-8"))
                    if validation.get("status") not in {"PASS", "BLOCKED"} or validation.get("status") != "PASS":
                        errors.append("cohort broad-family evidence validation did not pass")
                    if validation.get("manifest_sha256") != sha256(manifest_path):
                        errors.append("cohort broad-family evidence validation is stale")
                    selected_row = next(
                        (item for item in resolution_rows if float(item["resolution"]) == float(outcome["selected_resolution"])),
                        {},
                    )
                    if manifest.get("cluster_membership", {}).get("sha256") != selected_row.get("membership", {}).get("sha256"):
                        errors.append("cohort broad-family evidence is not bound to the selected-resolution membership")
                    table_ref = manifest.get("evidence_table", {})
                    table_path = path_at(manifest_path.parent, table_ref.get("path", ""))
                    if not table_path.is_file() or sha256(table_path) != table_ref.get("sha256"):
                        errors.append("cohort broad-family evidence table is missing or stale")
                    else:
                        family_rows = read_tsv(table_path)
                    profile_ref = manifest.get("biological_profile", {})
                    profile_path = path_at(manifest_path.parent, profile_ref.get("path", ""))
                    profile = json.loads(profile_path.read_text(encoding="utf-8"))
                    family_thresholds = profile.get("broad_family_evidence_contract", {}).get(
                        "absolute_family_support_thresholds", {}
                    )
                except (OSError, UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError):
                    errors.append("cohort broad-family evidence bindings are unreadable")
        decision_ref = outcome.get("lineage_decision_table")
        if not isinstance(decision_ref, dict):
            errors.append("cohort outcome schema 2.0 lacks lineage_decision_table")
        else:
            decision_path, decision_errors = validate_evidence_artifact(
                root, decision_ref, "cohort lineage decision table"
            )
            errors.extend(decision_errors)
            if decision_path and not decision_errors:
                decision_rows = read_tsv(decision_path)
                errors.extend(
                    validate_decision_rows(
                        decision_rows,
                        require_purity=True,
                        expected_clusters=selected_assignments,
                    )
                )
                if family_rows:
                    for decision_row in decision_rows:
                        try:
                            supported = absolute_supported_families(
                                family_rows,
                                cluster=decision_row.get("cluster", ""),
                                candidate_id=decision_row.get("candidate_id", ""),
                                thresholds=family_thresholds,
                            )
                        except (KeyError, TypeError, ValueError):
                            errors.append("cohort profile lacks valid absolute family support thresholds")
                            break
                        observed = split_families(decision_row.get("positive_families", ""))
                        if supported != observed:
                            errors.append(
                                f"subcluster {decision_row.get('cluster')} candidate {decision_row.get('candidate_id')} "
                                f"positive families were not derived from the bound absolute-expression matrix"
                            )
                else:
                    errors.append("cohort lineage decisions cannot be checked against absolute family evidence")
    outcome_members: set[str] = set()
    for row in outcome["subcluster_outcomes"]:
        members, member_errors = membership_ids(root, row["membership"], f"subcluster {row['subcluster_id']} outcome")
        errors.extend(member_errors)
        overlap = outcome_members & members
        if overlap:
            errors.append(f"subcluster outcomes overlap for {len(overlap)} observations")
        outcome_members |= members
        action = row["outcome"]
        if action in {"parent_return", "fine_return", "cross_lineage_return"} and not str(row.get("target_broad_label", "")).strip():
            errors.append(f"subcluster {row['subcluster_id']} return lacks target broad label")
        if action == "fine_return" and not str(row.get("target_fine_label", "")).strip():
            errors.append(f"subcluster {row['subcluster_id']} fine return lacks target fine label")
        if action == "targeted_successor" and not str(row.get("target_cohort_id", "")).strip():
            errors.append(f"subcluster {row['subcluster_id']} targeted successor lacks target cohort ID")
        if evidence_v2 and action in {"parent_return", "fine_return", "cross_lineage_return"}:
            subcluster = str(row["subcluster_id"])
            source_members = selected_assignments.get(subcluster, set())
            scope = row.get("return_scope")
            for field in ("target_candidate_id", "confidence", "return_scope"):
                if not str(row.get(field, "")).strip():
                    errors.append(f"subcluster {subcluster} v2 biological return lacks {field}")
            if scope not in {"whole_subcluster", "supported_subset"}:
                errors.append(f"subcluster {subcluster} v2 biological return has invalid return_scope")
            if not source_members:
                errors.append(f"subcluster {subcluster} return is absent from the selected-resolution partition")
            elif scope == "whole_subcluster" and members != source_members:
                errors.append(f"subcluster {subcluster} whole-subcluster return does not equal its selected membership")
            elif scope == "supported_subset":
                if not members < source_members:
                    errors.append(f"subcluster {subcluster} supported-subset return is not a strict subset of its selected membership")
                errors.extend(_validate_subset_purity(root, row, members))
            candidate_rows = decision_cluster_rows(decision_rows, subcluster)
            target_candidate = str(row.get("target_candidate_id", ""))
            candidate_by_id = {item.get("candidate_id", ""): item for item in candidate_rows}
            target_row = candidate_by_id.get(target_candidate)
            if not target_row:
                errors.append(f"subcluster {subcluster} target candidate is absent from its decision table")
            else:
                computed = rank_decisions(candidate_rows)
                if computed["winner"] != target_candidate or float(computed["winning_margin"]) <= 0:
                    errors.append(
                        f"subcluster {subcluster} return target is not the positive-margin numerical winner "
                        f"(expected {computed['winner']})"
                    )
                if str(target_row.get("candidate_broad_lineage", "")) != str(row.get("target_broad_label", "")):
                    errors.append(f"subcluster {subcluster} target broad label differs from the decision table")
                if not decision_eligible(target_row):
                    errors.append(f"subcluster {subcluster} return target lacks two families or has unresolved contradictions")
                if scope == "whole_subcluster" and not purity_passes(target_row):
                    errors.append(f"subcluster {subcluster} whole-subcluster return lacks an observation-level purity pass")
                if row.get("confidence") not in {"moderate", "high"}:
                    errors.append(f"subcluster {subcluster} biological return is below moderate confidence")
    if outcome_members != query:
        errors.append("subcluster outcomes do not form the exact query membership partition")
    homogeneous = outcome["homogeneous_parent_confirmed"]
    if homogeneous != (outcome["terminal_outcome"] == "homogeneous_parent_confirmed"):
        errors.append("homogeneous_parent_confirmed boolean and terminal_outcome disagree")
    if homogeneous and any(row["outcome"] != "parent_return" for row in outcome["subcluster_outcomes"]):
        errors.append("homogeneous_parent_confirmed may contain only parent returns")
    if outcome["cohort_type"] == "broad_class_recluster" and outcome["question_mode"] != "broad_purity_audit":
        errors.append("broad-class outcome must use broad_purity_audit")
    if outcome["cohort_type"] != "broad_class_recluster" and outcome["question_mode"] != "targeted_mixture":
        errors.append("targeted outcome must use targeted_mixture")
    if registry_path and registry_path.is_file():
        matches = [row for row in read_tsv(registry_path) if row.get("cohort_id") == outcome["cohort_id"]]
        if len(matches) != 1:
            errors.append("cohort outcome must match exactly one registry row")
        else:
            registry = matches[0]
            if registry.get("membership_sha256") != outcome["query_membership"]["sha256"]:
                errors.append("cohort outcome query hash differs from registry")
            if registry.get("question_mode") != outcome["question_mode"]:
                errors.append("cohort outcome question mode differs from registry")
            if registry.get("outcome_artifact_sha256") and registry["outcome_artifact_sha256"] != sha256(outcome_path):
                errors.append("cohort outcome artifact hash differs from registry")
    return {
        "status": "PASS" if not errors else "FAIL",
        "artifact": str(outcome_path),
        "artifact_sha256": sha256(outcome_path),
        "cohort_id": outcome["cohort_id"],
        "terminal_outcome": outcome["terminal_outcome"],
        "n_query": len(query),
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("outcome", type=Path)
    parser.add_argument("--registry", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = validate(args.project_root.resolve(), args.outcome.resolve(), args.registry)
    write_result(args.out, result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
