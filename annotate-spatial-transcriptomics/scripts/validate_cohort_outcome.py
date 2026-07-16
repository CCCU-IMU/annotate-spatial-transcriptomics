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


def validate(root: Path, outcome_path: Path, registry_path: Path | None = None) -> dict:
    schema = Path(__file__).resolve().parents[1] / "schemas/cohort_outcome.schema.json"
    outcome, errors = validate_json_against_schema(outcome_path, schema)
    if errors:
        return {"status": "FAIL", "artifact": str(outcome_path), "errors": errors}
    query, query_errors = membership_ids(root, outcome["query_membership"], "query membership")
    errors.extend(query_errors)
    grid = [float(value) for value in outcome["candidate_grid"]]
    resolution_rows = outcome["resolutions"]
    observed = [float(row["resolution"]) for row in resolution_rows]
    if sorted(observed) != sorted(grid) or len(observed) != len(grid):
        errors.append("resolution evidence does not cover the complete candidate grid exactly once")
    if float(outcome["selected_resolution"]) not in grid:
        errors.append("selected resolution is outside the candidate grid")
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
