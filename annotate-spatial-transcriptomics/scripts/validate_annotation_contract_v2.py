#!/usr/bin/env python3
"""Fail closed when a v2 annotation contract or any bound artifact is stale."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

from controller_thresholds import load_controller_thresholds
from evidence_schema_lib import sha256, validate_artifact_ref, validate_json_against_schema


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("contract", type=Path)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    schema = Path(__file__).resolve().parents[1] / "schemas/annotation_contract.schema.json"
    contract, errors = validate_json_against_schema(args.contract, schema)
    root = args.contract.parent
    resolved = {}
    for key in (
        "project_config", "workflow_profile", "biological_profile",
        "candidate_catalog", "threshold_registry",
    ):
        path, artifact_errors = validate_artifact_ref(root, contract.get(key, {}), key)
        errors.extend(artifact_errors)
        if path:
            resolved[key] = path
    scope = contract.get("input_scope", {})
    scope_paths: dict[str, Path] = {}
    for key in ("full_object", "analysis_set", "excluded_initial_qc"):
        path, artifact_errors = validate_artifact_ref(root, scope.get(key, {}), f"input scope {key}")
        errors.extend(artifact_errors)
        if path:
            scope_paths[key] = path
    scope_ids: dict[str, set[str]] = {}
    for key in ("analysis_set", "excluded_initial_qc"):
        path = scope_paths.get(key)
        if not path:
            continue
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            values = [str(row.get("cell_id", "")) for row in reader]
            if not reader.fieldnames or "cell_id" not in reader.fieldnames or "" in values or len(values) != len(set(values)):
                errors.append(f"input scope {key} is not a unique cell_id membership")
            scope_ids[key] = set(values)
    if scope_ids.get("analysis_set", set()) & scope_ids.get("excluded_initial_qc", set()):
        errors.append("analysis_set and excluded_initial_qc overlap")
    if len(scope_ids.get("analysis_set", set())) != int(scope.get("analysis_set_n", -1)):
        errors.append("analysis_set count differs from contract")
    if len(scope_ids.get("excluded_initial_qc", set())) != int(scope.get("excluded_initial_qc_n", -1)):
        errors.append("excluded_initial_qc count differs from contract")
    if scope.get("biological_unresolved_is_initial_qc") is not False:
        errors.append("biological unresolved members cannot be frozen as initial QC")
    controller = contract.get("canonical_lineage_controller", {})
    script_dir = Path(__file__).resolve().parent
    selector = script_dir / "select_lineage_resolution.py"
    selector_record = controller.get("resolution_selector", {})
    if (
        Path(str(selector_record.get("path", ""))).resolve() != selector.resolve()
        or not selector.is_file()
        or selector_record.get("sha256") != sha256(selector)
    ):
        errors.append("annotation contract does not bind the canonical resolution selector")
    builder = script_dir / "build_resolution_grid_evidence.py"
    builder_record = controller.get("resolution_evidence_builder", {})
    if (
        Path(str(builder_record.get("path", ""))).resolve() != builder.resolve()
        or not builder.is_file()
        or builder_record.get("sha256") != sha256(builder)
    ):
        errors.append("annotation contract does not bind canonical resolution evidence")
    discovery = script_dir / "discover_unmodeled_lineages.py"
    discovery_record = controller.get("unmodeled_discovery", {})
    if (
        Path(str(discovery_record.get("path", ""))).resolve() != discovery.resolve()
        or not discovery.is_file()
        or discovery_record.get("sha256") != sha256(discovery)
    ):
        errors.append("annotation contract does not bind canonical unmodeled-lineage discovery")
    cohort_contract = controller.get("cohort_reclustering", {})
    for key, name in (
        ("entrypoint", "run_seurat_cohort_recluster.R"),
        ("implementation", "run_seurat_cohort_recluster_impl.R"),
    ):
        expected = script_dir / name
        record = cohort_contract.get(key, {})
        if (
            Path(str(record.get("path", ""))).resolve() != expected.resolve()
            or not expected.is_file()
            or record.get("sha256") != sha256(expected)
        ):
            errors.append(f"annotation contract does not bind canonical cohort {key}")
    if (
        cohort_contract.get("sct_input_boundary")
        != "project_local_non_sct_raw_counts"
        or cohort_contract.get("clustering_path") != "SCT_PCA_SNN_Leiden"
    ):
        errors.append("annotation contract does not bind the cohort expression boundary")
    subset_policy = controller.get("subset_policy", {})
    if not (
        subset_policy.get("identity_core_component_required") is True
        and subset_policy.get("generic_support_cell_transitive_bridging") is False
        and subset_policy.get("fine_candidate_parent_lock") is True
        and subset_policy.get("fine_candidate_discriminator_seed_required") is True
    ):
        errors.append("annotation contract lacks the v2.2 identity-core subset policy")
    if subset_policy.get("activation_scope") != "second_round_mixed_subcluster_only":
        errors.append("candidate-local subset logic is not limited to second-round mixed subclusters")
    remainder_policy = controller.get("remainder_policy", {})
    if not (
        remainder_policy.get("allowed_scopes")
        == ["local_mixed_subcluster", "post_merge_reconciliation"]
        and remainder_policy.get("first_round_forbidden") is True
        and remainder_policy.get("unselected_member_is_qc") is False
    ):
        errors.append("annotation contract permits invalid first-round remainder closure")
    expected_phases = [
        "whole_tissue_partition", "cluster_cohort_recluster",
        "local_mixed_subcluster_split", "merge_and_freeze_broad",
        "atlas_and_completeness_review", "materialize_final_release",
    ]
    if controller.get("phase_order") != expected_phases:
        errors.append("annotation contract phase order differs from v2.2 architecture")
    if controller.get("phase_authority") != {
        "whole_tissue_partition": "provisional_only_no_release_membership",
        "cluster_cohort_recluster": "candidate_only_no_release_membership",
        "local_mixed_subcluster_split": "candidate_only_no_release_membership",
        "merge_and_freeze_broad": "formal_broad_freeze",
        "atlas_and_completeness_review": "unlabeled_broad_rescue_only",
        "materialize_final_release": "final_broad_fine_state_release",
    }:
        errors.append("annotation contract phase authority differs from v2.2 architecture")
    for name in (
        "run_observation_lineage_scoring.R",
        "derive_candidate_local_subsets.R",
        "close_exact_remainders.py",
        "build_whole_tissue_cohort_plan.py",
        "adjudicate_second_round_subclusters.py",
        "build_candidate_context_evidence.py",
        "merge_and_freeze_broad_membership.py",
        "route_global_atlas_v2.py",
        "validate_global_atlas_v2.py",
        "apply_post_merge_atlas_routing.py",
        "review_post_merge_unresolved_components.py",
        "audit_post_merge_completeness.py",
        "validate_sheep_ovary_biological_quality.py",
        "apply_sheep_ovary_follicle_roi_repair.py",
        "screen_rare_cell_programs.R",
        "screen_spatial_foci.py",
        "materialize_oocyte_cluster_membership.py",
        "apply_cell_id_membership_patch.py",
        "materialize_parent_locked_fine_proposals.py",
        "materialize_final_release_v2_2.py",
        "evaluate_annotation_robustness.py",
        "run_lineage_controller.py",
    ):
        expected = script_dir / name
        record = controller.get("scripts", {}).get(name, {})
        try:
            observed = Path(str(record.get("path", ""))).resolve()
        except (OSError, RuntimeError):
            observed = Path()
        if observed != expected.resolve():
            errors.append(f"canonical lineage controller path differs for {name}")
        elif not expected.is_file() or record.get("sha256") != sha256(expected):
            errors.append(f"canonical lineage controller hash differs for {name}")
    for dependency_name in (
        "controller_thresholds.py", "lineage_controller_lib.py"
    ):
        dependency = script_dir / dependency_name
        dependency_record = controller.get("dependencies", {}).get(
            dependency_name, {}
        )
        if (
            Path(str(dependency_record.get("path", ""))).resolve()
            != dependency.resolve()
            or not dependency.is_file()
            or dependency_record.get("sha256") != sha256(dependency)
        ):
            errors.append(
                f"annotation contract does not bind {dependency_name}"
            )
    if "threshold_registry" in resolved:
        try:
            thresholds = load_controller_thresholds(resolved["threshold_registry"])
            expected_writeback = thresholds["observation_writeback_policy"]
            expected_scoring = thresholds["scoring_policy"]
            expected_local = thresholds["local_subset_policy"]
            contract_writeback = contract.get(
                "observation_writeback", {}
            ).get("policy", {})
            scoring = controller.get("scoring_policy", {})
            for key in (
                "direct_weight", "local_weight", "anti_weight",
                "family_active_threshold", "local_gene_detection_fraction",
            ):
                if scoring.get(key) != expected_scoring.get(key):
                    errors.append(
                        f"controller scoring policy differs from threshold registry: {key}"
                    )
            subset = controller.get("subset_policy", {})
            expected_subset = {
                "subset_validation_supported_fraction": contract_writeback[
                    "supported_subset_min_lineage_supported_fraction"
                ],
                "subset_validation_competitor_margin": contract_writeback[
                    "supported_subset_min_purity_margin"
                ],
                "maximum_contradiction_fraction": contract_writeback[
                    "maximum_contradiction_fraction"
                ],
                "maximum_second_subset_rounds": expected_local[
                    "maximum_second_subset_rounds"
                ],
            }
            for key, value in expected_subset.items():
                if subset.get(key) != value:
                    errors.append(
                        f"controller subset policy differs from threshold registry: {key}"
                    )
            if contract_writeback != expected_writeback:
                project_policy = {}
                if "project_config" in resolved:
                    project_policy = json.loads(
                        resolved["project_config"].read_text(encoding="utf-8")
                    ).get("observation_writeback_policy", {})
                if contract_writeback != project_policy:
                    errors.append(
                        "observation-writeback policy is neither the canonical registry nor an explicit project override"
                    )
        except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
            errors.append(f"invalid controller threshold registry: {exc}")
    if "project_config" in resolved:
        project = json.loads(resolved["project_config"].read_text(encoding="utf-8"))
        for key in ("project_id", "sample_id", "framework_version"):
            if project.get(key) != contract.get(key):
                errors.append(f"project config disagrees with contract: {key}")
        if project.get("observation_subset_writeback_required") is True:
            bound = contract.get("observation_writeback", {})
            if bound.get("subset_writeback_required") is not True:
                errors.append("annotation contract does not bind required observation-level subset writeback")
            if bound.get("policy") != project.get("observation_writeback_policy"):
                errors.append("annotation contract observation-writeback policy differs from project config")
    snapshot_registry = Path(contract.get("selected_input_snapshot", {}).get("registry_path", ""))
    if snapshot_registry.is_file():
        with snapshot_registry.open(newline="", encoding="utf-8") as handle:
            rows = [row for row in csv.DictReader(handle, delimiter="\t") if row.get("snapshot_id") == contract.get("selected_input_snapshot", {}).get("snapshot_id")]
        if len(rows) != 1 or rows[0].get("sha256") != contract.get("selected_input_snapshot", {}).get("sha256"):
            errors.append("selected input snapshot is missing, duplicated or changed")
        elif Path(rows[0].get("path", "")).resolve() != Path(
            contract.get("selected_input_snapshot", {}).get("path", "")
        ).resolve():
            errors.append("selected input snapshot path differs from registry")
    else:
        errors.append("selected input snapshot registry is missing")
    whole = contract.get("whole_tissue_partition", {})
    query = contract.get("query_reclustering", {})
    if (
        whole.get("selection_endpoint") != "stable_cohort_partition"
        or whole.get("formal_membership_authority") is not False
    ):
        errors.append("whole-tissue phase is not provisional-only")
    partition_grid_path, partition_grid_errors = validate_artifact_ref(
        root, whole.get("partition_grid", {}), "whole-tissue partition grid"
    )
    errors.extend(partition_grid_errors)
    if not (
        query.get("cohort_unit") == "initial_cluster"
        and query.get("normalization_path") == "raw_counts_SCTv2_PCA_SNN_Leiden"
        and query.get("full_catalog_scan_required") is True
        and query.get("provisional_broad_blinded_until_score_freeze") is True
        and query.get("resolution_stability_neighbors")
        == "nearest_lower_and_higher_else_two_nearest"
    ):
        errors.append("query reclustering is not initial-cluster, raw-count and label-blind")
    writeback = contract.get("observation_writeback", {})
    if (
        writeback.get("first_round_formal_writeback") is not False
        or writeback.get("local_mixed_subcluster_only") is not True
    ):
        errors.append("observation writeback is not restricted to local second-round mixtures")
    role_policy = contract.get("artifact_role_policy", {})
    selected_record = contract.get("selected_input_snapshot", {})
    if (
        role_policy.get("runtime_allowed_roles")
        != ["runtime_input", "external_reference"]
        or role_policy.get("failed_diagnostic_runtime_forbidden") is not True
        or selected_record.get("artifact_role") != "runtime_input"
    ):
        errors.append("artifact role policy does not exclude failed diagnostics")
    diagnostic_record = role_policy.get("diagnostic_registry")
    if diagnostic_record:
        diagnostic_path, diagnostic_errors = validate_artifact_ref(
            root, diagnostic_record, "failed diagnostic registry"
        )
        errors.extend(diagnostic_errors)
        if diagnostic_path:
            with diagnostic_path.open(newline="", encoding="utf-8") as handle:
                failed = [
                    Path(row.get("artifact_root", "")).resolve()
                    for row in csv.DictReader(handle, delimiter="\t")
                    if row.get("artifact_role") == "failed_diagnostic"
                    and row.get("artifact_root")
                ]
            selected_path = Path(selected_record.get("path", "")).resolve()
            if any(selected_path == path or path in selected_path.parents for path in failed):
                errors.append("selected runtime input is inside failed_diagnostic artifacts")
    if whole.get("method") == "BANKSY" and whole.get("candidate_grid_source") not in {"bound_upstream_input", "fresh_project_computation"}:
        errors.append("BANKSY whole-tissue selection requires a declared fresh or bound-upstream grid source")
    if whole.get("candidate_grid_source") == "bound_upstream_input" or whole.get("method") == "BANKSY":
        grid_path, grid_errors = validate_artifact_ref(root, whole.get("grid_artifact", {}), "whole-tissue grid artifact")
        errors.extend(grid_errors)
        if grid_path:
            try:
                grid_payload = json.loads(grid_path.read_text(encoding="utf-8"))
                grid = grid_payload.get("candidate_resolutions", grid_payload.get("resolutions", []))
                if sorted({float(value) for value in grid}) != whole.get("candidate_resolutions"):
                    errors.append("whole-tissue resolutions differ from the bound upstream grid artifact")
            except (OSError, UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError):
                errors.append("whole-tissue grid artifact is not valid resolution JSON")
    if whole.get("candidate_grid_source") == "fresh_project_computation" and "biological_profile" in resolved:
        biological = json.loads(resolved["biological_profile"].read_text(encoding="utf-8"))
        fresh = biological.get("resolution_policy", {}).get("fresh_sct_banksy_whole_tissue_candidate_resolutions", [])
        if fresh and sorted({float(value) for value in fresh}) != whole.get("candidate_resolutions"):
            errors.append("fresh whole-tissue resolutions differ from the bound biological profile")
    if whole.get("method") == "BANKSY" and whole.get("candidate_resolutions") == query.get("candidate_resolutions"):
        # Equality is possible, but it must not arise by silently substituting
        # the query grid. The bound-upstream source above is the decisive gate.
        pass
    result = {
        "status": "PASS" if not errors else "BLOCKED",
        "schema_version": "2.0",
        "contract": str(args.contract.resolve()),
        "contract_sha256": sha256(args.contract),
        "errors": errors,
    }
    out = args.out or args.contract.parent.parent / "provenance/annotation_contract_validation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
