#!/usr/bin/env python3
"""Freeze the v2 project/profile/input/resolution contract before annotation."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
from pathlib import Path

from evidence_schema_lib import sha256
from controller_thresholds import REGISTRY_PATH, load_controller_thresholds
from lineage_decision_lib import observation_writeback_policy


def artifact(path: Path, role: str = "runtime_input") -> dict[str, str]:
    if not path.is_file():
        raise SystemExit(f"missing contract artifact: {path}")
    return {
        "path": str(path.resolve()), "sha256": sha256(path),
        "artifact_role": role,
    }


def freeze_artifact(
    source: Path, root: Path, role: str, contract_id: str = ""
) -> dict[str, str]:
    """Copy mutable installed profiles into the immutable project contract."""
    if not source.is_file():
        raise SystemExit(f"missing contract artifact: {source}")
    suffix = "".join(source.suffixes) or ".txt"
    profile_root = root / "config/contract_profiles"
    if contract_id:
        profile_root = profile_root / contract_id
    destination = profile_root / f"{role}{suffix}"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    if sha256(destination) != sha256(source):
        raise SystemExit(f"failed to freeze {role} into the project")
    return artifact(destination)


def resolutions(value: str) -> list[float]:
    try:
        result = sorted({float(x) for x in value.replace(";", ",").split(",") if x.strip()})
    except ValueError as exc:
        raise SystemExit(f"invalid resolution grid: {exc}")
    if len(result) < 3:
        raise SystemExit("v2 whole-tissue contract requires at least three candidate resolutions")
    return result


def grid_artifact(path: Path, expected: list[float]) -> dict[str, str]:
    record = artifact(path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        observed = payload.get("candidate_resolutions", payload.get("resolutions", []))
        observed = sorted({float(value) for value in observed})
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError):
        raise SystemExit("--whole-tissue-grid-artifact must be JSON with candidate_resolutions/resolutions")
    if observed != expected:
        raise SystemExit("whole-tissue grid differs from the bound upstream grid artifact")
    return record


def membership_ids(path: Path, allow_empty: bool = False) -> set[str]:
    import gzip

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "cell_id" not in reader.fieldnames:
            raise SystemExit(f"membership lacks cell_id: {path}")
        values = [str(row.get("cell_id", "")) for row in reader]
    if "" in values or len(values) != len(set(values)) or (not values and not allow_empty):
        raise SystemExit(f"membership has empty/duplicate observations: {path}")
    return set(values)


def validate_partition_grid(
    path: Path, expected_resolutions: list[float], analysis_ids: set[str]
) -> None:
    import gzip

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        required = {"cell_id", "resolution", "cluster"}
        if not required <= fields:
            raise SystemExit(
                "whole-tissue partition grid lacks: "
                + ", ".join(sorted(required - fields))
            )
        forbidden = {
            "provisional_broad", "initial_broad_label", "broad_label",
            "fine_label", "final_broad_label", "final_fine_label",
            "celltype", "cell_type", "annotation", "historical_label",
            "repair_membership", "atlas_label",
        }
        leaked = sorted(forbidden & {field.lower() for field in fields})
        if leaked:
            raise SystemExit(
                "whole-tissue partition grid exposes annotation columns: "
                + ", ".join(leaked)
            )
        coverage: dict[float, set[str]] = {}
        for row in reader:
            try:
                resolution = float(str(row.get("resolution", "")))
            except ValueError as exc:
                raise SystemExit("whole-tissue partition grid has invalid resolution") from exc
            cell = str(row.get("cell_id", ""))
            cluster = str(row.get("cluster", ""))
            if not cell or not cluster:
                raise SystemExit("whole-tissue partition grid contains an incomplete row")
            if cell in coverage.setdefault(resolution, set()):
                raise SystemExit(
                    "whole-tissue partition grid duplicates an observation within a resolution"
                )
            coverage[resolution].add(cell)
    if sorted(coverage) != expected_resolutions:
        raise SystemExit("whole-tissue partition grid resolutions differ from the contract")
    if any(ids != analysis_ids for ids in coverage.values()):
        raise SystemExit(
            "every whole-tissue partition resolution must exactly cover analysis_set"
        )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("project_root", type=Path)
    ap.add_argument("--workflow-profile", required=True, type=Path)
    ap.add_argument("--biological-profile", required=True, type=Path)
    ap.add_argument("--candidate-catalog", required=True, type=Path)
    ap.add_argument(
        "--context-evidence",
        type=Path,
        help=(
            "optional exogenous biological-context evidence; this only "
            "permits evaluation of context-gated candidates and cannot "
            "assign observations"
        ),
    )
    ap.add_argument("--snapshot-id", required=True)
    ap.add_argument("--analysis-membership", required=True, type=Path)
    ap.add_argument("--excluded-initial-qc", required=True, type=Path)
    ap.add_argument("--whole-tissue-method", choices=["BANKSY", "Seurat", "Scanpy", "external"], required=True)
    ap.add_argument("--whole-tissue-grid", required=True)
    ap.add_argument("--whole-tissue-grid-artifact", type=Path)
    ap.add_argument("--whole-tissue-partitions", required=True, type=Path)
    ap.add_argument("--grid-source", choices=["bound_upstream_input", "fresh_project_computation"], required=True)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--seed", type=int, default=2200)
    ap.add_argument(
        "--contract-id",
        default="",
        help="optional immutable profile namespace for a new blind attempt",
    )
    args = ap.parse_args()
    if args.contract_id and not re.fullmatch(r"[A-Za-z0-9._-]+", args.contract_id):
        raise SystemExit("--contract-id must be one safe path component")
    root = args.project_root.resolve()
    project_path = root / "config/project.json"
    project = json.loads(project_path.read_text(encoding="utf-8"))
    if project.get("framework_version") != "2.0.0":
        raise SystemExit("project framework_version must be 2.0.0")
    snapshot_registry = root / "state/input_snapshot_registry.tsv"
    with snapshot_registry.open(newline="", encoding="utf-8") as handle:
        matches = [row for row in csv.DictReader(handle, delimiter="\t") if row.get("snapshot_id") == args.snapshot_id and row.get("status") in {"frozen", "validated", "active"}]
    if len(matches) != 1 or len(matches[0].get("sha256", "")) != 64:
        raise SystemExit("selected snapshot is not one unique frozen/validated registry row")
    diagnostic_registry = root / "provenance/failed_diagnostic_artifact_registry.tsv"
    failed_roots: list[Path] = []
    if diagnostic_registry.is_file():
        with diagnostic_registry.open(newline="", encoding="utf-8") as handle:
            failed_rows = list(csv.DictReader(handle, delimiter="\t"))
        failed_roots = [
            Path(row["artifact_root"]).resolve()
            for row in failed_rows
            if row.get("artifact_role") == "failed_diagnostic"
            and row.get("artifact_root")
        ]
    selected_path = Path(matches[0].get("path", "")).resolve()
    if not selected_path.is_file() or sha256(selected_path) != matches[0]["sha256"]:
        raise SystemExit("selected runtime snapshot is missing or stale")
    if any(selected_path == failed or failed in selected_path.parents for failed in failed_roots):
        raise SystemExit("selected input snapshot is inside failed_diagnostic artifacts")
    biological = json.loads(args.biological_profile.read_text(encoding="utf-8"))
    controller_thresholds = load_controller_thresholds()
    scoring_policy = controller_thresholds["scoring_policy"]
    local_subset_policy = controller_thresholds["local_subset_policy"]
    writeback_policy = observation_writeback_policy(project)
    analysis_ids = membership_ids(args.analysis_membership)
    excluded_ids = membership_ids(args.excluded_initial_qc, allow_empty=True)
    if analysis_ids & excluded_ids:
        raise SystemExit("analysis_set and excluded_initial_qc overlap")
    whole_grid = resolutions(args.whole_tissue_grid)
    validate_partition_grid(args.whole_tissue_partitions, whole_grid, analysis_ids)
    query_grid = biological.get("resolution_policy", {}).get("query_reclustering_candidate_resolutions", [])
    if len(query_grid) < 3:
        raise SystemExit("biological profile lacks the v2 query-reclustering grid")
    if args.grid_source == "bound_upstream_input" and not args.whole_tissue_grid_artifact:
        raise SystemExit("bound_upstream_input requires --whole-tissue-grid-artifact")
    if args.whole_tissue_method == "BANKSY" and not args.whole_tissue_grid_artifact:
        raise SystemExit("BANKSY whole-tissue computation requires --whole-tissue-grid-artifact")
    fresh_grid = biological.get("resolution_policy", {}).get("fresh_sct_banksy_whole_tissue_candidate_resolutions", [])
    if args.grid_source == "fresh_project_computation" and fresh_grid and whole_grid != sorted({float(value) for value in fresh_grid}):
        raise SystemExit("fresh whole-tissue grid differs from the bound biological profile")
    contract = {
        "schema_version": "2.0",
        "framework_version": "2.0.0",
        "skill_release_version": "2.2.0",
        "project_id": project["project_id"],
        "sample_id": project["sample_id"],
        "project_config": artifact(project_path),
        "workflow_profile": freeze_artifact(
            args.workflow_profile, root, "workflow_profile", args.contract_id
        ),
        "biological_profile": freeze_artifact(
            args.biological_profile, root, "biological_profile", args.contract_id
        ),
        "candidate_catalog": freeze_artifact(
            args.candidate_catalog, root, "candidate_catalog", args.contract_id
        ),
        "threshold_registry": freeze_artifact(
            REGISTRY_PATH, root, "controller_threshold_registry", args.contract_id
        ),
        "candidate_context_evidence": (
            freeze_artifact(
                args.context_evidence,
                root,
                "candidate_context_evidence",
                args.contract_id,
            )
            if args.context_evidence else None
        ),
        "selected_input_snapshot": {
            "registry_path": str(snapshot_registry.resolve()),
            "snapshot_id": args.snapshot_id,
            "sha256": matches[0]["sha256"],
            "path": str(selected_path),
            "artifact_role": "runtime_input",
        },
        "input_scope": {
            "full_object": {
                "path": str(selected_path),
                "sha256": matches[0]["sha256"],
                "artifact_role": "runtime_input",
            },
            "analysis_set": artifact(args.analysis_membership),
            "excluded_initial_qc": artifact(args.excluded_initial_qc),
            "analysis_set_n": len(analysis_ids),
            "excluded_initial_qc_n": len(excluded_ids),
            "biological_unresolved_is_initial_qc": False,
        },
        "artifact_role_policy": {
            "runtime_allowed_roles": ["runtime_input", "external_reference"],
            "failed_diagnostic_runtime_forbidden": True,
            "diagnostic_registry": (
                artifact(diagnostic_registry, "failed_diagnostic_registry")
                if diagnostic_registry.is_file() else None
            ),
        },
        "expression_ancestry_policy": {
            "project_local_query_evidence": True,
            "cross_project_expression_requires_external_reference": True,
        },
        "whole_tissue_partition": {
            "method": args.whole_tissue_method,
            "candidate_grid_source": args.grid_source,
            "candidate_resolutions": whole_grid,
            "selection_endpoint": "stable_cohort_partition",
            "formal_membership_authority": False,
            "partition_grid": artifact(args.whole_tissue_partitions),
        },
        "query_reclustering": {
            "candidate_resolutions": query_grid,
            "separate_from_upstream_banksy_grid": True,
            "cohort_unit": "initial_cluster",
            "normalization_path": "raw_counts_SCTv2_PCA_SNN_Leiden",
            "full_catalog_scan_required": True,
            "provisional_broad_blinded_until_score_freeze": True,
            "resolution_stability_neighbors": "nearest_lower_and_higher_else_two_nearest",
        },
        "broad_family_evidence": {
            "required": True,
            "feature_scope": "full_feature",
            "complete_cartesian_product": True,
        },
        "canonical_lineage_controller": {
            "controller_version": "2.2.0",
            "formal_release_requires_canonical_chain": True,
            "historical_labels_blinded_until_membership_freeze": True,
            "phase_order": [
                "whole_tissue_partition",
                "cluster_cohort_recluster",
                "local_mixed_subcluster_split",
                "merge_and_freeze_broad",
                "atlas_and_completeness_review",
                "materialize_final_release",
            ],
            "phase_authority": {
                "whole_tissue_partition": "provisional_only_no_release_membership",
                "cluster_cohort_recluster": "candidate_only_no_release_membership",
                "local_mixed_subcluster_split": "candidate_only_no_release_membership",
                "merge_and_freeze_broad": "formal_broad_freeze",
                "atlas_and_completeness_review": "unlabeled_broad_rescue_only",
                "materialize_final_release": "final_broad_fine_state_release",
            },
            "random_seed": args.seed,
            "resolution_selector": artifact(
                Path(__file__).resolve().parent / "select_lineage_resolution.py"
            ),
            "resolution_evidence_builder": artifact(
                Path(__file__).resolve().parent / "build_resolution_grid_evidence.py"
            ),
            "unmodeled_discovery": artifact(
                Path(__file__).resolve().parent / "discover_unmodeled_lineages.py"
            ),
            "cohort_reclustering": {
                "entrypoint": artifact(
                    Path(__file__).resolve().parent
                    / "run_seurat_cohort_recluster.R"
                ),
                "implementation": artifact(
                    Path(__file__).resolve().parent
                    / "run_seurat_cohort_recluster_impl.R"
                ),
                "sct_input_boundary": "project_local_non_sct_raw_counts",
                "clustering_path": "SCT_PCA_SNN_Leiden",
            },
            "scripts": {
                name: artifact(Path(__file__).resolve().parent / name)
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
                )
            },
            "dependencies": {
                name: artifact(Path(__file__).resolve().parent / name)
                for name in (
                    "controller_thresholds.py", "lineage_controller_lib.py"
                )
            },
            "scoring_policy": {
                "gene_scale": "query_nonzero_q95_capped",
                "family_aggregation": "mean_top_two_available_genes",
                "direct_weight": scoring_policy["direct_weight"],
                "local_weight": scoring_policy["local_weight"],
                "anti_weight": scoring_policy["anti_weight"],
                "family_active_threshold": scoring_policy[
                    "family_active_threshold"
                ],
                "local_gene_detection_fraction": scoring_policy[
                    "local_gene_detection_fraction"
                ],
                "hard_contradiction": (
                    "coherent_multigene_multifamily_direct_anti_only"
                ),
            },
            "subset_policy": {
                "candidate_local_independent_proposals": True,
                "aggregate_winner_can_veto": False,
                "identity_core_component_required": True,
                "generic_support_cell_transitive_bridging": False,
                "fine_candidate_parent_lock": True,
                "fine_candidate_discriminator_seed_required": True,
                "sparse_tail_inheritance_scope": (
                    "validated_high_purity_expression_subcluster_only"
                ),
                "subset_validation_supported_fraction": writeback_policy[
                    "supported_subset_min_lineage_supported_fraction"
                ],
                "subset_validation_competitor_margin": writeback_policy[
                    "supported_subset_min_purity_margin"
                ],
                "maximum_contradiction_fraction": writeback_policy[
                    "maximum_contradiction_fraction"
                ],
                "maximum_second_subset_rounds": local_subset_policy[
                    "maximum_second_subset_rounds"
                ],
                "activation_scope": "second_round_mixed_subcluster_only",
            },
            "remainder_policy": {
                "allowed_scopes": [
                    "local_mixed_subcluster", "post_merge_reconciliation"
                ],
                "first_round_forbidden": True,
                "unselected_member_is_qc": False,
            },
        },
        "observation_writeback": {
            "first_round_formal_writeback": False,
            "local_mixed_subcluster_only": True,
            "subset_writeback_required": project.get("observation_subset_writeback_required") is True,
            "whole_subcluster_purity_evidence_required": project.get("whole_subcluster_purity_evidence_required") is True,
            "terminal_return_purity_audit_required": project.get("terminal_return_purity_audit_required") is True,
            "raw_two_family_writeback_audit_required": project.get("raw_two_family_writeback_audit_required") is True,
            "complete_fine_candidate_audit_required": project.get("complete_fine_candidate_audit_required") is True,
            "fine_writeback_broad_lock_required": project.get("fine_writeback_broad_lock_required") is True,
            "final_fine_state_validation_required": project.get("final_fine_state_validation_required") is True,
            "policy": writeback_policy,
        },
        "atlas_routing": {
            "authoritative_router": "route_global_atlas_v2.py",
            "mapping_scope": "complete_analysis_set",
            "writeback_scope": "unlabeled_after_broad_merge_only",
            "fine_anchor_eligible": False,
        },
        "release_taxonomy": {
            "vascular_parent": "Vascular-associated",
            "hierarchy_validation_required": True,
        },
    }
    if args.whole_tissue_grid_artifact:
        contract["whole_tissue_partition"]["grid_artifact"] = grid_artifact(args.whole_tissue_grid_artifact, whole_grid)
    out = args.out or root / "config/annotation_contract.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(contract, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(contract, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
