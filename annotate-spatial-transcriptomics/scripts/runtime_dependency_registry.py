#!/usr/bin/env python3
"""Single source of truth for release-authorized runtime dependencies.

The annotation contract builder, validator and controller must all import this
module.  Keeping independent script-name tuples in those three locations was a
source of false contract failures and, more importantly, allowed a project to
run a helper that was never frozen into its contract.
"""

from __future__ import annotations


REGISTRY_VERSION = "2.5.3"

PHASE_ORDER = (
    "whole_tissue_partition",
    "cluster_cohort_recluster",
    "local_mixed_subcluster_split",
    "merge_and_freeze_broad",
    "atlas_and_completeness_review",
    "materialize_final_release",
)

CANONICAL_SCRIPTS = (
    "bootstrap_sct_banksy_project.py",
    "freeze_sct_banksy_input.R",
    "validate_phase_runtime.py",
    "plan_cohort_resources.py",
    "run_observation_lineage_scoring.R",
    "derive_candidate_local_subsets.R",
    "close_exact_remainders.py",
    "build_whole_tissue_cohort_plan.py",
    "adjudicate_second_round_subclusters.py",
    "audit_local_split_workload.py",
    "build_candidate_context_evidence.py",
    "merge_and_freeze_broad_membership.py",
    "export_seurat_fixed_atlas_features.R",
    "map_fixed_atlas_projection.py",
    "build_program_anchor_membership.R",
    "calibrate_tiered_mapping_thresholds.py",
    "bind_atlas_routing_mapping.py",
    "route_global_atlas_v2.py",
    "validate_global_atlas_v2.py",
    "apply_post_merge_atlas_routing.py",
    "review_post_merge_unresolved_components.py",
    "audit_post_merge_completeness.py",
    "audit_catalog_wide_lineage_challengers.py",
    "build_cell_type_review_marker_manifest.py",
    "export_cell_type_review_counts.R",
    "build_zero_census_direct_challengers.py",
    "build_broad_cell_type_review_evidence.py",
    "export_broad_cell_type_review_pseudobulk.R",
    "summarize_broad_cell_type_review_pseudobulk.py",
    "build_broad_cell_type_review_packet_index.py",
    "manage_cell_type_review_queue.py",
    "validate_catalog_wide_lineage_review_decisions.py",
    "apply_catalog_wide_lineage_review.py",
    "manage_membership_transform_chain.py",
    "validate_fixed_atlas_bundle.py",
    "validate_sheep_ovary_biological_quality.py",
    "apply_sheep_ovary_follicle_roi_repair.py",
    "screen_rare_cell_programs.R",
    "screen_spatial_foci.py",
    "materialize_oocyte_cluster_membership.py",
    "apply_cell_id_membership_patch.py",
    "materialize_parent_locked_fine_proposals.py",
    "materialize_final_release_v2_2.py",
    "materialize_final_deliverables.py",
    "build_release_evidence_tables.py",
    "write_frozen_annotations_to_seurat.R",
    "build_annotation_maps.R",
    "build_marker_dotplots.R",
    "build_spatial_gene_maps.R",
    "run_final_label_deg.R",
    "build_frozen_review_report.py",
    "evaluate_annotation_robustness.py",
    "run_lineage_controller.py",
)

CANONICAL_DEPENDENCIES = (
    "controller_runtime_state.py",
    "controller_thresholds.py",
    "lineage_controller_lib.py",
    "membership_transform_lib.py",
    "runtime_dependency_registry.py",
)


def registry_document() -> dict[str, object]:
    return {
        "registry_version": REGISTRY_VERSION,
        "phase_order": list(PHASE_ORDER),
        "scripts": list(CANONICAL_SCRIPTS),
        "dependencies": list(CANONICAL_DEPENDENCIES),
    }
