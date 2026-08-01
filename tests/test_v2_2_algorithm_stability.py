from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "annotate-spatial-transcriptomics/scripts"
sys.path.insert(0, str(SCRIPTS))

from lineage_controller_lib import (  # noqa: E402
    aggregate_program_supported,
    canonical_cluster_challenger,
    candidate_can_release,
    candidate_core_seed,
    catalog_candidates,
    choose_group_parent,
    deterministic_membership_hash,
    dominant_generic_remainder_group,
    effective_broad_writeback_strategy,
    fine_audit_complete,
    group_candidate_detected,
    hard_contradiction,
    independent_group_program,
    local_split_worthy_group_program,
    pairwise_separable_identity_components,
    rare_group_program_watch,
    index_scores,
    rank_supported,
    specific_component_embedded_in_generic_parent,
    supported_seed,
    validate_subset,
    validate_canonical_identity_component,
    validate_unmodeled_programs,
)
from build_resolution_grid_evidence import candidate_detected  # noqa: E402
from close_exact_remainders import (  # noqa: E402
    provenance_subset,
    remainder_candidate_programs,
    validated_residual_component_candidates,
)


CANDIDATES = {
    "granulosa": {
        "candidate_id": "granulosa",
        "candidate_role": "broad",
        "release_broad_label": "Granulosa",
        "specificity_priority": 80,
    },
    "stromal_mesenchymal": {
        "candidate_id": "stromal_mesenchymal",
        "candidate_role": "broad",
        "release_broad_label": "Stromal/mesenchymal",
        "specificity_priority": 10,
    },
    "smooth_muscle": {
        "candidate_id": "smooth_muscle",
        "candidate_role": "broad",
        "release_broad_label": "Smooth muscle",
        "specificity_priority": 80,
    },
    "pericyte_mural": {
        "candidate_id": "pericyte_mural",
        "candidate_role": "broad",
        "release_broad_label": "Pericyte/mural",
        "release_fine_label": "",
        "specificity_priority": 90,
    },
    "vascular_endothelial": {
        "candidate_id": "vascular_endothelial",
        "candidate_role": "broad",
        "release_broad_label": "Endothelial",
        "specificity_priority": 90,
    },
    "lymphatic_endothelial": {
        "candidate_id": "lymphatic_endothelial",
        "candidate_role": "fine",
        "release_broad_label": "Endothelial",
        "release_fine_label": "Lymphatic endothelial",
        "parent_broad_label": "Endothelial",
        "specificity_priority": 95,
    },
    "epithelial_mesothelial": {
        "candidate_id": "epithelial_mesothelial",
        "candidate_role": "broad",
        "release_broad_label": "Epithelial/mesothelial",
        "specificity_priority": 90,
    },
    "theca_steroidogenic": {
        "candidate_id": "theca_steroidogenic",
        "candidate_role": "broad",
        "release_broad_label": "Theca",
        "specificity_priority": 90,
    },
    "oocyte": {
        "candidate_id": "oocyte",
        "candidate_role": "broad",
        "release_broad_label": "Oocyte",
        "specificity_priority": 100,
        "writeback_strategy": "canonical_cluster_membership",
    },
    "exploratory_neural": {
        "candidate_id": "exploratory_neural",
        "candidate_role": "exploratory",
        "release_broad_label": "Neural candidate",
        "specificity_priority": 95,
    },
}


def row(cell: str, candidate: str, supported: bool, score: float | None = None) -> dict[str, str]:
    value = {
        "cell_id": cell,
        "source_boundary": "b1",
        "source_cluster": "c1",
        "candidate_id": candidate,
        "program_score": str(score if score is not None else (0.4 if supported else 0.0)),
        "normalized_evidence": str(score if score is not None else (0.8 if supported else 0.0)),
        "positive_family_count": "2" if supported else "0",
        "positive_families": "family_a;family_b" if supported else "",
        "positive_gene_count": "4" if supported else "0",
        "family_coherent": "true" if supported else "false",
        "identity_core_direct": "true" if supported else "false",
        "split_discriminator_direct": "true" if supported else "false",
        "release_family_coherent": "true" if supported else "false",
        "candidate_seed": "true" if supported else "false",
        "direct_signal": "0.4" if supported else "0",
        "local_signal": "0.4" if supported else "0",
        "direct_anti_gene_count": "0",
        "direct_anti_family_count": "0",
        "hard_contradiction": "false",
        "ambient_suspect": "false",
        "specificity_priority": str(CANDIDATES[candidate]["specificity_priority"]),
        "cross_resolution_support_count": "3" if supported else "0",
        "local_seed_fraction": "0.8" if supported else "0",
    }
    return value


def score_fixture(assignments: dict[str, str], candidates: list[str]) -> tuple[
    list[dict[str, str]], dict[tuple[str, str], dict[str, str]], list[str]
]:
    rows = []
    for cell, target in assignments.items():
        for candidate in candidates:
            rows.append(row(cell, candidate, candidate == target))
    index, _, universe = index_scores(rows)
    return rows, index, universe


class V22AlgorithmStabilityTests(unittest.TestCase):
    def test_pairwise_mixedness_requires_material_exclusive_direct_cores(self) -> None:
        cells = [f"c{i}" for i in range(20)]
        score_index = {}
        for cell in cells:
            for candidate_id in ("smooth_muscle", "pericyte_mural"):
                supported = (
                    cell in {f"c{i}" for i in range(10)}
                    if candidate_id == "smooth_muscle"
                    else cell in {f"c{i}" for i in range(10, 20)}
                )
                score_index[(cell, candidate_id)] = row(
                    cell, candidate_id, supported
                )
                score_index[(cell, candidate_id)].update({
                    "neighbor_1_boundary": "b1_n1",
                    "neighbor_1_cluster": (
                        "smooth" if int(cell[1:]) < 10 else "pericyte"
                    ),
                    "neighbor_2_boundary": "b1_n2",
                    "neighbor_2_cluster": (
                        "smooth" if int(cell[1:]) < 10 else "pericyte"
                    ),
                })
        result = pairwise_separable_identity_components(
            cells, "smooth_muscle", "pericyte_mural",
            score_index, CANDIDATES,
        )
        self.assertTrue(result["separable"])
        self.assertEqual(result["left_only_n"], 10)
        self.assertEqual(result["right_only_n"], 10)

    def test_nested_shared_program_is_not_pairwise_mixed(self) -> None:
        cells = [f"c{i}" for i in range(20)]
        score_index = {}
        for cell in cells:
            score_index[(cell, "smooth_muscle")] = row(
                cell, "smooth_muscle", True
            )
            score_index[(cell, "pericyte_mural")] = row(
                cell, "pericyte_mural", int(cell[1:]) < 12
            )
        result = pairwise_separable_identity_components(
            cells, "smooth_muscle", "pericyte_mural",
            score_index, CANDIDATES,
        )
        self.assertFalse(result["separable"])
        self.assertEqual(result["right_only_n"], 0)
        self.assertEqual(result["reason"], "coexpressed_or_nested_direct_identity")

    def test_specific_direct_component_can_split_only_inside_generic_parent(self) -> None:
        cells = [f"c{i}" for i in range(20)]
        score_index = {
            (cell, "smooth_muscle"): row(
                cell, "smooth_muscle", int(cell[1:]) < 7
            )
            for cell in cells
        }
        for cell in cells:
            score_index[(cell, "smooth_muscle")].update({
                "neighbor_1_boundary": "b1_n1",
                "neighbor_1_cluster": (
                    "smooth" if int(cell[1:]) < 7 else "remainder"
                ),
                "neighbor_2_boundary": "b1_n2",
                "neighbor_2_cluster": (
                    "smooth" if int(cell[1:]) < 7 else "remainder"
                ),
            })
        result = specific_component_embedded_in_generic_parent(
            cells, "smooth_muscle", score_index, CANDIDATES
        )
        self.assertTrue(result["separable"])
        self.assertEqual(result["direct_identity_n"], 7)
        self.assertEqual(result["complement_n"], 13)

    def test_canonical_cluster_challenger_survives_diluted_required_family_and_somatic_anti(self) -> None:
        aggregate = {
            "candidate_id": "oocyte",
            "available_positive_family_count": "3",
            "group_positive_family_supported_count": "2",
            "group_required_positive_families_pass": "false",
            "observation_identity_core_fraction": "0.46",
            "observation_identity_core_direct_fraction": "0.44",
            "observation_release_family_coherent_fraction": "0.031",
            "positive_marker_detection_fraction": "0.78",
            "marker_deg_log2fc_mean": "2.77",
            "anti_marker_deg_log2fc_mean": "-0.20",
            "hard_contradiction_fraction": "0.93",
        }
        self.assertTrue(canonical_cluster_challenger(aggregate, CANDIDATES["oocyte"]))
        self.assertTrue(group_candidate_detected(aggregate, CANDIDATES["oocyte"]))
        self.assertTrue(independent_group_program(aggregate, CANDIDATES["oocyte"]))

    def test_weak_ordinary_program_cannot_count_as_canonical_oocyte(self) -> None:
        aggregate = {
            "candidate_id": "oocyte",
            "available_positive_family_count": "3",
            "group_positive_family_supported_count": "2",
            "group_required_positive_families_pass": "true",
            "observation_identity_core_fraction": "0.20",
            "observation_identity_core_direct_fraction": "0.05",
            "observation_release_family_coherent_fraction": "0.01",
            "positive_marker_detection_fraction": "0.50",
            "mean_program_score": "0.20",
            "marker_deg_log2fc_mean": "0.80",
            "anti_marker_deg_log2fc_mean": "0.0",
        }
        self.assertFalse(canonical_cluster_challenger(
            aggregate, CANDIDATES["oocyte"]
        ))
        self.assertFalse(group_candidate_detected(
            aggregate, CANDIDATES["oocyte"]
        ))

    def test_fine_subset_cannot_reconstruct_broad_without_parent_identity(self) -> None:
        parent = {
            "candidate_id": "epithelial_mesothelial",
            "candidate_role": "broad",
            "release_broad_label": "Epithelial/mesothelial",
            "required_positive_families": ["family_a", "family_b"],
            "specificity_priority": 90,
        }
        fine = {
            "candidate_id": "epithelial_mesothelial_like",
            "candidate_role": "fine",
            "release_broad_label": "Epithelial/mesothelial",
            "release_fine_label": "Mesothelial-like",
            "context_evidence_candidate_id": "epithelial_mesothelial",
            "specificity_priority": 90,
        }
        catalog = {
            "epithelial_mesothelial": parent,
            "epithelial_mesothelial_like": fine,
        }
        rows = []
        members = [f"e{index}" for index in range(10)]
        for index, cell in enumerate(members):
            fine_row = row(cell, "epithelial_mesothelial", True)
            fine_row.update({
                "candidate_id": "epithelial_mesothelial_like",
                "positive_families": "parent_identity;fine_discriminator",
                "specificity_priority": "90",
            })
            parent_row = row(cell, "epithelial_mesothelial", True)
            parent_row.update({
                "positive_family_count": "1",
                "positive_families": "family_a",
            })
            rows.extend([parent_row, fine_row])
        score_index, _, candidate_ids = index_scores(rows)
        evidence = validate_subset(
            members,
            "epithelial_mesothelial_like",
            score_index,
            candidate_ids,
            catalog=catalog,
            release_level="broad",
        )
        self.assertEqual(evidence["status"], "FAIL")
        self.assertEqual(evidence["parent_identity_status"], "FAIL")
        self.assertEqual(
            evidence["parent_lineage_supported_fraction"], 1.0
        )

    def test_fine_identity_can_reconstruct_sparse_parent_without_declared_required_families(self) -> None:
        parent = {
            "candidate_id": "luteal_steroidogenic",
            "candidate_role": "broad",
            "release_broad_label": "Luteal",
            "specificity_priority": 90,
        }
        fine = {
            "candidate_id": "luteal_early",
            "candidate_role": "fine",
            "release_broad_label": "Luteal",
            "release_fine_label": "Early luteal",
            "context_evidence_candidate_id": "luteal_steroidogenic",
            "specificity_priority": 90,
        }
        catalog = {
            "luteal_steroidogenic": parent,
            "luteal_early": fine,
        }
        members = [f"l{index}" for index in range(10)]
        rows = []
        for index, cell in enumerate(members):
            fine_row = row(cell, "epithelial_mesothelial", True)
            fine_row.update({
                "candidate_id": "luteal_early",
                "positive_families": "parent_identity;fine_discriminator",
                "specificity_priority": "90",
            })
            parent_row = row(cell, "epithelial_mesothelial", index < 4)
            parent_row.update({
                "candidate_id": "luteal_steroidogenic",
                "positive_family_count": "1" if index < 4 else "0",
                "positive_families": "luteal_core" if index < 4 else "",
                "specificity_priority": "90",
            })
            rows.extend([parent_row, fine_row])
        score_index, _, candidate_ids = index_scores(rows)
        evidence = validate_subset(
            members,
            "luteal_early",
            score_index,
            candidate_ids,
            catalog=catalog,
            release_level="broad",
        )
        self.assertEqual(evidence["status"], "PASS")
        self.assertEqual(evidence["parent_identity_status"], "PASS")
        self.assertEqual(
            evidence["parent_lineage_supported_fraction"], 0.4
        )

    def test_generic_stromal_remainder_is_not_an_independent_competitor(self) -> None:
        aggregate = {
            "candidate_id": "stromal_mesenchymal",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "2",
            "group_required_positive_families_pass": "true",
            "observation_identity_core_fraction": "0.97",
            "observation_identity_core_direct_fraction": "0.95",
            "positive_marker_detection_fraction": "0.95",
            "mean_program_score": "0.4",
            "marker_deg_log2fc_mean": "0.0",
            "anti_marker_deg_log2fc_mean": "1.5",
            "positive_marker_pseudobulk_sum": "100",
            "anti_marker_pseudobulk_sum": "100",
            "cross_resolution_stable_fraction": "0.8",
            "hard_contradiction_fraction": "0.96",
        }
        self.assertTrue(group_candidate_detected(aggregate, CANDIDATES["stromal_mesenchymal"]))
        self.assertFalse(independent_group_program(aggregate, CANDIDATES["stromal_mesenchymal"]))
        self.assertFalse(local_split_worthy_group_program(
            aggregate, CANDIDATES["stromal_mesenchymal"]
        ))

    def test_mixture_contradiction_does_not_hide_local_split_program(self) -> None:
        aggregate = {
            "candidate_id": "smooth_muscle",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "2",
            "group_required_positive_families_pass": "true",
            "observation_seed_fraction": "0.24",
            "observation_identity_core_fraction": "0.36",
            "observation_identity_core_direct_fraction": "0.31",
            "observation_release_family_coherent_fraction": "0.30",
            "positive_marker_detection_fraction": "0.85",
            "mean_program_score": "0.28",
            "marker_deg_log2fc_mean": "0.45",
            "anti_marker_deg_log2fc_mean": "0.10",
            "positive_marker_pseudobulk_sum": "120",
            "anti_marker_pseudobulk_sum": "80",
            "cross_resolution_stable_fraction": "0.85",
            "hard_contradiction_fraction": "0.42",
        }
        self.assertTrue(group_candidate_detected(aggregate, CANDIDATES["smooth_muscle"]))
        self.assertTrue(local_split_worthy_group_program(
            aggregate, CANDIDATES["smooth_muscle"]
        ))

    def test_prevalent_program_without_subcluster_deg_triggers_local_separability_check(self) -> None:
        aggregate = {
            "candidate_id": "smooth_muscle",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "2",
            "group_required_positive_families_pass": "true",
            "observation_seed_fraction": "0.40",
            "observation_identity_core_fraction": "0.70",
            "observation_identity_core_direct_fraction": "0.65",
            "observation_release_family_coherent_fraction": "0.65",
            "positive_marker_detection_fraction": "0.90",
            "mean_program_score": "0.35",
            "marker_deg_log2fc_mean": "0.10",
            "anti_marker_deg_log2fc_mean": "0.10",
            "positive_marker_pseudobulk_sum": "200",
            "anti_marker_pseudobulk_sum": "100",
            "cross_resolution_stable_fraction": "0.90",
            "hard_contradiction_fraction": "0.30",
        }
        self.assertTrue(group_candidate_detected(
            aggregate, CANDIDATES["smooth_muscle"]
        ))
        self.assertTrue(local_split_worthy_group_program(
            aggregate, CANDIDATES["smooth_muscle"]
        ))

    def test_single_family_shared_signal_does_not_trigger_local_split(self) -> None:
        aggregate = {
            "candidate_id": "smooth_muscle",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "1",
            "group_required_positive_families_pass": "false",
            "observation_seed_fraction": "0.40",
            "observation_identity_core_fraction": "0.70",
            "observation_identity_core_direct_fraction": "0.65",
            "observation_release_family_coherent_fraction": "0.65",
            "positive_marker_detection_fraction": "0.90",
            "mean_program_score": "0.35",
            "marker_deg_log2fc_mean": "0.10",
            "anti_marker_deg_log2fc_mean": "0.10",
            "positive_marker_pseudobulk_sum": "200",
            "anti_marker_pseudobulk_sum": "100",
            "cross_resolution_stable_fraction": "0.90",
            "hard_contradiction_fraction": "0.30",
        }
        self.assertFalse(local_split_worthy_group_program(
            aggregate, CANDIDATES["smooth_muscle"]
        ))

    def test_subpercent_reproducible_program_is_watch_not_local_split(self) -> None:
        aggregate = {
            "candidate_id": "smooth_muscle",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "2",
            "group_required_positive_families_pass": "true",
            "observation_seed_fraction": "0.006",
            "observation_identity_core_fraction": "0.006",
            "observation_identity_core_direct_fraction": "0.006",
            "observation_release_family_coherent_fraction": "0.006",
            "positive_marker_detection_fraction": "0.10",
            "mean_program_score": "0.05",
            "marker_deg_log2fc_mean": "0.80",
            "anti_marker_deg_log2fc_mean": "0.10",
            "positive_marker_pseudobulk_sum": "20",
            "anti_marker_pseudobulk_sum": "1",
            "cross_resolution_stable_fraction": "0.80",
            "hard_contradiction_fraction": "0",
        }
        self.assertTrue(group_candidate_detected(
            aggregate, CANDIDATES["smooth_muscle"]
        ))
        self.assertTrue(rare_group_program_watch(
            aggregate, CANDIDATES["smooth_muscle"]
        ))
        self.assertFalse(local_split_worthy_group_program(
            aggregate, CANDIDATES["smooth_muscle"]
        ))

    def test_generic_parent_requires_group_coherence_but_not_low_aggregate_anti(self) -> None:
        aggregate = {
            "candidate_id": "stromal_mesenchymal",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "2",
            "group_required_positive_families_pass": "true",
            "observation_identity_core_fraction": "0.98",
            "observation_identity_core_direct_fraction": "0.96",
            "group_positive_family_mean_fraction": "0.95",
            "positive_marker_detection_fraction": "0.98",
            "mean_program_score": "0.40",
            "marker_deg_log2fc_mean": "0.05",
            "anti_marker_deg_log2fc_mean": "0.60",
            "positive_marker_pseudobulk_sum": "300",
            "anti_marker_pseudobulk_sum": "100",
            "cross_resolution_stable_fraction": "0.90",
            "spatial_group_connectivity_fraction": "0.80",
            "hard_contradiction_fraction": "0.55",
        }
        self.assertTrue(dominant_generic_remainder_group(
            aggregate, CANDIDATES["stromal_mesenchymal"]
        ))
        disconnected = dict(aggregate, spatial_group_connectivity_fraction="0")
        self.assertFalse(dominant_generic_remainder_group(
            disconnected, CANDIDATES["stromal_mesenchymal"]
        ))

    def test_fine_candidate_inherits_parent_never_expand_strategy(self) -> None:
        catalog = catalog_candidates({
            "candidate_boundaries": [{
                "candidate_id": "epithelial_mesothelial",
                "candidate_role": "broad",
                "release_broad_label": "Epithelial/mesothelial",
                "release_fine_label": "",
                "writeback_strategy": "candidate_local_component_never_parent_expansion",
            }],
            "machine_actionable_fine_candidate_catalog": {
                "epithelial_mesothelial": [{
                    "candidate_id": "epithelial_surface",
                    "release_label": "Surface epithelial",
                    "parent_release_label": "Epithelial/mesothelial",
                }]
            },
        })
        fine = catalog["epithelial_surface"]
        self.assertEqual(
            effective_broad_writeback_strategy(fine),
            "candidate_local_component_never_parent_expansion",
        )

    def test_canonical_component_validator_recomputes_bounded_exception(self) -> None:
        members = [f"o{i}" for i in range(10)]
        background = [f"s{i}" for i in range(20)]
        candidate = {
            **CANDIDATES["oocyte"],
            "required_positive_families": ["germline_identity", "maternal_ooplasm"],
            "seed_required_positive_families": ["germline_identity", "maternal_ooplasm"],
        }
        rows = []
        for cell in members:
            value = row(cell, "oocyte", True)
            value.update(
                positive_families="germline_identity;maternal_ooplasm",
                identity_core_coherent="true",
                hard_contradiction="true",
                direct_anti_gene_count="4",
                direct_anti_family_count="2",
            )
            rows.append(value)
        for cell in background:
            rows.append(row(cell, "oocyte", False))
        index, _, _ = index_scores(rows)
        aggregate = {
            "available_positive_family_count": "3",
            "group_positive_family_supported_count": "2",
            "observation_identity_core_fraction": "0.46",
            "observation_identity_core_direct_fraction": "0.44",
            "observation_release_family_coherent_fraction": "0.031",
            "positive_marker_detection_fraction": "0.78",
            "marker_deg_log2fc_mean": "2.77",
            "anti_marker_deg_log2fc_mean": "-0.20",
        }
        evidence = validate_canonical_identity_component(
            members, members + background, "oocyte", index, candidate, aggregate
        )
        self.assertEqual(evidence["status"], "PASS")
        self.assertGreater(evidence["contradiction_fraction"], 0.90)

    def test_sparse_tail_in_high_purity_group_inherits_broad(self) -> None:
        cells = [f"s{i}" for i in range(10)]
        rows, index, candidates = score_fixture(
            {cell: "smooth_muscle" if i < 8 else "" for i, cell in enumerate(cells)},
            ["smooth_muscle", "stromal_mesenchymal"],
        )
        parent, evidence = choose_group_parent(
            cells, candidates, index,
            {key: CANDIDATES[key] for key in candidates},
            preferred_parent="smooth_muscle",
        )
        self.assertEqual(parent, "smooth_muscle")
        self.assertEqual(evidence["status"], "PASS")

    def test_certified_generic_remainder_returns_broad_unless_specific_blocker_remains(self) -> None:
        cells = [f"r{i}" for i in range(10)]
        rows = []
        for i, cell in enumerate(cells):
            parent = row(cell, "stromal_mesenchymal", i < 5)
            if i < 5:
                parent.update(
                    hard_contradiction="true",
                    direct_anti_gene_count="2",
                    direct_anti_family_count="2",
                )
            rows.append(parent)
            rows.append(row(cell, "smooth_muscle", False))
        index, _, universe = index_scores(rows)
        catalog = {
            key: CANDIDATES[key]
            for key in ("stromal_mesenchymal", "smooth_muscle")
        }
        parent, evidence = choose_group_parent(
            cells,
            universe,
            index,
            catalog,
            preferred_parent="stromal_mesenchymal",
            blocker_candidate_ids={"stromal_mesenchymal"},
            parent_candidate_ids={"stromal_mesenchymal"},
            certified_generic_parent_ids={"stromal_mesenchymal"},
        )
        self.assertEqual(parent, "stromal_mesenchymal")
        self.assertEqual(evidence["reason"], "certified_generic_remainder_parent")

        for i, cell in enumerate(cells[:3]):
            competitor = index[(cell, "smooth_muscle")]
            competitor.update(row(cell, "smooth_muscle", True))
        parent, evidence = choose_group_parent(
            cells,
            universe,
            index,
            catalog,
            preferred_parent="stromal_mesenchymal",
            blocker_candidate_ids={"stromal_mesenchymal", "smooth_muscle"},
            parent_candidate_ids={"stromal_mesenchymal"},
            certified_generic_parent_ids={"stromal_mesenchymal"},
        )
        self.assertEqual(parent, "")
        self.assertEqual(evidence["reason"], "embedded_competing_program")

    def test_post_split_sparse_generic_remainder_uses_aggregate_backbone(self) -> None:
        cells = [f"g{i}" for i in range(10)]
        rows = []
        for i, cell in enumerate(cells):
            stromal = row(cell, "stromal_mesenchymal", i < 3)
            if i >= 3:
                stromal.update(
                    hard_contradiction="true",
                    direct_anti_gene_count="2",
                    direct_anti_family_count="2",
                )
            rows.append(stromal)
            rows.append(row(cell, "smooth_muscle", False))
        index, _, universe = index_scores(rows)
        catalog = {
            key: CANDIDATES[key]
            for key in ("stromal_mesenchymal", "smooth_muscle")
        }
        aggregate = {
            "candidate_id": "stromal_mesenchymal",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "2",
            "group_required_positive_families_pass": "true",
            "observation_identity_core_fraction": "0.98",
            "observation_identity_core_direct_fraction": "0.96",
            "group_positive_family_mean_fraction": "0.95",
            "positive_marker_detection_fraction": "0.98",
            "mean_program_score": "0.40",
            "marker_deg_log2fc_mean": "0.05",
            "anti_marker_deg_log2fc_mean": "0.60",
            "positive_marker_pseudobulk_sum": "300",
            "anti_marker_pseudobulk_sum": "100",
            "cross_resolution_stable_fraction": "0.90",
            "spatial_group_connectivity_fraction": "0.80",
        }
        blockers, parents, preferred, audit = remainder_candidate_programs(
            cells,
            "b1",
            "c1",
            universe,
            index,
            {("selected", "b1", "c1", "stromal_mesenchymal"): aggregate},
            catalog,
        )
        stromal_audit = next(
            item for item in audit
            if item["candidate_id"] == "stromal_mesenchymal"
        )
        self.assertAlmostEqual(stromal_audit["supported_fraction"], 0.30)
        self.assertTrue(stromal_audit["generic_remainder_supported"])
        self.assertEqual(blockers, {"stromal_mesenchymal"})
        self.assertEqual(parents, {"stromal_mesenchymal"})
        self.assertEqual(preferred, "stromal_mesenchymal")

    def test_post_split_dominant_generic_tail_returns_parent_without_pure_source_group(self) -> None:
        cells = [f"t{i}" for i in range(10)]
        rows = []
        for i, cell in enumerate(cells):
            stromal = row(cell, "stromal_mesenchymal", i < 8)
            if i >= 8:
                stromal.update(
                    hard_contradiction="true",
                    direct_anti_gene_count="2",
                    direct_anti_family_count="2",
                )
            rows.append(stromal)
            rows.append(row(cell, "smooth_muscle", False))
        index, _, universe = index_scores(rows)
        catalog = {
            key: CANDIDATES[key]
            for key in ("stromal_mesenchymal", "smooth_muscle")
        }
        mixed_source_aggregate = {
            "candidate_id": "stromal_mesenchymal",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "2",
            "group_required_positive_families_pass": "true",
            "observation_identity_core_fraction": "0.50",
            "observation_identity_core_direct_fraction": "0.50",
            "group_positive_family_mean_fraction": "0.50",
            "positive_marker_detection_fraction": "0.50",
            "mean_program_score": "0.25",
            "marker_deg_log2fc_mean": "0.10",
            "anti_marker_deg_log2fc_mean": "0.20",
            "cross_resolution_stable_fraction": "0.80",
            "spatial_group_connectivity_fraction": "0.80",
        }
        blockers, parents, preferred, audit = remainder_candidate_programs(
            cells,
            "b1",
            "c1",
            universe,
            index,
            {
                ("selected", "b1", "c1", "stromal_mesenchymal"):
                    mixed_source_aggregate,
            },
            catalog,
        )
        stromal_audit = next(
            item for item in audit
            if item["candidate_id"] == "stromal_mesenchymal"
        )
        self.assertFalse(
            stromal_audit["aggregate_generic_remainder_supported"]
        )
        self.assertTrue(
            stromal_audit["remainder_dominant_generic_supported"]
        )
        self.assertEqual(blockers, {"stromal_mesenchymal"})
        self.assertEqual(parents, {"stromal_mesenchymal"})
        self.assertEqual(preferred, "stromal_mesenchymal")

        parent, evidence = choose_group_parent(
            cells,
            universe,
            index,
            catalog,
            preferred_parent=preferred,
            blocker_candidate_ids=blockers,
            parent_candidate_ids=parents,
            certified_generic_parent_ids=parents,
        )
        self.assertEqual(parent, "stromal_mesenchymal")
        self.assertEqual(evidence["reason"], "certified_generic_remainder_parent")

    def test_aggregate_only_shared_program_does_not_block_generic_remainder(self) -> None:
        cells = [f"w{i}" for i in range(10)]
        rows = []
        for i, cell in enumerate(cells):
            rows.append(row(cell, "stromal_mesenchymal", i < 3))
            challenger = row(cell, "exploratory_neural", i < 4)
            if i >= 4 and i < 8:
                challenger.update(
                    hard_contradiction="true",
                    direct_anti_gene_count="2",
                    direct_anti_family_count="2",
                )
            rows.append(challenger)
        index, _, universe = index_scores(rows)
        catalog = {
            key: CANDIDATES[key]
            for key in ("stromal_mesenchymal", "exploratory_neural")
        }
        generic_aggregate = {
            "candidate_id": "stromal_mesenchymal",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "2",
            "group_required_positive_families_pass": "true",
            "observation_identity_core_fraction": "0.98",
            "observation_identity_core_direct_fraction": "0.96",
            "group_positive_family_mean_fraction": "0.95",
            "positive_marker_detection_fraction": "0.98",
            "mean_program_score": "0.40",
            "marker_deg_log2fc_mean": "0.05",
            "anti_marker_deg_log2fc_mean": "0.60",
            "positive_marker_pseudobulk_sum": "300",
            "anti_marker_pseudobulk_sum": "100",
            "cross_resolution_stable_fraction": "0.90",
            "spatial_group_connectivity_fraction": "0.80",
        }
        shared_aggregate = {
            "candidate_id": "exploratory_neural",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "2",
            "observation_identity_core_fraction": "0.40",
            "observation_identity_core_direct_fraction": "0.40",
            "observation_coherent_fraction": "0.40",
            "observation_seed_fraction": "0.40",
            "positive_marker_detection_fraction": "0.40",
            "mean_program_score": "0.20",
            "marker_deg_log2fc_mean": "0.80",
            "anti_marker_deg_log2fc_mean": "0.00",
        }
        blockers, parents, preferred, audit = remainder_candidate_programs(
            cells,
            "b1",
            "c1",
            universe,
            index,
            {
                ("selected", "b1", "c1", "stromal_mesenchymal"): generic_aggregate,
                ("selected", "b1", "c1", "exploratory_neural"): shared_aggregate,
            },
            catalog,
        )
        challenger_audit = next(
            item for item in audit
            if item["candidate_id"] == "exploratory_neural"
        )
        self.assertTrue(challenger_audit["aggregate_supported"])
        self.assertTrue(challenger_audit["aggregate_only_watch"])
        self.assertFalse(challenger_audit["credible_blocker"])
        self.assertEqual(blockers, {"stromal_mesenchymal"})
        self.assertEqual(parents, {"stromal_mesenchymal"})
        self.assertEqual(preferred, "stromal_mesenchymal")

        blockers, parents, preferred, audit = remainder_candidate_programs(
            cells,
            "b1",
            "c1",
            universe,
            index,
            {
                ("selected", "b1", "c1", "stromal_mesenchymal"): generic_aggregate,
                ("selected", "b1", "c1", "exploratory_neural"): shared_aggregate,
            },
            catalog,
            {"exploratory_neural"},
        )
        self.assertIn("exploratory_neural", blockers)
        self.assertNotIn("stromal_mesenchymal", parents)
        self.assertEqual(preferred, "")

    def test_residual_separable_specific_component_blocks_generic_closure(self) -> None:
        cells = [f"m{i}" for i in range(10)]
        rows = []
        for i, cell in enumerate(cells):
            rows.append(row(cell, "stromal_mesenchymal", i < 3))
            smooth = row(cell, "smooth_muscle", i < 3)
            if 3 <= i < 7:
                smooth.update(
                    hard_contradiction="true",
                    direct_anti_gene_count="2",
                    direct_anti_family_count="2",
                )
            rows.append(smooth)
        index, _, universe = index_scores(rows)
        catalog = {
            key: CANDIDATES[key]
            for key in ("stromal_mesenchymal", "smooth_muscle")
        }
        aggregate = {
            "candidate_id": "stromal_mesenchymal",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "2",
            "group_required_positive_families_pass": "true",
            "observation_identity_core_fraction": "0.98",
            "observation_identity_core_direct_fraction": "0.96",
            "group_positive_family_mean_fraction": "0.95",
            "positive_marker_detection_fraction": "0.98",
            "mean_program_score": "0.40",
            "marker_deg_log2fc_mean": "0.05",
            "anti_marker_deg_log2fc_mean": "0.60",
            "positive_marker_pseudobulk_sum": "300",
            "anti_marker_pseudobulk_sum": "100",
            "cross_resolution_stable_fraction": "0.90",
            "spatial_group_connectivity_fraction": "0.80",
        }
        blockers, parents, preferred, audit = remainder_candidate_programs(
            cells,
            "b1",
            "c1",
            universe,
            index,
            {("selected", "b1", "c1", "stromal_mesenchymal"): aggregate},
            catalog,
            {"smooth_muscle"},
        )
        smooth_audit = next(
            item for item in audit
            if item["candidate_id"] == "smooth_muscle"
        )
        stromal_audit = next(
            item for item in audit
            if item["candidate_id"] == "stromal_mesenchymal"
        )
        self.assertTrue(smooth_audit["residual_separable_component"])
        self.assertTrue(smooth_audit["credible_blocker"])
        self.assertFalse(stromal_audit["generic_remainder_supported"])
        self.assertIn("smooth_muscle", blockers)
        self.assertEqual(parents, set())
        self.assertEqual(preferred, "")

    def test_component_history_alone_cannot_block_exact_remainder(self) -> None:
        cells = [f"h{i}" for i in range(10)]
        rows = []
        for i, cell in enumerate(cells):
            rows.append(row(cell, "stromal_mesenchymal", False))
            smooth = row(cell, "smooth_muscle", i < 3)
            if 3 <= i < 7:
                smooth.update(
                    hard_contradiction="true",
                    direct_anti_gene_count="2",
                    direct_anti_family_count="2",
                )
            rows.append(smooth)
        index, _, universe = index_scores(rows)
        catalog = {
            key: CANDIDATES[key]
            for key in ("stromal_mesenchymal", "smooth_muscle")
        }
        candidates, ambiguous, audit = validated_residual_component_candidates(
            cells,
            {"smooth_muscle": set(cells)},
            "b1",
            "c1",
            universe,
            index,
            catalog,
        )
        self.assertEqual(candidates, set())
        self.assertEqual(ambiguous, set())
        self.assertEqual(audit[0]["status"], "FAIL")

        _, clean_index, clean_universe = score_fixture(
            {
                cell: "smooth_muscle" if i < 8 else ""
                for i, cell in enumerate(cells)
            },
            ["stromal_mesenchymal", "smooth_muscle"],
        )
        candidates, ambiguous, audit = validated_residual_component_candidates(
            cells,
            {"smooth_muscle": set(cells)},
            "b1",
            "c1",
            clean_universe,
            clean_index,
            catalog,
        )
        self.assertEqual(candidates, {"smooth_muscle"})
        self.assertEqual(ambiguous, set(cells))
        self.assertEqual(audit[0]["status"], "PASS")

    def test_fine_or_exploratory_program_blocks_but_never_becomes_parent(self) -> None:
        cells = [f"x{i}" for i in range(10)]
        rows = []
        for i, cell in enumerate(cells):
            rows.append(row(cell, "stromal_mesenchymal", i < 8))
            rows.append(row(cell, "exploratory_neural", i < 3))
        index, _, universe = index_scores(rows)
        catalog = {
            key: CANDIDATES[key]
            for key in ("stromal_mesenchymal", "exploratory_neural")
        }
        parent, evidence = choose_group_parent(
            cells,
            universe,
            index,
            catalog,
            preferred_parent="stromal_mesenchymal",
            blocker_candidate_ids={"stromal_mesenchymal", "exploratory_neural"},
            parent_candidate_ids={"stromal_mesenchymal"},
        )
        self.assertEqual(parent, "")
        self.assertEqual(evidence["status"], "UNRESOLVED")

        fine_cells = [f"f{i}" for i in range(5)]
        _, fine_index, fine_universe = score_fixture(
            {cell: "lymphatic_endothelial" for cell in fine_cells},
            ["lymphatic_endothelial"],
        )
        parent, evidence = choose_group_parent(
            fine_cells,
            fine_universe,
            fine_index,
            {"lymphatic_endothelial": CANDIDATES["lymphatic_endothelial"]},
            parent_candidate_ids={"lymphatic_endothelial"},
        )
        self.assertEqual(parent, "")
        self.assertEqual(evidence["reason"], "no_releasable_positive_program")

    def test_credible_mixed_blocker_is_not_hidden_by_group_contradiction(self) -> None:
        cells = [f"b{i}" for i in range(10)]
        rows = []
        for i, cell in enumerate(cells):
            rows.append(row(cell, "granulosa", i < 9))
            smooth = row(cell, "smooth_muscle", i < 8)
            if i >= 8:
                smooth["hard_contradiction"] = "true"
                smooth["direct_anti_gene_count"] = "2"
                smooth["direct_anti_family_count"] = "2"
            rows.append(smooth)
        index, _, universe = index_scores(rows)
        catalog = {
            key: CANDIDATES[key] for key in ("granulosa", "smooth_muscle")
        }
        parent, evidence = choose_group_parent(
            cells,
            universe,
            index,
            catalog,
            blocker_candidate_ids={"granulosa", "smooth_muscle"},
            parent_candidate_ids={"granulosa"},
        )
        self.assertEqual(parent, "")
        self.assertEqual(evidence["reason"], "embedded_competing_program")

    def test_sparse_oocyte_members_inherit_strong_canonical_highres_cluster(self) -> None:
        cells = [f"o{i}" for i in range(100)]
        rows = []
        for i, cell in enumerate(cells):
            target = row(cell, "oocyte", i < 10)
            if i >= 10:
                target["hard_contradiction"] = "true"
                target["direct_anti_gene_count"] = "2"
                target["direct_anti_family_count"] = "2"
            rows.append(target)
        index, _, universe = index_scores(rows)
        aggregate = {
            "observation_coherent_fraction": "0.893",
            "observation_seed_fraction": "0.037",
            "hard_contradiction_fraction": "0.939",
            "mean_program_score": "0.051",
            "positive_marker_detection_fraction": "0.954",
            "marker_deg_log2fc_mean": "4.03",
            "anti_marker_deg_log2fc_mean": "-0.50",
        }
        self.assertTrue(
            aggregate_program_supported(
                aggregate,
                CANDIDATES["oocyte"],
                {"status": "PASS"},
            )
        )
        evidence = validate_subset(
            cells,
            "oocyte",
            index,
            universe,
            catalog={"oocyte": CANDIDATES["oocyte"]},
            aggregate_evidence=aggregate,
        )
        self.assertEqual(evidence["status"], "PASS")
        contaminated_identity = row("o-core", "oocyte", True)
        contaminated_identity.update(
            positive_families="germline_identity;maternal_ooplasm",
            identity_core_coherent="true",
            hard_contradiction="true",
            direct_anti_gene_count="4",
            direct_anti_family_count="2",
            candidate_seed="false",
        )
        self.assertTrue(
            candidate_core_seed(
                contaminated_identity,
                {
                    **CANDIDATES["oocyte"],
                    "seed_required_positive_families": [
                        "germline_identity", "maternal_ooplasm"
                    ],
                },
            )
        )

    def test_r_subset_deriver_consumes_scorer_frozen_identity_core(self) -> None:
        text = (SCRIPTS / "derive_candidate_local_subsets.R").read_text()
        family_body = text.split("group_family_evidence <- function", 1)[1].split(
            "aggregate_row_for <- function", 1
        )[0]
        identity_body = text.split("identity_core_mask <- function", 1)[1].split(
            "minimum_identity_core_fraction <- function", 1
        )[0]
        self.assertIn("candidate_rows$identity_core_coherent", identity_body)
        self.assertIn("candidate_rows$candidate_seed", identity_body)
        self.assertIn("candidate_rows$program_score >= 0.02", identity_body)
        self.assertIn("candidate_rows$direct_signal > 0", identity_body)
        self.assertNotIn("seed_required_positive_families", identity_body)
        self.assertIn("required_positive_families", family_body)
        self.assertIn("or_values", family_body)

    def test_aggregate_subset_cannot_lower_identity_support_to_twenty_five_percent(self) -> None:
        cells = [f"s{i}" for i in range(20)]
        rows = [row(cell, "smooth_muscle", True) for cell in cells]
        index, _, universe = index_scores(rows)
        evidence = validate_subset(
            cells,
            "smooth_muscle",
            index,
            universe,
            catalog={"smooth_muscle": CANDIDATES["smooth_muscle"]},
            aggregate_evidence={
                "observation_coherent_fraction": "0.90",
                "observation_identity_core_fraction": "0.25",
                "hard_contradiction_fraction": "0.00",
                "mean_program_score": "0.20",
                "positive_marker_detection_fraction": "0.90",
                "marker_deg_log2fc_mean": "2.0",
                "anti_marker_deg_log2fc_mean": "0.0",
            },
        )
        self.assertEqual(evidence["status"], "FAIL")
        self.assertEqual(evidence["lineage_supported_fraction"], 0.25)

    def test_negative_deg_prevalence_cannot_wholesale_return_candidate(self) -> None:
        aggregate = {
            "observation_coherent_fraction": "0.90",
            "observation_seed_fraction": "0.80",
            "mean_program_score": "0.30",
            "positive_marker_detection_fraction": "0.95",
            "marker_deg_log2fc_mean": "-1.20",
            "anti_marker_deg_log2fc_mean": "-0.10",
        }
        self.assertFalse(
            aggregate_program_supported(
                aggregate,
                CANDIDATES["stromal_mesenchymal"],
                {"status": "PASS"},
            )
        )

    def test_one_anti_is_soft_but_multigene_direct_anti_is_hard(self) -> None:
        one = row("a", "smooth_muscle", True)
        one.update(hard_contradiction="true", direct_anti_gene_count="1")
        two = dict(one, direct_anti_gene_count="2")
        self.assertFalse(hard_contradiction(one))
        self.assertFalse(hard_contradiction(two))
        one_family = dict(two, direct_anti_family_count="1")
        two_families = dict(two, direct_anti_family_count="2")
        self.assertFalse(hard_contradiction(one_family))
        self.assertTrue(hard_contradiction(two_families))

    def test_sparse_observations_jointly_validate_two_families_at_group_level(self) -> None:
        cells = [f"g{i}" for i in range(10)]
        rows = []
        for i, cell in enumerate(cells):
            target = row(cell, "smooth_muscle", True)
            target.update(
                positive_family_count="1",
                positive_families="contractile_a" if i < 5 else "contractile_b",
                release_family_coherent="false",
            )
            competitor = row(cell, "stromal_mesenchymal", False)
            rows.extend([target, competitor])
        index, _, universe = index_scores(rows)
        evidence = validate_subset(
            cells,
            "smooth_muscle",
            index,
            universe,
            catalog={
                "smooth_muscle": CANDIDATES["smooth_muscle"],
                "stromal_mesenchymal": CANDIDATES["stromal_mesenchymal"],
            },
        )
        self.assertEqual(evidence["status"], "PASS")
        self.assertEqual(
            evidence["supported_families"], "contractile_a;contractile_b"
        )

    def test_mixed_cluster_anti_does_not_hide_candidate_local_program(self) -> None:
        evidence = {
            "observation_seed_fraction": "0.04",
            "observation_coherent_fraction": "0.20",
            "hard_contradiction_fraction": "0.80",
            "marker_deg_log2fc_mean": "1.2",
            "positive_marker_detection_fraction": "0.30",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "2",
            "mean_program_score": "0.20",
        }
        self.assertTrue(candidate_detected(evidence))

    def test_rare_coherent_program_survives_three_percent_cluster_dilution(self) -> None:
        evidence = {
            "observation_seed_fraction": "0.012",
            "observation_coherent_fraction": "0.80",
            "hard_contradiction_fraction": "0.70",
            "marker_deg_log2fc_mean": "2.0",
            "positive_marker_detection_fraction": "0.25",
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "2",
            "mean_program_score": "0.20",
        }
        self.assertTrue(candidate_detected(evidence))

    def test_local_only_one_family_signal_is_not_a_group_identity(self) -> None:
        evidence = {
            "available_positive_family_count": "2",
            "group_positive_family_supported_count": "1",
            "observation_identity_core_fraction": "0.80",
            "observation_identity_core_direct_fraction": "0",
            "observation_coherent_fraction": "1",
            "mean_program_score": "0.40",
            "positive_marker_detection_fraction": "0.90",
            "marker_deg_log2fc_mean": "1.50",
        }
        self.assertFalse(group_candidate_detected(evidence, {}))

    def test_missing_required_group_family_blocks_identity_release(self) -> None:
        evidence = {
            "available_positive_family_count": "3",
            "group_positive_family_supported_count": "2",
            "group_required_positive_families_pass": "false",
            "observation_identity_core_fraction": "0.80",
            "observation_identity_core_direct_fraction": "0.60",
            "mean_program_score": "0.40",
            "positive_marker_detection_fraction": "0.90",
            "marker_deg_log2fc_mean": "1.50",
        }
        self.assertFalse(group_candidate_detected(evidence, {}))

    def test_r_scorer_binds_seed_to_declared_identity_core(self) -> None:
        source = (SCRIPTS / "run_observation_lineage_scoring.R").read_text()
        self.assertIn(
            "catalog_by_id[[candidate_id]]$seed_required_positive_families",
            source,
        )
        seed_block = source.split(
            "# Candidate-local seeds are computed independently", 1
        )[1].split("candidate_local_seed_fraction", 1)[0]
        self.assertIn("candidate_identity_core_coherent", seed_block)
        self.assertNotIn("candidate_family_coherent &", seed_block)
        identity_block = source.split(
            "seed_families <-", 1
        )[1].split("required_families <-", 1)[0]
        self.assertIn(") >= 1L", identity_block)
        self.assertNotIn(") == length(seed_columns)", identity_block)

    def test_mixed_cluster_splits_smooth_vascular_and_stromal_remainder(self) -> None:
        assignments = {
            **{f"m{i}": "smooth_muscle" for i in range(4)},
            **{f"v{i}": "pericyte_mural" for i in range(4)},
            **{f"s{i}": "stromal_mesenchymal" for i in range(4)},
        }
        _, index, universe = score_fixture(
            assignments, ["smooth_muscle", "pericyte_mural", "stromal_mesenchymal"]
        )
        smooth = validate_subset(
            [f"m{i}" for i in range(4)], "smooth_muscle", index, universe
        )
        vascular = validate_subset(
            [f"v{i}" for i in range(4)], "pericyte_mural", index, universe
        )
        parent, evidence = choose_group_parent(
            [f"s{i}" for i in range(4)], universe, index,
            {key: CANDIDATES[key] for key in universe},
            preferred_parent="stromal_mesenchymal",
        )
        self.assertEqual(smooth["status"], "PASS")
        self.assertEqual(vascular["status"], "PASS")
        self.assertEqual(parent, "stromal_mesenchymal")
        self.assertEqual(evidence["status"], "PASS")

    def test_aggregate_winner_does_not_suppress_independent_candidates(self) -> None:
        for candidate in ("theca_steroidogenic", "epithelial_mesothelial", "oocyte"):
            alternative = row("x", candidate, True, 0.25)
            aggregate = row("x", "stromal_mesenchymal", True, 0.9)
            self.assertTrue(supported_seed(alternative))
            self.assertEqual(rank_supported([aggregate, alternative])[0]["candidate_id"],
                             "stromal_mesenchymal")

    def test_zero_signal_never_wins_over_supported_negative_adjusted_candidate(self) -> None:
        zero = row("x", "stromal_mesenchymal", False, 0.0)
        supported = row("x", "smooth_muscle", True, -0.1)
        ranked = rank_supported([zero, supported])
        self.assertEqual([item["candidate_id"] for item in ranked], ["smooth_muscle"])

    def test_small_spatial_epithelial_program_survives_large_stromal_cluster(self) -> None:
        assignments = {
            **{f"e{i}": "epithelial_mesothelial" for i in range(4)},
            **{f"s{i}": "stromal_mesenchymal" for i in range(96)},
        }
        _, index, universe = score_fixture(
            assignments, ["epithelial_mesothelial", "stromal_mesenchymal"]
        )
        evidence = validate_subset(
            [f"e{i}" for i in range(4)],
            "epithelial_mesothelial", index, universe,
        )
        self.assertEqual(evidence["status"], "PASS")

    def test_fine_candidate_does_not_compete_with_parent_at_broad_release(self) -> None:
        catalog = {
            "granulosa": {
                "candidate_id": "granulosa", "candidate_role": "broad",
                "release_broad_label": "Granulosa", "release_fine_label": "",
            },
            "granulosa_mural": {
                "candidate_id": "granulosa_mural", "candidate_role": "fine",
                "release_broad_label": "Granulosa",
                "release_fine_label": "Mural/estrogenic granulosa",
            },
        }
        rows = []
        for cell in ("a", "b", "c", "d"):
            parent = row(cell, "stromal_mesenchymal", True)
            parent["candidate_id"] = "granulosa"
            fine = dict(parent, candidate_id="granulosa_mural")
            rows.extend([parent, fine])
        index, _, universe = index_scores(rows)
        evidence = validate_subset(
            ["a", "b", "c", "d"], "granulosa", index, universe,
            catalog=catalog, release_level="broad",
        )
        self.assertEqual(evidence["status"], "PASS")
        self.assertEqual(evidence["strongest_competing_fraction"], 0)

    def test_fine_candidate_parent_signal_cannot_seed_fine_component(self) -> None:
        candidate = {
            "candidate_id": "granulosa_cumulus_like",
            "candidate_role": "fine",
            "release_broad_label": "Granulosa",
            "release_fine_label": "Cumulus-like granulosa",
            "seed_required_positive_families": ["fine_discriminator"],
        }
        parent_only = row("g", "smooth_muscle", True)
        parent_only.update(
            candidate_id="granulosa_cumulus_like",
            positive_families="parent_identity",
        )
        discriminator = dict(
            parent_only, positive_families="fine_discriminator"
        )
        self.assertFalse(candidate_core_seed(parent_only, candidate))
        self.assertTrue(candidate_core_seed(discriminator, candidate))

        installed = json.loads((
            ROOT
            / "annotate-spatial-transcriptomics/references/profiles/"
            "sheep_ovary_candidate_lineage_catalog.json"
        ).read_text(encoding="utf-8"))
        cumulus = installed["machine_actionable_fine_candidate_catalog"][
            "granulosa"
        ][0]
        self.assertEqual(
            cumulus["fine_positive_families"]["fine_discriminator_core"],
            ["HAS2", "PTX3", "TNFAIP6"],
        )
        self.assertNotIn(
            "VCAN",
            cumulus["fine_positive_families"]["fine_discriminator_core"],
        )

    def test_contractile_support_cannot_bridge_smooth_identity_cores(self) -> None:
        candidate = {
            **CANDIDATES["smooth_muscle"],
            "seed_required_positive_families": ["mature_contractile_core"],
        }
        support_only = row("s", "smooth_muscle", True)
        support_only["positive_families"] = "contractile_structure_support"
        mature = dict(
            support_only,
            positive_families="mature_contractile_core",
        )
        self.assertFalse(candidate_core_seed(support_only, candidate))
        self.assertTrue(candidate_core_seed(mature, candidate))

    def test_generated_fine_candidate_inherits_parent_lock_and_context(self) -> None:
        catalog = catalog_candidates({
            "candidate_boundaries": [{
                "candidate_id": "luteal",
                "candidate_role": "broad",
                "release_broad_label": "Luteal",
                "release_fine_label": "",
                "formal_context_evidence_required": True,
            }],
            "machine_actionable_fine_candidate_catalog": {
                "luteal": [{
                    "candidate_id": "luteal_late",
                    "release_label": "Late luteal",
                    "parent_release_label": "Luteal",
                    "profile_program": "lineages.luteal.late",
                }]
            },
        })
        fine = catalog["luteal_late"]
        self.assertEqual(
            fine["required_positive_families"],
            ["parent_identity", "fine_discriminator"],
        )
        self.assertEqual(
            fine["seed_required_positive_families"],
            ["fine_discriminator"],
        )
        self.assertTrue(fine["formal_context_evidence_required"])
        self.assertEqual(fine["context_evidence_candidate_id"], "luteal")
        self.assertFalse(candidate_can_release(fine))
        fine["_context_ok"] = True
        self.assertTrue(candidate_can_release(fine))

    def test_fine_overlap_common_parent_retains_contributing_provenance(self) -> None:
        proposals = [("subset_b", "fine_b"), ("subset_a", "fine_a")]
        self.assertEqual(
            provenance_subset(proposals, "common_broad_parent"),
            "subset_a",
        )

    def test_external_context_candidate_requires_bound_support(self) -> None:
        candidate = {
            "candidate_role": "broad",
            "release_broad_label": "Luteal",
            "formal_context_evidence_required": True,
        }
        self.assertFalse(candidate_can_release(candidate))
        candidate["_context_ok"] = True
        self.assertTrue(candidate_can_release(candidate))

    def test_stable_unmodeled_program_is_recorded_without_label(self) -> None:
        rows = [{
            "program_id": "p",
            "resolutions": "0.2;0.4",
            "spatially_coherent": "true",
            "coexpressed_gene_count": "4",
            "excluded_program_classes": "",
        }]
        accepted = validate_unmodeled_programs(rows)
        self.assertEqual(accepted, rows)

    def test_empty_fine_census_requires_complete_parent_candidate_audit(self) -> None:
        catalog = {
            "machine_actionable_fine_candidate_catalog": {
                "granulosa": [
                    {
                        "candidate_id": "cumulus",
                        "parent_release_label": "Granulosa",
                    },
                    {
                        "candidate_id": "mural_granulosa",
                        "parent_release_label": "Granulosa",
                    },
                ]
            }
        }
        complete, expected, observed = fine_audit_complete(
            catalog,
            {"Granulosa"},
            [{"parent_broad_label": "Granulosa",
              "candidate_id": "cumulus", "status": "refuted"}],
        )
        self.assertFalse(complete)
        self.assertEqual(len(expected), 2)
        self.assertEqual(len(observed), 1)
        complete, _, _ = fine_audit_complete(
            catalog,
            {"Granulosa"},
            [
                {"parent_broad_label": "Granulosa",
                 "candidate_id": "cumulus", "status": "refuted"},
                {"parent_broad_label": "Granulosa",
                 "candidate_id": "mural_granulosa", "status": "not_evaluable"},
                {"parent_broad_label": "Immune",
                 "candidate_id": "extra", "status": "not_evaluable"},
            ],
        )
        self.assertTrue(complete)

    def test_resolution_selector_uses_nearest_neighbors_at_grid_edge(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            evidence = root / "grid.tsv"
            fields = [
                "resolution", "selection_purpose", "complete_catalog_scanned", "zero_census_audited",
                "catalog_recall", "embedded_program_separation",
                "deg_antideg_coherence", "pseudobulk_coherence",
                "spatial_morphology_coherence", "adjacent_resolution_stability",
                "observation_support_coherence", "technical_fragmentation",
                "directly_resolved_observation_fraction", "mean_identity_margin",
                "mixed_observation_fraction", "unresolved_observation_fraction",
                "state_overfragmentation", "complexity",
            ]
            with evidence.open("w", encoding="utf-8") as handle:
                handle.write("\t".join(fields) + "\n")
                for resolution, score in (
                    (0.2, 0.7), (0.4, 0.8), (0.6, 0.9), (0.8, 1.0)
                ):
                    values = {
                        **{key: score for key in fields},
                        "resolution": resolution,
                        "selection_purpose": "cohort_identity_resolution",
                        "complete_catalog_scanned": "true",
                        "zero_census_audited": "true",
                        "technical_fragmentation": 0,
                        "mixed_observation_fraction": 0,
                        "unresolved_observation_fraction": 0,
                        "state_overfragmentation": 0,
                        "complexity": resolution,
                    }
                    handle.write(
                        "\t".join(str(values[key]) for key in fields) + "\n"
                    )
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "select_lineage_resolution.py"),
                "--grid-evidence", str(evidence),
                "--selection-purpose", "cohort_identity_resolution",
                "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            selected = json.loads(
                (root / "out/resolution_selection.json").read_text()
            )
            self.assertEqual(selected["selected_resolution"], 0.8)
            self.assertEqual(
                selected["selected_and_neighbors"], [0.8, 0.6, 0.4]
            )
            self.assertEqual(
                selected["neighbor_strategy"],
                "nearest_lower_and_higher_else_two_nearest",
            )

    def test_resolution_simplicity_tiebreak_requires_metricwise_equivalence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            evidence = root / "grid.tsv"
            metrics = [
                "catalog_recall", "embedded_program_separation",
                "deg_antideg_coherence", "pseudobulk_coherence",
                "spatial_morphology_coherence",
                "adjacent_resolution_stability",
                "observation_support_coherence",
                "directly_resolved_observation_fraction",
                "mean_identity_margin",
            ]
            fields = [
                "resolution", "selection_purpose", "complete_catalog_scanned",
                "zero_census_audited", *metrics,
                "technical_fragmentation", "state_overfragmentation",
                "mixed_observation_fraction", "unresolved_observation_fraction",
                "complexity",
            ]
            with evidence.open("w", encoding="utf-8") as handle:
                handle.write("\t".join(fields) + "\n")
                for resolution in (0.2, 0.4, 0.6, 0.8):
                    values = {metric: 0.8 for metric in metrics}
                    if resolution == 0.2:
                        values["catalog_recall"] = 0.824
                        values["embedded_program_separation"] = 0.76
                    row_values = {
                        **values,
                        "resolution": resolution,
                        "selection_purpose": "cohort_identity_resolution",
                        "complete_catalog_scanned": "true",
                        "zero_census_audited": "true",
                        "technical_fragmentation": 0,
                        "mixed_observation_fraction": 0,
                        "unresolved_observation_fraction": 0,
                        "state_overfragmentation": 0,
                        "complexity": resolution,
                    }
                    handle.write(
                        "\t".join(str(row_values[field]) for field in fields)
                        + "\n"
                    )
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "select_lineage_resolution.py"),
                "--grid-evidence", str(evidence),
                "--selection-purpose", "cohort_identity_resolution",
                "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            selected = json.loads(
                (root / "out/resolution_selection.json").read_text()
            )
            self.assertEqual(selected["selected_resolution"], 0.4)
            self.assertTrue(
                selected["lower_complexity_used_only_as_tiebreaker"]
            )

    def test_cohort_resolution_rejects_low_grid_with_all_cells_mixed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            evidence = root / "grid.tsv"
            metrics = [
                "catalog_recall", "embedded_program_separation",
                "deg_antideg_coherence", "pseudobulk_coherence",
                "spatial_morphology_coherence",
                "adjacent_resolution_stability",
                "observation_support_coherence",
                "directly_resolved_observation_fraction",
                "mean_identity_margin",
            ]
            fields = [
                "resolution", "selection_purpose", "complete_catalog_scanned",
                "zero_census_audited", *metrics,
                "technical_fragmentation", "state_overfragmentation",
                "mixed_observation_fraction", "unresolved_observation_fraction",
                "complexity",
            ]
            rows = []
            for resolution in (0.1, 0.2, 0.4):
                row = {metric: 0.80 for metric in metrics}
                row.update({
                    "resolution": resolution,
                    "selection_purpose": "cohort_identity_resolution",
                    "complete_catalog_scanned": "true",
                    "zero_census_audited": "true",
                    "technical_fragmentation": 0,
                    "state_overfragmentation": 0,
                    "complexity": {0.1: 2, 0.2: 4, 0.4: 7}[resolution],
                })
                if resolution == 0.1:
                    row.update({
                        "directly_resolved_observation_fraction": 0,
                        "mean_identity_margin": 0.03,
                        "mixed_observation_fraction": 1,
                        "unresolved_observation_fraction": 0,
                    })
                else:
                    row.update({
                        "directly_resolved_observation_fraction": (
                            0.55 if resolution == 0.2 else 0.85
                        ),
                        "mean_identity_margin": 0.35,
                        "mixed_observation_fraction": (
                            0.35 if resolution == 0.2 else 0.10
                        ),
                        "unresolved_observation_fraction": 0.05,
                    })
                rows.append(row)
            with evidence.open("w", encoding="utf-8") as handle:
                handle.write("\t".join(fields) + "\n")
                for row in rows:
                    handle.write("\t".join(str(row[field]) for field in fields) + "\n")
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "select_lineage_resolution.py"),
                "--grid-evidence", str(evidence),
                "--selection-purpose", "cohort_identity_resolution",
                "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            selected = json.loads(
                (root / "out/resolution_selection.json").read_text()
            )
            self.assertEqual(selected["selected_resolution"], 0.4)

    def test_cohort_resolution_prefers_resolvable_mixture_over_unresolved_fragmentation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            evidence = root / "grid.tsv"
            metrics = [
                "catalog_recall", "embedded_program_separation",
                "deg_antideg_coherence", "pseudobulk_coherence",
                "spatial_morphology_coherence",
                "adjacent_resolution_stability",
                "observation_support_coherence",
                "directly_resolved_observation_fraction",
                "mean_identity_margin",
            ]
            fields = [
                "resolution", "selection_purpose", "complete_catalog_scanned",
                "zero_census_audited", *metrics,
                "technical_fragmentation", "state_overfragmentation",
                "mixed_observation_fraction", "unresolved_observation_fraction",
                "complexity",
            ]
            rows = []
            for resolution, mixed, unresolved in (
                (0.1, 0.20, 0.00),
                (0.2, 0.00, 0.20),
                (0.3, 0.05, 0.20),
            ):
                value = {metric: 0.80 for metric in metrics}
                value.update({
                    "resolution": resolution,
                    "selection_purpose": "cohort_identity_resolution",
                    "complete_catalog_scanned": "true",
                    "zero_census_audited": "true",
                    "technical_fragmentation": 0,
                    "state_overfragmentation": 0,
                    "mixed_observation_fraction": mixed,
                    "unresolved_observation_fraction": unresolved,
                    "complexity": int(resolution * 10),
                })
                rows.append(value)
            with evidence.open("w", encoding="utf-8") as handle:
                handle.write("\t".join(fields) + "\n")
                for row_value in rows:
                    handle.write(
                        "\t".join(str(row_value[field]) for field in fields)
                        + "\n"
                    )
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "select_lineage_resolution.py"),
                "--grid-evidence", str(evidence),
                "--selection-purpose", "cohort_identity_resolution",
                "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            selected = json.loads(
                (root / "out/resolution_selection.json").read_text()
            )
            self.assertEqual(selected["selected_resolution"], 0.1)

    def test_project_local_scorer_cannot_enter_formal_release(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            custom = root / "custom.R"
            custom.write_text("quit(status=0)\n", encoding="utf-8")
            dummy = root / "dummy"
            dummy.write_text("x\n", encoding="utf-8")
            scripts = {}
            for name in (
                "run_observation_lineage_scoring.R",
                "derive_candidate_local_subsets.R",
                "close_exact_remainders.py",
                "build_whole_tissue_cohort_plan.py",
                "adjudicate_second_round_subclusters.py",
                "merge_and_freeze_broad_membership.py",
                "route_global_atlas_v2.py",
                "validate_global_atlas_v2.py",
                "apply_post_merge_atlas_routing.py",
                "audit_post_merge_completeness.py",
                "audit_catalog_wide_lineage_challengers.py",
                "validate_catalog_wide_lineage_review_decisions.py",
                "apply_catalog_wide_lineage_review.py",
                "materialize_parent_locked_fine_proposals.py",
                "materialize_final_release_v2_2.py",
                "evaluate_annotation_robustness.py",
                "run_lineage_controller.py",
            ):
                path = custom if name == "run_observation_lineage_scoring.R" else SCRIPTS / name
                scripts[name] = {
                    "path": str(path),
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    "artifact_role": "runtime_input",
                }
            dependency = SCRIPTS / "lineage_controller_lib.py"
            artifact = lambda path: {
                "path": str(path),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "artifact_role": "runtime_input",
            }
            contract = root / "contract.json"
            contract.write_text(json.dumps({
                "schema_version": "2.0",
                "observation_unit": "cellbin",
                "selected_input_snapshot": {
                    "path": str(dummy), "sha256": "0" * 64,
                    "artifact_role": "runtime_input",
                },
                "artifact_role_policy": {
                    "runtime_allowed_roles": ["runtime_input", "external_reference"],
                    "failed_diagnostic_runtime_forbidden": True,
                    "diagnostic_registry": None,
                },
                "workflow_profile": artifact(dummy),
                "biological_profile": artifact(dummy),
                "candidate_catalog": artifact(dummy),
                "release_taxonomy": {
                    "independent_vascular_lineages": [
                        "Endothelial", "Pericyte/mural", "Smooth muscle",
                    ],
                    "lymphatic_parent": "Endothelial",
                    "legacy_vascular_associated_release_forbidden": True,
                    "single_public_annotation_column": "final_cell_type",
                },
                "canonical_lineage_controller": {
                    "controller_version": "2.2.0",
                    "phase_order": [
                        "whole_tissue_partition", "cluster_cohort_recluster",
                        "local_mixed_subcluster_split", "merge_and_freeze_broad",
                        "atlas_and_completeness_review", "materialize_final_release",
                    ],
                    "scripts": scripts,
                    "dependencies": {
                        "lineage_controller_lib.py": artifact(dependency),
                    },
                },
            }), encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "run_lineage_controller.py"),
                "whole_tissue_partition", "--contract", str(contract),
                "--rds", str(dummy),
                "--partitions", str(dummy), "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertTrue(
                "installed canonical" in result.stderr
                or "stale runtime dependency registry" in result.stderr,
                result.stderr,
            )

    def test_contract_builder_binds_identity_core_and_parent_lock_policy(self) -> None:
        source = (SCRIPTS / "build_annotation_contract_v2.py").read_text(
            encoding="utf-8"
        )
        for token in (
            '"identity_core_component_required": True',
            '"generic_support_cell_transitive_bridging": False',
            '"fine_candidate_parent_lock": True',
            '"fine_candidate_discriminator_seed_required": True',
            '"validated_high_purity_expression_subcluster_only"',
            '"sct_input_boundary": "project_local_non_sct_raw_counts"',
            '"clustering_path": "SCT_PCA_SNN_Leiden"',
            '"--contract-id"',
        ):
            self.assertIn(token, source)

    def test_sct_cohort_reclustering_restarts_from_raw_count_assay(self) -> None:
        source = (
            SCRIPTS / "run_seurat_cohort_recluster_impl.R"
        ).read_text(encoding="utf-8")
        self.assertIn('value != "SCT" && "counts" %in% Layers', source)
        self.assertIn(
            "refusing to run SCTransform on the SCT assay", source
        )
        self.assertIn('preferred <- c("RNA", "Spatial")', source)

    def test_two_runs_with_same_input_contract_and_seed_have_identical_membership_hash(self) -> None:
        rows = [
            {
                "cell_id": "b", "source_boundary": "x", "source_cluster": "1",
                "candidate_id": "smooth_muscle", "final_state": "defined_broad_only",
                "final_broad_label": "Smooth muscle", "final_fine_label": "",
                "confidence": "high", "assignment_origin": "supported_subset",
                "qc_reason": "",
            },
            {
                "cell_id": "a", "source_boundary": "x", "source_cluster": "1",
                "candidate_id": "stromal_mesenchymal", "final_state": "defined_broad_only",
                "final_broad_label": "Stromal/mesenchymal", "final_fine_label": "",
                "confidence": "moderate", "assignment_origin": "exact_remainder_parent",
                "qc_reason": "",
            },
        ]
        def run_once(seed: int) -> tuple[list[dict[str, str]], str]:
            self.assertEqual(seed, 2200)
            membership = sorted(rows, key=lambda item: item["cell_id"])
            return membership, deterministic_membership_hash(membership)

        first_membership, first_hash = run_once(2200)
        second_membership, second_hash = run_once(2200)
        self.assertEqual(first_membership, second_membership)
        self.assertEqual(first_hash, second_hash)


if __name__ == "__main__":
    unittest.main()
