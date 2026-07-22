from __future__ import annotations

import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "annotate-spatial-transcriptomics"
SCRIPTS = SKILL / "scripts"
PROFILE = SKILL / "references/profiles/sheep_ovary.json"
sys.path.insert(0, str(SCRIPTS))

from validate_broad_class_completeness import validate as validate_completeness  # noqa: E402
from validate_cohort_outcome import validate as validate_cohort  # noqa: E402
from validate_prelabel_broad_evidence import validate as validate_prelabel  # noqa: E402
from lineage_decision_lib import absolute_supported_families  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ref(path: Path, n: int | None = None) -> dict:
    record = {"path": str(path), "sha256": sha(path)}
    if n is not None:
        record["n_observations"] = n
    return record


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    return path


def decision_rows(epithelial_score: float = 0.01, stromal_score: float = 0.8) -> list[dict]:
    base = {
        "cluster": "1", "anti_program_burden": "0", "contradiction_count": "0",
        "eligible": "true", "purity_status": "pass", "lineage_supported_fraction": "0.8",
        "strongest_competing_fraction": "0.1",
    }
    return [
        {**base, "candidate_id": "epithelial_mesothelial", "candidate_broad_lineage": "Epithelial/mesothelial",
         "program_score": epithelial_score, "positive_family_count": "2", "positive_families": "keratin;surface"},
        {**base, "candidate_id": "stromal_mesenchymal", "candidate_broad_lineage": "Stromal/mesenchymal",
         "program_score": stromal_score, "positive_family_count": "2", "positive_families": "ecm;ovarian_stroma"},
    ]


class V202ContractTests(unittest.TestCase):
    def test_absolute_family_threshold_blocks_sparse_epithelial_claim(self) -> None:
        rows = [
            {"cluster": "4", "candidate_id": "epithelial_mesothelial", "family_name": "keratin",
             "detected_marker_count": "3", "any_detection_fraction": "0.0046", "mean_normalized_expression": "0.0059"},
            {"cluster": "4", "candidate_id": "epithelial_mesothelial", "family_name": "surface",
             "detected_marker_count": "4", "any_detection_fraction": "0.0056", "mean_normalized_expression": "0.0056"},
        ]
        thresholds = json.loads(PROFILE.read_text())["broad_family_evidence_contract"]["absolute_family_support_thresholds"]
        self.assertEqual(
            absolute_supported_families(rows, cluster="4", candidate_id="epithelial_mesothelial", thresholds=thresholds),
            [],
        )

    def test_prelabel_winner_is_recomputed_from_scores(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); (root / "config").mkdir(); (root / "state").mkdir(); (root / "evidence").mkdir()
            (root / "config/project.json").write_text(json.dumps({"framework_version": "1.9.0"}))
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({"candidate_boundaries": [
                {"candidate_id": "epithelial_mesothelial", "review_required": True},
                {"candidate_id": "stromal_mesenchymal", "review_required": True},
            ]}))
            generic = root / "evidence/generic.json"; generic.write_text(json.dumps({"content": ["marker", "anti"]}))
            hypotheses = [
                {"candidate_id": "epithelial_mesothelial", "candidate_broad_lineage": "Epithelial/mesothelial",
                 "eligible": True, "program_score": 0.01, "positive_marker_family_count": 2,
                 "positive_marker_families": ["keratin", "surface"], "anti_program_burden": 0,
                 "contradiction_count": 0, "evidence_artifacts": [ref(generic)]},
                {"candidate_id": "stromal_mesenchymal", "candidate_broad_lineage": "Stromal/mesenchymal",
                 "eligible": True, "program_score": 0.8, "positive_marker_family_count": 2,
                 "positive_marker_families": ["ecm", "ovarian_stroma"], "anti_program_burden": 0,
                 "contradiction_count": 0, "evidence_artifacts": [ref(generic)]},
            ]
            freeze = root / "evidence/freeze.json"
            freeze.write_text(json.dumps({
                "schema_version": "1.0", "freeze_id": "f1", "source_cluster": "1", "n_observations": 2,
                "label_blind": True, "favored_marker_interpretation_loaded": False,
                "candidate_catalog_artifact": ref(catalog),
                "candidate_lineages": ["epithelial_mesothelial", "stromal_mesenchymal"],
                "lineage_hypotheses": hypotheses, "winner": "epithelial_mesothelial",
                "runner_up": "stromal_mesenchymal", "winning_margin": 0.1, "winner_confidence": "moderate",
                "positive_deg_artifact": ref(generic), "anti_deg_artifact": ref(generic),
                "technical_flags": {"low_information": False}, "frozen_before_initial_label": True,
            }))
            ledger = write_tsv(root / "state/cluster_decision_ledger.tsv", [
                "decision_id", "source_cluster", "n_observations", "prelabel_evidence_artifact",
                "prelabel_evidence_sha256", "prelabel_evidence_frozen", "prelabel_winner",
                "prelabel_runner_up", "prelabel_winning_margin", "initial_broad_label", "confidence", "supersedes",
            ], [{
                "decision_id": "d1", "source_cluster": "1", "n_observations": "2",
                "prelabel_evidence_artifact": str(freeze), "prelabel_evidence_sha256": sha(freeze),
                "prelabel_evidence_frozen": "true", "prelabel_winner": "epithelial_mesothelial",
                "prelabel_runner_up": "stromal_mesenchymal", "prelabel_winning_margin": "0.1",
                "initial_broad_label": "Epithelial/mesothelial", "confidence": "moderate", "supersedes": "",
            }])
            result = validate_prelabel(root, ledger)
            self.assertEqual(result["status"], "FAIL")
            self.assertTrue(any("expected stromal_mesenchymal" in error for error in result["errors"]), result["errors"])

    def make_v2_cohort(self, root: Path, target_candidate: str, target_label: str) -> Path:
        (root / "config").mkdir(); (root / "evidence").mkdir(); (root / "memberships").mkdir()
        (root / "config/project.json").write_text(json.dumps({"framework_version": "2.0.2"}))
        query = write_tsv(root / "memberships/query.tsv", ["cell_id"], [{"cell_id": "a"}, {"cell_id": "b"}])
        resolutions = []
        for resolution in (0.1, 0.2):
            membership = write_tsv(root / f"memberships/res{resolution}.tsv", ["cell_id", "cluster"], [
                {"cell_id": "a", "cluster": "1"}, {"cell_id": "b", "cluster": "1"},
            ])
            index = root / f"evidence/res{resolution}.json"; index.write_text(json.dumps({"programs": ["all_catalog"], "anti": ["reviewed"]}))
            resolutions.append({"resolution": resolution, "membership": ref(membership, 2), "cluster_count": 1, "evidence_index": ref(index)})
        full = root / "evidence/full.json"; full.write_text(json.dumps({"feature_count": 20000, "programs": ["all_catalog"]}))
        decisions = write_tsv(root / "evidence/lineage_decisions.tsv", list(decision_rows()[0]), decision_rows())
        family_fields = ["cluster", "candidate_id", "family_name", "detected_marker_count", "any_detection_fraction", "mean_normalized_expression"]
        family_rows = []
        for candidate_id, families in {
            "epithelial_mesothelial": ["keratin", "surface"],
            "stromal_mesenchymal": ["ecm", "ovarian_stroma"],
        }.items():
            for family in families:
                family_rows.append({"cluster": "1", "candidate_id": candidate_id, "family_name": family,
                                    "detected_marker_count": "2", "any_detection_fraction": "0.5",
                                    "mean_normalized_expression": "0.2"})
        family_table = write_tsv(root / "evidence/broad_family.tsv", family_fields, family_rows)
        family_manifest = root / "evidence/broad_family_manifest.json"
        family_manifest.write_text(json.dumps({
            "cluster_membership": resolutions[0]["membership"], "evidence_table": ref(family_table),
            "biological_profile": ref(PROFILE),
        }))
        family_validation = root / "evidence/broad_family_validation.json"
        family_validation.write_text(json.dumps({"status": "PASS", "manifest_sha256": sha(family_manifest), "evidence_rows": 4}))
        returned = write_tsv(root / "memberships/returned.tsv", ["cell_id"], [{"cell_id": "a"}, {"cell_id": "b"}])
        outcome = root / "evidence/outcome.json"
        outcome.write_text(json.dumps({
            "schema_version": "2.0", "cohort_id": "c1", "cohort_type": "broad_class_recluster",
            "question_mode": "broad_purity_audit", "query_membership": ref(query, 2),
            "candidate_grid": [0.1, 0.2], "resolutions": resolutions, "selected_resolution": 0.1,
            "lineage_decision_table": ref(decisions),
            "broad_family_evidence_manifest": ref(family_manifest),
            "broad_family_evidence_validation": ref(family_validation),
            "marker_evidence": {"positive_marker_families": ["ecm", "ovarian_stroma"],
                "anti_marker_results": {"strongest_competitor": "reviewed"},
                "full_feature_scope": {"status": "PASS", "feature_space": "full_feature", "artifact": ref(full)}},
            "adjacent_stability": [{"left": 0.1, "right": 0.2, "metric": "ARI", "value": 1.0}],
            "spatial_morphology": {"status": "reviewed", "pattern": "coherent"},
            "source_qc_composition": {"status": "reviewed", "n": 2},
            "selected_resolution_rationale": "The selected resolution preserves the dominant resident program and purity.",
            "rejected_alternatives": [{"resolution": 0.2, "reason": "No additional stable lineage separation was supported."}],
            "subcluster_outcomes": [{"subcluster_id": "1", "outcome": "parent_return",
                "membership": ref(returned, 2), "target_candidate_id": target_candidate,
                "target_broad_label": target_label, "target_fine_label": "", "confidence": "moderate",
                "return_scope": "whole_subcluster"}],
            "terminal_outcome": "homogeneous_parent_confirmed", "homogeneous_parent_confirmed": True,
        }))
        return outcome

    def test_cohort_rejects_whole_cluster_epithelial_override(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outcome = self.make_v2_cohort(root, "epithelial_mesothelial", "Epithelial/mesothelial")
            result = validate_cohort(root, outcome)
            self.assertEqual(result["status"], "FAIL")
            self.assertTrue(any("expected stromal_mesenchymal" in error for error in result["errors"]), result["errors"])

    def test_cohort_accepts_numerical_winner_with_purity(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outcome = self.make_v2_cohort(root, "stromal_mesenchymal", "Stromal/mesenchymal")
            result = validate_cohort(root, outcome)
            self.assertEqual(result["status"], "PASS", result["errors"])

    def test_present_completeness_fails_when_competitor_dominates(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); (root / "config").mkdir(); (root / "state").mkdir(); (root / "evidence").mkdir()
            (root / "config/project.json").write_text(json.dumps({"framework_version": "2.0.2"}))
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({"catalog_id": "c", "candidate_boundaries": [{
                "candidate_id": "epithelial_mesothelial", "review_required": True, "release_level": "default_broad_candidate",
            }]}))
            evidence = root / "evidence/epithelial.json"
            evidence.write_text(json.dumps({
                "schema_version": "2.0", "candidate_lineage": "epithelial_mesothelial", "final_n_observations": 100,
                "own_program": {"positive_family_count": 2, "family_detection_fractions": {"keratin": 0.005, "surface": 0.006}, "program_score": 0.01},
                "strongest_competing_program": {"candidate_lineage": "stromal_mesenchymal", "program_score": 0.8},
                "observation_level_purity": {"status": "PASS", "lineage_supported_fraction": 0.01, "strongest_competing_fraction": 0.8, "contradiction_fraction": 0.8},
                "spatial_morphology": {"status": "PASS", "rationale": "A broad internal diffuse component."},
            }))
            fields = ["review_id", "candidate_lineage", "final_n_observations", "census_status", "query_marker_program_review",
                      "deg_review", "spatial_morphology_review", "observation_level_review", "large_label_embedding_review",
                      "decision", "evidence_artifact", "evidence_artifact_sha256", "closure_rationale", "status", "supersedes"]
            write_tsv(root / "state/broad_class_completeness_registry.tsv", fields, [{
                "review_id": "r1", "candidate_lineage": "epithelial_mesothelial", "final_n_observations": "100",
                "census_status": "present", "query_marker_program_review": "reviewed", "deg_review": "reviewed",
                "spatial_morphology_review": "reviewed", "observation_level_review": "reviewed",
                "large_label_embedding_review": "reviewed", "decision": "supported", "evidence_artifact": str(evidence),
                "evidence_artifact_sha256": sha(evidence), "closure_rationale": "numerical audit", "status": "audited", "supersedes": "",
            }])
            result = validate_completeness(root, catalog)
            self.assertEqual(result["status"], "BLOCKED")
            self.assertTrue(any("competing program" in error for error in result["errors"]), result["errors"])

    def test_profile_routes_mural_contractile_cells_to_vascular_parent(self) -> None:
        profile = json.loads(PROFILE.read_text())
        smooth = profile["lineages"]["smooth_muscle"]
        self.assertTrue({"RGS5", "PDGFRB", "NOTCH3"} <= set(smooth["mural_exclusion_program"]))
        self.assertIn("Vascular-associated", smooth["safety"])
        epithelial = profile["lineages"]["epithelial_mesothelial"]
        self.assertIn("never permission to expand", epithelial["broad_recall_rule"])


if __name__ == "__main__":
    unittest.main()
