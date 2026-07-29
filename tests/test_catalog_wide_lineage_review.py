from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "annotate-spatial-transcriptomics"
SCRIPTS = PACKAGE / "scripts"
THRESHOLDS = PACKAGE / "references/controller_thresholds_v2_2.json"
HAS_RUNTIME = all(importlib.util.find_spec(name) for name in ("numpy", "pandas", "scipy"))
sys.path.insert(0, str(SCRIPTS))

from lineage_controller_lib import candidate_can_support_broad_review  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(rows[0])
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def artifact(path: Path) -> dict[str, str]:
    return {"path": str(path.resolve()), "sha256": sha(path)}


def write_evidence_packet_manifest(
    root: Path, review_manifest: Path, queue: list[dict[str, str]],
) -> tuple[Path, dict[str, str]]:
    packet_hashes: dict[str, str] = {}
    rows = []
    for row in queue:
        review_id = row["review_id"]
        packet_hashes[review_id] = hashlib.sha256(
            f"packet::{review_id}::{row['unit_signature']}".encode()
        ).hexdigest()
        present = row["review_mode"] == "broad_lineage_review"
        recall_question = 5 if row["target_broad_label"] == "Granulosa" else 0
        rows.append({
            "review_id": review_id, "review_mode": row["review_mode"],
            "target_broad_label": row["target_broad_label"],
            "unit_signature": row["unit_signature"],
            "evidence_packet_sha256": packet_hashes[review_id],
            "current_n": 5 if present else 0,
            "current_competitor_question_n": 0,
            "cross_type_over_recall_question_n": 0,
            "outside_recall_question_n": recall_question,
            "precision_evaluable": str(present).lower(),
            "recall_evaluable": "true", "molecular_evaluable": "true",
            "spatial_evaluable": "true",
        })
    packet_index = root / "evidence_packet_index.tsv"
    write_tsv(packet_index, rows)
    current_questions = root / "current_questions.tsv.gz"
    write_tsv(current_questions, [{
        "broad_label": "Stromal/mesenchymal", "cell_id": "c0",
    }])
    recall_questions = root / "recall_questions.tsv.gz"
    write_tsv(recall_questions, [{
        "broad_label": "Granulosa", "cell_id": f"c{index}",
    } for index in range(5)])
    broad_manifest = root / "broad_evidence_manifest.json"
    broad_manifest.write_text(json.dumps({
        "artifacts": {
            "current_member_questions": artifact(current_questions),
            "recall_membership": artifact(recall_questions),
        },
    }), encoding="utf-8")
    manifest = root / "evidence_packet_manifest.json"
    manifest.write_text(json.dumps({
        "status": "PASS",
        "artifact_role": "broad_cell_type_review_evidence_packet_index",
        "review_manifest": artifact(review_manifest),
        "broad_evidence_manifest": artifact(broad_manifest),
        "packet_index": artifact(packet_index),
    }), encoding="utf-8")
    return manifest, packet_hashes


def cluster_evidence(cluster: str, candidate: str, positive: bool) -> dict[str, object]:
    return {
        "resolution_role": "selected", "source_boundary": "cohort_1",
        "source_cluster": cluster, "candidate_id": candidate,
        "n_observations": 5, "available_positive_family_count": 2,
        "group_positive_family_supported_count": 2 if positive else 0,
        "group_required_positive_families_pass": str(positive).lower(),
        "observation_identity_core_fraction": 0.8 if positive else 0,
        "positive_marker_detection_fraction": 0.8 if positive else 0,
        "mean_program_score": 0.4 if positive else 0,
        "marker_deg_log2fc_mean": 1.5 if positive else 0,
        "anti_marker_deg_log2fc_mean": 0,
    }


@unittest.skipUnless(HAS_RUNTIME, "numpy/pandas/scipy runtime unavailable")
class CatalogWideLineageReviewFunctionalTests(unittest.TestCase):
    def test_group_level_program_creates_recall_watch_when_observation_scores_are_sparse(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            contract = root / "contract.json"
            contract.write_text("{}\n", encoding="utf-8")
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({"candidate_boundaries": [
                {
                    "candidate_id": "epithelial", "candidate_role": "broad",
                    "release_broad_label": "Epithelial/mesothelial",
                    "specificity_priority": 80,
                    "required_positive_families": ["epi_a", "epi_b"],
                },
                {
                    "candidate_id": "stromal_mesenchymal", "candidate_role": "broad",
                    "release_broad_label": "Stromal/mesenchymal",
                    "specificity_priority": 10,
                    "required_positive_families": ["str_a", "str_b"],
                    "writeback_strategy": "generic_exact_remainder_after_specific_lineages",
                },
            ]}), encoding="utf-8")
            membership = root / "membership.tsv"
            members = [{
                "cell_id": f"c{i}", "source_boundary": "cohort_1",
                "source_cluster": "mixed", "final_broad_label": "Stromal/mesenchymal",
                "final_state": "defined_broad_only", "candidate_id": "stromal_mesenchymal",
                "assignment_origin": "fixture", "final_fine_label": "",
            } for i in range(100)]
            write_tsv(membership, members)
            scores = []
            for i in range(100):
                for candidate in ("epithelial", "stromal_mesenchymal"):
                    stromal = candidate == "stromal_mesenchymal"
                    scores.append({
                        "cell_id": f"c{i}", "source_boundary": "cohort_1",
                        "source_cluster": "mixed", "candidate_id": candidate,
                        "normalized_evidence": 0.4 if stromal else 0.1,
                        "direct_signal": 0.2 if stromal else 0,
                        "program_score": 0.2 if stromal else 0,
                        "positive_family_count": 2 if stromal else 0,
                        "family_coherent": str(stromal).lower(),
                        "release_family_coherent": str(stromal).lower(),
                        "identity_core_coherent": str(stromal).lower(),
                        "identity_core_direct": str(stromal).lower(),
                        "hard_contradiction": "false", "technical_flag": "false",
                        "x": i, "y": 0,
                    })
            score_path = root / "scores.tsv.gz"
            write_tsv(score_path, scores)
            epithelial_program = cluster_evidence("mixed", "epithelial", True)
            epithelial_program["n_observations"] = 100
            epithelial_program["observation_identity_core_fraction"] = 0.05
            epithelial_program["observation_identity_core_direct_fraction"] = 0.02
            evidence_path = root / "evidence.tsv"
            write_tsv(evidence_path, [
                epithelial_program,
                cluster_evidence("mixed", "stromal_mesenchymal", True),
            ])
            authority = root / "authority.json"
            authority.write_text(json.dumps({
                "mode": "stage_authority", "phase": "atlas_and_completeness_review",
                "annotation_contract_sha256": sha(contract),
                "post_atlas_membership": artifact(membership),
                "candidate_catalog": artifact(catalog),
                "threshold_registry": artifact(THRESHOLDS),
                "observation_scores": [artifact(score_path)],
                "cluster_evidence": [artifact(evidence_path)],
            }), encoding="utf-8")
            out = root / "review"
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "audit_catalog_wide_lineage_challengers.py"),
                "--contract", str(contract), "--stage-authority", str(authority),
                "--membership", str(membership), "--catalog", str(catalog),
                "--threshold-registry", str(THRESHOLDS), "--scores", str(score_path),
                "--cluster-evidence", str(evidence_path), "--out", str(out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2, result.stderr or result.stdout)
            queue = read_tsv(out / "catalog_wide_lineage_review_queue.tsv")
            self.assertTrue(any(
                row["review_mode"] == "missing_broad_review"
                and row["target_broad_label"] == "Epithelial/mesothelial"
                for row in queue
            ))
            watches = read_tsv(out / "outside_label_group_watch.tsv")
            self.assertTrue(any(
                row["target_broad_label"] == "Epithelial/mesothelial"
                for row in watches
            ))
            components = read_tsv(out / "outside_label_recall_components.tsv")
            self.assertFalse(any(
                row["target_broad_label"] == "Epithelial/mesothelial"
                for row in components
            ))

    def test_recall_component_is_reviewed_applied_and_closed_next_round(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            contract = root / "contract.json"
            contract.write_text("{}\n", encoding="utf-8")
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({"candidate_boundaries": [
                {
                    "candidate_id": "granulosa", "candidate_role": "broad",
                    "release_broad_label": "Granulosa", "specificity_priority": 80,
                    "required_positive_families": ["gran_a", "gran_b"],
                },
                {
                    "candidate_id": "stromal_mesenchymal", "candidate_role": "broad",
                    "release_broad_label": "Stromal/mesenchymal", "specificity_priority": 10,
                    "required_positive_families": ["str_a", "str_b"],
                    "writeback_strategy": "generic_exact_remainder_after_specific_lineages",
                },
            ]}), encoding="utf-8")
            membership = root / "membership.tsv"
            members = []
            scores = []
            for index in range(10):
                cluster = "g" if index < 5 else "s"
                cell = f"c{index}"
                members.append({
                    "cell_id": cell, "source_boundary": "cohort_1",
                    "source_cluster": cluster,
                    "final_broad_label": "Stromal/mesenchymal",
                    "final_state": "defined_broad_only", "candidate_id": "stromal_mesenchymal",
                    "assignment_origin": "fixture", "final_fine_label": "",
                })
                for candidate in ("granulosa", "stromal_mesenchymal"):
                    positive = (
                        candidate == "granulosa" and cluster == "g"
                    ) or (
                        candidate == "stromal_mesenchymal" and cluster == "s"
                    )
                    scores.append({
                        "cell_id": cell, "source_boundary": "cohort_1",
                        "source_cluster": cluster, "candidate_id": candidate,
                        "normalized_evidence": 0.8 if positive else 0.05,
                        "direct_signal": 0.7 if positive else 0.02,
                        "program_score": 0.5 if positive else 0,
                        "positive_family_count": 2 if positive else 0,
                        "family_coherent": str(positive).lower(),
                        "release_family_coherent": str(positive).lower(),
                        "identity_core_coherent": str(positive).lower(),
                        "identity_core_direct": str(positive).lower(),
                        "hard_contradiction": "false", "technical_flag": "false",
                        "x": index if cluster == "g" else index + 20, "y": 0,
                    })
            write_tsv(membership, members)
            score_path = root / "scores.tsv.gz"
            write_tsv(score_path, scores)
            evidence_path = root / "evidence.tsv"
            write_tsv(evidence_path, [
                cluster_evidence("g", "granulosa", True),
                cluster_evidence("g", "stromal_mesenchymal", False),
                cluster_evidence("s", "granulosa", False),
                cluster_evidence("s", "stromal_mesenchymal", True),
            ])
            authority = root / "review_authority.json"
            authority.write_text(json.dumps({
                "mode": "stage_authority", "phase": "atlas_and_completeness_review",
                "annotation_contract_sha256": sha(contract),
                "post_atlas_membership": artifact(membership),
                "candidate_catalog": artifact(catalog),
                "threshold_registry": artifact(THRESHOLDS),
                "observation_scores": [artifact(score_path)],
                "cluster_evidence": [artifact(evidence_path)],
            }), encoding="utf-8")
            review_1 = root / "review_1"
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "audit_catalog_wide_lineage_challengers.py"),
                "--contract", str(contract), "--stage-authority", str(authority),
                "--membership", str(membership), "--catalog", str(catalog),
                "--threshold-registry", str(THRESHOLDS), "--scores", str(score_path),
                "--cluster-evidence", str(evidence_path), "--out", str(review_1),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2, result.stderr or result.stdout)
            queue = read_tsv(review_1 / "catalog_wide_lineage_review_queue.tsv")
            self.assertTrue(any(
                row["review_mode"] == "missing_broad_review"
                and row["target_broad_label"] == "Granulosa" for row in queue
            ))
            self.assertEqual(
                {row["target_broad_label"] for row in queue},
                {"Granulosa", "Stromal/mesenchymal"},
            )
            granulosa_patch = root / "granulosa_patch.tsv"
            write_tsv(granulosa_patch, [{
                "cell_id": f"c{index}", "new_broad_label": "Granulosa",
                "candidate_id": "granulosa",
            } for index in range(5)])
            decisions = root / "decisions.tsv"
            packet_manifest, packet_hashes = write_evidence_packet_manifest(
                root, review_1 / "catalog_wide_lineage_review_manifest.json", queue
            )
            granulosa_review_id = next(
                row["review_id"] for row in queue
                if row["target_broad_label"] == "Granulosa"
            )
            targeted_manifest = root / "granulosa_targeted_review.json"
            targeted_manifest.write_text(json.dumps({
                "status": "PASS",
                "stage": "per_broad_targeted_membership_evidence",
                "route_class": "canonical_per_broad_evidence_packet",
                "review_id": granulosa_review_id,
                "target_broad_label": "Granulosa",
                "evidence_packet_sha256": packet_hashes[granulosa_review_id],
                "source_membership": artifact(membership),
                "patch_membership": artifact(granulosa_patch),
            }), encoding="utf-8")
            decision_rows = []
            for row in queue:
                recall = row["target_broad_label"] == "Granulosa"
                decision_rows.append({
                    "review_id": row["review_id"], "review_mode": row["review_mode"],
                    "outcome": "apply_cell_type_membership_patch" if recall else "retain_current_cell_type",
                    "evidence_packet_sha256": packet_hashes[row["review_id"]],
                    "current_member_precision": "not_applicable" if recall else "supported",
                    "whole_query_recall": "under_recall_detected" if recall else "complete",
                    "spatial_consistency": "localized_issue" if recall else "consistent",
                    "molecular_support": "supported",
                    "proposed_broad_label": "",
                    "evidence_basis": "query raw-count marker DEG pseudobulk and spatial review",
                    "rationale": "Direct multigene identity and the bounded spatial component support this exact decision.",
                    "membership_path": str(granulosa_patch) if recall else "",
                    "membership_sha256": sha(granulosa_patch) if recall else "",
                    "targeted_review_manifest_path": (
                        str(targeted_manifest) if recall else ""
                    ),
                    "targeted_review_manifest_sha256": (
                        sha(targeted_manifest) if recall else ""
                    ),
                })
            write_tsv(decisions, decision_rows)
            validation_out = root / "validation"
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "validate_catalog_wide_lineage_review_decisions.py"),
                "--review-manifest", str(review_1 / "catalog_wide_lineage_review_manifest.json"),
                "--evidence-packet-manifest", str(packet_manifest),
                "--decisions", str(decisions), "--out", str(validation_out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            validation = validation_out / "catalog_wide_lineage_decision_validation.json"
            apply_authority = root / "apply_authority.json"
            apply_authority.write_text(json.dumps({
                "mode": "stage_authority", "phase": "atlas_and_completeness_review",
                "annotation_contract_sha256": sha(contract),
                "post_atlas_membership": artifact(membership),
                "catalog_review_manifest": artifact(review_1 / "catalog_wide_lineage_review_manifest.json"),
                "catalog_decision_validation": artifact(validation),
                "candidate_catalog": artifact(catalog),
            }), encoding="utf-8")
            apply_out = root / "apply"
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "apply_catalog_wide_lineage_review.py"),
                "--contract", str(contract), "--stage-authority", str(apply_authority),
                "--membership", str(membership),
                "--review-manifest", str(review_1 / "catalog_wide_lineage_review_manifest.json"),
                "--decision-validation", str(validation), "--catalog", str(catalog),
                "--out", str(apply_out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            revised = apply_out / "catalog_wide_reviewed_membership.tsv.gz"
            revised_rows = read_tsv(revised)
            self.assertEqual(
                {row["final_broad_label"] for row in revised_rows[:5]}, {"Granulosa"}
            )
            authority_2 = root / "review_authority_2.json"
            authority_2.write_text(json.dumps({
                "mode": "stage_authority", "phase": "atlas_and_completeness_review",
                "annotation_contract_sha256": sha(contract),
                "post_atlas_membership": artifact(revised),
                "candidate_catalog": artifact(catalog),
                "threshold_registry": artifact(THRESHOLDS),
                "observation_scores": [artifact(score_path)],
                "cluster_evidence": [artifact(evidence_path)],
            }), encoding="utf-8")
            review_2 = root / "review_2"
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "audit_catalog_wide_lineage_challengers.py"),
                "--contract", str(contract), "--stage-authority", str(authority_2),
                "--membership", str(revised), "--catalog", str(catalog),
                "--threshold-registry", str(THRESHOLDS), "--scores", str(score_path),
                "--cluster-evidence", str(evidence_path), "--round-index", "2",
                "--previous-review-manifest", str(review_1 / "catalog_wide_lineage_review_manifest.json"),
                "--prior-decision-validation", str(validation), "--out", str(review_2),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2, result.stderr or result.stdout)
            manifest = json.loads(
                (review_2 / "catalog_wide_lineage_review_manifest.json").read_text()
            )
            self.assertEqual(manifest["status"], "ITERATION_REQUIRED")
            self.assertEqual(manifest["review_queue_n"], 2)


class CatalogWideLineageReviewArchitectureTests(unittest.TestCase):
    def test_parent_locked_fine_program_cannot_reconstruct_broad(self) -> None:
        self.assertFalse(candidate_can_support_broad_review({
            "candidate_role": "fine", "release_broad_label": "Luteal",
            "parent_broad_label": "Luteal",
        }))
        self.assertTrue(candidate_can_support_broad_review({
            "candidate_role": "fine", "release_broad_label": "Endothelial",
            "parent_broad_label": "Endothelial",
            "parent_broad_reconstruction_allowed": True,
        }))

    def test_controller_requires_catalog_wide_review_before_release(self) -> None:
        controller = (SCRIPTS / "run_lineage_controller.py").read_text(encoding="utf-8")
        validator = (SCRIPTS / "validate_lineage_controller_release.py").read_text(encoding="utf-8")
        self.assertIn("run_catalog_wide_review_iterations", controller)
        self.assertIn("--lineage-review-decisions", controller)
        self.assertIn("catalog_wide_lineage_review_status", validator)
        self.assertIn("catalog_wide_double_sided_review", validator)
        audit = (SCRIPTS / "audit_catalog_wide_lineage_challengers.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("逐大类全样本复核", audit)

    def test_review_thresholds_are_centralized(self) -> None:
        thresholds = json.loads(THRESHOLDS.read_text(encoding="utf-8"))
        policy = thresholds["catalog_wide_lineage_review_policy"]
        self.assertEqual(policy["maximum_decision_rounds"], 2)
        self.assertEqual(policy["minimum_group_watch_identity_fraction"], 0.005)
        source = (SCRIPTS / "audit_catalog_wide_lineage_challengers.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("catalog_wide_lineage_review_policy", source)
        self.assertIn("whole_object_per_cell_classifier_used", source)


if __name__ == "__main__":
    unittest.main()
