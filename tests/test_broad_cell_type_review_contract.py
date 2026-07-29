from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "annotate-spatial-transcriptomics"
SCRIPTS = PACKAGE / "scripts"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


class BroadCellTypeReviewContractTests(unittest.TestCase):
    def test_marker_manifest_includes_complete_broad_families_not_subtypes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "build_cell_type_review_marker_manifest.py"),
                "--profile", str(PACKAGE / "references/profiles/sheep_ovary.json"),
                "--catalog", str(PACKAGE / "references/profiles/sheep_ovary_candidate_lineage_catalog.json"),
                "--out", str(out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            with (out / "cell_type_review_marker_manifest.tsv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            genes = {row["gene"] for row in rows}
            self.assertTrue({"KRT7", "MUC16", "MYH11", "CYP17A1", "ZP3"} <= genes)
            self.assertFalse(any("subtype" in row["family_id"].lower() for row in rows))
            self.assertFalse(any("state" in row["family_id"].lower() for row in rows))

    def test_present_broad_cannot_close_without_precision_recall_and_spatial_conclusions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue = root / "queue.tsv"
            write_tsv(queue, [{
                "review_id": "r1", "review_mode": "broad_lineage_review",
                "target_broad_label": "Granulosa", "unit_signature": "x",
            }])
            review = root / "review.json"
            review.write_text(json.dumps({
                "stage": "post_atlas_catalog_wide_lineage_review",
                "status": "ITERATION_REQUIRED",
                "catalog_wide_double_sided_review": True,
                "review_round": 1,
                "artifacts": {"review_queue": {"path": str(queue), "sha256": sha(queue)}},
            }), encoding="utf-8")
            packet_index = root / "packet_index.tsv"
            packet_sha = "a" * 64
            write_tsv(packet_index, [{
                "review_id": "r1", "review_mode": "broad_lineage_review",
                "target_broad_label": "Granulosa", "unit_signature": "x",
                "evidence_packet_sha256": packet_sha,
                "current_n": 10, "current_competitor_question_n": 0,
                "cross_type_over_recall_question_n": 0,
                "outside_recall_question_n": 0,
                "precision_evaluable": "true", "recall_evaluable": "true",
                "molecular_evaluable": "true", "spatial_evaluable": "true",
            }])
            packet_manifest = root / "packet_manifest.json"
            packet_manifest.write_text(json.dumps({
                "status": "PASS",
                "artifact_role": "broad_cell_type_review_evidence_packet_index",
                "review_manifest": {"path": str(review.resolve()), "sha256": sha(review)},
                "packet_index": {"path": str(packet_index.resolve()), "sha256": sha(packet_index)},
            }), encoding="utf-8")
            decisions = root / "decisions.tsv"
            write_tsv(decisions, [{
                "review_id": "r1", "review_mode": "broad_lineage_review",
                "outcome": "retain_current_cell_type", "proposed_broad_label": "",
                "evidence_packet_sha256": packet_sha,
                "evidence_basis": "query marker DEG pseudobulk and spatial review",
                "rationale": "The current cell type is retained without a complete conclusion record.",
                "membership_path": "", "membership_sha256": "",
            }])
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "validate_catalog_wide_lineage_review_decisions.py"),
                "--review-manifest", str(review), "--decisions", str(decisions),
                "--evidence-packet-manifest", str(packet_manifest),
                "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            validation = json.loads((root / "out/catalog_wide_lineage_decision_validation.json").read_text())
            self.assertEqual(validation["status"], "BLOCKED")
            self.assertTrue(any("current_member_precision" in error for error in validation["errors"]))

    def test_keywords_and_disposition_text_cannot_close_bound_challengers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue = root / "queue.tsv"
            write_tsv(queue, [{
                "review_id": "r1", "review_mode": "broad_lineage_review",
                "target_broad_label": "Luteal", "unit_signature": "u1",
            }])
            review = root / "review.json"
            review.write_text(json.dumps({
                "stage": "post_atlas_catalog_wide_lineage_review",
                "status": "ITERATION_REQUIRED",
                "catalog_wide_double_sided_review": True,
                "review_round": 1,
                "artifacts": {"review_queue": {
                    "path": str(queue), "sha256": sha(queue),
                }},
            }), encoding="utf-8")
            packet_index = root / "packet.tsv"
            packet_sha = "b" * 64
            write_tsv(packet_index, [{
                "review_id": "r1", "review_mode": "broad_lineage_review",
                "target_broad_label": "Luteal", "unit_signature": "u1",
                "evidence_packet_sha256": packet_sha, "current_n": 100,
                "current_competitor_question_n": 4,
                "cross_type_over_recall_question_n": 0,
                "outside_recall_question_n": 20,
                "precision_evaluable": "true", "recall_evaluable": "true",
                "molecular_evaluable": "true", "spatial_evaluable": "true",
                "ovary_spatial_status": "PASS",
                "oocyte_review_status": "PASS",
                "follicle_histology_status": "PASS",
            }])
            packet_manifest = root / "packet.json"
            packet_manifest.write_text(json.dumps({
                "status": "PASS",
                "artifact_role": "broad_cell_type_review_evidence_packet_index",
                "review_manifest": {"path": str(review.resolve()), "sha256": sha(review)},
                "packet_index": {"path": str(packet_index.resolve()), "sha256": sha(packet_index)},
            }), encoding="utf-8")
            decisions = root / "decisions.tsv"
            write_tsv(decisions, [{
                "review_id": "r1", "review_mode": "broad_lineage_review",
                "outcome": "retain_current_cell_type",
                "evidence_packet_sha256": packet_sha,
                "current_member_precision": "supported",
                "whole_query_recall": "complete",
                "molecular_support": "supported",
                "spatial_consistency": "consistent",
                "precision_challenger_disposition": "refuted_by_bound_evidence",
                "recall_challenger_disposition": "refuted_by_bound_evidence",
                "rationale": "All required query marker DEG pseudobulk and spatial words are present here.",
                "proposed_broad_label": "", "membership_path": "",
                "membership_sha256": "",
            }])
            result = subprocess.run([
                sys.executable,
                str(SCRIPTS / "validate_catalog_wide_lineage_review_decisions.py"),
                "--review-manifest", str(review),
                "--evidence-packet-manifest", str(packet_manifest),
                "--decisions", str(decisions), "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            validation = json.loads(
                (root / "out/catalog_wide_lineage_decision_validation.json").read_text()
            )
            self.assertTrue(any(
                "exact patch or targeted review" in error
                for error in validation["errors"]
            ))

    def test_patch_without_canonical_targeted_review_manifest_is_blocked(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            membership = root / "membership.tsv"
            write_tsv(membership, [{
                "cell_id": "c1", "final_broad_label": "",
            }])
            patch = root / "patch.tsv"
            write_tsv(patch, [{
                "cell_id": "c1", "new_broad_label": "Granulosa",
                "candidate_id": "granulosa",
            }])
            queue = root / "queue.tsv"
            write_tsv(queue, [{
                "review_id": "r1", "review_mode": "missing_broad_review",
                "target_broad_label": "Granulosa", "unit_signature": "u1",
            }])
            review = root / "review.json"
            review.write_text(json.dumps({
                "stage": "post_atlas_catalog_wide_lineage_review",
                "status": "ITERATION_REQUIRED",
                "catalog_wide_double_sided_review": True,
                "review_round": 1,
                "membership": {"path": str(membership), "sha256": sha(membership)},
                "artifacts": {"review_queue": {
                    "path": str(queue), "sha256": sha(queue),
                }},
            }))
            packet_index = root / "packet.tsv"
            packet_sha = "c" * 64
            write_tsv(packet_index, [{
                "review_id": "r1", "review_mode": "missing_broad_review",
                "target_broad_label": "Granulosa", "unit_signature": "u1",
                "evidence_packet_sha256": packet_sha, "current_n": 0,
                "current_precision_question_n": 0,
                "outside_recall_question_n": 1,
                "precision_evaluable": "false", "recall_evaluable": "true",
                "molecular_evaluable": "true", "spatial_evaluable": "true",
            }])
            packet_manifest = root / "packet.json"
            packet_manifest.write_text(json.dumps({
                "status": "PASS",
                "artifact_role": "broad_cell_type_review_evidence_packet_index",
                "review_manifest": {"path": str(review), "sha256": sha(review)},
                "packet_index": {
                    "path": str(packet_index), "sha256": sha(packet_index),
                },
            }))
            decisions = root / "decisions.tsv"
            write_tsv(decisions, [{
                "review_id": "r1", "review_mode": "missing_broad_review",
                "outcome": "apply_cell_type_membership_patch",
                "evidence_packet_sha256": packet_sha,
                "current_member_precision": "not_applicable",
                "whole_query_recall": "under_recall_detected",
                "molecular_support": "supported",
                "spatial_consistency": "localized_issue",
                "rationale": "A bounded query-derived component supports this exact recall patch.",
                "proposed_broad_label": "",
                "membership_path": str(patch),
                "membership_sha256": sha(patch),
                "targeted_review_manifest_path": "",
                "targeted_review_manifest_sha256": "",
            }])
            result = subprocess.run([
                sys.executable,
                str(SCRIPTS / "validate_catalog_wide_lineage_review_decisions.py"),
                "--review-manifest", str(review),
                "--evidence-packet-manifest", str(packet_manifest),
                "--decisions", str(decisions), "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            validation = json.loads(
                (root / "out/catalog_wide_lineage_decision_validation.json").read_text()
            )
            self.assertTrue(any(
                "targeted review manifest" in error
                for error in validation["errors"]
            ))


if __name__ == "__main__":
    unittest.main()
