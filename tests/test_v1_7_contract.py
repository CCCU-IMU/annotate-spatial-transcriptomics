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
SKILL = ROOT / "annotate-spatial-transcriptomics"
SCRIPTS = SKILL / "scripts"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, rows: list[dict], fields=None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


class V17ContractTests(unittest.TestCase):
    def test_new_project_enables_evidence_freeze_and_global_atlas(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "project"
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "init_annotation_project.py"),
                "--sample", "s1", "--input-root", temp, "--project-root", str(root),
                "--modality", "single-cell",
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            project = json.loads((root / "config/project.json").read_text())
            self.assertEqual(project["framework_version"], "1.8.0")
            self.assertEqual(project["global_atlas_mapping_scope"], "complete_analysis_set")
            cluster_header = (root / "state/cluster_decision_ledger.tsv").read_text().splitlines()[0]
            route_header = (root / "state/route_attempt_registry.tsv").read_text().splitlines()[0]
            self.assertIn("prelabel_evidence_artifact", cluster_header)
            self.assertIn("concordance_artifact", route_header)

    def test_prelabel_evidence_rejects_single_marker_anchoring(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "state").mkdir(); (root / "evidence").mkdir()
            positive = write_tsv(root / "evidence/positive.tsv", [{"gene": "A", "score": 2}, {"gene": "B", "score": 1}])
            anti = write_tsv(root / "evidence/anti.tsv", [{"gene": "X", "score": -1}])
            catalog = root / "evidence/candidates.json"
            catalog.write_text(json.dumps({"candidate_boundaries": [
                {"candidate_id": "granulosa", "review_required": True},
                {"candidate_id": "stromal_mesenchymal", "review_required": True},
            ]}))
            def artifact(path: Path) -> dict:
                return {"path": str(path), "sha256": sha(path)}
            document = {
                "schema_version": "1.0", "freeze_id": "freeze_1", "source_cluster": "0", "n_observations": 10,
                "label_blind": True, "favored_marker_interpretation_loaded": False,
                "candidate_catalog_artifact": artifact(catalog),
                "candidate_lineages": ["granulosa", "stromal_mesenchymal"],
                "lineage_hypotheses": [
                    {"candidate_id": "granulosa", "candidate_broad_lineage": "Granulosa", "eligible": True, "program_score": 2.0,
                     "positive_marker_family_count": 2, "positive_marker_families": ["AMH-HSD17B1", "SERPINE2-GSTA1"],
                     "anti_program_burden": 0.0, "contradiction_count": 0, "evidence_artifacts": [artifact(positive)]},
                    {"candidate_id": "stromal_mesenchymal", "candidate_broad_lineage": "Stromal/mesenchymal", "eligible": True, "program_score": 0.5,
                     "positive_marker_family_count": 1, "positive_marker_families": ["DCN-LUM"],
                     "anti_program_burden": 1.0, "contradiction_count": 1, "evidence_artifacts": [artifact(anti)]},
                ],
                "winner": "granulosa", "runner_up": "stromal_mesenchymal", "winning_margin": 1.5,
                "winner_confidence": "high", "positive_deg_artifact": artifact(positive), "anti_deg_artifact": artifact(anti),
                "technical_flags": {"low_depth": False, "doublet": False}, "frozen_before_initial_label": True,
            }
            evidence = root / "evidence/freeze.json"; evidence.write_text(json.dumps(document))
            ledger = write_tsv(root / "state/cluster_decision_ledger.tsv", [{
                "decision_id": "d1", "source_cluster": "0", "n_observations": 10,
                "prelabel_evidence_artifact": str(evidence), "prelabel_evidence_sha256": sha(evidence),
                "prelabel_evidence_frozen": "true", "prelabel_winner": "granulosa",
                "prelabel_runner_up": "stromal_mesenchymal", "prelabel_winning_margin": "1.5",
                "initial_broad_label": "Granulosa", "confidence": "high", "supersedes": "",
            }])
            passed = subprocess.run([sys.executable, str(SCRIPTS / "validate_prelabel_broad_evidence.py"), str(root)], capture_output=True, text=True)
            self.assertEqual(passed.returncode, 0, passed.stdout + passed.stderr)
            document["lineage_hypotheses"][0]["positive_marker_family_count"] = 1
            document["lineage_hypotheses"][0]["positive_marker_families"] = ["AMH"]
            evidence.write_text(json.dumps(document))
            rows = list(csv.DictReader(ledger.read_text().splitlines(), delimiter="\t")); rows[0]["prelabel_evidence_sha256"] = sha(evidence)
            write_tsv(ledger, rows)
            failed = subprocess.run([sys.executable, str(SCRIPTS / "validate_prelabel_broad_evidence.py"), str(root)], capture_output=True, text=True)
            self.assertEqual(failed.returncode, 2)
            self.assertIn("fewer than two marker families", failed.stdout)

    def test_all_cell_atlas_routes_qc_and_requires_discrepancy_closure(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); audit = root / "audit"
            ledger_rows = [
                {"cell_id": "q1", "analysis_scope": "analysis_set", "final_state": "qc_holdout", "final_broad_label": "", "final_fine_label": "", "source_cluster": "q"},
                {"cell_id": "q2", "analysis_scope": "analysis_set", "final_state": "qc_holdout", "final_broad_label": "", "final_fine_label": "", "source_cluster": "q"},
                {"cell_id": "a1", "analysis_scope": "analysis_set", "final_state": "defined_broad_only", "final_broad_label": "A", "final_fine_label": "", "source_cluster": "0"},
                {"cell_id": "a2", "analysis_scope": "analysis_set", "final_state": "defined_broad_only", "final_broad_label": "A", "final_fine_label": "", "source_cluster": "0"},
                {"cell_id": "a3", "analysis_scope": "analysis_set", "final_state": "defined_broad_only", "final_broad_label": "A", "final_fine_label": "", "source_cluster": "0"},
                {"cell_id": "u1", "analysis_scope": "analysis_set", "final_state": "defined_broad_only", "final_broad_label": "A", "final_fine_label": "", "source_cluster": "1"},
            ]
            mapping_rows = [
                {"cell_id": "q1", "predicted_label": "B", "mapping_tier": "moderate_only", "consensus_pass": "true", "out_of_distribution": "false", "ontology_conflict": "false"},
                {"cell_id": "q2", "predicted_label": "B", "mapping_tier": "low_reject", "consensus_pass": "false", "out_of_distribution": "false", "ontology_conflict": "false"},
                {"cell_id": "a1", "predicted_label": "B", "mapping_tier": "high", "consensus_pass": "true", "out_of_distribution": "false", "ontology_conflict": "false"},
                {"cell_id": "a2", "predicted_label": "B", "mapping_tier": "high", "consensus_pass": "true", "out_of_distribution": "false", "ontology_conflict": "false"},
                {"cell_id": "a3", "predicted_label": "A", "mapping_tier": "high", "consensus_pass": "true", "out_of_distribution": "false", "ontology_conflict": "false"},
                {"cell_id": "u1", "predicted_label": "A", "mapping_tier": "low_reject", "consensus_pass": "false", "out_of_distribution": "true", "ontology_conflict": "false"},
            ]
            ledger = write_tsv(root / "ledger.tsv", ledger_rows); mapping = write_tsv(root / "mapping.tsv", mapping_rows)
            built = subprocess.run([
                sys.executable, str(SCRIPTS / "build_global_atlas_concordance.py"),
                "--cell-ledger", str(ledger), "--atlas-mapping", str(mapping), "--out", str(audit),
                "--min-discordant-n", "2", "--min-discordant-fraction", "0.5",
                "--ood-min-n", "1", "--ood-min-fraction", "1.0",
            ], capture_output=True, text=True)
            self.assertEqual(built.returncode, 2, built.stdout + built.stderr)
            manifest = json.loads((audit / "global_atlas_concordance_manifest.json").read_text())
            self.assertEqual(manifest["comparison_counts"]["qc_writeback_candidate"], 1)
            self.assertEqual(manifest["review_queue_n"], 2)
            open_review = subprocess.run([sys.executable, str(SCRIPTS / "validate_global_atlas_concordance.py"), "--audit-root", str(audit)], capture_output=True, text=True)
            self.assertEqual(open_review.returncode, 2)
            evidence_source = root / "orthogonal.tsv"; write_tsv(evidence_source, [{"hypothesis": "A", "positive_programs": 2}, {"hypothesis": "B", "positive_programs": 4}])
            queue = list(csv.DictReader((audit / "discrepancy_review_queue.tsv").read_text().splitlines(), delimiter="\t"))
            decision_rows = []
            for row in queue:
                outcome = "supersede_broad" if row["trigger_type"] == "material_broad_disagreement" else "unknown_candidate"
                evidence = root / f"{row['review_id']}.json"
                evidence.write_text(json.dumps({
                    "schema_version": "1.0", "review_id": row["review_id"], "review_scope": "complete_cluster",
                    "primary_hypothesis": row["primary_broad"], "atlas_hypothesis": row["atlas_broad"], "outcome": outcome,
                    "query_full_feature_evidence": {"status": "PASS"},
                    "positive_marker_comparison": {"winner": "B"}, "anti_marker_comparison": {"contradictions": 0},
                    "pseudobulk_comparison": {"preferred": "B"}, "spatial_or_sample_consistency": {"status": "compatible"},
                    "technical_alternatives": {"doublet": "rejected", "ambient": "rejected"}, "atlas_only": False,
                    "rationale": "Independent full-feature positive and anti-program review closed this complete cluster.",
                    "evidence_artifacts": [{"path": str(evidence_source), "sha256": sha(evidence_source)}],
                }))
                decision_rows.append({
                    "review_id": row["review_id"], "outcome": outcome,
                    "rationale": "Independent full-feature positive and anti-program review closed this complete cluster.",
                    "evidence_artifact": str(evidence), "evidence_sha256": sha(evidence), "atlas_only": "false",
                })
            decisions = write_tsv(root / "decisions.tsv", decision_rows)
            closed = subprocess.run([sys.executable, str(SCRIPTS / "validate_global_atlas_concordance.py"), "--audit-root", str(audit), "--decisions", str(decisions)], capture_output=True, text=True)
            self.assertEqual(closed.returncode, 0, closed.stdout + closed.stderr)
            rows = list(csv.DictReader(decisions.read_text().splitlines(), delimiter="\t")); rows[0]["atlas_only"] = "true"; write_tsv(decisions, rows)
            unsafe = subprocess.run([sys.executable, str(SCRIPTS / "validate_global_atlas_concordance.py"), "--audit-root", str(audit), "--decisions", str(decisions)], capture_output=True, text=True)
            self.assertEqual(unsafe.returncode, 2)
            self.assertIn("Atlas-only", unsafe.stdout)

    def test_skill_contract_declares_efficient_open_set_mapping(self) -> None:
        text = (SKILL / "SKILL.md").read_text(encoding="utf-8")
        self.assertIn("complete `analysis_set`", text)
        self.assertIn("out-of-distribution", text)
        self.assertIn("never form a full query-by-reference distance matrix", text)

    def test_material_ontology_conflict_enters_review_queue(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ledger = write_tsv(root / "ledger.tsv", [
                {"cell_id": "a1", "analysis_scope": "analysis_set", "final_state": "defined_broad_only", "final_broad_label": "A", "final_fine_label": "", "source_cluster": "0"},
                {"cell_id": "a2", "analysis_scope": "analysis_set", "final_state": "defined_broad_only", "final_broad_label": "A", "final_fine_label": "", "source_cluster": "0"},
            ])
            mapping = write_tsv(root / "mapping.tsv", [
                {"cell_id": "a1", "predicted_label": "A", "mapping_tier": "high", "consensus_pass": "true", "out_of_distribution": "false", "ontology_conflict": "true"},
                {"cell_id": "a2", "predicted_label": "A", "mapping_tier": "high", "consensus_pass": "true", "out_of_distribution": "false", "ontology_conflict": "true"},
            ])
            audit = root / "audit"
            built = subprocess.run([
                sys.executable, str(SCRIPTS / "build_global_atlas_concordance.py"),
                "--cell-ledger", str(ledger), "--atlas-mapping", str(mapping), "--out", str(audit),
                "--min-discordant-n", "2", "--min-discordant-fraction", "1.0",
            ], capture_output=True, text=True)
            self.assertEqual(built.returncode, 2, built.stdout + built.stderr)
            queue = list(csv.DictReader((audit / "discrepancy_review_queue.tsv").read_text().splitlines(), delimiter="\t"))
            self.assertEqual(len(queue), 1)
            self.assertEqual(queue[0]["trigger_type"], "material_ontology_conflict")

    def test_atlas_index_contract_forbids_costly_defaults(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            feature = write_tsv(root / "feature.tsv", [{"gene": "A", "mean": 0}])
            representation = write_tsv(root / "reference.tsv", [{"cell_id": "r1", "pc1": 0.1}])
            crosswalk = write_tsv(root / "crosswalk.tsv", [{"source": "x", "broad": "A"}])
            index = root / "index.bin"; index.write_bytes(b"ann-index")
            def artifact(path: Path) -> dict:
                return {"path": str(path), "sha256": sha(path)}
            manifest = root / "atlas_index.json"
            document = {
                "schema_version": "1.0", "reference_id": "ref", "mapping_engine": "fixed_projection_ann",
                "n_reference": 1, "n_dimensions": 30, "feature_transform": artifact(feature),
                "reference_representation": artifact(representation), "neighbor_index": artifact(index),
                "broad_crosswalk": artifact(crosswalk), "dense_pairwise_matrix": False,
                "query_reference_joint_retraining": False, "whole_object_rctd": False,
                "reusable_across_queries": True,
            }
            manifest.write_text(json.dumps(document))
            passed = subprocess.run([sys.executable, str(SCRIPTS / "validate_atlas_index_manifest.py"), "--project-root", str(root), "--manifest", str(manifest)], capture_output=True, text=True)
            self.assertEqual(passed.returncode, 0, passed.stdout + passed.stderr)
            document["dense_pairwise_matrix"] = True; manifest.write_text(json.dumps(document))
            failed = subprocess.run([sys.executable, str(SCRIPTS / "validate_atlas_index_manifest.py"), "--project-root", str(root), "--manifest", str(manifest)], capture_output=True, text=True)
            self.assertEqual(failed.returncode, 2)


if __name__ == "__main__":
    unittest.main()
