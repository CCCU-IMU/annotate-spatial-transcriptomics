from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "annotate-spatial-transcriptomics"
SCRIPTS = SKILL / "scripts"


class SheepOvaryReleaseContract(unittest.TestCase):
    def test_profile_has_rfirst_fixed_parameters_and_atlas_policy(self) -> None:
        profile = json.loads((SKILL / "references/profiles/sheep_ovary_rfirst_profile.json").read_text())
        self.assertEqual(profile["workflow"]["preferred_backbone"], "seurat_r_first")
        fixed = profile["stereopy_cellbin_pped_contract"]
        self.assertEqual(fixed["entry_qc"]["minimum_counts"], 100)
        self.assertEqual(fixed["entry_qc"]["minimum_features"], 75)
        self.assertEqual(fixed["sctransform"]["method"], "glmGamPoi")
        self.assertEqual(fixed["sctransform"]["variable_features"], 3000)
        self.assertEqual(fixed["pca"], {"computed_components": 50, "neighbor_components": 30})
        self.assertEqual(fixed["neighbors"], {"k": 30, "method": "annoy", "trees": 50, "metric": "cosine"})
        self.assertEqual(fixed["clustering"]["candidate_resolutions"], [0.1, 0.2, 0.3, 0.4, 0.6])
        self.assertEqual(profile["external_reference_policy"]["primary_public_atlas"], "GSE233801")
        self.assertEqual(profile["external_reference_policy"]["calibration_origin"], "query_like_heldout_current_query_anchors")
        self.assertEqual(profile["external_reference_policy"]["accepted_calibrated_tiers"], ["high", "moderate"])
        self.assertEqual(profile["external_reference_policy"]["atlas_rescue_ceiling"], "broad_only")
        self.assertTrue(profile["release_policy"]["all_broad_rescues_enter_final_inclusive_deg_and_dotplots"])

    def test_resolver_automatically_selects_sheep_ovary_cellbin_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            context = root / "context.json"
            inspect = root / "inspect.json"
            output = root / "resolved.json"
            context.write_text(json.dumps({"species": "Ovis aries", "tissue": "ovary"}))
            inspect.write_text(json.dumps({"type": "Seurat", "n_features": 24000, "assays": ["Spatial"]}))
            subprocess.run(
                [sys.executable, str(SCRIPTS / "resolve_workflow_profile.py"), "--context", str(context), "--r-inspection", str(inspect), "--input-path", "/example/cellbin_PPed/sample.seurat.rds", "--out", str(output)],
                check=True,
                capture_output=True,
                text=True,
            )
            result = json.loads(output.read_text())
            self.assertEqual(result["preferred_backbone"], "seurat_r_first")
            self.assertTrue(result["fixed_cellbin_preprocessing_required"])
            self.assertEqual(result["primary_public_atlas"], "GSE233801")
            self.assertFalse(result["dotplot_only_reference_allows_cell_transfer"])

    def test_resolver_accepts_single_assay_scalar_from_r_jsonlite(self) -> None:
        """jsonlite auto_unbox emits a scalar when a Seurat object has one assay."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            context = root / "context.json"
            inspect = root / "inspect.json"
            output = root / "resolved.json"
            context.write_text(json.dumps({"species": "绵羊", "tissue": "卵巢"}))
            inspect.write_text(json.dumps({"type": "Seurat", "n_features": 31663, "assays": "Spatial"}))
            subprocess.run(
                [sys.executable, str(SCRIPTS / "resolve_workflow_profile.py"), "--context", str(context), "--r-inspection", str(inspect), "--input-path", "/example/cellbin_PPed/sample.seurat.rds", "--out", str(output)],
                check=True,
                capture_output=True,
                text=True,
            )
            result = json.loads(output.read_text())
            self.assertTrue(result["sheep_ovary_context"])
            self.assertTrue(result["fixed_cellbin_preprocessing_required"])

    def test_cellbin_runner_rejects_unrecorded_drift_and_missing_hashes(self) -> None:
        text = (SCRIPTS / "run_seurat_sct_preprocess.R").read_text()
        self.assertIn("allow-batch-exception", text)
        self.assertIn("batch-exception-reason", text)
        self.assertIn("digest is required for fail-closed", text)
        for token in ["v2", "glmGamPoi", "3000", "50000", "pca-npcs", "annoy-trees", "0.1,0.2,0.3,0.4,0.6"]:
            self.assertIn(token, text)

    def test_reference_self_calibration_and_legacy_writeback_are_blocked(self) -> None:
        bundle = (SCRIPTS / "build_depth_matched_atlas_bundle.py").read_text()
        tiered = (SCRIPTS / "calibrate_tiered_mapping_thresholds.py").read_text()
        legacy = (SCRIPTS / "calibrate_mapping_thresholds.py").read_text()
        materialize = (SCRIPTS / "materialize_adjudication_partition.py").read_text()
        terminal = (SCRIPTS / "commit_terminal_partitions.py").read_text()
        self.assertIn('"heldout_origin": "reference_self_classification"', bundle)
        self.assertIn('"eligible_for_query_rescue_calibration": False', bundle)
        self.assertIn("calibration-origin-manifest", tiered)
        self.assertIn("query_like_heldout_current_query_anchors", tiered)
        self.assertIn("anchor_ids_sha256", tiered)
        self.assertIn("target_ids_sha256", tiered)
        self.assertIn("Legacy one-tier calibration is disabled for writeback", legacy)
        self.assertIn("legacy or combined mapping-tier adjudication is not writeback eligible", materialize)
        self.assertIn("atlas/mapping broad returns lack observation-level tier proof", terminal)
        self.assertIn("RCTD-assisted returns lack observation-level confidence tier", terminal)

    def test_new_literature_enforces_shallow_cross_species_boundaries(self) -> None:
        text = (SKILL / "references/profiles/sheep_ovary_literature_2025_2026.md").read_text()
        for doi in ["10.1126/science.adx0659", "10.1002/advs.202517633", "10.1016/j.ajog.2024.05.046"]:
            self.assertIn(doi, text)
        self.assertIn("DLG2", text)
        self.assertIn("single-marker", text)
        self.assertIn("not a quota", text)
        case = (SKILL / "references/profiles/sheep_ovary_rfirst_case_reference.md").read_text()
        self.assertIn("reusable strategy trace", case)
        self.assertIn("inclusive broad membership", case)
        self.assertNotIn("D055", case)
        self.assertNotIn("/" + "share" + "/" + "org" + "/", case)


class CohortOrchestrationContract(unittest.TestCase):
    def test_one_sample_one_root_wave_assignment_and_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = root / "samples.tsv"
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["sample_id", "input_root"], delimiter="\t")
                writer.writeheader()
                for sample in ["S1", "S2", "S3"]:
                    writer.writerow({"sample_id": sample, "input_root": f"/inputs/{sample}"})
            cohort = root / "cohort"
            subprocess.run([sys.executable, str(SCRIPTS / "init_annotation_cohort.py"), "--manifest", str(manifest), "--cohort-root", str(cohort), "--max-active-workers", "2", "--cohort-id", "fixture"], check=True, capture_output=True, text=True)
            with (cohort / "control/worker_registry.tsv").open(newline="", encoding="utf-8") as handle:
                workers = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(workers), 3)
            self.assertEqual([row["wave"] for row in workers], ["1", "1", "2"])
            self.assertEqual(len({row["sample_root"] for row in workers}), 3)
            subprocess.run([sys.executable, str(SCRIPTS / "validate_cohort_state.py"), str(cohort)], check=True, capture_output=True, text=True)

    def test_active_worker_double_claim_requires_takeover(self) -> None:
        text = (SCRIPTS / "update_cohort_worker.py").read_text()
        self.assertIn("sample already has an active worker", text)
        self.assertIn("--takeover", text)
        orchestration = (SKILL / "references/multi-sample-agent-orchestration.md").read_text()
        self.assertIn("exactly one logical worker Agent to each sample", orchestration)
        self.assertIn("Parallelism never reduces evidence gates", orchestration)

    def test_sample_run_registry_blocks_duplicate_active_work(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            project = root / "project"
            subprocess.run([sys.executable, str(SCRIPTS / "init_annotation_project.py"), "--sample", "S1", "--input-root", str(root), "--project-root", str(project), "--modality", "spatial", "--observation-unit", "cellbin"], check=True, capture_output=True, text=True)
            common = [sys.executable, str(SCRIPTS / "update_run_registry.py"), str(project), "--work-key", "whole_sct", "--execution-fingerprint", "fixture-sha", "--stage", "preprocess", "--script", "runner.R", "--environment", "R", "--status", "prepared", "--output-root", str(project / "runs")]
            subprocess.run(common + ["--run-id", "run1"], check=True, capture_output=True, text=True)
            failed = subprocess.run(common + ["--run-id", "run2"], capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("active run already owns this work_key", failed.stderr + failed.stdout)

    def test_cohort_confirmation_is_hash_bound(self) -> None:
        request = (SCRIPTS / "request_cohort_confirmation.py").read_text()
        record = (SCRIPTS / "record_cohort_confirmation.py").read_text()
        validate = (SCRIPTS / "validate_cohort_state.py").read_text()
        self.assertIn("every sample must be SAMPLE_FROZEN", request)
        self.assertIn("cross-sample audit is absent or not PASS", request)
        self.assertIn("request_sha256", record)
        self.assertIn("cohort confirmation was invalidated", validate)


if __name__ == "__main__":
    unittest.main()
