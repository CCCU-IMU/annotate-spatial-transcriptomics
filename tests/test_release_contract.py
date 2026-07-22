from __future__ import annotations

import csv
import gzip
import importlib.util
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
        self.assertEqual(profile["workflow"]["preferred_backbone"], "seurat_r_first_context_adaptive_graph")
        fixed = profile["stereopy_cellbin_pped_contract"]
        self.assertEqual(fixed["entry_qc"]["minimum_counts"], 100)
        self.assertEqual(fixed["entry_qc"]["minimum_features"], 75)
        self.assertEqual(fixed["sctransform"]["method"], "glmGamPoi")
        self.assertEqual(fixed["sctransform"]["variable_features"], 4000)
        self.assertEqual(fixed["banksy"]["input"], "SCT scale.data Pearson residuals for 4000 SCT variable features")
        self.assertEqual(fixed["banksy"]["k_geom"], 30)
        self.assertEqual(fixed["banksy"]["lambda"], 0.2)
        self.assertEqual(fixed["banksy"]["k_neighbors"], 50)
        self.assertEqual(fixed["clustering"]["candidate_resolutions"], [0.2, 0.4, 0.6, 0.8])
        self.assertEqual(profile["external_reference_policy"]["primary_public_atlas"], "GSE233801")
        self.assertEqual(profile["external_reference_policy"]["calibration_origin"], "query_like_heldout_current_query_anchors")
        self.assertEqual(profile["external_reference_policy"]["accepted_calibrated_tiers"], ["high", "moderate_only"])
        self.assertEqual(profile["external_reference_policy"]["atlas_rescue_ceiling"], "broad_only")
        self.assertEqual(profile["external_reference_policy"]["atlas_mapping_scope"], "complete_analysis_set_once_after_terminal_qc_freeze")
        self.assertEqual(profile["external_reference_policy"]["atlas_direct_writeback_scope"], "complete_residual_qc_holdout_only")
        self.assertTrue(profile["external_reference_policy"]["open_set_ood_required"])
        self.assertTrue(profile["release_policy"]["all_accepted_broad_rescues_enter_final_deg_and_dotplots"])

    def test_sheep_immunoglobulin_aliases_enable_guarded_plasma_route(self) -> None:
        profile = json.loads((SKILL / "references/profiles/sheep_ovary.json").read_text())
        immune = profile["lineages"]["immune"]
        self.assertEqual(
            immune["sheep_immunoglobulin_aliases"],
            ["LOC101108817", "LOC101108781", "LOC121817142", "LOC114108841"],
        )
        for regulator in ["JCHAIN", "POU2AF1", "TENT5C", "MZB1"]:
            self.assertIn(regulator, immune["b_plasma"])
        self.assertIn("at least two sheep immunoglobulin loci", immune["antibody_secreting_alternative_gate"])
        self.assertIn("single immunoglobulin locus", immune["safety"])

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
            self.assertEqual(result["preferred_backbone"], "seurat_r_first_context_adaptive_graph")
            self.assertFalse(result["fixed_cellbin_preprocessing_required"])
            self.assertTrue(result["stereopy_cellbin_path_or_feature_hint"])
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
            self.assertFalse(result["fixed_cellbin_preprocessing_required"])

    def test_cellbin_runner_rejects_unrecorded_drift_and_missing_hashes(self) -> None:
        text = (SCRIPTS / "run_seurat_sct_preprocess.R").read_text()
        self.assertIn("allow-batch-exception", text)
        self.assertIn("batch-exception-reason", text)
        self.assertIn('for (package in c("glmGamPoi", "digest"))', text)
        for token in ["v2", "glmGamPoi", "4000", "50000", "banksy-k-geom", "banksy-lambda", "0.2,0.4,0.6,0.8"]:
            self.assertIn(token, text)
        self.assertIn("imported normalization/PCA/UMAP/clusters/labels ignored", text)

    def test_seurat_resolution_grids_expose_real_parallel_workers(self) -> None:
        whole = (SCRIPTS / "run_seurat_sct_preprocess.R").read_text()
        cohort = (SCRIPTS / "run_seurat_cohort_recluster_impl.R").read_text()
        for token in ["resolution-workers", "resolution-future-plan", "future::multicore", "resolution_workers_used", "umap_threads"]:
            self.assertIn(token, cohort)
        for token in ["analysis-threads", "clusterBanksy", "k_neighbors = k_neighbors", "n_threads = analysis_threads"]:
            self.assertIn(token, whole)
        self.assertNotIn('FindNeighbors(', whole)
        self.assertIn('by_cols<-c("cluster",as.character(sc))', cohort)
        self.assertNotIn('by=c("cluster",sc)', cohort)
        self.assertIn('resolution = resolutions', whole)
        self.assertIn('FindClusters(object=q[["COHORT_snn"]],algorithm=4,resolution=resolutions', cohort)
        self.assertIn('resolution_contract=="sheep_ovary"', cohort)
        self.assertIn('resolution below formal minimum', cohort)
        orchestration = (SKILL / "references/job-orchestration.md").read_text()
        self.assertIn("Mandatory CPU-to-parallelism contract", orchestration)
        self.assertIn("CPU/wall ratio", orchestration)

    def test_reference_self_calibration_and_legacy_writeback_are_blocked(self) -> None:
        bundle = (SCRIPTS / "build_depth_matched_atlas_bundle.py").read_text()
        tiered = (SCRIPTS / "calibrate_tiered_mapping_thresholds.py").read_text()
        legacy = (SCRIPTS / "calibrate_mapping_thresholds.py").read_text()
        materialize = (SCRIPTS / "legacy/materialize_adjudication_partition.py").read_text()
        terminal = (SCRIPTS / "legacy/commit_terminal_partitions.py").read_text()
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
        self.assertIn("strategy trace", case)
        self.assertIn("final broad membership", case)
        self.assertNotIn("D055", case)
        self.assertNotIn("/" + "share" + "/" + "org" + "/", case)

    def test_spatial_deg_requires_separate_verified_lognormalize_object(self) -> None:
        prepare = (SCRIPTS / "prepare_seurat_full_feature_validation.R").read_text()
        helper = (SCRIPTS / "seurat_validation_layer.R").read_text()
        initial = (SCRIPTS / "run_initial_cluster_evidence.R").read_text()
        final = (SCRIPTS / "run_final_label_deg.R").read_text()
        for token in [
            'normalization.method = "LogNormalize"',
            'scale_factor = scale_factor',
            'role = "full_feature_deg_marker_validation_only"',
            "analysis_set_sha256",
            "reductions_removed = TRUE",
            "clustering_eligible = FALSE",
        ]:
            self.assertIn(token, prepare)
        self.assertIn("sparse_exact_equal(counts, data)", helper)
        self.assertIn("hash_observation_ids(colnames(counts))", helper)
        self.assertIn("Validation manifest points to a different normalized object", helper)
        for evidence_runner in [initial, final]:
            self.assertIn("assert_seurat_validation_layer", evidence_runner)
            self.assertIn("validation-manifest", evidence_runner)
            self.assertIn("object_path = a$rds", evidence_runner)

    def test_microclusters_are_audited_but_never_removed_from_all_cell_scores(self) -> None:
        compare = (SCRIPTS / "compare_clusterings.py").read_text()
        for token in [
            "detect_cluster_separator",
            'endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz"))',
            "ARI_all_observations",
            "AMI_all_observations",
            "ARI_macro_restricted",
            "clustering_microcluster_audit.tsv",
            "clustering_pairwise_migration.tsv",
            '"included_in_all_cell_scores": True',
            "never delete, relabel, or omit from DEG/spatial evidence by size alone",
        ]:
            self.assertIn(token, compare)
        clustering_policy = (SKILL / "references/clustering-selection.md").read_text()
        self.assertIn("size alone never deletes, relabels or suppresses", clustering_policy)

    def test_cortical_or_edge_location_is_not_negative_oocyte_evidence(self) -> None:
        profile = json.loads((SKILL / "references/profiles/sheep_ovary.json").read_text())
        rule = profile["context_specific_identity_rules"]["oocyte"]["spatial_location_rule"].lower()
        for token in ["cortical", "peripheral", "section-edge", "never negative evidence"]:
            self.assertIn(token, rule)
        self.assertIn("location alone is not positive evidence", rule)
        context = (SKILL / "references/context-and-biology.md").read_text().lower()
        self.assertIn("small and primordial oocytes may occur in the ovarian cortex", context)

    @unittest.skipUnless(importlib.util.find_spec("pandas"), "pandas runtime unavailable")
    def test_oocyte_route_reclusters_full_targeted_cohort_not_only_strict_seeds(self) -> None:
        profile = json.loads((SKILL / "references/profiles/sheep_ovary.json").read_text())
        oocyte = profile["context_specific_identity_rules"]["oocyte"]
        self.assertFalse(oocyte["spatial_focus_hard_filter_for_targeted_cohort"])
        self.assertIn("all observations", oocyte["targeted_cohort_policy"])
        self.assertIn("never the final census", oocyte["strict_seed_role"])
        route = (SKILL / "references/profiles/sheep_ovary_oocyte_rfirst_route.md").read_text()
        self.assertIn("complete cohort", route)
        self.assertIn("Strict seeds/spatial foci", route)

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            screen = root / "screen.tsv"
            out = root / "out"
            fields = [
                "cell_id",
                "starting_marker_gate",
                "spatial_focus_supported",
                "spatial_focus_id",
                "total_oocyte_program_hits",
                "identity_core_hits",
                "modules_detected",
                "contradictory_somatic_hits",
            ]
            with screen.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                for index in range(12):
                    writer.writerow(
                        {
                            "cell_id": f"cell_{index}",
                            "starting_marker_gate": "TRUE",
                            "spatial_focus_supported": "TRUE" if index < 10 else "FALSE",
                            "spatial_focus_id": 1 if index < 5 else (2 if index < 10 else -1),
                            "total_oocyte_program_hits": 8 if index < 3 else 4,
                            "identity_core_hits": 5 if index < 3 else 2,
                            "modules_detected": 3,
                            "contradictory_somatic_hits": 0 if index < 3 else 2,
                        }
                    )
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS / "calibrate_rare_cell_candidates.py"),
                    "--screen",
                    str(screen),
                    "--out",
                    str(out),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            summary = json.loads((out / "calibrated_rare_focus_summary.json").read_text())
            self.assertEqual(summary["full_targeted_recluster_cohort"], 12)
            self.assertLess(summary["strict_seeds"], 12)
            self.assertEqual(
                summary["canonical_recluster_membership"],
                "oocyte_targeted_recluster_cohort.tsv.gz",
            )
            with gzip.open(out / "oocyte_targeted_recluster_cohort.tsv.gz", "rt") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 12)
            self.assertTrue(all(row["full_targeted_cohort_member"] == "True" for row in rows))


class CohortOrchestrationContract(unittest.TestCase):
    def test_new_project_uses_cohorts_and_no_persistent_pool_registries(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp) / "project"
            subprocess.run(
                [sys.executable, str(SCRIPTS / "init_annotation_project.py"), "--sample", "S1", "--input-root", temp, "--project-root", str(project), "--modality", "spatial"],
                check=True,
                capture_output=True,
                text=True,
            )
            config = json.loads((project / "config/project.json").read_text())
            self.assertEqual(config["routing_model"], "direct_cross_lineage_recluster_cohorts_global_atlas")
            self.assertFalse(config["persistent_biological_pools"])
            self.assertTrue((project / "state/recluster_cohort_registry.tsv").exists())
            self.assertTrue((project / "state/direct_return_registry.tsv").exists())
            self.assertFalse((project / "state/pool_registry.tsv").exists())
            self.assertFalse((project / "state/pool_snapshot_registry.tsv").exists())
            self.assertFalse((project / "state/branch_control_board.tsv").exists())

    def test_scheduler_names_expose_sample_stage_scope_and_attempt(self) -> None:
        generated = subprocess.run(
            [sys.executable, str(SCRIPTS / "scheduler_job_name.py"), "--sample", "C05297F1", "--stage-code", "40", "--scope", "stromal", "--attempt", "2"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        self.assertEqual(generated, "C05297F1__P40_COHORT_stromal__A02")
        subprocess.run(
            [sys.executable, str(SCRIPTS / "scheduler_job_name.py"), "--validate", generated],
            check=True,
            capture_output=True,
            text=True,
        )
        invalid = subprocess.run(
            [sys.executable, str(SCRIPTS / "scheduler_job_name.py"), "--validate", "C05297F1_sct_preprocess_v0"],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(invalid.returncode, 0)

    def test_submitted_run_requires_stage_readable_scheduler_name(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            project = root / "project"
            subprocess.run([sys.executable, str(SCRIPTS / "init_annotation_project.py"), "--sample", "S1", "--input-root", str(root), "--project-root", str(project), "--modality", "spatial"], check=True, capture_output=True, text=True)
            command = [sys.executable, str(SCRIPTS / "update_run_registry.py"), str(project), "--run-id", "run1", "--work-key", "whole_sct", "--execution-fingerprint", "fixture-sha", "--stage", "preprocess", "--script", "runner.R", "--environment", "R", "--status", "submitted", "--output-root", str(project / "runs")]
            missing = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(missing.returncode, 0)
            self.assertIn("require --job-name", missing.stderr + missing.stdout)
            subprocess.run(command + ["--job-name", "S1__P10_SCT__A01", "--job-id", "123"], check=True, capture_output=True, text=True)
            with (project / "state/run_registry.tsv").open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["scheduler_job_name"], "S1__P10_SCT__A01")

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
