from __future__ import annotations

import csv
import gzip
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
PROFILE = SKILL / "references/profiles/sheep_ovary.json"
WORKFLOW = SKILL / "references/profiles/sheep_ovary_rfirst_profile.json"
CATALOG = SKILL / "references/profiles/sheep_ovary_candidate_lineage_catalog.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(*args: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True)


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


class V2ContractTests(unittest.TestCase):
    def test_annotation_contract_separates_banksy_and_query_grids(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "config").mkdir(); (root / "state").mkdir(); (root / "provenance").mkdir()
            (root / "config/project.json").write_text(json.dumps({
                "framework_version": "2.0.0", "project_id": "p", "sample_id": "s"
            }))
            write_tsv(root / "state/input_snapshot_registry.tsv",
                      ["snapshot_id", "sample_id", "path", "kind", "size_bytes", "sha256", "status", "created_at"],
                      [{"snapshot_id": "raw", "sample_id": "s", "path": "/input", "kind": "rds", "size_bytes": "1", "sha256": "a" * 64, "status": "frozen", "created_at": "now"}])
            grid = root / "banksy_grid.json"
            grid.write_text(json.dumps({"candidate_resolutions": [0.2, 0.4, 0.8]}))
            built = run(SCRIPTS / "build_annotation_contract_v2.py", root,
                        "--workflow-profile", WORKFLOW, "--biological-profile", PROFILE, "--candidate-catalog", CATALOG,
                        "--snapshot-id", "raw", "--whole-tissue-method", "BANKSY",
                        "--whole-tissue-grid", "0.2,0.4,0.8", "--grid-source", "bound_upstream_input",
                        "--whole-tissue-grid-artifact", grid)
            self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
            contract = json.loads((root / "config/annotation_contract.json").read_text())
            self.assertEqual(Path(contract["workflow_profile"]["path"]).parent, root / "config/contract_profiles")
            self.assertEqual(Path(contract["biological_profile"]["path"]).parent, root / "config/contract_profiles")
            self.assertEqual(Path(contract["candidate_catalog"]["path"]).parent, root / "config/contract_profiles")
            self.assertEqual(contract["whole_tissue_partition"]["candidate_resolutions"], [0.2, 0.4, 0.8])
            self.assertEqual(contract["whole_tissue_partition"]["grid_artifact"]["sha256"], sha(grid))
            self.assertEqual(contract["query_reclustering"]["candidate_resolutions"], [0.1, 0.2, 0.3, 0.4, 0.6])
            validated = run(SCRIPTS / "validate_annotation_contract_v2.py", root / "config/annotation_contract.json")
            self.assertEqual(validated.returncode, 0, validated.stdout + validated.stderr)
            self.assertTrue((root / "provenance/annotation_contract_validation.json").is_file())

    def test_complete_catalog_has_two_explicit_positive_families(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / "families.json"
            result = run(SCRIPTS / "validate_broad_marker_family_contract.py", "--profile", PROFILE, "--catalog", CATALOG, "--out", out)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            payload = json.loads(out.read_text())
            self.assertEqual(payload["review_required_candidates_checked"], 14)
            self.assertTrue(all(row["positive_family_n"] >= 2 for row in payload["candidates"]))
            self.assertEqual(payload["profile_sha256"], sha(PROFILE))

    def test_v2_taxonomy_rejects_legacy_broad_and_qc_fine_label(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ledger = root / "ledger.tsv"
            write_tsv(ledger, ["broad_label", "fine_label", "annotation_status", "recluster_cohort_id"], [
                {"broad_label": "Vascular/endothelial", "fine_label": "Blood endothelial", "annotation_status": "defined_fine", "recluster_cohort_id": "v"},
                {"broad_label": "Stromal/mesenchymal", "fine_label": "Low-information stromal", "annotation_status": "defined_fine", "recluster_cohort_id": "s"},
            ])
            result = run(SCRIPTS / "audit_release_taxonomy.py", ledger, "--profile", PROFILE, "--require-v2", "--out", root / "audit.json")
            self.assertEqual(result.returncode, 2)
            self.assertIn("legacy broad alias", result.stdout)
            self.assertIn("non-biological/QC semantics", result.stdout)

    def test_project_result_audit_distinguishes_candidate_from_release(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "config").mkdir(); (root / "results/x").mkdir(parents=True); (root / "provenance").mkdir(); (root / "state").mkdir()
            (root / "config/project.json").write_text(json.dumps({"framework_version": "1.9.2", "project_id": "p"}))
            ledger = root / "results/x/final_candidate_cell_ledger.tsv.gz"
            write_tsv(ledger, ["cell_id", "analysis_scope", "final_state", "final_broad_label", "final_fine_label"], [
                {"cell_id": "a", "analysis_scope": "analysis_set", "final_state": "defined", "final_broad_label": "Pericyte/mural", "final_fine_label": "RGS5 mural"},
                {"cell_id": "b", "analysis_scope": "analysis_set", "final_state": "qc_holdout", "final_broad_label": "", "final_fine_label": ""},
            ])
            result = run(SCRIPTS / "audit_project_results_v2.py", root, "--out", root / "audit.json")
            self.assertEqual(result.returncode, 2)
            payload = json.loads((root / "audit.json").read_text())
            self.assertEqual(payload["ledger_role"], "uncommitted_candidate")
            self.assertIn("latest_annotation_ledger_is_uncommitted_candidate", payload["blockers"])
            self.assertIn("legacy_or_child_label_used_as_release_broad", payload["blockers"])

    def test_banksy_selection_must_match_bound_upstream_grid(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            carry = root / "signals.tsv"; carry.write_text("signal\nwatch\n")
            reviews = {key: True for key in (
                "full_catalog_lineage_scan", "default_broad_recall_review", "large_cluster_purity_review",
                "zero_census_review", "deg_marker_coherence_review", "spatial_morphology_review",
                "adjacent_resolution_migration_review", "technical_fragmentation_review")}
            decision = root / "decision.json"
            decision.write_text(json.dumps({
                "method": "BANKSY", "question_mode": "whole_tissue_broad_annotation", "selected_resolution": 0.8,
                "carryforward_signal_artifact": str(carry), "carryforward_signal_sha256": sha(carry),
                "candidates": [{"resolution": value, "selection_basis": "broad_recall", "selection_rationale": "full catalog evidence", **reviews} for value in (0.2, 0.4, 0.8)],
            }))
            contract = root / "contract.json"
            contract.write_text(json.dumps({"whole_tissue_partition": {"method": "BANKSY", "candidate_grid_source": "bound_upstream_input", "candidate_resolutions": [0.1, 0.2, 0.4, 0.8]}}))
            result = run(SCRIPTS / "validate_banksy_broad_resolution_selection.py", decision, "--annotation-contract", contract)
            self.assertEqual(result.returncode, 2)
            self.assertIn("differ from the v2 annotation contract", result.stdout)

    def test_authoritative_atlas_router_requires_classwise_calibration(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ledger = root / "ledger.tsv"
            write_tsv(ledger, ["cell_id", "analysis_scope", "final_state", "final_broad_label", "final_fine_label", "source_cluster"], [
                {"cell_id": "a", "analysis_scope": "analysis_set", "final_state": "qc_holdout", "final_broad_label": "", "final_fine_label": "", "source_cluster": "1"},
                {"cell_id": "b", "analysis_scope": "analysis_set", "final_state": "qc_holdout", "final_broad_label": "", "final_fine_label": "", "source_cluster": "1"},
            ])
            mapping = root / "mapping.tsv.gz"
            write_tsv(mapping, ["cell_id", "predicted_label", "mapping_tier", "out_of_distribution", "ontology_conflict"], [
                {"cell_id": "a", "predicted_label": "Stromal", "mapping_tier": "moderate_only", "out_of_distribution": "false", "ontology_conflict": "false"},
                {"cell_id": "b", "predicted_label": "Immune", "mapping_tier": "high", "out_of_distribution": "false", "ontology_conflict": "false"},
            ])
            origin = root / "origin.json"; origin.write_text(json.dumps({
                "status": "PASS", "heldout_origin": "query_like_heldout_current_query_anchors",
                "reference_self_classification": False, "anchor_target_overlap": 0,
            }))
            cumulative = root / "heldout.tsv"
            write_tsv(cumulative, ["predicted_label", "cumulative_tier", "validation_n", "validation_precision", "target_precision"], [
                {"predicted_label": "Stromal", "cumulative_tier": "moderate_or_higher", "validation_n": 40, "validation_precision": 0.95, "target_precision": 0.9},
                {"predicted_label": "Immune", "cumulative_tier": "moderate_or_higher", "validation_n": 10, "validation_precision": 1.0, "target_precision": 0.9},
            ])
            target_mapping = root / "target_mapping.tsv.gz"
            write_tsv(target_mapping, ["cell_id", "predicted_label", "mapping_tier", "out_of_distribution", "ontology_conflict"], [])
            placeholders = {}
            for name, path in {"query_mapping": target_mapping, "heldout_cumulative_validation": cumulative}.items():
                placeholders[name] = {"path": str(path), "sha256": sha(path)}
            calibration_source = root / "calibration_source.json"
            calibration_source.write_text(json.dumps({
                "schema_version": "2.0", "status": "CALIBRATED_TIERED_EVIDENCE_ONLY",
                "heldout_origin": "query_like_heldout_current_query_anchors",
                "calibration_origin_manifest": str(origin), "calibration_origin_manifest_sha256": sha(origin),
                "artifacts": placeholders,
            }))
            calibration = root / "calibration.json"
            bound = run(SCRIPTS / "bind_atlas_routing_mapping.py", "--calibration-manifest", calibration_source,
                        "--heldout-mapping", mapping, "--combined-mapping", mapping, "--out", calibration)
            self.assertEqual(bound.returncode, 0, bound.stdout + bound.stderr)
            result = run(SCRIPTS / "route_global_atlas_v2.py", "--cell-ledger", ledger, "--atlas-mapping", mapping,
                         "--calibration-manifest", calibration, "--workflow-profile", WORKFLOW, "--out", root / "out")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with gzip.open(root / "out/atlas_state_routing.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
                rows = {row["cell_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["a"]["atlas_state_route"], "direct_qc_broad_return")
            self.assertEqual(rows["a"]["proposed_broad_label"], "Stromal/mesenchymal")
            self.assertEqual(rows["b"]["atlas_state_route"], "retain_qc")

    def test_atlas_material_review_uses_current_broad_group_and_has_review_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ledger = root / "ledger.tsv"
            ledger_rows = []
            mapping_rows = []
            for index in range(10):
                primary = "Granulosa" if index < 2 else "Stromal/mesenchymal"
                ledger_rows.append({
                    "cell_id": f"c{index}", "analysis_scope": "analysis_set", "final_state": "defined_broad_only",
                    "final_broad_label": primary, "final_fine_label": "", "source_cluster": "mixed_1",
                })
                mapping_rows.append({
                    "cell_id": f"c{index}", "predicted_label": "Stromal", "mapping_tier": "high",
                    "out_of_distribution": "true" if index < 2 else "false", "ontology_conflict": "false",
                })
            write_tsv(ledger, list(ledger_rows[0]), ledger_rows)
            mapping = root / "mapping.tsv.gz"
            write_tsv(mapping, list(mapping_rows[0]), mapping_rows)
            origin = root / "origin.json"
            origin.write_text(json.dumps({
                "status": "PASS", "heldout_origin": "query_like_heldout_current_query_anchors",
                "reference_self_classification": False, "anchor_target_overlap": 0,
            }))
            cumulative = root / "heldout.tsv"
            write_tsv(cumulative, ["predicted_label", "cumulative_tier", "validation_n", "validation_precision"], [{
                "predicted_label": "Stromal", "cumulative_tier": "moderate_or_higher",
                "validation_n": 40, "validation_precision": 0.95,
            }])
            target_mapping = root / "target_mapping.tsv.gz"
            write_tsv(target_mapping, list(mapping_rows[0]), [])
            calibration_source = root / "calibration_source.json"
            calibration_source.write_text(json.dumps({
                "schema_version": "2.0", "status": "CALIBRATED_TIERED_EVIDENCE_ONLY",
                "heldout_origin": "query_like_heldout_current_query_anchors",
                "calibration_origin_manifest": str(origin), "calibration_origin_manifest_sha256": sha(origin),
                "artifacts": {
                    "query_mapping": {"path": str(target_mapping), "sha256": sha(target_mapping)},
                    "heldout_cumulative_validation": {"path": str(cumulative), "sha256": sha(cumulative)},
                },
            }))
            calibration = root / "calibration.json"
            bound = run(SCRIPTS / "bind_atlas_routing_mapping.py", "--calibration-manifest", calibration_source,
                        "--heldout-mapping", mapping, "--combined-mapping", mapping, "--out", calibration)
            self.assertEqual(bound.returncode, 0, bound.stdout + bound.stderr)
            result = run(
                SCRIPTS / "route_global_atlas_v2.py", "--cell-ledger", ledger, "--atlas-mapping", mapping,
                "--calibration-manifest", calibration, "--workflow-profile", WORKFLOW, "--out", root / "out",
                "--min-discordant-n", 2, "--min-discordant-fraction", 0.9,
            )
            self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
            with (root / "out/atlas_discrepancy_review_queue.tsv").open(newline="", encoding="utf-8") as handle:
                queue = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(queue), 1)
            self.assertEqual(queue[0]["trigger_reason"], "out_of_distribution")
            self.assertEqual(queue[0]["primary_group_n"], "2")
            self.assertEqual(queue[0]["cluster_n"], "10")
            with gzip.open(root / "out/atlas_state_routing.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
                routes = list(csv.DictReader(handle, delimiter="\t"))
            reviewed = [row for row in routes if row["review_required"] == "true"]
            self.assertEqual({row["cell_id"] for row in reviewed}, {"c0", "c1"})
            self.assertTrue(all(row["review_id"] == queue[0]["review_id"] for row in reviewed))

    def test_query_only_cohort_validator_accepts_zero_anchors(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "tables").mkdir(); (root / "figures").mkdir()
            for name in ("RUN_COMPLETE.tsv", "sessionInfo.txt", "cohort_reclustered_query_seurat.rds"):
                (root / name).write_text("ok\n")
            write_tsv(root / "run_manifest.tsv", ["parameter", "value"], [
                {"parameter": "anchor_assisted", "value": "false"},
                {"parameter": "query_only_graph_umap_deg", "value": "true"},
                {"parameter": "n_query_analyzed", "value": "2"},
                {"parameter": "n_anchors_analyzed", "value": "0"},
            ])
            write_tsv(root / "tables/analyzed_membership.tsv.gz", ["cell_id"], [{"cell_id": "a"}, {"cell_id": "b"}])
            write_tsv(root / "tables/cluster_source_state_composition.tsv", ["cluster", "n"], [{"cluster": "0", "n": "2"}])
            write_tsv(root / "tables/cluster_QC_summary.tsv", ["cluster", "n"], [{"cluster": "0", "n": "2"}])
            write_tsv(root / "tables/framework_res0p1_clusters.tsv", ["cell_id", "cluster"], [{"cell_id": "a", "cluster": "0"}, {"cell_id": "b", "cluster": "0"}])
            write_tsv(root / "tables/framework_res0p1_DEG_all.tsv", ["gene", "cluster"], [{"gene": "G", "cluster": "0"}])
            write_tsv(root / "tables/framework_res0p1_DEG_top100.tsv", ["gene", "cluster"], [{"gene": "G", "cluster": "0"}])
            (root / "figures/framework_res0p1_UMAP.png").write_bytes(b"png")
            (root / "figures/framework_res0p1_UMAP.pdf").write_bytes(b"pdf")
            result = run(SCRIPTS / "validate_cohort_recluster_output.py", root,
                         "--expected-query", 2, "--expected-anchors", 0, "--resolutions", "0.1")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_oocyte_context_validator_reads_gzip_memberships(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            canonical = root / "canonical.tsv.gz"; context = root / "context.tsv.gz"
            write_tsv(canonical, ["cell_id"], [{"cell_id": "a"}])
            write_tsv(context, ["cell_id"], [{"cell_id": "a"}, {"cell_id": "b"}])
            outcomes = root / "outcomes.tsv"
            write_tsv(outcomes, ["cell_id", "final_broad_label"], [
                {"cell_id": "a", "final_broad_label": "Oocyte"},
                {"cell_id": "b", "final_broad_label": "Granulosa"},
            ])
            out = root / "validation.json"
            result = run(SCRIPTS / "validate_oocyte_context_boundary.py",
                         "--canonical-membership", canonical, "--context-membership", context,
                         "--outcomes", outcomes, "--out", out)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_analysis_scope_policy_is_hash_bound(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); membership = root / "scope.tsv.gz"
            write_tsv(membership, ["cell_id", "analysis_scope"], [
                {"cell_id": "a", "analysis_scope": "analysis_set"},
                {"cell_id": "b", "analysis_scope": "excluded_initial_qc"},
            ])
            result = run(SCRIPTS / "build_analysis_scope_policy.py", root, "--membership", membership)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            payload = json.loads((root / "provenance/analysis_scope_policy.json").read_text())
            self.assertEqual(payload["analysis_set_n"], 1)
            self.assertEqual(payload["excluded_initial_qc_n"], 1)
            self.assertEqual(payload["membership_sha256"], sha(membership))

    def test_preflight_rejects_invalid_scheduler_stage_before_submission(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / "job.py"; source.write_text("print('ok')\n")
            result = run(SCRIPTS / "preflight_generated_job.py", source,
                         "--scheduler-job-name", "S1__P74_COMPLETION__A01")
            self.assertEqual(result.returncode, 2)
            self.assertIn("scheduler_job_name", result.stdout)

    def test_preflight_checks_the_declared_python_environment(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / "job.py"; source.write_text("print('ok')\n")
            result = run(SCRIPTS / "preflight_generated_job.py", source,
                         "--python", sys.executable, "--python-imports", "module_that_cannot_exist_v201")
            self.assertEqual(result.returncode, 2)
            self.assertIn("python_imports", result.stdout)


if __name__ == "__main__":
    unittest.main()
