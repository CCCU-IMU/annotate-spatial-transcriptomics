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
SKILL = ROOT / "annotate-spatial-transcriptomics"
SCRIPTS = SKILL / "scripts"
BIOLOGICAL = SKILL / "references/profiles/sheep_ovary.json"
WORKFLOW = SKILL / "references/profiles/sheep_ovary_rfirst_profile.json"
CATALOG = SKILL / "references/profiles/sheep_ovary_candidate_lineage_catalog.json"


class V150Contract(unittest.TestCase):
    def test_explicit_same_batch_preset_binds_process_not_sample_answers(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            context = root / "context.json"
            inspection = root / "inspection.json"
            output = root / "active_strategy_preset.json"
            context.write_text(json.dumps({"species": "Ovis aries", "tissue": "ovary"}))
            inspection.write_text(json.dumps({"type": "Seurat", "assays": ["Spatial"], "n_features": 25000}))
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "resolve_workflow_profile.py"),
                "--context", str(context), "--r-inspection", str(inspection),
                "--input-path", "/example/cellbin_PPed/sample.seurat.rds",
                "--strategy-preset", "sheep_ovary_same_batch_rfirst", "--out", str(output),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            resolved = json.loads(output.read_text())
            self.assertEqual(resolved["strategy_preset_status"], "ACTIVE")
            self.assertEqual(resolved["strategy_preset_preprocessing_mode"], "fixed_verified_same_batch_contract")
            preset = resolved["strategy_preset"]
            self.assertIn("route_c_complete_qc_anchor_recluster_then_residual_only_calibrated_atlas", preset["phase_order"])
            self.assertIn("selected_pool_resolution", preset["always_adaptive"])
            self.assertIn("cluster_number_to_label_mapping", preset["forbidden_transfer_from_reference_case"])
            self.assertIn("disabled", preset["generic_rare_cell_route"])
            self.assertIn("candidate_lineage_catalog_sha256", resolved["strategy_preset_bindings"])
            self.assertTrue(resolved["strategy_preset_bindings"]["preset_sha256"])

    def test_open_world_lineage_audit_keeps_sample_specific_candidates_visible(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "project"
            (root / "state").mkdir(parents=True)
            ledger = root / "state/cluster_decision_ledger.tsv"
            ledger.write_text("decision_id\tsource_cluster\nD1\t0\n")
            audit = root / "open_world.json"
            scaffold = subprocess.run([
                sys.executable, str(SCRIPTS / "init_open_world_lineage_audit.py"),
                str(root), "--out", str(audit),
            ], capture_output=True, text=True)
            self.assertEqual(scaffold.returncode, 0, scaffold.stdout + scaffold.stderr)
            catalog = json.loads(CATALOG.read_text())
            required = [
                item["candidate_id"] for item in catalog["candidate_boundaries"]
                if item.get("review_required") is True
            ]
            self.assertGreaterEqual(len(required), 12)
            audit_data = json.loads(audit.read_text())
            reviews = audit_data["candidate_reviews"]
            self.assertEqual({item["candidate_id"] for item in reviews}, set(required))
            for item in reviews:
                item.update({
                    "outcome": "not_supported",
                    "evidence_summary": "full-gene marker and anti-marker review found no coherent program",
                    "action": "record negative audit",
                })
            audit_data["additional_candidates"] = [{"candidate_id": "sample_specific_candidate", "action": "create Route A pool"}]
            audit.write_text(json.dumps(audit_data))
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "validate_open_world_lineage_audit.py"),
                str(root), "--audit", str(audit),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            validation = json.loads((root / "provenance/open_world_lineage_audit_validation.json").read_text())
            self.assertEqual(validation["candidate_catalog_sha256"], hashlib.sha256(CATALOG.read_bytes()).hexdigest())
            incomplete = json.loads(audit.read_text())
            incomplete["candidate_reviews"] = incomplete["candidate_reviews"][:-1]
            audit.write_text(json.dumps(incomplete))
            missing = subprocess.run([
                sys.executable, str(SCRIPTS / "validate_open_world_lineage_audit.py"),
                str(root), "--audit", str(audit),
            ], capture_output=True, text=True)
            self.assertNotEqual(missing.returncode, 0)
            incomplete["candidate_reviews"] = reviews
            audit.write_text(json.dumps(incomplete))
            ledger.write_text("decision_id\tsource_cluster\nD1\t0\nD2\t1\n")
            stale = subprocess.run([
                sys.executable, str(SCRIPTS / "validate_open_world_lineage_audit.py"),
                str(root), "--audit", str(audit),
            ], capture_output=True, text=True)
            self.assertNotEqual(stale.returncode, 0)

    def test_master_quality_review_cannot_run_at_broad_only_stage(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "project"
            (root / "provenance").mkdir(parents=True)
            (root / "state").mkdir()
            (root / "provenance/completion_gate.json").write_text(json.dumps({"status": "PASS"}))
            result = subprocess.run(
                [sys.executable, str(SCRIPTS / "request_master_quality_review.py"), str(root)],
                capture_output=True, text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("required final artifacts are missing", result.stdout + result.stderr)

    def test_profile_roles_and_normal_resolution_grid_fail_closed(self) -> None:
        biological = json.loads(BIOLOGICAL.read_text())
        workflow = json.loads(WORKFLOW.read_text())
        self.assertEqual(biological["profile_role"], "biological_evidence")
        self.assertEqual(workflow["profile_role"], "workflow_preprocessing")
        expected = [0.1, 0.2, 0.3, 0.4, 0.6]
        self.assertEqual(biological["resolution_policy"]["formal_candidate_resolutions"], expected)
        self.assertEqual(workflow["pool_reclustering_contract"]["candidate_resolutions"], expected)
        good = subprocess.run(
            [sys.executable, str(SCRIPTS / "validate_resolution_grid.py"), "--workflow-profile", str(WORKFLOW), "--scope", "pool", "--resolutions", "0.1,0.2,0.3,0.4,0.6"],
            capture_output=True, text=True,
        )
        self.assertEqual(good.returncode, 0, good.stdout + good.stderr)
        for grid in ["0.01,0.02,0.05,0.1,0.2", "0.1,0.2,0.4,0.6"]:
            bad = subprocess.run(
                [sys.executable, str(SCRIPTS / "validate_resolution_grid.py"), "--workflow-profile", str(WORKFLOW), "--scope", "pool", "--resolutions", grid],
                capture_output=True, text=True,
            )
            self.assertNotEqual(bad.returncode, 0)

    def test_workflow_profile_cannot_satisfy_biological_profile_binding(self) -> None:
        result = subprocess.run(
            [sys.executable, str(SCRIPTS / "validate_profile_role.py"), str(WORKFLOW), "--expected-role", "biological_evidence"],
            capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("expected 'biological_evidence'", result.stdout + result.stderr)

    def test_anchor_route_rejects_ultralow_or_unrecorded_candidate_grid(self) -> None:
        sys.path.insert(0, str(SCRIPTS))
        try:
            from multiroute_lib import valid_attempt
        finally:
            sys.path.pop(0)
        profile = json.loads(BIOLOGICAL.read_text())
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            validation = root / "validation.json"
            validation.write_text("{}\n")
            membership = root / "membership.tsv"
            membership.write_text("cell_id\nc1\nc2\n")
            base = {
                "route_class": "biological_anchor_recluster", "status": "validated",
                "validation_artifact": str(validation), "applicability": "applicable",
                "n_query": "2", "n_anchors": "20", "query_only_graph": "true",
                "query_membership_artifact": str(membership),
                "query_membership_sha256": hashlib.sha256(membership.read_bytes()).hexdigest(),
                "selected_resolution": "0.2",
                "candidate_resolutions": "0.1,0.2,0.3,0.4,0.6",
            }
            self.assertTrue(valid_attempt(root, base, profile)[0])
            bad = dict(base, candidate_resolutions="0.01,0.05,0.1,0.2", selected_resolution="0.05")
            ok, errors = valid_attempt(root, bad, profile)
            self.assertFalse(ok)
            self.assertTrue(any("formal candidate_resolutions" in error for error in errors))

    def test_atlas_is_residual_qc_only_and_requires_qc_anchor_prerequisite(self) -> None:
        profile = json.loads(BIOLOGICAL.read_text())
        workflow = json.loads(WORKFLOW.read_text())
        self.assertEqual(profile["multi_route_policy"]["atlas_rescue_scope"], "residual_post_qc_anchor_holdout_only")
        self.assertTrue(profile["multi_route_policy"]["atlas_for_defined_broad_or_fine_forbidden"])
        self.assertTrue(workflow["external_reference_policy"]["forbid_full_object_or_biological_pool_atlas_classification"])
        sys.path.insert(0, str(SCRIPTS))
        try:
            from multiroute_lib import valid_attempt
        finally:
            sys.path.pop(0)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            validation = root / "validation.json"; validation.write_text("{}\n")
            calibration = root / "calibration.json"; calibration.write_text("{}\n")
            membership = root / "residual_qc.tsv"; membership.write_text("cell_id\nq1\n")
            attempt = {
                "route_class": "qc_atlas_review", "status": "validated", "applicability": "applicable",
                "validation_artifact": str(validation), "depth_matched_validation": "true",
                "observed_density_spatial_prior": "true",
                "calibration_origin": "query_like_heldout_current_query_anchors",
                "calibration_manifest": str(calibration), "n_query": "1",
                "query_membership_artifact": str(membership),
                "query_membership_sha256": hashlib.sha256(membership.read_bytes()).hexdigest(),
                "pool_snapshot_id": "residual_qc_pool_v1",
                "parent_pool_snapshot_id": "complete_qc_pool_v1",
                "source_state": "qc_holdout_residual_after_anchor",
            }
            ok, errors = valid_attempt(root, attempt, profile)
            self.assertFalse(ok)
            self.assertTrue(any("prerequisite QC-anchor" in error for error in errors))
        planner = (SCRIPTS / "plan_next_iteration.py").read_text()
        self.assertLess(
            planner.index('"qc_anchor_recluster" in gap_routes'),
            planner.index('"qc_atlas_review" in gap_routes'),
        )
        self.assertIn('required_route":"build_final_annotation"', planner)

    def test_rctd_medium_low_cannot_name_atlas_as_direct_successor(self) -> None:
        sys.path.insert(0, str(SCRIPTS))
        try:
            from multiroute_lib import valid_attempt
        finally:
            sys.path.pop(0)
        profile = json.loads(BIOLOGICAL.read_text())
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            validation = root / "rctd_validation.json"; validation.write_text("{}\n")
            reroute = root / "rctd_qc.tsv"; reroute.write_text("cell_id\nq1\nq2\n")
            base = {
                "route_class": "interface_deconvolution_review", "status": "validated",
                "validation_artifact": str(validation), "applicability": "applicable",
                "depth_matched_validation": "true", "n_query": "2", "n_rerouted": "2",
                "n_defined_fine": "0", "n_defined_broad_only": "0",
                "rctd_extreme_n": "0", "rctd_high_n": "0", "rctd_medium_low_n": "2",
                "rctd_fine_return_n": "0", "rctd_broad_return_n": "0",
                "reroute_membership_artifact": str(reroute),
                "reroute_membership_sha256": hashlib.sha256(reroute.read_bytes()).hexdigest(),
                "fallback_route_attempt_id": "future_qc_anchor",
                "successor_state": "atlas",
            }
            ok, errors = valid_attempt(root, base, profile)
            self.assertFalse(ok)
            self.assertTrue(any("successor_state=qc_holdout" in error for error in errors))

    def test_atlas_query_must_equal_qc_anchor_residual_membership(self) -> None:
        sys.path.insert(0, str(SCRIPTS))
        try:
            from multiroute_lib import validate_qc_atlas_prerequisite
        finally:
            sys.path.pop(0)
        profile = json.loads(BIOLOGICAL.read_text())
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            validation = root / "validation.json"; validation.write_text("{}\n")
            calibration = root / "calibration.json"; calibration.write_text("{}\n")
            full_qc = root / "full_qc.tsv"; full_qc.write_text("cell_id\nq1\nq2\nq3\n")
            residual = root / "residual_qc.tsv"; residual.write_text("cell_id\nq2\nq3\n")
            outcome = root / "qc_outcome.tsv"; outcome.write_text("cell_id\noutcome\nq1\tdefined_broad_only\nq2\tqc_holdout\nq3\tqc_holdout\n")
            qca = {
                "route_attempt_id": "qc_anchor_1", "route_class": "qc_anchor_recluster",
                "status": "validated", "applicability": "applicable", "source_state": "qc_holdout",
                "validation_artifact": str(validation), "outcome_artifact": str(outcome),
                "pool_snapshot_id": "qc_full_v1", "n_query": "3", "n_anchors": "20",
                "query_only_graph": "true", "query_membership_artifact": str(full_qc),
                "query_membership_sha256": hashlib.sha256(full_qc.read_bytes()).hexdigest(),
                "residual_qc_membership_artifact": str(residual),
                "residual_qc_membership_sha256": hashlib.sha256(residual.read_bytes()).hexdigest(),
                "n_qc_retained": "2", "candidate_resolutions": "0.1,0.2,0.3,0.4,0.6",
                "selected_resolution": "0.2", "created_at": "2026-07-15T01:00:00+00:00",
            }
            atlas = {
                "route_attempt_id": "atlas_1", "route_class": "qc_atlas_review",
                "status": "validated", "applicability": "applicable",
                "source_state": "qc_holdout_residual_after_anchor", "n_query": "2",
                "validation_artifact": str(validation), "depth_matched_validation": "true",
                "observed_density_spatial_prior": "true",
                "calibration_origin": "query_like_heldout_current_query_anchors",
                "calibration_manifest": str(calibration),
                "query_membership_artifact": str(residual),
                "query_membership_sha256": hashlib.sha256(residual.read_bytes()).hexdigest(),
                "pool_snapshot_id": "qc_residual_v1", "parent_pool_snapshot_id": "qc_full_v1",
                "prerequisite_route_attempt_id": "qc_anchor_1",
                "prerequisite_outcome_sha256": hashlib.sha256(outcome.read_bytes()).hexdigest(),
                "prerequisite_residual_qc_sha256": hashlib.sha256(residual.read_bytes()).hexdigest(),
                "created_at": "2026-07-15T02:00:00+00:00",
            }
            self.assertTrue(validate_qc_atlas_prerequisite(root, atlas, [qca], profile)[0])
            wrong = dict(atlas, query_membership_artifact=str(full_qc), query_membership_sha256=hashlib.sha256(full_qc.read_bytes()).hexdigest(), n_query="3")
            ok, errors = validate_qc_atlas_prerequisite(root, wrong, [qca], profile)
            self.assertFalse(ok)
            self.assertTrue(any("residual-QC" in error for error in errors))

    def test_single_final_annotation_thresholds(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "project"
            subprocess.run(
                [sys.executable, str(SCRIPTS / "init_annotation_project.py"), "--sample", "S1", "--input-root", temp, "--project-root", str(root), "--modality", "spatial", "--observation-unit", "cellbin"],
                check=True, capture_output=True, text=True,
            )
            ledger = root / "state/cell_ledger.tsv.gz"
            fields = ["sample_id", "cell_id", "analysis_scope", "state", "broad_label", "fine_label", "confidence", "fine_anchor_eligible"]
            fixtures = [
                ["S1", "high_fine", "analysis_set", "defined_fine", "Granulosa", "Mural granulosa", "high", "true"],
                ["S1", "moderate_fine", "analysis_set", "defined_fine", "Granulosa", "Cumulus", "moderate", "true"],
                ["S1", "rescued_broad", "analysis_set", "defined_broad_only", "Stromal/mesenchymal", "", "high", "false"],
                ["S1", "low_broad", "analysis_set", "defined_broad_only", "Immune", "", "low", "false"],
                ["S1", "excluded", "excluded_initial_qc", "excluded_initial_qc", "", "", "high", "false"],
            ]
            with gzip.open(ledger, "wt", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t"); writer.writerow(fields); writer.writerows(fixtures)
            subprocess.run(
                [sys.executable, str(SCRIPTS / "build_final_annotation.py"), str(root), "--cell-ledger", str(ledger), "--out", str(ledger), "--sample", "S1", "--version", "v001"],
                check=True, capture_output=True, text=True,
            )
            with gzip.open(ledger, "rt", newline="", encoding="utf-8") as handle:
                rows = {row["cell_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["high_fine"]["final_fine_label"], "Mural granulosa")
            self.assertEqual(rows["moderate_fine"]["final_state"], "defined_broad_only")
            self.assertEqual(rows["moderate_fine"]["final_fine_label"], "")
            self.assertEqual(rows["rescued_broad"]["final_broad_label"], "Stromal/mesenchymal")
            self.assertEqual(rows["rescued_broad"]["final_fine_label"], "")
            self.assertEqual(rows["low_broad"]["final_state"], "pending_review")
            self.assertEqual(rows["low_broad"]["final_broad_label"], "")
            self.assertEqual(rows["excluded"]["final_state"], "excluded_initial_qc")
            with (root / "state/annotation_view_registry.tsv").open(newline="", encoding="utf-8") as handle:
                registry = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["view"] for row in registry], ["final"])

    def test_incident_and_generated_job_preflights(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            good = root / "good.sh"; good.write_text("#!/usr/bin/env bash\necho ok\n")
            bad = root / "bad.sh"; bad.write_text("#!/usr/bin/env bash\nif true; then\n")
            params = root / "params.json"; params.write_text('{"resolutions":[0.1,0.2,0.3,0.4,0.6]}\n')
            ok = subprocess.run([sys.executable, str(SCRIPTS / "preflight_generated_job.py"), str(good), str(params)], capture_output=True, text=True)
            self.assertEqual(ok.returncode, 0, ok.stdout + ok.stderr)
            failed = subprocess.run([sys.executable, str(SCRIPTS / "preflight_generated_job.py"), str(bad)], capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            incident = root / "incidents.tsv"
            fields = [
                "incident_id", "scheduler_job_id", "failure_class", "failure_stage", "symptom",
                "root_cause", "failure_boundary", "accepted_prior_artifacts", "repair_action",
                "repair_verification", "state_mutated", "biological_labels_changed",
                "skill_prevention_candidate", "regression_test_candidate", "status", "evidence_paths",
            ]
            with incident.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader()
                row = {field: "recorded" for field in fields}; row.update({"incident_id": "INC-1", "status": "open"}); writer.writerow(row)
            open_result = subprocess.run([sys.executable, str(SCRIPTS / "validate_incident_registry.py"), str(incident)], capture_output=True, text=True)
            self.assertNotEqual(open_result.returncode, 0)

    def test_release_surface_has_one_final_annotation(self) -> None:
        skill_text = (SKILL / "SKILL.md").read_text()
        report_text = (SCRIPTS / "build_report.py").read_text()
        audit_text = (SCRIPTS / "audit_release.py").read_text()
        self.assertIn("one final annotation", skill_text.lower())
        self.assertNotIn("id='views'", report_text)
        self.assertIn("id='final'", report_text)
        self.assertIn("id='final'", audit_text)

    def test_broad_marker_panel_does_not_require_synthetic_subtypes(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            broad = root / "broad.tsv"
            broad.write_text(
                "label\tgene\tanalysis_view\tp_val_adj\tavg_log2FC\tpct.1\tpct.2\n"
                "Granulosa\tFSHR\tfinal\t0.001\t2\t0.8\t0.1\n"
                "Stromal/mesenchymal\tDCN\tfinal\t0.001\t2\t0.8\t0.1\n"
            )
            canonical = root / "canonical.tsv"
            canonical.write_text(
                "level\tpanel\tmarker_group\tgene\n"
                "broad\tcanonical\tGranulosa\tFSHR\n"
                "broad\tcanonical\tStromal/mesenchymal\tDCN\n"
            )
            out = root / "panel.tsv"
            result = subprocess.run(
                [sys.executable, str(SCRIPTS / "build_marker_panel_from_deg.py"), "--broad-deg", str(broad), "--canonical", str(canonical), "--out", str(out)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue(out.is_file())

    def test_report_metadata_never_fabricates_broad_only_subtypes(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ledger = root / "final.tsv"
            ledger.write_text(
                "cell_id\tanalysis_scope\tfinal_state\tfinal_broad_label\tfinal_fine_label\t"
                "final_confidence\tfinal_assignment_tier\tfinal_broad_eligible\tfinal_fine_eligible\t"
                "fine_anchor_eligible\troute\n"
                "c1\tanalysis_set\tdefined_broad_only\tGranulosa\t\thigh\thigh_broad\ttrue\tfalse\tfalse\tqc_atlas_review\n"
                "c2\tanalysis_set\tdefined_fine\tGranulosa\tMural granulosa\thigh\thigh_fine\ttrue\ttrue\ttrue\tbiological_anchor_recluster\n"
                "c3\tanalysis_set\tqc_holdout\t\t\tlow\tretained_or_unresolved\tfalse\tfalse\tfalse\tqc_atlas_review\n"
            )
            out = root / "report.tsv"
            result = subprocess.run(
                [sys.executable, str(SCRIPTS / "prepare_report_metadata.py"), "--cell-ledger", str(ledger), "--out", str(out)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with out.open(newline="", encoding="utf-8") as handle:
                rows = {row["cell_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["c1"]["primary_broad_label"], "Granulosa")
            self.assertEqual(rows["c1"]["primary_subtype_label"], "")
            self.assertEqual(rows["c1"]["subtype_display"], "")
            self.assertEqual(rows["c3"]["retained_state_display"], "QC holdout")
            self.assertNotIn("Broad only:", out.read_text())

    def test_lightweight_review_is_required_before_user_confirmation(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "project"
            for name in ["config", "state", "provenance", "tables", "review/confirmation/assets"]:
                (root / name).mkdir(parents=True, exist_ok=True)
            (root / "config/project.json").write_text(json.dumps({"sample_id": "S1", "decision_version": "v001"}))
            cell_ledger = root / "state/cell_ledger.tsv.gz"
            with gzip.open(cell_ledger, "wt", encoding="utf-8", newline="") as handle:
                handle.write(
                    "cell_id\tanalysis_scope\tstate\tfinal_state\tfinal_broad_label\tfinal_fine_label\n"
                    "c1\tanalysis_set\tdefined_broad_only\tdefined_broad_only\tGranulosa\t\n"
                    "c2\tanalysis_set\tdefined_fine\tdefined_fine\tGranulosa\tMural granulosa\n"
                    "c3\tanalysis_set\tqc_holdout\tqc_holdout\t\t\n"
                )
            cluster_ledger = root / "state/cluster_decision_ledger.tsv"
            cluster_ledger.write_text("decision_id\tsample_id\nD1\tS1\n")
            completion = root / "provenance/completion_gate.json"
            completion.write_text(json.dumps({"status": "PASS"}))
            (root / "provenance/multiroute_audit.json").write_text(json.dumps({"status": "PASS"}))
            (root / "provenance/state_validation.json").write_text(json.dumps({"status": "PASS"}))
            (root / "provenance/iteration_plan.json").write_text(json.dumps({"status": "READY_FOR_COMPLETION_AUDIT"}))
            taxonomy = root / "provenance/release_taxonomy_audit.json"
            taxonomy.write_text(json.dumps({
                "pass": True, "metadata_sha256": hashlib.sha256(cell_ledger.read_bytes()).hexdigest(),
                "biological_broad_census": {"Granulosa": 2}, "retained_state_census": {"qc_holdout": 1},
            }))
            support = root / "state/annotation_support_registry.tsv"
            support.write_text(
                "support_id\tsample_id\tlabel_level\tbroad_label\tfine_label\tn_observations\tconfidence\tpositive_marker_evidence\tanti_marker_evidence\tresolution_evidence\tspatial_evidence\tliterature_context\troute_summary\tsource_decision_ids\tvalidation_artifacts\tstatus\tsupersedes\tcreated_at\n"
                "B1\tS1\tbroad\tGranulosa\t\t2\thigh\tFSHR;INHA\tlow DCN\tstable res0.2-0.4\tfollicular wall\tadult sheep ovary\tRoute A\tD1\tevidence.tsv\tvalidated\t\t2026-01-01\n"
                "F1\tS1\tfine\tGranulosa\tMural granulosa\t1\thigh\tIHH;HSD17B1\tlow oocyte core\tstable res0.4\tfollicular layer\tadult sheep ovary\tRoute A\tD1\tevidence.tsv\tvalidated\t\t2026-01-01\n"
            )
            (root / "state/next_action_queue.tsv").write_text("action_id\n")
            (root / "state/pool_registry.tsv").write_text("pool_id\tstatus\n")
            (root / "state/route_attempt_registry.tsv").write_text("route_attempt_id\tstatus\n")
            (root / "state/run_registry.tsv").write_text("run_id\tstatus\n")
            (root / "state/annotation_view_registry.tsv").write_text("view\tstatus\nfinal\tvalidated\n")
            (root / "tables/final_annotation_census.tsv").write_text("broad_label\tn_observations\nGranulosa\t2\n")
            assets = root / "review/confirmation/assets"
            spatial = assets / "broad_spatial_review.png"; spatial.write_bytes(b"\x89PNG\r\n\x1a\nspatial")
            dotplot = assets / "broad_canonical_marker_dotplot_review.png"; dotplot.write_bytes(b"\x89PNG\r\n\x1a\ndotplot")
            source = assets / "dotplot.tsv"; source.write_text("gene\tlabel\nFSHR\tGranulosa\n")
            palette = assets / "palette.tsv"; palette.write_text("label\tcolor\nGranulosa\t#0072B2\n")
            manifest = {
                "status": "PASS", "label_column": "primary_broad_label",
                "spatial_png": str(spatial.relative_to(root)), "spatial_png_sha256": hashlib.sha256(spatial.read_bytes()).hexdigest(),
                "dotplot_png": str(dotplot.relative_to(root)), "dotplot_png_sha256": hashlib.sha256(dotplot.read_bytes()).hexdigest(),
                "dotplot_source": str(source.relative_to(root)), "dotplot_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "palette_tsv": str(palette.relative_to(root)), "palette_tsv_sha256": hashlib.sha256(palette.read_bytes()).hexdigest(),
            }
            asset_manifest = assets / "review_asset_manifest.json"
            asset_manifest.write_text(json.dumps(manifest))
            premature = subprocess.run(
                [sys.executable, str(SCRIPTS / "request_final_annotation_confirmation.py"), str(root)],
                capture_output=True, text=True,
            )
            self.assertNotEqual(premature.returncode, 0)
            premature_review = subprocess.run(
                [sys.executable, str(SCRIPTS / "build_confirmation_review.py"), str(root)],
                capture_output=True, text=True,
            )
            self.assertNotEqual(premature_review.returncode, 0)
            master_request = subprocess.run(
                [sys.executable, str(SCRIPTS / "request_master_quality_review.py"), str(root)],
                capture_output=True, text=True,
            )
            self.assertEqual(master_request.returncode, 0, master_request.stdout + master_request.stderr)
            assessment = root / "provenance/master_quality_assessment.json"
            assessment.write_text(json.dumps({
                "reviewer_role": "main_conversation_agent",
                "reviewer_id": "root-test",
                "reviewed_after_all_annotation_complete": True,
                "verdict": "PASS",
                "comparison_to_reference": "Biological evidence quality is comparable to the validated R-first reference.",
                "rationale": "Broad identity, spatial morphology and rare-lineage safeguards are biologically coherent.",
                "biological_concerns": ["minor broad-only remainder retained"],
                "requested_revisions": [],
                "checklist": {
                    "broad_annotation_reasonableness": {"status": "PASS", "note": "coherent broad label"},
                    "marker_antimarker_and_spatial_support": {"status": "PASS", "note": "marker and spatial support agree"},
                    "context_specific_and_confounded_lineage_safety": {"status": "CONCERN", "note": "no unsafe context-specific call; retain caution"},
                    "comparable_to_validated_rfirst_reference": {"status": "PASS", "note": "comparable evidence depth"},
                },
            }))
            master_record = subprocess.run(
                [sys.executable, str(SCRIPTS / "record_master_quality_approval.py"), str(root), "--assessment", str(assessment)],
                capture_output=True, text=True,
            )
            self.assertEqual(master_record.returncode, 0, master_record.stdout + master_record.stderr)
            result = subprocess.run(
                [sys.executable, str(SCRIPTS / "build_confirmation_review.py"), str(root)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            report = root / "review/confirmation/index.html"
            report_text = report.read_text()
            self.assertIn("大类注释支持原因", report_text)
            self.assertIn("主 Agent 注释质量审批", report_text)
            self.assertIn("data:image/png;base64,", report_text)
            request = subprocess.run(
                [sys.executable, str(SCRIPTS / "request_final_annotation_confirmation.py"), str(root)],
                capture_output=True, text=True,
            )
            self.assertEqual(request.returncode, 0, request.stdout + request.stderr)
            frozen = json.loads((root / "provenance/final_annotation_confirmation_request.json").read_text())
            self.assertEqual(frozen["confirmation_review_report_sha256"], hashlib.sha256(report.read_bytes()).hexdigest())
            self.assertEqual(frozen["annotation_support_registry_sha256"], hashlib.sha256(support.read_bytes()).hexdigest())
            confirmed = subprocess.run(
                [sys.executable, str(SCRIPTS / "record_final_annotation_confirmation.py"), str(root), "--confirmed-by", "test_user", "--user-message", "approved"],
                capture_output=True, text=True,
            )
            self.assertEqual(confirmed.returncode, 0, confirmed.stdout + confirmed.stderr)
            confirmation = json.loads((root / "state/final_annotation_confirmation.json").read_text())
            self.assertEqual(confirmation["confirmation_review_report_sha256"], frozen["confirmation_review_report_sha256"])
            self.assertEqual(confirmation["master_quality_approval_sha256"], frozen["master_quality_approval_sha256"])
            autopilot = (SCRIPTS / "autopilot_status.py").read_text()
            self.assertLess(autopilot.index('"master_quality_approval"'), autopilot.index('"confirmation_review"'))
            self.assertLess(autopilot.index('"confirmation_review"'), autopilot.index('"user_confirmation"'))


if __name__ == "__main__":
    unittest.main()
