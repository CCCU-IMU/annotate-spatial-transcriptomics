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
CATALOG = SKILL / "references/profiles/sheep_ovary_candidate_lineage_catalog.json"


def run(*args: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class V19ContractTests(unittest.TestCase):
    def init_project(self, temp: str) -> Path:
        project = Path(temp) / "project"
        result = run(SCRIPTS / "init_annotation_project.py", "--sample", "s1", "--input-root", temp,
                     "--project-root", project, "--modality", "spatial", "--observation-unit", "cellbin")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return project

    def test_project_enables_new_fail_closed_gates(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            project = self.init_project(temp)
            config = json.loads((project / "config/project.json").read_text())
            self.assertEqual(config["framework_version"], "1.10.0")
            self.assertTrue(config["project_input_boundary_validation_required"])
            self.assertTrue(config["broad_class_completeness_review_required"])
            self.assertTrue((project / "state/derived_expression_registry.tsv").is_file())
            self.assertTrue((project / "state/broad_class_completeness_registry.tsv").is_file())

    def test_cross_project_query_expression_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            project = self.init_project(temp)
            artifact = Path(temp) / "derived.rds"
            artifact.write_bytes(b"fixture")
            result = run(SCRIPTS / "register_derived_expression_artifact.py", project,
                         "--artifact-id", "a1", "--sample", "s1", "--path", artifact,
                         "--artifact-kind", "seurat_rds", "--purpose", "marker_validation",
                         "--parent-snapshot-id", "raw", "--raw-counts-sha256", "r" * 64,
                         "--analysis-set-sha256", "a" * 64, "--source-project-id", "project_A")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("cross-project", result.stderr + result.stdout)

    def test_large_input_snapshot_accepts_matching_scheduler_hash_record(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            project = self.init_project(temp)
            raw = Path(temp) / "raw.rds"; raw.write_bytes(b"raw")
            record = Path(temp) / "raw.sha256"; record.write_text(f"{digest(raw)}  {raw.resolve()}\n")
            result = run(SCRIPTS / "register_input_snapshot.py", project, "--snapshot-id", "raw", "--sample", "s1",
                         "--path", raw, "--kind", "raw_counts", "--sha256-file", record)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_atlas_channel_can_register_cross_project_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            project = self.init_project(temp)
            raw = Path(temp) / "raw.rds"; raw.write_bytes(b"raw")
            run(SCRIPTS / "register_input_snapshot.py", project, "--snapshot-id", "raw", "--sample", "s1",
                "--path", raw, "--kind", "raw_counts")
            artifact = Path(temp) / "atlas.rds"; artifact.write_bytes(b"atlas")
            result = run(SCRIPTS / "register_derived_expression_artifact.py", project,
                         "--artifact-id", "atlas", "--sample", "s1", "--path", artifact,
                         "--artifact-kind", "reference_rds", "--purpose", "atlas", "--parent-snapshot-id", "raw",
                         "--raw-counts-sha256", "r" * 64, "--analysis-set-sha256", "a" * 64,
                         "--source-project-id", "reference_project", "--external-reference")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            validated = run(SCRIPTS / "validate_project_input_boundary.py", project)
            self.assertEqual(validated.returncode, 0, validated.stdout + validated.stderr)

    def test_zero_census_default_lineages_fail_without_query_review(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            project = self.init_project(temp)
            result = run(SCRIPTS / "validate_broad_class_completeness.py", project, "--catalog", CATALOG)
            self.assertEqual(result.returncode, 2)
            payload = json.loads(result.stdout)
            self.assertTrue(any("default broad candidates" in error for error in payload["errors"]))

    def test_single_positive_family_cannot_be_closed_as_absent(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            project = self.init_project(temp)
            catalog = Path(temp) / "catalog.json"
            catalog.write_text(json.dumps({"candidate_boundaries": [{"candidate_id": "epithelial_mesothelial"}]}))
            universe = Path(temp) / "universe.tsv"
            universe.write_text("resolution\tcluster\n0.1\t1\n0.2\t1\n0.3\t1\n")
            unexplained = Path(temp) / "unexplained.tsv"; unexplained.write_text("status\nreviewed\n")
            boundary_path = project / "state/lineage_signal_boundary_registry.tsv"
            with boundary_path.open(newline="", encoding="utf-8") as handle:
                fields = csv.DictReader(handle, delimiter="\t").fieldnames or []
            boundary = {field: "" for field in fields}
            boundary.update({"boundary_id": "b1", "sample_id": "s1", "boundary_type": "whole_tissue",
                             "source_run_id": "r1", "n_observations": "100", "analysis_fraction": "0.01",
                             "candidate_resolutions": "0.1,0.2,0.3", "selected_resolution": "0.1",
                             "audited_resolutions": "0.1,0.2,0.3", "cluster_universe_artifact": str(universe),
                             "cluster_universe_sha256": digest(universe), "candidate_catalog": str(catalog),
                             "candidate_catalog_sha256": digest(catalog), "unexplained_program_audit": "true",
                             "unexplained_program_artifact": str(unexplained), "unexplained_program_sha256": digest(unexplained),
                             "status": "audited"})
            with boundary_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader(); writer.writerow(boundary)
            signal_path = project / "state/lineage_signal_registry.tsv"
            with signal_path.open(newline="", encoding="utf-8") as handle:
                signal_fields = csv.DictReader(handle, delimiter="\t").fieldnames or []
            rows = []
            for index, resolution in enumerate(("0.1", "0.2", "0.3"), start=1):
                row = {field: "" for field in signal_fields}
                row.update({"signal_id": f"s{index}", "sample_id": "s1", "boundary_id": "b1",
                            "resolution": resolution, "cluster": "1", "candidate_lineage": "epithelial_mesothelial",
                            "candidate_source": "catalog", "signal_status": "absent", "review_status": "resolved"})
                if resolution == "0.1":
                    row.update({"positive_family_count": "1", "positive_families": "keratin",
                                "positive_genes": "KRT19", "required_action": "carry_forward"})
                rows.append(row)
            with signal_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, signal_fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
            result = run(SCRIPTS / "validate_lineage_signal_coverage.py", project)
            self.assertEqual(result.returncode, 2)
            self.assertIn("positive lineage evidence cannot be recorded as absent", result.stdout)

    def test_graph_cluster_increase_alone_is_not_rescue(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            table = Path(temp) / "graphs.tsv"
            table.write_text("graph_id\trole\tk\tn_clusters\tcore_support_separation\tanti_clearance\tspatial_coherence\n"
                             "g1\tdefault\t30\t2\t0.20\t0.30\t0.80\n"
                             "g2\tlocal\t10\t10\t0.22\t0.29\t0.80\n")
            out = Path(temp) / "result.json"
            result = run(SCRIPTS / "evaluate_graph_sensitivity.py", table, "--out", out)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            payload = json.loads(out.read_text())
            self.assertFalse(payload["rescued"])
            self.assertEqual(payload["decision"], "not_rescued_do_not_name_small_clusters")

    def test_spatially_coherent_four_percent_lineage_is_retained(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            table = Path(temp) / "scores.tsv"
            lines = ["cell_id\tparent_label\tcandidate_lineage\tprogram_score\tcore_hits\tsupport_hits\tanti_hits\tx\ty"]
            for i in range(1000):
                positive = i < 40
                lines.append(f"c{i}\tStromal\tepithelial_mesothelial\t{10 if positive else 0}\t{2 if positive else 0}\t{2 if positive else 0}\t0\t{i % 20 if positive else 1000+i}\t{i // 20 if positive else 1000}")
            table.write_text("\n".join(lines) + "\n")
            out = Path(temp) / "embedded"
            result = run(SCRIPTS / "screen_embedded_lineage_components.py", table, "--out", out,
                         "--radius", "1.5", "--min-component", "25", "--positive-quantile", "0.95")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            manifest = json.loads((out / "embedded_lineage_manifest.json").read_text())
            self.assertGreaterEqual(manifest["n_embedded_candidates"], 1)

    def test_oocyte_context_cannot_expand_oocyte_census(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            canonical = Path(temp) / "canonical.tsv"; canonical.write_text("cell_id\na\n")
            context = Path(temp) / "context.tsv"; context.write_text("cell_id\na\nb\n")
            outcomes = Path(temp) / "outcomes.tsv"; outcomes.write_text("cell_id\tfinal_broad_label\na\tOocyte\nb\tOocyte\n")
            result = run(SCRIPTS / "validate_oocyte_context_boundary.py", "--canonical-membership", canonical,
                         "--context-membership", context, "--outcomes", outcomes, "--out", Path(temp) / "audit.json")
            self.assertEqual(result.returncode, 2)
            self.assertIn("context-only", result.stdout)

    def test_banksy_cluster_count_selection_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            carry = Path(temp) / "signals.tsv"; carry.write_text("signal\nwatch\n")
            reviews = {key: True for key in (
                "full_catalog_lineage_scan", "default_broad_recall_review", "large_cluster_purity_review",
                "zero_census_review", "deg_marker_coherence_review", "spatial_morphology_review",
                "adjacent_resolution_migration_review", "technical_fragmentation_review")}
            decision = {"method": "BANKSY", "question_mode": "whole_tissue_broad_annotation",
                        "selected_resolution": "0.4", "carryforward_signal_artifact": str(carry),
                        "carryforward_signal_sha256": digest(carry),
                        "candidates": [{"resolution": "0.4", **reviews, "selection_basis": "cluster_count_only",
                                        "selection_rationale": "target number"}]}
            path = Path(temp) / "decision.json"; path.write_text(json.dumps(decision))
            result = run(SCRIPTS / "validate_banksy_broad_resolution_selection.py", path)
            self.assertEqual(result.returncode, 2)
            self.assertIn("cluster count", result.stdout)

    def test_cross_lineage_reconstruction_remains_open_world(self) -> None:
        text = (SKILL / "references/direct-lineage-controller.md").read_text(encoding="utf-8")
        self.assertIn("must not narrow the candidates", text)
        self.assertIn("Direct cross-lineage return", text)

    def test_reviewed_mapping_preserves_prelabel_freeze_fields(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            project = self.init_project(temp)
            evidence = Path(temp) / "cluster_1.json"; evidence.write_text("{}\n")
            mapping = Path(temp) / "mapping.tsv"
            mapping.write_text(
                "source_cluster\tbroad_label\tfine_label\tstate\tconfidence\tevidence_status\troute\tfine_anchor_eligible\tnext_action\tclosed\tprelabel_evidence_artifact\tprelabel_evidence_sha256\tprelabel_evidence_frozen\tprelabel_winner\tprelabel_runner_up\tprelabel_winning_margin\n"
                f"1\tGranulosa\t\tdefined\thigh\tsupported\tbroad_class_recluster\tfalse\trecluster\tfalse\t{evidence}\t{digest(evidence)}\ttrue\tGranulosa\tTheca\t0.42\n"
            )
            counts = Path(temp) / "counts.tsv"; counts.write_text("cluster\tn_observations\n1\t100\n")
            result = run(SCRIPTS / "commit_reviewed_mapping.py", project, "--mapping", mapping,
                         "--counts", counts, "--selected-run", "res0p8", "--method", "BANKSY",
                         "--sample", "s1", "--selection-rationale", "query evidence")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (project / "state/cluster_decision_ledger.tsv").open(newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["prelabel_evidence_frozen"], "true")
            self.assertEqual(row["prelabel_winner"], "Granulosa")
            self.assertEqual(row["prelabel_runner_up"], "Theca")
            self.assertEqual(row["prelabel_winning_margin"], "0.42")
            self.assertEqual(row["state"], "defined_broad_only")

    def test_state_validator_understands_frozen_full_scope_membership(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            project = self.init_project(temp)
            membership = project / "membership.tsv"
            membership.write_text(
                "sample_id\tcell_id\tanalysis_scope\n"
                "s1\ta\tanalysis_set\n"
                "s1\tb\tanalysis_set\n"
                "s1\tc\texcluded_initial_qc\n"
            )
            policy = {
                "status": "PASS", "membership_path": str(membership),
                "membership_sha256": digest(membership), "full_object_count": 3,
                "analysis_set_count": 2, "excluded_initial_qc_count": 1,
            }
            (project / "provenance/analysis_scope_policy.json").write_text(json.dumps(policy))
            result = run(SCRIPTS / "validate_state.py", project)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
