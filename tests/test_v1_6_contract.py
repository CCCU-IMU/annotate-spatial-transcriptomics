from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "annotate-spatial-transcriptomics"
SCRIPTS = SKILL / "scripts"
PROFILE = SKILL / "references/profiles/sheep_ovary.json"
GRID = [0.1, 0.2, 0.3, 0.4, 0.6]
sys.path.insert(0, str(SCRIPTS))

from audit_annotation_membership_partition import audit as audit_partition  # noqa: E402
from evaluate_evidence_ablation import main as _ablation_main  # noqa: F401,E402
from validate_benchmark_isolation import validate as validate_isolation  # noqa: E402
from validate_cohort_outcome import validate as validate_cohort  # noqa: E402
from validate_direct_return_evidence import validate as validate_return  # noqa: E402
from validate_annotation_support_registry import validate as validate_support  # noqa: E402
from dependency_manifest import build as build_dependencies, stale as dependency_stale  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    return path


def ref(path: Path, n: int | None = None) -> dict:
    value = {"path": str(path), "sha256": sha(path)}
    if n is not None:
        value["n_observations"] = n
    return value


def evidence_files(root: Path) -> tuple[Path, Path]:
    full = root / "evidence/full_feature.json"
    full.parent.mkdir(parents=True, exist_ok=True)
    full.write_text(json.dumps({"feature_count": 22000, "programs": ["identity", "anti-lineage"]}))
    independent = root / "evidence/biological_evidence.json"
    independent.write_text(json.dumps({"marker_families": ["FSHR-INHA", "GJA1"], "spatial_objects": 3}))
    return full, independent


def make_cohort_outcome(
    root: Path,
    query_ids: list[str],
    subclusters: list[tuple[str, str, list[str], str, str]],
    *,
    cohort_id: str = "broad_granulosa_v001",
    cohort_type: str = "broad_class_recluster",
    question_mode: str = "broad_purity_audit",
    homogeneous: bool = False,
) -> Path:
    query = write_tsv(root / f"memberships/{cohort_id}_query.tsv", [{"cell_id": value} for value in query_ids])
    full, _ = evidence_files(root)
    resolutions = []
    for value in GRID:
        token = str(value).replace(".", "p")
        membership = write_tsv(
            root / f"memberships/{cohort_id}_res{token}.tsv",
            [{"cell_id": cell_id, "cluster": "0"} for cell_id in query_ids],
        )
        index = root / f"evidence/{cohort_id}_res{token}.json"
        index.write_text(json.dumps({"resolution": value, "marker_programs": ["granulosa"], "anti_marker_hits": 0}))
        resolutions.append({"resolution": value, "membership": ref(membership, len(query_ids)), "cluster_count": 1, "evidence_index": ref(index)})
    outcomes = []
    for subcluster_id, outcome, ids, broad, fine in subclusters:
        membership = write_tsv(root / f"memberships/{cohort_id}_{subcluster_id}.tsv", [{"cell_id": value} for value in ids])
        row = {"subcluster_id": subcluster_id, "outcome": outcome, "membership": ref(membership, len(ids)), "target_broad_label": broad, "target_fine_label": fine}
        outcomes.append(row)
    document = {
        "schema_version": "1.0", "cohort_id": cohort_id, "cohort_type": cohort_type,
        "question_mode": question_mode, "query_membership": ref(query, len(query_ids)),
        "candidate_grid": GRID, "resolutions": resolutions, "selected_resolution": 0.3,
        "marker_evidence": {
            "positive_marker_families": ["FSHR-INHA", "GJA1-FST"],
            "anti_marker_results": {"oocyte_core": "absent", "immune": "absent"},
            "full_feature_scope": {"status": "PASS", "feature_space": "full_feature", "artifact": ref(full)},
        },
        "adjacent_stability": [
            {"left": left, "right": right, "metric": "ARI", "value": 1.0}
            for left, right in zip(GRID, GRID[1:])
        ],
        "spatial_morphology": {"pattern": "coherent follicular layer", "object_count": 3},
        "source_qc_composition": {"source_fraction": 1.0, "median_features": 820},
        "selected_resolution_rationale": "Resolution 0.3 best preserves the coherent lineage program without technical fragmentation.",
        "rejected_alternatives": [
            {"resolution": value, "reason": "Equivalent or weaker biological separation and stability evidence."}
            for value in GRID if value != 0.3
        ],
        "subcluster_outcomes": outcomes,
        "terminal_outcome": "homogeneous_parent_confirmed" if homogeneous else "subclusters_adjudicated",
        "homogeneous_parent_confirmed": homogeneous,
    }
    path = root / f"evidence/{cohort_id}_outcome.json"
    path.write_text(json.dumps(document, indent=2))
    return path


def make_return_evidence(root: Path, return_id: str, cohort_id: str, subcluster: str, membership: Path, broad: str, fine: str = "", confidence: str = "high") -> Path:
    full, independent = evidence_files(root)
    with membership.open(newline="", encoding="utf-8") as handle:
        n = len(list(csv.DictReader(handle, delimiter="\t")))
    document = {
        "schema_version": "1.0", "return_id": return_id, "source_cohort_id": cohort_id,
        "source_subcluster": subcluster, "target_broad_label": broad, "target_fine_label": fine,
        "membership": ref(membership, n), "positive_marker_families": ["FSHR-INHA", "GJA1"],
        "anti_marker_results": {"competing_lineages": "not_supported"},
        "spatial_evidence": {"pattern": "coherent follicular layer", "object_count": 3},
        "full_feature_scope": {"status": "PASS", "feature_space": "full_feature", "artifact": ref(full)},
        "confidence": confidence, "alternative_hypotheses": ["stromal contamination rejected"],
        "evidence_artifacts": [ref(independent)],
    }
    path = root / f"evidence/{return_id}.json"; path.write_text(json.dumps(document, indent=2)); return path


class V160ActiveContract(unittest.TestCase):
    def test_expensive_asset_staleness_uses_content_hash_not_mtime(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); dependency = root / "input.tsv"; target = root / "report.html"
            dependency.write_text("a\n"); target.write_text("report\n")
            build_dependencies(target, [dependency])
            original = dependency.stat().st_mtime_ns
            os.utime(dependency, ns=(original + 1_000_000, original + 1_000_000))
            self.assertFalse(dependency_stale(target, [dependency]))
            dependency.write_text("b\n"); os.utime(dependency, ns=(original, original))
            self.assertTrue(dependency_stale(target, [dependency]))

    def test_full_feature_gate_reads_loaded_profile_not_context_hint(self) -> None:
        text = (SCRIPTS / "check_completion_gate.py").read_text(encoding="utf-8")
        self.assertIn('requires_full_feature=bool(profile.get("final_validation"))', text)
        self.assertNotIn('context.get("profile")', text)

    def test_new_rctd_evidence_uses_canonical_confidence_tiers(self) -> None:
        runner = (SCRIPTS / "run_calibrated_rctd_review.R").read_text(encoding="utf-8")
        registry = (SCRIPTS / "register_route_attempt.py").read_text(encoding="utf-8")
        for token in ['levels=c("high","moderate","low")', "rctd_moderate_n", "rctd_low_n"]:
            self.assertIn(token, runner + registry)
        self.assertNotIn('rctd_confidence_tier:="medium_low"', runner)

    def test_fixed_grid_requires_verified_active_workflow_and_preset_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); (root / "config").mkdir(); (root / "provenance").mkdir()
            (root / "config/project.json").write_text(json.dumps({"strategy_preset_requested": "sheep_ovary_same_batch_rfirst"}))
            context = root / "context.json"; context.write_text(json.dumps({"species": "Ovis aries", "tissue": "ovary"}))
            inspection = root / "inspection.json"; inspection.write_text(json.dumps({"type": "Seurat", "assays": ["Spatial"], "n_features": 22000}))
            hints_only = subprocess.run([
                sys.executable, str(SCRIPTS / "resolve_workflow_profile.py"), "--context", str(context), "--r-inspection", str(inspection),
                "--input-path", "/portable/cellbin_PPed/example.rds", "--out", str(root / "hints_only.json"),
            ], capture_output=True, text=True)
            self.assertEqual(hints_only.returncode, 0, hints_only.stdout + hints_only.stderr)
            self.assertFalse(json.loads((root / "hints_only.json").read_text())["fixed_cellbin_preprocessing_required"])
            manifests = {
                "conversion.json": {"status": "PASS", "source_platform": "StereoPy", "cellbin_pped_conversion": True, "same_batch": True},
                "layers.json": {"status": "PASS", "raw_counts_verified": True},
                "coordinates.json": {"status": "PASS", "cell_coordinate_ids_match": True},
                "markers.json": {"status": "PASS", "coverage_sufficient": True},
            }
            for name, value in manifests.items(): (root / name).write_text(json.dumps(value))
            verified_without_preset = subprocess.run([
                sys.executable, str(SCRIPTS / "resolve_workflow_profile.py"), "--context", str(context), "--r-inspection", str(inspection),
                "--conversion-manifest", str(root / "conversion.json"), "--layer-audit", str(root / "layers.json"),
                "--coordinate-audit", str(root / "coordinates.json"), "--marker-coverage-audit", str(root / "markers.json"),
                "--out", str(root / "verified_without_preset.json"),
            ], capture_output=True, text=True)
            self.assertEqual(verified_without_preset.returncode, 0, verified_without_preset.stdout + verified_without_preset.stderr)
            self.assertFalse(json.loads((root / "verified_without_preset.json").read_text())["fixed_cellbin_preprocessing_required"])
            resolved = subprocess.run([
                sys.executable, str(SCRIPTS / "resolve_workflow_profile.py"), "--context", str(context), "--r-inspection", str(inspection),
                "--input-path", "/portable/example.rds", "--conversion-manifest", str(root / "conversion.json"),
                "--layer-audit", str(root / "layers.json"), "--coordinate-audit", str(root / "coordinates.json"),
                "--marker-coverage-audit", str(root / "markers.json"), "--strategy-preset", "sheep_ovary_same_batch_rfirst",
                "--out", str(root / "config/active_strategy_preset.json"), "--workflow-record-out", str(root / "config/active_workflow_profile.json"),
                "--input-contract-out", str(root / "provenance/input_contract_validation.json"),
            ], capture_output=True, text=True)
            self.assertEqual(resolved.returncode, 0, resolved.stdout + resolved.stderr)
            self.assertEqual(
                json.loads((root / "config/active_strategy_preset.json").read_text())["preferred_backbone"],
                "seurat_sct_banksy_whole_tissue_r_evidence",
            )
            grid = subprocess.run([
                sys.executable, str(SCRIPTS / "validate_resolution_grid.py"), "--project-root", str(root),
                "--workflow-profile", str(SKILL / "references/profiles/sheep_ovary_rfirst_profile.json"),
                "--scope", "cohort", "--resolutions", "0.1,0.2,0.3,0.4,0.6",
            ], capture_output=True, text=True)
            self.assertEqual(grid.returncode, 0, grid.stdout + grid.stderr)
            whole_grid = subprocess.run([
                sys.executable, str(SCRIPTS / "validate_resolution_grid.py"), "--project-root", str(root),
                "--workflow-profile", str(SKILL / "references/profiles/sheep_ovary_rfirst_profile.json"),
                "--scope", "whole_tissue", "--resolutions", "0.2,0.4,0.6,0.8",
                "--out", str(root / "provenance/whole_tissue_resolution_grid_validation.json"),
            ], capture_output=True, text=True)
            self.assertEqual(whole_grid.returncode, 0, whole_grid.stdout + whole_grid.stderr)

    def test_homogeneous_parent_is_success_after_complete_grid(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outcome = make_cohort_outcome(root, ["c1", "c2"], [("0", "parent_return", ["c1", "c2"], "Granulosa", "")], homogeneous=True)
            result = validate_cohort(root, outcome)
            self.assertEqual(result["status"], "PASS", result["errors"])
            self.assertEqual(result["terminal_outcome"], "homogeneous_parent_confirmed")

    def test_true_subtype_and_mixed_lineage_have_content_validated_returns(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outcome = make_cohort_outcome(
                root, ["c1", "c2"],
                [("0", "fine_return", ["c1"], "Granulosa", "Mural granulosa"), ("1", "cross_lineage_return", ["c2"], "Stromal/mesenchymal", "")],
            )
            self.assertEqual(validate_cohort(root, outcome)["status"], "PASS")
            fine_membership = root / "memberships/broad_granulosa_v001_0.tsv"
            evidence = make_return_evidence(root, "return_fine", "broad_granulosa_v001", "0", fine_membership, "Granulosa", "Mural granulosa")
            self.assertEqual(validate_return(root, evidence)["status"], "PASS")

    def test_state_only_split_can_merge_every_subcluster_to_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outcome = make_cohort_outcome(
                root, ["c1", "c2"],
                [("ecm_rich", "parent_return", ["c1"], "Stromal/mesenchymal", ""), ("stress_ribosomal", "parent_return", ["c2"], "Stromal/mesenchymal", "")],
            )
            result = validate_cohort(root, outcome)
            self.assertEqual(result["status"], "PASS", result["errors"])

    def test_empty_or_status_only_evidence_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outcome = make_cohort_outcome(root, ["c1"], [("0", "parent_return", ["c1"], "Granulosa", "")], homogeneous=True)
            document = json.loads(outcome.read_text())
            evidence_path = Path(document["resolutions"][0]["evidence_index"]["path"])
            evidence_path.write_text('{"status":"PASS"}')
            document["resolutions"][0]["evidence_index"]["sha256"] = sha(evidence_path)
            outcome.write_text(json.dumps(document))
            self.assertEqual(validate_cohort(root, outcome)["status"], "FAIL")

    def test_support_registry_requires_content_and_exact_label_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); (root / "state").mkdir(); full, independent = evidence_files(root)
            ledger = write_tsv(root / "state/cell_ledger.tsv", [{
                "cell_id": "c1", "analysis_scope": "analysis_set", "final_broad_label": "Granulosa",
                "final_fine_label": "", "final_confidence": "moderate", "fine_anchor_eligible": "false",
            }])
            support = {
                "schema_version": "1.0", "support_id": "support_granulosa", "label_level": "broad",
                "broad_label": "Granulosa", "fine_label": "", "confidence": "moderate",
                "positive_marker_evidence": {"families": ["FSHR-INHA", "GJA1-FST"]},
                "anti_marker_evidence": {"oocyte_core": "not_supported"},
                "full_feature_evidence": {"status": "PASS", "feature_space": "full_feature", "artifact": ref(full)},
                "resolution_evidence": {"selected": 0.3, "adjacent_ari": 0.9},
                "spatial_evidence": {"pattern": "follicular layer", "object_count": 3},
                "provenance_evidence": {"cohort_ids": ["broad_granulosa_v001"]},
                "alternative_hypotheses": ["stromal lineage rejected"], "evidence_artifacts": [ref(independent)],
            }
            artifact = root / "evidence/support_granulosa.json"; artifact.write_text(json.dumps(support))
            registry = write_tsv(root / "state/annotation_support_registry.tsv", [{
                "support_id": "support_granulosa", "label_level": "broad", "broad_label": "Granulosa",
                "fine_label": "", "confidence": "moderate", "support_artifact": str(artifact),
                "support_artifact_sha256": sha(artifact), "status": "validated",
            }])
            self.assertEqual(validate_support(root, registry, ledger)["status"], "PASS")
            support["spatial_evidence"] = {"status": "PASS"}; artifact.write_text(json.dumps(support))
            rows = list(csv.DictReader(registry.read_text().splitlines(), delimiter="\t")); rows[0]["support_artifact_sha256"] = sha(artifact)
            write_tsv(registry, rows, list(rows[0]))
            self.assertEqual(validate_support(root, registry, ledger)["status"], "FAIL")

    def test_partition_rejects_overlap_omission_and_atlas_query_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self._make_partition_project(root)
            self.assertEqual(audit_partition(root)["status"], "PASS", audit_partition(root)["errors"])

            route_path = root / "state/route_attempt_registry.tsv"
            original_route = route_path.read_text()
            wrong = write_tsv(root / "memberships/wrong_atlas.tsv", [{"cell_id": "c3"}])
            rows = list(csv.DictReader(original_route.splitlines(), delimiter="\t"))
            rows[0]["query_membership_artifact"] = str(wrong); rows[0]["query_membership_sha256"] = sha(wrong); rows[0]["n_query"] = "1"
            write_tsv(route_path, rows, list(rows[0]))
            self.assertEqual(audit_partition(root)["status"], "FAIL")
            route_path.write_text(original_route)

            direct_path = root / "state/direct_return_registry.tsv"
            with direct_path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            duplicate = dict(rows[0]); duplicate["return_id"] = "duplicate"
            write_tsv(direct_path, rows + [duplicate], list(rows[0]))
            self.assertEqual(audit_partition(root)["status"], "FAIL")
            write_tsv(direct_path, rows, list(rows[0]))

            cohort_path = root / "evidence/broad_granulosa_v001_outcome.json"
            document = json.loads(cohort_path.read_text()); document["subcluster_outcomes"][0]["membership"] = ref(write_tsv(root / "memberships/omitted.tsv", [{"cell_id": "c1"}]), 1)
            cohort_path.write_text(json.dumps(document)); self.assertEqual(audit_partition(root)["status"], "FAIL")

    def _make_partition_project(self, root: Path) -> None:
        (root / "state").mkdir(parents=True); (root / "config").mkdir(); (root / "provenance").mkdir()
        cells = [
            {"cell_id": "c1", "analysis_scope": "analysis_set", "initial_broad_label": "Granulosa", "final_broad_label": "Granulosa", "final_fine_label": "", "final_state": "defined_broad_only", "state": "defined_broad_only", "assignment_mode": "parent_broad_direct", "final_confidence": "high", "fine_anchor_eligible": "false"},
            {"cell_id": "c2", "analysis_scope": "analysis_set", "initial_broad_label": "Granulosa", "final_broad_label": "Granulosa", "final_fine_label": "", "final_state": "defined_broad_only", "state": "defined_broad_only", "assignment_mode": "parent_broad_direct", "final_confidence": "high", "fine_anchor_eligible": "false"},
            {"cell_id": "c3", "analysis_scope": "analysis_set", "initial_broad_label": "", "final_broad_label": "Stromal/mesenchymal", "final_fine_label": "", "final_state": "defined_broad_only", "state": "qc_holdout", "assignment_mode": "atlas_broad_rescue", "final_confidence": "moderate", "fine_anchor_eligible": "false"},
            {"cell_id": "c4", "analysis_scope": "analysis_set", "initial_broad_label": "", "final_broad_label": "", "final_fine_label": "", "final_state": "qc_reject", "state": "qc_holdout", "assignment_mode": "terminal_qc_reject", "final_confidence": "low", "fine_anchor_eligible": "false"},
        ]
        write_tsv(root / "state/cell_ledger.tsv", cells)
        outcome = make_cohort_outcome(root, ["c1", "c2"], [("0", "parent_return", ["c1", "c2"], "Granulosa", "")], homogeneous=True)
        query = root / "memberships/broad_granulosa_v001_query.tsv"
        cohort = [{"cohort_id": "broad_granulosa_v001", "cohort_type": "broad_class_recluster", "question_mode": "broad_purity_audit", "source_broad_label": "Granulosa", "membership_path": str(query), "membership_sha256": sha(query), "n_observations": 2, "candidate_resolutions": "0.1,0.2,0.3,0.4,0.6", "selected_resolution": "0.3", "terminal_outcome": "homogeneous_parent_confirmed", "outcome_artifact": str(outcome), "outcome_artifact_sha256": sha(outcome), "status": "validated_done"}]
        write_tsv(root / "state/recluster_cohort_registry.tsv", cohort)
        membership = root / "memberships/broad_granulosa_v001_0.tsv"
        evidence = make_return_evidence(root, "return_parent", "broad_granulosa_v001", "0", membership, "Granulosa")
        direct = [{"return_id": "return_parent", "source_cohort_id": "broad_granulosa_v001", "source_cluster": "0", "membership_path": str(membership), "membership_sha256": sha(membership), "n_observations": 2, "target_broad_label": "Granulosa", "target_fine_label": "", "confidence": "high", "assignment_mode": "parent_broad_direct", "evidence_artifact": str(evidence), "evidence_artifact_sha256": sha(evidence), "fine_anchor_eligible": "false", "status": "validated_done"}]
        write_tsv(root / "state/direct_return_registry.tsv", direct)
        atlas_query = write_tsv(root / "memberships/atlas_query.tsv", [{"cell_id": "c3"}, {"cell_id": "c4"}])
        accepted = write_tsv(root / "memberships/atlas_accepted.tsv", [{"cell_id": "c3"}]); rejected = write_tsv(root / "memberships/atlas_rejected.tsv", [{"cell_id": "c4"}])
        route = [{"route_attempt_id": "atlas_v001", "route_class": "residual_qc_atlas_review", "status": "validated_done", "query_membership_artifact": str(atlas_query), "query_membership_sha256": sha(atlas_query), "n_query": 2, "accepted_membership_artifact": str(accepted), "accepted_membership_sha256": sha(accepted), "accepted_membership_n_observations": 1, "rejected_membership_artifact": str(rejected), "rejected_membership_sha256": sha(rejected), "rejected_membership_n_observations": 1}]
        write_tsv(root / "state/route_attempt_registry.tsv", route)

    @unittest.skipUnless(importlib.util.find_spec("numpy") and importlib.util.find_spec("pandas") and importlib.util.find_spec("sklearn"), "ranking dependencies unavailable")
    def test_ranker_does_not_penalize_one_cluster_broad_cohort_and_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); data = root / "data"; data.mkdir()
            for value in GRID:
                token = str(value).replace(".", "p")
                write_tsv(data / f"clusters_res{token}.tsv", [{"cell_id": "c1", "cluster": "0"}, {"cell_id": "c2", "cluster": "0"}])
                write_tsv(data / f"deg_res{token}.tsv", [{"cluster": "0", "gene": gene} for gene in ["FSHR", "INHA", "GJA1", "FST"]])
            outputs = []
            for name in ("out1", "out2"):
                out = root / name
                result = subprocess.run([
                    sys.executable, str(SCRIPTS / "rank_cohort_resolutions.py"),
                    "--cluster-glob", str(data / "clusters_*.tsv"), "--deg-glob", str(data / "deg_*.tsv"),
                    "--profile", str(PROFILE), "--lineages", "granulosa", "--out", str(out),
                    "--question-mode", "broad_purity_audit", "--expected-grid", "0.1,0.2,0.3,0.4,0.6",
                ], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr); outputs.append(out)
            summary = json.loads((outputs[0] / "cohort_resolution_ranking.json").read_text())
            self.assertFalse(summary["under_split_penalty_applied"]); self.assertTrue(summary["homogeneous_parent_confirmed_is_success"])
            self.assertEqual((outputs[0] / "cohort_resolution_ranking.tsv").read_bytes(), (outputs[1] / "cohort_resolution_ranking.tsv").read_bytes())

    def test_benchmark_leakage_metrics_and_ablation(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); (root / "state").mkdir(); (root / "config").mkdir()
            config = root / "benchmark.json"; config.write_text((SKILL / "assets/benchmark_config_template.json").read_text())
            (root / "state/route_attempt_registry.tsv").write_text("reference_id\nGSE233801\n")
            self.assertEqual(validate_isolation(root, config)["status"], "FAIL")
            (root / "state/route_attempt_registry.tsv").write_text("reference_id\nindependent_reference\n")
            frozen = write_tsv(root / "frozen.tsv", [{"cell_id": f"c{i}", "sample_id": "s1" if i < 3 else "s2", "final_broad_label": label} for i, label in enumerate(["A", "A", "B", "B", "B", "QC"])])
            author = write_tsv(root / "author.tsv", [{"cell_id": f"c{i}", "sample_id": "s1" if i < 3 else "s2", "author_broad_label": label} for i, label in enumerate(["A", "A", "B", "B", "A", "B"])])
            benchmark = subprocess.run([sys.executable, str(SCRIPTS / "run_biological_benchmark.py"), str(root), "--config", str(config), "--frozen-annotation", str(frozen), "--author-labels", str(author), "--rare-labels", "Oocyte", "--out", str(root / "benchmark_out.json")], capture_output=True, text=True)
            self.assertEqual(benchmark.returncode, 0, benchmark.stdout + benchmark.stderr)
            result = json.loads((root / "benchmark_out.json").read_text()); self.assertIn("macro_f1", result); self.assertIn("cross_sample_stability", result); self.assertIn("Oocyte", result["rare_type_false_positives"])
            fixture = root / "ablation.json"; fixture.write_text(json.dumps({"truth": ["A", "A", "B", "B"], "predictions": {"full_evidence": ["A", "A", "B", "B"], "without_anti_marker": ["A", "B", "B", "B"], "without_spatial": ["A", "A", "A", "B"], "without_broad_cohort": ["A", "B", "A", "B"]}}))
            ablation = subprocess.run([sys.executable, str(SCRIPTS / "evaluate_evidence_ablation.py"), str(fixture)], capture_output=True, text=True)
            self.assertEqual(ablation.returncode, 0, ablation.stdout + ablation.stderr)


if __name__ == "__main__":
    unittest.main()
