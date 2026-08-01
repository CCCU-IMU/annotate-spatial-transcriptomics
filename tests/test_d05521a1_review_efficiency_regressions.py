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
HAS_RUNTIME = all(
    importlib.util.find_spec(name) for name in ("numpy", "pandas", "scipy")
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()), "sha256": sha(path),
        "n_bytes": path.stat().st_size,
    }


def write_tsv(
    path: Path, rows: list[dict[str, object]], fields: list[str] | None = None,
) -> None:
    fields = fields or list(rows[0])
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class D05521A1ReviewEfficiencyArchitectureTests(unittest.TestCase):
    def test_controller_reuses_static_counts_and_skips_noop_transform(self) -> None:
        source = (SCRIPTS / "run_lineage_controller.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("prepare_catalog_review_static_evidence", source)
        self.assertIn("catalog_wide_static_evidence", source)
        self.assertIn(
            'if int(apply_doc.get("n_changed_observations", 0)) == 0:',
            source,
        )

    def test_runtime_resource_classes_are_bound(self) -> None:
        thresholds = json.loads(THRESHOLDS.read_text(encoding="utf-8"))
        policy = thresholds["runtime_resource_policy"]
        self.assertEqual(policy["state_manifest_and_preflight_cpus"], 4)
        self.assertEqual(policy["review_evidence_default_cpus"], 8)
        self.assertEqual(policy["review_evidence_max_cpus"], 16)
        self.assertEqual(policy["heavy_clustering_default_cpus"], 64)

    def test_zero_change_decision_reuses_source_membership(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            contract = root / "contract.json"
            contract.write_text("{}\n", encoding="utf-8")
            membership = root / "membership.tsv"
            write_tsv(membership, [{
                "cell_id": "c1", "source_boundary": "cohort_1",
                "source_cluster": "1", "final_broad_label": "Granulosa",
                "final_state": "defined_broad_only", "candidate_id": "granulosa",
            }])
            signature = hashlib.sha256(b"granulosa-scope").hexdigest()
            queue = root / "queue.tsv"
            write_tsv(queue, [{
                "review_id": "review_1", "review_mode": "broad_lineage_review",
                "target_broad_label": "Granulosa", "unit_signature": signature,
            }])
            components = root / "components.tsv.gz"
            write_tsv(components, [], fields=[
                "component_id", "cell_id", "candidate_id",
            ])
            scope = root / "scope.tsv.gz"
            write_tsv(scope, [{
                "review_id": "review_1", "cell_id": "c1",
            }])
            review = root / "review.json"
            review.write_text(json.dumps({
                "review_round": 1, "membership": artifact(membership),
                "artifacts": {
                    "review_queue": artifact(queue),
                    "outside_label_recall_component_membership": artifact(components),
                    "broad_lineage_review_scope_membership": artifact(scope),
                },
            }), encoding="utf-8")
            current_questions = root / "current_questions.tsv.gz"
            recall_questions = root / "recall_questions.tsv.gz"
            write_tsv(current_questions, [], fields=[
                "broad_label", "cell_id",
            ])
            write_tsv(recall_questions, [], fields=[
                "broad_label", "cell_id", "component_id",
            ])
            broad_evidence = root / "broad_evidence.json"
            broad_evidence.write_text(json.dumps({
                "artifacts": {
                    "current_member_questions": artifact(current_questions),
                    "recall_membership": artifact(recall_questions),
                },
            }), encoding="utf-8")
            packet = root / "packet.json"
            packet.write_text(json.dumps({
                "broad_evidence_manifest": artifact(broad_evidence),
            }), encoding="utf-8")
            state = root / "state.json"
            state.write_text(json.dumps({
                "artifact_role": "sequential_cell_type_review_state",
                "active_review_n": 1,
                "formal_batch_closure_forbidden": True,
                "active_cell_type_review": {
                    "review_id": "review_1", "target_broad_label": "Granulosa",
                },
            }), encoding="utf-8")
            decisions = root / "decisions.tsv"
            write_tsv(decisions, [{
                "review_id": "review_1",
                "outcome": "retain_current_cell_type",
                "current_member_precision": "supported",
                "whole_query_recall": "adequate",
                "spatial_consistency": "supported",
                "molecular_support": "supported",
                "evidence_basis": "bound multichannel evidence",
                "rationale": "No membership correction is supported.",
            }])
            validation = root / "validation.json"
            validation.write_text(json.dumps({
                "status": "PASS", "formal_batch_closure_performed": False,
                "decision_n": 1, "active_review_id": "review_1",
                "active_cell_type": "Granulosa",
                "review_state": artifact(state),
                "evidence_packet_manifest": artifact(packet),
                "review_manifest": artifact(review),
                "validated_decisions": artifact(decisions),
            }), encoding="utf-8")
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({"candidate_boundaries": [{
                "candidate_id": "granulosa", "candidate_role": "broad",
                "release_broad_label": "Granulosa",
            }]}), encoding="utf-8")
            authority = root / "authority.json"
            authority.write_text(json.dumps({
                "mode": "stage_authority",
                "phase": "atlas_and_completeness_review",
                "annotation_contract_sha256": sha(contract),
                "post_atlas_membership": artifact(membership),
                "catalog_review_manifest": artifact(review),
                "cell_type_review_state": artifact(state),
                "catalog_decision_validation": artifact(validation),
                "candidate_catalog": artifact(catalog),
            }), encoding="utf-8")
            out = root / "apply"
            result = subprocess.run([
                sys.executable,
                str(SCRIPTS / "apply_catalog_wide_lineage_review.py"),
                "--contract", str(contract),
                "--stage-authority", str(authority),
                "--membership", str(membership),
                "--review-manifest", str(review),
                "--review-state", str(state),
                "--decision-validation", str(validation),
                "--catalog", str(catalog), "--out", str(out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            manifest = json.loads(
                (out / "catalog_wide_lineage_review_apply_manifest.json")
                .read_text(encoding="utf-8")
            )
            self.assertFalse(manifest["membership_rewritten"])
            self.assertEqual(manifest["n_changed_observations"], 0)
            self.assertEqual(
                Path(manifest["membership"]["path"]).resolve(),
                membership.resolve(),
            )
            self.assertFalse(
                (out / "catalog_wide_reviewed_membership.tsv.gz").exists()
            )

    def test_manual_adjudication_closes_exact_unit_even_after_reaudit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            membership = root / "membership.tsv"
            write_tsv(membership, [{
                "cell_id": "c1", "final_broad_label": "Granulosa",
                "source_boundary": "cohort_1", "source_cluster": "1",
            }])
            queue = root / "queue.tsv"
            write_tsv(queue, [], fields=[
                "review_id", "review_mode", "target_broad_label",
                "unit_signature",
            ])
            signature = hashlib.sha256(b"exact-scope").hexdigest()
            summary = root / "summary.tsv"
            write_tsv(summary, [{
                "broad_label": "Granulosa",
                "review_mode": "broad_lineage_review",
                "unit_signature": signature,
            }])
            review = root / "review.json"
            review.write_text(json.dumps({
                "stage": "post_atlas_catalog_wide_lineage_review",
                "membership": artifact(membership),
                "artifacts": {
                    "review_queue": artifact(queue),
                    "broad_lineage_review_summary": artifact(summary),
                },
            }), encoding="utf-8")
            blocked_review = root / "blocked_review.json"
            blocked_review.write_text("{}\n", encoding="utf-8")
            ledger = root / "blocked_ledger.tsv"
            write_tsv(ledger, [{
                "review_mode": "broad_lineage_review",
                "target_broad_label": "Granulosa",
                "unit_signature": signature,
                "status": "blocked_maximum_decisions",
                "decision_count": 2,
            }])
            blocked_state = root / "blocked_state.json"
            blocked_state.write_text(json.dumps({
                "status": "BLOCKED",
                "membership": artifact(membership),
                "review_manifest": artifact(blocked_review),
                "task_ledger": artifact(ledger),
            }), encoding="utf-8")
            support_a = root / "support_a.txt"
            support_b = root / "support_b.txt"
            support_a.write_text("molecular\n", encoding="utf-8")
            support_b.write_text("spatial\n", encoding="utf-8")
            manual = root / "manual.json"
            manual.write_text(json.dumps({
                "schema_version": "1.0", "status": "PASS",
                "artifact_role": "user_authorized_manual_biological_adjudication",
                "review_mode": "broad_lineage_review",
                "target_broad_label": "Granulosa",
                "unit_signature": signature,
                "outcome": "retain_current_cell_type",
                "membership_changed": False,
                "counts_as_automatic_decision_round": False,
                "membership": artifact(membership),
                "blocked_review_state": artifact(blocked_state),
                "blocked_review_manifest": artifact(blocked_review),
                "user_authorization": {
                    "explicitly_confirmed": True,
                    "verbatim_text": "保留当前颗粒细胞判定。",
                },
                "five_conclusions": {
                    "current_member_precision": "supported",
                    "whole_query_recall": "adequate",
                    "molecular_identity": "supported",
                    "whole_section_spatial_consistency": "supported",
                    "literature_boundary_consistency": "supported",
                },
                "supporting_artifacts": [artifact(support_a), artifact(support_b)],
            }), encoding="utf-8")
            out = root / "state"
            result = subprocess.run([
                sys.executable,
                str(SCRIPTS / "manage_cell_type_review_queue.py"),
                "--review-manifest", str(review),
                "--manual-biological-adjudication", str(manual),
                "--maximum-decisions-per-cell-type", "2",
                "--out", str(out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            state = json.loads(
                (out / "cell_type_review_state.json").read_text(encoding="utf-8")
            )
            self.assertEqual(state["status"], "COMPLETE")
            self.assertEqual(state["closed_review_n"], 1)
            ledger_rows = read_tsv(out / "cell_type_review_task_ledger.tsv")
            self.assertEqual(
                ledger_rows[0]["status"], "closed_by_manual_adjudication"
            )


@unittest.skipUnless(HAS_RUNTIME, "numpy/pandas/scipy runtime unavailable")
class D05521A1ReviewEfficiencyFunctionalTests(unittest.TestCase):
    def test_zero_census_direct_multifamily_program_cannot_disappear(self) -> None:
        import numpy as np
        from scipy.io import mmwrite
        from scipy.sparse import csc_matrix

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            count_root = root / "counts"
            count_root.mkdir()
            genes = ["PECAM1", "CDH5", "KDR", "FLT1"]
            cells = [f"c{i}" for i in range(6)]
            matrix = np.zeros((4, 6), dtype=int)
            matrix[:, :3] = 1
            matrix_path = count_root / "cell_type_review_marker_counts.mtx"
            mmwrite(matrix_path, csc_matrix(matrix))
            with matrix_path.open("rb") as source, gzip.open(
                str(matrix_path) + ".gz", "wb"
            ) as target:
                target.write(source.read())
            matrix_path.unlink()
            write_tsv(count_root / "cell_type_review_gene_map.tsv", [
                {"requested_gene": gene, "matched_feature": gene, "status": "matched"}
                for gene in genes
            ])
            write_tsv(count_root / "cell_type_review_cells.tsv", [
                {"cell_index": i + 1, "cell_id": cell}
                for i, cell in enumerate(cells)
            ])
            count_manifest = count_root / "cell_type_review_count_export_manifest.json"
            count_manifest.write_text(json.dumps({
                "artifact_role": "query_raw_count_cell_type_review_export",
                "assay_ancestry": "project_local_non_SCT_raw_counts",
                "raw_count_assay": "RNA",
            }), encoding="utf-8")
            marker = root / "marker.tsv"
            write_tsv(marker, [
                {"candidate_id": "endothelial", "broad_label": "Endothelial",
                 "evidence_role": "positive_family", "family_id": family,
                 "gene": gene}
                for family, genes_in_family in {
                    "endothelial_junction_backbone": ["PECAM1", "CDH5"],
                    "angiovascular_support": ["KDR", "FLT1"],
                }.items() for gene in genes_in_family
            ])
            membership = root / "membership.tsv"
            write_tsv(membership, [{
                "cell_id": cell, "source_boundary": "cohort_1",
                "source_cluster": "vascular_watch" if i < 3 else "stromal",
                "final_broad_label": "Stromal/mesenchymal",
            } for i, cell in enumerate(cells)])
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({"candidate_boundaries": [{
                "candidate_id": "endothelial", "candidate_role": "broad",
                "release_broad_label": "Endothelial",
                "required_positive_families": [
                    "endothelial_junction_backbone", "angiovascular_support"
                ],
                "seed_required_positive_families": [
                    "endothelial_junction_backbone"
                ],
            }]}), encoding="utf-8")
            out = root / "out"
            result = subprocess.run([
                sys.executable,
                str(SCRIPTS / "build_zero_census_direct_challengers.py"),
                "--count-export", str(count_root),
                "--marker-manifest", str(marker),
                "--membership", str(membership),
                "--catalog", str(catalog),
                "--threshold-registry", str(THRESHOLDS),
                "--out", str(out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            summary = read_tsv(
                out / "zero_census_direct_multifamily_challengers.tsv"
            )
            self.assertEqual(summary[0]["status"], "review_required")
            self.assertEqual(summary[0]["direct_multifamily_n"], "3")
            self.assertEqual(len(read_tsv(
                out / "zero_census_direct_multifamily_membership.tsv.gz"
            )), 3)

    def test_small_monotonic_subtraction_does_not_reopen_closed_type(self) -> None:
        sys.path.insert(0, str(SCRIPTS))
        from lineage_controller_lib import deterministic_membership_hash

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            contract = root / "contract.json"
            contract.write_text("{}\n", encoding="utf-8")
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({"candidate_boundaries": [
                {"candidate_id": "granulosa", "candidate_role": "broad",
                 "release_broad_label": "Granulosa", "specificity_priority": 80,
                 "required_positive_families": ["g1", "g2"]},
                {"candidate_id": "smooth", "candidate_role": "broad",
                 "release_broad_label": "Smooth muscle", "specificity_priority": 80,
                 "required_positive_families": ["s1", "s2"]},
            ]}), encoding="utf-8")
            members = []
            scores = []
            for i in range(20):
                broad = "Granulosa" if i < 10 else "Smooth muscle"
                cluster = "g" if i < 10 else "s"
                members.append({
                    "cell_id": f"c{i}", "source_boundary": "cohort_1",
                    "source_cluster": cluster, "final_broad_label": broad,
                    "final_state": "defined_broad_only", "candidate_id": "",
                })
                for candidate, label in (
                    ("granulosa", "Granulosa"), ("smooth", "Smooth muscle")
                ):
                    positive = broad == label
                    scores.append({
                        "cell_id": f"c{i}", "source_boundary": "cohort_1",
                        "source_cluster": cluster, "candidate_id": candidate,
                        "normalized_evidence": 0.8 if positive else 0.02,
                        "direct_signal": 0.7 if positive else 0.01,
                        "program_score": 0.5 if positive else 0,
                        "positive_family_count": 2 if positive else 0,
                        "family_coherent": str(positive).lower(),
                        "release_family_coherent": str(positive).lower(),
                        "identity_core_coherent": str(positive).lower(),
                        "identity_core_direct": str(positive).lower(),
                        "hard_contradiction": "false", "technical_flag": "false",
                        "x": i, "y": 0,
                    })
            membership_1 = root / "membership_1.tsv"
            score_path = root / "scores.tsv.gz"
            write_tsv(membership_1, members)
            write_tsv(score_path, scores)
            evidence = root / "evidence.tsv"
            evidence_rows = []
            for cluster, candidate in (("g", "granulosa"), ("s", "smooth")):
                evidence_rows.append({
                    "resolution_role": "selected",
                    "source_boundary": "cohort_1", "source_cluster": cluster,
                    "candidate_id": candidate, "n_observations": 10,
                    "available_positive_family_count": 2,
                    "group_positive_family_supported_count": 2,
                    "group_required_positive_families_pass": "true",
                    "observation_identity_core_fraction": 0.9,
                    "observation_identity_core_direct_fraction": 0.8,
                    "positive_marker_detection_fraction": 0.9,
                    "mean_program_score": 0.5,
                    "marker_deg_log2fc_mean": 1.5,
                    "anti_marker_deg_log2fc_mean": 0,
                    "positive_marker_pseudobulk_sum": 8,
                    "anti_marker_pseudobulk_sum": 0,
                    "cross_resolution_stable_fraction": 0.9,
                })
            write_tsv(evidence, evidence_rows)

            def authority(path: Path, membership: Path) -> None:
                path.write_text(json.dumps({
                    "mode": "stage_authority",
                    "phase": "atlas_and_completeness_review",
                    "annotation_contract_sha256": sha(contract),
                    "post_atlas_membership": artifact(membership),
                    "candidate_catalog": artifact(catalog),
                    "threshold_registry": artifact(THRESHOLDS),
                    "observation_scores": [artifact(score_path)],
                    "cluster_evidence": [artifact(evidence)],
                }), encoding="utf-8")

            auth_1 = root / "auth_1.json"
            authority(auth_1, membership_1)
            out_1 = root / "review_1"
            first = subprocess.run([
                sys.executable, str(SCRIPTS / "audit_catalog_wide_lineage_challengers.py"),
                "--contract", str(contract), "--stage-authority", str(auth_1),
                "--membership", str(membership_1), "--catalog", str(catalog),
                "--threshold-registry", str(THRESHOLDS), "--scores", str(score_path),
                "--cluster-evidence", str(evidence), "--out", str(out_1),
            ], capture_output=True, text=True)
            self.assertEqual(first.returncode, 2, first.stderr or first.stdout)
            queue_1 = read_tsv(out_1 / "catalog_wide_lineage_review_queue.tsv")
            granulosa_task = next(
                row for row in queue_1 if row["target_broad_label"] == "Granulosa"
            )
            decisions = root / "validated.tsv"
            write_tsv(decisions, [{
                "review_id": granulosa_task["review_id"],
                "outcome": "retain_current_cell_type",
            }])
            validation = root / "validation.json"
            validation.write_text(json.dumps({
                "status": "PASS",
                "review_manifest": artifact(
                    out_1 / "catalog_wide_lineage_review_manifest.json"
                ),
                "validated_decisions": artifact(decisions),
            }), encoding="utf-8")

            members_2 = [dict(row) for row in members]
            members_2[0]["final_broad_label"] = "Smooth muscle"
            membership_2 = root / "membership_2.tsv"
            write_tsv(membership_2, members_2)
            auth_2 = root / "auth_2.json"
            authority(auth_2, membership_2)
            out_2 = root / "review_2"
            second = subprocess.run([
                sys.executable, str(SCRIPTS / "audit_catalog_wide_lineage_challengers.py"),
                "--contract", str(contract), "--stage-authority", str(auth_2),
                "--membership", str(membership_2), "--catalog", str(catalog),
                "--threshold-registry", str(THRESHOLDS), "--scores", str(score_path),
                "--cluster-evidence", str(evidence),
                "--prior-decision-validation", str(validation),
                "--out", str(out_2),
            ], capture_output=True, text=True)
            self.assertEqual(second.returncode, 2, second.stderr or second.stdout)
            summary = read_tsv(out_2 / "broad_lineage_review_summary.tsv")
            granulosa = next(
                row for row in summary if row["broad_label"] == "Granulosa"
            )
            self.assertEqual(
                granulosa["status"], "closed_after_monotonic_subtraction"
            )
            self.assertEqual(granulosa["monotonic_removed_n"], "1")
            queue_2 = read_tsv(out_2 / "catalog_wide_lineage_review_queue.tsv")
            self.assertNotIn(
                "Granulosa", {row["target_broad_label"] for row in queue_2}
            )


if __name__ == "__main__":
    unittest.main()
