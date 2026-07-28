from __future__ import annotations

import csv
import gzip
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
sys.path.insert(0, str(SCRIPTS))

from lineage_controller_lib import apply_candidate_context  # noqa: E402


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def catalog() -> dict:
    return {
        "candidate_boundaries": [
            {
                "candidate_id": "luteal_steroidogenic",
                "candidate_role": "broad",
                "release_broad_label": "Luteal",
                "release_fine_label": "",
                "formal_context_evidence_required": True,
                "required_positive_families": ["luteal_core", "cl_identity"],
                "writeback_strategy": "context_gated_candidate_local",
            },
            {
                "candidate_id": "stromal_mesenchymal",
                "candidate_role": "broad",
                "release_broad_label": "Stromal/mesenchymal",
                "release_fine_label": "",
                "required_positive_families": ["stromal_a", "stromal_b"],
                "writeback_strategy": "generic_exact_remainder_after_specific_lineages",
            },
        ],
        "machine_actionable_fine_candidate_catalog": {
            "luteal": [
                {
                    "candidate_id": "late_luteal",
                    "release_label": "Late luteal",
                    "parent_release_label": "Luteal",
                    "profile_program": "lineages.luteal.late",
                }
            ]
        },
    }


def evidence(candidate_id: str, positive: bool) -> dict[str, object]:
    return {
        "resolution_role": "selected",
        "source_boundary": "cohort_1",
        "source_cluster": "subcluster_1",
        "candidate_id": candidate_id,
        "n_observations": 1,
        "available_positive_family_count": 2,
        "group_positive_family_supported_count": 2 if positive else 0,
        "group_required_positive_families_pass": str(positive).lower(),
        "observation_identity_core_fraction": 0.5 if positive else 0,
        "positive_marker_detection_fraction": 0.8 if positive else 0,
        "mean_program_score": 0.4 if positive else 0,
        "marker_deg_log2fc_mean": 1.5 if positive else 0,
        "anti_marker_deg_log2fc_mean": 0,
    }


class ContextGatedCandidateIntegrityTests(unittest.TestCase):
    def test_invalid_or_foreign_context_rows_fail_closed(self) -> None:
        candidates = {
            "luteal_steroidogenic": {
                "candidate_id": "luteal_steroidogenic",
                "candidate_role": "broad",
                "release_broad_label": "Luteal",
                "formal_context_evidence_required": True,
            }
        }
        with self.assertRaises(ValueError):
            apply_candidate_context(candidates, [{
                "candidate_id": "luteal_steroidogenic",
                "status": "maybe",
            }])
        with self.assertRaises(ValueError):
            apply_candidate_context(candidates, [{
                "candidate_id": "not_in_catalog",
                "status": "supported",
            }])

    def run_completeness(
        self, root: Path, *, context_status: str, final_broad: str
    ) -> subprocess.CompletedProcess[str]:
        catalog_path = root / "catalog.json"
        catalog_path.write_text(json.dumps(catalog()), encoding="utf-8")
        context_path = root / "context.tsv"
        write_tsv(context_path, [{
            "candidate_id": "luteal_steroidogenic",
            "status": context_status,
            "reason": "fixture stage permission",
        }])
        membership = root / "membership.tsv"
        write_tsv(membership, [{
            "cell_id": "cell_1",
            "final_broad_label": final_broad,
            "final_state": "defined_broad_only",
            "assignment_origin": "post_merge_atlas_unlabeled_broad_rescue",
            "candidate_id": "",
            "source_boundary": "cohort_1",
            "source_cluster": "subcluster_1",
        }])
        cluster_evidence = root / "cluster_evidence.tsv"
        write_tsv(cluster_evidence, [
            evidence("luteal_steroidogenic", True),
            evidence("late_luteal", True),
            evidence("stromal_mesenchymal", False),
        ])
        return subprocess.run([
            sys.executable, str(SCRIPTS / "audit_post_merge_completeness.py"),
            "--membership", str(membership),
            "--catalog", str(catalog_path),
            "--context-evidence", str(context_path),
            "--cluster-evidence", str(cluster_evidence),
            "--out", str(root / "audit"),
        ], capture_output=True, text=True)

    def test_not_evaluable_luteal_program_does_not_block_zero_census(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self.run_completeness(
                root,
                context_status="not_evaluable",
                final_broad="Stromal/mesenchymal",
            )
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            rows = read_tsv(root / "audit/broad_completeness_audit.tsv")
            luteal = next(row for row in rows if row["broad_label"] == "Luteal")
            self.assertEqual(luteal["status"], "not_evaluable")
            self.assertEqual(luteal["positive_program_n"], "0")
            self.assertIn("luteal_steroidogenic", luteal["context_not_evaluable_candidate_ids"])
            self.assertIn("late_luteal", luteal["context_not_evaluable_candidate_ids"])

    def test_context_ineligible_luteal_membership_is_blocked(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self.run_completeness(
                root, context_status="not_evaluable", final_broad="Luteal"
            )
            self.assertEqual(result.returncode, 2)
            manifest = json.loads(
                (root / "audit/post_merge_completeness_manifest.json").read_text()
            )
            self.assertEqual(manifest["status"], "BLOCKED")
            self.assertTrue(any(
                "context-ineligible" in error for error in manifest["errors"]
            ))

    def test_context_ineligible_supported_fine_proposal_is_not_materialized(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            catalog_path = root / "catalog.json"
            catalog_path.write_text(json.dumps(catalog()), encoding="utf-8")
            context_path = root / "context.tsv"
            write_tsv(context_path, [{
                "candidate_id": "luteal_steroidogenic",
                "status": "not_evaluable",
                "reason": "fixture stage permission",
            }])
            membership = root / "membership.tsv"
            write_tsv(membership, [{
                "cell_id": "cell_1", "source_boundary": "cohort_1",
                "source_cluster": "subcluster_1", "final_broad_label": "Luteal",
            }])
            fine_audit = root / "fine.tsv"
            write_tsv(fine_audit, [{
                "status": "supported", "release_candidate": "true",
                "candidate_id": "late_luteal", "parent_broad_label": "Luteal",
                "lineage_supported_fraction": 0.9,
                "strongest_competing_fraction": 0.1,
                "contradiction_fraction": 0,
                "cohort_id": "cohort_1", "subcluster_id": "subcluster_1",
            }])
            result = subprocess.run([
                sys.executable,
                str(SCRIPTS / "materialize_parent_locked_fine_proposals.py"),
                "--membership", str(membership), "--catalog", str(catalog_path),
                "--context-evidence", str(context_path),
                "--threshold-registry", str(THRESHOLDS),
                "--fine-audit", str(fine_audit), "--out", str(root / "fine_out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            self.assertEqual(
                read_tsv(root / "fine_out/parent_locked_fine_assignments.tsv.gz"), []
            )
            manifest = json.loads(
                (root / "fine_out/parent_locked_fine_manifest.json").read_text()
            )
            self.assertEqual(
                manifest["n_context_ineligible_supported_sources_skipped"], 1
            )

    def test_atlas_writeback_rejects_context_ineligible_broad(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            catalog_path = root / "catalog.json"
            catalog_path.write_text(json.dumps(catalog()), encoding="utf-8")
            context_path = root / "context.tsv"
            write_tsv(context_path, [{
                "candidate_id": "luteal_steroidogenic",
                "status": "not_evaluable", "reason": "fixture",
            }])
            frozen = root / "frozen.tsv"
            write_tsv(frozen, [{
                "cell_id": "cell_1", "final_broad_label": "",
                "final_state": "unresolved_biological",
            }])
            routing = root / "routing.tsv"
            write_tsv(routing, [{
                "cell_id": "cell_1", "proposed_broad_label": "Luteal",
                "atlas_state_route": "direct_unlabeled_broad_return",
                "atlas_broad": "Luteal", "atlas_tier": "high", "review_id": "",
            }])
            validation = root / "validation.json"
            validation.write_text(json.dumps({"status": "PASS"}), encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "apply_post_merge_atlas_routing.py"),
                "--frozen-broad", str(frozen), "--routing", str(routing),
                "--atlas-validation", str(validation), "--catalog", str(catalog_path),
                "--context-evidence", str(context_path), "--out", str(root / "atlas"),
            ], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("context-ineligible", result.stderr + result.stdout)

    def test_every_post_scoring_release_outlet_is_context_aware(self) -> None:
        required = {
            "build_resolution_grid_evidence.py": "apply_candidate_context",
            "adjudicate_second_round_subclusters.py": "apply_candidate_context",
            "close_exact_remainders.py": "apply_candidate_context",
            "merge_and_freeze_broad_membership.py": "apply_candidate_context",
            "review_post_merge_unresolved_components.py": "apply_candidate_context",
            "audit_post_merge_completeness.py": "apply_candidate_context",
            "route_global_atlas_v2.py": "apply_candidate_context",
            "apply_post_merge_atlas_routing.py": "apply_candidate_context",
            "audit_catalog_wide_lineage_challengers.py": "apply_candidate_context",
            "apply_catalog_wide_lineage_review.py": "apply_candidate_context",
            "apply_sheep_ovary_follicle_roi_repair.py": "apply_candidate_context",
            "materialize_oocyte_cluster_membership.py": "apply_candidate_context",
            "apply_cell_id_membership_patch.py": "apply_candidate_context",
            "materialize_parent_locked_fine_proposals.py": "apply_candidate_context",
            "materialize_final_release_v2_2.py": "apply_candidate_context",
        }
        for name, token in required.items():
            source = (SCRIPTS / name).read_text(encoding="utf-8")
            self.assertIn(token, source, name)
            self.assertIn("context-evidence", source, name)


if __name__ == "__main__":
    unittest.main()
