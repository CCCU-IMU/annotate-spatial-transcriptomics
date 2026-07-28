from __future__ import annotations

import json
import hashlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT / "annotate-spatial-transcriptomics" / "scripts"
    / "validate_sheep_ovary_biological_quality.py"
)


class SheepOvaryBiologicalQualityReviewTest(unittest.TestCase):
    candidates = {
        "granulosa": "Granulosa",
        "oocyte": "Oocyte",
        "theca_steroidogenic": "Theca",
        "vascular_endothelial": "Vascular-associated",
        "pericyte_mural": "Vascular-associated",
        "smooth_muscle": "Smooth muscle",
        "stromal_mesenchymal": "Stromal/mesenchymal",
    }

    @staticmethod
    def ring(prefix: str, radius: float, n: int, source_cluster: str):
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
        return [
            {
                "cell_id": f"{prefix}_{index:03d}",
                "x": radius * np.cos(angle),
                "y": radius * np.sin(angle),
                "source_cluster": source_cluster,
            }
            for index, angle in enumerate(angles)
        ]

    def build_fixture(self, root: Path, under_recall_theca: bool = False):
        cells = []
        identity = {}
        broad = {}
        for candidate, label, radius, n in (
            ("granulosa", "Granulosa", 20.0, 96),
            ("theca_steroidogenic", "Theca", 23.0, 96),
            ("vascular_endothelial", "Vascular-associated", 25.0, 48),
            ("smooth_muscle", "Smooth muscle", 28.0, 96),
            ("stromal_mesenchymal", "Stromal/mesenchymal", 31.0, 96),
        ):
            rows = self.ring(candidate, radius, n, candidate)
            cells.extend(rows)
            for row in rows:
                identity[row["cell_id"]] = candidate
                broad[row["cell_id"]] = label
        oocytes = self.ring("oocyte", 1.5, 12, "canonical_oocyte")
        cells.extend(oocytes)
        for row in oocytes:
            identity[row["cell_id"]] = "oocyte"
            broad[row["cell_id"]] = "Oocyte"
        if under_recall_theca:
            for cell_id, candidate in identity.items():
                if candidate == "theca_steroidogenic":
                    broad[cell_id] = "Stromal/mesenchymal"

        membership = root / "membership.tsv.gz"
        pd.DataFrame([
            {"cell_id": row["cell_id"], "final_broad_label": broad[row["cell_id"]]}
            for row in cells
        ]).to_csv(membership, sep="\t", index=False, compression="gzip")

        score_rows = []
        for row in cells:
            for candidate, release_broad in self.candidates.items():
                supported = identity[row["cell_id"]] == candidate
                score_rows.append({
                    "cell_id": row["cell_id"],
                    "source_boundary": "selected_res0.2",
                    "source_cluster": row["source_cluster"],
                    "candidate_id": candidate,
                    "release_broad_label": release_broad,
                    "normalized_evidence": 0.90 if supported else 0.05,
                    "program_score": 0.90 if supported else 0.05,
                    "family_coherent": supported,
                    "identity_core_coherent": supported,
                    "identity_core_direct": supported,
                    "release_family_coherent": supported,
                    "hard_contradiction": False,
                    "candidate_seed": supported,
                    "technical_flag": False,
                    "x": row["x"], "y": row["y"],
                })
        scores = root / "scores.tsv.gz"
        pd.DataFrame(score_rows).to_csv(
            scores, sep="\t", index=False, compression="gzip",
        )
        catalog = root / "catalog.json"
        catalog.write_text(json.dumps({
            "candidate_boundaries": [
                {"candidate_id": candidate, "release_broad_label": label}
                for candidate, label in self.candidates.items()
            ]
        }), encoding="utf-8")
        return membership, scores, catalog

    def run_review(self, root: Path, under_recall_theca: bool = False):
        membership, scores, catalog = self.build_fixture(root, under_recall_theca)
        out = root / "review"
        result = subprocess.run([
            sys.executable, str(SCRIPT),
            "--membership", str(membership),
            "--scores", str(scores),
            "--catalog", str(catalog),
            "--out", str(out),
        ], capture_output=True, text=True, check=False)
        manifest = json.loads(
            (out / "sheep_ovary_biological_quality_review.json").read_text()
        )
        return result, manifest, out

    def test_concordant_antral_follicle_passes_all_endpoints(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, manifest, _ = self.run_review(Path(tmp))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(manifest["status"], "PASS")
            self.assertEqual(
                manifest["quality_endpoints"]["oocyte_annotation_quality"]["status"],
                "PASS",
            )
            follicle = manifest["quality_endpoints"]["follicle_roi_histology"]
            self.assertEqual(follicle["status"], "PASS")
            self.assertGreaterEqual(follicle["antral_roi_n"], 1)
            self.assertFalse(manifest["formal_membership_written"])

    def test_coherent_theca_program_in_stromal_remainder_requires_iteration(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, manifest, out = self.run_review(Path(tmp), True)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertEqual(manifest["status"], "ITERATION_REQUIRED")
            actions = pd.read_csv(
                out / "biological_quality_next_actions.tsv", sep="\t",
            )
            self.assertIn(
                "theca_steroidogenic_coherent_program_under_recalled",
                set(actions.issue_code),
            )
            repaired_membership = pd.read_csv(
                manifest["membership"]["path"], sep="\t",
            )
            self.assertNotIn("Theca", set(repaired_membership.final_broad_label))
            self.assertFalse(manifest["formal_membership_written"])

    def test_missing_vascular_layer_is_detected_without_theca_mislabeling(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            membership, scores, catalog = self.build_fixture(root)
            frame = pd.read_csv(scores, sep="\t")
            mask = frame.candidate_id.isin({
                "vascular_endothelial", "pericyte_mural",
            })
            for column in (
                "family_coherent", "identity_core_coherent",
                "identity_core_direct", "release_family_coherent",
                "candidate_seed",
            ):
                frame.loc[mask, column] = False
            frame.loc[mask, ["normalized_evidence", "program_score"]] = 0.05
            frame.to_csv(scores, sep="\t", index=False, compression="gzip")
            out = root / "review"
            result = subprocess.run([
                sys.executable, str(SCRIPT), "--membership", str(membership),
                "--scores", str(scores), "--catalog", str(catalog),
                "--out", str(out),
            ], capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 2, result.stderr)
            actions = pd.read_csv(
                out / "biological_quality_next_actions.tsv", sep="\t",
            )
            self.assertIn(
                "vascular_interna_program_not_resolved",
                set(actions.issue_code),
            )
            self.assertNotIn(
                "theca_interna_label_under_recall",
                set(actions.issue_code),
            )

    def test_structural_perifollicular_program_cannot_substitute_for_theca_interna(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            membership, scores, catalog = self.build_fixture(root)
            frame = pd.read_csv(scores, sep="\t")
            steroid = frame.candidate_id == "theca_steroidogenic"
            for column in (
                "family_coherent", "identity_core_coherent",
                "identity_core_direct", "release_family_coherent",
                "candidate_seed",
            ):
                frame.loc[steroid, column] = False
            frame.loc[steroid, ["normalized_evidence", "program_score"]] = 0.05
            structural = frame.loc[frame.candidate_id == "theca_steroidogenic"].copy()
            structural["candidate_id"] = "theca_structural_perifollicular"
            structural["release_broad_label"] = ""
            is_theca_ring = structural.source_cluster == "theca_steroidogenic"
            for column in (
                "family_coherent", "identity_core_coherent",
                "identity_core_direct", "release_family_coherent",
                "candidate_seed",
            ):
                structural[column] = is_theca_ring
            structural["normalized_evidence"] = np.where(is_theca_ring, 0.95, 0.05)
            structural["program_score"] = np.where(is_theca_ring, 0.95, 0.05)
            frame = pd.concat([frame, structural], ignore_index=True)
            frame.to_csv(scores, sep="\t", index=False, compression="gzip")
            document = json.loads(catalog.read_text())
            document["candidate_boundaries"].append({
                "candidate_id": "theca_structural_perifollicular",
                "release_broad_label": "",
            })
            catalog.write_text(json.dumps(document), encoding="utf-8")
            out = root / "review"
            result = subprocess.run([
                sys.executable, str(SCRIPT), "--membership", str(membership),
                "--scores", str(scores), "--catalog", str(catalog),
                "--out", str(out),
            ], capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 2, result.stderr)
            actions = pd.read_csv(out / "biological_quality_next_actions.tsv", sep="\t")
            self.assertIn("theca_interna_program_not_resolved", set(actions.issue_code))

    def test_no_follicle_and_no_oocyte_are_not_evaluable_not_failures(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cells = self.ring("stromal", 20.0, 96, "stromal")
            membership = root / "membership.tsv"
            pd.DataFrame([
                {"cell_id": row["cell_id"], "final_broad_label": "Stromal/mesenchymal"}
                for row in cells
            ]).to_csv(membership, sep="\t", index=False)
            rows = []
            for row in cells:
                for candidate, label in self.candidates.items():
                    supported = candidate == "stromal_mesenchymal"
                    rows.append({
                        "cell_id": row["cell_id"],
                        "source_boundary": "selected_res0.1",
                        "source_cluster": "stromal",
                        "candidate_id": candidate,
                        "release_broad_label": label,
                        "normalized_evidence": 0.9 if supported else 0.05,
                        "program_score": 0.9 if supported else 0.05,
                        "family_coherent": supported,
                        "identity_core_coherent": supported,
                        "identity_core_direct": supported,
                        "release_family_coherent": supported,
                        "hard_contradiction": False,
                        "candidate_seed": supported,
                        "technical_flag": False,
                        "x": row["x"], "y": row["y"],
                    })
            scores = root / "scores.tsv"
            pd.DataFrame(rows).to_csv(scores, sep="\t", index=False)
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({
                "candidate_boundaries": [
                    {"candidate_id": candidate, "release_broad_label": label}
                    for candidate, label in self.candidates.items()
                ]
            }), encoding="utf-8")
            out = root / "review"
            result = subprocess.run([
                sys.executable, str(SCRIPT), "--membership", str(membership),
                "--scores", str(scores), "--catalog", str(catalog),
                "--out", str(out),
            ], capture_output=True, text=True, check=False)
            manifest = json.loads(
                (out / "sheep_ovary_biological_quality_review.json").read_text()
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(manifest["status"], "PASS")
            self.assertEqual(
                manifest["quality_endpoints"]["follicle_roi_histology"]["status"],
                "NOT_EVALUABLE",
            )

    def test_exact_canonical_oocyte_review_supersedes_stale_ordinary_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            membership, scores, catalog = self.build_fixture(root)
            frame = pd.read_csv(scores, sep="\t")
            stale = frame.candidate_id.eq("oocyte")
            for column in (
                "family_coherent", "identity_core_coherent",
                "identity_core_direct", "release_family_coherent",
                "candidate_seed",
            ):
                frame.loc[stale, column] = False
            frame.loc[stale, "hard_contradiction"] = True
            frame.loc[stale, ["normalized_evidence", "program_score"]] = 0.01
            frame.to_csv(scores, sep="\t", index=False, compression="gzip")
            digest = hashlib.sha256(membership.read_bytes()).hexdigest()
            canonical = root / "canonical_oocyte_review.json"
            canonical.write_text(json.dumps({
                "status": "FROZEN_OOCYTE_MEMBERSHIP",
                "membership_path": str(membership),
                "membership_sha256": digest,
                "n_canonical_cluster_cellbins": 14,
                "n_final_oocyte_cellbins": 12,
                "n_direct_hard_somatic_contradiction_retained_in_resident_broad": 2,
                "n_putative_oocyte_objects": 4,
                "selected_resolution": "oocyte_res0p1",
                "cross_resolution_jaccard": 0.95,
                "spatial_location_used_for_admission": False,
                "zona_only_admission_forbidden": True,
                "independent_non_zona_deg_gene_n": 6,
            }), encoding="utf-8")
            out = root / "review"
            result = subprocess.run([
                sys.executable, str(SCRIPT), "--membership", str(membership),
                "--scores", str(scores), "--catalog", str(catalog),
                "--canonical-oocyte-review", str(canonical),
                "--out", str(out),
            ], capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = json.loads(
                (out / "sheep_ovary_biological_quality_review.json").read_text()
            )
            oocyte = manifest["quality_endpoints"]["oocyte_annotation_quality"]
            self.assertEqual(oocyte["status"], "PASS")
            self.assertEqual(oocyte["canonical_supported_group_n"], 1)

    def test_molecularly_supported_theca_is_not_failed_for_global_fragmentation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            membership, scores, catalog = self.build_fixture(root)
            frame = pd.read_csv(scores, sep="\t")
            theca_cells = sorted(
                frame.loc[
                    frame.source_cluster.eq("theca_steroidogenic"), "cell_id"
                ].unique()
            )
            coordinates = {
                cell: (1000.0 * index, 1000.0 * (index % 7))
                for index, cell in enumerate(theca_cells)
            }
            for cell, (x, y) in coordinates.items():
                mask = frame.cell_id.eq(cell)
                frame.loc[mask, ["x", "y"]] = [x, y]
            frame.to_csv(scores, sep="\t", index=False, compression="gzip")
            out = root / "review"
            subprocess.run([
                sys.executable, str(SCRIPT), "--membership", str(membership),
                "--scores", str(scores), "--catalog", str(catalog),
                "--out", str(out),
            ], capture_output=True, text=True, check=False)
            broad = pd.read_csv(out / "broad_spatial_localization_review.tsv", sep="\t")
            theca = broad.loc[broad.broad_label.eq("Theca")].iloc[0]
            self.assertEqual(theca.status, "PASS")
            self.assertEqual(theca.identity_supported_fraction, 1.0)

    def test_targeted_repair_cannot_erase_previously_confirmed_antral_roi(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result, manifest, first_out = self.run_review(root)
            self.assertEqual(result.returncode, 0, result.stderr)
            membership = Path(manifest["membership"]["path"])
            frame = pd.read_csv(membership, sep="\t")
            frame.loc[frame.final_broad_label == "Granulosa", "final_broad_label"] = "Stromal/mesenchymal"
            erased = root / "erased_membership.tsv.gz"
            frame.to_csv(erased, sep="\t", index=False, compression="gzip")
            scores = root / "scores.tsv.gz"
            catalog = root / "catalog.json"
            out = root / "erased_review"
            second = subprocess.run([
                sys.executable, str(SCRIPT), "--membership", str(erased),
                "--scores", str(scores), "--catalog", str(catalog),
                "--expected-roi-review", str(first_out / "follicle_roi_histology_review.tsv"),
                "--out", str(out),
            ], capture_output=True, text=True, check=False)
            self.assertEqual(second.returncode, 2, second.stderr)
            actions = pd.read_csv(out / "biological_quality_next_actions.tsv", sep="\t")
            self.assertIn(
                "previously_detected_large_antral_roi_lost_after_repair",
                set(actions.issue_code),
            )


if __name__ == "__main__":
    unittest.main()
