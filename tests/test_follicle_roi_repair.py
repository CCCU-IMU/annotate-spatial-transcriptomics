from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "annotate-spatial-transcriptomics/scripts/apply_sheep_ovary_follicle_roi_repair.py"
CANDIDATES = [
    "theca_steroidogenic", "vascular_endothelial", "pericyte_mural",
    "lymphatic_endothelial", "smooth_muscle", "stromal_mesenchymal",
]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


class FollicleRoiRepairTest(unittest.TestCase):
    def build_fixture(self, root: Path, tie: bool = False, raw_assay: str = "RNA") -> list[str]:
        contract = root / "contract.json"
        contract.write_text("{}\n")
        catalog = root / "catalog.json"
        release = {
            "theca_steroidogenic": "Theca",
            "vascular_endothelial": "Vascular-associated",
            "pericyte_mural": "Vascular-associated",
            "lymphatic_endothelial": "Vascular-associated",
            "smooth_muscle": "Smooth muscle",
            "stromal_mesenchymal": "Stromal/mesenchymal",
        }
        catalog.write_text(json.dumps({"candidate_boundaries": [
            {
                "candidate_id": candidate, "candidate_role": "broad",
                "release_broad_label": broad,
            }
            for candidate, broad in release.items()
        ]}))
        membership = root / "membership.tsv"
        rows = []
        labels = {
            "a": "Stromal/mesenchymal", "b": "Stromal/mesenchymal",
            "c": "Stromal/mesenchymal", "d": "Smooth muscle",
            "e": "Oocyte", "f": "Stromal/mesenchymal",
        }
        for index, (cell, label) in enumerate(labels.items()):
            rows.append({
                "cell_id": cell, "x": index, "y": index,
                "final_broad_label": label, "final_broad_confidence": "moderate",
                "final_fine_label": "", "final_fine_confidence": "",
                "final_state": "defined_broad_only", "qc_reason": "",
                "assignment_origin": "second_round", "broad_freeze_source": "second_round",
            })
        write_tsv(membership, rows)

        roi_membership = root / "roi_membership.tsv"
        write_tsv(roi_membership, [
            {"cell_id": cell, "follicle_roi_id": "F001", "histological_shell": "outer_wall"}
            for cell in "abcd"
        ])
        layers = root / "layers.tsv"
        write_tsv(layers, [
            {"follicle_roi_id": "F001", "layer_name": name, "status": "ITERATION_REQUIRED"}
            for name in (
                "theca_interna", "vascular_interna",
                "outer_nonvascular_contractile", "outer_stromal_background",
            )
        ])
        roi_review = root / "roi_review.tsv"
        write_tsv(roi_review, [{
            "follicle_roi_id": "F001", "follicle_stage_geometry": "large_antral_candidate",
            "cavity_structure_status": "PASS",
        }])
        actions = root / "actions.tsv"
        write_tsv(actions, [{
            "endpoint": "follicle_roi_histology", "scope_id": "F001",
            "issue_code": "complete_follicle_wall_hierarchy_review_required",
            "detail": "fixture", "recommended_action": "bounded raw-count review",
        }])
        quality = root / "quality.json"
        quality.write_text(json.dumps({
            "status": "ITERATION_REQUIRED",
            "required_next_actions": {"path": str(actions), "sha256": sha(actions)},
            "quality_endpoints": {"follicle_roi_histology": {
                "roi_membership": {"path": str(roi_membership), "sha256": sha(roi_membership)},
                "layer_hierarchy": {"path": str(layers), "sha256": sha(layers)},
                "roi_review": {"path": str(roi_review), "sha256": sha(roi_review)},
            }},
        }))

        base_scores = root / "base_scores.tsv"
        repair_scores = root / "repair_scores.tsv"

        def score_rows(cells: str) -> list[dict]:
            winners = {
                "a": ("smooth_muscle", 0.90),
                "b": ("vascular_endothelial", 0.80),
                "c": ("theca_steroidogenic", 0.85),
                "d": ("stromal_mesenchymal", 0.70),
            }
            result = []
            for cell in cells:
                for candidate in CANDIDATES:
                    evidence = 0.02
                    eligible = False
                    if cell in winners and candidate == winners[cell][0]:
                        evidence = winners[cell][1]
                        eligible = True
                    if tie and cell == "a" and candidate == "vascular_endothelial":
                        evidence = 0.88
                        eligible = True
                    result.append({
                        "cell_id": cell, "candidate_id": candidate,
                        "normalized_evidence": evidence,
                        "family_coherent": eligible,
                        "identity_core_direct": eligible,
                        "release_family_coherent": eligible,
                        "hard_contradiction": False, "technical_flag": False,
                        "positive_families": "identity_core" if eligible else "",
                    })
            return result

        write_tsv(base_scores, score_rows("abcdef"))
        write_tsv(repair_scores, score_rows("abcd"))
        ancestry = root / "ancestry.json"
        ancestry.write_text(json.dumps({
            "status": "PASS", "raw_count_assay": raw_assay,
            "clustering_path": "raw_counts_SCTv2_PCA_SNN_Leiden_follicle_ROI",
        }))
        authority = root / "authority.json"
        authority.write_text(json.dumps({
            "mode": "stage_authority", "phase": "atlas_and_completeness_review",
            "annotation_contract_sha256": sha(contract),
            "pre_repair_membership": {"path": str(membership), "sha256": sha(membership)},
            "pre_repair_biological_quality": {"path": str(quality), "sha256": sha(quality)},
            "base_scores": [{"path": str(base_scores), "sha256": sha(base_scores)}],
            "repair_scores": [{"path": str(repair_scores), "sha256": sha(repair_scores)}],
            "repair_ancestry": [{"path": str(ancestry), "sha256": sha(ancestry)}],
            "candidate_catalog": {"path": str(catalog), "sha256": sha(catalog)},
        }))
        out = root / "out"
        return [
            sys.executable, str(SCRIPT), "--contract", str(contract),
            "--stage-authority", str(authority), "--membership", str(membership),
            "--quality-review", str(quality), "--catalog", str(catalog),
            "--base-scores", str(base_scores),
            "--repair-score", f"F001={repair_scores}",
            "--repair-ancestry", f"F001={ancestry}", "--out", str(out),
        ]

    def test_specific_identities_precede_stromal_remainder(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            subprocess.run(self.build_fixture(root), check=True, capture_output=True, text=True)
            result = pd.read_csv(root / "out/post_follicle_roi_repair_membership.tsv.gz", sep="\t", dtype=str).fillna("").set_index("cell_id")
            self.assertEqual(result.at["a", "final_broad_label"], "Smooth muscle")
            self.assertEqual(result.at["b", "final_broad_label"], "Vascular-associated")
            self.assertEqual(result.at["c", "final_broad_label"], "Theca")
            self.assertEqual(result.at["d", "final_broad_label"], "Stromal/mesenchymal")
            self.assertEqual(result.at["e", "final_broad_label"], "Oocyte")
            self.assertEqual(result.at["f", "final_broad_label"], "Stromal/mesenchymal")

    def test_specific_specific_tie_preserves_existing_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            subprocess.run(self.build_fixture(root, tie=True), check=True, capture_output=True, text=True)
            result = pd.read_csv(root / "out/post_follicle_roi_repair_membership.tsv.gz", sep="\t", dtype=str).fillna("").set_index("cell_id")
            self.assertEqual(result.at["a", "final_broad_label"], "Stromal/mesenchymal")

    def test_sct_corrected_assay_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proc = subprocess.run(self.build_fixture(root, raw_assay="SCT"), capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("non-SCT raw-count", proc.stderr + proc.stdout)


if __name__ == "__main__":
    unittest.main()
