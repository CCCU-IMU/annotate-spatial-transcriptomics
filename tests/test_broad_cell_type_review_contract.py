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
PACKAGE = ROOT / "annotate-spatial-transcriptomics"
SCRIPTS = PACKAGE / "scripts"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


class BroadCellTypeReviewContractTests(unittest.TestCase):
    def test_marker_manifest_includes_complete_broad_families_not_subtypes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "build_cell_type_review_marker_manifest.py"),
                "--profile", str(PACKAGE / "references/profiles/sheep_ovary.json"),
                "--catalog", str(PACKAGE / "references/profiles/sheep_ovary_candidate_lineage_catalog.json"),
                "--out", str(out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            with (out / "cell_type_review_marker_manifest.tsv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            genes = {row["gene"] for row in rows}
            self.assertTrue({"KRT7", "MUC16", "MYH11", "CYP17A1", "ZP3"} <= genes)
            self.assertFalse(any("subtype" in row["family_id"].lower() for row in rows))
            self.assertFalse(any("state" in row["family_id"].lower() for row in rows))

    def test_present_broad_cannot_close_without_precision_recall_and_spatial_conclusions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue = root / "queue.tsv"
            write_tsv(queue, [{
                "review_id": "r1", "review_mode": "broad_lineage_review",
                "target_broad_label": "Granulosa", "unit_signature": "x",
            }])
            review = root / "review.json"
            review.write_text(json.dumps({
                "stage": "post_atlas_catalog_wide_lineage_review",
                "status": "ITERATION_REQUIRED",
                "catalog_wide_double_sided_review": True,
                "review_round": 1,
                "artifacts": {"review_queue": {"path": str(queue), "sha256": sha(queue)}},
            }), encoding="utf-8")
            decisions = root / "decisions.tsv"
            write_tsv(decisions, [{
                "review_id": "r1", "review_mode": "broad_lineage_review",
                "outcome": "retain_current_cell_type", "proposed_broad_label": "",
                "evidence_basis": "query marker DEG pseudobulk and spatial review",
                "rationale": "The current cell type is retained without a complete conclusion record.",
                "membership_path": "", "membership_sha256": "",
            }])
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "validate_catalog_wide_lineage_review_decisions.py"),
                "--review-manifest", str(review), "--decisions", str(decisions),
                "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            validation = json.loads((root / "out/catalog_wide_lineage_decision_validation.json").read_text())
            self.assertEqual(validation["status"], "BLOCKED")
            self.assertTrue(any("current_member_precision" in error for error in validation["errors"]))


if __name__ == "__main__":
    unittest.main()
