from __future__ import annotations

import csv
import gzip
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "annotate-spatial-transcriptomics"
BUILD = SKILL / "scripts/build_reliable_fine_annotation.py"
VALIDATE = SKILL / "scripts/validate_final_fine_annotation.py"
CATALOG = SKILL / "references/profiles/sheep_ovary_candidate_lineage_catalog.json"


FIELDS = [
    "cell_id", "sample_id", "decision_id", "generation", "state", "broad_label", "fine_label",
    "confidence", "evidence_status", "validation_status", "validation_artifact", "validation_feature_scope",
    "source_run_id", "recluster_cohort_id", "route", "assignment_mode", "fine_anchor_eligible",
    "decision_version", "supersedes", "closed", "next_action", "created_at", "final_state",
    "final_broad_label", "final_fine_label", "fine_confidence", "final_fine_confidence",
    "final_fine_eligible",
]


def row(cell: str, broad: str, fine: str = "", state: str = "defined_broad_only") -> dict[str, str]:
    value = {field: "" for field in FIELDS}
    value.update(
        cell_id=cell,
        sample_id="S1",
        decision_id=f"old_{cell}",
        generation="1",
        state=state,
        broad_label=broad,
        fine_label=fine,
        confidence="high",
        fine_anchor_eligible="true" if fine else "false",
        decision_version="2.0.4",
        closed="true",
        final_state=state,
        final_broad_label=broad,
        final_fine_label=fine,
        fine_confidence="high" if fine else "",
        final_fine_confidence="high" if fine else "",
        final_fine_eligible="true" if fine else "false",
    )
    return value


class FineAnnotationTests(unittest.TestCase):
    def test_writer_preserves_broad_and_repairs_fine_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ledger = root / "ledger.tsv.gz"
            with gzip.open(ledger, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerow(row("g1", "Granulosa", "Cumulus-like granulosa", "defined_broad_only"))
                writer.writerow(row("i1", "Immune"))
                writer.writerow(row("q1", "" , "", "qc_holdout"))
            assignments = root / "assignments.tsv"
            assignments.write_text(
                "cell_id\tparent_broad_label\tfinal_fine_label\tconfidence\tsource_partition\tassignment_mode\n"
                "i1\tImmune\tB/plasma\thigh\timmune_c2\tstable_parent_cohort_fine_direct\n",
                encoding="utf-8",
            )
            evidence = root / "evidence.tsv"
            evidence.write_text("status\tPASS\n", encoding="utf-8")
            out = root / "final.tsv.gz"
            manifest = root / "manifest.json"
            subprocess.run(
                [
                    "python", str(BUILD), "--ledger", str(ledger), "--assignments", str(assignments),
                    "--catalog", str(CATALOG), "--out", str(out), "--manifest", str(manifest),
                    "--source-run-id", "S1_P61_FINE_A02", "--validation-artifact", str(evidence),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            with gzip.open(out, "rt", encoding="utf-8", newline="") as handle:
                rows = {item["cell_id"]: item for item in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["g1"]["final_state"], "defined_fine")
            self.assertEqual(rows["i1"]["final_fine_label"], "B/plasma")
            self.assertEqual(rows["i1"]["final_broad_label"], "Immune")
            self.assertEqual(rows["q1"]["final_state"], "qc_holdout")
            result = subprocess.run(
                ["python", str(VALIDATE), str(out), "--catalog", str(CATALOG)],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(json.loads(result.stdout)["status"], "PASS")
            self.assertFalse(json.loads(manifest.read_text(encoding="utf-8"))["broad_labels_modified"])

    def test_writer_rejects_cross_parent_assignment(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ledger = root / "ledger.tsv.gz"
            with gzip.open(ledger, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerow(row("x1", "Immune"))
            assignments = root / "assignments.tsv"
            assignments.write_text(
                "cell_id\tparent_broad_label\tfinal_fine_label\tconfidence\n"
                "x1\tGranulosa\tCumulus-like granulosa\thigh\n",
                encoding="utf-8",
            )
            evidence = root / "evidence.tsv"
            evidence.write_text("status\tPASS\n", encoding="utf-8")
            result = subprocess.run(
                [
                    "python", str(BUILD), "--ledger", str(ledger), "--assignments", str(assignments),
                    "--catalog", str(CATALOG), "--out", str(root / "out.tsv.gz"),
                    "--manifest", str(root / "manifest.json"), "--source-run-id", "run",
                    "--validation-artifact", str(evidence),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("crosses locked broad parent", result.stderr)


if __name__ == "__main__":
    unittest.main()
