from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "annotate-spatial-transcriptomics"
SCRIPTS = SKILL / "scripts"
sys.path.insert(0, str(SCRIPTS))

from project_context import resolve_context_path  # noqa: E402


class V203RuntimeContractTests(unittest.TestCase):
    def test_context_alias_is_resolved_without_duplicate_files(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "config").mkdir()
            alias = root / "config/context.json"
            alias.write_text("{}\n", encoding="utf-8")
            self.assertEqual(resolve_context_path(root), alias)
            canonical = root / "config/biological_context.json"
            canonical.write_text("{}\n", encoding="utf-8")
            self.assertEqual(resolve_context_path(root), canonical)

    def test_residual_qc_validator_autodetects_v2_final_state(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ledger = root / "ledger.tsv"
            with ledger.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=["cell_id", "final_state", "qc_reason"], delimiter="\t"
                )
                writer.writeheader()
                writer.writerows([
                    {"cell_id": "a", "final_state": "defined_broad_only", "qc_reason": ""},
                    {"cell_id": "b", "final_state": "qc_holdout", "qc_reason": "low_information"},
                ])
            output = root / "validation.json"
            subprocess.run(
                [
                    sys.executable, str(SCRIPTS / "validate_residual_qc_audit.py"),
                    "--ledger", str(ledger), "--fraction-trigger", "1.0",
                    "--count-trigger", "999", "--out", str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            result = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(result["state_column"], "final_state")
            self.assertEqual(result["qc_reason_column"], "qc_reason")
            self.assertEqual(result["residual_qc_n"], 1)

    def test_report_and_spatial_assets_expose_existing_evidence(self) -> None:
        report = (SCRIPTS / "build_report.py").read_text(encoding="utf-8")
        for token in [
            "annotation_support_registry.tsv", "final_broad_DEG_top100.tsv",
            "final_subtype_DEG_top100.tsv", "positive_marker_evidence", "anti_marker_evidence",
        ]:
            self.assertIn(token, report)
        spatial = (SCRIPTS / "build_spatial_gene_maps.R").read_text(encoding="utf-8")
        self.assertIn('scope!="all_analysis_set_observations"', spatial)
        self.assertIn("expected-observations", spatial)

    def test_controller_accepts_zero_count_route_partitions_and_terminal_returns(self) -> None:
        route = (SCRIPTS / "register_route_attempt.py").read_text(encoding="utf-8")
        state = (SCRIPTS / "validate_state.py").read_text(encoding="utf-8")
        membership = (SCRIPTS / "audit_annotation_membership_partition.py").read_text(encoding="utf-8")
        self.assertIn("allow_empty=expected_count == 0", route)
        self.assertIn("active_cell_decision_ids", state)
        self.assertIn("terminal_residual_qc_freeze", state)
        self.assertIn("terminal_qc | seen_direct", membership)


if __name__ == "__main__":
    unittest.main()
