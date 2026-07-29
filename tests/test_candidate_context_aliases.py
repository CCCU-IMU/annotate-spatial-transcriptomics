from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "annotate-spatial-transcriptomics/scripts/build_candidate_context_evidence.py"


class CandidateContextAliasTests(unittest.TestCase):
    def test_stage_replicate_suffix_is_normalized_without_identity_authority(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "context.tsv"
            result = subprocess.run([
                sys.executable, str(SCRIPT), "--species", "绵羊",
                "--tissue", "卵巢", "--reproductive-stage", "发情后期2",
                "--out", str(out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            with out.open(encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["status"], "supported")
            self.assertEqual(row["evidence_scope"], "evaluation_permission_only")
            self.assertEqual(row["identity_writeback_authority"], "false")


if __name__ == "__main__":
    unittest.main()
