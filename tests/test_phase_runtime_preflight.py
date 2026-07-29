from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "annotate-spatial-transcriptomics/scripts/validate_phase_runtime.py"


class PhaseRuntimePreflightTests(unittest.TestCase):
    def test_exact_python_and_declared_semantic_readers_pass(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            module = root / "module.py"
            module.write_text("value = 1\n", encoding="utf-8")
            document = root / "document.json"
            document.write_text('{"ok": true}\n', encoding="utf-8")
            table = root / "table.tsv"
            table.write_text("cell_id\ncell_1\n", encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(SCRIPT), "--python", sys.executable,
                "--python-script", str(module),
                "--semantic-input", f"document:json:{document}",
                "--semantic-input", f"membership:tsv:{table}",
                "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            manifest = json.loads((root / "out/phase_runtime_preflight.json").read_text())
            self.assertEqual(manifest["status"], "PASS")

    def test_missing_dependency_and_wrong_reader_fail_before_phase(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = root / "context.tsv"
            table.write_text("key\tvalue\nstage\tdiestrus\n", encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(SCRIPT), "--python", sys.executable,
                "--python-import", "module_that_must_not_exist_2200",
                "--semantic-input", f"context:json:{table}",
                "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            manifest = json.loads((root / "out/phase_runtime_preflight.json").read_text())
            self.assertEqual(manifest["status"], "BLOCKED")
            self.assertTrue(any("lacks" in error for error in manifest["errors"]))


if __name__ == "__main__":
    unittest.main()
