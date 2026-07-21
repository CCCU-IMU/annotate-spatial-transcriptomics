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


class V18ContractTests(unittest.TestCase):
    def test_new_project_enables_continuous_lineage_memory(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp) / "project"
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "init_annotation_project.py"),
                "--sample", "s1", "--input-root", temp, "--project-root", str(project),
                "--modality", "spatial", "--observation-unit", "cellbin",
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            config = json.loads((project / "config/project.json").read_text())
            self.assertEqual(config["framework_version"], "2.0.0")
            self.assertTrue(config["continuous_open_world_lineage_scan_required"])
            self.assertEqual(config["lineage_signal_resolution_policy"], "selected_plus_two_higher_available")
            self.assertTrue((project / "state/lineage_signal_boundary_registry.tsv").is_file())
            self.assertTrue((project / "state/lineage_signal_registry.tsv").is_file())

    def test_empty_signal_memory_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp) / "project"
            subprocess.run([
                sys.executable, str(SCRIPTS / "init_annotation_project.py"),
                "--sample", "s1", "--input-root", temp, "--project-root", str(project),
                "--modality", "single-cell",
            ], check=True, capture_output=True, text=True)
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "validate_lineage_signal_coverage.py"), str(project),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            payload = json.loads(result.stdout)
            self.assertEqual(payload["status"], "BLOCKED")
            self.assertTrue(any("whole_tissue" in error for error in payload["errors"]))

    def test_skill_forbids_parent_filtered_subcluster_search(self) -> None:
        controller = (SKILL / "references/direct-lineage-controller.md").read_text(encoding="utf-8")
        skill = (SKILL / "SKILL.md").read_text(encoding="utf-8")
        self.assertIn("must not narrow the candidates", controller)
        self.assertIn("parent broad label is provenance, never a search-space restriction", skill)
        self.assertIn("watch and carry forward", skill)

    def test_signal_registry_contains_explicit_memory_and_closure_fields(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp) / "project"
            subprocess.run([
                sys.executable, str(SCRIPTS / "init_annotation_project.py"),
                "--sample", "s1", "--input-root", temp, "--project-root", str(project),
                "--modality", "single-cell",
            ], check=True, capture_output=True, text=True)
            with (project / "state/lineage_signal_registry.tsv").open(newline="", encoding="utf-8") as handle:
                fields = csv.DictReader(handle, delimiter="\t").fieldnames or []
            for field in ("signal_status", "required_action", "review_status", "resolution_outcome", "closure_rationale"):
                self.assertIn(field, fields)


if __name__ == "__main__":
    unittest.main()
