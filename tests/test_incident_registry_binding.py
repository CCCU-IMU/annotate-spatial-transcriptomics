from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "annotate-spatial-transcriptomics/scripts/validate_incident_registry.py"
FIELDS = [
    "incident_id", "scheduler_job_id", "failure_class", "failure_stage",
    "symptom", "root_cause", "failure_boundary", "accepted_prior_artifacts",
    "repair_action", "repair_verification", "state_mutated",
    "biological_labels_changed", "skill_prevention_candidate",
    "regression_test_candidate", "status", "evidence_paths",
]


def write_registry(path: Path, incident_ids: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        for incident_id in incident_ids:
            writer.writerow({
                "incident_id": incident_id, "scheduler_job_id": "1",
                "failure_class": "runtime", "failure_stage": "test",
                "symptom": "failure", "root_cause": "fixture",
                "failure_boundary": "test only",
                "accepted_prior_artifacts": "none",
                "repair_action": "fixed", "repair_verification": "pass",
                "state_mutated": "false", "biological_labels_changed": "false",
                "skill_prevention_candidate": "true",
                "regression_test_candidate": "true", "status": "closed",
                "evidence_paths": "fixture",
            })


class IncidentRegistryBindingTests(unittest.TestCase):
    def test_appending_incident_invalidates_saved_validation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            registry = root / "incident_registry.tsv"
            validation = root / "validation.json"
            write_registry(registry, ["i1"])
            first = subprocess.run([
                sys.executable, str(SCRIPT), str(registry), "--out", str(validation),
            ], capture_output=True, text=True)
            self.assertEqual(first.returncode, 0, first.stderr or first.stdout)
            saved = json.loads(validation.read_text())
            self.assertEqual(saved["registry"]["row_count"], 1)
            write_registry(registry, ["i1", "i2"])
            second = subprocess.run([
                sys.executable, str(SCRIPT), str(registry),
                "--verify-existing", str(validation),
            ], capture_output=True, text=True)
            self.assertEqual(second.returncode, 2)
            result = json.loads(second.stdout)
            self.assertIn("stale", " ".join(result["errors"]))


if __name__ == "__main__":
    unittest.main()
