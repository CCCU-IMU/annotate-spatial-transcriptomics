from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "annotate-spatial-transcriptomics/scripts/plan_cohort_resources.py"


class CohortResourcePlanTests(unittest.TestCase):
    def test_large_cohort_forces_one_resolution_worker(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rds = root / "input.rds"
            rds.write_bytes(b"0" * 1024)
            membership = root / "membership.tsv"
            with membership.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["cell_id"], delimiter="\t")
                writer.writeheader()
                writer.writerows({"cell_id": f"c{i}"} for i in range(20_000))
            out = root / "out"
            result = subprocess.run([
                sys.executable, str(SCRIPT), "--rds", str(rds),
                "--membership", str(membership), "--resolution-workers", "5",
                "--out", str(out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            plan = json.loads((out / "cohort_resource_plan.json").read_text())
            self.assertEqual(plan["effective_resolution_workers"], 1)
            self.assertEqual(plan["recommended_cpu"], 64)

    def test_declared_undersized_memory_blocks_before_recluster(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rds = root / "input.rds"
            rds.write_bytes(b"0" * 1024)
            membership = root / "membership.tsv"
            membership.write_text("cell_id\nc1\n", encoding="utf-8")
            out = root / "out"
            result = subprocess.run([
                sys.executable, str(SCRIPT), "--rds", str(rds),
                "--membership", str(membership), "--allocated-memory-gb", "1",
                "--out", str(out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertEqual(
                json.loads((out / "cohort_resource_plan.json").read_text())["status"],
                "BLOCKED",
            )


if __name__ == "__main__":
    unittest.main()
