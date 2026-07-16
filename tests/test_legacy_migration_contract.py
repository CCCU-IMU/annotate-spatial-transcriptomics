from __future__ import annotations

import importlib.util
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "annotate-spatial-transcriptomics"


class LegacyMigrationContract(unittest.TestCase):
    def test_retired_controller_is_isolated_from_active_scripts(self) -> None:
        self.assertTrue((SKILL / "references/legacy/multi-route-controller.md").is_file())
        self.assertTrue((SKILL / "scripts/legacy/multiroute_lib.py").is_file())
        self.assertFalse((SKILL / "scripts/multiroute_lib.py").exists())
        self.assertFalse((SKILL / "scripts/audit_multiroute_state.py").exists())

    def test_three_view_helper_requires_explicit_migration_mode(self) -> None:
        path = SKILL / "scripts/legacy/build_annotation_views.py"
        text = path.read_text(encoding="utf-8")
        self.assertIn("--legacy-migration-only", text)
        self.assertIn("strict/inclusive/display generation is retired", text)

    def test_legacy_confidence_is_readable_only_through_migration_helper(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "confidence_lib", SKILL / "scripts/confidence_lib.py"
        )
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        self.assertEqual(module.canonical("medium_high", migration_mode=True), "moderate")
        with self.assertRaises(ValueError):
            module.canonical("medium_high", migration_mode=False)

    def test_v160_rctd_tiers_migrate_once_to_canonical_values(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); (root / "state").mkdir(); (root / "provenance").mkdir()
            route = root / "state/route_attempt_registry.tsv"
            fields = ["route_attempt_id", "rctd_extreme_n", "rctd_high_n", "rctd_medium_low_n"]
            with route.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader(); writer.writerow({"route_attempt_id": "r1", "rctd_extreme_n": 2, "rctd_high_n": 3, "rctd_medium_low_n": 5})
            result = subprocess.run([sys.executable, str(SKILL / "scripts/migrate_project_v1_6_0_to_v1_6_1.py"), str(root)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with route.open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual((row["rctd_high_n"], row["rctd_moderate_n"], row["rctd_low_n"]), ("2", "3", "5"))
            self.assertNotIn("rctd_extreme_n", row)


if __name__ == "__main__":
    unittest.main()
