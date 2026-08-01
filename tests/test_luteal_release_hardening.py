from __future__ import annotations

import importlib.util
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "annotate-spatial-transcriptomics/scripts"
sys.path.insert(0, str(SCRIPTS))

from lineage_controller_lib import (  # noqa: E402
    candidate_allows_dominant_whole_subcluster,
    candidate_specific_group_release_pass,
    group_release_supported_fraction,
)


class LutealReleaseHardeningTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        catalog = json.loads((
            ROOT
            / "annotate-spatial-transcriptomics/references/profiles/"
            "sheep_ovary_candidate_lineage_catalog.json"
        ).read_text(encoding="utf-8"))
        cls.luteal = next(
            row for row in catalog["candidate_boundaries"]
            if row["candidate_id"] == "luteal_steroidogenic"
        )

    def test_shared_steroidogenesis_and_stability_do_not_release_luteal(self) -> None:
        row = {
            "observation_seed_fraction": "0.99",
            "observation_identity_core_fraction": "0.88",
            "observation_identity_core_direct_fraction": "0.88",
            "observation_required_positive_joint_direct_fraction": "0.30",
            "observation_split_discriminator_direct_fraction": "0.33",
            "marker_deg_log2fc_mean": "-0.02",
            "anti_marker_deg_log2fc_mean": "0",
            "hard_contradiction_fraction": "0",
            "cross_resolution_stable_fraction": "1",
        }
        self.assertEqual(
            group_release_supported_fraction(row, self.luteal), 0.30
        )
        self.assertFalse(candidate_specific_group_release_pass(row, self.luteal))
        self.assertFalse(
            candidate_allows_dominant_whole_subcluster(self.luteal)
        )

    def test_direct_joint_query_enriched_luteal_program_can_pass(self) -> None:
        row = {
            "observation_required_positive_joint_direct_fraction": "0.62",
            "observation_split_discriminator_direct_fraction": "0.51",
            "marker_deg_log2fc_mean": "0.90",
            "anti_marker_deg_log2fc_mean": "0.05",
            "hard_contradiction_fraction": "0.01",
        }
        self.assertTrue(candidate_specific_group_release_pass(row, self.luteal))

    def test_scorer_emits_direct_joint_evidence(self) -> None:
        source = (SCRIPTS / "run_observation_lineage_scoring.R").read_text(
            encoding="utf-8"
        )
        self.assertIn("required_family_minimum_direct_genes", source)
        self.assertIn("required_positive_families_joint_direct", source)
        self.assertIn(
            "observation_required_positive_joint_direct_fraction", source
        )

    def test_biological_validator_treats_luteal_as_restricted(self) -> None:
        module_path = SCRIPTS / "validate_sheep_ovary_biological_quality.py"
        spec = importlib.util.spec_from_file_location("sheep_quality", module_path)
        assert spec and spec.loader
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertIn("Luteal", module.RESTRICTED_BROADS)
        self.assertIn(
            "required_positive_families_joint_direct",
            module.OPTIONAL_BOOLEAN_COLUMNS,
        )


if __name__ == "__main__":
    unittest.main()
