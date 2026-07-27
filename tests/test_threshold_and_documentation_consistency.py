from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "annotate-spatial-transcriptomics"
SCRIPTS = SKILL / "scripts"
REGISTRY = SKILL / "references/controller_thresholds_v2_2.json"
PROFILE = SKILL / "references/profiles/sheep_ovary.json"
sys.path.insert(0, str(SCRIPTS))

from controller_thresholds import load_controller_thresholds  # noqa: E402
from lineage_decision_lib import (  # noqa: E402
    DEFAULT_OBSERVATION_WRITEBACK_POLICY,
)


class ThresholdAndDocumentationConsistencyTests(unittest.TestCase):
    def test_release_critical_defaults_have_one_registry(self) -> None:
        thresholds = load_controller_thresholds(REGISTRY)
        self.assertEqual(
            DEFAULT_OBSERVATION_WRITEBACK_POLICY,
            thresholds["observation_writeback_policy"],
        )
        self.assertAlmostEqual(
            thresholds["scoring_policy"]["direct_weight"]
            + thresholds["scoring_policy"]["local_weight"],
            1.0,
        )
        for purpose in (
            "whole_tissue_cohort_partition", "cohort_identity_resolution"
        ):
            total = sum(
                thresholds["resolution_selection"][purpose][
                    "metric_weights"
                ].values()
            )
            self.assertGreater(total, 0)
            self.assertLessEqual(total, 1.0)

    def test_biological_profile_references_controller_registry(self) -> None:
        profile = json.loads(PROFILE.read_text(encoding="utf-8"))
        evidence = profile["broad_family_evidence_contract"]
        self.assertEqual(
            evidence["controller_threshold_registry"],
            "../controller_thresholds_v2_2.json",
        )
        subset = evidence["observation_subset_writeback"]
        self.assertEqual(
            subset["release_threshold_source"],
            "../controller_thresholds_v2_2.json#observation_writeback_policy",
        )
        duplicated = {
            key for key in subset
            if key.startswith("whole_subcluster_")
            or key.startswith("supported_subset_")
            or key == "maximum_contradiction_fraction"
        }
        self.assertEqual(duplicated, set())

    def test_python_and_r_runtime_paths_bind_the_registry(self) -> None:
        controller = (SCRIPTS / "run_lineage_controller.py").read_text(
            encoding="utf-8"
        )
        scorer = (SCRIPTS / "run_observation_lineage_scoring.R").read_text(
            encoding="utf-8"
        )
        subset = (SCRIPTS / "derive_candidate_local_subsets.R").read_text(
            encoding="utf-8"
        )
        selector = (SCRIPTS / "select_lineage_resolution.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('paths["threshold_registry"]', controller)
        for text in (scorer, subset, selector):
            self.assertIn("threshold-registry", text)
        for text in (scorer, subset):
            self.assertIn("controller_thresholds_v2_2.json", text)
        self.assertIn("load_controller_thresholds", selector)

        registry_consumers = (
            "build_resolution_grid_evidence.py",
            "close_exact_remainders.py",
            "materialize_parent_locked_fine_proposals.py",
            "materialize_final_release_v2_2.py",
            "validate_lineage_controller_release.py",
            "validate_residual_qc_audit.py",
        )
        for name in registry_consumers:
            source = (SCRIPTS / name).read_text(encoding="utf-8")
            self.assertIn("load_controller_thresholds", source, name)
        for name in (
            "materialize_final_release_v2_2.py",
            "validate_lineage_controller_release.py",
        ):
            source = (SCRIPTS / name).read_text(encoding="utf-8")
            self.assertNotIn("qc_n >= 50000", source, name)
            self.assertNotIn("qc_fraction >= 0.10", source, name)

    def test_readme_has_only_the_active_staged_architecture(self) -> None:
        text = (ROOT / "README.md").read_text(encoding="utf-8")
        for stale in (
            "每个初步大类独立投递一次 broad_class_recluster cohort",
            "中高置信：直接写入初步大类",
            "SCT规模 | 3,000 variable features",
        ):
            self.assertNotIn(stale, text)
        self.assertIn("4,000 variable features", text)
        self.assertIn("one-initial-cluster : one-cohort", text)
        whole_command = text.split("全组织统一前处理可直接运行：", 1)[1].split(
            "脚本会输出：", 1
        )[0]
        self.assertIn("--resolutions 0.2,0.4,0.6,0.8", whole_command)
        self.assertNotIn("--resolution-workers", whole_command)
        self.assertNotIn("--resolution-future-plan", whole_command)


if __name__ == "__main__":
    unittest.main()
