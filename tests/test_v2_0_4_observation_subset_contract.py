from __future__ import annotations

import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "annotate-spatial-transcriptomics/scripts"
sys.path.insert(0, str(SCRIPTS))

from validate_broad_class_completeness import validate  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_registry(path: Path, evidence: Path) -> None:
    fields = [
        "review_id", "candidate_lineage", "final_n_observations", "census_status",
        "query_marker_program_review", "deg_review", "spatial_morphology_review",
        "observation_level_review", "large_label_embedding_review", "decision",
        "evidence_artifact", "evidence_artifact_sha256", "closure_rationale", "status", "supersedes",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({
            "review_id": "r1", "candidate_lineage": "granulosa", "final_n_observations": "100",
            "census_status": "present", "query_marker_program_review": "reviewed", "deg_review": "reviewed",
            "spatial_morphology_review": "reviewed", "observation_level_review": "reviewed",
            "large_label_embedding_review": "reviewed", "decision": "supported",
            "evidence_artifact": str(evidence), "evidence_artifact_sha256": sha(evidence),
            "closure_rationale": "complete numerical source and fine-candidate audit",
            "status": "audited", "supersedes": "",
        })


class V204ObservationSubsetContractTests(unittest.TestCase):
    def fixture(self, root: Path) -> tuple[Path, Path]:
        (root / "config").mkdir(); (root / "state").mkdir(); (root / "evidence").mkdir()
        (root / "config/project.json").write_text(json.dumps({
            "framework_version": "2.0.0", "numerical_broad_completeness_evidence_required": True,
            "terminal_return_purity_audit_required": True, "complete_fine_candidate_audit_required": True,
            "observation_writeback_policy": {
                "whole_subcluster_min_lineage_supported_fraction": 0.5,
                "whole_subcluster_min_purity_margin": 0.2,
                "supported_subset_min_lineage_supported_fraction": 0.7,
                "supported_subset_min_purity_margin": 0.3,
                "present_label_min_lineage_supported_fraction": 0.5,
                "present_label_min_purity_margin": 0.2,
                "maximum_contradiction_fraction": 0.05,
            },
        }))
        catalog = root / "catalog.json"
        catalog.write_text(json.dumps({
            "catalog_id": "test", "candidate_boundaries": [{
                "candidate_id": "granulosa", "review_required": True, "release_level": "default_broad_candidate",
            }],
            "machine_actionable_fine_candidate_catalog": {"granulosa": [
                {"candidate_id": "granulosa_cumulus_like"},
                {"candidate_id": "granulosa_mural_estrogenic"},
            ]},
        }))
        evidence = root / "evidence/granulosa.json"
        evidence.write_text(json.dumps({
            "schema_version": "2.0", "candidate_lineage": "granulosa", "final_n_observations": 100,
            "own_program": {"positive_family_count": 2, "family_detection_fractions": {"identity": 0.9, "support": 0.8}, "program_score": 2.0},
            "strongest_competing_program": {"candidate_lineage": "stromal_mesenchymal", "program_score": 0.5},
            "observation_level_purity": {"status": "PASS", "lineage_supported_fraction": 0.9, "strongest_competing_fraction": 0.1, "contradiction_fraction": 0.01},
            "spatial_morphology": {"status": "PASS", "rationale": "Coherent follicular-wall spatial morphology."},
            "source_writeback_audits": [{
                "return_id": "ret1", "n_observations": 100, "membership_sha256": "a" * 64,
                "return_scope": "supported_subset", "status": "PASS",
                "lineage_supported_fraction": 0.9, "strongest_competing_fraction": 0.1,
                "contradiction_fraction": 0.01,
            }],
            "fine_candidate_audit": [
                {"candidate_id": "granulosa_cumulus_like", "status": "refuted", "evidence_channels": ["full_feature", "spatial"]},
                {"candidate_id": "granulosa_mural_estrogenic", "status": "supported", "evidence_channels": ["full_feature", "cross_resolution"]},
            ],
        }))
        write_registry(root / "state/broad_class_completeness_registry.tsv", evidence)
        return catalog, evidence

    def test_complete_source_and_fine_audits_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); catalog, _ = self.fixture(root)
            result = validate(root, catalog)
            self.assertEqual(result["status"], "PASS", result["errors"])

    def test_low_purity_source_return_and_unsearched_fine_candidate_block(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); catalog, evidence = self.fixture(root)
            payload = json.loads(evidence.read_text())
            payload["source_writeback_audits"][0]["lineage_supported_fraction"] = 0.2
            payload["fine_candidate_audit"] = payload["fine_candidate_audit"][:1]
            evidence.write_text(json.dumps(payload))
            write_registry(root / "state/broad_class_completeness_registry.tsv", evidence)
            result = validate(root, catalog)
            self.assertEqual(result["status"], "BLOCKED")
            self.assertTrue(any("scope-specific purity" in error for error in result["errors"]), result["errors"])
            self.assertTrue(any("complete parent-specific catalog" in error for error in result["errors"]), result["errors"])


if __name__ == "__main__":
    unittest.main()
