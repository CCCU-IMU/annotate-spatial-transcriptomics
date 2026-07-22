from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "annotate-spatial-transcriptomics"
SCRIPTS = SKILL / "scripts"
PROFILE = SKILL / "references/profiles/sheep_ovary.json"
CATALOG = SKILL / "references/profiles/sheep_ovary_candidate_lineage_catalog.json"


def run(*args: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True)


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class V110ContractTests(unittest.TestCase):
    def test_sheep_default_broad_candidates_have_two_explicit_families(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / "family_contract.json"
            result = run(
                SCRIPTS / "validate_broad_marker_family_contract.py",
                "--profile", PROFILE, "--catalog", CATALOG, "--out", out,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            payload = json.loads(out.read_text())
            self.assertTrue(all(row["positive_family_n"] >= 2 for row in payload["candidates"]))

    def test_single_family_default_profile_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            profile = json.loads(PROFILE.read_text())
            profile["lineages"]["stromal"]["positive_families"] = {"all_markers": ["DCN", "LUM"]}
            broken = Path(temp) / "profile.json"
            broken.write_text(json.dumps(profile))
            result = run(
                SCRIPTS / "validate_broad_marker_family_contract.py",
                "--profile", broken, "--catalog", CATALOG, "--out", Path(temp) / "audit.json",
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("stromal_mesenchymal", result.stdout)

    def test_atlas_router_returns_unlabeled_qc_but_reviews_defined_disagreement(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            table = root / "concordance.tsv"
            table.write_text(
                "cell_id\tprimary_state\tprimary_broad_label\tatlas_broad_label\tatlas_tier\tatlas_ood\tontology_conflict\tmarker_audit\n"
                "a\tqc_holdout\t\tGranulosa\tmoderate_only\tfalse\tfalse\tfail\n"
                "b\tdefined\tStromal/mesenchymal\tGranulosa\thigh\tfalse\tfalse\tpass\n"
                "c\tqc_holdout\t\tOocyte\thigh\tfalse\tfalse\tpass\n"
                "d\tqc_holdout\t\tImmune\thigh\ttrue\tfalse\tpass\n"
            )
            config = root / "config.json"
            config.write_text(json.dumps({"atlas_direct_return_excluded_labels": ["Oocyte", "Theca", "Epithelial/mesothelial"]}))
            result = run(SCRIPTS / "route_calibrated_atlas_by_primary_state.py",
                         "--concordance", table, "--config", config, "--out", root / "out")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            rows = {row["cell_id"]: row for row in read_gzip_tsv(root / "out/atlas_state_routing.tsv.gz")}
            self.assertEqual(rows["a"]["atlas_state_route"], "direct_qc_broad_return")
            self.assertEqual(rows["a"]["proposed_broad_label"], "Granulosa")
            self.assertEqual(rows["b"]["atlas_state_route"], "defined_label_disagreement_review")
            self.assertEqual(rows["b"]["proposed_broad_label"], "Stromal/mesenchymal")
            self.assertEqual(rows["c"]["atlas_state_route"], "retain_qc")
            self.assertEqual(rows["d"]["atlas_state_route"], "retain_qc")

    def test_global_concordance_uses_state_aware_qc_policy(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ledger = root / "ledger.tsv"
            ledger.write_text(
                "cell_id\tanalysis_scope\tfinal_state\tfinal_broad_label\tfinal_fine_label\tsource_cluster\n"
                "a\tanalysis_set\tqc_holdout\t\t\t1\n"
                "b\tanalysis_set\tqc_holdout\t\t\t1\n"
                "c\tanalysis_set\tdefined_broad_only\tGranulosa\t\t2\n"
            )
            mapping = root / "mapping.tsv"
            mapping.write_text(
                "cell_id\tpredicted_label\tmapping_tier\tconsensus_pass\tout_of_distribution\tontology_conflict\n"
                "a\tGranulosa\tmoderate_only\tfalse\tfalse\tfalse\n"
                "b\tOocyte\thigh\ttrue\tfalse\tfalse\n"
                "c\tGranulosa\thigh\ttrue\tfalse\tfalse\n"
            )
            out = root / "audit"
            built = run(SCRIPTS / "build_global_atlas_concordance.py", "--cell-ledger", ledger,
                        "--atlas-mapping", mapping, "--out", out,
                        "--atlas-direct-return-excluded-label", "Oocyte")
            self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
            rows = {row["cell_id"]: row for row in read_gzip_tsv(out / "all_cell_atlas_concordance.tsv.gz")}
            self.assertEqual(rows["a"]["comparison_status"], "qc_writeback_candidate")
            self.assertEqual(rows["b"]["comparison_status"], "qc_reject")
            validated = run(SCRIPTS / "validate_global_atlas_concordance.py", "--audit-root", out)
            self.assertEqual(validated.returncode, 0, validated.stdout + validated.stderr)

    def test_passing_oocyte_cluster_materializes_all_canonical_members(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            canonical = root / "canonical.tsv"
            canonical.write_text(
                "cell_id\trecluster_cluster\tstrict_seed\tspatial_object_id\n"
                "a\t2\ttrue\tobj1\n"
                "b\t2\tfalse\t\n"
                "c\t2\tfalse\t\n"
                "d\t3\tfalse\t\n"
            )
            passing = root / "passing.tsv"
            passing.write_text("recluster_cluster\tadjudication_status\n2\tpass\n")
            result = run(SCRIPTS / "materialize_oocyte_cluster_membership.py",
                         "--canonical-membership", canonical, "--passing-clusters", passing,
                         "--out", root / "out")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (root / "out/materialized_oocyte_membership.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["cell_id"] for row in rows}, {"a", "b", "c"})

            context = root / "context.tsv"
            context.write_text("cell_id\na\nb\nc\nd\ne\n")
            incomplete = root / "incomplete.tsv"
            incomplete.write_text("cell_id\tfinal_broad_label\na\tOocyte\n")
            validation = run(SCRIPTS / "validate_oocyte_context_boundary.py",
                             "--canonical-membership", canonical, "--context-membership", context,
                             "--outcomes", incomplete, "--passing-clusters", passing,
                             "--out", root / "boundary.json")
            self.assertEqual(validation.returncode, 2)
            self.assertIn("were omitted", validation.stdout)

    def test_large_residual_qc_requires_hash_bound_upstream_audit(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ledger = root / "ledger.tsv"
            ledger.write_text(
                "cell_id\tannotation_status\n" +
                "".join(f"q{i}\tqc_holdout\n" for i in range(11)) +
                "".join(f"d{i}\tdefined_broad_only\n" for i in range(89))
            )
            blocked = run(SCRIPTS / "validate_residual_qc_audit.py", "--ledger", ledger,
                          "--count-trigger", "50000", "--fraction-trigger", "0.10",
                          "--out", root / "blocked.json")
            self.assertEqual(blocked.returncode, 2)
            audit = root / "audit.json"
            audit.write_text(json.dumps({
                "status": "PASS", "cell_ledger_sha256": hashlib.sha256(ledger.read_bytes()).hexdigest(),
                "residual_qc_n": 11, "unresolved_query_lineage_signals": 0,
                "completed_reviews": [
                    "initial_broad_resolution_recall_review",
                    "selected_plus_two_higher_catalog_scan",
                    "large_label_embedded_program_review",
                    "atlas_tier_census_review",
                ],
            }))
            passed = run(SCRIPTS / "validate_residual_qc_audit.py", "--ledger", ledger,
                         "--audit", audit, "--count-trigger", "50000", "--fraction-trigger", "0.10",
                         "--out", root / "passed.json")
            self.assertEqual(passed.returncode, 0, passed.stdout + passed.stderr)

    def test_vascular_release_parent_contains_endothelial_and_pericyte_children(self) -> None:
        profile = json.loads(PROFILE.read_text())
        taxonomy = profile["release_taxonomy"]
        self.assertIn("Vascular-associated", taxonomy["default_biological_broad_classes"])
        self.assertNotIn("Pericyte/mural", taxonomy["evidence_dependent_standalone_classes"])
        self.assertEqual(
            set(taxonomy["vascular_hierarchy"]["fine_children"]),
            {"Blood endothelial", "Lymphatic endothelial", "Pericyte/mural"},
        )


if __name__ == "__main__":
    unittest.main()
