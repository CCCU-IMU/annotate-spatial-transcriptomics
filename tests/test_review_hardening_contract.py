from __future__ import annotations

import csv
import gzip
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


def run(script: str, *args: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPTS / script), *(str(value) for value in args)],
        text=True,
        capture_output=True,
        check=False,
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ReviewHardeningContractTests(unittest.TestCase):
    def test_membership_patch_joins_by_cell_id_and_preserves_base_order(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            base = root / "base.tsv"
            base.write_text(
                "cell_id\tfinal_broad_label\tstate_annotations\n"
                "c\tStromal/mesenchymal\t\n"
                "a\tQC/Unknown\t\n"
                "b\tTheca\told\n",
                encoding="utf-8",
            )
            proposal = root / "proposal.tsv"
            proposal.write_text(
                "cell_id\tfinal_broad_label\tstate_annotations\n"
                "b\tLuteal\tupdated\n"
                "a\tOocyte\t\n",
                encoding="utf-8",
            )
            output = root / "final.tsv.gz"
            context = root / "context.tsv"
            context.write_text(
                "candidate_id\tstatus\treason\n"
                "luteal_steroidogenic\tsupported\tfixture stage permits evaluation\n",
                encoding="utf-8",
            )
            result = run(
                "apply_cell_id_membership_patch.py",
                "--base-membership", base,
                "--proposal", proposal,
                "--update-column", "final_broad_label",
                "--update-column", "state_annotations",
                "--catalog", CATALOG,
                "--context-evidence", context,
                "--out", output,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            rows = read_tsv(output)
            self.assertEqual([row["cell_id"] for row in rows], ["c", "a", "b"])
            self.assertEqual([row["final_broad_label"] for row in rows], [
                "Stromal/mesenchymal", "Oocyte", "Luteal",
            ])
            self.assertEqual(rows[0]["state_annotations"], "")
            self.assertEqual(rows[2]["state_annotations"], "updated")
            manifest = json.loads(
                (root / "final.tsv.gz.cell_id_patch_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertTrue(manifest["base_order_preserved"])
            self.assertFalse(manifest["positional_assignment_used"])

    def test_membership_patch_rejects_duplicate_or_foreign_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            base = root / "base.tsv"
            base.write_text("cell_id\tlabel\na\tA\nb\tB\n", encoding="utf-8")
            duplicate = root / "duplicate.tsv"
            duplicate.write_text("cell_id\tlabel\na\tX\na\tY\n", encoding="utf-8")
            result = run(
                "apply_cell_id_membership_patch.py",
                "--base-membership", base, "--proposal", duplicate,
                "--update-column", "label", "--out", root / "dup.tsv",
            )
            self.assertNotEqual(result.returncode, 0)
            foreign = root / "foreign.tsv"
            foreign.write_text("cell_id\tlabel\nz\tX\n", encoding="utf-8")
            result = run(
                "apply_cell_id_membership_patch.py",
                "--base-membership", base, "--proposal", foreign,
                "--update-column", "label", "--out", root / "foreign_out.tsv",
            )
            self.assertNotEqual(result.returncode, 0)

    def test_oocyte_full_cluster_only_allows_typed_hard_exclusions(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            canonical = root / "canonical.tsv.gz"
            with gzip.open(canonical, "wt", encoding="utf-8") as handle:
                handle.write(
                    "cell_id\trecluster_cluster\tstrict_seed\n"
                    "a\t2\ttrue\n"
                    "b\t2\tfalse\n"
                    "c\t2\tfalse\n"
                    "d\t3\tfalse\n"
                )
            passing = root / "passing.tsv"
            passing.write_text(
                "recluster_cluster\tadjudication_status\n2\tpass\n",
                encoding="utf-8",
            )
            exclusions = root / "exclusions.tsv"
            exclusions.write_text(
                "cell_id\texclusion_class\n"
                "c\tdirect_multifamily_somatic_hard_contradiction\n",
                encoding="utf-8",
            )
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({"candidate_boundaries": [{
                "candidate_id": "oocyte", "candidate_role": "broad",
                "release_broad_label": "Oocyte",
            }]}), encoding="utf-8")
            result = run(
                "materialize_oocyte_cluster_membership.py",
                "--canonical-membership", canonical,
                "--passing-clusters", passing,
                "--catalog", catalog,
                "--explicit-exclusions", exclusions,
                "--out", root / "out",
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            rows = read_tsv(root / "out/materialized_oocyte_membership.tsv")
            self.assertEqual({row["cell_id"] for row in rows}, {"a", "b"})

    def test_sheep_ovary_profile_encodes_review_boundaries(self) -> None:
        profile = json.loads(PROFILE.read_text(encoding="utf-8"))
        catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
        by_id = {row["candidate_id"]: row for row in catalog["candidate_boundaries"]}

        luteal = by_id["luteal_steroidogenic"]
        self.assertEqual(luteal["release_broad_label"], "Luteal")
        self.assertEqual(luteal["required_positive_families"], [
            "luteal_steroidogenic_core", "corpus_luteum_identity",
        ])
        self.assertEqual(
            luteal["seed_required_positive_families"],
            ["corpus_luteum_identity"],
        )
        self.assertEqual(luteal["required_family_minimum_direct_genes"], 2)
        self.assertEqual(luteal["whole_subcluster_support_metric"], "required_joint_direct")
        self.assertFalse(
            luteal["whole_subcluster_release_policy"]["allow_dominant_identity_route"]
        )
        self.assertIn(
            "OXT",
            profile["lineages"]["luteal_steroidogenic"]
            ["positive_families"]["corpus_luteum_identity"],
        )

        theca = by_id["theca_steroidogenic"]
        self.assertEqual(theca["context_requirements"], [
            "steroidogenic/androgenic program",
        ])
        self.assertIn("post hoc", theca["anatomy_review_role"])

        oocyte = profile["context_specific_identity_rules"]["oocyte"]
        self.assertIn("exactly one", oocyte["zero_census_targeted_cohort_trigger"])
        self.assertIn("all observations", oocyte["targeted_cohort_policy"])


if __name__ == "__main__":
    unittest.main()
