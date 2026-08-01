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
CATALOG_PATH = (
    SKILL / "references/profiles/sheep_ovary_candidate_lineage_catalog.json"
)
sys.path.insert(0, str(SCRIPTS))

from lineage_controller_lib import (  # noqa: E402
    candidate_core_seed,
    catalog_candidates,
    hard_contradiction,
    index_scores,
    resolve_overlap,
    validate_subset,
)


CATALOG_DOCUMENT = json.loads(CATALOG_PATH.read_text(encoding="utf-8"))
CATALOG = catalog_candidates(CATALOG_DOCUMENT)


def score_row(
    cell: str,
    candidate_id: str,
    supported: bool,
    *,
    score: float = 0.90,
    families: list[str] | None = None,
) -> dict[str, str]:
    candidate = CATALOG[candidate_id]
    positive = (
        families
        if families is not None
        else list(candidate.get("required_positive_families", []))
    )
    if not supported:
        positive = []
    return {
        "cell_id": cell,
        "source_boundary": "cohort",
        "source_cluster": "mixed_0",
        "candidate_id": candidate_id,
        "program_score": str(score if supported else 0.0),
        "normalized_evidence": str(score if supported else 0.0),
        "positive_family_count": str(len(positive)),
        "positive_families": ";".join(positive),
        "positive_gene_count": str(2 * len(positive)),
        "family_coherent": str(bool(positive)).lower(),
        "identity_core_coherent": str(bool(positive)).lower(),
        "release_family_coherent": str(
            set(candidate.get("required_positive_families", [])) <= set(positive)
        ).lower(),
        "candidate_seed": str(bool(positive)).lower(),
        "direct_signal": "0.8" if supported else "0",
        "local_signal": "0.2" if supported else "0",
        "direct_anti_gene_count": "0",
        "direct_anti_family_count": "0",
        "hard_contradiction": "false",
        "specificity_priority": str(candidate.get("specificity_priority", 0)),
    }


class V22VascularTaxonomyTests(unittest.TestCase):
    def test_final_taxonomy_has_three_independent_vascular_boundary_broads(self) -> None:
        expected = {
            "vascular_endothelial": ("Endothelial", "vascular_endothelium"),
            "pericyte_mural": ("Pericyte/mural", "mural_mesenchyme"),
            "smooth_muscle": ("Smooth muscle", "contractile_mesenchyme"),
        }
        for candidate_id, (label, family) in expected.items():
            candidate = CATALOG[candidate_id]
            self.assertEqual(candidate["candidate_role"], "broad")
            self.assertEqual(candidate["release_broad_label"], label)
            self.assertEqual(candidate["family"], family)
            self.assertFalse(candidate.get("parent_broad_label"))
        self.assertIn(
            "Vascular-associated",
            CATALOG_DOCUMENT["taxonomy_policy"][
                "forbidden_runtime_release_labels"
            ],
        )
        self.assertNotIn(
            "Vascular-associated",
            {candidate.get("release_broad_label") for candidate in CATALOG.values()},
        )

    def test_pure_endothelial_pericyte_and_smooth_identity_cores_seed(self) -> None:
        for candidate_id in (
            "vascular_endothelial", "pericyte_mural", "smooth_muscle"
        ):
            candidate = CATALOG[candidate_id]
            row = score_row("pure", candidate_id, True)
            self.assertTrue(candidate_core_seed(row, candidate), candidate_id)
            self.assertEqual(
                set(row["positive_families"].split(";")),
                set(candidate["required_positive_families"]),
            )

    def test_acta2_tagln_only_cannot_seed_pericyte_or_smooth_muscle(self) -> None:
        pericyte = score_row(
            "contractile", "pericyte_mural", True,
            families=["mural_contractile_support"],
        )
        smooth = score_row(
            "contractile", "smooth_muscle", True,
            families=["contractile_structure_support"],
        )
        self.assertFalse(candidate_core_seed(pericyte, CATALOG["pericyte_mural"]))
        self.assertFalse(candidate_core_seed(smooth, CATALOG["smooth_muscle"]))

    def test_separable_endothelial_pericyte_mixture_releases_two_subsets(self) -> None:
        rows: list[dict[str, str]] = []
        endothelial = [f"e{i}" for i in range(6)]
        pericyte = [f"p{i}" for i in range(6)]
        for cell in endothelial + pericyte:
            rows.append(score_row(
                cell, "vascular_endothelial", cell in endothelial
            ))
            rows.append(score_row(cell, "pericyte_mural", cell in pericyte))
        index, _, universe = index_scores(rows)
        endothelial_result = validate_subset(
            endothelial, "vascular_endothelial", index, universe,
            catalog=CATALOG,
        )
        pericyte_result = validate_subset(
            pericyte, "pericyte_mural", index, universe,
            catalog=CATALOG,
        )
        self.assertEqual(endothelial_result["status"], "PASS")
        self.assertEqual(pericyte_result["status"], "PASS")
        self.assertEqual(
            resolve_overlap(
                endothelial[0],
                ["vascular_endothelial", "pericyte_mural"],
                index,
                CATALOG,
            )[0],
            "vascular_endothelial",
        )

    def test_inseparable_complete_endothelial_pericyte_cellbin_is_unresolved(self) -> None:
        rows = [
            score_row("mixed", "vascular_endothelial", True, score=0.90),
            score_row("mixed", "pericyte_mural", True, score=0.90),
        ]
        index, _, _ = index_scores(rows)
        chosen, reason = resolve_overlap(
            "mixed", ["vascular_endothelial", "pericyte_mural"],
            index, CATALOG,
        )
        self.assertEqual(chosen, "")
        self.assertEqual(reason, "unresolved_candidate_overlap")

    def test_lymphatic_identity_is_parent_locked_to_endothelial(self) -> None:
        candidate = CATALOG["lymphatic_endothelial"]
        self.assertEqual(candidate["candidate_role"], "fine")
        self.assertEqual(candidate["release_broad_label"], "Endothelial")
        self.assertEqual(candidate["parent_broad_label"], "Endothelial")
        self.assertEqual(
            candidate["release_fine_label"], "Lymphatic endothelial"
        )
        self.assertEqual(
            set(candidate["required_positive_families"]),
            {"endothelial_backbone", "lymphatic_identity_support"},
        )

    def test_theca_uses_observation_unit_specific_vascular_competitors(self) -> None:
        theca = CATALOG["theca_steroidogenic"]
        by_unit = theca["hard_anti_families_by_observation_unit"]
        self.assertNotIn("vascular_endothelial", by_unit["cellbin"])
        self.assertNotIn("pericyte_mural", by_unit["spot"])
        self.assertIn("vascular_endothelial", by_unit["cell"])
        self.assertIn("pericyte_mural", by_unit["nucleus"])

        sparse_neighbor = score_row("t", "theca_steroidogenic", True)
        sparse_neighbor.update(
            hard_contradiction="true",
            direct_anti_gene_count="1",
            direct_anti_family_count="1",
        )
        self.assertFalse(hard_contradiction(sparse_neighbor))
        full_direct_competitor = dict(
            sparse_neighbor,
            direct_anti_gene_count="2",
            direct_anti_family_count="2",
        )
        self.assertTrue(hard_contradiction(full_direct_competitor))

    def test_true_endothelial_identity_beats_sparse_theca_signal(self) -> None:
        rows = [
            score_row("v", "vascular_endothelial", True, score=0.90),
            score_row(
                "v", "theca_steroidogenic", True, score=0.15,
                families=["steroidogenic_core"],
            ),
        ]
        index, _, _ = index_scores(rows)
        chosen, _ = resolve_overlap(
            "v", ["vascular_endothelial", "theca_steroidogenic"],
            index, CATALOG,
        )
        self.assertEqual(chosen, "vascular_endothelial")

    def test_legacy_broad_only_vascular_is_never_blindly_migrated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ledger = root / "legacy.tsv.gz"
            with gzip.open(ledger, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    ["cell_id", "final_broad_label", "final_fine_label", "final_state"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerow({
                    "cell_id": "legacy_1",
                    "final_broad_label": "Vascular-associated",
                    "final_fine_label": "",
                    "final_state": "defined_broad_only",
                })
            output = root / "migrated.tsv.gz"
            result = subprocess.run([
                sys.executable,
                str(SCRIPTS / "migrate_release_taxonomy_v2.py"),
                "--ledger", str(ledger),
                "--profile", str(PROFILE),
                "--catalog", str(CATALOG_PATH),
                "--out", str(output),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with gzip.open(output, "rt", encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["final_broad_label"], "")
            self.assertEqual(row["final_fine_label"], "")
            self.assertEqual(row["final_state"], "unresolved_biological")
            self.assertEqual(row["final_cell_type"], "QC/Unknown")
            self.assertEqual(
                row["taxonomy_migration_status"],
                "requires_source_subcluster_readjudication",
            )

    def test_marker_families_encode_identity_not_shared_contractility(self) -> None:
        profile = json.loads(PROFILE.read_text(encoding="utf-8"))["lineages"]
        endothelial = profile["blood_endothelial"]["positive_families"]
        mural = profile["pericyte_mural"]["positive_families"]
        smooth = profile["smooth_muscle"]["positive_families"]
        theca = profile["theca_steroidogenic"]["positive_families"]
        self.assertTrue(
            {"PECAM1", "CDH5", "VWF", "CLDN5", "ESAM"}
            <= set(endothelial["endothelial_junction_backbone"])
        )
        self.assertTrue(
            {"RGS5", "PDGFRB", "CSPG4", "NOTCH3"}
            <= set(mural["mural_identity_backbone"])
        )
        self.assertTrue(
            {"MYH11", "CNN1", "ACTG2", "SMTN", "LMOD1"}
            <= set(smooth["mature_contractile_core"])
        )
        self.assertTrue(
            {"CYP11A1", "STAR", "HSD3B1"}
            <= set(theca["steroidogenic_core"])
        )
        self.assertTrue(
            {"CYP17A1", "INSL3", "ANPEP"}
            <= set(theca["theca_androgenic_identity"])
        )
        self.assertTrue(
            {"NR5A1", "FDX1", "FDXR", "POR", "CYB5A"}
            <= set(theca["theca_metabolic_support"])
        )

    def test_legacy_final_builders_have_no_v22_release_authority(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "config").mkdir()
            (root / "config/project.json").write_text(
                json.dumps({"canonical_lineage_controller_version": "2.2.0"}),
                encoding="utf-8",
            )
            final = subprocess.run([
                sys.executable, str(SCRIPTS / "build_final_annotation.py"),
                str(root), "--cell-ledger", str(root / "missing.tsv"),
                "--out", str(root / "out.tsv"), "--sample", "S",
            ], capture_output=True, text=True)
            report = subprocess.run([
                sys.executable, str(SCRIPTS / "build_report.py"), str(root),
            ], capture_output=True, text=True)
            self.assertNotEqual(final.returncode, 0)
            self.assertIn("no v2.2 release authority", final.stderr + final.stdout)
            self.assertNotEqual(report.returncode, 0)
            self.assertIn("retired for v2.2", report.stderr + report.stdout)

    def test_release_evidence_writes_single_public_cell_type_census(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            membership = root / "membership.tsv.gz"
            with gzip.open(membership, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    ["cell_id", "final_broad_label", "final_fine_label", "final_cell_type"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerows([
                    {"cell_id": "e", "final_broad_label": "Endothelial", "final_fine_label": "", "final_cell_type": "Endothelial"},
                    {"cell_id": "l", "final_broad_label": "Endothelial", "final_fine_label": "Lymphatic endothelial", "final_cell_type": "Lymphatic endothelial"},
                    {"cell_id": "q", "final_broad_label": "", "final_fine_label": "", "final_cell_type": "QC/Unknown"},
                ])
            out = root / "tables"
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "build_release_evidence_tables.py"),
                "--membership", str(membership), "--catalog", str(CATALOG_PATH),
                "--profile", str(PROFILE), "--out", str(out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (out / "final_cell_type_census.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                counts = {
                    row["final_cell_type"]: int(row["n_observations"])
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(
                counts,
                {"Endothelial": 1, "Lymphatic endothelial": 1, "QC/Unknown": 1},
            )

    def test_master_review_requires_final_cell_type_assets(self) -> None:
        source = (SCRIPTS / "request_master_quality_review.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('asset_manifest.get("label_column") != "final_cell_type"', source)
        self.assertNotIn('asset_manifest.get("label_column") != "primary_broad_label"', source)


if __name__ == "__main__":
    unittest.main()
