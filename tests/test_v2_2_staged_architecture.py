from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "annotate-spatial-transcriptomics/scripts"


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = fields or (list(rows[0]) if rows else ["cell_id"])
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(script: str, *args: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPTS / script), *map(str, args)],
        capture_output=True, text=True,
    )


CATALOG = {
    "candidate_boundaries": [
        {"candidate_id": "granulosa", "candidate_role": "broad", "release_broad_label": "Granulosa", "release_fine_label": ""},
        {"candidate_id": "stromal", "candidate_role": "broad", "release_broad_label": "Stromal/mesenchymal", "release_fine_label": ""},
        {"candidate_id": "smooth", "candidate_role": "broad", "release_broad_label": "Smooth muscle", "release_fine_label": ""},
        {"candidate_id": "pericyte", "candidate_role": "fine", "release_broad_label": "Vascular-associated", "release_fine_label": "Pericyte/mural"},
        {"candidate_id": "epithelial", "candidate_role": "broad", "release_broad_label": "Epithelial/mesothelial", "release_fine_label": ""},
        {"candidate_id": "oocyte", "candidate_role": "broad", "release_broad_label": "Oocyte", "release_fine_label": "", "writeback_strategy": "canonical_cluster_membership"},
        {"candidate_id": "granulosa_hypoxia", "candidate_role": "state", "release_broad_label": "", "release_fine_label": "", "release_state_label": "Hypoxia", "parent_broad_label": "Granulosa"},
    ]
}


def evidence(candidate: str, coherent: float = 0.0, seed: float = 0.0, deg: float = 0.0) -> dict:
    return {
        "resolution": "0.2", "resolution_role": "selected",
        "source_boundary": "cohort_1", "source_cluster": "0",
        "candidate_id": candidate, "n_observations": "20",
        "observation_seed_fraction": str(seed),
        "observation_identity_core_fraction": str(seed),
        "observation_identity_core_direct_fraction": str(seed),
        "observation_coherent_fraction": str(coherent),
        "observation_release_family_coherent_fraction": str(coherent),
        "hard_contradiction_fraction": "0",
        "mean_program_score": "0.4" if coherent else "0",
        "available_positive_gene_count": "4",
        "available_positive_family_count": "2",
        "group_positive_family_supported_count": "2" if coherent else "0",
        "group_positive_family_mean_fraction": str(coherent),
        "positive_marker_detection_fraction": "0.8" if coherent else "0",
        "positive_marker_pseudobulk_sum": "100" if coherent else "0",
        "marker_deg_log2fc_mean": str(deg),
        "anti_marker_detection_fraction": "0",
        "anti_marker_pseudobulk_sum": "0",
        "anti_marker_deg_log2fc_mean": "0",
        "spatial_local_support_fraction": "0.8" if coherent else "0",
        "spatial_group_connectivity_fraction": "0.8" if coherent else "0",
        "cross_resolution_stable_fraction": "0.8" if coherent else "0",
    }


class V22StagedArchitectureTests(unittest.TestCase):
    def test_interestrus_context_enables_evaluation_without_label_authority(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / "context.tsv"
            result = run(
                "build_candidate_context_evidence.py",
                "--species", "sheep",
                "--tissue", "ovary",
                "--reproductive-stage", "interestrus",
                "--out", out,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            rows = read_tsv(out)
            self.assertEqual(rows[0]["candidate_id"], "luteal_steroidogenic")
            self.assertEqual(rows[0]["status"], "supported")
            self.assertEqual(rows[0]["evidence_scope"], "evaluation_permission_only")
            self.assertEqual(rows[0]["identity_writeback_authority"], "false")

    def test_luteal_lineage_origin_programs_are_not_hard_anti(self) -> None:
        catalog = json.loads((
            ROOT
            / "annotate-spatial-transcriptomics/references/profiles/"
            "sheep_ovary_candidate_lineage_catalog.json"
        ).read_text(encoding="utf-8"))
        luteal = next(
            row for row in catalog["candidate_boundaries"]
            if row["candidate_id"] == "luteal_steroidogenic"
        )
        self.assertNotIn("granulosa", luteal["hard_anti_families"])
        self.assertNotIn("theca_steroidogenic", luteal["hard_anti_families"])
        self.assertTrue(luteal["formal_context_evidence_required"])

    def test_post_merge_seeded_component_uses_luteal_anatomical_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            contract = root / "contract.json"
            contract.write_text("{}\n", encoding="utf-8")
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({"candidate_boundaries": [
                {
                    "candidate_id": "theca", "candidate_role": "broad",
                    "release_broad_label": "Theca",
                    "writeback_strategy": "candidate_local_then_exact_remainder",
                    "specificity_priority": 90,
                    "contextual_parent_overrides": [{
                        "context_broad_label": "Luteal",
                        "writeback_broad_label": "Luteal",
                        "minimum_component_neighbor_fraction": 0.5,
                    }],
                },
                {
                    "candidate_id": "luteal", "candidate_role": "broad",
                    "release_broad_label": "Luteal",
                    "writeback_strategy": "context_gated_candidate_local",
                    "specificity_priority": 90,
                },
                {
                    "candidate_id": "stromal", "candidate_role": "broad",
                    "release_broad_label": "Stromal/mesenchymal",
                    "writeback_strategy": "generic_exact_remainder_after_specific_lineages",
                    "specificity_priority": 10,
                },
            ]}), encoding="utf-8")
            membership = root / "membership.tsv.gz"
            members = []
            for index in range(6):
                members.append({
                    "cell_id": f"u{index}", "source_boundary": "cohort_1",
                    "source_cluster": "1", "candidate_id": "",
                    "final_state": "unresolved_biological", "final_broad_label": "",
                    "final_fine_label": "", "confidence": "",
                    "assignment_origin": "exact_remainder_unresolved",
                    "unresolved_reason": "irreducible_lineage_overlap",
                    "broad_frozen": "true", "fine_anchor_eligible": "false",
                })
            for index in range(20):
                members.append({
                    "cell_id": f"l{index}", "source_boundary": "cohort_1",
                    "source_cluster": "2", "candidate_id": "luteal",
                    "final_state": "defined_broad_only",
                    "final_broad_label": "Luteal",
                    "final_fine_label": "", "confidence": "high",
                    "assignment_origin": "second_round_whole_subcluster",
                    "unresolved_reason": "", "broad_frozen": "true",
                    "fine_anchor_eligible": "false",
                })
            write_tsv(membership, members)
            scores = root / "scores.tsv.gz"
            score_rows = []
            for index, member in enumerate(members):
                x = index % 7
                y = index // 7
                for candidate, label, evidence_value, direct, core, contradiction in (
                    ("theca", "Theca", 0.8 if member["cell_id"].startswith("u") else 0.1, 0.8, True, False),
                    ("luteal", "Luteal", 0.2, 0.1, False, False),
                    ("stromal", "Stromal/mesenchymal", 0.1, 0.1, False, True),
                ):
                    score_rows.append({
                        "cell_id": member["cell_id"], "source_boundary": "cohort_1",
                        "source_cluster": member["source_cluster"],
                        "candidate_id": candidate, "candidate_role": "broad",
                        "release_broad_label": label, "direct_signal": direct,
                        "program_score": evidence_value, "normalized_evidence": evidence_value,
                        "positive_family_count": 2, "family_coherent": core,
                        "identity_core_coherent": core, "identity_core_direct": core,
                        "hard_contradiction": contradiction, "technical_flag": False,
                        "x": x, "y": y,
                    })
            write_tsv(scores, score_rows)
            cluster_evidence = root / "cluster_evidence.tsv"
            evidence_rows = []
            for candidate, coherent, seed, deg in (
                ("theca", 0.8, 0.8, 1.0),
                ("luteal", 0.8, 0.8, 1.0),
                ("stromal", 0.0, 0.0, 0.0),
            ):
                row = evidence(candidate, coherent=coherent, seed=seed, deg=deg)
                row["source_cluster"] = "1"
                evidence_rows.append(row)
            write_tsv(cluster_evidence, evidence_rows)
            authority = root / "authority.json"
            authority.write_text(json.dumps({
                "mode": "stage_authority",
                "phase": "atlas_and_completeness_review",
                "annotation_contract_sha256": hashlib.sha256(contract.read_bytes()).hexdigest(),
                "post_atlas_membership": {
                    "path": str(membership.resolve()),
                    "sha256": hashlib.sha256(membership.read_bytes()).hexdigest(),
                },
                "observation_scores": [{
                    "path": str(scores.resolve()),
                    "sha256": hashlib.sha256(scores.read_bytes()).hexdigest(),
                }],
                "cluster_evidence": [{
                    "path": str(cluster_evidence.resolve()),
                    "sha256": hashlib.sha256(cluster_evidence.read_bytes()).hexdigest(),
                }],
            }), encoding="utf-8")
            result = run(
                "review_post_merge_unresolved_components.py",
                "--contract", contract, "--stage-authority", authority,
                "--membership", membership, "--catalog", catalog,
                "--scores", scores, "--cluster-evidence", cluster_evidence,
                "--workers", 2, "--out", root / "out",
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            reviewed = read_tsv(root / "out/reviewed_post_atlas_broad_membership.tsv.gz")
            unresolved_rows = [row for row in reviewed if row["cell_id"].startswith("u")]
            self.assertEqual({row["final_broad_label"] for row in unresolved_rows}, {"Luteal"})
            self.assertTrue(all("contextual_parent_return" in row["assignment_origin"] for row in unresolved_rows))
            defined_rows = [row for row in reviewed if row["cell_id"].startswith("l")]
            self.assertEqual({row["final_broad_label"] for row in defined_rows}, {"Luteal"})

            # Anatomical neighborhood alone cannot relabel a Theca-origin
            # component as Luteal.  Remove selected-subcluster Luteal support
            # and verify that the same component remains biologically unresolved.
            evidence_rows[1] = evidence("luteal", coherent=0.0, seed=0.0, deg=0.0)
            evidence_rows[1]["source_cluster"] = "1"
            write_tsv(cluster_evidence, evidence_rows)
            authority_doc = json.loads(authority.read_text(encoding="utf-8"))
            authority_doc["cluster_evidence"][0]["sha256"] = hashlib.sha256(
                cluster_evidence.read_bytes()
            ).hexdigest()
            authority.write_text(json.dumps(authority_doc), encoding="utf-8")
            unsupported = run(
                "review_post_merge_unresolved_components.py",
                "--contract", contract, "--stage-authority", authority,
                "--membership", membership, "--catalog", catalog,
                "--scores", scores, "--cluster-evidence", cluster_evidence,
                "--workers", 2, "--out", root / "unsupported",
            )
            self.assertEqual(
                unsupported.returncode, 0, unsupported.stdout + unsupported.stderr
            )
            unsupported_rows = read_tsv(
                root / "unsupported/reviewed_post_atlas_broad_membership.tsv.gz"
            )
            unresolved_rows = [
                row for row in unsupported_rows if row["cell_id"].startswith("u")
            ]
            self.assertEqual({row["final_broad_label"] for row in unresolved_rows}, {""})
            context_audit = read_tsv(
                root / "unsupported/contextual_parent_override_audit.tsv"
            )
            self.assertEqual(
                {row["status"] for row in context_audit},
                {"PARENT_SOURCE_UNSUPPORTED"},
            )

    def test_local_split_inherits_bound_context_evidence(self) -> None:
        source = (SCRIPTS / "run_lineage_controller.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('context_record = trigger.get("context_evidence")', source)
        self.assertIn('subset_command.extend(["--context-evidence"', source)
        self.assertIn('closure_command.extend(["--context-evidence"', source)
        self.assertIn('contract.get("candidate_context_evidence")', source)

    def test_first_pass_builds_grid_evidence_inside_bound_controller(self) -> None:
        source = (SCRIPTS / "run_lineage_controller.py").read_text(encoding="utf-8")
        self.assertNotIn('whole.add_argument("--resolution-grid-evidence"', source)
        self.assertIn('grid_scoring = output / "00_full_grid_scoring"', source)
        self.assertIn("materialize_selected_grid_evidence", source)
        self.assertIn("WHOLE_TISSUE_FORK_WORKER_CAP = 64", source)
        self.assertNotIn('scoring = output / "03_provisional_scoring"', source)
        self.assertIn('paths["whole_tissue_partitions"]', source)

    def test_recluster_cache_fingerprint_excludes_annotation_semantics(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            membership = root / "membership.tsv"
            write_tsv(membership, [{"cell_id": "c1"}, {"cell_id": "c2"}])
            sys.path.insert(0, str(SCRIPTS))
            try:
                spec = importlib.util.spec_from_file_location(
                    "v22_controller_cache", SCRIPTS / "run_lineage_controller.py"
                )
                module = importlib.util.module_from_spec(spec)
                assert spec and spec.loader
                spec.loader.exec_module(module)
            finally:
                sys.path.pop(0)
            args = SimpleNamespace(
                membership=membership, cohort_id="initial_cluster__1",
                source_initial_cluster="1", resolution_contract="sheep_ovary",
                seed=2200,
            )
            contract = {
                "selected_input_snapshot": {"sha256": "a" * 64},
                "candidate_catalog": {"sha256": "catalog_semantics_are_excluded"},
            }
            paths = {
                "run_seurat_cohort_recluster.R": (
                    SCRIPTS / "run_seurat_cohort_recluster.R"
                ),
                "run_seurat_cohort_recluster_impl.R": (
                    SCRIPTS / "run_seurat_cohort_recluster_impl.R"
                ),
            }
            first, payload = module.recluster_cache_fingerprint(
                contract, paths, args, [0.1, 0.2, 0.3]
            )
            contract["candidate_catalog"]["sha256"] = "changed_semantics"
            semantic_only, _ = module.recluster_cache_fingerprint(
                contract, paths, args, [0.1, 0.2, 0.3]
            )
            self.assertEqual(first, semantic_only)
            self.assertNotIn("candidate_catalog", payload)
            args.seed = 2201
            changed_seed, _ = module.recluster_cache_fingerprint(
                contract, paths, args, [0.1, 0.2, 0.3]
            )
            self.assertNotEqual(first, changed_seed)

    def test_local_second_extraction_reuses_components_not_whole_object_cells(self) -> None:
        source = (SCRIPTS / "close_exact_remainders.py").read_text(encoding="utf-8")
        self.assertNotIn("def second_round_candidates", source)
        self.assertIn("residual_components", source)
        self.assertIn("candidate_local_spatial_component", source)
        self.assertIn("It cannot invent a whole-object per-cell", source)

    def test_r_local_split_has_candidate_specific_canonical_component_route(self) -> None:
        source = (SCRIPTS / "derive_candidate_local_subsets.R").read_text(encoding="utf-8")
        self.assertIn("canonical_cluster_challenger_supported", source)
        self.assertIn('"canonical_identity_component"', source)
        self.assertIn("spatial_fraction >= 0.40", source)
        self.assertIn("evidence$program_score_delta >= 0.05", source)
        self.assertIn("evidence$direct_signal_delta >= 0.05", source)

    def test_r_local_component_seeds_exclude_direct_hard_contradictions(self) -> None:
        source = (SCRIPTS / "derive_candidate_local_subsets.R").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            "core_mask <- identity_core_mask(candidate_scores, candidate_id) &",
            source,
        )
        self.assertIn("!candidate_scores$hard_block", source)
        self.assertIn(
            "Those observations\n      # remain in the exact remainder",
            source,
        )

    def test_r_fine_subset_requires_exact_parent_broad_identity(self) -> None:
        source = (SCRIPTS / "derive_candidate_local_subsets.R").read_text(
            encoding="utf-8"
        )
        self.assertIn("context_evidence_candidate_id", source)
        self.assertIn("parent_identity_status", source)
        self.assertIn("parent_family_floor_pass && parent_support >= 0.25", source)
        self.assertIn('parent_identity_status != "FAIL"', source)

    def test_canonical_zero_census_is_deferred_only_to_required_biological_review(self) -> None:
        audit = (SCRIPTS / "audit_post_merge_completeness.py").read_text(
            encoding="utf-8"
        )
        controller = (SCRIPTS / "run_lineage_controller.py").read_text(
            encoding="utf-8"
        )
        flag = "--defer-canonical-zero-to-biological-review"
        self.assertIn(flag, audit)
        self.assertIn("canonical_zero_census_deferred_to_required_object_level_", audit)
        self.assertIn("canonical_cluster_membership", audit)
        self.assertIn("if quality_required", controller)
        self.assertGreaterEqual(controller.count(flag), 2)

    def test_r_scorer_accepts_one_cluster_partition_without_fake_deg(self) -> None:
        source = (SCRIPTS / "run_observation_lineage_scoring.R").read_text(
            encoding="utf-8"
        )
        self.assertIn("deg[rest_n < 1L] <- 0", source)
        self.assertIn("if (nlevels(group) < 2L) next", source)

    def test_r_scorer_atomically_closes_empty_gzip_tables(self) -> None:
        source = (SCRIPTS / "run_observation_lineage_scoring.R").read_text(
            encoding="utf-8"
        )
        self.assertIn("write_tsv_atomic <- function", source)
        self.assertIn("compressed && nrow(value) == 0L", source)
        self.assertIn('connection <- gzfile(temporary, open = "wt")', source)
        self.assertIn('system2("gzip", c("-t", shQuote(temporary)))', source)
        self.assertIn("file.rename(temporary, path)", source)
        self.assertIn(
            "write_tsv_atomic(\n  program_table,\n",
            source,
        )

    def test_specific_split_program_precedes_generic_remainder(self) -> None:
        source = (SCRIPTS / "adjudicate_second_round_subclusters.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("specific_split_worthy", source)
        self.assertIn("generic_positive", source)
        self.assertLess(
            source.index("specific_split_worthy[0]"),
            source.index("generic_positive[0]"),
        )
        self.assertIn("elif winner and specific_split_worthy", source)

    def test_first_pass_reuses_selected_grid_evidence_without_observation_scores(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            grid = root / "grid"
            tables = grid / "tables"
            tables.mkdir(parents=True)
            rows = []
            for resolution, role in (("0.2", "selected"), ("0.4", "grid"), ("0.6", "grid")):
                item = evidence("granulosa", 0.8, 0.7, 1.2)
                item.update(resolution=resolution, resolution_role=role)
                rows.append(item)
            write_tsv(
                tables / "cluster_candidate_multichannel_evidence.tsv.gz", rows
            )
            write_tsv(
                tables / "resolution_deg_coexpression_programs.tsv.gz",
                [
                    {
                        "program_id": f"p{index}", "resolution": resolution,
                        "source_boundary": "whole", "source_cluster": "0",
                        "n_observations": "20", "genes": "A;B",
                        "coexpressed_gene_count": "2", "mean_top_log2fc": "1",
                        "mean_detection_difference": "0.5",
                        "catalog_marker_overlap_fraction": "0",
                        "spatially_coherent": "true",
                        "excluded_program_classes": "",
                        "candidate_status": "unmodeled_program_seed",
                    }
                    for index, resolution in enumerate(("0.2", "0.4", "0.6"), 1)
                ],
            )
            (grid / "observation_scoring_manifest.json").write_text(
                json.dumps({"status": "PASS"})
            )
            sys.path.insert(0, str(SCRIPTS))
            try:
                spec = importlib.util.spec_from_file_location(
                    "v22_controller_grid_view_test",
                    SCRIPTS / "run_lineage_controller.py",
                )
                module = importlib.util.module_from_spec(spec)
                assert spec and spec.loader
                spec.loader.exec_module(module)
                selected, manifest_path = module.materialize_selected_grid_evidence(
                    grid,
                    {
                        "selected_resolution": 0.4,
                        "selected_and_neighbors": [0.4, 0.2, 0.6],
                    },
                    root / "selected",
                )
            finally:
                sys.path.pop(0)
            selected_rows = read_tsv(selected)
            self.assertEqual({row["resolution"] for row in selected_rows}, {"0.4"})
            self.assertEqual({row["resolution_role"] for row in selected_rows}, {"selected"})
            manifest = json.loads(manifest_path.read_text())
            self.assertFalse(manifest["observation_scores_written"])
            self.assertFalse(
                (root / "selected/tables/observation_lineage_scores.tsv.gz").exists()
            )

    def test_resolution_evidence_is_weighted_by_subcluster_identity_purity(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            scoring = root / "scoring"
            tables = scoring / "tables"
            tables.mkdir(parents=True)
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps(CATALOG))
            candidates = [
                row["candidate_id"] for row in CATALOG["candidate_boundaries"]
            ]
            rows = []
            layouts = {
                "0.1": [("0", 100, {"granulosa": (0.80, 0.70, 1.0), "smooth": (0.70, 0.60, 1.0)})],
                "0.2": [("0", 100, {"granulosa": (0.80, 0.70, 1.0), "smooth": (0.35, 0.20, 0.7)})],
                "0.4": [
                    ("0", 50, {"granulosa": (0.90, 0.80, 1.5)}),
                    ("1", 50, {"smooth": (0.90, 0.80, 1.5)}),
                ],
            }
            for resolution, groups in layouts.items():
                for cluster, n_obs, supported in groups:
                    for candidate in candidates:
                        coherent, seed, deg = supported.get(
                            candidate, (0.0, 0.0, 0.0)
                        )
                        row = evidence(candidate, coherent, seed, deg)
                        row.update({
                            "resolution": resolution,
                            "resolution_role": "grid",
                            "source_boundary": "fixture",
                            "source_cluster": cluster,
                            "n_observations": str(n_obs),
                        })
                        rows.append(row)
            evidence_path = tables / "cluster_candidate_multichannel_evidence.tsv.gz"
            write_tsv(evidence_path, rows)
            scorer = SCRIPTS / "run_observation_lineage_scoring.R"
            (scoring / "observation_scoring_manifest.json").write_text(
                json.dumps({
                    "status": "PASS",
                    "scorer": {
                        "path": str(scorer),
                        "sha256": hashlib.sha256(scorer.read_bytes()).hexdigest(),
                    },
                    "candidate_universe": candidates,
                })
            )
            out = root / "out"
            result = run(
                "build_resolution_grid_evidence.py",
                "--scoring-output", scoring,
                "--catalog", catalog,
                "--selection-purpose", "cohort_identity_resolution",
                "--out", out,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            grid = {
                row["resolution"]: row
                for row in read_tsv(out / "resolution_grid_evidence.tsv")
            }
            self.assertEqual(float(grid["0.1"]["directly_resolved_observation_fraction"]), 0)
            self.assertEqual(float(grid["0.1"]["mixed_observation_fraction"]), 1)
            self.assertEqual(float(grid["0.4"]["directly_resolved_observation_fraction"]), 1)
            self.assertEqual(float(grid["0.4"]["mixed_observation_fraction"]), 0)

    def test_second_round_query_is_exactly_bound_to_first_round_cohort(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            contract = root / "contract.json"
            contract.write_text("{}")
            membership = root / "cohort.tsv.gz"
            write_tsv(membership, [{"cell_id": "c1"}, {"cell_id": "c2"}])
            plan = root / "plan.tsv"
            write_tsv(plan, [{
                "cohort_id": "initial_cluster__0",
                "source_initial_cluster": "0",
                "provisional_status": "unknown",
                "provisional_broad_after_score_freeze": "",
                "n_observations": "2",
                "membership_path": str(membership),
                "membership_sha256": hashlib.sha256(membership.read_bytes()).hexdigest(),
            }])
            whole = root / "whole.json"
            whole.write_text(json.dumps({
                "status": "PASS", "controller_version": "2.2.0",
                "phase": "whole_tissue_partition",
                "annotation_contract": {
                    "sha256": hashlib.sha256(contract.read_bytes()).hexdigest()
                },
                "cohort_plan": {
                    "path": str(plan),
                    "sha256": hashlib.sha256(plan.read_bytes()).hexdigest(),
                },
            }))
            sys.path.insert(0, str(SCRIPTS))
            try:
                spec = importlib.util.spec_from_file_location(
                    "v22_controller_test", SCRIPTS / "run_lineage_controller.py"
                )
                module = importlib.util.module_from_spec(spec)
                assert spec and spec.loader
                spec.loader.exec_module(module)
            finally:
                sys.path.pop(0)
            args = SimpleNamespace(
                whole_manifest=whole, contract=contract,
                cohort_id="initial_cluster__0", source_initial_cluster="0",
                provisional_status="unknown", provisional_broad="",
                membership=membership,
            )
            _, row = module.validate_cohort_plan_binding(args)
            self.assertEqual(row["membership_sha256"], hashlib.sha256(membership.read_bytes()).hexdigest())
            altered = root / "altered.tsv.gz"
            write_tsv(altered, [{"cell_id": "c1"}, {"cell_id": "c2"}, {"cell_id": "c3"}])
            args.membership = altered
            with self.assertRaisesRegex(RuntimeError, "exact bound initial-cluster"):
                module.validate_cohort_plan_binding(args)

    def test_tiny_initial_cluster_is_recorded_as_underpowered_not_qc(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            membership = root / "cohort.tsv.gz"
            write_tsv(membership, [{"cell_id": "c1"}, {"cell_id": "c2"}])
            plan = root / "plan.tsv"
            write_tsv(plan, [{
                "cohort_id": "initial_cluster__tiny",
                "source_initial_cluster": "tiny",
                "provisional_status": "unknown",
                "provisional_broad_after_score_freeze": "",
                "n_observations": "2",
                "membership_path": str(membership),
                "membership_sha256": hashlib.sha256(membership.read_bytes()).hexdigest(),
            }])
            whole = root / "whole.json"
            whole.write_text(json.dumps({
                "status": "PASS", "controller_version": "2.2.0",
                "phase": "whole_tissue_partition",
                "formal_membership_written": False,
                "release_authority_written": False,
                "cohort_plan": {
                    "path": str(plan),
                    "sha256": hashlib.sha256(plan.read_bytes()).hexdigest(),
                },
            }))
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps(CATALOG))
            selected_input = root / "input.rds"
            selected_input.write_bytes(b"fixture")
            sys.path.insert(0, str(SCRIPTS))
            try:
                spec = importlib.util.spec_from_file_location(
                    "v22_controller_underpowered", SCRIPTS / "run_lineage_controller.py"
                )
                module = importlib.util.module_from_spec(spec)
                assert spec and spec.loader
                spec.loader.exec_module(module)
            finally:
                sys.path.pop(0)
            args = SimpleNamespace(
                out=root / "out", cohort_id="initial_cluster__tiny",
                source_initial_cluster="tiny", provisional_broad="",
                membership=membership, whole_manifest=whole,
                context_evidence=None,
            )
            contract = {
                "selected_input_snapshot": {"sha256": hashlib.sha256(selected_input.read_bytes()).hexdigest()},
                "query_reclustering": {"candidate_resolutions": [0.1, 0.2, 0.3]},
            }
            result = module.materialize_underpowered_cohort(
                args, contract,
                {
                    "catalog": catalog, "selected_input": selected_input,
                    "run_observation_lineage_scoring.R": SCRIPTS / "run_observation_lineage_scoring.R",
                },
                {"c1", "c2"},
            )
            self.assertEqual(result["cohort_status"], "UNDERPOWERED_NOT_EVALUABLE")
            unresolved = read_tsv(Path(result["base_candidate_membership"]["path"]))
            self.assertEqual({row["proposed_state"] for row in unresolved}, {"unresolved_biological"})
            self.assertTrue(all(not row["proposed_broad_label"] for row in unresolved))
            outcome = json.loads((root / "out/cohort_outcome.json").read_text())
            self.assertEqual(outcome["terminal_outcome"], "underpowered_not_evaluable")
            validation = run(
                "validate_cohort_outcome.py", root,
                root / "out/cohort_outcome.json",
            )
            self.assertEqual(validation.returncode, 0, validation.stdout + validation.stderr)

    def test_calibrated_atlas_can_rescue_an_underpowered_cohort(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            membership = root / "post_atlas.tsv.gz"
            write_tsv(membership, [
                {
                    "cell_id": cell, "final_broad_label": "Granulosa",
                    "assignment_origin": "post_merge_atlas_unlabeled_broad_rescue",
                }
                for cell in ("c1", "c2")
            ])
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps(CATALOG))
            evidence_path = root / "underpowered_evidence.tsv.gz"
            write_tsv(evidence_path, [
                {
                    "resolution": "",
                    "resolution_role": "underpowered_not_evaluable",
                    "source_boundary": "tiny", "source_cluster": "underpowered",
                    "candidate_id": row["candidate_id"],
                    "evaluation_status": "underpowered_not_evaluable",
                }
                for row in CATALOG["candidate_boundaries"]
            ])
            fine = root / "fine.tsv"
            write_tsv(fine, [{
                "parent_broad_label": "Vascular-associated",
                "candidate_id": "pericyte", "status": "not_evaluable",
            }])
            accepted = root / "accepted.tsv"
            write_tsv(accepted, [], ["program_id", "status"])
            unmodeled = root / "unmodeled.json"
            unmodeled.write_text(json.dumps({
                "accepted_program_n": 0,
                "accepted_programs": {
                    "path": str(accepted),
                    "sha256": hashlib.sha256(accepted.read_bytes()).hexdigest(),
                },
            }))
            out = root / "audit"
            result = run(
                "audit_post_merge_completeness.py",
                "--membership", membership, "--catalog", catalog,
                "--cluster-evidence", evidence_path,
                "--fine-audit", fine, "--unmodeled", unmodeled,
                "--out", out,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            manifest = json.loads((out / "post_merge_completeness_manifest.json").read_text())
            self.assertEqual(manifest["status"], "PASS")
            self.assertEqual(manifest["underpowered_cohort_n"], 1)

    def test_project_local_candidate_membership_cannot_freeze_broad(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            contract = root / "contract.json"
            contract.write_text(json.dumps({
                "canonical_lineage_controller": {"controller_version": "2.2.0"},
                "selected_input_snapshot": {"sha256": "a" * 64},
            }))
            analysis = root / "analysis.tsv"
            write_tsv(analysis, [{"cell_id": "c1"}])
            candidate = root / "candidate.tsv.gz"
            write_tsv(candidate, [{
                "cell_id": "c1", "source_boundary": "custom",
                "source_cluster": "0", "candidate_id": "granulosa",
                "proposed_state": "broad_candidate",
                "proposed_broad_label": "Granulosa", "confidence": "high",
                "assignment_origin": "project_local_scorer", "unresolved_reason": "",
            }])
            source = root / "project_local_source.json"
            source.write_text(json.dumps({
                "status": "PASS", "controller_version": "2.2.0",
                "phase": "cluster_cohort_recluster",
                "formal_membership_written": False,
                "annotation_contract": {
                    "sha256": hashlib.sha256(contract.read_bytes()).hexdigest()
                },
                "base_candidate_membership": {
                    "path": str(candidate),
                    "sha256": hashlib.sha256(candidate.read_bytes()).hexdigest(),
                },
            }))
            authority = root / "authority.json"
            authority.write_text(json.dumps({
                "mode": "stage_authority", "phase": "merge_and_freeze_broad",
                "annotation_contract_sha256": hashlib.sha256(contract.read_bytes()).hexdigest(),
                "candidate_source_manifests": [{
                    "path": str(source),
                    "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                }],
            }))
            result = run(
                "merge_and_freeze_broad_membership.py",
                "--contract", contract, "--stage-authority", authority,
                "--analysis-membership", analysis,
                "--candidate-membership", candidate,
                "--candidate-source-manifest", source,
                "--out", root / "out",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("cohort adjudication", result.stdout + result.stderr)

    def test_reviewed_unmodeled_program_can_close_without_forced_label(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            membership = root / "post_atlas.tsv.gz"
            write_tsv(membership, [{
                "cell_id": "c1", "final_broad_label": "Granulosa",
                "assignment_origin": "post_merge_atlas_unlabeled_broad_rescue",
            }])
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps(CATALOG), encoding="utf-8")
            evidence_path = root / "evidence.tsv.gz"
            write_tsv(
                evidence_path,
                [evidence(row["candidate_id"]) for row in CATALOG["candidate_boundaries"]],
            )
            program_dir = root / "initial_cluster__1/05_unmodeled"
            program = program_dir / "unmodeled_lineage_candidates.tsv"
            write_tsv(program, [{
                "program_id": "unmodeled_lineage_candidate_001",
                "resolutions": "0.1;0.2", "genes": "A;B;C",
                "status": "Unmodeled lineage candidate",
            }])
            manifest_path = program_dir / "unmodeled_discovery_manifest.json"
            manifest_path.write_text(json.dumps({
                "accepted_program_n": 1,
                "accepted_programs": {
                    "path": str(program),
                    "sha256": hashlib.sha256(program.read_bytes()).hexdigest(),
                },
            }), encoding="utf-8")
            review = root / "unmodeled_review.tsv"
            write_tsv(review, [{
                "source_boundary": "initial_cluster__1",
                "program_id": "unmodeled_lineage_candidate_001",
                "outcome": "insufficient_identity_program",
                "catalog_candidate_id": "",
                "rationale": "Stable expression state lacks an independent multigene identity core.",
            }])
            result = run(
                "audit_post_merge_completeness.py",
                "--membership", membership, "--catalog", catalog,
                "--cluster-evidence", evidence_path,
                "--unmodeled", manifest_path,
                "--unmodeled-review", review,
                "--out", root / "out",
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            audit = read_tsv(root / "out/unmodeled_lineage_audit.tsv")
            self.assertEqual(
                audit[0]["review_outcome"], "insufficient_identity_program"
            )

    def adjudicate(self, supported: dict[str, tuple[float, float, float]]) -> tuple[dict, list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        catalog = root / "catalog.json"
        catalog.write_text(json.dumps(CATALOG))
        contract = root / "contract.json"
        contract.write_text(json.dumps({
            "observation_writeback": {"policy": {
                "whole_subcluster_min_raw_two_family_supported_fraction": 0.40,
                "whole_subcluster_min_raw_two_family_margin": 0.20,
                "maximum_contradiction_fraction": 0.05,
            }}
        }))
        cells = [f"c{i:02d}" for i in range(20)]
        partitions = root / "partitions.tsv"
        write_tsv(partitions, [
            {"cell_id": cell, "resolution": "0.2", "cluster": "0", "resolution_role": "selected"}
            for cell in cells
        ])
        evidence_rows = []
        for candidate in [row["candidate_id"] for row in CATALOG["candidate_boundaries"]]:
            coherent, seed, deg = supported.get(candidate, (0.0, 0.0, 0.0))
            evidence_rows.append(evidence(candidate, coherent, seed, deg))
        evidence_path = root / "evidence.tsv"
        write_tsv(evidence_path, evidence_rows)
        scores = root / "scores.tsv"
        write_tsv(scores, [
            {"cell_id": cell, "candidate_id": candidate}
            for cell in cells
            for candidate in [row["candidate_id"] for row in CATALOG["candidate_boundaries"]]
        ])
        out = root / "out"
        result = run(
            "adjudicate_second_round_subclusters.py",
            "--partitions", partitions, "--cluster-evidence", evidence_path,
            "--scores", scores, "--catalog", catalog, "--contract", contract,
            "--cohort-id", "cohort_1", "--source-initial-cluster", "7",
            "--provisional-status", "unknown", "--out", out,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return (
            json.loads((out / "second_round_adjudication_manifest.json").read_text()),
            read_tsv(out / "base_candidate_membership.tsv.gz"),
            read_tsv(out / "pending_local_split_membership.tsv.gz"),
            read_tsv(out / "fine_candidate_proposals.tsv"),
            read_tsv(out / "state_annotation_proposals.tsv"),
        )

    def test_granulosa_rich_subcluster_inherits_sparse_tail(self) -> None:
        manifest, base, pending, fine, _ = self.adjudicate({"granulosa": (0.85, 0.70, 1.5)})
        self.assertEqual(len(base), 20)
        self.assertTrue(all(row["proposed_broad_label"] == "Granulosa" for row in base))
        self.assertFalse(pending)
        self.assertFalse(any(row["status"] == "supported" for row in fine))
        self.assertEqual(manifest["formal_membership_written"], False)

    def test_minor_specific_program_can_trigger_local_check_despite_negative_group_deg(self) -> None:
        manifest, base, pending, _, _ = self.adjudicate({
            "stromal": (0.85, 0.75, 0.4),
            "smooth": (0.35, 0.35, -0.3),
        })
        self.assertFalse(base)
        self.assertEqual(len(pending), 20)
        self.assertEqual(manifest["status"], "LOCAL_SPLIT_REQUIRED")

    def test_state_is_released_only_for_a_high_purity_resolved_parent(self) -> None:
        _, base, pending, _, states = self.adjudicate({
            "granulosa": (0.85, 0.70, 1.5),
            "granulosa_hypoxia": (0.60, 0.50, 1.0),
        })
        self.assertEqual(len(base), 20)
        self.assertFalse(pending)
        self.assertEqual(len(states), 1)
        self.assertEqual(states[0]["state_annotation"], "Hypoxia")
        self.assertEqual(
            states[0]["assignment_scope"],
            "whole_high_purity_second_round_subcluster",
        )

    def test_stromal_smooth_vascular_mixture_requires_local_split(self) -> None:
        manifest, base, pending, fine, state = self.adjudicate({
            "stromal": (0.70, 0.60, 1.0),
            "smooth": (0.35, 0.20, 1.2),
            "pericyte": (0.30, 0.18, 1.1),
        })
        self.assertFalse(base)
        self.assertEqual(len(pending), 20)
        self.assertEqual(manifest["status"], "LOCAL_SPLIT_REQUIRED")
        pericyte = next(row for row in fine if row["candidate_id"] == "pericyte")
        self.assertEqual(pericyte["status"], "not_evaluable")
        self.assertEqual(pericyte["release_candidate"], "false")
        self.assertFalse(state)

    def test_low_fraction_epithelial_program_is_not_silently_diluted(self) -> None:
        _, base, pending, _, _ = self.adjudicate({
            "stromal": (0.80, 0.70, 1.0),
            "epithelial": (0.05, 0.01, 2.0),
        })
        self.assertFalse(base)
        self.assertEqual(len(pending), 20)

    def test_oocyte_canonical_subcluster_returns_whole_membership(self) -> None:
        _, base, pending, _, _ = self.adjudicate({"oocyte": (0.60, 0.10, 2.5)})
        self.assertEqual(len(base), 20)
        self.assertTrue(all(row["proposed_broad_label"] == "Oocyte" for row in base))
        self.assertFalse(pending)

    def test_first_round_creates_only_one_cluster_one_cohort_plan(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps(CATALOG))
            partitions = root / "partitions.tsv"
            rows = []
            for role, resolution in (("selected", 0.4), ("neighbor_1", 0.2), ("neighbor_2", 0.6)):
                for index in range(8):
                    rows.append({
                        "cell_id": f"c{index}", "boundary_id": "whole",
                        "resolution": resolution, "cluster": str(index // 4),
                        "resolution_role": role,
                    })
            write_tsv(partitions, rows)
            ev = root / "evidence.tsv"
            ev_rows = []
            for cluster in ("0", "1"):
                for candidate in [row["candidate_id"] for row in CATALOG["candidate_boundaries"]]:
                    item = evidence(candidate, 0.8 if candidate == ("granulosa" if cluster == "0" else "stromal") else 0, 0.6, 1.0)
                    item.update(source_boundary="whole", source_cluster=cluster, resolution_role="selected")
                    ev_rows.append(item)
            write_tsv(ev, ev_rows)
            out = root / "out"
            result = run(
                "build_whole_tissue_cohort_plan.py", "--partitions", partitions,
                "--cluster-evidence", ev, "--catalog", catalog, "--out", out,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            plan = read_tsv(out / "whole_tissue_cohort_plan.tsv")
            self.assertEqual(len(plan), 2)
            self.assertTrue(all(row["formal_label_written"] == "false" for row in plan))
            self.assertFalse((out / "release_membership.tsv.gz").exists())

    def test_stable_unmodeled_program_blocks_automatic_naming(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            programs = root / "programs.tsv"
            write_tsv(programs, [
                {"program_id": "p1", "resolution": "0.2", "genes": "A;B;C", "candidate_status": "unmodeled_program_seed", "spatially_coherent": "true", "catalog_marker_overlap_fraction": "0", "coexpressed_gene_count": "3"},
                {"program_id": "p2", "resolution": "0.4", "genes": "A;B;D", "candidate_status": "unmodeled_program_seed", "spatially_coherent": "true", "catalog_marker_overlap_fraction": "0", "coexpressed_gene_count": "3"},
            ])
            out = root / "out"
            result = run("discover_unmodeled_lineages.py", "--programs", programs, "--out", out)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            candidates = read_tsv(out / "unmodeled_lineage_candidates.tsv")
            self.assertEqual(len(candidates), 1)
            self.assertEqual(candidates[0]["status"], "Unmodeled lineage candidate")
            self.assertFalse(any("label" in key.lower() for key in candidates[0] if key != "catalog_match"))

    def test_recurrent_state_program_is_not_promoted_to_unmodeled_lineage(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            programs = root / "programs.tsv"
            write_tsv(programs, [
                {"program_id": "h1", "resolution": "0.2", "genes": "HIF1A;VEGFA;BNIP3;CA9", "candidate_status": "unmodeled_program_seed", "spatially_coherent": "true", "catalog_marker_overlap_fraction": "0", "coexpressed_gene_count": "4"},
                {"program_id": "h2", "resolution": "0.4", "genes": "HIF1A;VEGFA;BNIP3;EGLN3", "candidate_status": "unmodeled_program_seed", "spatially_coherent": "true", "catalog_marker_overlap_fraction": "0", "coexpressed_gene_count": "4"},
            ])
            out = root / "out"
            result = run("discover_unmodeled_lineages.py", "--programs", programs, "--out", out)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertFalse(read_tsv(out / "unmodeled_lineage_candidates.tsv"))
            excluded = read_tsv(out / "excluded_state_or_technical_programs.tsv")
            self.assertEqual({row["excluded_program_class"] for row in excluded}, {"hypoxia"})

    def test_known_multifamily_program_is_not_promoted_to_unmodeled_lineage(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            programs = root / "programs.tsv"
            write_tsv(programs, [
                {
                    "program_id": "i1", "resolution": "0.2",
                    "genes": "C1QA;CD74;PTPRC;X", "candidate_status": "unmodeled_program_seed",
                    "spatially_coherent": "true", "catalog_marker_overlap_fraction": "0.1",
                    "coexpressed_gene_count": "4", "best_catalog_candidate_id": "immune",
                    "best_catalog_overlap_gene_count": "3",
                    "best_catalog_overlap_family_count": "2",
                },
                {
                    "program_id": "i2", "resolution": "0.4",
                    "genes": "C1QA;CD74;PTPRC;Y", "candidate_status": "unmodeled_program_seed",
                    "spatially_coherent": "true", "catalog_marker_overlap_fraction": "0.1",
                    "coexpressed_gene_count": "4", "best_catalog_candidate_id": "immune",
                    "best_catalog_overlap_gene_count": "3",
                    "best_catalog_overlap_family_count": "2",
                },
            ])
            out = root / "out"
            result = run(
                "discover_unmodeled_lineages.py", "--programs", programs,
                "--out", out,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertFalse(read_tsv(out / "unmodeled_lineage_candidates.tsv"))
            excluded = read_tsv(out / "excluded_state_or_technical_programs.tsv")
            self.assertEqual(
                {row["excluded_program_class"] for row in excluded},
                {"modeled_catalog_program"},
            )

    def test_final_unresolved_members_become_typed_qc_only_at_final_stage(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            contract = root / "contract.json"
            contract.write_text("{}")
            membership = root / "post_atlas.tsv"
            rows = []
            for index in range(20):
                broad = "" if index == 0 else "Granulosa"
                rows.append({
                    "cell_id": f"c{index}", "analysis_scope": "analysis_set",
                    "source_boundary": "cohort", "source_cluster": "0",
                    "candidate_id": "granulosa" if broad else "",
                    "final_state": "defined_broad_only" if broad else "unresolved_biological",
                    "final_broad_label": broad, "final_fine_label": "",
                    "confidence": "high" if broad else "low",
                    "assignment_origin": "second_round", "qc_reason": "",
                    "unresolved_reason": "ambiguous_biological_program" if not broad else "",
                    "fine_anchor_eligible": "false",
                })
            write_tsv(membership, rows)
            completeness = root / "completeness.json"
            completeness.write_text(json.dumps({
                "status": "PASS",
                "membership": {
                    "path": str(membership),
                    "sha256": hashlib.sha256(membership.read_bytes()).hexdigest(),
                },
            }))
            atlas_validation = root / "atlas_validation.json"
            atlas_validation.write_text(json.dumps({"status": "PASS"}))
            review = root / "review.json"
            review.write_text(json.dumps({
                "status": "PASS", "phase": "atlas_and_completeness_review",
                "membership": {"path": str(membership), "sha256": hashlib.sha256(membership.read_bytes()).hexdigest()},
                "completeness": {
                    "path": str(completeness),
                    "sha256": hashlib.sha256(completeness.read_bytes()).hexdigest(),
                },
                "atlas_validation": {
                    "path": str(atlas_validation),
                    "sha256": hashlib.sha256(atlas_validation.read_bytes()).hexdigest(),
                },
            }))
            authority = root / "authority.json"
            authority.write_text(json.dumps({
                "mode": "stage_authority", "phase": "materialize_final_release",
                "annotation_contract_sha256": hashlib.sha256(contract.read_bytes()).hexdigest(),
                "post_atlas_membership": {
                    "path": str(membership),
                    "sha256": hashlib.sha256(membership.read_bytes()).hexdigest(),
                },
                "prerequisite_manifest": {
                    "path": str(review),
                    "sha256": hashlib.sha256(review.read_bytes()).hexdigest(),
                },
                "state_annotation_proposals": [],
            }))
            out = root / "out"
            result = run(
                "materialize_final_release_v2_2.py", "--contract", contract,
                "--stage-authority", authority, "--post-atlas-membership", membership,
                "--atlas-completeness-manifest", review, "--out", out,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            final = read_tsv(out / "final_release_membership.tsv.gz")
            self.assertEqual(final[0]["final_state"], "qc_holdout")
            self.assertEqual(final[0]["qc_reason"], "ambiguous_biological_program")

    def test_completeness_requires_source_linked_candidate_support(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({"candidate_boundaries": [{
                "candidate_id": "granulosa", "candidate_role": "broad",
                "release_broad_label": "Granulosa",
                "release_fine_label": "", "required_positive_families": [],
            }]}))
            membership = root / "membership.tsv"
            write_tsv(membership, [
                {
                    "cell_id": "c0", "source_boundary": "cohort_a",
                    "source_cluster": "0", "candidate_id": "granulosa",
                    "final_broad_label": "Granulosa",
                    "assignment_origin": "second_round_whole_subcluster",
                },
                {
                    "cell_id": "c1", "source_boundary": "cohort_b",
                    "source_cluster": "0", "candidate_id": "",
                    "final_broad_label": "", "assignment_origin": "",
                },
            ])
            evidence_path = root / "evidence.tsv"
            evidence_rows = []
            for boundary, core in (("cohort_a", 0.0), ("cohort_b", 0.8)):
                evidence_rows.append({
                    "resolution_role": "selected", "source_boundary": boundary,
                    "source_cluster": "0", "candidate_id": "granulosa",
                    "available_positive_family_count": "2",
                    "group_positive_family_supported_count": "2",
                    "group_required_positive_families_pass": "true",
                    "observation_identity_core_fraction": str(core),
                    "observation_identity_core_direct_fraction": str(core),
                    "positive_marker_detection_fraction": "0.8",
                    "mean_program_score": "0.2",
                    "marker_deg_log2fc_mean": "1.0",
                    "anti_marker_deg_log2fc_mean": "0.0",
                    "positive_marker_pseudobulk_sum": "10",
                    "anti_marker_pseudobulk_sum": "0",
                    "cross_resolution_stable_fraction": "0.8",
                })
            write_tsv(evidence_path, evidence_rows)
            result = run(
                "audit_post_merge_completeness.py",
                "--membership", membership, "--catalog", catalog,
                "--cluster-evidence", evidence_path, "--out", root / "out",
            )
            self.assertEqual(result.returncode, 2)
            manifest = json.loads(
                (root / "out/post_merge_completeness_manifest.json").read_text()
            )
            self.assertEqual(manifest["status"], "BLOCKED")
            audit = read_tsv(root / "out/broad_completeness_audit.tsv")
            self.assertEqual(audit[0]["unsupported_release_observation_n"], "1")

    def test_completeness_accepts_bound_local_subset_support(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            membership = root / "membership.tsv.gz"
            write_tsv(membership, [{
                "cell_id": "c1", "final_broad_label": "Smooth muscle",
                "candidate_id": "smooth", "source_boundary": "cohort_1",
                "source_cluster": "0", "assignment_origin": "supported_subset",
            }])
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps(CATALOG), encoding="utf-8")
            cluster_evidence = root / "cluster_evidence.tsv.gz"
            write_tsv(
                cluster_evidence,
                [evidence(row["candidate_id"]) for row in CATALOG["candidate_boundaries"]],
            )
            local_validation = root / "local_subset_validation.tsv"
            write_tsv(local_validation, [{
                "source_boundary": "cohort_1", "source_cluster": "0",
                "candidate_id": "smooth", "status": "PASS",
            }])
            result = run(
                "audit_post_merge_completeness.py",
                "--membership", membership, "--catalog", catalog,
                "--cluster-evidence", cluster_evidence,
                "--local-subset-validation", local_validation,
                "--out", root / "out",
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            manifest = json.loads(
                (root / "out/post_merge_completeness_manifest.json").read_text()
            )
            self.assertEqual(manifest["status"], "PASS")

    def test_completeness_accepts_bound_post_merge_component_support(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            membership = root / "membership.tsv.gz"
            write_tsv(membership, [{
                "cell_id": "c1", "final_broad_label": "Smooth muscle",
                "candidate_id": "smooth", "source_boundary": "cohort_1",
                "source_cluster": "0",
                "assignment_origin": (
                    "post_merge_unresolved_component_review__"
                    "specific_unique_component"
                ),
            }])
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps(CATALOG), encoding="utf-8")
            cluster_evidence = root / "cluster_evidence.tsv.gz"
            write_tsv(
                cluster_evidence,
                [evidence(row["candidate_id"]) for row in CATALOG["candidate_boundaries"]],
            )
            decisions = root / "candidate_component_decisions.tsv"
            write_tsv(decisions, [{
                "cell_id": "c1", "source_boundary": "cohort_1",
                "source_cluster": "0", "candidate_id": "smooth",
                "release_broad_label": "Smooth muscle",
                "decision": "specific_unique_component",
                "component_id": "component_1",
            }])
            review = root / "post_merge_unresolved_review_manifest.json"
            review.write_text(json.dumps({
                "status": "PASS",
                "stage": "post_merge_unresolved_component_review",
                "membership": {
                    "path": str(membership.resolve()),
                    "sha256": hashlib.sha256(membership.read_bytes()).hexdigest(),
                },
                "component_artifacts": {
                    "candidate_component_decisions.tsv": {
                        "path": str(decisions.resolve()),
                        "sha256": hashlib.sha256(decisions.read_bytes()).hexdigest(),
                    },
                },
            }), encoding="utf-8")
            result = run(
                "audit_post_merge_completeness.py",
                "--membership", membership, "--catalog", catalog,
                "--cluster-evidence", cluster_evidence,
                "--post-merge-review-manifest", review,
                "--out", root / "out",
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            manifest = json.loads(
                (root / "out/post_merge_completeness_manifest.json").read_text()
            )
            self.assertEqual(manifest["status"], "PASS")

    def test_robustness_evaluator_requires_exact_repeat_and_biological_equivalence(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            baseline = root / "baseline.tsv.gz"
            replicate = root / "replicate.tsv.gz"
            variant = root / "variant.tsv.gz"
            rows = [
                {
                    "cell_id": f"c{index:03d}",
                    "final_broad_label": "Granulosa" if index < 60 else "Stromal/mesenchymal",
                    "final_fine_label": "Mural/estrogenic granulosa" if index < 40 else "",
                    "final_fine_candidate_id": "granulosa_mural_estrogenic" if index < 40 else "",
                    "final_state": "defined_broad_only",
                }
                for index in range(100)
            ]
            write_tsv(baseline, rows)
            write_tsv(replicate, rows)
            perturbed = [dict(row) for row in rows]
            perturbed[58]["final_broad_label"] = "Stromal/mesenchymal"
            perturbed[60]["final_broad_label"] = "Granulosa"
            write_tsv(variant, perturbed)
            out = root / "robustness"
            result = run(
                "evaluate_annotation_robustness.py",
                "--baseline", baseline,
                "--deterministic-replicate", replicate,
                "--variant", variant,
                "--out", out,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            manifest = json.loads((out / "parameter_robustness_manifest.json").read_text())
            self.assertEqual(manifest["technical_determinism"], "PASS")
            self.assertEqual(manifest["status"], "PASS")


if __name__ == "__main__":
    unittest.main()
