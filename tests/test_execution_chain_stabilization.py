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
PACKAGE = ROOT / "annotate-spatial-transcriptomics"
SCRIPTS = PACKAGE / "scripts"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "sha256": sha(path),
        "n_bytes": path.stat().st_size,
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(rows[0])
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def membership_row(
    cell_id: str, broad: str, state: str, candidate: str, origin: str,
) -> dict[str, object]:
    return {
        "cell_id": cell_id,
        "source_boundary": "cohort_1",
        "source_cluster": "0",
        "final_state": state,
        "final_broad_label": broad,
        "final_fine_label": "",
        "final_cell_type": broad,
        "candidate_id": candidate,
        "state_annotations": "",
        "confidence": "high" if broad else "",
        "assignment_origin": origin,
        "qc_reason": "",
        "unresolved_reason": "" if broad else "unresolved",
    }


class ExecutionChainStabilizationTests(unittest.TestCase):
    def test_canonical_sct_banksy_bootstrap_builds_one_valid_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rds = root / "sample.rds"
            rds.write_bytes(b"mock-rds-for-bootstrap-contract")
            fake_rscript = root / "Rscript"
            fake_rscript.write_text(
                """#!/usr/bin/env python3
import csv,gzip,json,sys
from pathlib import Path
out=Path(sys.argv[-1]);out.mkdir(parents=True,exist_ok=True)
def write(path,fields,rows):
 op=gzip.open if str(path).endswith('.gz') else open
 with op(path,'wt',encoding='utf-8',newline='') as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter='\\t',lineterminator='\\n');w.writeheader();w.writerows(rows)
cells=['c1','c2']
write(out/'analysis_membership.tsv.gz',['cell_id','x','y','analysis_scope'],[{'cell_id':c,'x':i,'y':i,'analysis_scope':'analysis_set'} for i,c in enumerate(cells)])
write(out/'excluded_initial_qc.tsv.gz',['cell_id','x','y','analysis_scope','exclusion_reason'],[])
write(out/'analysis_scope.tsv.gz',['cell_id','x','y','analysis_scope'],[{'cell_id':c,'x':i,'y':i,'analysis_scope':'analysis_set'} for i,c in enumerate(cells)])
rows=[]
for res in (0.2,0.4,0.6):
 rows += [{'cell_id':c,'boundary_id':'whole_tissue','resolution':res,'cluster':str(i),'resolution_role':'grid'} for i,c in enumerate(cells)]
write(out/'partition_grid.tsv.gz',['cell_id','boundary_id','resolution','cluster','resolution_role'],rows)
(out/'whole_tissue_grid.json').write_text(json.dumps({'candidate_resolutions':[0.2,0.4,0.6]}))
(out/'input_audit_manifest.json').write_text(json.dumps({
 'status':'PASS','sample_id':sys.argv[-2],
 'input_sha256':Path(sys.argv[3]).read_text().split()[0],
 'raw_count_assay':'RNA','coordinate_columns':['x','y'],
 'scoring_exports_historical_columns':False,
 'expression_boundary':'project-local non-SCT raw counts'
}))
""",
                encoding="utf-8",
            )
            fake_rscript.chmod(0o755)
            project = root / "project"
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "bootstrap_sct_banksy_project.py"),
                "--sample", "s1", "--rds", str(rds),
                "--project-root", str(project), "--rscript", str(fake_rscript),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            manifest = json.loads(
                (project / "provenance/bootstrap_manifest.json").read_text()
            )
            self.assertEqual(manifest["status"], "PASS")
            contract = json.loads(
                (project / "config/annotation_contract.json").read_text()
            )
            self.assertEqual(contract["input_scope"]["analysis_set_n"], 2)
            self.assertEqual(
                contract["atlas_routing"]["bundle_id"],
                "sheep_ovary_GSE233801_split_wall_v2",
            )

    def test_membership_transform_chain_accepts_legal_nonfixed_order(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            base = root / "base.tsv.gz"
            write_tsv(base, [
                membership_row("c1", "", "unresolved_biological", "", "merge"),
                membership_row("c2", "Stromal/mesenchymal", "defined_broad_only", "stromal", "merge"),
            ])
            atlas = root / "atlas.tsv.gz"
            write_tsv(atlas, [
                membership_row("c1", "Granulosa", "defined_broad_only", "granulosa", "post_merge_atlas_unlabeled_broad_rescue"),
                membership_row("c2", "Stromal/mesenchymal", "defined_broad_only", "stromal", "merge"),
            ])
            atlas_manifest = root / "atlas.json"
            atlas_manifest.write_text(json.dumps({
                "status": "PASS", "stage": "atlas_and_completeness_review",
                "membership": artifact(atlas),
            }), encoding="utf-8")
            cell_review = root / "cell_review.tsv.gz"
            write_tsv(cell_review, [
                membership_row("c1", "Granulosa", "defined_broad_only", "granulosa", "post_merge_atlas_unlabeled_broad_rescue"),
                membership_row("c2", "Theca", "defined_broad_only", "theca", "catalog_wide_lineage_review_round_1"),
            ])
            cell_manifest = root / "cell_review.json"
            cell_manifest.write_text(json.dumps({
                "status": "PASS_REQUIRES_NEXT_REVIEW_ROUND",
                "stage": "catalog_wide_lineage_review_apply",
                "formal_batch_closure_performed": False,
                "active_cell_type": "Theca",
                "membership": artifact(cell_review),
            }), encoding="utf-8")
            roi = root / "roi.tsv.gz"
            write_tsv(roi, [
                membership_row("c1", "", "unresolved_biological", "", "follicle_roi_review"),
                membership_row("c2", "Theca", "defined_broad_only", "theca", "catalog_wide_lineage_review_round_1"),
            ])
            roi_manifest = root / "roi.json"
            roi_manifest.write_text(json.dumps({
                "status": "PENDING_POST_REPAIR_BIOLOGICAL_REVIEW",
                "stage": "follicle_roi_repair_apply",
                "repaired_membership": artifact(roi),
            }), encoding="utf-8")

            script = SCRIPTS / "manage_membership_transform_chain.py"
            init = root / "t0"
            result = subprocess.run([
                sys.executable, str(script), "init", "--membership", str(base),
                "--out", str(init),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            chain = init / "membership_transform_chain.json"
            transitions = [
                ("atlas_unlabeled_broad_rescue", base, atlas, atlas_manifest, ""),
                ("cell_type_review_patch", atlas, cell_review, cell_manifest, "Theca"),
                # ROI review is intentionally after a cell-type review. The
                # ledger validates provenance, not one brittle hard-coded order.
                ("follicle_roi_reconciliation", cell_review, roi, roi_manifest, ""),
            ]
            for index, (operation, source, target, evidence, cell_type) in enumerate(transitions, 1):
                out = root / f"t{index}"
                command = [
                    sys.executable, str(script), "append", "--chain", str(chain),
                    "--operation", operation, "--source", str(source),
                    "--result", str(target), "--evidence-manifest", str(evidence),
                    "--out", str(out),
                ]
                if cell_type:
                    command.extend(["--target-cell-type", cell_type])
                result = subprocess.run(command, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                chain = out / "membership_transform_chain.json"
            document = json.loads(chain.read_text())
            self.assertEqual(document["transform_n"], 3)
            self.assertEqual(
                [row["operation"] for row in document["transforms"]],
                [row[0] for row in transitions],
            )

    def test_identity_neutral_source_sync_cannot_change_a_label(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "source.tsv"
            result_membership = root / "result.tsv"
            write_tsv(source, [
                membership_row("c1", "Granulosa", "defined_broad_only", "granulosa", "source")
            ])
            write_tsv(result_membership, [
                membership_row("c1", "Theca", "defined_broad_only", "theca", "sync")
            ])
            evidence = root / "sync.json"
            evidence.write_text(json.dumps({"status": "PASS", "stage": "source_sync"}))
            script = SCRIPTS / "manage_membership_transform_chain.py"
            init = root / "init"
            subprocess.run([
                sys.executable, str(script), "init", "--membership", str(source),
                "--out", str(init),
            ], check=True, capture_output=True, text=True)
            attempted = subprocess.run([
                sys.executable, str(script), "append",
                "--chain", str(init / "membership_transform_chain.json"),
                "--operation", "source_unit_sync", "--source", str(source),
                "--result", str(result_membership),
                "--evidence-manifest", str(evidence), "--out", str(root / "bad"),
            ], capture_output=True, text=True)
            self.assertNotEqual(attempted.returncode, 0)
            self.assertIn("identity-neutral", attempted.stderr + attempted.stdout)

    def test_sequential_queue_exposes_exactly_one_active_cell_type(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            queue = root / "queue.tsv"
            rows = [
                {
                    "review_id": f"r{index}", "review_mode": "broad_lineage_review",
                    "target_broad_label": label,
                    "unit_signature": hashlib.sha256(label.encode()).hexdigest(),
                }
                for index, label in enumerate(
                    ["Granulosa", "Theca", "Endothelial"], 1
                )
            ]
            write_tsv(queue, rows)
            review = root / "review.json"
            review.write_text(json.dumps({
                "stage": "post_atlas_catalog_wide_lineage_review",
                "membership": artifact(queue),
                "artifacts": {"review_queue": artifact(queue)},
            }), encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(SCRIPTS / "manage_cell_type_review_queue.py"),
                "--review-manifest", str(review), "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
            state = json.loads((root / "out/cell_type_review_state.json").read_text())
            self.assertEqual(state["active_review_n"], 1)
            self.assertEqual(state["active_cell_type_review"]["target_broad_label"], "Granulosa")
            self.assertEqual(state["queued_review_n"], 2)
            self.assertTrue(state["formal_batch_closure_forbidden"])

    def test_fixed_atlas_uses_independent_split_wall_prototypes(self) -> None:
        descriptor = (
            PACKAGE
            / "references/atlases/sheep_ovary_GSE233801_split_wall_v2.json"
        )
        document = json.loads(descriptor.read_text())
        for label in (
            "Granulosa", "Stromal/mesenchymal", "Smooth muscle",
            "Endothelial", "Pericyte/mural", "Immune",
        ):
            self.assertEqual(
                document["release_broad_capability"][label], "supported"
            )
        self.assertEqual(
            document["release_broad_capability"]["Theca"], "challenge_only"
        )
        self.assertEqual(
            document["release_broad_capability"]["Epithelial/mesothelial"],
            "challenge_only",
        )
        self.assertFalse(document["legacy_source_labels_may_be_released"])

        crosswalk = read_tsv(
            PACKAGE
            / "references/atlases/GSE233801_res0p4_split_wall_crosswalk_v2.tsv"
        )
        by_cluster = {row["reference_cluster"]: row for row in crosswalk}
        self.assertEqual(by_cluster["8"]["framework_broad_label"], "Smooth muscle")
        self.assertEqual(by_cluster["12"]["framework_broad_label"], "Pericyte/mural")
        self.assertEqual(by_cluster["1"]["framework_broad_label"], "Endothelial")
        self.assertEqual(by_cluster["4"]["include_in_prototype"], "FALSE")

    def test_fixed_atlas_runtime_assets_are_distributed_with_skill(self) -> None:
        descriptor = (
            PACKAGE
            / "references/atlases/sheep_ovary_GSE233801_split_wall_v2.json"
        )
        document = json.loads(descriptor.read_text())
        asset_root = (
            PACKAGE
            / "references/atlases/sheep_ovary_GSE233801_split_wall_v2_assets"
        )
        self.assertEqual(
            set(document["asset_hashes"]),
            {
                "fixed_features.tsv",
                "feature_transform.joblib",
                "reference_prototypes.npz",
                "reference_heldout_predictions.tsv.gz",
                "reference_cluster_crosswalk.tsv",
                "reference_split.tsv",
            },
        )
        for name, expected in document["asset_hashes"].items():
            path = asset_root / name
            self.assertTrue(path.is_file(), name)
            self.assertEqual(sha(path), expected, name)

        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "atlas_validation.json"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS / "validate_fixed_atlas_bundle.py"),
                    "--bundle-manifest", str(descriptor),
                    "--asset-root", str(asset_root),
                    "--out", str(output),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(json.loads(output.read_text())["status"], "PASS")

    def test_merged_atlas_is_resume_only(self) -> None:
        descriptor = PACKAGE / "references/atlases/sheep_ovary_GSE233801_v1.json"
        document = json.loads(descriptor.read_text())
        self.assertTrue(document["legacy_resume_only"])
        self.assertFalse(document["new_routing_authority"])

    def test_runtime_dependency_registry_is_the_only_bound_script_list(self) -> None:
        sys.path.insert(0, str(SCRIPTS))
        from runtime_dependency_registry import CANONICAL_SCRIPTS  # noqa: PLC0415

        self.assertIn("manage_cell_type_review_queue.py", CANONICAL_SCRIPTS)
        self.assertIn("manage_membership_transform_chain.py", CANONICAL_SCRIPTS)
        self.assertIn("bootstrap_sct_banksy_project.py", CANONICAL_SCRIPTS)
        self.assertIn("freeze_sct_banksy_input.R", CANONICAL_SCRIPTS)
        for name in CANONICAL_SCRIPTS:
            self.assertTrue((SCRIPTS / name).is_file(), name)
        for consumer in (
            "build_annotation_contract_v2.py",
            "validate_annotation_contract_v2.py",
            "run_lineage_controller.py",
        ):
            source = (SCRIPTS / consumer).read_text(encoding="utf-8")
            self.assertIn("runtime_dependency_registry", source)

    def test_final_rds_writer_keeps_pigz_binary_stdout_out_of_r_memory(self) -> None:
        source = (SCRIPTS / "write_frozen_annotations_to_seurat.R").read_text(
            encoding="utf-8"
        )
        self.assertIn("stdout=compressed_tmp,stderr=pigz_error_tmp", source)
        self.assertNotIn("stdout=compressed_tmp,stderr=TRUE", source)

    def test_absolute_dotplot_can_bind_a_separate_raw_count_assay(self) -> None:
        source = (SCRIPTS / "build_marker_dotplots.R").read_text(
            encoding="utf-8"
        )
        self.assertIn('arg$`count-assay`', source)
        self.assertIn("safe_layer(count_assay, count_layer)", source)

    def test_controller_runtime_state_distinguishes_pause_done_and_failure(self) -> None:
        sys.path.insert(0, str(SCRIPTS))
        from controller_runtime_state import materialize_runtime_state  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            contract = root / "contract.json"
            contract.write_text("{}", encoding="utf-8")
            _, pause = materialize_runtime_state(
                root / "pause", "atlas_and_completeness_review", contract,
                result={
                    "status": "REVIEW_REQUIRED",
                    "required_progress_message": "现在开始对 Theca 进行专项复核。",
                },
            )
            self.assertTrue(json.loads(pause.read_text())["safe_scheduler_exit"])
            _, done = materialize_runtime_state(
                root / "done", "materialize_final_release", contract,
                result={"status": "PASS"},
            )
            self.assertEqual(json.loads(done.read_text())["status"], "DONE_PENDING_USER_REVIEW")
            _, failed = materialize_runtime_state(
                root / "failed", "cluster_cohort_recluster", contract,
                error=RuntimeError("boom"),
            )
            self.assertFalse(json.loads(failed.read_text())["safe_scheduler_exit"])

    def test_final_release_materialization_is_an_ordered_transform(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "source.tsv.gz"
            source_rows = [
                membership_row("c1", "", "unresolved_biological", "", "review"),
                membership_row("c2", "Granulosa", "defined_broad_only", "granulosa", "review"),
            ]
            source_rows[0]["final_cell_type"] = ""
            source_rows[1]["final_cell_type"] = ""
            write_tsv(source, source_rows)
            final = root / "final.tsv.gz"
            final_rows = [dict(row) for row in source_rows]
            final_rows[0].update({
                "final_state": "qc_holdout", "final_cell_type": "QC/Unknown",
                "assignment_origin": "final_typed_residual_qc",
            })
            final_rows[1]["final_cell_type"] = "Granulosa"
            write_tsv(final, final_rows)
            evidence = root / "final_release_manifest.json"
            evidence.write_text(json.dumps({
                "stage": "materialize_final_release", "status": "PASS",
                "membership": artifact(final),
            }), encoding="utf-8")
            manager = SCRIPTS / "manage_membership_transform_chain.py"
            subprocess.run([
                sys.executable, str(manager), "init", "--membership", str(source),
                "--out", str(root / "init"),
            ], check=True, capture_output=True, text=True)
            result = subprocess.run([
                sys.executable, str(manager), "append",
                "--chain", str(root / "init/membership_transform_chain.json"),
                "--operation", "final_release_materialization",
                "--source", str(source), "--result", str(final),
                "--evidence-manifest", str(evidence), "--out", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            chain = json.loads(
                (root / "out/membership_transform_chain.json").read_text()
            )
            self.assertEqual(chain["transforms"][-1]["operation"], "final_release_materialization")

    def test_final_deliverables_resume_without_repeating_completed_steps(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            membership = root / "membership.tsv.gz"
            write_tsv(membership, [
                membership_row("c1", "Granulosa", "defined_broad_only", "granulosa", "review")
            ])
            excluded = root / "excluded.tsv"
            write_tsv(excluded, [{"cell_id": "x"}])
            selected = root / "input.rds"
            selected.write_bytes(b"mock")
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({
                "candidate_boundaries": [{
                    "candidate_id": "granulosa", "candidate_role": "broad",
                    "release_broad_label": "Granulosa", "specificity_priority": 1,
                }]
            }), encoding="utf-8")
            profile = root / "profile.json"
            profile.write_text("{}", encoding="utf-8")
            prerequisite = root / "atlas_review.json"
            prerequisite.write_text(json.dumps({"status": "PASS"}), encoding="utf-8")
            release = root / "release.json"
            release.write_text(json.dumps({
                "status": "PASS", "n_analysis_set": 1, "residual_qc_n": 0,
                "residual_qc_fraction": 0.0,
                "final_cell_type_census": {"Granulosa": 1},
                "state_census": {},
                "membership": {
                    **artifact(membership), "semantic_sha256": "a" * 64,
                },
                "atlas_completeness_review": artifact(prerequisite),
            }), encoding="utf-8")
            contract = root / "contract.json"
            contract.write_text(json.dumps({
                "sample_id": "s1", "observation_unit": "cellbin",
                "selected_input_snapshot": artifact(selected),
                "input_scope": {
                    "analysis_set": artifact(membership),
                    "excluded_initial_qc": artifact(excluded),
                },
                "biological_profile": artifact(profile),
                "candidate_catalog": artifact(catalog),
                "canonical_lineage_controller": {"random_seed": 2200},
            }), encoding="utf-8")
            calls = root / "calls.txt"
            fake = root / "Rscript"
            fake.write_text(
                """#!/usr/bin/env python3
import csv,json,sys
from pathlib import Path
script=Path(sys.argv[1]).name; args=sys.argv[2:]
def val(x): return Path(args[args.index(x)+1])
with open(sys.argv[0]+'.calls','a') as h: h.write(script+'\\n')
out=val('--out'); out.parent.mkdir(parents=True,exist_ok=True)
if script=='write_frozen_annotations_to_seurat.R':
 out.write_bytes(b'annotated'); m=val('--writer-manifest');m.parent.mkdir(parents=True,exist_ok=True);m.write_text(json.dumps({'status':'PASS'}))
elif script=='build_annotation_maps.R':
 (out/'figures').mkdir(parents=True,exist_ok=True);(out/'tables').mkdir(parents=True,exist_ok=True);(out/'figures/final_cell_type_spatial.png').write_bytes(b'png');(out/'tables/spatial_node_asset_index.tsv').write_text('level\\tlabel\\tpng\\ncell_type\\tGranulosa\\t\\n')
elif script=='build_marker_dotplots.R':
 out.mkdir(parents=True,exist_ok=True);(out/'marker_dotplot_asset_index.tsv').write_text('level\\tpanel\\tpng\\tabsolute_png\\ncell_type\\tcanonical\\t\\t\\n')
elif script=='build_spatial_gene_maps.R':
 (out/'tables').mkdir(parents=True,exist_ok=True);(out/'tables/spatial_gene_group_asset_index.tsv').write_text('marker_group\\tavailable_genes\\tpng\\nGranulosa\\t\\t\\n')
elif script=='run_final_label_deg.R':
 (out/'tables').mkdir(parents=True,exist_ok=True);(out/'tables/cell_type_DEG_one_vs_rest_all.tsv').write_text('label\\tgene\\tavg_log2FC\\tpct_expressed_absolute\\nGranulosa\\tFOXL2\\t1\\t50\\n')
""",
                encoding="utf-8",
            )
            fake.chmod(0o755)
            # The fake writes its call log next to its executable.
            command = [
                sys.executable, str(SCRIPTS / "materialize_final_deliverables.py"),
                "--contract", str(contract), "--membership", str(membership),
                "--release-manifest", str(release), "--rscript", str(fake),
                "--out", str(root / "delivery"),
            ]
            first = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            call_log = Path(str(fake) + ".calls")
            first_calls = call_log.read_text()
            second = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
            self.assertEqual(call_log.read_text(), first_calls)
            checkpoint = json.loads(
                (root / "delivery/final_deliverables_checkpoint.json").read_text()
            )
            self.assertEqual(checkpoint["status"], "PASS")
            self.assertEqual(checkpoint["public_annotation_column"], "final_cell_type")


if __name__ == "__main__":
    unittest.main()
