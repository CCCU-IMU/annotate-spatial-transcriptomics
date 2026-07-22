#!/usr/bin/env python3
from __future__ import annotations

import ast
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "annotate-spatial-transcriptomics"
TEXT_SUFFIXES = {".md", ".py", ".R", ".r", ".json", ".yaml", ".yml", ".tsv", ".txt", ".sh"}
BANNED = (
    "/" + "share" + "/" + "org" + "/",
    "bgi" + "_" + "zhangsch",
    "bgi" + "_" + "baiyy",
    "bgi" + "_" + "xinq",
    "D055" + "22A3",
)


def fail(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)


def main() -> int:
    problems: list[str] = []
    required = [
        ROOT / "README.md",
        ROOT / "requirements-ci.txt",
        ROOT / "install.sh",
        ROOT / "LICENSE",
        SKILL / "SKILL.md",
        SKILL / "references" / "direct-lineage-controller.md",
        SKILL / "references" / "legacy" / "multi-route-controller.md",
        SKILL / "references" / "matched-single-cell-reference.md",
        SKILL / "references" / "multi-sample-agent-orchestration.md",
        SKILL / "references" / "profiles" / "sheep_ovary_rfirst_profile.json",
        SKILL / "references" / "profiles" / "sheep_ovary_same_batch_rfirst_preset.json",
        SKILL / "references" / "profiles" / "sheep_ovary_candidate_lineage_catalog.json",
        SKILL / "references" / "profiles" / "sheep_ovary_literature_2025_2026.md",
        SKILL / "references" / "profiles" / "sheep_ovary_rfirst_case_reference.md",
        SKILL / "references" / "profiles" / "sheep_ovary_standard_workflow.md",
        SKILL / "references" / "efficient-operation.md",
        SKILL / "references" / "r-first-workflow.md",
        SKILL / "references" / "taxonomy-and-cohort-design.md",
        SKILL / "references" / "report-contract.md",
        SKILL / "scripts" / "audit_release_taxonomy.py",
        SKILL / "scripts" / "validate_matched_reference_crosswalk.py",
        SKILL / "scripts" / "resolve_workflow_profile.py",
        SKILL / "scripts" / "init_annotation_cohort.py",
        SKILL / "scripts" / "validate_cohort_state.py",
        SKILL / "scripts" / "migrate_project_v1_3_to_v1_4.py",
        SKILL / "scripts" / "migrate_project_v1_4_to_v1_5.py",
        SKILL / "scripts" / "migrate_project_v1_6_0_to_v1_6_1.py",
        SKILL / "scripts" / "migrate_project_v1_6_1_to_v1_7_0.py",
        SKILL / "scripts" / "migrate_project_v1_7_to_v1_8.py",
        SKILL / "scripts" / "validate_lineage_signal_coverage.py",
        SKILL / "scripts" / "register_derived_expression_artifact.py",
        SKILL / "scripts" / "validate_project_input_boundary.py",
        SKILL / "scripts" / "validate_banksy_broad_resolution_selection.py",
        SKILL / "scripts" / "validate_broad_class_completeness.py",
        SKILL / "scripts" / "screen_embedded_lineage_components.py",
        SKILL / "scripts" / "evaluate_graph_sensitivity.py",
        SKILL / "scripts" / "validate_oocyte_context_boundary.py",
        SKILL / "scripts" / "build_final_annotation.py",
        SKILL / "scripts" / "validate_profile_role.py",
        SKILL / "scripts" / "validate_resolution_grid.py",
        SKILL / "scripts" / "rank_cohort_resolutions.py",
        SKILL / "scripts" / "validate_cohort_outcome.py",
        SKILL / "scripts" / "validate_direct_return_evidence.py",
        SKILL / "scripts" / "validate_annotation_support_registry.py",
        SKILL / "scripts" / "validate_prelabel_broad_evidence.py",
        SKILL / "scripts" / "lineage_decision_lib.py",
        SKILL / "scripts" / "build_global_atlas_concordance.py",
        SKILL / "scripts" / "route_global_atlas_v2.py",
        SKILL / "scripts" / "validate_global_atlas_v2.py",
        SKILL / "scripts" / "build_annotation_contract_v2.py",
        SKILL / "scripts" / "validate_annotation_contract_v2.py",
        SKILL / "scripts" / "run_broad_family_evidence.R",
        SKILL / "scripts" / "validate_broad_family_evidence.py",
        SKILL / "scripts" / "audit_project_results_v2.py",
        SKILL / "scripts" / "migrate_project_v1_10_to_v2.py",
        SKILL / "scripts" / "migrate_release_taxonomy_v2.py",
        SKILL / "scripts" / "validate_global_atlas_concordance.py",
        SKILL / "scripts" / "validate_atlas_index_manifest.py",
        SKILL / "scripts" / "audit_annotation_membership_partition.py",
        SKILL / "scripts" / "controller_step.py",
        SKILL / "scripts" / "validate_benchmark_isolation.py",
        SKILL / "scripts" / "run_biological_benchmark.py",
        SKILL / "schemas" / "cohort_outcome.schema.json",
        SKILL / "schemas" / "direct_return_evidence.schema.json",
        SKILL / "schemas" / "annotation_support.schema.json",
        SKILL / "schemas" / "prelabel_broad_evidence.schema.json",
        SKILL / "schemas" / "atlas_index_manifest.schema.json",
        SKILL / "schemas" / "atlas_discrepancy_decision.schema.json",
        SKILL / "schemas" / "annotation_contract.schema.json",
        SKILL / "schemas" / "broad_family_evidence.schema.json",
        SKILL / "schemas" / "residual_qc_audit.schema.json",
        SKILL / "scripts" / "validate_direct_lineage_workflow.py",
        SKILL / "scripts" / "register_recluster_cohort.py",
        SKILL / "scripts" / "register_direct_return.py",
        SKILL / "scripts" / "preflight_generated_job.py",
        SKILL / "scripts" / "validate_incident_registry.py",
        SKILL / "scripts" / "request_cohort_confirmation.py",
        SKILL / "scripts" / "record_cohort_confirmation.py",
        SKILL / "scripts" / "master_quality_lib.py",
        SKILL / "scripts" / "request_master_quality_review.py",
        SKILL / "scripts" / "record_master_quality_approval.py",
        SKILL / "scripts" / "validate_open_world_lineage_audit.py",
        SKILL / "scripts" / "init_open_world_lineage_audit.py",
        SKILL / "scripts" / "prepare_seurat_full_feature_validation.R",
        SKILL / "scripts" / "seurat_validation_layer.R",
        SKILL / "scripts" / "run_initial_cluster_evidence.R",
        SKILL / "scripts" / "run_final_label_deg.R",
        SKILL / "scripts" / "compare_clusterings.py",
        SKILL / "scripts" / "scheduler_job_name.py",
        SKILL / "scripts" / "build_confirmation_review.py",
        SKILL / "assets" / "matched_reference_crosswalk_template.tsv",
    ]
    for path in required:
        if not path.is_file():
            problems.append(f"missing required file: {path.relative_to(ROOT)}")

    for path in sorted(ROOT.rglob("*")):
        if (
            not path.is_file()
            or ".git" in path.parts
            or ".release_extract" in path.parts
            or "dist" in path.parts
        ):
            continue
        if "__pycache__" in path.parts or path.suffix == ".pyc":
            problems.append(f"generated cache file present: {path.relative_to(ROOT)}")
            continue
        if path.suffix in TEXT_SUFFIXES or path.name in {"LICENSE", "VERSION"}:
            try:
                text = path.read_text(encoding="utf-8")
            except UnicodeDecodeError as exc:
                problems.append(f"not UTF-8 text: {path.relative_to(ROOT)}: {exc}")
                continue
            for token in BANNED:
                if token in text:
                    problems.append(f"private/sample token {token!r} in {path.relative_to(ROOT)}")
            if path.suffix == ".py":
                try:
                    ast.parse(text, filename=str(path))
                except SyntaxError as exc:
                    problems.append(f"Python syntax: {path.relative_to(ROOT)}: {exc}")
            elif path.suffix == ".json":
                try:
                    json.loads(text)
                except json.JSONDecodeError as exc:
                    problems.append(f"JSON syntax: {path.relative_to(ROOT)}: {exc}")

    regression = ROOT / "tests/test_release_contract.py"
    if not regression.is_file():
        problems.append("missing release-contract regression tests")
    if not (ROOT / "tests/test_v1_6_contract.py").is_file():
        problems.append("missing v1.6 active-contract regression tests")
    if not (ROOT / "tests/test_v1_7_contract.py").is_file():
        problems.append("missing v1.7 active-contract regression tests")
    if not (ROOT / "tests/test_v1_8_contract.py").is_file():
        problems.append("missing v1.8 active-contract regression tests")
    if not (ROOT / "tests/test_v1_9_contract.py").is_file():
        problems.append("missing v1.9 broad-completeness regression tests")
    if not (ROOT / "tests/test_legacy_migration_contract.py").is_file():
        problems.append("missing legacy migration regression tests")
    if not (ROOT / "tests/test_v2_contract.py").is_file():
        problems.append("missing v2 contract regression tests")
    if not (ROOT / "tests/test_v2_0_2_contract.py").is_file():
        problems.append("missing v2.0.2 biological-writeback regression tests")

    version = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    release_workflow = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    if f"当前版本：`{version}`" not in readme:
        problems.append("README current version differs from VERSION")
    if f"VERSION: v{version}" not in release_workflow:
        problems.append("release workflow version differs from VERSION")

    if problems:
        for problem in problems:
            fail(problem)
        return 1
    print("PASS: repository structure, portability scan, Python syntax and JSON validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
