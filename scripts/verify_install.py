#!/usr/bin/env python3
from __future__ import annotations

import ast
import json
import sys
from pathlib import Path


REQUIRED = (
    "SKILL.md",
    "agents/openai.yaml",
    "references/context-and-biology.md",
    "references/iterative-controller.md",
    "references/direct-lineage-controller.md",
    "references/legacy/multi-route-controller.md",
    "references/matched-single-cell-reference.md",
    "references/multi-sample-agent-orchestration.md",
    "references/profiles/sheep_ovary_rfirst_profile.json",
    "references/profiles/sheep_ovary_same_batch_rfirst_preset.json",
    "references/profiles/sheep_ovary_candidate_lineage_catalog.json",
    "references/profiles/sheep_ovary_literature_2025_2026.md",
    "references/profiles/sheep_ovary_rfirst_case_reference.md",
    "references/r-first-workflow.md",
    "references/taxonomy-and-cohort-design.md",
    "references/quality-standard.md",
    "references/report-contract.md",
    "references/state-schema.md",
    "schemas/cohort_outcome.schema.json",
    "schemas/direct_return_evidence.schema.json",
    "schemas/annotation_support.schema.json",
    "schemas/prelabel_broad_evidence.schema.json",
    "schemas/atlas_index_manifest.schema.json",
    "schemas/atlas_discrepancy_decision.schema.json",
    "schemas/annotation_contract.schema.json",
    "schemas/broad_family_evidence.schema.json",
    "schemas/residual_qc_audit.schema.json",
    "scripts/discover_inputs.py",
    "scripts/init_annotation_project.py",
    "scripts/autopilot_status.py",
    "scripts/plan_next_iteration.py",
    "scripts/check_completion_gate.py",
    "scripts/validate_direct_lineage_workflow.py",
    "scripts/rank_cohort_resolutions.py",
    "scripts/validate_cohort_outcome.py",
    "scripts/validate_direct_return_evidence.py",
    "scripts/validate_annotation_support_registry.py",
    "scripts/validate_prelabel_broad_evidence.py",
    "scripts/build_global_atlas_concordance.py",
    "scripts/route_global_atlas_v2.py",
    "scripts/validate_global_atlas_v2.py",
    "scripts/build_annotation_contract_v2.py",
    "scripts/validate_annotation_contract_v2.py",
    "scripts/run_broad_family_evidence.R",
    "scripts/validate_broad_family_evidence.py",
    "scripts/audit_project_results_v2.py",
    "scripts/migrate_project_v1_10_to_v2.py",
    "scripts/migrate_release_taxonomy_v2.py",
    "scripts/validate_global_atlas_concordance.py",
    "scripts/validate_atlas_index_manifest.py",
    "scripts/audit_annotation_membership_partition.py",
    "scripts/controller_step.py",
    "scripts/register_recluster_cohort.py",
    "scripts/register_direct_return.py",
    "scripts/audit_release_taxonomy.py",
    "scripts/validate_matched_reference_crosswalk.py",
    "scripts/resolve_workflow_profile.py",
    "scripts/init_annotation_cohort.py",
    "scripts/validate_cohort_state.py",
    "scripts/migrate_project_v1_3_to_v1_4.py",
    "scripts/migrate_project_v1_6_0_to_v1_6_1.py",
    "scripts/migrate_project_v1_6_1_to_v1_7_0.py",
    "scripts/request_cohort_confirmation.py",
    "scripts/record_cohort_confirmation.py",
    "scripts/master_quality_lib.py",
    "scripts/request_master_quality_review.py",
    "scripts/record_master_quality_approval.py",
    "scripts/validate_open_world_lineage_audit.py",
    "scripts/init_open_world_lineage_audit.py",
    "scripts/prepare_seurat_full_feature_validation.R",
    "scripts/seurat_validation_layer.R",
    "scripts/run_initial_cluster_evidence.R",
    "scripts/run_final_label_deg.R",
    "scripts/compare_clusterings.py",
    "scripts/scheduler_job_name.py",
    "scripts/request_final_annotation_confirmation.py",
    "scripts/build_confirmation_review.py",
    "scripts/build_report.py",
    "scripts/audit_release.py",
)


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "annotate-spatial-transcriptomics").resolve()
    missing = [name for name in REQUIRED if not (root / name).is_file()]
    if missing:
        print("Installation is incomplete; missing:", file=sys.stderr)
        for name in missing:
            print(f"  - {name}", file=sys.stderr)
        return 1

    failures: list[str] = []
    for path in sorted((root / "scripts").glob("*.py")):
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except (SyntaxError, UnicodeDecodeError) as exc:
            failures.append(f"{path.relative_to(root)}: {exc}")
    for path in sorted(root.rglob("*.json")):
        try:
            json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            failures.append(f"{path.relative_to(root)}: {exc}")
    if failures:
        print("Installation validation failed:", file=sys.stderr)
        for item in failures:
            print(f"  - {item}", file=sys.stderr)
        return 1

    print(f"PASS: {root}")
    print("The Skill structure, Python syntax and JSON assets are valid.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
