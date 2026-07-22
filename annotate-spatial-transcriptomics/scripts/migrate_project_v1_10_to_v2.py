#!/usr/bin/env python3
"""Upgrade project controls to v2 without grandfathering biological results."""

from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


REQUIRED_REVALIDATION = [
    "build_and_validate_annotation_contract_v2",
    "rerun_full_feature_broad_family_evidence_matrix",
    "reselect_or_revalidate_banksy_whole_tissue_resolution_against_bound_upstream_grid",
    "reroute_all_cell_atlas_with_route_global_atlas_v2",
    "migrate_and_validate_release_broad_fine_hierarchy",
    "rerun_typed_residual_qc_audit_if_triggered",
    "run_read_only_project_results_audit",
    "rerun_completion_and_master_quality_gates",
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("project_root", type=Path)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    root = args.project_root.resolve()
    path = root / "config/project.json"
    project = json.loads(path.read_text(encoding="utf-8"))
    source = str(project.get("framework_version", ""))
    if source != "1.10.0":
        raise SystemExit(f"migration requires framework_version 1.10.0, found {source or '<empty>'}")
    result = {
        "status": "BLOCKED_PENDING_V2_REVALIDATION",
        "source_framework_version": source,
        "target_framework_version": "2.0.0",
        "existing_biological_labels_grandfathered": False,
        "existing_completion_gate_grandfathered": False,
        "required_revalidation": REQUIRED_REVALIDATION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    if args.dry_run:
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    backup = root / "provenance/migrations/v1_10_to_v2"
    backup.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, backup / "project.v1.10.0.json")
    project.update({
        "framework_version": "2.0.0",
        "annotation_contract_required": True,
        "broad_family_evidence_matrix_required": True,
        "atlas_authoritative_router": "route_global_atlas_v2.py",
        "release_taxonomy_hierarchy_validation_required": True,
        "project_results_readonly_audit_required": True,
        "migration_status": "blocked_pending_v2_revalidation",
    })
    path.write_text(json.dumps(project, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (backup / "migration_manifest.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    completion = root / "provenance/completion_gate.json"
    if completion.is_file():
        shutil.copy2(completion, backup / "completion_gate.v1.10.0.json")
        completion.unlink()
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
