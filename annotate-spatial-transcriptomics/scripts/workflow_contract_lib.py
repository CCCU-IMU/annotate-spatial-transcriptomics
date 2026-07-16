#!/usr/bin/env python3
"""Resolve active preprocessing/resolution contracts from verified provenance, never path hints."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


REQUIRED_INPUT_CHECKS = {
    "object_inspection_passed", "counts_data_layer_audit_passed",
    "profile_marker_coverage_passed", "cell_coordinate_consistency_passed",
    "conversion_provenance_verified",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def active_workflow_contract(root: Path) -> dict[str, Any]:
    project = load_json(root / "config/project.json")
    workflow = load_json(root / "config/active_workflow_profile.json")
    preset = load_json(root / "config/active_strategy_preset.json")
    provenance = load_json(root / "provenance/input_contract_validation.json")
    checks = provenance.get("checks", {}) if isinstance(provenance.get("checks"), dict) else {}
    checks_pass = all(checks.get(key) is True for key in REQUIRED_INPUT_CHECKS)
    workflow_profile_path = Path(str(workflow.get("workflow_profile", "")))
    workflow_binding_valid = (
        workflow_profile_path.is_file()
        and bool(workflow.get("workflow_profile_sha256"))
        and sha256(workflow_profile_path) == workflow.get("workflow_profile_sha256")
    )
    preset_bindings = preset.get("strategy_preset_bindings", {}) if isinstance(preset.get("strategy_preset_bindings"), dict) else {}
    preset_path = Path(str(preset_bindings.get("preset", "")))
    preset_binding_valid = (
        preset_path.is_file()
        and bool(preset_bindings.get("preset_sha256"))
        and sha256(preset_path) == preset_bindings.get("preset_sha256")
    )
    same_batch = (
        workflow.get("status") == "ACTIVE"
        and workflow_binding_valid
        and workflow.get("workflow_profile_id") == "seurat_stereopy_cellbin_same_batch_rfirst"
        and preset.get("strategy_preset_status") == "ACTIVE"
        and preset.get("strategy_preset_id") == project.get("strategy_preset_requested") == "sheep_ovary_same_batch_rfirst"
        and preset_binding_valid
        and provenance.get("status") == "PASS"
        and provenance.get("verified_stereopy_cellbin_same_batch") is True
        and checks_pass
    )
    grid = workflow.get("candidate_resolution_grid", [])
    if same_batch and grid != [0.1, 0.2, 0.3, 0.4, 0.6]:
        same_batch = False
    bound_profile = load_json(workflow_profile_path) if workflow_binding_valid else {}
    whole_tissue_grid = []
    cohort_grid = []
    if workflow_binding_valid:
        if same_batch:
            whole_tissue_grid = bound_profile.get("stereopy_cellbin_pped_contract", {}).get("clustering", {}).get("candidate_resolutions", [])
        else:
            whole_tissue_grid = bound_profile.get("whole_tissue_contract", {}).get("candidate_resolutions", [])
        cohort_grid = bound_profile.get("cohort_reclustering_contract", {}).get("candidate_resolutions", [])
    return {
        "same_batch_stereopy_cellbin_rfirst": same_batch,
        "candidate_resolution_grid": grid if isinstance(grid, list) else [],
        "workflow_profile_id": workflow.get("workflow_profile_id"),
        "workflow_profile": str(workflow_profile_path) if workflow_profile_path else "",
        "workflow_profile_sha256": workflow.get("workflow_profile_sha256"),
        "whole_tissue_resolution_grid": whole_tissue_grid if isinstance(whole_tissue_grid, list) else [],
        "cohort_resolution_grid": cohort_grid if isinstance(cohort_grid, list) else [],
        "strategy_preset_id": preset.get("strategy_preset_id"),
        "input_contract_status": provenance.get("status"),
        "input_contract_checks_complete": checks_pass,
        "workflow_profile_binding_valid": workflow_binding_valid,
        "strategy_preset_binding_valid": preset_binding_valid,
        "path_and_feature_hints_are_nonbinding": True,
    }
