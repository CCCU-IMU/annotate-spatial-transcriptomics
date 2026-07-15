#!/usr/bin/env python3
"""Shared hash-bound validation for the post-annotation master quality gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


BOUND_ARTIFACTS = (
    "cell_ledger",
    "cluster_ledger",
    "completion_gate",
    "release_taxonomy_audit",
    "multiroute_audit",
    "state_validation",
    "iteration_plan",
    "next_action_queue",
    "pool_registry",
    "route_attempt_registry",
    "run_registry",
    "final_annotation_census",
    "annotation_view_registry",
    "annotation_support_registry",
    "review_asset_manifest",
    "quality_reference",
)


CHECKLIST_IDS = (
    "broad_annotation_reasonableness",
    "marker_antimarker_and_spatial_support",
    "context_specific_and_confounded_lineage_safety",
    "comparable_to_validated_rfirst_reference",
)

OPTIONAL_BOUND_ARTIFACTS = (
    "strategy_preset_record",
    "open_world_lineage_audit",
    "open_world_lineage_audit_source",
    "candidate_lineage_catalog",
    "open_world_biological_profile",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def validate_bound_record(root: Path, record: dict, artifacts=BOUND_ARTIFACTS) -> list[str]:
    errors: list[str] = []
    selected = list(artifacts) + [key for key in OPTIONAL_BOUND_ARTIFACTS if record.get(key)]
    for key in selected:
        value = str(record.get(key, ""))
        expected = str(record.get(f"{key}_sha256", ""))
        if not value or not expected:
            errors.append(f"missing bound artifact field: {key}")
            continue
        path = resolve(root, value)
        if not path.is_file():
            errors.append(f"bound artifact is missing: {key}")
        elif sha256(path) != expected:
            errors.append(f"bound artifact changed after master review: {key}")
    return errors


def validate_master_approval(root: Path, approval_path: Path | None = None) -> tuple[bool, list[str], dict]:
    path = approval_path or root / "state/master_quality_approval.json"
    if not path.is_file():
        return False, ["master quality approval is missing"], {}
    try:
        approval = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False, ["master quality approval is unreadable"], {}
    errors: list[str] = []
    if approval.get("status") != "PASS":
        errors.append("master quality approval is not PASS")
    if approval.get("reviewer_role") != "main_conversation_agent":
        errors.append("quality approval was not recorded by the main conversation Agent")
    if approval.get("reviewed_after_all_annotation_complete") is not True:
        errors.append("quality approval did not verify the post-completion stage")
    request_path = resolve(root, str(approval.get("request_path", "")))
    if not request_path.is_file():
        errors.append("master quality review request is missing")
    elif sha256(request_path) != approval.get("request_sha256"):
        errors.append("master quality review request changed after approval")
    errors.extend(validate_bound_record(root, approval))
    checklist = approval.get("checklist", {})
    for item in CHECKLIST_IDS:
        entry = checklist.get(item, {}) if isinstance(checklist, dict) else {}
        if entry.get("status") not in {"PASS", "CONCERN"}:
            errors.append(f"master quality checklist is incomplete: {item}")
    return not errors, errors, approval
