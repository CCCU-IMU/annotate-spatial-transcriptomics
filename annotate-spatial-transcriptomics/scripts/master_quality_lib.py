#!/usr/bin/env python3
"""Shared hash-bound validation for the post-annotation master quality gate."""

from __future__ import annotations

import hashlib
import json
import csv
import io
from datetime import datetime
from pathlib import Path


BOUND_ARTIFACTS = (
    "cell_ledger",
    "cluster_ledger",
    "completion_gate",
    "release_taxonomy_audit",
    "annotation_workflow_audit",
    "state_validation",
    "iteration_plan",
    "next_action_queue",
    "recluster_cohort_registry",
    "direct_return_registry",
    "route_attempt_registry",
    "run_registry",
    "final_annotation_census",
    "annotation_view_registry",
    "annotation_support_registry",
    "annotation_support_validation",
    "annotation_membership_partition_audit",
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
    # The biological approval freezes the run registry as it existed at review
    # time. Post-confirmation release-only runs are expected to append to that
    # registry and must not invalidate the biological snapshot. Verify the
    # exact reviewed prefix byte-for-byte, then permit only terminal release
    # stages after the recorded review timestamp.
    errors.extend(validate_bound_record(root, approval, artifacts=tuple(key for key in BOUND_ARTIFACTS if key != "run_registry")))
    run_path = resolve(root, str(approval.get("run_registry", "")))
    expected_run_hash = str(approval.get("run_registry_sha256", ""))
    if not run_path.is_file():
        errors.append("bound artifact is missing: run_registry")
    elif sha256(run_path) != expected_run_hash:
        try:
            reviewed_at = datetime.fromisoformat(str(approval.get("reviewed_at", "")))
            with run_path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fields = list(reader.fieldnames or [])
                rows = list(reader)
            reviewed_rows = []
            appended_rows = []
            for row in rows:
                started = datetime.fromisoformat(row.get("started_at", ""))
                (reviewed_rows if started <= reviewed_at else appended_rows).append(row)
            reviewed_hashes = set()
            for lineterminator in ("\r\n", "\n"):
                buffer = io.StringIO(newline="")
                writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator=lineterminator)
                writer.writeheader(); writer.writerows(reviewed_rows)
                reviewed_hashes.add(hashlib.sha256(buffer.getvalue().encode("utf-8")).hexdigest())
            terminal = {"validated_done", "failed_preserved", "cancelled_preserved"}
            allowed_prefixes = ("final_", "release_", "report_")
            if expected_run_hash not in reviewed_hashes:
                errors.append("bound artifact changed after master review: run_registry reviewed prefix")
            if not appended_rows or any(not row.get("stage", "").startswith(allowed_prefixes) or row.get("status") not in terminal for row in appended_rows):
                errors.append("run_registry contains non-terminal or non-release rows after master review")
        except (ValueError, OSError, csv.Error):
            errors.append("bound artifact changed after master review: run_registry")
    checklist = approval.get("checklist", {})
    for item in CHECKLIST_IDS:
        entry = checklist.get(item, {}) if isinstance(checklist, dict) else {}
        if entry.get("status") not in {"PASS", "CONCERN"}:
            errors.append(f"master quality checklist is incomplete: {item}")
    return not errors, errors, approval
