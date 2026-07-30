#!/usr/bin/env python3
"""Materialize deterministic controller pause, resume and failure state."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from evidence_schema_lib import sha256


CONTROLLED_STATUSES = {
    "PASS", "REVIEW_REQUIRED", "ITERATION_REQUIRED",
    "PENDING_USER_REVIEW_HIGH_UNRESOLVED",
}


def atomic_json(path: Path, document: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(document, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        Path(name).replace(path)
    finally:
        Path(name).unlink(missing_ok=True)


def controller_status(phase: str, result_status: str) -> str:
    if result_status == "PENDING_USER_REVIEW_HIGH_UNRESOLVED":
        return "REVIEW_REQUIRED"
    if phase == "materialize_final_release" and result_status == "PASS":
        return "DONE_PENDING_USER_REVIEW"
    if result_status in CONTROLLED_STATUSES:
        return result_status
    return "FAILED_RUNTIME"


def default_next_action(phase: str, status: str) -> str:
    if status == "PASS":
        return "advance_to_next_canonical_phase"
    if status == "DONE_PENDING_USER_REVIEW":
        return "review_and_approve_frozen_annotation_report"
    if status == "REVIEW_REQUIRED":
        return "complete_the_single_controller_bound_review_action_then_resume"
    if status == "ITERATION_REQUIRED":
        return "complete_the_bounded_biological_iteration_then_resume"
    return "inspect_runtime_failure_and_resume_from_last_frozen_checkpoint"


def resume_token(document: dict[str, object]) -> str:
    stable = {
        key: document.get(key)
        for key in (
            "phase", "controller_status", "annotation_contract_sha256",
            "membership_sha256", "review_state_sha256", "next_action",
        )
    }
    return hashlib.sha256(
        json.dumps(stable, ensure_ascii=False, sort_keys=True).encode("utf-8")
    ).hexdigest()


def materialize_runtime_state(
    out: Path,
    phase: str,
    contract: Path,
    result: dict[str, object] | None = None,
    controller_manifest: Path | None = None,
    error: BaseException | None = None,
) -> tuple[Path, Path]:
    result = result or {}
    result_status = str(result.get("status", "FAILED_RUNTIME"))
    status = "FAILED_RUNTIME" if error else controller_status(phase, result_status)
    membership = result.get("membership", {})
    review_state = result.get("cell_type_review_state", {})
    next_action = str(
        result.get("required_next_action", "")
        or result.get("required_progress_message", "")
        or default_next_action(phase, status)
    )
    now = datetime.now(timezone.utc).isoformat()
    state: dict[str, object] = {
        "schema_version": "1.0",
        "phase": phase,
        "controller_status": status,
        "phase_result_status": result_status,
        "scheduler_exit_semantics": (
            "successful_controlled_pause_or_completion"
            if status != "FAILED_RUNTIME" else "runtime_failure"
        ),
        "annotation_contract": str(contract.resolve()),
        "annotation_contract_sha256": sha256(contract) if contract.is_file() else "",
        "membership_sha256": (
            str(membership.get("sha256", ""))
            if isinstance(membership, dict) else ""
        ),
        "review_state_sha256": (
            str(review_state.get("sha256", ""))
            if isinstance(review_state, dict) else ""
        ),
        "next_action": next_action,
        "updated_at_utc": now,
    }
    if controller_manifest and controller_manifest.is_file():
        state["controller_manifest"] = {
            "path": str(controller_manifest.resolve()),
            "sha256": sha256(controller_manifest),
        }
    if error:
        state["failure"] = {
            "exception_type": type(error).__name__,
            "message": str(error),
        }
    state["resume_token"] = resume_token(state)
    current = out / "current_stage.json"
    next_manifest = out / "next_action_manifest.json"
    atomic_json(current, state)
    atomic_json(next_manifest, {
        "schema_version": "1.0",
        "status": status,
        "phase": phase,
        "next_action": next_action,
        "resume_token": state["resume_token"],
        "safe_scheduler_exit": status != "FAILED_RUNTIME",
        "current_stage": {
            "path": str(current.resolve()), "sha256": sha256(current),
        },
    })
    return current, next_manifest
