#!/usr/bin/env python3
"""Materialize the strictly serial, one-cell-type-at-a-time review state."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from evidence_schema_lib import sha256, validate_json_against_schema
from lineage_controller_lib import read_tsv, write_tsv
from membership_transform_lib import deterministic_membership_hash


CLOSING_OUTCOMES = {
    "retain_current_cell_type",
    "confirm_absent_or_not_evaluable",
    "apply_cell_type_membership_patch",
}


def artifact(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise SystemExit(f"cell-type review artifact is missing: {path}")
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "n_bytes": path.stat().st_size,
    }


def bound_artifact(document: dict, *keys: str) -> Path:
    record: object = document
    for key in keys:
        if not isinstance(record, dict):
            raise SystemExit("cell-type review artifact record is malformed")
        record = record.get(key, {})
    if not isinstance(record, dict):
        raise SystemExit("cell-type review artifact record is malformed")
    path = Path(str(record.get("path", "")))
    if not path.is_file() or record.get("sha256") != sha256(path):
        raise SystemExit("cell-type review artifact is missing or stale")
    return path


def bound_record(record: object, label: str) -> Path:
    if not isinstance(record, dict):
        raise SystemExit(f"manual adjudication {label} artifact is malformed")
    path = Path(str(record.get("path", "")))
    if not path.is_file() or record.get("sha256") != sha256(path):
        raise SystemExit(f"manual adjudication {label} artifact is missing or stale")
    return path


def manual_closed_keys(
    adjudication_paths: list[Path], review: dict, queue: list[dict[str, str]],
    maximum_decisions: int,
) -> tuple[set[tuple[str, str, str]], list[dict]]:
    """Validate a legacy/user-authorized non-mutating blocked-scope closure."""
    closed: set[tuple[str, str, str]] = set()
    records: list[dict] = []
    if not adjudication_paths:
        return closed, records
    summary_path = bound_artifact(
        review, "artifacts", "broad_lineage_review_summary"
    )
    review_units = {
        (
            str(row.get("review_mode", "")),
            str(row.get("broad_label", "")),
            str(row.get("unit_signature", "")),
        )
        for row in read_tsv(summary_path)
    }
    current_membership = review.get("membership", {})
    for path in adjudication_paths:
        document = json.loads(path.read_text(encoding="utf-8"))
        if (
            document.get("schema_version") != "1.0"
            or document.get("status") != "PASS"
            or document.get("artifact_role")
            != "user_authorized_manual_biological_adjudication"
            or document.get("outcome") not in {
                "retain_current_cell_type", "confirm_absent_or_not_evaluable"
            }
            or document.get("membership_changed") is not False
            or document.get("counts_as_automatic_decision_round") is not False
        ):
            raise SystemExit("manual adjudication is not a canonical non-mutating closure")
        key = (
            str(document.get("review_mode", "")),
            str(document.get("target_broad_label", "")),
            str(document.get("unit_signature", "")),
        )
        # A canonical re-audit may already have removed a manually closed unit
        # from the open queue.  The broad summary is the authoritative census
        # of exact review units, while the queue contains open units only.
        if not all(key):
            raise SystemExit("manual adjudication lacks an exact scope")
        membership_path = bound_record(document.get("membership"), "membership")
        current_membership_path = Path(str(current_membership.get("path", "")))
        if (
            not current_membership_path.is_file()
            or current_membership.get("sha256") != sha256(current_membership_path)
        ):
            raise SystemExit("current review membership is missing or stale")
        same_physical_membership = (
            current_membership.get("sha256") == sha256(membership_path)
            and current_membership_path.resolve() == membership_path.resolve()
        )
        if not same_physical_membership:
            target_label = str(document.get("target_broad_label", ""))
            current_target = [
                row for row in read_tsv(current_membership_path)
                if str(row.get("final_broad_label", row.get("broad_label", "")))
                == target_label
            ]
            adjudicated_target = [
                row for row in read_tsv(membership_path)
                if str(row.get("final_broad_label", row.get("broad_label", "")))
                == target_label
            ]
            if (
                deterministic_membership_hash(current_target)
                != deterministic_membership_hash(adjudicated_target)
            ):
                raise SystemExit(
                    "manual adjudication target membership differs from current review"
                )
        state_path = bound_record(document.get("blocked_review_state"), "blocked state")
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("status") != "BLOCKED":
            raise SystemExit("manual adjudication source state is not BLOCKED")
        state_membership = state.get("membership", {})
        state_membership_path = bound_record(
            state_membership, "blocked-state membership"
        )
        target_label = str(document.get("target_broad_label", ""))
        state_target = [
            row for row in read_tsv(state_membership_path)
            if str(row.get("final_broad_label", row.get("broad_label", "")))
            == target_label
        ]
        adjudicated_target = [
            row for row in read_tsv(membership_path)
            if str(row.get("final_broad_label", row.get("broad_label", "")))
            == target_label
        ]
        if deterministic_membership_hash(state_target) != deterministic_membership_hash(
            adjudicated_target
        ):
            raise SystemExit(
                "manual adjudication source state has different target membership"
            )
        ledger_path = bound_artifact(state, "task_ledger")
        blocked_rows = [
            row for row in read_tsv(ledger_path)
            if (
                str(row.get("review_mode", "")),
                str(row.get("target_broad_label", "")),
                str(row.get("unit_signature", "")),
            ) == key
            and str(row.get("status", "")) == "blocked_maximum_decisions"
            and int(float(row.get("decision_count", 0) or 0)) >= maximum_decisions
        ]
        if len(blocked_rows) != 1:
            raise SystemExit("manual adjudication lacks one exact exhausted blocked task")
        blocked_review_path = bound_record(
            document.get("blocked_review_manifest"), "blocked review"
        )
        state_review = state.get("review_manifest", {})
        if (
            state_review.get("sha256") != sha256(blocked_review_path)
            or Path(str(state_review.get("path", ""))).resolve()
            != blocked_review_path.resolve()
        ):
            raise SystemExit("manual adjudication is not bound to its blocked review")
        authorization = document.get("user_authorization", {})
        if (
            not isinstance(authorization, dict)
            or authorization.get("explicitly_confirmed") is not True
            or not str(authorization.get("verbatim_text", "")).strip()
        ):
            raise SystemExit("manual adjudication lacks explicit user authorization")
        conclusions = document.get("five_conclusions", {})
        required = {
            "current_member_precision", "whole_query_recall", "molecular_identity",
            "whole_section_spatial_consistency", "literature_boundary_consistency",
        }
        if (
            not isinstance(conclusions, dict)
            or set(conclusions) != required
            or any(not str(conclusions[name]).strip() for name in required)
        ):
            raise SystemExit("manual adjudication lacks the five required conclusions")
        supporting = document.get("supporting_artifacts", [])
        if not isinstance(supporting, list) or len(supporting) < 2:
            raise SystemExit("manual adjudication lacks bound supporting evidence")
        for index, record in enumerate(supporting, 1):
            bound_record(record, f"supporting evidence {index}")
        applies_to_current_scope = key in review_units
        if key in closed:
            raise SystemExit("duplicate manual adjudication for one exact scope")
        if applies_to_current_scope:
            closed.add(key)
        records.append({
            **artifact(path),
            "review_mode": key[0],
            "target_broad_label": key[1],
            "unit_signature": key[2],
            "applies_to_current_scope": applies_to_current_scope,
        })
    return closed, records


def closed_keys(
    validation_paths: list[Path],
) -> tuple[set[str], Counter[str], list[dict]]:
    """Return broad labels that completed their one authoritative review.

    Closure is deliberately label-scoped rather than membership-signature
    scoped.  A later review may transfer observations into or out of an
    already reviewed type, but that bookkeeping change must not schedule a
    second full-query biological review for the closed type.
    """
    closed: set[str] = set()
    decision_counts: Counter[str] = Counter()
    records: list[dict] = []
    for path in validation_paths:
        validation = json.loads(path.read_text(encoding="utf-8"))
        if validation.get("status") != "PASS":
            raise SystemExit("prior cell-type decision validation is not PASS")
        decisions_path = bound_artifact(validation, "validated_decisions")
        review_path = bound_artifact(validation, "review_manifest")
        review = json.loads(review_path.read_text(encoding="utf-8"))
        queue_path = bound_artifact(review, "artifacts", "review_queue")
        queue = {str(row.get("review_id", "")): row for row in read_tsv(queue_path)}
        decisions = read_tsv(decisions_path)
        if len(decisions) != 1:
            raise SystemExit("one validation may close only one cell type")
        decision = decisions[0]
        queued = queue.get(str(decision.get("review_id", "")), {})
        target = str(queued.get("target_broad_label", ""))
        if not target:
            raise SystemExit("prior decision validation lacks a target cell type")
        decision_counts[target] += 1
        if str(decision.get("outcome", "")) in CLOSING_OUTCOMES:
            closed.add(target)
        records.append(artifact(path))
    return closed, decision_counts, records


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--review-manifest", required=True, type=Path)
    ap.add_argument("--previous-state", type=Path)
    ap.add_argument(
        "--prior-decision-validation", action="append", type=Path, default=[]
    )
    ap.add_argument(
        "--manual-biological-adjudication", action="append", type=Path, default=[]
    )
    ap.add_argument("--maximum-decisions-per-cell-type", type=int, default=1)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    review = json.loads(args.review_manifest.read_text(encoding="utf-8"))
    if review.get("stage") != "post_atlas_catalog_wide_lineage_review":
        raise SystemExit("cell-type review state requires the canonical whole-query audit")
    queue_path = bound_artifact(review, "artifacts", "review_queue")
    queue = read_tsv(queue_path)
    if len({str(row.get("target_broad_label", "")) for row in queue}) != len(queue):
        raise SystemExit("cell-type review queue must contain at most one task per type")
    if args.maximum_decisions_per_cell_type < 1:
        raise SystemExit("maximum decisions per cell type must be positive")
    closed, decision_counts, validation_records = closed_keys(
        args.prior_decision_validation
    )
    manual_closed, manual_records = manual_closed_keys(
        args.manual_biological_adjudication, review, queue,
        args.maximum_decisions_per_cell_type,
    )
    closed.update(key[1] for key in manual_closed)

    previous_by_label: dict[str, dict[str, str]] = {}
    previous_record = None
    if args.previous_state:
        previous = json.loads(args.previous_state.read_text(encoding="utf-8"))
        if previous.get("artifact_role") != "sequential_cell_type_review_state":
            raise SystemExit("previous cell-type review state is not canonical")
        previous_tasks_path = bound_artifact(previous, "task_ledger")
        previous_by_label = {
            str(row.get("target_broad_label", "")): row
            for row in read_tsv(previous_tasks_path)
        }
        previous_record = artifact(args.previous_state)

    open_rows: list[dict[str, str]] = []
    tasks: list[dict[str, object]] = []
    for order, row in enumerate(queue, 1):
        key = (
            str(row.get("review_mode", "")),
            str(row.get("target_broad_label", "")),
            str(row.get("unit_signature", "")),
        )
        if key[1] in closed and key not in manual_closed:
            tasks.append({
                "queue_order": order,
                "review_id": row.get("review_id", ""),
                "review_mode": key[0],
                "target_broad_label": key[1],
                "unit_signature": key[2],
                "status": "closed",
                "decision_count": decision_counts[key[1]],
                "required_progress_message": f"现在开始对 {key[1]} 进行专项复核。",
            })
            continue
        if key in manual_closed:
            tasks.append({
                "queue_order": order,
                "review_id": row.get("review_id", ""),
                "review_mode": key[0],
                "target_broad_label": key[1],
                "unit_signature": key[2],
                "status": "closed_by_manual_adjudication",
                "decision_count": decision_counts[key[1]],
                "required_progress_message": f"现在开始对 {key[1]} 进行专项复核。",
            })
            continue
        exhausted = (
            decision_counts[key[1]] >= args.maximum_decisions_per_cell_type
        )
        if not exhausted:
            open_rows.append(row)
        tasks.append({
            "queue_order": order,
            "review_id": row.get("review_id", ""),
            "review_mode": key[0],
            "target_broad_label": key[1],
            "unit_signature": key[2],
            "status": (
                "blocked_maximum_decisions"
                if exhausted else "queued"
            ),
            "decision_count": decision_counts[key[1]],
            "required_progress_message": f"现在开始对 {key[1]} 进行专项复核。",
        })

    # Manually closed units are absent from the new open queue by design.  Add
    # them explicitly so the final ledger remains a complete biological audit.
    task_keys = {
        (
            str(row.get("review_mode", "")),
            str(row.get("target_broad_label", "")),
            str(row.get("unit_signature", "")),
        )
        for row in tasks
    }
    for record in manual_records:
        if record.get("applies_to_current_scope") is not True:
            continue
        key = (
            str(record.get("review_mode", "")),
            str(record.get("target_broad_label", "")),
            str(record.get("unit_signature", "")),
        )
        if key in task_keys:
            continue
        tasks.append({
            "queue_order": 0,
            "review_id": "",
            "review_mode": key[0],
            "target_broad_label": key[1],
            "unit_signature": key[2],
            "status": "closed_by_manual_adjudication",
            "decision_count": decision_counts[key[1]],
            "required_progress_message": f"现在开始对 {key[1]} 进行专项复核。",
        })

    blocked = any(
        row["status"] == "blocked_maximum_decisions" for row in tasks
    )
    active = None if blocked else next(
        (row for row in tasks if row["status"] == "queued"),
        None,
    )
    if active is not None:
        active["status"] = "active"
    if sum(row["status"] == "active" for row in tasks) > 1:
        raise SystemExit("more than one cell type became active")

    # Preserve exact closed targets in the ledger for an auditable census.
    open_labels = {str(row["target_broad_label"]) for row in tasks}
    for label, previous in sorted(previous_by_label.items()):
        if label in open_labels or str(previous.get("status", "")) != "closed":
            continue
        tasks.append({
            "queue_order": int(float(previous.get("queue_order", 0) or 0)),
            "review_id": previous.get("review_id", ""),
            "review_mode": previous.get("review_mode", ""),
            "target_broad_label": label,
            "unit_signature": previous.get("unit_signature", ""),
            "status": "closed",
            "decision_count": decision_counts[label],
            "required_progress_message": previous.get(
                "required_progress_message", f"现在开始对 {label} 进行专项复核。"
            ),
        })

    args.out.mkdir(parents=True, exist_ok=True)
    ledger_path = args.out / "cell_type_review_task_ledger.tsv"
    write_tsv(ledger_path, tasks, fields=[
        "queue_order", "review_id", "review_mode", "target_broad_label",
        "unit_signature", "status", "decision_count",
        "required_progress_message",
    ])
    membership = review.get("membership", {})
    state = {
        "schema_version": "1.0",
        "status": (
            "BLOCKED"
            if blocked
            else "COMPLETE" if active is None else "REVIEW_REQUIRED"
        ),
        "artifact_role": "sequential_cell_type_review_state",
        "user_facing_stage_name": "逐大类全样本复核",
        "review_manifest": artifact(args.review_manifest),
        "review_queue": artifact(queue_path),
        "membership": membership,
        "previous_state": previous_record,
        "prior_decision_validations": validation_records,
        "manual_biological_adjudications": manual_records,
        "task_ledger": artifact(ledger_path),
        "active_cell_type_review": (
            {
                "review_id": active["review_id"],
                "target_broad_label": active["target_broad_label"],
                "unit_signature": active["unit_signature"],
                "required_progress_message": active["required_progress_message"],
            }
            if active is not None else None
        ),
        "active_review_n": 0 if active is None else 1,
        "queued_review_n": sum(row["status"] == "queued" for row in tasks),
        "closed_review_n": sum(
            row["status"] in {"closed", "closed_by_manual_adjudication"}
            for row in tasks
        ),
        "blocked_review_n": sum(
            row["status"] == "blocked_maximum_decisions" for row in tasks
        ),
        "maximum_decisions_per_cell_type": args.maximum_decisions_per_cell_type,
        "single_pass_no_reopen": True,
        "membership_changes_do_not_reopen_closed_types": True,
        "closed_types_accept_bounded_incoming_writeback": True,
        "closed_cell_type_labels": sorted(closed),
        "formal_batch_closure_forbidden": True,
        "next_action": (
            "manual_biological_adjudication_required"
            if blocked
            else "all_cell_types_closed" if active is None
            else active["required_progress_message"]
        ),
    }
    state_path = args.out / "cell_type_review_state.json"
    state_path.write_text(
        json.dumps(state, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    schema = (
        Path(__file__).resolve().parents[1]
        / "schemas/cell_type_review_state.schema.json"
    )
    _, errors = validate_json_against_schema(state_path, schema)
    if errors:
        raise SystemExit(
            "cell-type review state violates schema: " + "; ".join(errors)
        )
    print(json.dumps(state, ensure_ascii=False, indent=2))
    return 0 if state["status"] == "COMPLETE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
