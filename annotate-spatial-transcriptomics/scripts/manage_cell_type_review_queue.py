#!/usr/bin/env python3
"""Materialize the strictly serial, one-cell-type-at-a-time review state."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from evidence_schema_lib import sha256, validate_json_against_schema
from lineage_controller_lib import read_tsv, write_tsv


CLOSING_OUTCOMES = {
    "retain_current_cell_type",
    "confirm_absent_or_not_evaluable",
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


def closed_keys(
    validation_paths: list[Path],
) -> tuple[set[tuple[str, str, str]], Counter[str], list[dict]]:
    closed: set[tuple[str, str, str]] = set()
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
            key = (
                str(queued.get("review_mode", "")),
                str(queued.get("target_broad_label", "")),
                str(queued.get("unit_signature", "")),
            )
            if not all(key):
                raise SystemExit("closed cell-type decision lacks a stable scope")
            closed.add(key)
        records.append(artifact(path))
    return closed, decision_counts, records


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--review-manifest", required=True, type=Path)
    ap.add_argument("--previous-state", type=Path)
    ap.add_argument(
        "--prior-decision-validation", action="append", type=Path, default=[]
    )
    ap.add_argument("--maximum-decisions-per-cell-type", type=int, default=2)
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
        if key in closed:
            continue
        exhausted = (
            decision_counts[key[1]] >= args.maximum_decisions_per_cell_type
        )
        if not exhausted:
            open_rows.append(row)
        previous_task = previous_by_label.get(key[1], {})
        reopened = bool(
            previous_task
            and str(previous_task.get("unit_signature", "")) != key[2]
            and str(previous_task.get("status", "")) == "closed"
        )
        tasks.append({
            "queue_order": order,
            "review_id": row.get("review_id", ""),
            "review_mode": key[0],
            "target_broad_label": key[1],
            "unit_signature": key[2],
            "status": (
                "blocked_maximum_decisions"
                if exhausted else "reopened" if reopened else "queued"
            ),
            "decision_count": decision_counts[key[1]],
            "required_progress_message": f"现在开始对 {key[1]} 进行专项复核。",
        })

    blocked = any(
        row["status"] == "blocked_maximum_decisions" for row in tasks
    )
    active = None if blocked else next(
        (row for row in tasks if row["status"] in {"queued", "reopened"}),
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
        "queued_review_n": sum(row["status"] in {"queued", "reopened"} for row in tasks),
        "closed_review_n": sum(row["status"] == "closed" for row in tasks),
        "blocked_review_n": sum(
            row["status"] == "blocked_maximum_decisions" for row in tasks
        ),
        "maximum_decisions_per_cell_type": args.maximum_decisions_per_cell_type,
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
