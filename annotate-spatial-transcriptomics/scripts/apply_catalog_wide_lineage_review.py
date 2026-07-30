#!/usr/bin/env python3
"""Apply one validated decision for the single active cell-type review."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from evidence_schema_lib import sha256
from lineage_controller_lib import (
    apply_candidate_context, candidate_can_release, catalog_candidates,
    deterministic_cell_id_set_hash, deterministic_membership_hash,
    read_tsv, write_tsv,
)


def artifact(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()), "sha256": sha256(path),
        "n_bytes": path.stat().st_size,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", required=True, type=Path)
    ap.add_argument("--stage-authority", required=True, type=Path)
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--review-manifest", required=True, type=Path)
    ap.add_argument("--review-state", required=True, type=Path)
    ap.add_argument("--decision-validation", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    authority = json.loads(args.stage_authority.read_text(encoding="utf-8"))
    if (
        authority.get("mode") != "stage_authority"
        or authority.get("phase") != "atlas_and_completeness_review"
        or authority.get("annotation_contract_sha256") != sha256(args.contract)
    ):
        raise SystemExit("stage authority does not permit catalog-wide review apply")
    review = json.loads(args.review_manifest.read_text(encoding="utf-8"))
    state = json.loads(args.review_state.read_text(encoding="utf-8"))
    validation = json.loads(args.decision_validation.read_text(encoding="utf-8"))
    active = state.get("active_cell_type_review") or {}
    if (
        validation.get("status") != "PASS"
        or validation.get("formal_batch_closure_performed") is not False
        or int(validation.get("decision_n", -1)) != 1
        or state.get("artifact_role") != "sequential_cell_type_review_state"
        or int(state.get("active_review_n", -1)) != 1
        or state.get("formal_batch_closure_forbidden") is not True
        or validation.get("active_review_id") != active.get("review_id")
        or validation.get("active_cell_type") != active.get("target_broad_label")
        or Path(str(validation.get("review_state", {}).get("path", ""))).resolve()
        != args.review_state.resolve()
        or validation.get("review_state", {}).get("sha256") != sha256(args.review_state)
    ):
        raise SystemExit("catalog-wide review decisions are not fully resolved")
    packet_record = validation.get("evidence_packet_manifest", {})
    packet_path = Path(str(packet_record.get("path", "")))
    if (
        not packet_path.is_file()
        or packet_record.get("sha256") != sha256(packet_path)
    ):
        raise SystemExit("catalog-wide decision evidence packets are missing or stale")
    packet_manifest = json.loads(packet_path.read_text(encoding="utf-8"))
    broad_record = packet_manifest.get("broad_evidence_manifest", {})
    broad_manifest_path = Path(str(broad_record.get("path", "")))
    if (
        not broad_manifest_path.is_file()
        or broad_record.get("sha256") != sha256(broad_manifest_path)
    ):
        raise SystemExit("catalog-wide decision lacks bound broad evidence")
    broad_manifest = json.loads(
        broad_manifest_path.read_text(encoding="utf-8")
    )
    question_record = broad_manifest.get("artifacts", {}).get(
        "current_member_questions", {}
    )
    recall_record = broad_manifest.get("artifacts", {}).get(
        "recall_membership", {}
    )
    question_path = Path(str(question_record.get("path", "")))
    recall_path = Path(str(recall_record.get("path", "")))
    if (
        not question_path.is_file()
        or question_record.get("sha256") != sha256(question_path)
        or not recall_path.is_file()
        or recall_record.get("sha256") != sha256(recall_path)
    ):
        raise SystemExit("catalog-wide exact challenger memberships are missing or stale")
    precision_question_ids: dict[str, set[str]] = {}
    for row in read_tsv(question_path):
        precision_question_ids.setdefault(
            str(row.get("broad_label", "")), set()
        ).add(str(row.get("cell_id", "")))
    recall_question_ids: dict[str, set[str]] = {}
    for row in read_tsv(recall_path):
        recall_question_ids.setdefault(
            str(row.get("broad_label", "")), set()
        ).add(str(row.get("cell_id", "")))
    authority_records = {
        "post_atlas_membership": args.membership,
        "catalog_review_manifest": args.review_manifest,
        "cell_type_review_state": args.review_state,
        "catalog_decision_validation": args.decision_validation,
        "candidate_catalog": args.catalog,
    }
    if args.context_evidence:
        authority_records["context_evidence"] = args.context_evidence
    for key, path in authority_records.items():
        record = authority.get(key, {})
        if (
            Path(str(record.get("path", ""))).resolve() != path.resolve()
            or record.get("sha256") != sha256(path)
        ):
            raise SystemExit(f"catalog-wide apply authority differs for {key}")
    if (
        Path(str(review.get("membership", {}).get("path", ""))).resolve()
        != args.membership.resolve()
        or review.get("membership", {}).get("sha256") != sha256(args.membership)
        or Path(str(validation.get("review_manifest", {}).get("path", ""))).resolve()
        != args.review_manifest.resolve()
        or validation.get("review_manifest", {}).get("sha256")
        != sha256(args.review_manifest)
    ):
        raise SystemExit("catalog-wide decisions do not bind the supplied membership/review")
    decision_path = Path(str(validation.get("validated_decisions", {}).get("path", "")))
    if not decision_path.is_file() or validation.get("validated_decisions", {}).get("sha256") != sha256(decision_path):
        raise SystemExit("validated catalog-wide decisions are missing or stale")
    queue_path = Path(str(review.get("artifacts", {}).get("review_queue", {}).get("path", "")))
    component_membership_path = Path(str(
        review.get("artifacts", {}).get("outside_label_recall_component_membership", {}).get("path", "")
    ))
    type_scope_path = Path(str(
        review.get("artifacts", {}).get("broad_lineage_review_scope_membership", {}).get("path", "")
    ))
    for path, record_name in (
        (queue_path, "review_queue"),
        (component_membership_path, "outside_label_recall_component_membership"),
        (type_scope_path, "broad_lineage_review_scope_membership"),
    ):
        record = review.get("artifacts", {}).get(record_name, {})
        if not path.is_file() or record.get("sha256") != sha256(path):
            raise SystemExit(f"catalog-wide review artifact is missing or stale: {record_name}")

    candidates = catalog_candidates(json.loads(args.catalog.read_text(encoding="utf-8")))
    context_summary = apply_candidate_context(
        candidates, read_tsv(args.context_evidence) if args.context_evidence else []
    )
    eligible_broad = {
        str(candidate.get("release_broad_label", ""))
        for candidate in candidates.values()
        if str(candidate.get("candidate_role", "")) == "broad"
        and candidate_can_release(candidate)
    }
    queue = {str(row.get("review_id", "")): row for row in read_tsv(queue_path)}
    decisions = read_tsv(decision_path)
    if (
        len(decisions) != 1
        or str(decisions[0].get("review_id", "")) != str(active.get("review_id", ""))
    ):
        raise SystemExit("apply may modify only the single active cell type")
    component_members: dict[str, list[dict[str, str]]] = {}
    for row in read_tsv(component_membership_path):
        component_members.setdefault(str(row.get("component_id", "")), []).append(row)
    type_scope_members: dict[str, set[str]] = {}
    for row in read_tsv(type_scope_path):
        type_scope_members.setdefault(str(row.get("review_id", "")), set()).add(
            str(row.get("cell_id", ""))
        )
    source = read_tsv(args.membership)
    by_cell = {str(row.get("cell_id", "")): dict(row) for row in source}
    if not by_cell or "" in by_cell or len(by_cell) != len(source):
        raise SystemExit("catalog-wide review membership is empty or duplicated")

    proposed: dict[str, tuple[str, str, str, str]] = {}
    decision_audit: list[dict[str, object]] = []
    for decision in decisions:
        review_id = str(decision.get("review_id", ""))
        queued = queue[review_id]
        outcome = str(decision.get("outcome", ""))
        affected: list[tuple[str, str, str]] = []
        if outcome == "accept_recall_component":
            target = str(queued.get("target_broad_label", ""))
            if target not in eligible_broad:
                raise SystemExit(f"{review_id}: target broad is context-ineligible")
            for row in component_members.get(str(queued.get("component_id", "")), []):
                affected.append((
                    str(row.get("cell_id", "")), target,
                    str(row.get("candidate_id", "")),
                ))
        elif outcome in {
            "apply_targeted_membership_patch", "apply_cell_type_membership_patch",
        }:
            patch_path = Path(str(decision.get("membership_path", "")))
            for row in read_tsv(patch_path):
                cell = str(row.get("cell_id", ""))
                target = str(row.get("new_broad_label", ""))
                candidate_id = str(row.get("candidate_id", ""))
                if target and target not in eligible_broad:
                    raise SystemExit(f"{review_id}: patch target is context-ineligible")
                candidate = candidates.get(candidate_id, {})
                if target and (
                    not candidate_can_release(candidate)
                    or str(candidate.get("release_broad_label", "")) != target
                ):
                    raise SystemExit(
                        f"{review_id}: patch candidate cannot release the target broad"
                    )
                affected.append((cell, target, candidate_id))
        for cell, target, candidate_id in affected:
            if cell not in by_cell:
                raise SystemExit(f"{review_id}: patch contains a foreign cell")
            review_mode = str(queued.get("review_mode", ""))
            review_target = str(queued.get("target_broad_label", ""))
            current_label = str(by_cell[cell].get("final_broad_label", ""))
            if review_mode in {"broad_lineage_review", "missing_broad_review"}:
                # One cell-type review can remove doubtful current members to
                # an evidenced competitor, or recall outside members only to
                # the reviewed target.  It cannot use a whole-query scope to
                # rewrite unrelated lineages.
                if current_label != review_target and target != review_target:
                    raise SystemExit(
                        f"{review_id}: cell-type patch assigns an unrelated outside member"
                    )
                if current_label == review_target and target != current_label:
                    if cell not in precision_question_ids.get(review_target, set()):
                        raise SystemExit(
                            f"{review_id}: current-member patch is outside the bound precision questions"
                        )
                elif current_label != review_target and target == review_target:
                    if cell not in recall_question_ids.get(review_target, set()):
                        raise SystemExit(
                            f"{review_id}: recall patch is outside the bound whole-query questions"
                        )
            if review_mode == "outside_label_recall":
                allowed = {
                    str(row.get("cell_id", ""))
                    for row in component_members.get(str(queued.get("component_id", "")), [])
                }
            elif review_mode in {"broad_lineage_review", "missing_broad_review"}:
                allowed = type_scope_members.get(review_id, set())
            elif review_mode == "outside_label_group_watch":
                allowed = {
                    key for key, row in by_cell.items()
                    if str(row.get("source_boundary", "")) == str(queued.get("source_boundary", ""))
                    and str(row.get("source_cluster", "")) == str(queued.get("source_cluster", ""))
                }
            else:
                allowed = {
                    key for key, row in by_cell.items()
                    if str(row.get("source_boundary", "")) == str(queued.get("source_boundary", ""))
                    and str(row.get("source_cluster", "")) == str(queued.get("source_cluster", ""))
                    and str(row.get("final_broad_label", "")) == str(queued.get("target_broad_label", ""))
                }
            if cell not in allowed:
                raise SystemExit(f"{review_id}: patch escaped the bounded review unit")
            previous = proposed.get(cell)
            value = (target, candidate_id, review_id, outcome)
            if previous and previous[:2] != value[:2]:
                raise SystemExit("catalog-wide decisions assign competing labels to one cell")
            proposed[cell] = value
        decision_audit.append({
            "review_id": review_id,
            "review_mode": queued.get("review_mode", ""),
            "outcome": outcome,
            "target_broad_label": queued.get("target_broad_label", ""),
            "n_affected": len(affected),
            "current_member_precision": decision.get("current_member_precision", ""),
            "whole_query_recall": decision.get("whole_query_recall", ""),
            "spatial_consistency": decision.get("spatial_consistency", ""),
            "molecular_support": decision.get("molecular_support", ""),
            "evidence_basis": decision.get("evidence_basis", ""),
            "rationale": decision.get("rationale", ""),
        })

    changes: list[dict[str, object]] = []
    output: list[dict[str, object]] = []
    round_index = int(review.get("review_round", 0))
    for cell in sorted(by_cell):
        row = dict(by_cell[cell])
        proposal = proposed.get(cell)
        if proposal:
            new_broad, candidate_id, review_id, outcome = proposal
            old_broad = str(row.get("final_broad_label", ""))
            if new_broad != old_broad:
                row["final_broad_label"] = new_broad
                row["final_state"] = (
                    "defined_broad_only" if new_broad else "unresolved_biological"
                )
                row["candidate_id"] = candidate_id if new_broad else ""
                row["confidence"] = "moderate" if new_broad else "low"
                row["assignment_origin"] = (
                    f"catalog_wide_lineage_review_round_{round_index}__{outcome}"
                )
                row["unresolved_reason"] = (
                    "" if new_broad else "catalog_wide_review_unresolved"
                )
                for column in (
                    "final_fine_label", "final_fine_confidence",
                    "final_fine_candidate_id", "final_fine_assignment_source",
                ):
                    if column in row:
                        row[column] = ""
                if "fine_anchor_eligible" in row:
                    row["fine_anchor_eligible"] = "false"
                changes.append({
                    "cell_id": cell, "old_broad_label": old_broad,
                    "new_broad_label": new_broad, "candidate_id": candidate_id,
                    "review_id": review_id, "outcome": outcome,
                })
        output.append(row)

    args.out.mkdir(parents=True, exist_ok=True)
    membership_path = args.out / "catalog_wide_reviewed_membership.tsv.gz"
    changes_path = args.out / "catalog_wide_review_changes.tsv"
    audit_path = args.out / "catalog_wide_review_decision_audit.tsv"
    write_tsv(membership_path, output)
    write_tsv(changes_path, changes, fields=[
        "cell_id", "old_broad_label", "new_broad_label", "candidate_id",
        "review_id", "outcome",
    ])
    write_tsv(audit_path, decision_audit)
    transitions = Counter(
        (str(row["old_broad_label"]), str(row["new_broad_label"]))
        for row in changes
    )
    manifest = {
        "schema_version": "2.2",
        "status": "PASS_REQUIRES_NEXT_REVIEW_ROUND",
        "stage": "catalog_wide_lineage_review_apply",
        "active_cell_type": active.get("target_broad_label", ""),
        "active_review_id": active.get("review_id", ""),
        "formal_batch_closure_performed": False,
        "review_round": round_index,
        "formal_membership_written": False,
        "membership": artifact(membership_path),
        "source_membership": artifact(args.membership),
        "annotation_contract": artifact(args.contract),
        "stage_authority": artifact(args.stage_authority),
        "review_manifest": artifact(args.review_manifest),
        "review_state": artifact(args.review_state),
        "decision_validation": artifact(args.decision_validation),
        "candidate_catalog": artifact(args.catalog),
        "context_evidence": artifact(args.context_evidence) if args.context_evidence else None,
        "context_release_eligibility": context_summary,
        "changes": artifact(changes_path),
        "decision_audit": artifact(audit_path),
        "n_changed_observations": len(changes),
        "changed_transitions": {
            f"{old} -> {new}": count
            for (old, new), count in sorted(transitions.items())
        },
        "membership_semantic_sha256": deterministic_membership_hash(output),
        "membership_transform": {
            "operation": "cell_type_review_patch",
            "target_cell_type": active.get("target_broad_label", ""),
            "source_physical_sha256": sha256(args.membership),
            "source_semantic_sha256": deterministic_membership_hash(source),
            "result_physical_sha256": sha256(membership_path),
            "result_semantic_sha256": deterministic_membership_hash(output),
            "source_cell_id_set_sha256": deterministic_cell_id_set_hash(source),
            "result_cell_id_set_sha256": deterministic_cell_id_set_hash(output),
            "delta_physical_sha256": sha256(changes_path),
            "changed_observation_n": len(changes),
        },
        "next_required_action": (
            "rerun_same_cell_type_if_membership_changed_otherwise_activate_next_cell_type"
        ),
    }
    manifest_path = args.out / "catalog_wide_lineage_review_apply_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
