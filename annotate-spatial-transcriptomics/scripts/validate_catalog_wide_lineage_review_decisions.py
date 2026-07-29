#!/usr/bin/env python3
"""Validate exact, evidence-backed decisions for one catalog-wide review round."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from evidence_schema_lib import sha256
from lineage_controller_lib import read_tsv, write_tsv


ALLOWED = {
    "broad_lineage_review": {
        "retain_current_cell_type",
        "apply_cell_type_membership_patch",
        "targeted_raw_count_review_required",
    },
    "missing_broad_review": {
        "confirm_absent_or_not_evaluable",
        "apply_cell_type_membership_patch",
        "targeted_raw_count_review_required",
    },
    "present_label_precision": {
        "retain_current_label",
        "apply_targeted_membership_patch",
        "targeted_raw_count_review_required",
    },
    "outside_label_recall": {
        "accept_recall_component",
        "reject_shared_or_ambient",
        "retain_current_parent",
        "apply_targeted_membership_patch",
        "targeted_raw_count_review_required",
    },
    "outside_label_group_watch": {
        "reject_shared_or_ambient",
        "retain_current_parent",
        "apply_targeted_membership_patch",
        "targeted_raw_count_review_required",
    },
}

CELL_TYPE_PRECISION = {"supported", "over_recall_detected", "not_applicable", "not_evaluable"}
CELL_TYPE_RECALL = {"complete", "under_recall_detected", "confirmed_absent", "not_evaluable"}
CELL_TYPE_SPATIAL = {"consistent", "localized_issue", "inconsistent", "not_evaluable"}
CELL_TYPE_MOLECULAR = {"supported", "mixed", "unsupported", "not_evaluable"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--review-manifest", required=True, type=Path)
    ap.add_argument("--evidence-packet-manifest", required=True, type=Path)
    ap.add_argument("--decisions", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    review = json.loads(args.review_manifest.read_text(encoding="utf-8"))
    if (
        review.get("stage") != "post_atlas_catalog_wide_lineage_review"
        or review.get("status") != "ITERATION_REQUIRED"
        or review.get("catalog_wide_double_sided_review") is not True
    ):
        raise SystemExit("catalog-wide review is not an open canonical round")
    queue_record = review.get("artifacts", {}).get("review_queue", {})
    queue_path = Path(str(queue_record.get("path", "")))
    if not queue_path.is_file() or queue_record.get("sha256") != sha256(queue_path):
        raise SystemExit("catalog-wide review queue is missing or stale")
    queue_rows = read_tsv(queue_path)
    queue = {str(row.get("review_id", "")): row for row in queue_rows}
    packet_manifest = json.loads(
        args.evidence_packet_manifest.read_text(encoding="utf-8")
    )
    packet_review = packet_manifest.get("review_manifest", {})
    if (
        packet_manifest.get("status") != "PASS"
        or packet_manifest.get("artifact_role")
        != "broad_cell_type_review_evidence_packet_index"
        or Path(str(packet_review.get("path", ""))).resolve()
        != args.review_manifest.resolve()
        or packet_review.get("sha256") != sha256(args.review_manifest)
    ):
        raise SystemExit("broad review evidence packets do not bind this open review")
    packet_record = packet_manifest.get("packet_index", {})
    packet_path = Path(str(packet_record.get("path", "")))
    if (
        not packet_path.is_file()
        or packet_record.get("sha256") != sha256(packet_path)
    ):
        raise SystemExit("broad review evidence packet index is missing or stale")
    packet_rows = read_tsv(packet_path)
    packets = {str(row.get("review_id", "")): row for row in packet_rows}
    if (
        "" in packets or len(packets) != len(packet_rows)
        or set(packets) != set(queue)
    ):
        raise SystemExit("broad review evidence packets do not cover the queue exactly")
    decisions = read_tsv(args.decisions)
    decision_by_id = {str(row.get("review_id", "")): row for row in decisions}
    errors: list[str] = []
    if (
        "" in queue or len(queue) != len(queue_rows)
        or "" in decision_by_id or len(decision_by_id) != len(decisions)
        or set(queue) != set(decision_by_id)
    ):
        errors.append("decisions must cover the review queue exactly once")
    pending = 0
    for review_id, queued in queue.items():
        decision = decision_by_id.get(review_id, {})
        mode = str(queued.get("review_mode", ""))
        outcome = str(decision.get("outcome", ""))
        if outcome not in ALLOWED.get(mode, set()):
            errors.append(f"{review_id}: outcome is invalid for {mode}")
            continue
        if str(decision.get("review_mode", "")) != mode:
            errors.append(f"{review_id}: decision review_mode differs from queue")
        packet = packets.get(review_id, {})
        if (
            str(packet.get("review_mode", "")) != mode
            or str(packet.get("target_broad_label", ""))
            != str(queued.get("target_broad_label", ""))
            or str(packet.get("unit_signature", ""))
            != str(queued.get("unit_signature", ""))
            or str(decision.get("evidence_packet_sha256", ""))
            != str(packet.get("evidence_packet_sha256", ""))
        ):
            errors.append(f"{review_id}: decision does not bind its canonical evidence packet")
        rationale = str(decision.get("rationale", "")).strip()
        if len(rationale) < 20:
            errors.append(f"{review_id}: biological rationale is too short")
        if mode in {"broad_lineage_review", "missing_broad_review"}:
            precision = str(decision.get("current_member_precision", "")).strip()
            recall = str(decision.get("whole_query_recall", "")).strip()
            spatial = str(decision.get("spatial_consistency", "")).strip()
            molecular = str(decision.get("molecular_support", "")).strip()
            if precision not in CELL_TYPE_PRECISION:
                errors.append(f"{review_id}: current_member_precision is missing or invalid")
            if recall not in CELL_TYPE_RECALL:
                errors.append(f"{review_id}: whole_query_recall is missing or invalid")
            if spatial not in CELL_TYPE_SPATIAL:
                errors.append(f"{review_id}: spatial_consistency is missing or invalid")
            if molecular not in CELL_TYPE_MOLECULAR:
                errors.append(f"{review_id}: molecular_support is missing or invalid")
            precision_evaluable = str(packet.get("precision_evaluable", "")).lower() == "true"
            recall_evaluable = str(packet.get("recall_evaluable", "")).lower() == "true"
            molecular_evaluable = str(packet.get("molecular_evaluable", "")).lower() == "true"
            spatial_evaluable = str(packet.get("spatial_evaluable", "")).lower() == "true"
            if precision in {"supported", "over_recall_detected"} and not precision_evaluable:
                errors.append(f"{review_id}: precision conclusion lacks evaluable bound evidence")
            if recall in {"complete", "under_recall_detected", "confirmed_absent"} and not recall_evaluable:
                errors.append(f"{review_id}: recall conclusion lacks evaluable bound evidence")
            if molecular != "not_evaluable" and not molecular_evaluable:
                errors.append(f"{review_id}: molecular conclusion lacks evaluable bound evidence")
            if spatial != "not_evaluable" and not spatial_evaluable:
                errors.append(f"{review_id}: spatial conclusion lacks evaluable bound evidence")
            if mode == "broad_lineage_review":
                if outcome != "targeted_raw_count_review_required" and precision in {"not_applicable", "not_evaluable"}:
                    errors.append(f"{review_id}: an annotated broad requires current-member precision review")
                if outcome != "targeted_raw_count_review_required" and recall in {"confirmed_absent", "not_evaluable"}:
                    errors.append(f"{review_id}: an annotated broad requires whole-query recall review")
                if outcome == "retain_current_cell_type" and (
                    precision != "supported"
                    or recall != "complete"
                    or spatial != "consistent"
                    or molecular != "supported"
                ):
                    errors.append(
                        f"{review_id}: retaining a broad requires supported precision, complete recall, "
                        "consistent spatial distribution and supported molecular identity"
                    )
                if str(packet.get("current_precision_question_n", "")).strip():
                    current_question_n = int(float(
                        packet.get("current_precision_question_n", 0) or 0
                    ))
                else:
                    current_question_n = int(float(
                        packet.get("current_competitor_question_n", 0) or 0
                    )) + int(float(
                        packet.get("cross_type_over_recall_question_n", 0) or 0
                    ))
                recall_question_n = int(float(
                    packet.get("outside_recall_question_n", 0) or 0
                ))
                if (
                    outcome == "retain_current_cell_type"
                    and current_question_n > 0
                ):
                    errors.append(
                        f"{review_id}: current-member challengers require an exact patch or targeted review"
                    )
                if (
                    outcome == "retain_current_cell_type"
                    and recall_question_n > 0
                ):
                    errors.append(
                        f"{review_id}: whole-query recall challengers require an exact patch or targeted review"
                    )
                ovary_spatial_status = str(
                    packet.get("ovary_spatial_status", "not_applicable")
                )
                oocyte_review_status = str(
                    packet.get("oocyte_review_status", "not_applicable")
                )
                follicle_status = str(
                    packet.get("follicle_histology_status", "not_applicable")
                )
                if (
                    outcome == "retain_current_cell_type"
                    and ovary_spatial_status not in {"PASS", "not_applicable"}
                ):
                    errors.append(
                        f"{review_id}: sheep-ovary broad spatial review is not closed"
                    )
                if (
                    str(queued.get("target_broad_label", "")) == "Oocyte"
                    and outcome == "retain_current_cell_type"
                    and (
                        oocyte_review_status != "PASS"
                        or str(packet.get("oocyte_canonical_bound", "")).lower()
                        != "true"
                    )
                ):
                    errors.append(
                        f"{review_id}: Oocyte cannot close without the canonical quality review"
                    )
                follicle_related = {
                    "Granulosa", "Theca", "Vascular-associated",
                    "Smooth muscle", "Stromal/mesenchymal",
                }
                if (
                    str(queued.get("target_broad_label", ""))
                    in follicle_related
                    and outcome == "retain_current_cell_type"
                    and follicle_status not in {"PASS", "not_applicable"}
                ):
                    errors.append(
                        f"{review_id}: follicle ROI histology remains open for this lineage"
                    )
                if outcome == "apply_cell_type_membership_patch" and not (
                    precision == "over_recall_detected"
                    or recall == "under_recall_detected"
                    or spatial in {"localized_issue", "inconsistent"}
                    or molecular in {"mixed", "unsupported"}
                ):
                    errors.append(f"{review_id}: a cell-type patch lacks a stated broad-level defect")
            else:
                if precision != "not_applicable":
                    errors.append(f"{review_id}: missing-broad review precision must be not_applicable")
                if outcome == "confirm_absent_or_not_evaluable" and recall not in {
                    "confirmed_absent", "not_evaluable"
                }:
                    errors.append(f"{review_id}: absent broad was not actually closed by recall review")
        proposed = str(decision.get("proposed_broad_label", ""))
        if outcome == "accept_recall_component" and proposed != str(
            queued.get("target_broad_label", "")
        ):
            errors.append(f"{review_id}: accepted recall changed its target broad")
        if mode == "outside_label_group_watch" and outcome == "accept_recall_component":
            errors.append(f"{review_id}: a group watch cannot assign a whole source group")
        patch_path = str(decision.get("membership_path", "")).strip()
        patch_hash = str(decision.get("membership_sha256", "")).strip()
        if outcome in {
            "apply_targeted_membership_patch", "apply_cell_type_membership_patch",
        }:
            path = Path(patch_path)
            if not path.is_file() or patch_hash != sha256(path):
                errors.append(f"{review_id}: targeted membership patch is missing or stale")
        elif patch_path or patch_hash:
            errors.append(f"{review_id}: non-patch decision carries membership authority")
        if outcome == "targeted_raw_count_review_required":
            pending += 1

    args.out.mkdir(parents=True, exist_ok=True)
    validated_path = args.out / "validated_catalog_wide_lineage_decisions.tsv"
    write_tsv(validated_path, decisions)
    status = "PASS" if not errors and not pending else "ITERATION_REQUIRED" if not errors else "BLOCKED"
    manifest = {
        "schema_version": "2.2",
        "status": status,
        "stage": "catalog_wide_lineage_review_decision_validation",
        "review_round": review.get("review_round"),
        "review_manifest": {
            "path": str(args.review_manifest.resolve()),
            "sha256": sha256(args.review_manifest),
        },
        "review_queue": {
            "path": str(queue_path.resolve()), "sha256": sha256(queue_path),
        },
        "evidence_packet_manifest": {
            "path": str(args.evidence_packet_manifest.resolve()),
            "sha256": sha256(args.evidence_packet_manifest),
        },
        "evidence_packet_index": {
            "path": str(packet_path.resolve()), "sha256": sha256(packet_path),
        },
        "source_decisions": {
            "path": str(args.decisions.resolve()), "sha256": sha256(args.decisions),
        },
        "validated_decisions": {
            "path": str(validated_path.resolve()), "sha256": sha256(validated_path),
        },
        "decision_n": len(decisions),
        "pending_targeted_review_n": pending,
        "errors": errors,
    }
    manifest_path = args.out / "catalog_wide_lineage_decision_validation.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
