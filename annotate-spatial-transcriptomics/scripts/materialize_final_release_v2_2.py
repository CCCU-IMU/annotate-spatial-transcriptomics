#!/usr/bin/env python3
"""Materialize final v2.2 broad/fine/state/QC membership after all reviews pass."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from controller_thresholds import load_controller_thresholds
from evidence_schema_lib import sha256
from lineage_controller_lib import (
    apply_candidate_context, candidate_can_release, catalog_candidates,
    deterministic_membership_hash, read_tsv, write_tsv,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", required=True, type=Path)
    ap.add_argument("--stage-authority", required=True, type=Path)
    ap.add_argument("--post-atlas-membership", required=True, type=Path)
    ap.add_argument("--atlas-completeness-manifest", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument("--fine-assignments", type=Path)
    ap.add_argument("--state-proposals", action="append", type=Path, default=[])
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    authority = json.loads(args.stage_authority.read_text(encoding="utf-8"))
    if (
        authority.get("mode") != "stage_authority"
        or authority.get("phase") != "materialize_final_release"
        or authority.get("annotation_contract_sha256") != sha256(args.contract)
    ):
        raise SystemExit("stage authority does not permit final materialization")
    threshold_record = contract.get("threshold_registry", {})
    threshold_path = Path(str(threshold_record.get("path", "")))
    if not threshold_path.is_absolute():
        threshold_path = (args.contract.parent / threshold_path).resolve()
    authority_threshold = authority.get("threshold_registry", {})
    if (
        not threshold_path.is_file()
        or threshold_record.get("sha256") != sha256(threshold_path)
        or Path(str(authority_threshold.get("path", ""))).resolve()
        != threshold_path.resolve()
        or authority_threshold.get("sha256") != sha256(threshold_path)
    ):
        raise SystemExit("final release threshold registry is missing or stale")
    thresholds = load_controller_thresholds(threshold_path)
    completion_policy = thresholds["completion_policy"]
    state_policy = thresholds["state_release_policy"]
    contradiction_ceiling = thresholds[
        "observation_writeback_policy"
    ]["maximum_contradiction_fraction"]
    catalog_record = contract.get("candidate_catalog", {})
    if (
        Path(str(catalog_record.get("path", ""))).resolve()
        != args.catalog.resolve()
        or catalog_record.get("sha256") != sha256(args.catalog)
    ):
        raise SystemExit("final candidate catalog differs from annotation contract")
    bound_context = contract.get("candidate_context_evidence")
    authority_context = authority.get("context_evidence")
    if args.context_evidence:
        if (
            not bound_context
            or Path(str(bound_context.get("path", ""))).resolve()
            != args.context_evidence.resolve()
            or bound_context.get("sha256") != sha256(args.context_evidence)
            or not authority_context
            or Path(str(authority_context.get("path", ""))).resolve()
            != args.context_evidence.resolve()
            or authority_context.get("sha256") != sha256(args.context_evidence)
        ):
            raise SystemExit("final context evidence is missing, stale or unbound")
    elif bound_context or authority_context:
        raise SystemExit("final release omitted contract-bound context evidence")
    candidates = catalog_candidates(
        json.loads(args.catalog.read_text(encoding="utf-8"))
    )
    context_summary = apply_candidate_context(
        candidates,
        read_tsv(args.context_evidence) if args.context_evidence else [],
    )
    eligible_broad_labels = {
        str(candidate.get("release_broad_label", ""))
        for candidate in candidates.values()
        if str(candidate.get("candidate_role", "")) == "broad"
        and candidate_can_release(candidate)
    }
    authority_membership = authority.get("post_atlas_membership", {})
    if (
        Path(str(authority_membership.get("path", ""))).resolve()
        != args.post_atlas_membership.resolve()
        or authority_membership.get("sha256")
        != sha256(args.post_atlas_membership)
    ):
        raise SystemExit("post-Atlas membership differs from final stage authority")
    authority_review = authority.get("prerequisite_manifest", {})
    if (
        Path(str(authority_review.get("path", ""))).resolve()
        != args.atlas_completeness_manifest.resolve()
        or authority_review.get("sha256")
        != sha256(args.atlas_completeness_manifest)
    ):
        raise SystemExit("Atlas/completeness review differs from final stage authority")
    review = json.loads(
        args.atlas_completeness_manifest.read_text(encoding="utf-8")
    )
    if review.get("status") != "PASS" or review.get("phase") != "atlas_and_completeness_review":
        raise SystemExit("Atlas/completeness review is not PASS")
    membership_record = review.get("membership", {})
    if (
        Path(str(membership_record.get("path", ""))).resolve()
        != args.post_atlas_membership.resolve()
        or membership_record.get("sha256") != sha256(args.post_atlas_membership)
    ):
        raise SystemExit("post-Atlas membership differs from the reviewed membership")
    completeness_record = review.get("completeness", {})
    completeness_path = Path(str(completeness_record.get("path", "")))
    if (
        not completeness_path.is_file()
        or completeness_record.get("sha256") != sha256(completeness_path)
    ):
        raise SystemExit("post-merge completeness audit is missing or stale")
    completeness = json.loads(completeness_path.read_text(encoding="utf-8"))
    if completeness.get("status") != "PASS":
        raise SystemExit("post-merge biological completeness audit is not PASS")
    reviewed_membership = completeness.get("membership", {})
    if (
        Path(str(reviewed_membership.get("path", ""))).resolve()
        != args.post_atlas_membership.resolve()
        or reviewed_membership.get("sha256") != sha256(args.post_atlas_membership)
    ):
        raise SystemExit("completeness audit does not bind post-Atlas membership")
    atlas_record = review.get("atlas_validation", {})
    atlas_path = Path(str(atlas_record.get("path", "")))
    if not atlas_path.is_file() or atlas_record.get("sha256") != sha256(atlas_path):
        raise SystemExit("Atlas validation is missing or stale")
    if json.loads(atlas_path.read_text(encoding="utf-8")).get("status") != "PASS":
        raise SystemExit("Atlas validation is not PASS")

    rows = read_tsv(args.post_atlas_membership)
    ids = [str(row.get("cell_id", "")) for row in rows]
    if not ids or "" in ids or len(ids) != len(set(ids)):
        raise SystemExit("post-Atlas membership must contain unique cell IDs")

    fine: dict[str, dict[str, str]] = {}
    if args.fine_assignments:
        authority_fine = authority.get("fine_assignments", {})
        if (
            Path(str(authority_fine.get("path", ""))).resolve()
            != args.fine_assignments.resolve()
            or authority_fine.get("sha256") != sha256(args.fine_assignments)
        ):
            raise SystemExit("fine assignments differ from final stage authority")
        for row in read_tsv(args.fine_assignments):
            cell = str(row.get("cell_id", ""))
            if (
                not cell or cell in fine
                or row.get("confidence") != "high"
                or not row.get("parent_broad_label")
                or not row.get("final_fine_label")
            ):
                raise SystemExit("fine assignments are duplicated or not release-grade")
            candidate = candidates.get(str(row.get("candidate_id", "")), {})
            if (
                not candidate_can_release(candidate)
                or str(candidate.get("release_broad_label", ""))
                != str(row.get("parent_broad_label", ""))
                or str(candidate.get("release_fine_label", ""))
                != str(row.get("final_fine_label", ""))
            ):
                raise SystemExit(
                    "fine assignment targets a context-ineligible or mismatched candidate"
                )
            fine[cell] = row

    state_by_group: dict[tuple[str, str], set[str]] = defaultdict(set)
    authority_states = {
        (Path(str(record.get("path", ""))).resolve(), str(record.get("sha256", "")))
        for record in authority.get("state_annotation_proposals", [])
    }
    supplied_states = {
        (path.resolve(), sha256(path)) for path in args.state_proposals
    }
    if authority_states != supplied_states:
        raise SystemExit("state proposals differ from final stage authority")
    for path in args.state_proposals:
        for row in read_tsv(path):
            key = (str(row.get("cohort_id", "")), str(row.get("subcluster_id", "")))
            value = str(row.get("state_annotation", ""))
            if value and (
                row.get("assignment_scope")
                != "whole_high_purity_second_round_subcluster"
                or float(row.get("lineage_supported_fraction", 0) or 0)
                < state_policy["minimum_parent_lineage_supported_fraction"]
                or float(row.get("contradiction_fraction", 1) or 1)
                > contradiction_ceiling
            ):
                raise SystemExit(
                    "state proposal is not a high-purity parent-resolved subcluster state"
                )
            if all(key) and value:
                state_by_group[key].add(value)

    final: list[dict[str, object]] = []
    for source in sorted(rows, key=lambda row: str(row.get("cell_id", ""))):
        row = dict(source)
        cell = str(row["cell_id"])
        broad = str(row.get("final_broad_label", ""))
        if broad and broad not in eligible_broad_labels:
            raise SystemExit(
                f"final broad label is not eligible under bound context: {broad}"
            )
        assignment = fine.get(cell)
        if assignment:
            if broad != assignment.get("parent_broad_label"):
                raise SystemExit("fine assignment crosses the frozen broad parent")
            row["final_fine_label"] = assignment["final_fine_label"]
            row["final_state"] = "defined_fine"
            row["final_fine_confidence"] = "high"
            row["final_fine_eligible"] = "true"
            row["fine_anchor_eligible"] = "true"
            row["final_fine_candidate_id"] = assignment.get("candidate_id", "")
        else:
            row["final_fine_label"] = ""
            row["final_fine_confidence"] = ""
            row["final_fine_eligible"] = "false"
            row["fine_anchor_eligible"] = "false"
            row["final_fine_candidate_id"] = ""
        states = state_by_group.get(
            (str(row.get("source_boundary", "")), str(row.get("source_cluster", ""))),
            set(),
        )
        row["state_annotations"] = ";".join(sorted(states))
        if not broad:
            if row.get("final_state") != "unresolved_biological":
                raise SystemExit("unlabeled final member lacks unresolved biological state")
            row["final_state"] = "qc_holdout"
            row["qc_reason"] = (
                row.get("unresolved_reason")
                or "unresolved_after_second_round_and_post_merge_atlas"
            )
            row["assignment_origin"] = "final_typed_residual_qc"
            row["confidence"] = "low"
        elif not assignment:
            row["final_state"] = "defined_broad_only"
            row["qc_reason"] = ""
        final.append(row)

    qc_n = sum(row["final_state"] == "qc_holdout" for row in final)
    qc_fraction = qc_n / len(final)
    if (
        qc_n >= int(completion_policy["residual_qc_count_trigger"])
        or qc_fraction >= float(
            completion_policy["residual_qc_fraction_trigger"]
        )
    ):
        raise SystemExit(
            "residual QC reaches the bound completion threshold; return to the implicated "
            "second-round cohort or post-merge unresolved review"
        )

    broad_census = Counter(str(row.get("final_broad_label", "")) for row in final if row.get("final_broad_label"))
    fine_census = Counter(str(row.get("final_fine_label", "")) for row in final if row.get("final_fine_label"))
    state_census = Counter(
        state for row in final
        for state in str(row.get("state_annotations", "")).split(";") if state
    )
    qc_reason_census = Counter(str(row.get("qc_reason", "")) for row in final if row.get("qc_reason"))
    args.out.mkdir(parents=True, exist_ok=True)
    membership_path = args.out / "final_release_membership.tsv.gz"
    write_tsv(membership_path, final)
    manifest = {
        "status": "PASS",
        "schema_version": "2.2",
        "phase": "materialize_final_release",
        "controller_version": "2.2.0",
        "threshold_registry": {
            "path": str(threshold_path.resolve()),
            "sha256": sha256(threshold_path),
        },
        "n_analysis_set": len(final),
        "residual_qc_n": qc_n,
        "residual_qc_fraction": qc_fraction,
        "broad_census": dict(sorted(broad_census.items())),
        "fine_census": dict(sorted(fine_census.items())),
        "state_census": dict(sorted(state_census.items())),
        "qc_reason_census": dict(sorted(qc_reason_census.items())),
        "membership": {
            "path": str(membership_path.resolve()),
            "sha256": sha256(membership_path),
            "semantic_sha256": deterministic_membership_hash(final),
        },
        "atlas_completeness_review": {
            "path": str(args.atlas_completeness_manifest.resolve()),
            "sha256": sha256(args.atlas_completeness_manifest),
        },
        "formal_broad_fine_state_qc_membership_written": True,
        "context_evidence": (
            {
                "path": str(args.context_evidence.resolve()),
                "sha256": sha256(args.context_evidence),
            }
            if args.context_evidence else None
        ),
        "context_release_eligibility": context_summary,
    }
    manifest_path = args.out / "final_release_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
