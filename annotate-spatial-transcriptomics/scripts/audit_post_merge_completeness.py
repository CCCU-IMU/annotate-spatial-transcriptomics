#!/usr/bin/env python3
"""Audit broad/fine census only after second-round broad merge and Atlas rescue."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from evidence_schema_lib import sha256
from lineage_controller_lib import (
    apply_candidate_context, candidate_can_release, catalog_candidates,
    deterministic_cell_id_set_hash, deterministic_membership_hash,
    group_candidate_detected, independent_group_program, read_tsv, write_tsv,
)
from membership_transform_lib import load_and_validate_chain


def detected(row: dict[str, str], candidate: dict) -> bool:
    return group_candidate_detected(row, candidate)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument(
        "--cluster-evidence", required=True, action="append", type=Path,
        help="repeat once for every completed second-round cohort",
    )
    ap.add_argument("--fine-audit", action="append", type=Path, default=[])
    ap.add_argument("--unmodeled", action="append", type=Path, default=[])
    ap.add_argument("--unmodeled-review", type=Path)
    ap.add_argument(
        "--local-subset-validation", action="append", type=Path, default=[],
    )
    ap.add_argument(
        "--local-remainder-audit", action="append", type=Path, default=[],
    )
    ap.add_argument(
        "--post-merge-review-manifest", type=Path,
        help=(
            "canonical post-merge component review whose accepted cell-level "
            "decisions can provide source support for minority local programs"
        ),
    )
    ap.add_argument(
        "--follicle-roi-repair-manifest", type=Path,
        help=(
            "canonical bounded follicle-ROI repair inserted after the post-merge "
            "unresolved review and before any per-broad review patches"
        ),
    )
    ap.add_argument(
        "--catalog-wide-review-manifest", action="append", type=Path, default=[],
        help="validated catalog-wide review apply manifests in chronological order",
    )
    ap.add_argument(
        "--catalog-wide-review-summary", type=Path,
        help=(
            "current canonical PASS one-pass catalog-wide review summary; "
            "binds exact zero-census absence/refutation decisions"
        ),
    )
    ap.add_argument(
        "--membership-transform-chain", type=Path,
        help=(
            "canonical ordered transform chain ending at --membership; preferred "
            "over stage-specific provenance arguments"
        ),
    )
    ap.add_argument(
        "--defer-canonical-zero-to-biological-review",
        action="store_true",
        help=(
            "permit a zero-census canonical-cluster candidate to remain "
            "not-evaluable until the contract-required object-level "
            "biological-quality review; never use without that downstream review"
        ),
    )
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    membership = read_tsv(args.membership)
    ids = [str(row.get("cell_id", "")) for row in membership]
    if not ids or "" in ids or len(ids) != len(set(ids)):
        raise SystemExit("post-merge membership must contain unique nonempty cell_id")
    broad_census = Counter(
        str(row.get("final_broad_label", "")) for row in membership
        if str(row.get("final_broad_label", ""))
    )
    atlas_rescue_census = Counter(
        str(row.get("final_broad_label", "")) for row in membership
        if str(row.get("final_broad_label", ""))
        and str(row.get("assignment_origin", ""))
        == "post_merge_atlas_unlabeled_broad_rescue"
    )

    catalog_doc = json.loads(args.catalog.read_text(encoding="utf-8"))
    candidates = catalog_candidates(catalog_doc)
    context_summary = apply_candidate_context(
        candidates,
        read_tsv(args.context_evidence) if args.context_evidence else [],
    )
    candidate_universe = set(candidates)
    broad_candidates: dict[str, set[str]] = defaultdict(set)
    broad_release_candidates: dict[str, set[str]] = defaultdict(set)
    fine_candidates: dict[str, set[str]] = defaultdict(set)
    for candidate_id, candidate in candidates.items():
        role = str(candidate.get("candidate_role", ""))
        broad = str(candidate.get("release_broad_label", ""))
        if role in {"broad", "fine"} and broad:
            broad_candidates[broad].add(candidate_id)
        if role == "broad" and broad:
            broad_release_candidates[broad].add(candidate_id)
        if role == "fine" and broad:
            fine_candidates[broad].add(candidate_id)

    evidence_by_candidate: dict[str, list[dict[str, str]]] = defaultdict(list)
    evidence_by_group_candidate: dict[
        tuple[str, str, str], dict[str, str]
    ] = {}
    cohort_boundaries: set[str] = set()
    underpowered_cohort_n = 0
    for path in args.cluster_evidence:
        all_rows = read_tsv(path)
        rows = [
            row for row in all_rows
            if row.get("resolution_role") == "selected"
        ]
        if not rows:
            if all_rows and all(
                row.get("evaluation_status") == "underpowered_not_evaluable"
                for row in all_rows
            ):
                underpowered_cohort_n += 1
                cohort_boundaries.update(
                    str(row.get("source_boundary", "")) for row in all_rows
                    if str(row.get("source_boundary", ""))
                )
                continue
            raise SystemExit(f"selected second-round evidence is empty: {path}")
        by_group: dict[tuple[str, str], set[str]] = defaultdict(set)
        for row in rows:
            boundary = str(row.get("source_boundary", ""))
            cluster = str(row.get("source_cluster", ""))
            candidate_id = str(row.get("candidate_id", ""))
            if not boundary or not cluster:
                raise SystemExit("second-round evidence lacks boundary/cluster")
            by_group[(boundary, cluster)].add(candidate_id)
            evidence_by_candidate[candidate_id].append(row)
            key = (boundary, cluster, candidate_id)
            if key in evidence_by_group_candidate:
                raise SystemExit(
                    "selected second-round evidence duplicates a group x candidate"
                )
            evidence_by_group_candidate[key] = row
            cohort_boundaries.add(boundary)
        incomplete = [group for group, observed in by_group.items() if observed != candidate_universe]
        if incomplete:
            raise SystemExit(
                f"second-round full-catalog scan is incomplete for {len(incomplete)} subclusters"
            )

    errors: list[str] = []
    supported_local_subsets: set[tuple[str, str, str]] = set()
    for path in args.local_subset_validation:
        for row in read_tsv(path):
            if str(row.get("status", "")) != "PASS":
                continue
            key = (
                str(row.get("source_boundary", "")),
                str(row.get("source_cluster", "")),
                str(row.get("candidate_id", "")),
            )
            if all(key):
                supported_local_subsets.add(key)
    supported_local_remainders: set[tuple[str, str, str]] = set()
    for path in args.local_remainder_audit:
        for row in read_tsv(path):
            if str(row.get("status", "")) != "PASS":
                continue
            key = (
                str(row.get("source_boundary", "")),
                str(row.get("source_cluster", "")),
                str(row.get("candidate_id", "")),
            )
            if all(key):
                supported_local_remainders.add(key)
    supported_transform_cells: set[tuple[str, str, str]] = set()
    supported_transform_programs: set[tuple[str, str, str]] = set()
    transform_operation_census: Counter[str] = Counter()
    transform_chain = None
    if args.membership_transform_chain:
        transform_chain = load_and_validate_chain(
            args.membership_transform_chain, args.membership.resolve()
        )
        for entry in transform_chain.get("transforms", []):
            operation = str(entry.get("operation", ""))
            transform_operation_census[operation] += 1
            if operation in {
                "atlas_unlabeled_broad_rescue", "source_unit_sync",
                "final_release_materialization",
            }:
                continue
            delta_path = Path(str(entry.get("delta", {}).get("path", "")))
            transform_id = str(entry.get("transform_id", ""))
            for row in read_tsv(delta_path):
                cell_id = str(row.get("cell_id", ""))
                candidate_id = str(row.get("candidate_id", ""))
                broad = str(row.get("new_broad_label", ""))
                if cell_id and candidate_id and broad:
                    supported_transform_cells.add(
                        (cell_id, candidate_id, broad)
                    )
                    supported_transform_programs.add(
                        (transform_id, candidate_id, broad)
                    )
    supported_post_merge_cells: set[tuple[str, str, str]] = set()
    supported_post_merge_programs: set[tuple[str, str, str]] = set()
    post_merge_review_membership_path: Path | None = None
    if args.post_merge_review_manifest:
        review = json.loads(
            args.post_merge_review_manifest.read_text(encoding="utf-8")
        )
        if (
            review.get("status") != "PASS"
            or review.get("stage") != "post_merge_unresolved_component_review"
        ):
            raise SystemExit("post-merge component review is not canonical PASS")
        membership_record = review.get("membership", {})
        post_merge_review_membership_path = Path(
            str(membership_record.get("path", ""))
        )
        if (
            not post_merge_review_membership_path.is_file()
            or membership_record.get("sha256")
            != sha256(post_merge_review_membership_path)
        ):
            raise SystemExit("post-merge review membership is missing or stale")
        if (
            not args.catalog_wide_review_manifest
            and post_merge_review_membership_path.resolve()
            != args.membership.resolve()
        ):
            raise SystemExit("post-merge review does not bind audited membership")
        decision_record = review.get("component_artifacts", {}).get(
            "candidate_component_decisions.tsv", {}
        )
        decision_path = Path(str(decision_record.get("path", "")))
        if (
            not decision_path.is_file()
            or decision_record.get("sha256") != sha256(decision_path)
        ):
            raise SystemExit("post-merge component decisions are missing or stale")
        for row in read_tsv(decision_path):
            decision = str(row.get("decision", ""))
            cell_id = str(row.get("cell_id", ""))
            candidate_id = str(row.get("candidate_id", ""))
            broad = str(row.get("release_broad_label", ""))
            if (
                decision
                and not decision.startswith("remain_unresolved")
                and cell_id and candidate_id and broad
            ):
                if not candidate_can_release(candidates.get(candidate_id, {})):
                    errors.append(
                        f"{candidate_id}: post-merge review released a context-ineligible candidate"
                    )
                    continue
                supported_post_merge_cells.add((cell_id, candidate_id, broad))
                supported_post_merge_programs.add((
                    str(row.get("component_id", "")), candidate_id, broad
                ))
    supported_follicle_roi_cells: set[tuple[str, str, str]] = set()
    supported_follicle_roi_programs: set[tuple[str, str, str]] = set()
    follicle_roi_membership_path: Path | None = None
    if args.follicle_roi_repair_manifest:
        repair = json.loads(
            args.follicle_roi_repair_manifest.read_text(encoding="utf-8")
        )
        if (
            repair.get("stage") != "follicle_roi_repair_apply"
            or repair.get("status")
            != "PENDING_POST_REPAIR_BIOLOGICAL_REVIEW"
        ):
            raise SystemExit("follicle ROI repair manifest is not canonical")
        authority_record = repair.get("stage_authority", {})
        authority_path = Path(str(authority_record.get("path", "")))
        if (
            not authority_path.is_file()
            or authority_record.get("sha256") != sha256(authority_path)
        ):
            raise SystemExit("follicle ROI repair stage authority is missing or stale")
        authority = json.loads(authority_path.read_text(encoding="utf-8"))
        writer = Path(__file__).resolve().parent / "apply_sheep_ovary_follicle_roi_repair.py"
        writer_record = authority.get("scripts", {}).get(writer.name, {})
        if (
            authority.get("mode") != "stage_authority"
            or authority.get("phase") != "atlas_and_completeness_review"
            or Path(str(writer_record.get("path", ""))).resolve()
            != writer.resolve()
            or writer_record.get("sha256") != sha256(writer)
        ):
            raise SystemExit("follicle ROI repair lacks canonical writer authority")
        source_record = repair.get("pre_repair_membership", {})
        result_record = repair.get("repaired_membership", {})
        source_path = Path(str(source_record.get("path", "")))
        follicle_roi_membership_path = Path(str(result_record.get("path", "")))
        if (
            not source_path.is_file()
            or source_record.get("sha256") != sha256(source_path)
            or not follicle_roi_membership_path.is_file()
            or result_record.get("sha256") != sha256(follicle_roi_membership_path)
            or (
                post_merge_review_membership_path is not None
                and source_path.resolve()
                != post_merge_review_membership_path.resolve()
            )
        ):
            raise SystemExit("follicle ROI repair membership chain is missing or stale")
        changes_record = repair.get("changes", {})
        changes_path = Path(str(changes_record.get("path", "")))
        if (
            not changes_path.is_file()
            or changes_record.get("sha256") != sha256(changes_path)
        ):
            raise SystemExit("follicle ROI repair changes are missing or stale")
        repaired_rows = {
            str(row.get("cell_id", "")): row
            for row in read_tsv(follicle_roi_membership_path)
        }
        for row in read_tsv(changes_path):
            cell_id = str(row.get("cell_id", ""))
            candidate_id = str(row.get("candidate_id", ""))
            broad = str(row.get("new_broad_label", ""))
            result_row = repaired_rows.get(cell_id, {})
            if (
                not cell_id or not candidate_id or not broad
                or str(result_row.get("final_broad_label", "")) != broad
                or str(result_row.get("assignment_origin", ""))
                != "follicle_roi_raw_count_direct_identity_repair"
            ):
                raise SystemExit("follicle ROI repair delta differs from membership")
            supported_follicle_roi_cells.add((cell_id, candidate_id, broad))
            supported_follicle_roi_programs.add((
                str(row.get("follicle_roi_id", "")), candidate_id, broad
            ))
    supported_catalog_review_cells: set[tuple[str, str, str]] = set()
    supported_catalog_review_programs: set[tuple[str, str, str]] = set()
    previous_membership_path: Path | None = (
        follicle_roi_membership_path or post_merge_review_membership_path
    )
    for index, path in enumerate(args.catalog_wide_review_manifest):
        manifest = json.loads(path.read_text(encoding="utf-8"))
        if (
            manifest.get("status") != "PASS_REQUIRES_NEXT_REVIEW_ROUND"
            or manifest.get("stage") != "catalog_wide_lineage_review_apply"
        ):
            raise SystemExit("catalog-wide review apply manifest is not canonical PASS")
        authority_record = manifest.get("stage_authority", {})
        authority_path = Path(str(authority_record.get("path", "")))
        if (
            not authority_path.is_file()
            or authority_record.get("sha256") != sha256(authority_path)
        ):
            raise SystemExit("catalog-wide review stage authority is missing or stale")
        authority = json.loads(authority_path.read_text(encoding="utf-8"))
        writer = Path(__file__).resolve().parent / "apply_catalog_wide_lineage_review.py"
        writer_record = authority.get("scripts", {}).get(writer.name, {})
        if (
            authority.get("mode") != "stage_authority"
            or authority.get("phase") != "atlas_and_completeness_review"
            or Path(str(writer_record.get("path", ""))).resolve() != writer.resolve()
            or writer_record.get("sha256") != sha256(writer)
        ):
            raise SystemExit("catalog-wide review was not written by the canonical controller")
        source_record = manifest.get("source_membership", {})
        result_record = manifest.get("membership", {})
        source_path = Path(str(source_record.get("path", "")))
        result_path = Path(str(result_record.get("path", "")))
        if (
            not source_path.is_file() or source_record.get("sha256") != sha256(source_path)
            or not result_path.is_file() or result_record.get("sha256") != sha256(result_path)
            or (previous_membership_path is not None and source_path.resolve() != previous_membership_path.resolve())
        ):
            raise SystemExit("catalog-wide review membership chain is missing or stale")
        source_rows = read_tsv(source_path)
        result_rows = read_tsv(result_path)
        transform = manifest.get("membership_transform", {})
        if (
            transform.get("operation") != "catalog_wide_exact_cell_id_patch"
            or transform.get("source_physical_sha256") != sha256(source_path)
            or transform.get("result_physical_sha256") != sha256(result_path)
            or transform.get("source_semantic_sha256")
            != deterministic_membership_hash(source_rows)
            or transform.get("result_semantic_sha256")
            != deterministic_membership_hash(result_rows)
            or transform.get("source_cell_id_set_sha256")
            != deterministic_cell_id_set_hash(source_rows)
            or transform.get("result_cell_id_set_sha256")
            != deterministic_cell_id_set_hash(result_rows)
            or transform.get("source_cell_id_set_sha256")
            != transform.get("result_cell_id_set_sha256")
        ):
            raise SystemExit("catalog-wide membership transform ledger is invalid")
        previous_membership_path = result_path
        changes_record = manifest.get("changes", {})
        changes_path = Path(str(changes_record.get("path", "")))
        if not changes_path.is_file() or changes_record.get("sha256") != sha256(changes_path):
            raise SystemExit("catalog-wide review changes are missing or stale")
        if (
            transform.get("delta_physical_sha256") != sha256(changes_path)
            or int(transform.get("changed_observation_n", -1))
            != len(read_tsv(changes_path))
        ):
            raise SystemExit("catalog-wide membership delta differs from transform ledger")
        for row in read_tsv(changes_path):
            cell_id = str(row.get("cell_id", ""))
            candidate_id = str(row.get("candidate_id", ""))
            broad = str(row.get("new_broad_label", ""))
            if cell_id and candidate_id and broad:
                supported_catalog_review_cells.add((cell_id, candidate_id, broad))
                supported_catalog_review_programs.add((
                    str(row.get("review_id", "")), candidate_id, broad
                ))
    if previous_membership_path is not None and previous_membership_path.resolve() != args.membership.resolve():
        raise SystemExit("last catalog-wide review manifest does not bind audited membership")

    exact_refuted_zero_census: set[str] = set()
    catalog_wide_review_summary_record = None
    if args.catalog_wide_review_summary:
        summary_path = args.catalog_wide_review_summary.resolve()
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if (
            summary.get("status") != "PASS"
            or summary.get("stage") != "post_atlas_catalog_wide_lineage_review"
            or int(summary.get("review_queue_n", -1)) != 0
        ):
            raise SystemExit(
                "catalog-wide review summary is not a closed canonical PASS"
            )
        summary_membership = summary.get("membership", {})
        summary_membership_path = Path(str(summary_membership.get("path", "")))
        if (
            not summary_membership_path.is_file()
            or summary_membership.get("sha256") != sha256(summary_membership_path)
            or summary_membership_path.resolve() != args.membership.resolve()
        ):
            raise SystemExit(
                "catalog-wide review summary does not bind audited membership"
            )
        matrix_record = summary.get("artifacts", {}).get(
            "lineage_review_matrix", {}
        )
        matrix_path = Path(str(matrix_record.get("path", "")))
        if (
            not matrix_path.is_file()
            or matrix_record.get("sha256") != sha256(matrix_path)
        ):
            raise SystemExit("catalog-wide lineage review matrix is missing or stale")
        current_exact_absence = {
            str(row.get("broad_label", ""))
            for row in read_tsv(matrix_path)
            if str(row.get("status", "")) in {
                "closed_after_single_cell_type_review",
                "supported_after_exact_cell_type_review",
            }
            and int(row.get("final_n_observations", "0") or 0) == 0
            and int(row.get("precision_review_source_group_n", "0") or 0) == 0
            and int(row.get("recall_challenger_component_n", "0") or 0) == 0
            and int(row.get("recall_challenger_observation_n", "0") or 0) == 0
            and int(row.get("recall_group_watch_n", "0") or 0) == 0
        }
        validated_exact_absence: set[str] = set()
        for record in summary.get("prior_decision_validations", []):
            validation_path = Path(str(record.get("path", "")))
            if (
                not validation_path.is_file()
                or record.get("sha256") != sha256(validation_path)
            ):
                raise SystemExit(
                    "catalog-wide prior decision validation is missing or stale"
                )
            validation = json.loads(validation_path.read_text(encoding="utf-8"))
            if (
                validation.get("status") != "PASS"
                or validation.get("stage")
                != "catalog_wide_lineage_review_decision_validation"
            ):
                raise SystemExit(
                    "catalog-wide prior decision validation is not canonical PASS"
                )
            cell_type = str(validation.get("active_cell_type", ""))
            decisions_record = validation.get("validated_decisions", {})
            decisions_path = Path(str(decisions_record.get("path", "")))
            if (
                not decisions_path.is_file()
                or decisions_record.get("sha256") != sha256(decisions_path)
            ):
                raise SystemExit(
                    "catalog-wide validated decisions are missing or stale"
                )
            if any(
                str(row.get("outcome", ""))
                == "confirm_absent_or_not_evaluable"
                and str(row.get("whole_query_recall", ""))
                == "confirmed_absent"
                and str(row.get("outside_recall_challenger_resolution", ""))
                == "refuted_by_multichannel_evidence"
                for row in read_tsv(decisions_path)
            ):
                validated_exact_absence.add(cell_type)
        exact_refuted_zero_census = current_exact_absence & validated_exact_absence
        catalog_wide_review_summary_record = {
            "path": str(summary_path),
            "sha256": sha256(summary_path),
        }
    unsupported_release_members: Counter[str] = Counter()
    unsupported_release_groups: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for row in membership:
        broad = str(row.get("final_broad_label", ""))
        if not broad or str(row.get("assignment_origin", "")) == (
            "post_merge_atlas_unlabeled_broad_rescue"
        ):
            continue
        candidate_id = str(row.get("candidate_id", ""))
        boundary = str(row.get("source_boundary", ""))
        cluster = str(row.get("source_cluster", ""))
        candidate = candidates.get(candidate_id)
        evidence = evidence_by_group_candidate.get(
            (boundary, cluster, candidate_id)
        )
        key = (boundary, cluster, candidate_id)
        origin = str(row.get("assignment_origin", ""))
        source_supported = bool(
            evidence and detected(evidence, candidate)
        ) if candidate else False
        if origin in {"supported_subset", "supported_subset_round2"}:
            source_supported = key in supported_local_subsets
        elif origin == "exact_remainder_parent":
            source_supported = key in supported_local_remainders
        elif origin.startswith("post_merge_unresolved_component_review__"):
            source_supported = (
                str(row.get("cell_id", "")), candidate_id, broad
            ) in supported_post_merge_cells
        elif origin == "follicle_roi_raw_count_direct_identity_repair":
            source_supported = (
                str(row.get("cell_id", "")), candidate_id, broad
            ) in supported_follicle_roi_cells
        elif origin.startswith("catalog_wide_lineage_review_round_"):
            source_supported = (
                str(row.get("cell_id", "")), candidate_id, broad
            ) in supported_catalog_review_cells
        if (
            str(row.get("cell_id", "")), candidate_id, broad
        ) in supported_transform_cells:
            source_supported = True
        supported = bool(
            candidate
            and candidate_can_release(candidate)
            and str(candidate.get("release_broad_label", "")) == broad
            and source_supported
        )
        if not supported:
            unsupported_release_members[broad] += 1
            unsupported_release_groups[broad].add((boundary, cluster))

    broad_rows: list[dict[str, object]] = []
    for broad, candidate_ids in sorted(broad_candidates.items()):
        eligible_candidate_ids = {
            candidate_id for candidate_id in candidate_ids
            if candidate_can_release(candidates[candidate_id])
        }
        eligible_broad_release_ids = {
            candidate_id for candidate_id in broad_release_candidates[broad]
            if candidate_can_release(candidates[candidate_id])
        }
        context_ineligible_ids = candidate_ids - eligible_candidate_ids
        evidence = [
            row for candidate_id in eligible_candidate_ids
            for row in evidence_by_candidate.get(candidate_id, [])
        ]
        # A visible shared program is sufficient to keep a candidate in the
        # descriptive audit, but a zero census conflicts only with an
        # independent identity program.  Use the same separability rule as the
        # exact catalog-wide reviewer.
        positive = [
            row for row in evidence
            if independent_group_program(
                row, candidates[str(row.get("candidate_id", ""))]
            )
        ]
        local_positive_n = sum(
            candidate_id in eligible_candidate_ids
            for _, _, candidate_id in (
                supported_local_subsets | supported_local_remainders
            )
        )
        local_positive_n += sum(
            candidate_id in eligible_candidate_ids and component_broad == broad
            for _, candidate_id, component_broad in supported_post_merge_programs
        )
        local_positive_n += sum(
            candidate_id in eligible_candidate_ids and component_broad == broad
            for _, candidate_id, component_broad in supported_follicle_roi_programs
        )
        local_positive_n += sum(
            candidate_id in eligible_candidate_ids and component_broad == broad
            for _, candidate_id, component_broad in supported_catalog_review_programs
        )
        local_positive_n += sum(
            candidate_id in eligible_candidate_ids and component_broad == broad
            for _, candidate_id, component_broad in supported_transform_programs
        )
        positive_program_n = len(positive) + local_positive_n
        n_final = broad_census[broad]
        unsupported_n = unsupported_release_members[broad]
        if n_final and not eligible_broad_release_ids:
            status = "blocked"
            rationale = "released_broad_is_not_evaluable_under_bound_context"
            errors.append(
                f"{broad}: final membership exists although every broad release "
                "candidate is context-ineligible"
            )
        elif n_final and positive_program_n and not unsupported_n:
            status = "supported"
            rationale = "each_non_atlas_source_group_matches_its_multichannel_candidate_program"
        elif n_final and atlas_rescue_census[broad] == n_final:
            status = "supported"
            rationale = "calibrated_post_merge_atlas_rescued_underpowered_members"
        elif n_final:
            status = "blocked"
            rationale = "released_broad_contains_source_groups_without_matching_multichannel_support"
            errors.append(
                f"{broad}: {unsupported_n or n_final} released observations lack "
                "source-linked second-round support"
            )
        elif not eligible_broad_release_ids and context_ineligible_ids:
            status = "not_evaluable"
            rationale = "bound_context_does_not_permit_this_stage_dependent_lineage_evaluation"
        elif positive_program_n and broad in exact_refuted_zero_census:
            status = "refuted"
            rationale = (
                "positive_group_program_refuted_by_current_exact_cell_type_review"
            )
        elif positive_program_n:
            canonical_zero = any(
                str(candidates[candidate_id].get("writeback_strategy", ""))
                == "canonical_cluster_membership"
                for candidate_id in candidate_ids
            )
            if (
                canonical_zero
                and args.defer_canonical_zero_to_biological_review
            ):
                status = "not_evaluable"
                rationale = (
                    "canonical_zero_census_deferred_to_required_object_level_"
                    "biological_review"
                )
            else:
                status = "blocked"
                rationale = "positive_program_remains_without_broad_reconstruction"
                errors.append(
                    f"{broad}: zero census conflicts with a positive second-round program"
                )
        elif underpowered_cohort_n:
            status = "not_evaluable"
            rationale = "complete_evaluable_scan_negative_but_tiny_cohort_remained_underpowered"
        else:
            status = "refuted"
            rationale = "complete_full_catalog_scan_found_no_coherent_program"
        broad_rows.append({
            "broad_label": broad,
            "candidate_ids": ";".join(sorted(candidate_ids)),
            "release_eligible_candidate_ids": ";".join(
                sorted(eligible_candidate_ids)
            ),
            "context_not_evaluable_candidate_ids": ";".join(
                sorted(context_ineligible_ids)
            ),
            "final_n_observations": n_final,
            "selected_subcluster_evidence_n": len(evidence),
            "positive_program_n": positive_program_n,
            "local_supported_program_n": local_positive_n,
            "unsupported_release_observation_n": unsupported_n,
            "unsupported_release_group_n": len(unsupported_release_groups[broad]),
            "status": status,
            "rationale": rationale,
        })

    fine_status: dict[tuple[str, str], list[str]] = defaultdict(list)
    for path in args.fine_audit:
        for row in read_tsv(path):
            key = (
                str(row.get("parent_broad_label", "")),
                str(row.get("candidate_id", "")),
            )
            if all(key):
                fine_status[key].append(str(row.get("status", "")))
    fine_rows: list[dict[str, object]] = []
    for parent in sorted(label for label in broad_census if label in fine_candidates):
        for candidate_id in sorted(fine_candidates[parent]):
            statuses = fine_status.get((parent, candidate_id), [])
            if not candidate_can_release(candidates[candidate_id]):
                status = "not_evaluable"
            elif not statuses:
                status = "missing"
                errors.append(f"{parent}/{candidate_id}: fine candidate was not audited")
            elif "supported" in statuses:
                status = "supported"
            elif set(statuses) <= {"refuted", "not_evaluable"}:
                status = "not_evaluable" if set(statuses) == {"not_evaluable"} else "refuted"
            else:
                status = "invalid"
                errors.append(f"{parent}/{candidate_id}: fine audit has invalid statuses")
            fine_rows.append({
                "parent_broad_label": parent,
                "candidate_id": candidate_id,
                "status": status,
                "evaluated_subcluster_n": len(statuses),
                "context_status": candidates[candidate_id].get(
                    "_context_status", "not_evaluable"
                ),
            })

    unmodeled_rows: list[dict[str, str]] = []
    for path in args.unmodeled:
        manifest = json.loads(path.read_text(encoding="utf-8"))
        record = manifest.get("accepted_programs", {})
        program_path = Path(str(record.get("path", "")))
        if program_path.is_file() and record.get("sha256") == sha256(program_path):
            source_boundary = next(
                (
                    part for part in program_path.parts
                    if part.startswith("initial_cluster__")
                ),
                "",
            )
            for row in read_tsv(program_path):
                unmodeled_rows.append({
                    **row, "source_boundary": source_boundary,
                })
        elif int(manifest.get("accepted_program_n", 0) or 0):
            errors.append("unmodeled discovery manifest contains stale accepted programs")
    unmodeled_decisions = (
        read_tsv(args.unmodeled_review) if args.unmodeled_review else []
    )
    decision_by_key = {
        (str(row.get("source_boundary", "")), str(row.get("program_id", ""))): row
        for row in unmodeled_decisions
    }
    expected_unmodeled_keys = {
        (row["source_boundary"], str(row.get("program_id", "")))
        for row in unmodeled_rows
    }
    if len(decision_by_key) != len(unmodeled_decisions):
        errors.append("unmodeled program review contains duplicate or empty keys")
    if set(decision_by_key) != expected_unmodeled_keys:
        if unmodeled_rows:
            errors.append(
                "stable unmodeled lineage candidates require exact catalog review before final release"
            )
        elif unmodeled_decisions:
            errors.append("unmodeled decisions exist although no stable program was found")
    allowed_unmodeled_outcomes = {
        "catalog_program", "state_or_technical_program",
        "insufficient_identity_program", "novel_lineage_candidate",
    }
    for row in unmodeled_rows:
        decision = decision_by_key.get(
            (row["source_boundary"], str(row.get("program_id", ""))), {}
        )
        outcome = str(decision.get("outcome", ""))
        candidate_id = str(decision.get("catalog_candidate_id", ""))
        rationale = str(decision.get("rationale", "")).strip()
        row["review_outcome"] = outcome
        row["review_catalog_candidate_id"] = candidate_id
        row["review_rationale"] = rationale
        if not decision:
            continue
        if outcome not in allowed_unmodeled_outcomes:
            errors.append(
                f"{row['source_boundary']}/{row.get('program_id', '')}: invalid unmodeled review outcome"
            )
        elif len(rationale) < 20:
            errors.append(
                f"{row['source_boundary']}/{row.get('program_id', '')}: unmodeled review rationale is too short"
            )
        elif outcome == "catalog_program" and candidate_id not in candidate_universe:
            errors.append(
                f"{row['source_boundary']}/{row.get('program_id', '')}: catalog review names an unknown candidate"
            )
        elif outcome == "novel_lineage_candidate":
            errors.append(
                f"{row['source_boundary']}/{row.get('program_id', '')}: novel lineage requires catalog extension and cohort rerun"
            )

    args.out.mkdir(parents=True, exist_ok=True)
    broad_path = args.out / "broad_completeness_audit.tsv"
    fine_path = args.out / "fine_candidate_census_audit.tsv"
    unmodeled_path = args.out / "unmodeled_lineage_audit.tsv"
    write_tsv(broad_path, broad_rows)
    write_tsv(
        fine_path, fine_rows,
        fields=[
            "parent_broad_label", "candidate_id", "status",
            "evaluated_subcluster_n", "context_status",
        ],
    )
    write_tsv(unmodeled_path, unmodeled_rows)
    manifest = {
        "status": "PASS" if not errors else "BLOCKED",
        "schema_version": "2.2",
        "stage": "post_merge_completeness",
        "n_analysis_set": len(membership),
        "n_second_round_cohorts": len(cohort_boundaries),
        "full_catalog_scan": True,
        "broad_census": dict(sorted(broad_census.items())),
        "unmodeled_candidate_n": len(unmodeled_rows),
        "unmodeled_review_n": len(unmodeled_decisions),
        "underpowered_cohort_n": underpowered_cohort_n,
        "supported_post_merge_component_observation_n": len(
            supported_post_merge_cells
        ),
        "supported_follicle_roi_repair_observation_n": len(
            supported_follicle_roi_cells
        ),
        "supported_catalog_wide_review_observation_n": len(
            supported_catalog_review_cells
        ),
        "supported_membership_transform_observation_n": len(
            supported_transform_cells
        ),
        "membership_transform_operation_census": dict(
            sorted(transform_operation_census.items())
        ),
        "membership_transform_chain": (
            {
                "path": str(args.membership_transform_chain.resolve()),
                "sha256": sha256(args.membership_transform_chain),
            }
            if args.membership_transform_chain else None
        ),
        "catalog_wide_review_apply_manifests": [
            {"path": str(path.resolve()), "sha256": sha256(path)}
            for path in args.catalog_wide_review_manifest
        ],
        "catalog_wide_review_summary": catalog_wide_review_summary_record,
        "exact_refuted_zero_census_broad_labels": sorted(
            exact_refuted_zero_census
        ),
        "follicle_roi_repair_manifest": (
            {
                "path": str(args.follicle_roi_repair_manifest.resolve()),
                "sha256": sha256(args.follicle_roi_repair_manifest),
            }
            if args.follicle_roi_repair_manifest else None
        ),
        "membership": {
            "path": str(args.membership.resolve()),
            "sha256": sha256(args.membership),
        },
        "candidate_catalog": {
            "path": str(args.catalog.resolve()),
            "sha256": sha256(args.catalog),
        },
        "context_evidence": (
            {
                "path": str(args.context_evidence.resolve()),
                "sha256": sha256(args.context_evidence),
            }
            if args.context_evidence else None
        ),
        "context_release_eligibility": context_summary,
        "broad_audit": {"path": str(broad_path.resolve()), "sha256": sha256(broad_path)},
        "fine_audit": {"path": str(fine_path.resolve()), "sha256": sha256(fine_path)},
        "unmodeled_audit": {"path": str(unmodeled_path.resolve()), "sha256": sha256(unmodeled_path)},
        "errors": errors,
    }
    manifest_path = args.out / "post_merge_completeness_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
