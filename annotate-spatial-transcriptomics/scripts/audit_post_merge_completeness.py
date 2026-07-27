#!/usr/bin/env python3
"""Audit broad/fine census only after second-round broad merge and Atlas rescue."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from evidence_schema_lib import sha256
from lineage_controller_lib import (
    catalog_candidates, group_candidate_detected, read_tsv, write_tsv,
)


def detected(row: dict[str, str], candidate: dict) -> bool:
    return group_candidate_detected(row, candidate)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
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
    candidate_universe = set(candidates)
    broad_candidates: dict[str, set[str]] = defaultdict(set)
    fine_candidates: dict[str, set[str]] = defaultdict(set)
    for candidate_id, candidate in candidates.items():
        role = str(candidate.get("candidate_role", ""))
        broad = str(candidate.get("release_broad_label", ""))
        if role in {"broad", "fine"} and broad:
            broad_candidates[broad].add(candidate_id)
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
    supported_post_merge_cells: set[tuple[str, str, str]] = set()
    supported_post_merge_programs: set[tuple[str, str, str]] = set()
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
        if (
            Path(str(membership_record.get("path", ""))).resolve()
            != args.membership.resolve()
            or membership_record.get("sha256") != sha256(args.membership)
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
                supported_post_merge_cells.add((cell_id, candidate_id, broad))
                supported_post_merge_programs.add((
                    str(row.get("component_id", "")), candidate_id, broad
                ))
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
        supported = bool(
            candidate
            and str(candidate.get("release_broad_label", "")) == broad
            and source_supported
        )
        if not supported:
            unsupported_release_members[broad] += 1
            unsupported_release_groups[broad].add((boundary, cluster))

    broad_rows: list[dict[str, object]] = []
    for broad, candidate_ids in sorted(broad_candidates.items()):
        evidence = [
            row for candidate_id in candidate_ids
            for row in evidence_by_candidate.get(candidate_id, [])
        ]
        positive = [
            row for row in evidence
            if detected(row, candidates[str(row.get("candidate_id", ""))])
        ]
        local_positive_n = sum(
            candidate_id in candidate_ids
            for _, _, candidate_id in (
                supported_local_subsets | supported_local_remainders
            )
        )
        local_positive_n += sum(
            candidate_id in candidate_ids and component_broad == broad
            for _, candidate_id, component_broad in supported_post_merge_programs
        )
        positive_program_n = len(positive) + local_positive_n
        n_final = broad_census[broad]
        unsupported_n = unsupported_release_members[broad]
        if n_final and positive_program_n and not unsupported_n:
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
            if not statuses:
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
            "evaluated_subcluster_n",
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
        "membership": {
            "path": str(args.membership.resolve()),
            "sha256": sha256(args.membership),
        },
        "candidate_catalog": {
            "path": str(args.catalog.resolve()),
            "sha256": sha256(args.catalog),
        },
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
