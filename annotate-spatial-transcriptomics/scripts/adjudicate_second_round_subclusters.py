#!/usr/bin/env python3
"""Adjudicate second-round Leiden subclusters without freezing release labels."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from evidence_schema_lib import sha256
from lineage_controller_lib import (
    candidate_can_release, canonical_cluster_challenger, catalog_candidates,
    dominant_generic_remainder_group,
    effective_broad_writeback_strategy, group_candidate_detected,
    group_candidate_score, group_identity_core_fraction,
    group_identity_core_direct_fraction,
    local_split_worthy_group_program,
    group_orthogonal_support_count, group_release_supported_fraction,
    number, read_tsv, write_tsv,
)


def truth(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "pass"}


def detected(row: dict[str, str], candidate: dict | None = None) -> bool:
    return group_candidate_detected(row, candidate)


def orthogonal_support_count(
    row: dict[str, str], candidate: dict | None = None
) -> int:
    return group_orthogonal_support_count(row, candidate)


def candidate_score(row: dict[str, str], candidate: dict | None = None) -> float:
    return group_candidate_score(row, candidate)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--partitions", required=True, type=Path)
    ap.add_argument("--cluster-evidence", required=True, type=Path)
    ap.add_argument("--scores", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--contract", required=True, type=Path)
    ap.add_argument("--cohort-id", required=True)
    ap.add_argument("--source-initial-cluster", required=True)
    ap.add_argument("--provisional-status", required=True, choices=["provisional_broad", "mixed", "unknown"])
    ap.add_argument("--provisional-broad", default="")
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    partition_rows = read_tsv(args.partitions)
    selected = [row for row in partition_rows if row.get("resolution_role") == "selected"]
    if not selected:
        raise SystemExit("selected second-round partition is missing")
    members: dict[str, list[str]] = defaultdict(list)
    for row in selected:
        members[str(row.get("cluster", ""))].append(str(row.get("cell_id", "")))
    if any(not cluster or not all(cells) for cluster, cells in members.items()):
        raise SystemExit("second-round selected partition is incomplete")

    catalog_doc = json.loads(args.catalog.read_text(encoding="utf-8"))
    candidates = catalog_candidates(catalog_doc)
    context_supported: set[str] = set()
    if args.context_evidence:
        context_supported = {
            str(row.get("candidate_id", ""))
            for row in read_tsv(args.context_evidence)
            if str(row.get("status", "")).strip().lower()
            in {"supported", "pass"}
        }
    for candidate_id, candidate in candidates.items():
        context_id = str(
            candidate.get("context_evidence_candidate_id", "") or candidate_id
        )
        candidate["_context_ok"] = (
            candidate_id in context_supported or context_id in context_supported
        )
    evidence_rows = [
        row for row in read_tsv(args.cluster_evidence)
        if row.get("resolution_role") == "selected"
    ]
    by_cluster: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in evidence_rows:
        by_cluster[str(row.get("source_cluster", ""))].append(row)
    score_rows = read_tsv(args.scores)
    candidate_universe = sorted(candidates)
    per_cell_candidates: dict[str, set[str]] = defaultdict(set)
    for row in score_rows:
        per_cell_candidates[str(row.get("cell_id", ""))].add(str(row.get("candidate_id", "")))
    selected_ids = {cell for cells in members.values() for cell in cells}
    if (
        set(per_cell_candidates) != selected_ids
        or any(value != set(candidate_universe) for value in per_cell_candidates.values())
    ):
        raise SystemExit("second-round scorer did not expose the complete candidate catalog")
    if any(
        {str(row.get("candidate_id", "")) for row in by_cluster.get(cluster, [])}
        != set(candidate_universe)
        for cluster in members
    ):
        raise SystemExit("second-round group evidence does not cover every subcluster x candidate")

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    policy = contract.get("observation_writeback", {}).get("policy", {})
    required_policy = {
        "whole_subcluster_min_raw_two_family_supported_fraction",
        "whole_subcluster_min_raw_two_family_margin",
        "maximum_contradiction_fraction",
        "whole_subcluster_dominant_seed_fraction",
        "whole_subcluster_dominant_direct_core_fraction",
        "whole_subcluster_dominant_identity_core_fraction",
        "whole_subcluster_dominant_max_contradiction_fraction",
    }
    if not required_policy.issubset(policy):
        raise SystemExit(
            "annotation contract lacks canonical observation-writeback thresholds"
        )
    support_floor = number(
        policy["whole_subcluster_min_raw_two_family_supported_fraction"]
    )
    margin_floor = number(policy["whole_subcluster_min_raw_two_family_margin"])
    contradiction_ceiling = number(policy["maximum_contradiction_fraction"])
    dominant_seed_floor = number(
        policy["whole_subcluster_dominant_seed_fraction"]
    )
    dominant_direct_floor = number(
        policy["whole_subcluster_dominant_direct_core_fraction"]
    )
    dominant_core_floor = number(
        policy["whole_subcluster_dominant_identity_core_fraction"]
    )
    dominant_contradiction_ceiling = number(
        policy["whole_subcluster_dominant_max_contradiction_fraction"]
    )

    outcome_rows: list[dict[str, object]] = []
    candidate_membership: list[dict[str, object]] = []
    pending_membership: list[dict[str, object]] = []
    fine_rows: list[dict[str, object]] = []
    state_rows: list[dict[str, object]] = []
    for cluster in sorted(members, key=lambda value: (number(value, float("inf")), value)):
        rows = by_cluster.get(cluster, [])
        broad_by_label: dict[str, list[tuple[float, str, dict[str, str], dict]]] = defaultdict(list)
        for row in rows:
            candidate = candidates.get(str(row.get("candidate_id", "")), {})
            role = str(candidate.get("candidate_role", ""))
            if role not in {"broad", "fine"} or not candidate_can_release(candidate):
                continue
            label = str(candidate.get("release_broad_label", ""))
            if label:
                broad_by_label[label].append(
                    (candidate_score(row, candidate), label, row, candidate)
                )
        broad_rows = [
            sorted(items, key=lambda item: (-item[0], str(item[2].get("candidate_id", ""))))[0]
            for _, items in sorted(broad_by_label.items())
        ]
        broad_rows.sort(key=lambda item: (-item[0], str(item[2].get("candidate_id", ""))))
        positive = [
            item for item in broad_rows if detected(item[2], item[3])
        ]
        specific_positive = [
            item for item in positive
            if str(item[3].get("candidate_id", "")) != "stromal_mesenchymal"
        ]
        specific_split_worthy = [
            item for item in specific_positive
            if local_split_worthy_group_program(item[2], item[3])
        ]
        generic_positive = [
            item for item in positive
            if str(item[3].get("candidate_id", "")) == "stromal_mesenchymal"
        ]
        winner = (
            specific_split_worthy[0]
            if specific_split_worthy
            else generic_positive[0]
            if generic_positive
            else specific_positive[0]
            if specific_positive
            else positive[0]
            if positive
            else None
        )
        competing = [item for item in positive if item is not winner]
        target_raw = (
            group_release_supported_fraction(winner[2], winner[3])
            if winner else 0.0
        )
        contradiction = number(winner[2].get("hard_contradiction_fraction")) if winner else 1.0
        split_worthy_competitors = []
        for item in competing:
            row, candidate = item[2], item[3]
            if local_split_worthy_group_program(row, candidate):
                split_worthy_competitors.append(item)
        runner_raw = max(
            (
                group_release_supported_fraction(item[2], item[3])
                for item in split_worthy_competitors
            ),
            default=0.0,
        )
        whole_strategy_ok = bool(
            winner
            and effective_broad_writeback_strategy(winner[3])
            != "candidate_local_component_never_parent_expansion"
        )
        strict_whole_pass = bool(
            winner
            and whole_strategy_ok
            and target_raw >= support_floor
            and target_raw - runner_raw >= margin_floor
            and contradiction <= contradiction_ceiling
            and orthogonal_support_count(winner[2], winner[3]) >= 3
            and not split_worthy_competitors
        )
        dominant_whole_pass = bool(
            winner
            and whole_strategy_ok
            and str(winner[3].get("candidate_id", ""))
            not in {"stromal_mesenchymal"}
            and not canonical_cluster_challenger(winner[2], winner[3])
            and number(winner[2].get("observation_seed_fraction"))
            >= dominant_seed_floor
            and group_identity_core_direct_fraction(winner[2])
            >= dominant_direct_floor
            and group_identity_core_fraction(winner[2]) >= dominant_core_floor
            and contradiction <= dominant_contradiction_ceiling
            and orthogonal_support_count(winner[2], winner[3]) >= 3
            and not split_worthy_competitors
        )
        generic_remainder_whole_pass = bool(
            winner
            and whole_strategy_ok
            and dominant_generic_remainder_group(winner[2], winner[3])
            and not split_worthy_competitors
            and not specific_split_worthy
        )
        whole_pass = (
            strict_whole_pass
            or dominant_whole_pass
            or generic_remainder_whole_pass
        )
        if whole_pass and not split_worthy_competitors:
            target_label = winner[1]
            target_candidate = str(winner[2].get("candidate_id", ""))
            if args.provisional_status == "unknown" or not args.provisional_broad:
                outcome = "missing_broad_reconstruction"
            elif target_label == args.provisional_broad:
                outcome = "parent_return"
            else:
                outcome = "cross_lineage_return"
            for cell in sorted(members[cluster]):
                candidate_membership.append({
                    "cell_id": cell,
                    "source_boundary": args.cohort_id,
                    "source_cluster": cluster,
                    "candidate_id": target_candidate,
                    "proposed_state": "broad_candidate",
                    "proposed_broad_label": target_label,
                    "confidence": "high",
                    "assignment_origin": "second_round_whole_subcluster",
                    "unresolved_reason": "",
                })
            local_split = False
            whole_subcluster_route = (
                "strict_purity"
                if strict_whole_pass
                else "dominant_identity"
                if dominant_whole_pass
                else "dominant_generic_remainder"
            )
        elif winner and specific_split_worthy:
            outcome = "local_split_required"
            target_label = ""
            target_candidate = ""
            local_split = True
            whole_subcluster_route = "local_split_specific_programs"
            for cell in sorted(members[cluster]):
                pending_membership.append({
                    "cell_id": cell,
                    "source_boundary": args.cohort_id,
                    "source_cluster": cluster,
                    "pending_reason": "competing_separable_lineage_programs",
                })
        else:
            outcome = "unresolved_biological"
            target_label = ""
            target_candidate = ""
            local_split = False
            whole_subcluster_route = "unresolved"
            for cell in sorted(members[cluster]):
                candidate_membership.append({
                    "cell_id": cell,
                    "source_boundary": args.cohort_id,
                    "source_cluster": cluster,
                    "candidate_id": "",
                    "proposed_state": "unresolved_biological",
                    "proposed_broad_label": "",
                    "confidence": "low",
                    "assignment_origin": "second_round_unresolved",
                    "unresolved_reason": "insufficient_coherent_group_identity",
                })
        outcome_rows.append({
            "cohort_id": args.cohort_id,
            "source_initial_cluster": args.source_initial_cluster,
            "subcluster_id": cluster,
            "n_observations": len(members[cluster]),
            "outcome": outcome,
            "target_candidate_id": target_candidate,
            "target_broad_label": target_label,
            "provisional_broad_after_score_freeze": args.provisional_broad,
            "local_split_required": str(local_split).lower(),
            "competing_candidate_ids": ";".join(str(item[2].get("candidate_id", "")) for item in positive),
            "target_raw_supported_fraction": target_raw,
            "strongest_competitor_raw_fraction": runner_raw,
            "contradiction_fraction": contradiction,
            "whole_subcluster_route": whole_subcluster_route,
        })
        fine_by_parent: dict[str, list[tuple[float, dict[str, str], dict]]] = defaultdict(list)
        for row in rows:
            candidate = candidates.get(str(row.get("candidate_id", "")), {})
            role = str(candidate.get("candidate_role", ""))
            if role == "fine":
                parent = str(candidate.get("release_broad_label", ""))
                if parent:
                    fine_by_parent[parent].append(
                        (candidate_score(row, candidate), row, candidate)
                    )
            elif (
                role == "state" and detected(row, candidate)
                and whole_pass and not local_split and bool(target_label)
                and (
                    not str(candidate.get("parent_broad_label", ""))
                    or str(candidate.get("parent_broad_label", ""))
                    == target_label
                )
                and group_identity_core_fraction(row)
                >= support_floor
                and number(row.get("hard_contradiction_fraction"))
                <= contradiction_ceiling
                and orthogonal_support_count(row, candidate) >= 3
            ):
                state_rows.append({
                    "cohort_id": args.cohort_id,
                    "subcluster_id": cluster,
                    "candidate_id": row.get("candidate_id", ""),
                    "state_annotation": candidate.get("release_state_label", "") or row.get("candidate_id", ""),
                    "lineage_supported_fraction": group_identity_core_fraction(row),
                    "contradiction_fraction": row.get("hard_contradiction_fraction", ""),
                    "assignment_scope": "whole_high_purity_second_round_subcluster",
                })
        for parent, items in sorted(fine_by_parent.items()):
            ordered = sorted(
                items,
                key=lambda item: (-item[0], str(item[1].get("candidate_id", ""))),
            )
            evaluable = [
                item for item in ordered
                if number(item[1].get("available_positive_family_count")) >= 2
                and candidate_can_release(item[2])
            ]
            fine_positive = [
                item for item in evaluable if detected(item[1], item[2])
            ]
            fine_winner = fine_positive[0] if fine_positive else None
            fine_runner_raw = max(
                (
                    group_release_supported_fraction(item[1], item[2])
                    for item in fine_positive[1:]
                ),
                default=0.0,
            )
            fine_target_raw = (
                group_release_supported_fraction(
                    fine_winner[1], fine_winner[2]
                )
                if fine_winner else 0.0
            )
            fine_pass = bool(
                fine_winner
                and whole_pass
                and not local_split
                and parent == target_label
                and fine_target_raw >= support_floor
                and fine_target_raw - fine_runner_raw >= margin_floor
                and number(fine_winner[1].get("hard_contradiction_fraction"))
                <= contradiction_ceiling
                and number(
                    fine_winner[1].get(
                        "observation_release_family_coherent_fraction"
                    )
                ) >= 0.03
                and orthogonal_support_count(
                    fine_winner[1], fine_winner[2]
                ) >= 3
            )
            winner_id = str(fine_winner[1].get("candidate_id", "")) if fine_winner else ""
            for _, row, candidate in ordered:
                candidate_id = str(row.get("candidate_id", ""))
                if not whole_pass or local_split or parent != target_label:
                    status = "not_evaluable"
                    reason = "broad_parent_not_stably_resolved_for_this_subcluster"
                elif number(row.get("available_positive_family_count")) < 2:
                    status = "not_evaluable"
                    reason = "fewer_than_two_marker_families_available"
                elif fine_pass and candidate_id == winner_id:
                    status = "supported"
                    reason = "stable_parent_locked_fine_program"
                else:
                    status = "refuted"
                    reason = "not_the_stable_parent_specific_fine_winner"
                fine_rows.append({
                    "cohort_id": args.cohort_id,
                    "subcluster_id": cluster,
                    "candidate_id": candidate_id,
                    "parent_broad_label": parent,
                    "proposed_fine_label": candidate.get("release_fine_label", ""),
                    "status": status,
                    "release_candidate": str(status == "supported").lower(),
                    "lineage_supported_fraction": group_identity_core_fraction(row),
                    "strongest_competing_fraction": fine_runner_raw if candidate_id == winner_id else fine_target_raw,
                    "contradiction_fraction": row.get("hard_contradiction_fraction", ""),
                    "reason": reason,
                })

    args.out.mkdir(parents=True, exist_ok=True)
    outcome_path = args.out / "second_round_subcluster_outcomes.tsv"
    write_tsv(outcome_path, outcome_rows)
    candidate_path = args.out / "base_candidate_membership.tsv.gz"
    write_tsv(
        candidate_path, candidate_membership,
        fields=[
            "cell_id", "source_boundary", "source_cluster", "candidate_id",
            "proposed_state", "proposed_broad_label", "confidence",
            "assignment_origin", "unresolved_reason",
        ],
    )
    pending_path = args.out / "pending_local_split_membership.tsv.gz"
    write_tsv(
        pending_path, pending_membership,
        fields=["cell_id", "source_boundary", "source_cluster", "pending_reason"],
    )
    write_tsv(
        args.out / "fine_candidate_proposals.tsv", fine_rows,
        fields=[
            "cohort_id", "subcluster_id", "candidate_id",
            "parent_broad_label", "proposed_fine_label", "status",
            "release_candidate", "lineage_supported_fraction",
            "strongest_competing_fraction", "contradiction_fraction",
            "reason",
        ],
    )
    write_tsv(
        args.out / "state_annotation_proposals.tsv", state_rows,
        fields=[
            "cohort_id", "subcluster_id", "candidate_id", "state_annotation",
            "lineage_supported_fraction", "contradiction_fraction",
            "assignment_scope",
        ],
    )
    manifest = {
        "status": "LOCAL_SPLIT_REQUIRED" if pending_membership else "PASS",
        "schema_version": "2.2",
        "stage": "cluster_cohort_recluster",
        "cohort_id": args.cohort_id,
        "source_initial_cluster": args.source_initial_cluster,
        "formal_membership_written": False,
        "full_catalog_scan": True,
        "provisional_broad_visible_during_scoring": False,
        "n_observations": len(selected),
        "n_subclusters": len(members),
        "n_pending_local_split": len(pending_membership),
        "outcomes": {"path": str(outcome_path.resolve()), "sha256": sha256(outcome_path)},
        "base_candidate_membership": {"path": str(candidate_path.resolve()), "sha256": sha256(candidate_path)},
        "pending_local_split_membership": {"path": str(pending_path.resolve()), "sha256": sha256(pending_path)},
    }
    (args.out / "second_round_adjudication_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
