#!/usr/bin/env python3
"""Close exact local remainders as proposals; never write release membership."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from controller_thresholds import (
    load_controller_thresholds, observation_writeback_defaults,
)
from lineage_controller_lib import (
    apply_candidate_context,
    GENERIC_REMAINDER_IDS,
    aggregate_program_supported,
    aggregate_score,
    assignment_row,
    candidate_anchor,
    candidate_can_release,
    catalog_candidates,
    choose_group_parent,
    deterministic_candidate_membership_hash,
    dominant_generic_remainder_group,
    group_identity_core_fraction,
    hard_contradiction,
    group_family_support,
    index_scores,
    number,
    read_tsv,
    resolve_overlap,
    sha256,
    truth,
    validate_canonical_identity_component,
    validate_subset,
    validate_unmodeled_programs,
    write_manifest,
    write_tsv,
)


CANONICAL_CHAIN = (
    "run_observation_lineage_scoring.R",
    "derive_candidate_local_subsets.R",
    "close_exact_remainders.py",
    "run_lineage_controller.py",
)


def artifact_ok(record: dict, path: Path) -> bool:
    try:
        return (
            Path(str(record["path"])).resolve() == path.resolve()
            and str(record["sha256"]) == sha256(path)
        )
    except (KeyError, OSError):
        return False


def validate_stage_authority(
    authority_path: Path, contract_path: Path, expected_phase: str
) -> tuple[dict, dict]:
    authority = json.loads(authority_path.read_text(encoding="utf-8"))
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    errors: list[str] = []
    if authority.get("mode") != "stage_authority":
        errors.append("authority is not stage_authority")
    if authority.get("phase") != expected_phase:
        errors.append("authority does not permit the requested remainder phase")
    if authority.get("annotation_contract_sha256") != sha256(contract_path):
        errors.append("stage authority is stale for the annotation contract")
    controller = contract.get("canonical_lineage_controller", {})
    if controller.get("controller_version") != "2.2.0":
        errors.append("annotation contract does not bind lineage controller v2.2.0")
    scripts = controller.get("scripts", {})
    script_dir = Path(__file__).resolve().parent
    for name in CANONICAL_CHAIN:
        path = script_dir / name
        if not artifact_ok(scripts.get(name, {}), path):
            errors.append(f"formal chain does not bind the installed canonical {name}")
        if not artifact_ok(authority.get("scripts", {}).get(name, {}), path):
            errors.append(f"stage authority does not bind the installed canonical {name}")
    for dependency_name in (
        "controller_thresholds.py", "lineage_controller_lib.py"
    ):
        dependency = script_dir / dependency_name
        if not artifact_ok(
            controller.get("dependencies", {}).get(dependency_name, {}),
            dependency,
        ):
            errors.append(f"formal chain does not bind {dependency_name}")
        if not artifact_ok(
            authority.get("dependencies", {}).get(dependency_name, {}),
            dependency,
        ):
            errors.append(f"stage authority does not bind {dependency_name}")
    threshold_path = Path(str(contract.get("threshold_registry", {}).get("path", "")))
    if not threshold_path.is_absolute():
        threshold_path = (contract_path.parent / threshold_path).resolve()
    if not artifact_ok(contract.get("threshold_registry", {}), threshold_path):
        errors.append("formal chain does not bind the controller threshold registry")
    if not artifact_ok(authority.get("threshold_registry", {}), threshold_path):
        errors.append("stage authority does not bind the controller threshold registry")
    if errors:
        raise RuntimeError("; ".join(errors))
    return authority, contract


def load_proposals(
    membership_path: Path,
    evidence_path: Path,
) -> tuple[dict[str, list[str]], dict[str, dict[str, str]]]:
    memberships: dict[str, list[str]] = defaultdict(list)
    for row in read_tsv(membership_path):
        memberships[str(row["subset_id"])].append(str(row["cell_id"]))
    evidence = {str(row["subset_id"]): row for row in read_tsv(evidence_path)}
    if set(memberships) != set(evidence):
        raise RuntimeError("subset membership/evidence universes differ")
    return {key: sorted(set(value)) for key, value in memberships.items()}, evidence


def provenance_subset(
    proposals: list[tuple[str, str]], chosen_candidate: str
) -> str:
    matching = sorted(
        subset_id
        for subset_id, candidate_id in proposals
        if candidate_id == chosen_candidate
    )
    if matching:
        return matching[0]
    contributing = sorted(subset_id for subset_id, _ in proposals)
    if not contributing:
        raise ValueError("cannot trace an accepted assignment without proposals")
    return contributing[0]


def validated_residual_component_candidates(
    exact_remainder: list[str],
    validated_component_members: dict[str, set[str]],
    source_boundary: str,
    source_cluster: str,
    candidate_ids: list[str],
    score_index: dict[tuple[str, str], dict[str, str]],
    catalog: dict[str, dict],
    *,
    release_level: str = "broad",
    minimum_component_members: int = 5,
) -> tuple[set[str], set[str], list[dict[str, object]]]:
    """Audit whether a component's post-overlap tail remains separable.

    Membership in a previously validated expanded component is insufficient:
    only the exact residual intersection can preserve blocker authority.
    """
    exact_remainder_set = set(exact_remainder)
    separable: set[str] = set()
    ambiguous_members: set[str] = set()
    audit: list[dict[str, object]] = []
    for candidate, component_members in sorted(validated_component_members.items()):
        residual_members = sorted(exact_remainder_set.intersection(component_members))
        if len(residual_members) < minimum_component_members:
            continue
        validation = validate_subset(
            residual_members,
            candidate,
            score_index,
            candidate_ids,
            catalog=catalog,
            release_level=release_level,
        )
        audit.append({
            "subset_id": (
                f"{source_boundary}__{source_cluster}__{candidate}"
                "__post_overlap_separability"
            ),
            "source_boundary": source_boundary,
            "source_cluster": source_cluster,
            "candidate_id": candidate,
            "n_observations": len(residual_members),
            "round": "post_overlap_audit",
            **validation,
        })
        if validation["status"] == "PASS":
            separable.add(candidate)
            ambiguous_members.update(residual_members)
    return separable, ambiguous_members, audit


def remainder_candidate_programs(
    members: list[str],
    source_boundary: str,
    source_cluster: str,
    candidate_ids: list[str],
    score_index: dict[tuple[str, str], dict[str, str]],
    cluster_evidence: dict[tuple[str, str, str, str], dict[str, str]],
    catalog: dict[str, dict],
    residual_separable_candidate_ids: set[str] | None = None,
    *,
    writeback_policy: dict[str, float] | None = None,
) -> tuple[set[str], set[str], str, list[dict[str, object]]]:
    """Recompute the immutable-score evidence of one exact remainder.

    Broad/fine identity programs can block a coarse parent when they retain
    exact low-contradiction support or a validated local component. Aggregate-
    only shared/exploratory programs remain watches after both bounded local
    extraction rounds. State programs are recorded separately and do not
    replace or veto broad identity. Only a formal broad candidate with strong
    remainder-level support can itself become parent.
    """
    writeback_policy = writeback_policy or observation_writeback_defaults()
    residual_separable_candidate_ids = set(
        residual_separable_candidate_ids or set()
    ).intersection(candidate_ids)
    residual_specific_candidates = (
        residual_separable_candidate_ids - GENERIC_REMAINDER_IDS
    )
    records: list[dict[str, object]] = []
    for candidate_id in candidate_ids:
        rows = [score_index[(cell, candidate_id)] for cell in members]
        family = group_family_support(rows, catalog[candidate_id])
        supported_fraction = sum(candidate_anchor(row) for row in rows) / len(rows)
        contradiction_fraction = (
            sum(hard_contradiction(row) for row in rows) / len(rows)
        )
        supported_program_scores = [
            number(row.get("program_score"))
            for row in rows if candidate_anchor(row)
        ]
        mean_supported_program = (
            sum(supported_program_scores) / len(supported_program_scores)
            if supported_program_scores else 0.0
        )
        aggregate = cluster_evidence.get(
            ("selected", source_boundary, source_cluster, candidate_id)
        )
        aggregate_supported = bool(
            aggregate
            and aggregate_program_supported(
                aggregate, catalog[candidate_id], family
            )
        )
        exact_program_supported = (
            family["status"] == "PASS"
            and supported_fraction >= writeback_policy[
                "whole_subcluster_min_lineage_supported_fraction"
            ]
            and mean_supported_program >= 0.02
            and contradiction_fraction <= writeback_policy[
                "maximum_contradiction_fraction"
            ]
        )
        aggregate_generic_supported = bool(
            candidate_id in GENERIC_REMAINDER_IDS
            and aggregate
            and family["status"] == "PASS"
            and dominant_generic_remainder_group(
                aggregate, catalog[candidate_id]
            )
            # Specific identities have already received two bounded local
            # extraction attempts.  The exact tail may therefore contain
            # sparse/dropout observations whose generic parent program is
            # visible in fewer observations than in the source subcluster.
            # The dominant aggregate backbone, not a whole-object per-cell
            # classifier, certifies inheritance.  A still-separable specific
            # component blocks this route below.
            and supported_fraction >= writeback_policy[
                "whole_subcluster_min_lineage_supported_fraction"
            ]
            and mean_supported_program >= 0.02
            and not residual_specific_candidates
        )
        # After both bounded specific-lineage extractions, the exact tail can
        # itself become the evidentiary group.  A dominant multigene generic
        # program in that tail should return to its broad parent even when the
        # pre-split source subcluster was mixed and therefore failed the
        # whole-subcluster purity rule.  Candidate-specific residual
        # components remain the blocker; the generic candidate's anti signal
        # is not reused as a whole-tail veto.
        remainder_dominant_generic_supported = bool(
            candidate_id in GENERIC_REMAINDER_IDS
            and family["status"] == "PASS"
            and supported_fraction >= writeback_policy[
                "supported_subset_min_lineage_supported_fraction"
            ]
            and mean_supported_program >= 0.15
            and not residual_specific_candidates
        )
        generic_supported = bool(
            aggregate_generic_supported
            or remainder_dominant_generic_supported
        )
        candidate = catalog[candidate_id]
        candidate_role = str(candidate.get("candidate_role", ""))
        residual_separable = candidate_id in residual_separable_candidate_ids
        credible = bool(
            candidate_role != "state"
            and (
                generic_supported
                or residual_separable
                or (
                    exact_program_supported
                    and candidate_role in {"broad", "fine"}
                )
            )
        )
        formal_broad = (
            candidate_can_release(candidate)
            and str(candidate.get("candidate_role", "")) == "broad"
            and not str(candidate.get("release_fine_label", ""))
        )
        parent_eligible = bool(
            credible
            and formal_broad
            and family["status"] == "PASS"
            and (
                (
                    generic_supported
                    and max(
                        supported_fraction,
                        group_identity_core_fraction(aggregate),
                    ) >= writeback_policy[
                        "supported_subset_min_lineage_supported_fraction"
                    ]
                )
                or (
                    not generic_supported
                    and supported_fraction >= writeback_policy[
                        "supported_subset_min_lineage_supported_fraction"
                    ]
                    and contradiction_fraction <= writeback_policy[
                        "maximum_contradiction_fraction"
                    ]
                )
            )
            and (
                generic_supported
                or aggregate_supported
                or mean_supported_program >= 0.04
            )
        )
        records.append({
            "candidate_id": candidate_id,
            "candidate_role": candidate_role,
            "release_broad_label": str(candidate.get("release_broad_label", "")),
            "supported_fraction": supported_fraction,
            "effective_parent_supported_fraction": (
                max(
                    supported_fraction,
                    group_identity_core_fraction(aggregate),
                )
                if generic_supported and aggregate
                else supported_fraction
            ),
            "contradiction_fraction": contradiction_fraction,
            "mean_supported_program": mean_supported_program,
            "family_support_status": family["status"],
            "aggregate_supported": aggregate_supported,
            "aggregate_only_watch": bool(
                aggregate_supported
                and not exact_program_supported
                and not residual_separable
                and not generic_supported
            ),
            "aggregate_score": aggregate_score(aggregate) if aggregate else float("-inf"),
            "generic_remainder_supported": generic_supported,
            "aggregate_generic_remainder_supported": aggregate_generic_supported,
            "remainder_dominant_generic_supported": (
                remainder_dominant_generic_supported
            ),
            "residual_separable_component": residual_separable,
            "credible_blocker": credible,
            "parent_eligible": parent_eligible,
        })
    credible_ids = {
        str(row["candidate_id"]) for row in records if row["credible_blocker"]
    }
    parent_ids = {
        str(row["candidate_id"]) for row in records if row["parent_eligible"]
    }
    preferred = ""
    if parent_ids:
        ranked = sorted(
            (row for row in records if row["candidate_id"] in parent_ids),
            key=lambda row: (
                -number(row["effective_parent_supported_fraction"]),
                -number(row["aggregate_score"], -1e9),
                str(row["candidate_id"]),
            ),
        )
        winner = ranked[0]
        winner_broad = str(winner["release_broad_label"])
        runner_fraction = max(
            (
                number(row["supported_fraction"])
                for row in records
                if row["credible_blocker"]
                and str(row["release_broad_label"]) != winner_broad
            ),
            default=0.0,
        )
        if (
            number(winner["effective_parent_supported_fraction"])
            >= writeback_policy[
                "supported_subset_min_lineage_supported_fraction"
            ]
            and number(winner["effective_parent_supported_fraction"])
            - runner_fraction >= writeback_policy[
                "supported_subset_min_purity_margin"
            ]
        ):
            preferred = str(winner["candidate_id"])
    return credible_ids, parent_ids, preferred, records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--cluster-evidence", required=True, type=Path)
    parser.add_argument("--subset-membership", required=True, type=Path)
    parser.add_argument("--subset-evidence", required=True, type=Path)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--stage-authority", required=True, type=Path)
    parser.add_argument(
        "--scope", required=True,
        choices=["local_mixed_subcluster", "post_merge_reconciliation"],
    )
    parser.add_argument("--source-boundary")
    parser.add_argument("--source-cluster")
    parser.add_argument("--unmodeled-programs", type=Path)
    parser.add_argument("--context-evidence", type=Path)
    parser.add_argument("--release-level", choices=["broad"], default="broad")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    expected_phase = (
        "local_mixed_subcluster_split"
        if args.scope == "local_mixed_subcluster"
        else "merge_and_freeze_broad"
    )
    if args.scope == "local_mixed_subcluster" and (
        not args.source_boundary or not args.source_cluster
    ):
        raise RuntimeError("local remainder closure requires one source subcluster")
    authority, contract = validate_stage_authority(
        args.stage_authority.resolve(), args.contract.resolve(), expected_phase
    )
    threshold_registry_path = Path(str(contract["threshold_registry"]["path"]))
    if not threshold_registry_path.is_absolute():
        threshold_registry_path = (
            args.contract.parent / threshold_registry_path
        ).resolve()
    threshold_registry = load_controller_thresholds(threshold_registry_path)
    writeback_policy = contract.get("observation_writeback", {}).get("policy", {})
    required_writeback = set(
        threshold_registry["observation_writeback_policy"]
    )
    if not required_writeback.issubset(writeback_policy):
        raise RuntimeError("contract lacks the complete observation-writeback policy")
    local_subset_policy = threshold_registry["local_subset_policy"]
    minimum_component_members = int(
        local_subset_policy["minimum_component_members"]
    )
    maximum_second_subset_rounds = int(
        local_subset_policy["maximum_second_subset_rounds"]
    )
    catalog_document = json.loads(args.catalog.read_text(encoding="utf-8"))
    catalog = catalog_candidates(catalog_document)
    bound_context = authority.get("context_evidence")
    if args.context_evidence:
        if (
            not bound_context
            or Path(str(bound_context.get("path", ""))).resolve()
            != args.context_evidence.resolve()
            or bound_context.get("sha256") != sha256(args.context_evidence)
        ):
            raise RuntimeError("context evidence differs from stage authority")
    elif bound_context:
        raise RuntimeError("stage authority binds context evidence but none was supplied")
    apply_candidate_context(
        catalog,
        read_tsv(args.context_evidence) if args.context_evidence else [],
    )
    # Freeze the broad phase first.  Fine labels are materialized only after
    # every exact remainder has closed and only inside the candidate's frozen
    # broad parent.  The fine candidate_id remains provenance in broad phase.
    release_catalog = {
        candidate_id: {**candidate, "release_fine_label": ""}
        for candidate_id, candidate in catalog.items()
    }
    if contract.get("candidate_catalog", {}).get("sha256") != sha256(args.catalog):
        raise RuntimeError("candidate catalog differs from the annotation contract")

    score_rows = read_tsv(args.scores)
    if args.scope == "local_mixed_subcluster":
        score_rows = [
            row for row in score_rows
            if row.get("source_boundary") == args.source_boundary
            and row.get("source_cluster") == args.source_cluster
        ]
        if not score_rows:
            raise RuntimeError("requested local mixed subcluster is absent from scores")
    score_index, groups, candidate_ids = index_scores(score_rows)
    if set(candidate_ids) != set(catalog):
        raise RuntimeError("score candidate universe differs from the bound catalog")
    memberships, proposal_evidence = load_proposals(
        args.subset_membership, args.subset_evidence
    )
    cluster_evidence_rows = read_tsv(args.cluster_evidence)
    if args.scope == "local_mixed_subcluster":
        cluster_evidence_rows = [
            row for row in cluster_evidence_rows
            if row.get("source_boundary") == args.source_boundary
            and row.get("source_cluster") == args.source_cluster
        ]
    cluster_evidence = {
        (
            row.get("resolution_role", ""),
            row.get("source_boundary", ""),
            row.get("source_cluster", ""),
            row.get("candidate_id", ""),
        ): row
        for row in cluster_evidence_rows
    }
    bound_cluster_evidence = authority.get("cluster_evidence", {})
    if not artifact_ok(bound_cluster_evidence, args.cluster_evidence):
        raise RuntimeError("cluster evidence differs from stage authority")

    accepted_by_cell: dict[str, tuple[str, str, str]] = {}
    subset_audit: list[dict[str, object]] = []
    proposals_by_cell: dict[str, list[tuple[str, str]]] = defaultdict(list)
    validated_component_members: dict[str, set[str]] = defaultdict(set)
    for subset_id, members in memberships.items():
        evidence = proposal_evidence[subset_id]
        candidate = str(evidence.get("candidate_id", ""))
        source_boundary = str(evidence.get("source_boundary", ""))
        source_cluster = str(evidence.get("source_cluster", ""))
        source_members = set(groups.get((source_boundary, source_cluster), []))
        if str(evidence.get("status", "")).upper() != "PASS":
            validation = {
                "status": "FAIL",
                "reason": "preliminary_subset_validation_failed",
            }
        elif candidate not in catalog or not candidate_can_release(catalog[candidate]):
            validation = {"status": "FAIL", "reason": "candidate_not_release_eligible"}
        elif not set(members).issubset(source_members):
            validation = {"status": "FAIL", "reason": "membership_outside_source"}
        else:
            proposal_scope = str(evidence.get("proposal_scope", ""))
            aggregate_evidence = None
            if proposal_scope == "whole_subcluster":
                aggregate_evidence = cluster_evidence.get(
                    ("selected", source_boundary, source_cluster, candidate)
                )
            elif proposal_scope == "neighboring_resolution_expression_subcluster":
                aggregate_evidence = cluster_evidence.get(
                    (
                        str(evidence.get("source_resolution_role", "")),
                        source_boundary,
                        str(evidence.get("source_resolution_cluster", "")),
                        candidate,
                    )
                )
            if proposal_scope == "canonical_identity_component":
                aggregate_evidence = cluster_evidence.get(
                    ("selected", source_boundary, source_cluster, candidate)
                )
                validation = validate_canonical_identity_component(
                    members,
                    sorted(source_members),
                    candidate,
                    score_index,
                    catalog[candidate],
                    aggregate_evidence,
                )
            elif (
                proposal_scope in {
                    "whole_subcluster",
                    "neighboring_resolution_expression_subcluster",
                }
                and aggregate_evidence is None
            ):
                validation = {
                    "status": "FAIL",
                    "reason": "aggregate_multichannel_evidence_missing",
                }
            else:
                validation = validate_subset(
                    members, candidate, score_index, candidate_ids,
                    catalog=catalog, release_level=args.release_level,
                    aggregate_evidence=aggregate_evidence,
                )
        subset_audit.append({
            "subset_id": subset_id,
            "source_boundary": source_boundary,
            "source_cluster": source_cluster,
            "candidate_id": candidate,
            "n_observations": len(members),
            **validation,
        })
        if validation["status"] == "PASS":
            validated_component_members[candidate].update(members)
            for cell in members:
                proposals_by_cell[cell].append((subset_id, candidate))

    # Resolve overlapping accepted proposals from scores, never proposal order.
    for cell in sorted(proposals_by_cell):
        proposals = proposals_by_cell[cell]
        candidates = sorted({candidate for _, candidate in proposals})
        if len(candidates) == 1:
            chosen = candidates[0]
        else:
            chosen, _ = resolve_overlap(
                cell, candidates, score_index, catalog,
                release_level=args.release_level,
            )
        if chosen:
            # Fine-candidate overlap may deliberately collapse to their
            # supported common broad parent. That parent is not itself one of
            # the proposal tuples, so retain the first contributing validated
            # subset as provenance instead of indexing an empty match.
            subset = provenance_subset(proposals, chosen)
            proposal_scope = str(
                proposal_evidence[subset].get("proposal_scope", "")
            )
            accepted_by_cell[cell] = (
                chosen,
                (
                    "whole_subcluster"
                    if proposal_scope == "whole_subcluster"
                    else "supported_subset"
                ),
                subset,
            )

    assignments: list[dict[str, object]] = []
    remainder_audit: list[dict[str, object]] = []
    for (source_boundary, source_cluster), source_members in sorted(groups.items()):
        assigned = {
            cell for cell in source_members if cell in accepted_by_cell
        }
        for cell in sorted(assigned):
            candidate, origin, _ = accepted_by_cell[cell]
            assignments.append(assignment_row(
                cell, source_boundary, source_cluster, candidate, release_catalog[candidate],
                origin, "high",
            ))

        remainder = sorted(set(source_members) - assigned)
        round_two_assigned: dict[str, str] = {}
        round_two_by_cell: dict[str, list[str]] = defaultdict(list)
        if remainder and maximum_second_subset_rounds > 0:
            # The final local extraction may only reconsider precomputed
            # candidate-local spatial components after accepted memberships
            # have been removed. It cannot invent a whole-object per-cell
            # classifier from the remainder.
            remainder_set = set(remainder)
            residual_components = []
            for subset_id, proposed_members in memberships.items():
                evidence = proposal_evidence[subset_id]
                if str(evidence.get("status", "")).upper() == "PASS":
                    continue
                if str(evidence.get("proposal_scope", "")) not in {
                    "candidate_local_spatial_component",
                    "canonical_identity_component",
                }:
                    continue
                residual_members = sorted(set(proposed_members) & remainder_set)
                if len(residual_members) < minimum_component_members:
                    continue
                residual_components.append(
                    (subset_id, str(evidence.get("candidate_id", "")), residual_members)
                )
            residual_components.sort(
                key=lambda item: (
                    item[1] == "stromal_mesenchymal",
                    -int(catalog.get(item[1], {}).get("specificity_priority", 0)),
                    item[1],
                    item[0],
                )
            )
            for subset_id, candidate, residual_members in residual_components:
                if candidate not in catalog or not candidate_can_release(catalog[candidate]):
                    continue
                validation = validate_subset(
                    residual_members, candidate, score_index, candidate_ids,
                    catalog=catalog, release_level=args.release_level,
                )
                subset_audit.append({
                    "subset_id": (
                        f"{subset_id}__remainder_round2"
                    ),
                    "source_boundary": source_boundary,
                    "source_cluster": source_cluster,
                    "candidate_id": candidate,
                    "n_observations": len(residual_members),
                    "round": 2,
                    **validation,
                })
                if validation["status"] == "PASS":
                    validated_component_members[candidate].update(residual_members)
                    for cell in residual_members:
                        round_two_by_cell[cell].append(candidate)
            for cell, candidates in sorted(round_two_by_cell.items()):
                if len(candidates) == 1:
                    chosen = candidates[0]
                else:
                    chosen, _ = resolve_overlap(
                        cell, candidates, score_index, catalog,
                        release_level=args.release_level,
                    )
                if chosen:
                    round_two_assigned[cell] = chosen
            for cell, candidate in sorted(round_two_assigned.items()):
                assignments.append(assignment_row(
                    cell, source_boundary, source_cluster, candidate,
                    release_catalog[candidate], "supported_subset_round2", "high",
                ))

        exact_remainder = sorted(set(remainder) - set(round_two_assigned))
        if exact_remainder:
            # A validated component can leave only mutually ambiguous or
            # low-information expansion members in the exact remainder.  The
            # historical component alone is not proof that this residual tail
            # remains separable.  Revalidate each candidate's exact residual
            # intersection; only a still release-grade subset may block the
            # generic parent.  This is an ambiguity audit, not a third subset
            # extraction or a whole-object per-cell classifier.
            (
                residual_separable_candidates,
                ambiguous_component_members,
                post_overlap_audit,
            ) = validated_residual_component_candidates(
                exact_remainder,
                validated_component_members,
                source_boundary,
                source_cluster,
                candidate_ids,
                score_index,
                catalog,
                release_level=args.release_level,
                minimum_component_members=minimum_component_members,
            )
            subset_audit.extend(post_overlap_audit)
            if ambiguous_component_members:
                remainder_audit.append({
                    "source_boundary": source_boundary,
                    "source_cluster": source_cluster,
                    "closure_partition": "post_overlap_ambiguous_core",
                    "n_source_exact_remainder": len(exact_remainder),
                    "n_exact_remainder": len(ambiguous_component_members),
                    "selected_parent_candidate": "",
                    "credible_blocker_candidates": ";".join(
                        sorted(residual_separable_candidates)
                    ),
                    "residual_separable_candidates": ";".join(
                        sorted(residual_separable_candidates)
                    ),
                    "eligible_parent_candidates": "",
                    "preferred_parent_candidate": "",
                    "candidate_program_audit_json": "[]",
                    "status": "UNRESOLVED",
                    "reason": "residual_component_overlap_after_two_extractions",
                })
                for cell in sorted(ambiguous_component_members):
                    assignments.append(assignment_row(
                        cell,
                        source_boundary,
                        source_cluster,
                        "",
                        None,
                        "post_overlap_ambiguous_remainder",
                        "low",
                        qc_reason="irreducible_lineage_overlap",
                    ))

            parent_remainder = sorted(
                set(exact_remainder) - ambiguous_component_members
            )
            if parent_remainder:
                (
                    credible_blockers,
                    parent_candidates,
                    preferred,
                    remainder_programs,
                ) = remainder_candidate_programs(
                    parent_remainder,
                    source_boundary,
                    source_cluster,
                    candidate_ids,
                    score_index,
                    cluster_evidence,
                    catalog,
                    residual_separable_candidates,
                    writeback_policy=writeback_policy,
                )
                parent, parent_evidence = choose_group_parent(
                    parent_remainder,
                    candidate_ids,
                    score_index,
                    catalog,
                    preferred_parent=preferred,
                    release_level=args.release_level,
                    blocker_candidate_ids=credible_blockers,
                    parent_candidate_ids=parent_candidates,
                    certified_generic_parent_ids=(
                        parent_candidates.intersection(GENERIC_REMAINDER_IDS)
                    ),
                )
                remainder_audit.append({
                    "source_boundary": source_boundary,
                    "source_cluster": source_cluster,
                    "closure_partition": "post_overlap_parent_tail",
                    "n_source_exact_remainder": len(exact_remainder),
                    "n_exact_remainder": len(parent_remainder),
                    "selected_parent_candidate": parent,
                    "credible_blocker_candidates": ";".join(
                        sorted(credible_blockers)
                    ),
                    "residual_separable_candidates": "",
                    "eligible_parent_candidates": ";".join(
                        sorted(parent_candidates)
                    ),
                    "preferred_parent_candidate": preferred,
                    "candidate_program_audit_json": json.dumps(
                        remainder_programs, ensure_ascii=False, sort_keys=True
                    ),
                    **parent_evidence,
                })
                if parent:
                    for cell in parent_remainder:
                        assignments.append(assignment_row(
                            cell,
                            source_boundary,
                            source_cluster,
                            parent,
                            release_catalog[parent],
                            "exact_remainder_parent",
                            "moderate",
                        ))
                else:
                    for cell in parent_remainder:
                        rows = [
                            score_index[(cell, candidate)]
                            for candidate in candidate_ids
                        ]
                        positive = any(
                            int(number(row.get("positive_family_count"), 0)) > 0
                            for row in rows
                        )
                        technical = any(
                            truth(row.get("technical_flag")) for row in rows
                        )
                        reason = (
                            "technical_or_low_rna" if technical
                            else "irreducible_lineage_overlap" if positive
                            else "low_information"
                        )
                        assignments.append(assignment_row(
                            cell,
                            source_boundary,
                            source_cluster,
                            "",
                            None,
                            "exact_remainder_qc",
                            "low",
                            qc_reason=reason,
                        ))

    unmodeled_rows = (
        validate_unmodeled_programs(read_tsv(args.unmodeled_programs))
        if args.unmodeled_programs and args.unmodeled_programs.is_file()
        else []
    )
    assignment_ids = [str(row["cell_id"]) for row in assignments]
    if len(assignment_ids) != len(set(assignment_ids)):
        raise RuntimeError("release membership is not mutually exclusive")
    expected_ids = {cell for members in groups.values() for cell in members}
    if set(assignment_ids) != expected_ids:
        raise RuntimeError("release membership does not exactly close the source universe")

    assignments = sorted(assignments, key=lambda row: str(row["cell_id"]))
    proposal_rows: list[dict[str, object]] = []
    for row in assignments:
        broad = str(row.get("final_broad_label", ""))
        unresolved = str(row.get("qc_reason", ""))
        proposal_rows.append({
            "cell_id": row.get("cell_id", ""),
            "source_boundary": row.get("source_boundary", ""),
            "source_cluster": row.get("source_cluster", ""),
            "candidate_id": row.get("candidate_id", ""),
            "proposed_state": "broad_candidate" if broad else "unresolved_biological",
            "proposed_broad_label": broad,
            "confidence": row.get("confidence", ""),
            "assignment_origin": (
                "exact_remainder_unresolved"
                if row.get("assignment_origin") == "exact_remainder_qc"
                else row.get("assignment_origin", "")
            ),
            "unresolved_reason": unresolved,
        })
    membership_hash = deterministic_candidate_membership_hash(proposal_rows)
    args.out.mkdir(parents=True, exist_ok=True)
    membership_path = args.out / "candidate_membership.tsv.gz"
    write_tsv(membership_path, proposal_rows)
    write_tsv(args.out / "subset_validation.tsv", subset_audit)
    write_tsv(args.out / "exact_remainder_audit.tsv", remainder_audit)
    write_tsv(
        args.out / "unmodeled_lineage_candidates.tsv",
        [
            {
                **row,
                "release_status": "Unmodeled lineage candidate",
                "formal_label_written": "false",
            }
            for row in unmodeled_rows
        ],
        fields=(
            list(unmodeled_rows[0]) + ["release_status", "formal_label_written"]
            if unmodeled_rows else
            ["program_id", "release_status", "formal_label_written"]
        ),
    )
    unresolved_n = sum(
        str(row["proposed_state"]) == "unresolved_biological"
        for row in proposal_rows
    )
    manifest = {
        "status": "PASS",
        "schema_version": "2.2",
        "controller_version": "2.2.0",
        "stage": expected_phase,
        "scope": args.scope,
        "formal_membership_written": False,
        "release_level": "broad_candidate_only",
        "annotation_contract": {
            "path": str(args.contract.resolve()),
            "sha256": sha256(args.contract),
        },
        "stage_authority": {
            "path": str(args.stage_authority.resolve()),
            "sha256": sha256(args.stage_authority),
        },
        "scores": {"path": str(args.scores.resolve()), "sha256": sha256(args.scores)},
        "cluster_evidence": {
            "path": str(args.cluster_evidence.resolve()),
            "sha256": sha256(args.cluster_evidence),
        },
        "subset_membership": {
            "path": str(args.subset_membership.resolve()),
            "sha256": sha256(args.subset_membership),
        },
        "candidate_membership": {
            "path": str(membership_path.resolve()),
            "sha256": sha256(membership_path),
            "semantic_sha256": membership_hash,
            "n_observations": len(proposal_rows),
        },
        "unresolved_biological_n": unresolved_n,
        "unresolved_biological_fraction": (
            unresolved_n / len(proposal_rows) if proposal_rows else 0.0
        ),
        "second_subset_round_limit": maximum_second_subset_rounds,
        "unmodeled_lineage_candidate_n": len(unmodeled_rows),
        "context_evidence": (
            {
                "path": str(args.context_evidence.resolve()),
                "sha256": sha256(args.context_evidence),
            }
            if args.context_evidence else None
        ),
        "immutable_observation_scores": True,
        "historical_labels_visible_during_scoring": False,
    }
    write_manifest(args.out / "exact_remainder_closure_manifest.json", manifest)
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
