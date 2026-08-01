#!/usr/bin/env python3
"""Audit every context-evaluable broad lineage as one biological review unit.

This is a post-Atlas reviewer, not a whole-object classifier.  One annotated
broad cell type produces one precision + recall + spatial review question.
Second-round source subclusters, direct-seeded spatial components and group
watches are internal evidence and bounded patch regions; they are not separate
user-facing biological decisions.  The audit never writes labels.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from controller_thresholds import load_controller_thresholds
from evidence_schema_lib import sha256, validate_json_against_schema
from lineage_controller_lib import (
    GENERIC_REMAINDER_IDS,
    apply_candidate_context,
    candidate_can_release,
    candidate_can_support_broad_review,
    catalog_candidates,
    deterministic_cell_id_set_hash,
    deterministic_membership_hash,
    group_candidate_detected,
    group_identity_core_direct_fraction,
    group_identity_core_fraction,
    independent_group_program,
    number,
    read_tsv,
)


REQUIRED_SCORE_COLUMNS = {
    "cell_id", "source_boundary", "source_cluster", "candidate_id",
    "normalized_evidence", "direct_signal", "program_score",
    "positive_family_count", "family_coherent", "identity_core_coherent",
    "identity_core_direct", "hard_contradiction", "technical_flag", "x", "y",
}

PRECISION_COLUMNS = [
    "broad_label", "source_boundary", "source_cluster", "n_current_label",
    "source_group_n", "current_label_fraction", "group_program_detected",
    "group_program_candidate_ids", "observation_supported_n",
    "effective_supported_n", "effective_supported_fraction",
    "sparse_inheritance_used", "strongest_competing_broad",
    "strongest_competing_direct_n", "strongest_competing_direct_fraction",
    "status", "reason", "unit_signature",
]
COMPONENT_COLUMNS = [
    "component_id", "review_round", "target_broad_label", "source_boundary",
    "source_cluster", "n_observations", "n_direct_seeds",
    "direct_seed_fraction", "median_normalized_evidence",
    "median_direct_signal", "median_broad_evidence_margin",
    "group_program_detected", "group_program_candidate_ids",
    "candidate_id_census", "current_broad_census", "status", "reason",
    "unit_signature",
]
GROUP_WATCH_COLUMNS = [
    "review_round", "target_broad_label", "source_boundary",
    "source_cluster", "source_group_n", "current_target_n",
    "current_target_fraction", "expected_identity_core_fraction",
    "expected_identity_core_direct_fraction", "expected_missing_n",
    "identity_fraction_deficit", "strongest_candidate_id",
    "group_program_candidate_ids", "status", "reason", "unit_signature",
]


def artifact(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "n_bytes": path.stat().st_size,
    }


def unit_signature(mode: str, broad: str, cells: list[str]) -> str:
    payload = "\n".join([mode, broad, *sorted(cells)]).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def prior_closed_units(paths: list[Path]) -> set[tuple[str, str, str]]:
    closed: set[tuple[str, str, str]] = set()
    for path in paths:
        validation = json.loads(path.read_text(encoding="utf-8"))
        if validation.get("status") != "PASS":
            raise SystemExit("prior catalog-wide decision validation is not PASS")
        review_path = Path(str(validation.get("review_manifest", {}).get("path", "")))
        decision_path = Path(str(validation.get("validated_decisions", {}).get("path", "")))
        if (
            not review_path.is_file()
            or validation.get("review_manifest", {}).get("sha256") != sha256(review_path)
            or not decision_path.is_file()
            or validation.get("validated_decisions", {}).get("sha256") != sha256(decision_path)
        ):
            raise SystemExit("prior catalog-wide decision validation is stale")
        review = json.loads(review_path.read_text(encoding="utf-8"))
        queue_path = Path(str(review.get("artifacts", {}).get("review_queue", {}).get("path", "")))
        if not queue_path.is_file() or review.get("artifacts", {}).get("review_queue", {}).get("sha256") != sha256(queue_path):
            raise SystemExit("prior catalog-wide review queue is stale")
        queue = {str(row.get("review_id", "")): row for row in read_tsv(queue_path)}
        for decision in read_tsv(decision_path):
            outcome = str(decision.get("outcome", ""))
            if outcome not in {
                "retain_current_label", "reject_shared_or_ambient",
                "retain_current_parent",
                "retain_current_cell_type", "confirm_absent_or_not_evaluable",
            }:
                continue
            queued = queue.get(str(decision.get("review_id", "")), {})
            key = (
                str(queued.get("review_mode", "")),
                str(queued.get("target_broad_label", "")),
                str(queued.get("unit_signature", "")),
            )
            if all(key):
                closed.add(key)
    return closed


def prior_broad_scopes(paths: list[Path]) -> dict[str, dict[str, object]]:
    """Recover the exact member/recall/watch scope of a closed broad review."""
    result: dict[str, dict[str, object]] = {}
    for path in paths:
        validation = json.loads(path.read_text(encoding="utf-8"))
        if validation.get("status") != "PASS":
            continue
        review_path = Path(str(validation.get("review_manifest", {}).get("path", "")))
        decisions_path = Path(str(validation.get("validated_decisions", {}).get("path", "")))
        if not review_path.is_file() or not decisions_path.is_file():
            continue
        review = json.loads(review_path.read_text(encoding="utf-8"))
        queue_record = review.get("artifacts", {}).get("review_queue", {})
        scope_record = review.get("artifacts", {}).get(
            "broad_lineage_review_scope_membership", {}
        )
        queue_path = Path(str(queue_record.get("path", "")))
        scope_path = Path(str(scope_record.get("path", "")))
        if (
            not queue_path.is_file()
            or queue_record.get("sha256") != sha256(queue_path)
            or not scope_path.is_file()
            or scope_record.get("sha256") != sha256(scope_path)
        ):
            continue
        queue = {str(row.get("review_id", "")): row for row in read_tsv(queue_path)}
        scope_by_review: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in read_tsv(scope_path):
            scope_by_review[str(row.get("review_id", ""))].append(row)
        for decision in read_tsv(decisions_path):
            if str(decision.get("outcome", "")) not in {
                "retain_current_label", "reject_shared_or_ambient",
                "retain_current_parent", "retain_current_cell_type",
                "confirm_absent_or_not_evaluable",
            }:
                continue
            queued = queue.get(str(decision.get("review_id", "")), {})
            broad = str(queued.get("target_broad_label", ""))
            if str(queued.get("review_mode", "")) != "broad_lineage_review" or not broad:
                continue
            rows = scope_by_review.get(str(decision.get("review_id", "")), [])
            role_ids: dict[str, set[str]] = defaultdict(set)
            for row in rows:
                cell = str(row.get("cell_id", ""))
                for role in str(row.get("scope_roles", "")).split(";"):
                    if role and cell:
                        role_ids[role].add(cell)
            result[broad] = {
                "unit_signature": str(queued.get("unit_signature", "")),
                "current_ids": role_ids["current_label"],
                "component_ids": role_ids["direct_recall_component"],
                "watch_ids": role_ids["group_watch_source"],
                "zero_direct_ids": role_ids[
                    "zero_census_direct_multifamily_challenger"
                ],
            }
    return result


def manual_closed_units(
    paths: list[Path], membership: pd.DataFrame,
) -> tuple[set[tuple[str, str, str]], list[dict[str, object]]]:
    """Bind non-mutating user adjudication to the exact current target set."""
    closed: set[tuple[str, str, str]] = set()
    records: list[dict[str, object]] = []
    current_rows = membership.to_dict("records")
    for path in paths:
        document = json.loads(path.read_text(encoding="utf-8"))
        if (
            document.get("status") != "PASS"
            or document.get("artifact_role")
            != "user_authorized_manual_biological_adjudication"
            or document.get("membership_changed") is not False
            or document.get("counts_as_automatic_decision_round") is not False
            or document.get("outcome") not in {
                "retain_current_cell_type", "confirm_absent_or_not_evaluable"
            }
        ):
            raise SystemExit("manual adjudication is not a canonical closure")
        broad = str(document.get("target_broad_label", ""))
        key = (
            str(document.get("review_mode", "")), broad,
            str(document.get("unit_signature", "")),
        )
        adjudicated_path = Path(str(document.get("membership", {}).get("path", "")))
        if (
            not all(key)
            or not adjudicated_path.is_file()
            or document.get("membership", {}).get("sha256") != sha256(adjudicated_path)
        ):
            raise SystemExit("manual adjudication membership or exact scope is stale")
        current_target = [
            row for row in current_rows
            if str(row.get("final_broad_label", row.get("broad_label", ""))) == broad
        ]
        adjudicated_target = [
            row for row in read_tsv(adjudicated_path)
            if str(row.get("final_broad_label", row.get("broad_label", ""))) == broad
        ]
        if deterministic_membership_hash(current_target) != deterministic_membership_hash(
            adjudicated_target
        ):
            raise SystemExit("manual adjudication target membership differs from current review")
        if key in closed:
            raise SystemExit("duplicate manual adjudication for one exact scope")
        closed.add(key)
        records.append(artifact(path))
    return closed, records


def load_zero_census_challengers(
    path: Path | None, membership: pd.DataFrame,
) -> tuple[dict[str, set[str]], dict[str, dict[str, str]], dict | None]:
    if path is None:
        return {}, {}, None
    document = json.loads(path.read_text(encoding="utf-8"))
    if (
        document.get("status") != "PASS"
        or document.get("artifact_role")
        != "query_derived_zero_census_challenger"
        or document.get("formal_membership_written") is not False
        or document.get("membership_cell_id_set_sha256")
        != deterministic_cell_id_set_hash(membership.to_dict("records"))
    ):
        raise SystemExit("zero-census direct challenger manifest is incompatible")
    artifacts = document.get("artifacts", {})
    summary_path = Path(str(artifacts.get("summary", {}).get("path", "")))
    members_path = Path(str(artifacts.get("membership", {}).get("path", "")))
    for label, source in (("summary", summary_path), ("membership", members_path)):
        if not source.is_file() or artifacts.get(label, {}).get("sha256") != sha256(source):
            raise SystemExit(f"zero-census challenger {label} is missing or stale")
    summary = {
        str(row.get("broad_label", "")): row for row in read_tsv(summary_path)
    }
    direct: dict[str, set[str]] = defaultdict(set)
    eligible = {
        broad for broad, row in summary.items()
        if str(row.get("status", "")) == "review_required"
    }
    for row in read_tsv(members_path):
        broad = str(row.get("broad_label", ""))
        if broad in eligible:
            direct[broad].add(str(row.get("cell_id", "")))
    return dict(direct), summary, artifact(path)


def validate_authority(args: argparse.Namespace) -> dict:
    authority = json.loads(args.stage_authority.read_text(encoding="utf-8"))
    if (
        authority.get("mode") != "stage_authority"
        or authority.get("phase") != "atlas_and_completeness_review"
        or authority.get("annotation_contract_sha256") != sha256(args.contract)
    ):
        raise SystemExit("stage authority does not permit catalog-wide review")
    records = {
        "post_atlas_membership": args.membership,
        "candidate_catalog": args.catalog,
        "threshold_registry": args.threshold_registry,
    }
    if args.context_evidence:
        records["context_evidence"] = args.context_evidence
    if args.zero_census_direct_challenger_manifest:
        records["zero_census_direct_challenger"] = (
            args.zero_census_direct_challenger_manifest
        )
    for key, path in records.items():
        record = authority.get(key, {})
        if (
            Path(str(record.get("path", ""))).resolve() != path.resolve()
            or record.get("sha256") != sha256(path)
        ):
            raise SystemExit(f"catalog-wide review authority differs for {key}")
    for key, supplied in (
        ("observation_scores", args.scores),
        ("cluster_evidence", args.cluster_evidence),
    ):
        bound = {
            (Path(str(row.get("path", ""))).resolve(), str(row.get("sha256", "")))
            for row in authority.get(key, [])
        }
        observed = {(path.resolve(), sha256(path)) for path in supplied}
        if bound != observed:
            raise SystemExit(f"catalog-wide review authority differs for {key}")
    return authority


def boolean(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(
        {"1", "true", "yes", "pass", "passed", "t"}
    )


def load_scores(
    paths: list[Path], candidates: dict[str, dict], workers: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    allowed = set(candidates)

    def read_one(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
        header = set(pd.read_csv(path, sep="\t", nrows=0).columns)
        if missing := REQUIRED_SCORE_COLUMNS.difference(header):
            raise SystemExit(f"observation score lacks {sorted(missing)}: {path}")
        usecols = sorted(REQUIRED_SCORE_COLUMNS | ({"release_family_coherent"} & header))
        frame = pd.read_csv(
            path, sep="\t", usecols=usecols, dtype={"cell_id": str},
            low_memory=False,
        )
        coordinates = frame[["cell_id", "x", "y"]].drop_duplicates("cell_id")
        frame = frame.loc[frame.candidate_id.astype(str).isin(allowed)].copy()
        return frame, coordinates

    maximum_workers = max(1, min(int(workers), len(paths)))
    if maximum_workers == 1:
        loaded = [read_one(path) for path in paths]
    else:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=maximum_workers
        ) as executor:
            loaded = list(executor.map(read_one, paths))
    frames = [item[0] for item in loaded]
    coordinate_frames = [item[1] for item in loaded]
    if not frames:
        raise SystemExit("catalog-wide review requires observation scores")
    scores = pd.concat(frames, ignore_index=True)
    scores["cell_id"] = scores.cell_id.astype(str)
    scores["candidate_id"] = scores.candidate_id.astype(str)
    if scores.duplicated(["cell_id", "candidate_id"]).any():
        raise SystemExit("observation scores duplicate cell_id x candidate_id")
    coords = pd.concat(coordinate_frames, ignore_index=True)
    coords["cell_id"] = coords.cell_id.astype(str)
    coordinate_conflict = coords.groupby("cell_id")[["x", "y"]].nunique().max(axis=1).gt(1)
    if coordinate_conflict.any():
        raise SystemExit("observation score files disagree on cell coordinates")
    coords = coords.drop_duplicates("cell_id")
    for column in (
        "normalized_evidence", "direct_signal", "program_score",
        "positive_family_count", "x", "y",
    ):
        scores[column] = pd.to_numeric(scores[column], errors="coerce").fillna(0.0)
    for column in (
        "family_coherent", "identity_core_coherent", "identity_core_direct",
        "hard_contradiction", "technical_flag",
    ):
        scores[column] = boolean(scores[column])
    if "release_family_coherent" in scores:
        scores["release_family_coherent"] = boolean(scores.release_family_coherent)
    else:
        scores["release_family_coherent"] = scores.family_coherent
    return scores, coords


def add_broad_margins(scores: pd.DataFrame) -> pd.DataFrame:
    best = (
        scores.sort_values(
            ["cell_id", "release_broad_label", "normalized_evidence", "direct_signal", "specificity_priority", "candidate_id"],
            ascending=[True, True, False, False, False, True], kind="mergesort",
        )
        .drop_duplicates(["cell_id", "release_broad_label"])
    )
    ranked = best.sort_values(
        ["cell_id", "normalized_evidence", "direct_signal", "specificity_priority", "release_broad_label"],
        ascending=[True, False, False, False, True], kind="mergesort",
    ).copy()
    ranked["broad_rank"] = ranked.groupby("cell_id", sort=False).cumcount() + 1
    first = ranked.loc[
        ranked.broad_rank.eq(1),
        ["cell_id", "release_broad_label", "normalized_evidence"],
    ].rename(columns={
        "release_broad_label": "top_broad_label",
        "normalized_evidence": "top_broad_evidence",
    })
    second = ranked.loc[
        ranked.broad_rank.eq(2), ["cell_id", "normalized_evidence"]
    ].rename(columns={"normalized_evidence": "second_broad_evidence"})
    best = best.merge(first, on="cell_id", how="left").merge(
        second, on="cell_id", how="left"
    )
    competitor = np.where(
        best.release_broad_label.eq(best.top_broad_label),
        best.second_broad_evidence.fillna(0.0),
        best.top_broad_evidence.fillna(0.0),
    )
    best["broad_evidence_margin"] = best.normalized_evidence - competitor
    return best


def graph_components(neighbors: np.ndarray, nodes: set[int]) -> list[list[int]]:
    seen: set[int] = set()
    result: list[list[int]] = []
    for start in sorted(nodes):
        if start in seen:
            continue
        seen.add(start)
        stack = [start]
        members: list[int] = []
        while stack:
            current = stack.pop()
            members.append(current)
            for other in np.atleast_1d(neighbors[current]):
                other = int(other)
                if other in nodes and other not in seen:
                    seen.add(other)
                    stack.append(other)
        result.append(members)
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", required=True, type=Path)
    ap.add_argument("--stage-authority", required=True, type=Path)
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--threshold-registry", required=True, type=Path)
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument("--scores", required=True, action="append", type=Path)
    ap.add_argument("--cluster-evidence", required=True, action="append", type=Path)
    ap.add_argument("--round-index", type=int, default=1)
    ap.add_argument("--previous-review-manifest", type=Path)
    ap.add_argument(
        "--prior-decision-validation", action="append", type=Path, default=[]
    )
    ap.add_argument("--zero-census-direct-challenger-manifest", type=Path)
    ap.add_argument(
        "--manual-biological-adjudication", action="append", type=Path,
        default=[],
    )
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    validate_authority(args)
    thresholds = load_controller_thresholds(args.threshold_registry)
    policy = thresholds["catalog_wide_lineage_review_policy"]
    if args.round_index < 1:
        raise SystemExit("catalog-wide review round index must be >=1")
    prior_closed = prior_closed_units(args.prior_decision_validation)
    prior_scopes = prior_broad_scopes(args.prior_decision_validation)

    membership = pd.read_csv(
        args.membership, sep="\t", dtype=str, low_memory=False
    ).fillna("")
    if membership.empty or membership.cell_id.eq("").any() or membership.cell_id.duplicated().any():
        raise SystemExit("review membership must contain unique nonempty cell_id")
    for column in ("source_boundary", "source_cluster", "final_broad_label"):
        if column not in membership:
            raise SystemExit(f"review membership lacks {column}")
        membership[column] = membership[column].astype(str)
    manual_closed, manual_records = manual_closed_units(
        args.manual_biological_adjudication, membership
    )
    prior_closed.update(manual_closed)
    zero_direct_by_broad, zero_direct_summary, zero_direct_record = (
        load_zero_census_challengers(
            args.zero_census_direct_challenger_manifest, membership
        )
    )

    candidates = catalog_candidates(json.loads(args.catalog.read_text(encoding="utf-8")))
    context_summary = apply_candidate_context(
        candidates, read_tsv(args.context_evidence) if args.context_evidence else []
    )
    release_candidates = {
        candidate_id: candidate for candidate_id, candidate in candidates.items()
        if candidate_can_support_broad_review(candidate)
        and str(candidate.get("release_broad_label", ""))
        and candidate_can_release(candidate)
    }
    broad_release_labels = sorted({
        str(candidate.get("release_broad_label", ""))
        for candidate in candidates.values()
        if str(candidate.get("candidate_role", "")) == "broad"
    })
    eligible_broad_labels = {
        str(candidate.get("release_broad_label", ""))
        for candidate in candidates.values()
        if str(candidate.get("candidate_role", "")) == "broad"
        and candidate_can_release(candidate)
    }
    broad_candidate_ids: dict[str, set[str]] = defaultdict(set)
    for candidate_id, candidate in release_candidates.items():
        broad_candidate_ids[str(candidate["release_broad_label"])].add(candidate_id)

    scores, coords = load_scores(
        args.scores, release_candidates, max(1, args.workers)
    )
    if set(coords.cell_id) != set(membership.cell_id):
        raise SystemExit("review score coordinates do not exactly cover membership")
    scores["release_broad_label"] = scores.candidate_id.map(
        lambda value: str(release_candidates[value]["release_broad_label"])
    )
    scores["specificity_priority"] = scores.candidate_id.map(
        lambda value: int(release_candidates[value].get("specificity_priority", 0))
    )
    best = add_broad_margins(scores)
    best = best.merge(
        membership[["cell_id", "source_boundary", "source_cluster", "final_broad_label"]],
        on="cell_id", how="left", suffixes=("", "_membership"), validate="many_to_one",
    )
    for column in ("source_boundary", "source_cluster"):
        alternate = f"{column}_membership"
        if alternate in best and not best[column].astype(str).equals(best[alternate].astype(str)):
            raise SystemExit("score source partition differs from reviewed membership")

    group_detected: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    group_detected_rows: dict[
        tuple[str, str, str], list[tuple[str, dict[str, str]]]
    ] = defaultdict(list)
    for path in args.cluster_evidence:
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if "resolution_role" in frame:
            frame = frame.loc[frame.resolution_role.eq("selected")]
        for row in frame.to_dict("records"):
            candidate_id = str(row.get("candidate_id", ""))
            candidate = release_candidates.get(candidate_id)
            if not candidate or not group_candidate_detected(row, candidate):
                continue
            key = (
                str(row.get("source_boundary", "")),
                str(row.get("source_cluster", "")),
                str(candidate.get("release_broad_label", "")),
            )
            group_detected[key].add(candidate_id)
            if independent_group_program(row, candidate):
                group_detected_rows[key].append((candidate_id, row))

    min_group_n = int(policy["minimum_current_label_observations_per_source_group"])
    min_precision_support = float(policy["minimum_present_label_supported_fraction"])
    min_competitor = float(policy["minimum_present_label_competitor_fraction"])
    support_margin = float(policy["minimum_recall_support_evidence_margin"])
    seed_margin = float(policy["minimum_recall_seed_evidence_margin"])
    min_component = int(policy["minimum_recall_component_members"])
    min_seed_n = int(policy["minimum_recall_direct_seed_members"])
    min_seed_fraction = float(policy["minimum_recall_direct_seed_fraction"])
    knn_k = int(policy["spatial_knn_k"])
    min_group_watch_fraction = float(
        policy["minimum_group_watch_identity_fraction"]
    )
    min_group_watch_deficit = float(
        policy["minimum_group_watch_fraction_deficit"]
    )
    min_group_watch_missing_n = int(
        policy["minimum_group_watch_expected_missing_observations"]
    )

    best["observation_support"] = (
        best.identity_core_coherent
        & best.release_family_coherent
        & best.positive_family_count.ge(2)
        & ~best.hard_contradiction
        & ~best.technical_flag
        & best.program_score.gt(0)
        & best.broad_evidence_margin.ge(support_margin)
    )
    best["direct_seed"] = (
        best.observation_support
        & best.identity_core_direct
        & best.broad_evidence_margin.ge(seed_margin)
    )

    precision_rows: list[dict[str, object]] = []
    precision_queue: list[dict[str, object]] = []
    for broad in broad_release_labels:
        current = membership.loc[membership.final_broad_label.eq(broad)]
        if current.empty or broad not in eligible_broad_labels:
            continue
        target = best.loc[best.release_broad_label.eq(broad)].set_index("cell_id")
        current_group_sizes = current.groupby(
            ["source_boundary", "source_cluster"], sort=True
        ).size()
        full_group_sizes = membership.groupby(
            ["source_boundary", "source_cluster"], sort=True
        ).size()
        for (boundary, cluster), n_current in current_group_sizes.items():
            ids = set(current.loc[
                current.source_boundary.eq(boundary)
                & current.source_cluster.eq(cluster), "cell_id"
            ])
            target_group = target.loc[target.index.intersection(ids)]
            observation_supported = set(
                target_group.index[target_group.observation_support]
            )
            group_program = bool(group_detected.get((boundary, cluster, broad)))
            current_fraction = n_current / int(full_group_sizes.loc[(boundary, cluster)])
            inherited = group_program and current_fraction >= 0.50
            effective_supported_n = (
                n_current if inherited else len(observation_supported)
            )
            supported_fraction = effective_supported_n / n_current
            competitors = best.loc[
                best.cell_id.isin(ids)
                & best.direct_seed
                & ~best.release_broad_label.eq(broad)
            ]
            competitor_counts = Counter(competitors.release_broad_label)
            strongest_competitor, strongest_n = (
                sorted(competitor_counts.items(), key=lambda item: (-item[1], item[0]))[0]
                if competitor_counts else ("", 0)
            )
            competitor_fraction = strongest_n / n_current
            status = "supported"
            reasons: list[str] = []
            if n_current >= min_group_n and supported_fraction < min_precision_support:
                status = "review_required"
                reasons.append("current_label_lacks_target_identity_support")
            if n_current >= min_group_n and competitor_fraction >= min_competitor:
                status = "review_required"
                reasons.append("current_label_contains_competing_direct_identity")
            row = {
                "broad_label": broad,
                "source_boundary": boundary,
                "source_cluster": cluster,
                "n_current_label": int(n_current),
                "source_group_n": int(full_group_sizes.loc[(boundary, cluster)]),
                "current_label_fraction": current_fraction,
                "group_program_detected": str(group_program).lower(),
                "group_program_candidate_ids": ";".join(sorted(group_detected.get((boundary, cluster, broad), set()))),
                "observation_supported_n": len(observation_supported),
                "effective_supported_n": effective_supported_n,
                "effective_supported_fraction": supported_fraction,
                "sparse_inheritance_used": str(inherited).lower(),
                "strongest_competing_broad": strongest_competitor,
                "strongest_competing_direct_n": strongest_n,
                "strongest_competing_direct_fraction": competitor_fraction,
                "status": status,
                "reason": ";".join(reasons) or "source_group_or_observation_identity_supports_current_label",
                "unit_signature": unit_signature(
                    "present_label_precision", broad, sorted(ids)
                ),
            }
            if (
                status == "review_required"
                and (
                    "present_label_precision", broad,
                    str(row["unit_signature"]),
                ) in prior_closed
            ):
                row["status"] = "supported_after_prior_biological_review"
                row["reason"] = "prior_exact_unit_review_retained_current_label"
            precision_rows.append(row)
            if row["status"] == "review_required":
                precision_queue.append(row)

    coordinate_index = coords.set_index("cell_id")
    recall_components: list[dict[str, object]] = []
    recall_members: list[dict[str, object]] = []
    component_counter = 0
    generic_recall = bool(policy["generic_remainder_recall_components_enabled"])
    for (boundary, cluster), source_members in membership.groupby(
        ["source_boundary", "source_cluster"], sort=True
    ):
        group_ids = list(source_members.cell_id)
        group_coords = coordinate_index.loc[group_ids, ["x", "y"]].astype(float)
        if len(group_ids) > 1:
            k = min(knn_k, len(group_ids) - 1)
            tree = cKDTree(group_coords.to_numpy())
            neighbors = tree.query(
                group_coords.to_numpy(), k=k + 1, workers=max(1, args.workers)
            )[1][:, 1:]
        else:
            neighbors = np.empty((len(group_ids), 0), dtype=int)
        index_by_cell = {cell: index for index, cell in enumerate(group_ids)}
        source_best = best.loc[
            best.source_boundary.eq(boundary) & best.source_cluster.eq(cluster)
        ]
        for broad in sorted(eligible_broad_labels):
            candidate_ids = broad_candidate_ids.get(broad, set())
            if not generic_recall and candidate_ids and candidate_ids <= GENERIC_REMAINDER_IDS:
                continue
            challenger = source_best.loc[
                source_best.release_broad_label.eq(broad)
                & ~source_best.final_broad_label.eq(broad)
                & source_best.observation_support
            ]
            if challenger.empty:
                continue
            seed_ids = set(challenger.loc[challenger.direct_seed, "cell_id"])
            nodes = {index_by_cell[cell] for cell in challenger.cell_id}
            for component in graph_components(neighbors, nodes):
                ids = [group_ids[index] for index in component]
                n_seed = len(set(ids) & seed_ids)
                seed_fraction = n_seed / len(ids)
                if (
                    len(ids) < min_component
                    or n_seed < min_seed_n
                    or seed_fraction < min_seed_fraction
                ):
                    continue
                frame = challenger.loc[challenger.cell_id.isin(ids)]
                component_counter += 1
                component_id = f"catalog_recall_r{args.round_index}__{component_counter:05d}"
                current_census = Counter(frame.final_broad_label.replace("", "unresolved_biological"))
                candidate_census = Counter(frame.candidate_id)
                group_program = bool(group_detected.get((boundary, cluster, broad)))
                record = {
                    "component_id": component_id,
                    "review_round": args.round_index,
                    "target_broad_label": broad,
                    "source_boundary": boundary,
                    "source_cluster": cluster,
                    "n_observations": len(ids),
                    "n_direct_seeds": n_seed,
                    "direct_seed_fraction": seed_fraction,
                    "median_normalized_evidence": float(frame.normalized_evidence.median()),
                    "median_direct_signal": float(frame.direct_signal.median()),
                    "median_broad_evidence_margin": float(frame.broad_evidence_margin.median()),
                    "group_program_detected": str(group_program).lower(),
                    "group_program_candidate_ids": ";".join(sorted(group_detected.get((boundary, cluster, broad), set()))),
                    "candidate_id_census": ";".join(f"{key}:{value}" for key, value in sorted(candidate_census.items())),
                    "current_broad_census": ";".join(f"{key}:{value}" for key, value in sorted(current_census.items())),
                    "status": "review_required",
                    "reason": "coherent_direct_seeded_candidate_component_outside_current_broad",
                    "unit_signature": unit_signature(
                        "outside_label_recall", broad, ids
                    ),
                }
                if (
                    "outside_label_recall", broad,
                    str(record["unit_signature"]),
                ) in prior_closed:
                    record["status"] = "refuted_after_prior_biological_review"
                    record["reason"] = "prior_exact_component_review_rejected_or_retained_parent"
                recall_components.append(record)
                for cell in sorted(ids):
                    source = frame.loc[frame.cell_id.eq(cell)].sort_values(
                        ["normalized_evidence", "direct_signal", "candidate_id"],
                        ascending=[False, False, True], kind="mergesort",
                    ).iloc[0]
                    recall_members.append({
                        "component_id": component_id,
                        "cell_id": cell,
                        "target_broad_label": broad,
                        "current_broad_label": source.final_broad_label,
                        "candidate_id": source.candidate_id,
                        "normalized_evidence": source.normalized_evidence,
                        "direct_signal": source.direct_signal,
                        "broad_evidence_margin": source.broad_evidence_margin,
                        "direct_seed": str(bool(source.direct_seed)).lower(),
                    })

    # A direct-seeded spatial component is intentionally strict.  It must not
    # be the only recall path: a real minority lineage can be visible in
    # subcluster DEG/pseudobulk and cross-resolution evidence while old
    # observation scores remain sparse.  Such cases become bounded source-
    # subcluster watches.  They authorize raw-count review only, never a whole-
    # group label assignment.
    active_component_keys = {
        (
            str(row["source_boundary"]), str(row["source_cluster"]),
            str(row["target_broad_label"]),
        )
        for row in recall_components if row["status"] == "review_required"
    }
    group_watch_rows: list[dict[str, object]] = []
    full_group_sizes = membership.groupby(
        ["source_boundary", "source_cluster"], sort=True
    ).size()
    current_target_sizes = membership.loc[
        membership.final_broad_label.ne("")
    ].groupby(
        ["source_boundary", "source_cluster", "final_broad_label"], sort=True
    ).size()
    for key in sorted(group_detected_rows):
        boundary, cluster, broad = key
        if broad not in eligible_broad_labels or key in active_component_keys:
            continue
        candidate_ids = broad_candidate_ids.get(broad, set())
        if not generic_recall and candidate_ids and candidate_ids <= GENERIC_REMAINDER_IDS:
            continue
        source_group_n = int(full_group_sizes.get((boundary, cluster), 0))
        if source_group_n < min_group_n:
            continue
        evidence = sorted(
            group_detected_rows[key],
            key=lambda item: (
                -group_identity_core_fraction(item[1]),
                -group_identity_core_direct_fraction(item[1]),
                -number(item[1].get("marker_deg_log2fc_mean")),
                item[0],
            ),
        )[0]
        strongest_candidate_id, strongest_row = evidence
        expected_fraction = group_identity_core_fraction(strongest_row)
        expected_direct_fraction = group_identity_core_direct_fraction(strongest_row)
        current_n = int(current_target_sizes.get((boundary, cluster, broad), 0))
        current_fraction = current_n / source_group_n
        deficit = expected_fraction - current_fraction
        expected_missing_n = max(
            0, int(round(expected_fraction * source_group_n)) - current_n
        )
        if (
            expected_fraction < min_group_watch_fraction
            or deficit < min_group_watch_deficit
            or expected_missing_n < min_group_watch_missing_n
        ):
            continue
        group_ids = sorted(membership.loc[
            membership.source_boundary.eq(boundary)
            & membership.source_cluster.eq(cluster), "cell_id"
        ])
        record = {
            "review_round": args.round_index,
            "target_broad_label": broad,
            "source_boundary": boundary,
            "source_cluster": cluster,
            "source_group_n": source_group_n,
            "current_target_n": current_n,
            "current_target_fraction": current_fraction,
            "expected_identity_core_fraction": expected_fraction,
            "expected_identity_core_direct_fraction": expected_direct_fraction,
            "expected_missing_n": expected_missing_n,
            "identity_fraction_deficit": deficit,
            "strongest_candidate_id": strongest_candidate_id,
            "group_program_candidate_ids": ";".join(sorted(group_detected[key])),
            "status": "review_required",
            "reason": "subcluster_program_exceeds_current_target_membership_without_direct_seeded_component",
            "unit_signature": unit_signature(
                "outside_label_group_watch", broad, group_ids
            ),
        }
        if (
            "outside_label_group_watch", broad,
            str(record["unit_signature"]),
        ) in prior_closed:
            record["status"] = "refuted_after_prior_biological_review"
            record["reason"] = "prior_exact_source_group_review_rejected_or_retained_parent"
        group_watch_rows.append(record)

    recall_by_broad = Counter(
        row["target_broad_label"] for row in recall_components
        if row["status"] == "review_required"
    )
    recall_n_by_broad = Counter()
    for row in recall_components:
        if row["status"] == "review_required":
            recall_n_by_broad[str(row["target_broad_label"])] += int(row["n_observations"])
    group_watch_by_broad = Counter(
        row["target_broad_label"] for row in group_watch_rows
        if row["status"] == "review_required"
    )
    precision_review_by_broad = Counter(row["broad_label"] for row in precision_queue)
    broad_census = Counter(membership.final_broad_label)
    component_members_by_broad: dict[str, set[str]] = defaultdict(set)
    component_id_to_broad = {
        str(row["component_id"]): str(row["target_broad_label"])
        for row in recall_components if row["status"] == "review_required"
    }
    for row in recall_members:
        broad = component_id_to_broad.get(str(row["component_id"]), "")
        if broad:
            component_members_by_broad[broad].add(str(row["cell_id"]))
    watch_groups_by_broad: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for row in group_watch_rows:
        if row["status"] == "review_required":
            watch_groups_by_broad[str(row["target_broad_label"])].add((
                str(row["source_boundary"]), str(row["source_cluster"]),
            ))

    # The user-facing and decision-level unit is one broad lineage.  Source
    # subclusters, direct components and group watches remain child evidence
    # and exact patch bounds; they are never exposed as hundreds of separate
    # biological decisions.
    type_review_rows: list[dict[str, object]] = []
    review_scope_rows: list[dict[str, object]] = []
    queue_rows: list[dict[str, object]] = []
    review_index = 0
    type_status_by_broad: dict[str, str] = {}
    current_broad_by_cell = dict(zip(
        membership.cell_id.astype(str), membership.final_broad_label.astype(str)
    ))
    all_query_ids = sorted(current_broad_by_cell)
    global_membership_semantic_sha256 = deterministic_membership_hash(
        membership.to_dict("records")
    )
    for broad in broad_release_labels:
        if broad not in eligible_broad_labels:
            type_status_by_broad[broad] = "not_evaluable"
            continue
        current_ids = set(membership.loc[
            membership.final_broad_label.eq(broad), "cell_id"
        ].astype(str))
        component_ids = set(component_members_by_broad.get(broad, set()))
        zero_direct_ids = (
            set(zero_direct_by_broad.get(broad, set())) if not current_ids else set()
        )
        watch_ids: set[str] = set()
        for boundary, cluster in watch_groups_by_broad.get(broad, set()):
            watch_ids.update(membership.loc[
                membership.source_boundary.eq(boundary)
                & membership.source_cluster.eq(cluster), "cell_id"
            ].astype(str))
        child_issue_n = (
            precision_review_by_broad[broad] + recall_by_broad[broad]
            + group_watch_by_broad[broad]
        )
        if not current_ids and not component_ids and not watch_ids and not zero_direct_ids:
            type_status_by_broad[broad] = "supported_zero_census_no_query_challenger"
            continue
        review_mode = (
            "broad_lineage_review" if current_ids else "missing_broad_review"
        )
        signature_tokens = (
            [f"current:{cell}" for cell in sorted(current_ids)]
            + [f"component:{cell}" for cell in sorted(component_ids)]
            + [f"zero_direct:{cell}" for cell in sorted(zero_direct_ids)]
            + [
                f"watch_group:{boundary}:{cluster}"
                for boundary, cluster in sorted(watch_groups_by_broad.get(broad, set()))
            ]
        )
        signature = unit_signature(review_mode, broad, signature_tokens)
        exact_key = (review_mode, broad, signature)
        closed = exact_key in prior_closed
        manual_closed_scope = exact_key in manual_closed
        monotonic_subtraction = False
        monotonic_removed_n = 0
        monotonic_removed_fraction = 0.0
        previous_scope = prior_scopes.get(broad, {})
        previous_current = set(previous_scope.get("current_ids", set()))
        if (
            not closed
            and review_mode == "broad_lineage_review"
            and current_ids
            and previous_current
            and current_ids < previous_current
            and component_ids.issubset(set(previous_scope.get("component_ids", set())))
            and watch_ids.issubset(set(previous_scope.get("watch_ids", set())))
            and zero_direct_ids.issubset(
                set(previous_scope.get("zero_direct_ids", set()))
            )
        ):
            monotonic_removed_n = len(previous_current - current_ids)
            monotonic_removed_fraction = monotonic_removed_n / len(previous_current)
            monotonic_subtraction = monotonic_removed_fraction <= float(
                policy[
                    "maximum_monotonic_subtraction_fraction_without_full_reopen"
                ]
            )
        # The review still scans the whole query, but closure is invalidated by
        # a change to this target's current members, recall components or source
        # watches.  An unrelated lineage patch must not reopen every already
        # closed type.  Patches that move cells into or out of this target alter
        # these tokens and therefore reopen it deterministically.
        open_review = not closed and not monotonic_subtraction
        review_id = ""
        if open_review:
            review_index += 1
            review_id = f"catalog_review_r{args.round_index}__{review_index:05d}"
            queue_rows.append({
                "review_id": review_id,
                "review_round": args.round_index,
                "review_mode": review_mode,
                "target_broad_label": broad,
                "source_boundary": "*",
                "source_cluster": "*",
                "component_id": "",
                "n_observations": len(all_query_ids),
                "reason": (
                    f"complete_cell_type_review:precision_groups={precision_review_by_broad[broad]};"
                    f"recall_components={recall_by_broad[broad]};"
                    f"group_watches={group_watch_by_broad[broad]};"
                    f"zero_census_direct={len(zero_direct_ids)}"
                ),
                "unit_signature": signature,
                "required_review": "query_raw_count_marker_families;target_vs_outside_DEG_pseudobulk;pairwise_competitors;whole_section_spatial_distribution;over_recall;under_recall;targeted_membership_if_supported",
            })
            for cell in all_query_ids:
                roles: list[str] = []
                if cell in current_ids:
                    roles.append("current_label")
                else:
                    roles.append("whole_query_recall_scan")
                if cell in component_ids:
                    roles.append("direct_recall_component")
                if cell in watch_ids:
                    roles.append("group_watch_source")
                if cell in zero_direct_ids:
                    roles.append("zero_census_direct_multifamily_challenger")
                review_scope_rows.append({
                    "review_id": review_id,
                    "cell_id": cell,
                    "target_broad_label": broad,
                    "scope_roles": ";".join(roles),
                    "current_broad_label": current_broad_by_cell[cell],
                })
        status = (
            "review_required" if open_review
            else "closed_by_manual_adjudication" if manual_closed_scope
            else "supported_after_exact_cell_type_review" if closed
            else "closed_after_monotonic_subtraction" if monotonic_subtraction
            else "supported_after_membership_change_reaudit"
        )
        type_status_by_broad[broad] = status
        type_review_rows.append({
            "broad_label": broad,
            "review_mode": review_mode,
            "current_label_n": len(current_ids),
            "precision_child_source_group_n": precision_review_by_broad[broad],
            "recall_child_component_n": recall_by_broad[broad],
            "recall_child_group_watch_n": group_watch_by_broad[broad],
            "zero_census_direct_challenger_n": len(zero_direct_ids),
            "monotonic_removed_n": monotonic_removed_n,
            "monotonic_removed_fraction": monotonic_removed_fraction,
            "review_scope_n": len(all_query_ids),
            "status": status,
            "review_id": review_id,
            "unit_signature": signature,
        })

    matrix_rows: list[dict[str, object]] = []
    for broad in broad_release_labels:
        context_candidates = [
            candidate for candidate in candidates.values()
            if str(candidate.get("candidate_role", "")) == "broad"
            and str(candidate.get("release_broad_label", "")) == broad
        ]
        evaluable = broad in eligible_broad_labels
        matrix_rows.append({
            "broad_label": broad,
            "context_status": (
                "evaluable" if evaluable else ";".join(sorted({
                    str(candidate.get("_context_status", "not_evaluable"))
                    for candidate in context_candidates
                }))
            ),
            "final_n_observations": int(broad_census[broad]),
            "precision_source_group_n": sum(
                row["broad_label"] == broad for row in precision_rows
            ),
            "precision_review_source_group_n": precision_review_by_broad[broad],
            "recall_challenger_component_n": recall_by_broad[broad],
            "recall_challenger_observation_n": recall_n_by_broad[broad],
            "recall_group_watch_n": group_watch_by_broad[broad],
            "zero_census_direct_challenger_n": len(
                zero_direct_by_broad.get(broad, set())
            ),
            "primary_review_unit": "broad_lineage",
            "status": type_status_by_broad.get(
                broad, "not_evaluable" if not evaluable else "supported"
            ),
        })

    args.out.mkdir(parents=True, exist_ok=True)
    matrix_path = args.out / "catalog_wide_lineage_review_matrix.tsv"
    precision_path = args.out / "present_label_precision_audit.tsv"
    component_path = args.out / "outside_label_recall_components.tsv"
    component_membership_path = args.out / "outside_label_recall_component_membership.tsv.gz"
    group_watch_path = args.out / "outside_label_group_watch.tsv"
    type_review_path = args.out / "broad_lineage_review_summary.tsv"
    review_scope_path = args.out / "broad_lineage_review_scope_membership.tsv.gz"
    queue_path = args.out / "catalog_wide_lineage_review_queue.tsv"
    pd.DataFrame(matrix_rows).to_csv(matrix_path, sep="\t", index=False)
    pd.DataFrame(precision_rows, columns=PRECISION_COLUMNS).to_csv(
        precision_path, sep="\t", index=False
    )
    pd.DataFrame(recall_components, columns=COMPONENT_COLUMNS).to_csv(
        component_path, sep="\t", index=False
    )
    pd.DataFrame(recall_members, columns=[
        "component_id", "cell_id", "target_broad_label", "current_broad_label",
        "candidate_id", "normalized_evidence", "direct_signal",
        "broad_evidence_margin", "direct_seed",
    ]).to_csv(component_membership_path, sep="\t", index=False, compression="gzip")
    pd.DataFrame(group_watch_rows, columns=GROUP_WATCH_COLUMNS).to_csv(
        group_watch_path, sep="\t", index=False
    )
    pd.DataFrame(type_review_rows).to_csv(
        type_review_path, sep="\t", index=False
    )
    pd.DataFrame(review_scope_rows, columns=[
        "review_id", "cell_id", "target_broad_label", "scope_roles",
        "current_broad_label",
    ]).to_csv(review_scope_path, sep="\t", index=False, compression="gzip")
    pd.DataFrame(queue_rows, columns=[
        "review_id", "review_round", "review_mode", "target_broad_label",
        "source_boundary", "source_cluster", "component_id", "n_observations",
        "reason", "unit_signature", "required_review",
    ]).to_csv(queue_path, sep="\t", index=False)
    previous = None
    if args.previous_review_manifest:
        previous = artifact(args.previous_review_manifest)
    manifest = {
        "schema_version": "2.2",
        "status": "ITERATION_REQUIRED" if queue_rows else "PASS",
        "stage": "post_atlas_catalog_wide_lineage_review",
        "user_facing_stage_name": "逐大类全样本复核",
        "artifact_role": "biological_quality_review",
        "review_round": args.round_index,
        "maximum_decision_rounds": int(policy["maximum_decision_rounds"]),
        "formal_membership_written": False,
        "catalog_wide_double_sided_review": True,
        "broad_lineage_is_primary_review_unit": True,
        "source_subclusters_are_internal_evidence_only": True,
        "whole_object_per_cell_classifier_used": False,
        "geometry_can_trigger_review_but_cannot_assign_labels": True,
        "historical_labels_used_as_identity_evidence": False,
        "membership": artifact(args.membership),
        "membership_semantic_sha256": global_membership_semantic_sha256,
        "annotation_contract": artifact(args.contract),
        "stage_authority": artifact(args.stage_authority),
        "candidate_catalog": artifact(args.catalog),
        "threshold_registry": artifact(args.threshold_registry),
        "context_evidence": artifact(args.context_evidence) if args.context_evidence else None,
        "context_release_eligibility": context_summary,
        "observation_scores": [artifact(path) for path in args.scores],
        "cluster_evidence": [artifact(path) for path in args.cluster_evidence],
        "previous_review_manifest": previous,
        "prior_decision_validations": [
            artifact(path) for path in args.prior_decision_validation
        ],
        "manual_biological_adjudications": manual_records,
        "zero_census_direct_challenger": zero_direct_record,
        "eligible_broad_n": len(eligible_broad_labels),
        "context_not_evaluable_broad_n": len(set(broad_release_labels) - eligible_broad_labels),
        "precision_review_source_group_n": len(precision_queue),
        "recall_challenger_component_n": sum(
            row["status"] == "review_required" for row in recall_components
        ),
        "recall_group_watch_n": sum(
            row["status"] == "review_required" for row in group_watch_rows
        ),
        "zero_census_direct_challenger_broad_n": sum(
            str(row.get("status", "")) == "review_required"
            for row in zero_direct_summary.values()
        ),
        "cell_type_review_n": len(type_review_rows),
        "review_queue_n": len(queue_rows),
        "artifacts": {
            "lineage_review_matrix": artifact(matrix_path),
            "present_label_precision_audit": artifact(precision_path),
            "outside_label_recall_components": artifact(component_path),
            "outside_label_recall_component_membership": artifact(component_membership_path),
            "outside_label_group_watch": artifact(group_watch_path),
            "broad_lineage_review_summary": artifact(type_review_path),
            "broad_lineage_review_scope_membership": artifact(review_scope_path),
            "review_queue": artifact(queue_path),
        },
    }
    manifest_path = args.out / "catalog_wide_lineage_review_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    schema = Path(__file__).resolve().parents[1] / "schemas/catalog_wide_lineage_review.schema.json"
    _, schema_errors = validate_json_against_schema(manifest_path, schema)
    if schema_errors:
        raise SystemExit(
            "catalog-wide review manifest violates its schema: "
            + "; ".join(schema_errors)
        )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0 if not queue_rows else 2


if __name__ == "__main__":
    raise SystemExit(main())
