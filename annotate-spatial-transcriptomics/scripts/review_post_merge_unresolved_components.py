#!/usr/bin/env python3
"""Resolve only bounded post-merge biological remainders.

The second-round Leiden subcluster remains the primary annotation unit.  This
stage is allowed to act only on observations that remain unresolved after
second-round local splitting and Atlas routing.  It forms candidate-specific
spatial components from strict direct multi-family seeds, expands through
coherent low-RNA support members, resolves overlaps without candidate-order
writeback, and processes generic remainder candidates last.
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from evidence_schema_lib import sha256
from lineage_controller_lib import (
    deterministic_membership_hash,
    group_candidate_detected,
    group_candidate_score,
)


def artifact(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "n_bytes": path.stat().st_size,
    }


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"cell_id": str}, **kwargs)


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes", "t"})


def validate_authority(args: argparse.Namespace) -> dict:
    authority = json.loads(args.stage_authority.read_text(encoding="utf-8"))
    if (
        authority.get("mode") != "stage_authority"
        or authority.get("phase") != "atlas_and_completeness_review"
        or authority.get("annotation_contract_sha256") != sha256(args.contract)
    ):
        raise SystemExit("stage authority does not permit post-merge unresolved review")
    record = authority.get("post_atlas_membership", {})
    if (
        Path(str(record.get("path", ""))).resolve() != args.membership.resolve()
        or record.get("sha256") != sha256(args.membership)
    ):
        raise SystemExit("post-Atlas membership differs from stage authority")
    bound_scores = {
        (Path(str(row.get("path", ""))).resolve(), str(row.get("sha256", "")))
        for row in authority.get("observation_scores", [])
    }
    supplied_scores = {(path.resolve(), sha256(path)) for path in args.scores}
    if not supplied_scores or bound_scores != supplied_scores:
        raise SystemExit("observation scores differ from stage authority")
    bound_evidence = {
        (Path(str(row.get("path", ""))).resolve(), str(row.get("sha256", "")))
        for row in authority.get("cluster_evidence", [])
    }
    supplied_evidence = {
        (path.resolve(), sha256(path)) for path in args.cluster_evidence
    }
    if not supplied_evidence or bound_evidence != supplied_evidence:
        raise SystemExit("cluster evidence differs from stage authority")
    return authority


def candidate_catalog(path: Path) -> tuple[dict[str, dict], set[str]]:
    catalog = json.loads(path.read_text(encoding="utf-8"))
    candidates: dict[str, dict] = {}
    generic: set[str] = set()
    for row in catalog.get("candidate_boundaries", []):
        strategy = str(row.get("writeback_strategy", ""))
        if (
            row.get("candidate_role") not in {"broad", "fine"}
            or not row.get("release_broad_label")
            or strategy.startswith("watch_only")
            or strategy == "canonical_cluster_membership"
        ):
            continue
        candidates[str(row["candidate_id"])] = row
        if strategy == "generic_exact_remainder_after_specific_lineages":
            generic.add(str(row["candidate_id"]))
    if not candidates or not generic:
        raise SystemExit("candidate catalog lacks releasable or generic remainder candidates")
    return candidates, generic


def top_two_by_cell(frame: pd.DataFrame, value: str, prefix: str) -> pd.DataFrame:
    ranked = frame.sort_values(
        ["cell_id", value, "direct_signal", "specificity_priority", "release_broad_label"],
        ascending=[True, False, False, False, True],
        kind="mergesort",
    ).copy()
    ranked["rank"] = ranked.groupby("cell_id", sort=False).cumcount() + 1
    first = ranked.loc[ranked["rank"].eq(1), ["cell_id", "release_broad_label", value]].rename(
        columns={"release_broad_label": f"{prefix}_top_label", value: f"{prefix}_top_value"}
    )
    second = ranked.loc[ranked["rank"].eq(2), ["cell_id", value]].rename(
        columns={value: f"{prefix}_second_value"}
    )
    result = first.merge(second, on="cell_id", how="left")
    result[f"{prefix}_second_value"] = result[f"{prefix}_second_value"].fillna(0.0)
    return result


def add_candidate_margins(scores: pd.DataFrame, generic_ids: set[str]) -> pd.DataFrame:
    best_by_broad = (
        scores.sort_values(
            ["cell_id", "release_broad_label", "normalized_evidence", "direct_signal", "specificity_priority"],
            ascending=[True, True, False, False, False], kind="mergesort",
        )
        .drop_duplicates(["cell_id", "release_broad_label"])
    )
    best_by_broad["is_generic"] = best_by_broad.candidate_id.isin(generic_ids)
    specific = best_by_broad.loc[~best_by_broad.is_generic]
    all_ev = top_two_by_cell(best_by_broad, "normalized_evidence", "all_ev")
    all_direct = top_two_by_cell(best_by_broad, "direct_signal", "all_direct")
    specific_ev = top_two_by_cell(specific, "normalized_evidence", "specific_ev")
    specific_direct = top_two_by_cell(specific, "direct_signal", "specific_direct")
    scores = scores.merge(all_ev, on="cell_id", how="left").merge(
        all_direct, on="cell_id", how="left"
    ).merge(specific_ev, on="cell_id", how="left").merge(
        specific_direct, on="cell_id", how="left"
    )
    is_generic = scores.candidate_id.isin(generic_ids)
    for value, prefix in (("normalized_evidence", "ev"), ("direct_signal", "direct")):
        all_top_label = scores[f"all_{prefix}_top_label"]
        all_competitor = np.where(
            scores.release_broad_label.eq(all_top_label),
            scores[f"all_{prefix}_second_value"],
            scores[f"all_{prefix}_top_value"],
        )
        specific_top_label = scores[f"specific_{prefix}_top_label"]
        specific_competitor = np.where(
            scores.release_broad_label.eq(specific_top_label),
            scores[f"specific_{prefix}_second_value"],
            scores[f"specific_{prefix}_top_value"],
        )
        competitor = np.where(is_generic, all_competitor, specific_competitor)
        scores[f"candidate_{prefix}_margin"] = scores[value] - np.nan_to_num(competitor, nan=0.0)
    return scores


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


def form_candidate_components(
    scores: pd.DataFrame,
    candidates: dict[str, dict],
    generic_ids: set[str],
    workers: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    proposal_rows: list[dict] = []
    component_rows: list[dict] = []
    for (boundary, cluster), group in scores.groupby(["source_boundary", "source_cluster"], sort=True):
        full = group.drop_duplicates("cell_id")[["cell_id", "x", "y"]].reset_index(drop=True)
        index_by_id = dict(zip(full.cell_id, full.index))
        if len(full) > 1:
            k = min(12, len(full) - 1)
            tree = cKDTree(full[["x", "y"]].to_numpy(float))
            neighbors = tree.query(full[["x", "y"]].to_numpy(float), k=k + 1, workers=workers)[1][:, 1:]
        else:
            neighbors = np.empty((len(full), 0), dtype=int)
        for candidate_id, candidate in group.groupby("candidate_id", sort=True):
            support = candidate.loc[
                candidate.identity_core_coherent
                & candidate.family_coherent
                & candidate.positive_family_count.ge(2)
                & ~candidate.hard_contradiction
                & ~candidate.technical_flag
                & candidate.program_score.gt(0)
                & candidate.candidate_ev_margin.ge(-0.05)
            ]
            seed = candidate.loc[
                candidate.identity_core_coherent
                & candidate.identity_core_direct
                & candidate.family_coherent
                & candidate.positive_family_count.ge(2)
                & ~candidate.hard_contradiction
                & ~candidate.technical_flag
                & candidate.program_score.gt(0)
                & candidate.candidate_ev_margin.ge(0.10)
                & candidate.candidate_direct_margin.ge(0.03)
            ]
            nodes = {index_by_id[cell] for cell in support.cell_id if cell in index_by_id}
            seed_ids = set(seed.cell_id)
            component_number = 0
            for members in graph_components(neighbors, nodes):
                member_ids = set(full.loc[members, "cell_id"])
                n_seed = len(member_ids & seed_ids)
                seed_fraction = n_seed / len(members)
                if len(members) < 5 or n_seed < 3 or seed_fraction < 0.10:
                    continue
                component_number += 1
                component_id = f"{boundary}__{cluster}__{candidate_id}__{component_number}"
                comp = candidate.loc[candidate.cell_id.isin(member_ids)]
                component_rows.append({
                    "component_id": component_id,
                    "source_boundary": boundary,
                    "source_cluster": cluster,
                    "candidate_id": candidate_id,
                    "release_broad_label": candidates[candidate_id]["release_broad_label"],
                    "specificity_priority": int(candidates[candidate_id].get("specificity_priority", 0)),
                    "is_generic": candidate_id in generic_ids,
                    "n_observations": len(comp),
                    "n_direct_seeds": n_seed,
                    "seed_fraction": seed_fraction,
                    "median_program_score": float(comp.program_score.median()),
                    "median_candidate_margin": float(comp.candidate_ev_margin.median()),
                    "median_direct_signal": float(comp.direct_signal.median()),
                })
                for row in comp.itertuples(index=False):
                    proposal_rows.append({
                        "cell_id": row.cell_id,
                        "source_boundary": boundary,
                        "source_cluster": cluster,
                        "candidate_id": candidate_id,
                        "release_broad_label": candidates[candidate_id]["release_broad_label"],
                        "specificity_priority": int(candidates[candidate_id].get("specificity_priority", 0)),
                        "is_generic": candidate_id in generic_ids,
                        "component_id": component_id,
                        "normalized_evidence": row.normalized_evidence,
                        "direct_signal": row.direct_signal,
                        "program_score": row.program_score,
                    })
    return pd.DataFrame(proposal_rows), pd.DataFrame(component_rows)


def resolve_proposal_overlaps(proposals: pd.DataFrame) -> pd.DataFrame:
    if proposals.empty:
        return pd.DataFrame(columns=[
            "cell_id", "source_boundary", "source_cluster", "decision",
            "release_broad_label", "candidate_id", "component_id",
            "evidence_margin",
        ])
    broad = (
        proposals.sort_values(
            ["cell_id", "release_broad_label", "normalized_evidence", "direct_signal", "specificity_priority"],
            ascending=[True, True, False, False, False], kind="mergesort",
        )
        .drop_duplicates(["cell_id", "release_broad_label"])
    )
    decisions: list[dict] = []
    for cell_id, frame in broad.groupby("cell_id", sort=True):
        source_boundary = str(frame.iloc[0].source_boundary)
        source_cluster = str(frame.iloc[0].source_cluster)
        specific = frame.loc[~frame.is_generic].sort_values(
            ["normalized_evidence", "direct_signal", "specificity_priority", "release_broad_label"],
            ascending=[False, False, False, True], kind="mergesort",
        )
        generic = frame.loc[frame.is_generic]
        if len(specific) == 1:
            chosen, decision, margin = specific.iloc[0], "specific_unique_component", np.nan
        elif len(specific) > 1:
            margin = float(specific.iloc[0].normalized_evidence - specific.iloc[1].normalized_evidence)
            if margin < 0.10:
                decisions.append({"cell_id": cell_id, "source_boundary": source_boundary, "source_cluster": source_cluster, "decision": "remain_unresolved_specific_overlap", "release_broad_label": "", "candidate_id": "", "component_id": "", "evidence_margin": margin})
                continue
            chosen, decision = specific.iloc[0], "specific_overlap_resolved"
        elif len(generic):
            chosen, decision, margin = generic.iloc[0], "generic_after_specific_absent", np.nan
        else:
            continue
        decisions.append({
            "cell_id": cell_id,
            "source_boundary": source_boundary,
            "source_cluster": source_cluster,
            "decision": decision,
            "release_broad_label": chosen.release_broad_label,
            "candidate_id": chosen.candidate_id,
            "component_id": chosen.component_id,
            "evidence_margin": margin,
        })
    return pd.DataFrame(decisions)


def source_supported_parent_candidates(
    paths: list[Path], candidates: dict[str, dict]
) -> dict[tuple[str, str, str], str]:
    """Choose one detected candidate that can authorize each parent writeback.

    Spatial neighborhood establishes anatomical context, not identity.  A
    contextual parent can therefore receive a component only when the same
    selected second-round subcluster independently contains a detected program
    for that parent.  The chosen candidate replaces the lineage-of-origin
    challenger in release provenance so downstream completeness can verify the
    exact candidate x source-group support.
    """
    supported: dict[tuple[str, str, str], list[tuple[int, float, str]]] = {}
    for path in paths:
        frame = read_tsv(path).fillna("")
        if "resolution_role" not in frame or "candidate_id" not in frame:
            raise SystemExit(f"cluster evidence is malformed: {path}")
        frame = frame.loc[frame.resolution_role.eq("selected")]
        for row in frame.to_dict("records"):
            candidate_id = str(row.get("candidate_id", ""))
            candidate = candidates.get(candidate_id)
            if not candidate or not group_candidate_detected(row, candidate):
                continue
            boundary = str(row.get("source_boundary", ""))
            cluster = str(row.get("source_cluster", ""))
            broad = str(candidate.get("release_broad_label", ""))
            if not boundary or not cluster or not broad:
                continue
            broad_role = int(str(candidate.get("candidate_role", "")) == "broad")
            supported.setdefault((boundary, cluster, broad), []).append(
                (broad_role, group_candidate_score(row, candidate), candidate_id)
            )
    result: dict[tuple[str, str, str], str] = {}
    for key, rows in supported.items():
        rows.sort(key=lambda item: (-item[0], -item[1], item[2]))
        result[key] = rows[0][2]
    return result


def apply_contextual_parent_overrides(
    decisions: pd.DataFrame,
    full: pd.DataFrame,
    candidates: dict[str, dict],
    source_supported_parents: dict[tuple[str, str, str], str],
    workers: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    accepted = decisions.loc[
        ~decisions.decision.astype(str).str.startswith("remain_unresolved")
    ].copy()
    if accepted.empty:
        return decisions, pd.DataFrame()
    override_labels = sorted({
        str(rule.get("context_broad_label", ""))
        for row in candidates.values()
        for rule in row.get("contextual_parent_overrides", [])
        if rule.get("context_broad_label")
    })
    if not override_labels:
        return decisions, pd.DataFrame()
    coords = full.set_index("cell_id")[["x", "y"]]
    accepted = accepted.join(coords, on="cell_id")
    tree = cKDTree(full[["x", "y"]].to_numpy(float))
    k = min(41, len(full))
    neighbor_indices = tree.query(accepted[["x", "y"]].to_numpy(float), k=k, workers=workers)[1][:, 1:]
    labels = full.final_broad_label.to_numpy(str)
    context_rows = []
    for row_index, indices in enumerate(neighbor_indices):
        defined = labels[np.asarray(indices, dtype=int)]
        defined = defined[defined != ""]
        record = {"cell_id": accepted.iloc[row_index].cell_id, "defined_neighbor_n": len(defined)}
        for label in override_labels:
            record[f"context_fraction__{label}"] = float(np.mean(defined == label)) if len(defined) else 0.0
        context_rows.append(record)
    context = pd.DataFrame(context_rows)
    accepted = accepted.merge(context, on="cell_id", validate="one_to_one")
    override_audit: list[dict] = []
    for component_id, component in accepted.groupby("component_id", sort=True):
        candidate_id = str(component.iloc[0].candidate_id)
        for rule in candidates.get(candidate_id, {}).get("contextual_parent_overrides", []):
            context_label = str(rule.get("context_broad_label", ""))
            target_label = str(rule.get("writeback_broad_label", context_label))
            minimum = float(rule.get("minimum_component_neighbor_fraction", 0.50))
            fractions = component[f"context_fraction__{context_label}"]
            passing_fraction = float(np.mean(fractions >= minimum))
            median_fraction = float(fractions.median())
            context_pass = passing_fraction >= minimum and median_fraction >= minimum
            boundary = str(component.iloc[0].source_boundary)
            cluster = str(component.iloc[0].source_cluster)
            support_candidate_id = source_supported_parents.get(
                (boundary, cluster, target_label), ""
            )
            apply_override = context_pass and bool(support_candidate_id)
            if apply_override:
                status = "PARENT_OVERRIDE"
            elif context_pass:
                status = "PARENT_SOURCE_UNSUPPORTED"
            else:
                status = "KEEP_CANDIDATE"
            override_audit.append({
                "component_id": component_id,
                "candidate_id": candidate_id,
                "source_boundary": boundary,
                "source_cluster": cluster,
                "original_broad_label": str(component.iloc[0].release_broad_label),
                "context_broad_label": context_label,
                "writeback_broad_label": target_label,
                "n_observations": len(component),
                "member_context_pass_fraction": passing_fraction,
                "median_neighbor_context_fraction": median_fraction,
                "parent_source_supported": bool(support_candidate_id),
                "parent_support_candidate_id": support_candidate_id,
                "status": status,
            })
            mask = decisions.component_id.eq(component_id) & ~decisions.decision.astype(str).str.startswith("remain_unresolved")
            if apply_override:
                decisions.loc[mask, "release_broad_label"] = target_label
                decisions.loc[mask, "candidate_id"] = support_candidate_id
                decisions.loc[mask, "decision"] = "contextual_parent_return"
                break
            if context_pass:
                decisions.loc[mask, "release_broad_label"] = ""
                decisions.loc[mask, "candidate_id"] = ""
                decisions.loc[mask, "decision"] = "remain_unresolved_context_parent_unsupported"
                break
    return decisions, pd.DataFrame(override_audit)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", required=True, type=Path)
    ap.add_argument("--stage-authority", required=True, type=Path)
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--scores", action="append", type=Path, default=[])
    ap.add_argument("--cluster-evidence", action="append", type=Path, default=[])
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    validate_authority(args)
    candidates, generic_ids = candidate_catalog(args.catalog)
    source_supported_parents = source_supported_parent_candidates(
        args.cluster_evidence, candidates
    )
    membership = read_tsv(args.membership).fillna("")
    membership.cell_id = membership.cell_id.astype(str)
    if membership.cell_id.duplicated().any():
        raise SystemExit("post-Atlas membership duplicates cell IDs")
    unresolved = membership.loc[membership.final_state.eq("unresolved_biological")]
    unresolved_ids = set(unresolved.cell_id)
    required = {
        "cell_id", "source_boundary", "source_cluster", "candidate_id",
        "candidate_role", "release_broad_label", "direct_signal", "program_score",
        "normalized_evidence", "positive_family_count", "family_coherent",
        "identity_core_coherent", "identity_core_direct", "hard_contradiction",
        "technical_flag", "x", "y",
    }
    score_frames, coordinate_frames = [], []
    for path in args.scores:
        header = set(pd.read_csv(path, sep="\t", nrows=0).columns)
        if missing := required.difference(header):
            raise SystemExit(f"observation score lacks {sorted(missing)}: {path}")
        frame = read_tsv(path, usecols=sorted(required))
        coordinate_frames.append(frame.drop_duplicates("cell_id")[["cell_id", "x", "y"]])
        score_frames.append(frame.loc[frame.cell_id.isin(unresolved_ids) & frame.candidate_id.isin(candidates)].copy())
    scores = pd.concat(score_frames, ignore_index=True) if score_frames else pd.DataFrame()
    full_coords = pd.concat(coordinate_frames, ignore_index=True).drop_duplicates("cell_id")
    if set(full_coords.cell_id) != set(membership.cell_id):
        raise SystemExit("observation score coordinates do not cover post-Atlas membership")
    full = membership.merge(full_coords, on="cell_id", validate="one_to_one")
    if not scores.empty:
        for column in ["direct_signal", "program_score", "normalized_evidence", "positive_family_count", "x", "y"]:
            scores[column] = pd.to_numeric(scores[column], errors="coerce").fillna(0)
        for column in ["family_coherent", "identity_core_coherent", "identity_core_direct", "hard_contradiction", "technical_flag"]:
            scores[column] = as_bool(scores[column])
        scores["specificity_priority"] = scores.candidate_id.map(lambda x: int(candidates[x].get("specificity_priority", 0)))
        scores = add_candidate_margins(scores, generic_ids)
        proposals, components = form_candidate_components(scores, candidates, generic_ids, args.workers)
        decisions = resolve_proposal_overlaps(proposals)
        decisions, context_audit = apply_contextual_parent_overrides(
            decisions, full, candidates, source_supported_parents, args.workers
        )
    else:
        proposals, components, decisions, context_audit = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    accepted = decisions.loc[
        ~decisions.decision.astype(str).str.startswith("remain_unresolved")
    ].copy() if len(decisions) else decisions
    accepted_map = accepted.set_index("cell_id").to_dict("index") if len(accepted) else {}
    original_defined = membership.final_broad_label.ne("")
    output_membership = membership.copy()
    for index, row in output_membership.loc[~original_defined].iterrows():
        decision = accepted_map.get(str(row.cell_id))
        if not decision:
            continue
        output_membership.at[index, "candidate_id"] = decision["candidate_id"]
        output_membership.at[index, "final_state"] = "defined_broad_only"
        output_membership.at[index, "final_broad_label"] = decision["release_broad_label"]
        output_membership.at[index, "confidence"] = "moderate"
        output_membership.at[index, "assignment_origin"] = "post_merge_unresolved_component_review__" + decision["decision"]
        output_membership.at[index, "unresolved_reason"] = ""
        output_membership.at[index, "broad_frozen"] = "true"
        output_membership.at[index, "fine_anchor_eligible"] = "false"
    if not output_membership.loc[original_defined, "final_broad_label"].equals(membership.loc[original_defined, "final_broad_label"]):
        raise SystemExit("post-merge unresolved review altered a previously defined broad label")
    args.out.mkdir(parents=True, exist_ok=True)
    membership_path = args.out / "reviewed_post_atlas_broad_membership.tsv.gz"
    with gzip.open(membership_path, "wt", encoding="utf-8", newline="") as handle:
        output_membership.to_csv(handle, sep="\t", index=False)
    written_artifacts: dict[str, dict[str, object]] = {}
    for frame, name, compressed in (
        (proposals, "candidate_component_proposals.tsv.gz", True),
        (decisions, "candidate_component_decisions.tsv", False),
        (components, "candidate_component_summary.tsv", False),
        (context_audit, "contextual_parent_override_audit.tsv", False),
    ):
        path = args.out / name
        if compressed:
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                frame.to_csv(handle, sep="\t", index=False)
        else:
            frame.to_csv(path, sep="\t", index=False)
        written_artifacts[name] = artifact(path)
    remaining = int(output_membership.final_broad_label.eq("").sum())
    census = Counter(output_membership.final_broad_label.replace("", "unresolved_biological"))
    manifest = {
        "status": "PASS",
        "stage": "post_merge_unresolved_component_review",
        "formal_membership_written": True,
        "scope": "post_atlas_unresolved_biological_only",
        "membership": artifact(membership_path),
        "component_artifacts": written_artifacts,
        "membership_semantic_hash": deterministic_membership_hash(
            output_membership.to_dict("records")
        ),
        "n_unresolved_input": int(len(unresolved)),
        "n_broad_returns": int(len(accepted)),
        "n_remaining_unresolved": remaining,
        "remaining_unresolved_fraction": remaining / len(output_membership),
        "broad_return_census": dict(Counter(accepted.release_broad_label)) if len(accepted) else {},
        "final_broad_census": dict(census),
        "policy": "strict direct multi-family seeds; coherent local expansion; specific before generic; unresolved overlaps preserved; anatomical parent override requires matching selected-subcluster parent identity support",
    }
    (args.out / "post_merge_unresolved_review_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
