#!/usr/bin/env python3
"""Calibrate and adjudicate broad-only rescue from independent evidence channels.

The input must contain query-like held-out anchors and query observations in one
table.  No expression matrix is read and no labels are written to a project
ledger.  The output is an auditable proposal that must still be committed with
the framework's atomic state tools.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


TRUE_VALUES = {"1", "true", "yes", "y", "pass", "passed", "accepted", "accept"}


def truth(value: object) -> bool:
    return str(value).strip().lower() in TRUE_VALUES


def safe(value: object) -> str:
    return str(value).strip() if value is not None else ""


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, low_memory=False)


def channel_accepts(frame: pd.DataFrame, channel: dict) -> pd.Series:
    column = channel.get("accept_column")
    if not column:
        return frame[channel["label_column"]].astype(str).str.strip().ne("")
    values = channel.get("accepted_values")
    if values is None:
        return frame[column].map(truth)
    accepted = {str(value).strip().lower() for value in values}
    return frame[column].astype(str).str.strip().str.lower().isin(accepted)


def choose_candidate(row: pd.Series, channels: list[dict]) -> tuple[str, list[str]]:
    votes: dict[str, list[tuple[int, str]]] = {}
    for order, channel in enumerate(channels):
        if not truth(row[f"__accept__{channel['name']}"]):
            continue
        label = safe(row[channel["label_column"]])
        if not label:
            continue
        priority = int(channel.get("priority", order))
        votes.setdefault(label, []).append((priority, channel["name"]))
    if not votes:
        return "", []
    ranked = sorted(
        votes,
        key=lambda label: (-len(votes[label]), min(item[0] for item in votes[label]), label),
    )
    label = ranked[0]
    names = [name for _, name in sorted(votes[label])]
    return label, names


def calibrate_threshold(
    heldout: pd.DataFrame,
    minimum_votes: int,
    target_precision: float,
    minimum_selected: int,
) -> dict | None:
    if heldout.empty:
        return None
    candidates = []
    # Calibrate only at support counts actually represented in the held-out
    # anchors.  A lower unobserved threshold would silently extrapolate to a
    # query evidence pattern that was never precision-tested.
    observed_thresholds = sorted(
        int(value) for value in heldout["support_channel_n"].unique() if int(value) >= minimum_votes
    )
    for threshold in observed_thresholds:
        selected = heldout[heldout["joint_gate_pass"] & heldout["support_channel_n"].ge(threshold)]
        if len(selected) < minimum_selected:
            continue
        precision = float(selected["candidate_correct"].mean())
        if precision >= target_precision:
            candidates.append(
                {
                    "minimum_support_channels": threshold,
                    "n_selected": int(len(selected)),
                    "n_correct": int(selected["candidate_correct"].sum()),
                    "precision": precision,
                }
            )
    if not candidates:
        return None
    return sorted(candidates, key=lambda item: (-item["n_selected"], item["minimum_support_channels"]))[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    config = json.loads(args.config.read_text(encoding="utf-8"))
    frame = read_table(args.evidence)
    id_column = config.get("id_column", "cell_id")
    split_column = config.get("split_column", "calibration_split")
    truth_column = config.get("truth_column", "truth_label")
    route_column = config.get("route_column")
    heldout_value = str(config.get("heldout_value", "heldout"))
    query_value = str(config.get("query_value", "query"))
    channels = config.get("channels", [])
    if len(channels) < 2:
        raise SystemExit("at least two independent evidence channels are required")
    names = [channel.get("name", "") for channel in channels]
    if any(not name for name in names) or len(names) != len(set(names)):
        raise SystemExit("channel names must be nonempty and unique")

    required = {id_column, split_column, truth_column}
    if route_column:
        required.add(route_column)
    for channel in channels:
        required.add(channel["label_column"])
        if channel.get("accept_column"):
            required.add(channel["accept_column"])
    required.update(config.get("required_boolean_columns", []))
    required.update(config.get("required_heldout_boolean_columns", []))
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise SystemExit(f"evidence table lacks required columns: {missing}")
    if frame[id_column].eq("").any() or frame[id_column].duplicated().any():
        raise SystemExit("observation IDs must be nonempty and globally unique")
    allowed_splits = {heldout_value, query_value}
    unknown_splits = sorted(set(frame[split_column]).difference(allowed_splits))
    if unknown_splits:
        raise SystemExit(f"unknown calibration split values: {unknown_splits}")
    heldout_mask = frame[split_column].eq(heldout_value)
    query_mask = frame[split_column].eq(query_value)
    if not heldout_mask.any() or not query_mask.any():
        raise SystemExit("both held-out anchors and query observations are required")
    if frame.loc[heldout_mask, truth_column].eq("").any():
        raise SystemExit("held-out anchors require truth labels")
    for column in config.get("required_heldout_boolean_columns", []):
        if not frame.loc[heldout_mask, column].map(truth).all():
            raise SystemExit(f"held-out anchors fail required query-like audit column: {column}")

    for channel in channels:
        frame[f"__accept__{channel['name']}"] = channel_accepts(frame, channel)
    candidates = frame.apply(lambda row: choose_candidate(row, channels), axis=1)
    frame["candidate_label"] = [item[0] for item in candidates]
    frame["supporting_channels"] = [";".join(item[1]) for item in candidates]
    frame["support_channel_n"] = [len(item[1]) for item in candidates]
    independent = {channel["name"] for channel in channels if channel.get("role", "independent") != "primary"}
    frame["independent_support_n"] = [
        sum(name in independent for name in item[1]) for item in candidates
    ]

    required_channels = set(config.get("required_channels", []))
    unknown = required_channels.difference(names)
    if unknown:
        raise SystemExit(f"required_channels are not configured: {sorted(unknown)}")
    any_groups = [set(group) for group in config.get("any_of_channel_groups", [])]
    unknown_any = set().union(*any_groups).difference(names) if any_groups else set()
    if unknown_any:
        raise SystemExit(f"any_of_channel_groups contain unknown channels: {sorted(unknown_any)}")
    required_boolean_columns = config.get("required_boolean_columns", [])

    gate_values = []
    gate_reasons = []
    for _, row in frame.iterrows():
        support = set(filter(None, safe(row["supporting_channels"]).split(";")))
        reasons = []
        if not safe(row["candidate_label"]):
            reasons.append("no_consensus_label")
        missing_required = sorted(required_channels.difference(support))
        if missing_required:
            reasons.append("missing_required_channels:" + ",".join(missing_required))
        for number, group in enumerate(any_groups, 1):
            if not support.intersection(group):
                reasons.append(f"missing_any_channel_group_{number}:" + ",".join(sorted(group)))
        failed_booleans = [column for column in required_boolean_columns if not truth(row[column])]
        if failed_booleans:
            reasons.append("failed_boolean_gates:" + ",".join(failed_booleans))
        gate_values.append(not reasons)
        gate_reasons.append("pass" if not reasons else ";".join(reasons))
    frame["joint_gate_pass"] = gate_values
    frame["joint_gate_reason"] = gate_reasons
    frame["candidate_correct"] = frame["candidate_label"].eq(frame[truth_column])

    group_input = list(config.get("calibration_group_columns", []))
    for column in group_input:
        if column not in frame:
            raise SystemExit(f"calibration group column is missing: {column}")
    group_columns = group_input + ["candidate_label"]
    moderate_target = float(config.get("moderate_target_precision", 0.90))
    high_target = float(config.get("high_target_precision", 0.95))
    moderate_min_votes = int(config.get("moderate_min_support_channels", 2))
    high_min_votes = int(config.get("high_min_support_channels", 3))
    moderate_min_selected = int(config.get("moderate_min_heldout_selected", 20))
    high_min_selected = int(config.get("high_min_heldout_selected", 20))
    if not (0 < moderate_target <= high_target <= 1):
        raise SystemExit("precision targets must satisfy 0 < moderate <= high <= 1")
    if not (1 <= moderate_min_votes <= high_min_votes <= len(channels)):
        raise SystemExit("support-channel minima must satisfy 1 <= moderate <= high <= n_channels")

    threshold_rows = []
    heldout = frame.loc[heldout_mask].copy()
    for key, group in heldout.groupby(group_columns, dropna=False, sort=True):
        key = key if isinstance(key, tuple) else (key,)
        labels = dict(zip(group_columns, key))
        if not labels["candidate_label"]:
            continue
        moderate = calibrate_threshold(group, moderate_min_votes, moderate_target, moderate_min_selected)
        high = calibrate_threshold(group, high_min_votes, high_target, high_min_selected)
        if high and (not moderate or moderate["minimum_support_channels"] > high["minimum_support_channels"]):
            # A high-precision cumulative set always qualifies for the lower
            # moderate target.  Reuse it so high can never exist outside the
            # moderate-or-higher gate.
            high_selected = group[
                group["joint_gate_pass"]
                & group["support_channel_n"].ge(high["minimum_support_channels"])
            ]
            if len(high_selected) >= moderate_min_selected and float(high_selected["candidate_correct"].mean()) >= moderate_target:
                moderate = {
                    "minimum_support_channels": high["minimum_support_channels"],
                    "n_selected": int(len(high_selected)),
                    "n_correct": int(high_selected["candidate_correct"].sum()),
                    "precision": float(high_selected["candidate_correct"].mean()),
                }
        row = dict(labels)
        for tier, result in (("moderate_or_higher", moderate), ("high", high)):
            row[f"{tier}_minimum_support_channels"] = "" if result is None else result["minimum_support_channels"]
            row[f"{tier}_n_selected"] = 0 if result is None else result["n_selected"]
            row[f"{tier}_n_correct"] = 0 if result is None else result["n_correct"]
            row[f"{tier}_precision"] = "" if result is None else result["precision"]
        threshold_rows.append(row)
    thresholds = pd.DataFrame(threshold_rows)
    if thresholds.empty:
        raise SystemExit("no calibration group reached the held-out precision/support contract")

    threshold_key = group_columns
    query = frame.loc[query_mask].copy()
    query = query.merge(thresholds, on=threshold_key, how="left", validate="many_to_one")
    moderate_threshold = pd.to_numeric(
        query["moderate_or_higher_minimum_support_channels"], errors="coerce"
    )
    high_threshold = pd.to_numeric(query["high_minimum_support_channels"], errors="coerce")
    query["meets_moderate_or_higher"] = (
        query["joint_gate_pass"]
        & moderate_threshold.notna()
        & query["support_channel_n"].ge(moderate_threshold)
    )
    query["meets_high"] = (
        query["joint_gate_pass"]
        & high_threshold.notna()
        & query["support_channel_n"].ge(high_threshold)
    )
    if (query["meets_high"] & ~query["meets_moderate_or_higher"]).any():
        raise SystemExit("nested-tier invariant failed: a high row misses moderate-or-higher")
    query["consensus_tier"] = "low-reject"
    query.loc[query["meets_moderate_or_higher"], "consensus_tier"] = "moderate-only"
    query.loc[query["meets_high"], "consensus_tier"] = "high"
    query["writeback_eligible"] = query["meets_moderate_or_higher"]
    query["proposed_state"] = config.get("reject_state", "qc_holdout")
    query["proposed_broad_label"] = ""
    eligible = query["writeback_eligible"]
    query.loc[eligible, "proposed_state"] = "defined_broad_only"
    query.loc[eligible, "proposed_broad_label"] = query.loc[eligible, "candidate_label"]
    query["proposed_fine_label"] = ""
    query["fine_anchor_eligible"] = False
    query["writeback_status"] = "proposal_only_requires_atomic_commit"

    internal_columns = [column for column in frame if column.startswith("__accept__")]
    heldout_out = frame.loc[heldout_mask].drop(columns=internal_columns)
    query_out = query.drop(columns=internal_columns)
    heldout_out.to_csv(args.out / "heldout_multichannel_evidence.tsv", sep="\t", index=False)
    thresholds.to_csv(args.out / "multichannel_thresholds.tsv", sep="\t", index=False)
    query_out.to_csv(
        args.out / "calibrated_multichannel_broad_rescue.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    tier_counts = query["consensus_tier"].value_counts().to_dict()
    label_counts = (
        query.loc[eligible, "proposed_broad_label"].value_counts().rename_axis("broad_label").reset_index(name="n")
    )
    label_counts.to_csv(args.out / "eligible_broad_label_counts.tsv", sep="\t", index=False)
    manifest = {
        "status": "CALIBRATED_MULTICHANNEL_PROPOSAL_ONLY",
        "n_heldout": int(heldout_mask.sum()),
        "n_query": int(query_mask.sum()),
        "n_channels": len(channels),
        "channels": names,
        "calibration_group_columns": group_columns,
        "moderate_target_precision": moderate_target,
        "high_target_precision": high_target,
        "query_high": int(tier_counts.get("high", 0)),
        "query_moderate_only": int(tier_counts.get("moderate-only", 0)),
        "query_moderate_or_higher": int(query["meets_moderate_or_higher"].sum()),
        "query_low_reject": int(tier_counts.get("low-reject", 0)),
        "nested_tier_invariant": bool((~query["meets_high"] | query["meets_moderate_or_higher"]).all()),
        "broad_only": True,
        "fine_anchor_eligible": False,
        "ledger_writeback_performed": False,
        "warning": "This artifact is a calibrated broad-only proposal. Commit only after biological review and the framework's atomic state checks.",
    }
    (args.out / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
