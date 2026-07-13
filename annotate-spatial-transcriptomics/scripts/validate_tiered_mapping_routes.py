#!/usr/bin/env python3
"""Fail-closed validation for mutually exclusive route-level mapping tiers."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "route"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tiered-root", required=True, type=Path)
    parser.add_argument("--source-mapping-root", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--route-col", required=True)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--include-routes", default="")
    parser.add_argument("--state-col", default="")
    parser.add_argument("--include-states", default="")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    membership = pd.read_csv(
        args.membership, sep="\t", dtype={args.cell_id_col: str}, low_memory=False
    )
    errors: list[str] = []
    for column in (args.cell_id_col, args.route_col):
        if column not in membership:
            raise SystemExit(f"membership lacks {column}")
    if not membership[args.cell_id_col].is_unique:
        errors.append("membership IDs are duplicated")
    routes = parse_csv(args.include_routes) or membership[args.route_col].astype(str).drop_duplicates().tolist()
    membership = membership[membership[args.route_col].astype(str).isin(routes)].copy()
    if args.include_states:
        if not args.state_col or args.state_col not in membership:
            raise SystemExit("state filtering requested but state column is absent")
        membership = membership[
            membership[args.state_col].astype(str).isin(parse_csv(args.include_states))
        ].copy()

    route_rows: list[dict] = []
    count_rows: list[pd.DataFrame] = []
    threshold_rows: list[pd.DataFrame] = []
    validation_rows: list[pd.DataFrame] = []
    cumulative_validation_rows: list[pd.DataFrame] = []
    observed_ids: list[str] = []
    for route in routes:
        expected_ids = membership.loc[
            membership[args.route_col].astype(str).eq(route), args.cell_id_col
        ].astype(str).tolist()
        route_dir = args.tiered_root / safe_name(route)
        source_dir = args.source_mapping_root / safe_name(route)
        route_errors: list[str] = []
        try:
            manifest = json.loads((route_dir / "calibration_summary.json").read_text())
            source_manifest = json.loads(
                (source_dir / "query_mapping" / "mapping_manifest.json").read_text()
            )
            mapping = pd.read_csv(
                route_dir / "calibrated_query_mapping.tsv.gz",
                sep="\t",
                dtype={args.cell_id_col: str},
            )
            thresholds = pd.read_csv(route_dir / "mapping_thresholds.tsv", sep="\t")
            heldout_validation = pd.read_csv(
                route_dir / "heldout_tier_validation.tsv", sep="\t"
            )
            cumulative_validation = pd.read_csv(
                route_dir / "heldout_cumulative_validation.tsv", sep="\t"
            )
        except (FileNotFoundError, OSError, ValueError) as exc:
            errors.append(f"{route}: cannot load tiered artifacts: {exc}")
            continue

        if manifest.get("status") != "CALIBRATED_TIERED_EVIDENCE_ONLY":
            route_errors.append("tiered calibration manifest status is invalid")
        if source_manifest.get("membership_exact_order") is not True:
            route_errors.append("source query membership order contract failed")
        if source_manifest.get("prediction_rows_equal_membership") is not True:
            route_errors.append("source mapping row-count contract failed")
        required = {
            args.cell_id_col,
            "predicted_label",
            "confidence",
            "margin",
            "mapping_tier",
            "mapping_status",
            "fine_anchor_eligible",
            "meets_moderate_or_higher",
        }
        missing = sorted(required.difference(mapping.columns))
        if missing:
            route_errors.append(f"mapping lacks columns {missing}")
        if len(mapping) != len(expected_ids):
            route_errors.append("mapping row count differs from frozen route")
        elif mapping[args.cell_id_col].astype(str).tolist() != expected_ids:
            route_errors.append("mapping IDs/order differ from frozen route")
        if not mapping[args.cell_id_col].is_unique:
            route_errors.append("mapping IDs are duplicated")
        observed_ids.extend(mapping[args.cell_id_col].astype(str).tolist())

        allowed_tiers = {"high", "moderate", "low_reject"}
        unknown_tiers = sorted(set(mapping["mapping_tier"].astype(str)).difference(allowed_tiers))
        if unknown_tiers:
            route_errors.append(f"unknown mapping tiers: {unknown_tiers}")
        expected_status = {
            "high": "calibrated_high_broad_candidate",
            "moderate": "calibrated_moderate_broad_candidate",
            "low_reject": "rejected_to_qc_or_review",
        }
        bad_status = mapping["mapping_status"].astype(str).ne(
            mapping["mapping_tier"].astype(str).map(expected_status)
        )
        if bad_status.any():
            route_errors.append("mapping tier/status contract failed")
        if mapping["fine_anchor_eligible"].astype(str).str.lower().isin(["true", "1"]).any():
            route_errors.append("mapping candidates were incorrectly made fine anchors")
        expected_moderate_gate = mapping["mapping_tier"].astype(str).isin(
            ["high", "moderate"]
        )
        observed_moderate_gate = mapping["meets_moderate_or_higher"].astype(str).str.lower().isin(
            ["true", "1"]
        )
        if not expected_moderate_gate.equals(observed_moderate_gate):
            route_errors.append("moderate-or-higher cumulative gate disagrees with tier labels")
        manifest_total = sum(
            int(manifest.get(key, -1))
            for key in ("query_high", "query_moderate", "query_low_reject")
        )
        if manifest_total != len(mapping) or int(manifest.get("query_n", -1)) != len(mapping):
            route_errors.append("manifest tier counts do not partition query")

        required_threshold = {
            "predicted_label",
            "tier",
            "target_precision",
            "minimum_support",
            "calibration_precision",
            "calibration_support",
        }
        if not required_threshold.issubset(thresholds.columns):
            route_errors.append("threshold table lacks required columns")
        elif len(thresholds):
            if (thresholds["calibration_precision"] < thresholds["target_precision"] - 1e-12).any():
                route_errors.append("one or more thresholds miss their target precision")
            if (thresholds["calibration_support"] < thresholds["minimum_support"]).any():
                route_errors.append("one or more thresholds miss their minimum support")
            copy = thresholds.copy()
            copy.insert(0, "route", route)
            threshold_rows.append(copy)

        for tier, target_key, support_key in (
            ("high", "high_target_precision", "high_min_support"),
            (
                "moderate_or_higher",
                "moderate_target_precision",
                "moderate_min_support",
            ),
        ):
            selected = cumulative_validation[
                cumulative_validation["cumulative_tier"].astype(str).eq(tier)
            ]
            if len(selected):
                if (selected["validation_precision"] < float(manifest[target_key]) - 1e-12).any():
                    route_errors.append(f"held-out cumulative {tier} tier misses target precision")
                if (selected["validation_n"] < int(manifest[support_key])).any():
                    route_errors.append(f"held-out cumulative {tier} tier misses minimum support")
        if int(manifest.get("query_moderate_or_higher", -1)) != int(
            expected_moderate_gate.sum()
        ):
            route_errors.append("manifest moderate-or-higher count disagrees with mapping")
        if int(manifest.get("query_moderate_only", -1)) != int(
            mapping["mapping_tier"].astype(str).eq("moderate").sum()
        ):
            route_errors.append("manifest moderate-only count disagrees with mapping")
        validation_copy = heldout_validation.copy()
        validation_copy.insert(0, "route", route)
        validation_rows.append(validation_copy)
        cumulative_copy = cumulative_validation.copy()
        cumulative_copy.insert(0, "route", route)
        cumulative_validation_rows.append(cumulative_copy)

        counts = (
            mapping.groupby(["predicted_label", "mapping_tier", "mapping_status"], dropna=False)
            .size()
            .rename("n")
            .reset_index()
        )
        counts.insert(0, "route", route)
        count_rows.append(counts)
        tier_counts = mapping["mapping_tier"].value_counts()
        route_rows.append(
            {
                "route": route,
                "expected_n": len(expected_ids),
                "mapping_n": len(mapping),
                "high_n": int(tier_counts.get("high", 0)),
                "moderate_n": int(tier_counts.get("moderate", 0)),
                "moderate_only_n": int(tier_counts.get("moderate", 0)),
                "moderate_or_higher_n": int(
                    tier_counts.get("high", 0) + tier_counts.get("moderate", 0)
                ),
                "low_reject_n": int(tier_counts.get("low_reject", 0)),
                "labels_with_high_threshold": manifest.get("labels_with_high_threshold"),
                "labels_with_moderate_threshold": manifest.get("labels_with_moderate_threshold"),
                "route_status": "PASS" if not route_errors else "FAIL",
                "route_errors": "; ".join(route_errors),
            }
        )
        errors.extend(f"{route}: {message}" for message in route_errors)

    expected_ids = membership[args.cell_id_col].astype(str).tolist()
    if len(observed_ids) != len(expected_ids):
        errors.append("tiered mapping union row count differs from frozen membership")
    if len(observed_ids) != len(set(observed_ids)):
        errors.append("tiered mapping routes overlap or contain duplicate IDs")
    if set(observed_ids) != set(expected_ids):
        errors.append("tiered mapping union differs from frozen membership")

    pd.DataFrame(route_rows).to_csv(args.out / "route_tier_summary.tsv", sep="\t", index=False)
    for frames, name, columns in (
        (count_rows, "route_label_tier_counts.tsv", ["route", "predicted_label", "mapping_tier", "mapping_status", "n"]),
        (threshold_rows, "route_tier_thresholds.tsv", ["route", "predicted_label", "tier"]),
        (validation_rows, "heldout_tier_validation.tsv", ["route", "predicted_label", "mapping_tier", "validation_n", "validation_precision"]),
        (cumulative_validation_rows, "heldout_cumulative_validation.tsv", ["route", "predicted_label", "cumulative_tier", "validation_n", "validation_precision", "target_precision"]),
    ):
        output = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns)
        output.to_csv(args.out / name, sep="\t", index=False)
    manifest = {
        "status": "PASS" if not errors else "FAIL",
        "tiered_root": str(args.tiered_root.resolve()),
        "source_mapping_root": str(args.source_mapping_root.resolve()),
        "membership": str(args.membership.resolve()),
        "n_routes": len(routes),
        "expected_union_n": len(expected_ids),
        "observed_union_n": len(observed_ids),
        "observed_unique_n": len(set(observed_ids)),
        "exact_union": set(observed_ids) == set(expected_ids),
        "routes_disjoint": len(observed_ids) == len(set(observed_ids)),
        "fine_anchor_eligible": False,
        "errors": errors,
        "warning": "PASS validates nested cumulative high and moderate-or-higher broad candidates. The mutually exclusive moderate tier means moderate-only; independent biological and spatial adjudication is still mandatory.",
    }
    (args.out / "validation_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
