#!/usr/bin/env python3
"""Fail-closed validation for route-bounded, calibrated mapping evidence."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "route"


def read_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping-root", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--bundle-root", type=Path)
    parser.add_argument("--route-col", required=True)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--include-routes", default="")
    parser.add_argument("--state-col", default="")
    parser.add_argument("--include-states", default="")
    parser.add_argument("--minimum-calibration-support", type=int, default=20)
    parser.add_argument("--minimum-target-precision", type=float, default=0.9)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    membership = pd.read_csv(
        args.membership, sep="\t", dtype={args.cell_id_col: str}, low_memory=False
    )
    errors: list[str] = []
    for column in (args.cell_id_col, args.route_col):
        if column not in membership:
            errors.append(f"membership lacks {column}")
    if errors:
        raise SystemExit("; ".join(errors))
    if not membership[args.cell_id_col].is_unique:
        errors.append("membership cell IDs are duplicated")

    requested_routes = parse_csv(args.include_routes)
    if requested_routes:
        membership = membership[membership[args.route_col].astype(str).isin(requested_routes)].copy()
    if args.include_states:
        if not args.state_col or args.state_col not in membership:
            errors.append("state filtering requested but state column is absent")
        else:
            membership = membership[
                membership[args.state_col].astype(str).isin(parse_csv(args.include_states))
            ].copy()
    routes = requested_routes or membership[args.route_col].astype(str).drop_duplicates().tolist()
    missing_routes = sorted(set(routes).difference(set(membership[args.route_col].astype(str))))
    if missing_routes:
        errors.append(f"requested routes absent from membership: {missing_routes}")

    route_rows: list[dict] = []
    label_rows: list[pd.DataFrame] = []
    threshold_rows: list[pd.DataFrame] = []
    confusion_rows: list[pd.DataFrame] = []
    observed_ids: list[str] = []
    for route in routes:
        route_dir = args.mapping_root / safe_name(route)
        expected_ids = membership.loc[
            membership[args.route_col].astype(str).eq(route), args.cell_id_col
        ].astype(str).tolist()
        try:
            query_manifest = read_json(route_dir / "query_mapping" / "mapping_manifest.json")
            heldout_manifest = read_json(route_dir / "heldout_mapping" / "mapping_manifest.json")
            calibration_summary = read_json(route_dir / "calibration" / "calibration_summary.json")
            mapping = pd.read_csv(
                route_dir / "calibration" / "calibrated_query_mapping.tsv.gz",
                sep="\t",
                dtype={args.cell_id_col: str},
            )
            thresholds = pd.read_csv(
                route_dir / "calibration" / "mapping_thresholds.tsv", sep="\t"
            )
        except (FileNotFoundError, OSError, ValueError) as exc:
            errors.append(f"{route}: cannot load mapping artifacts: {exc}")
            continue

        route_errors: list[str] = []
        required_mapping = {
            args.cell_id_col,
            "predicted_label",
            "confidence",
            "margin",
            "mapping_status",
        }
        missing_mapping_cols = sorted(required_mapping.difference(mapping.columns))
        if missing_mapping_cols:
            route_errors.append(f"mapping lacks columns {missing_mapping_cols}")
        if query_manifest.get("status") != "EVIDENCE_ONLY":
            route_errors.append("query manifest is not EVIDENCE_ONLY")
        if heldout_manifest.get("status") != "EVIDENCE_ONLY":
            route_errors.append("held-out manifest is not EVIDENCE_ONLY")
        if calibration_summary.get("status") != "CALIBRATED_EVIDENCE_ONLY":
            route_errors.append("calibration manifest is not CALIBRATED_EVIDENCE_ONLY")
        if query_manifest.get("membership_exact_order") is not True:
            route_errors.append("query membership exact-order contract failed")
        if query_manifest.get("prediction_rows_equal_membership") is not True:
            route_errors.append("prediction/membership row-count contract failed")
        if heldout_manifest.get("prediction_rows_equal_membership") is not True:
            route_errors.append("held-out prediction row-count contract failed")
        if int(query_manifest.get("n_query", -1)) != len(expected_ids):
            route_errors.append("query manifest row count differs from frozen route")
        if int(query_manifest.get("membership_n", -1)) != len(expected_ids):
            route_errors.append("query manifest membership count differs from frozen route")
        if len(mapping) != len(expected_ids):
            route_errors.append("calibrated mapping row count differs from frozen route")
        if args.cell_id_col in mapping:
            ids = mapping[args.cell_id_col].astype(str).tolist()
            if ids != expected_ids:
                route_errors.append("calibrated mapping IDs/order differ from frozen route")
            if not mapping[args.cell_id_col].is_unique:
                route_errors.append("calibrated mapping IDs are duplicated")
            observed_ids.extend(ids)

        target_precision = float(calibration_summary.get("target_precision", np.nan))
        if not np.isfinite(target_precision) or target_precision < args.minimum_target_precision:
            route_errors.append(
                f"target precision {target_precision} is below {args.minimum_target_precision}"
            )
        required_thresholds = {
            "predicted_label",
            "confidence_threshold",
            "margin_threshold",
            "calibration_precision",
            "calibration_support",
        }
        missing_threshold_cols = sorted(required_thresholds.difference(thresholds.columns))
        if missing_threshold_cols:
            route_errors.append(f"threshold table lacks columns {missing_threshold_cols}")
        elif len(thresholds):
            if (thresholds["calibration_precision"] < target_precision - 1e-12).any():
                route_errors.append("one or more label thresholds miss target precision")
            if (thresholds["calibration_support"] < args.minimum_calibration_support).any():
                route_errors.append("one or more label thresholds miss minimum support")
        if len(thresholds):
            threshold_copy = thresholds.copy()
            threshold_copy.insert(0, "route", route)
            threshold_rows.append(threshold_copy)

        if args.bundle_root:
            try:
                heldout_predictions = pd.read_csv(
                    route_dir / "heldout_mapping" / "anchor_knn_predictions.tsv.gz",
                    sep="\t",
                    dtype={args.cell_id_col: str},
                )
                heldout_truth = pd.read_csv(
                    args.bundle_root / safe_name(route) / "heldout_truth.tsv",
                    sep="\t",
                    dtype={args.cell_id_col: str},
                )
                if not heldout_predictions[args.cell_id_col].is_unique:
                    route_errors.append("held-out prediction IDs are duplicated")
                if not heldout_truth[args.cell_id_col].is_unique:
                    route_errors.append("held-out truth IDs are duplicated")
                heldout = heldout_predictions.merge(
                    heldout_truth, on=args.cell_id_col, how="inner", validate="one_to_one"
                )
                if len(heldout) != len(heldout_predictions) or len(heldout) != len(heldout_truth):
                    route_errors.append("held-out prediction/truth membership mismatch")
                confusion = (
                    heldout.groupby(["true_label", "predicted_label"], dropna=False)
                    .size()
                    .rename("n")
                    .reset_index()
                )
                confusion.insert(0, "route", route)
                confusion_rows.append(confusion)
            except (FileNotFoundError, OSError, ValueError, KeyError) as exc:
                route_errors.append(f"cannot validate held-out truth: {exc}")

        allowed_status = {
            "calibrated_medium_high_broad_only",
            "rejected_to_qc_or_review",
        }
        if "mapping_status" in mapping:
            unknown_status = sorted(set(mapping["mapping_status"].astype(str)).difference(allowed_status))
            if unknown_status:
                route_errors.append(f"unknown mapping statuses: {unknown_status}")
            accepted = mapping["mapping_status"].eq("calibrated_medium_high_broad_only")
            if accepted.any():
                check = mapping.loc[accepted]
                invalid = ~(
                    check["confidence"].ge(check["confidence_threshold"])
                    & check["margin"].ge(check["margin_threshold"])
                )
                if invalid.any():
                    route_errors.append("accepted mappings violate calibrated thresholds")
            accepted_n = int(accepted.sum())
        else:
            accepted_n = 0

        route_rows.append(
            {
                "route": route,
                "expected_n": len(expected_ids),
                "mapping_n": len(mapping),
                "accepted_n": accepted_n,
                "rejected_n": len(mapping) - accepted_n,
                "target_precision": target_precision,
                "shared_genes": query_manifest.get("shared_genes"),
                "labels_with_thresholds": len(thresholds),
                "route_status": "PASS" if not route_errors else "FAIL",
                "route_errors": "; ".join(route_errors),
            }
        )
        if route_errors:
            errors.extend(f"{route}: {message}" for message in route_errors)
        if required_mapping.issubset(mapping.columns):
            counts = (
                mapping.groupby(["predicted_label", "mapping_status"], dropna=False)
                .size()
                .rename("n")
                .reset_index()
            )
            counts.insert(0, "route", route)
            label_rows.append(counts)

    expected_union = membership[args.cell_id_col].astype(str).tolist()
    if len(observed_ids) != len(expected_union):
        errors.append("mapping union row count differs from expected membership union")
    if len(observed_ids) != len(set(observed_ids)):
        errors.append("mapping routes overlap or contain duplicated IDs")
    if set(observed_ids) != set(expected_union):
        errors.append("mapping union does not equal expected membership union")

    route_summary = pd.DataFrame(route_rows)
    route_summary.to_csv(args.out / "route_validation_summary.tsv", sep="\t", index=False)
    if label_rows:
        pd.concat(label_rows, ignore_index=True).to_csv(
            args.out / "route_label_status_counts.tsv", sep="\t", index=False
        )
    else:
        pd.DataFrame(columns=["route", "predicted_label", "mapping_status", "n"]).to_csv(
            args.out / "route_label_status_counts.tsv", sep="\t", index=False
        )
    if threshold_rows:
        pd.concat(threshold_rows, ignore_index=True).to_csv(
            args.out / "route_mapping_thresholds.tsv", sep="\t", index=False
        )
    else:
        pd.DataFrame(
            columns=[
                "route",
                "predicted_label",
                "confidence_threshold",
                "margin_threshold",
                "calibration_precision",
                "calibration_coverage",
                "calibration_support",
            ]
        ).to_csv(args.out / "route_mapping_thresholds.tsv", sep="\t", index=False)
    if confusion_rows:
        pd.concat(confusion_rows, ignore_index=True).to_csv(
            args.out / "heldout_confusion_counts.tsv", sep="\t", index=False
        )
    else:
        pd.DataFrame(columns=["route", "true_label", "predicted_label", "n"]).to_csv(
            args.out / "heldout_confusion_counts.tsv", sep="\t", index=False
        )

    manifest = {
        "status": "PASS" if not errors else "FAIL",
        "mapping_root": str(args.mapping_root.resolve()),
        "membership": str(args.membership.resolve()),
        "bundle_root": str(args.bundle_root.resolve()) if args.bundle_root else None,
        "routes": routes,
        "n_routes": len(routes),
        "expected_union_n": len(expected_union),
        "observed_union_n": len(observed_ids),
        "observed_unique_n": len(set(observed_ids)),
        "exact_union": set(observed_ids) == set(expected_union),
        "routes_disjoint": len(observed_ids) == len(set(observed_ids)),
        "fine_anchor_eligible": False,
        "errors": errors,
        "warning": "PASS validates mapping evidence boundaries only; it never authorizes biological writeback without independent marker and spatial adjudication.",
    }
    (args.out / "validation_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
