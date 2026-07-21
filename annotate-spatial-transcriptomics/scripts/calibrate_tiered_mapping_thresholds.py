#!/usr/bin/env python3
"""Calibrate mutually exclusive high/moderate broad-only mapping tiers."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def hash_ids(values: pd.Series) -> str:
    payload = "\n".join(sorted(values.astype(str).tolist())) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def select_threshold(
    frame: pd.DataFrame,
    truth_col: str,
    target_precision: float,
    minimum_support: int,
    required_indices: pd.Index | None = None,
) -> dict | None:
    if frame.empty:
        return None
    confidence_grid = np.unique(
        np.quantile(frame["confidence"], np.linspace(0, 1, min(41, len(frame))))
    )
    margin_grid = np.unique(
        np.quantile(frame["margin"], np.linspace(0, 1, min(41, len(frame))))
    )
    best = None
    for confidence in confidence_grid:
        for margin in margin_grid:
            selected = frame[
                frame["confidence"].ge(confidence) & frame["margin"].ge(margin)
            ]
            if required_indices is not None and not required_indices.isin(selected.index).all():
                continue
            support = len(selected)
            if support < minimum_support:
                continue
            precision = float(
                selected["predicted_label"].eq(selected[truth_col]).mean()
            )
            if precision < target_precision:
                continue
            candidate = {
                "confidence_threshold": float(confidence),
                "margin_threshold": float(margin),
                "calibration_precision": precision,
                "calibration_support": support,
                "calibration_coverage": support / len(frame),
            }
            if best is None or (
                candidate["calibration_support"],
                candidate["calibration_precision"],
                -candidate["confidence_threshold"],
                -candidate["margin_threshold"],
            ) > (
                best["calibration_support"],
                best["calibration_precision"],
                -best["confidence_threshold"],
                -best["margin_threshold"],
            ):
                best = candidate
    return best


def passes(frame: pd.DataFrame, threshold: dict | None) -> pd.Series:
    if threshold is None:
        return pd.Series(False, index=frame.index)
    return frame["confidence"].ge(threshold["confidence_threshold"]) & frame[
        "margin"
    ].ge(threshold["margin_threshold"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--query-predictions", required=True, type=Path)
    parser.add_argument("--calibration-origin-manifest", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--truth-col", default="true_label")
    parser.add_argument("--high-target-precision", type=float, default=0.95)
    parser.add_argument("--moderate-target-precision", type=float, default=0.90)
    parser.add_argument("--high-min-support", type=int, default=20)
    parser.add_argument("--moderate-min-support", type=int, default=20)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.moderate_target_precision > args.high_target_precision:
        raise SystemExit("moderate target precision cannot exceed high target precision")

    origin = json.loads(args.calibration_origin_manifest.read_text(encoding="utf-8"))
    if origin.get("status") != "PASS":
        raise SystemExit("calibration origin manifest is not PASS")
    if origin.get("heldout_origin") != "query_like_heldout_current_query_anchors":
        raise SystemExit("final rescue calibration requires query-like held-out current-query anchors")
    if origin.get("reference_self_classification") is not False:
        raise SystemExit("reference self-classification cannot calibrate query rescue")
    if int(origin.get("anchor_target_overlap", -1)) != 0:
        raise SystemExit("held-out anchors overlap the target query")

    predictions = pd.read_csv(
        args.predictions, sep="\t", dtype={args.cell_id_col: str}
    )
    truth = pd.read_csv(args.truth, sep="\t", dtype={args.cell_id_col: str})
    query = pd.read_csv(
        args.query_predictions, sep="\t", dtype={args.cell_id_col: str}
    )
    if int(origin.get("n_anchors", -1)) != len(predictions) or hash_ids(predictions[args.cell_id_col]) != origin.get("anchor_ids_sha256"):
        raise SystemExit("held-out prediction membership differs from the origin manifest")
    if int(origin.get("n_target", -1)) != len(query) or hash_ids(query[args.cell_id_col]) != origin.get("target_ids_sha256"):
        raise SystemExit("query prediction membership differs from the origin manifest")
    required_prediction = {
        args.cell_id_col,
        "predicted_label",
        "confidence",
        "margin",
    }
    if not required_prediction.issubset(predictions.columns) or not required_prediction.issubset(
        query.columns
    ):
        raise SystemExit("prediction table lacks required columns")
    if args.truth_col not in truth or args.cell_id_col not in truth:
        raise SystemExit("truth table lacks required columns")
    if not predictions[args.cell_id_col].is_unique or not truth[args.cell_id_col].is_unique:
        raise SystemExit("held-out prediction/truth IDs are duplicated")
    heldout = predictions.merge(
        truth, on=args.cell_id_col, how="inner", validate="one_to_one"
    )
    if len(heldout) != len(predictions) or len(heldout) != len(truth):
        raise SystemExit("held-out prediction/truth membership mismatch")

    threshold_rows: list[dict] = []
    heldout["mapping_tier"] = "low_reject"
    for label in sorted(heldout["predicted_label"].astype(str).unique()):
        group = heldout[heldout["predicted_label"].astype(str).eq(label)].copy()
        high = select_threshold(
            group,
            args.truth_col,
            args.high_target_precision,
            args.high_min_support,
        )
        high_pass = passes(group, high)
        if high is not None:
            threshold_rows.append(
                {
                    "predicted_label": label,
                    "tier": "high",
                    "target_precision": args.high_target_precision,
                    "minimum_support": args.high_min_support,
                    **high,
                }
            )
            heldout.loc[group.index[high_pass], "mapping_tier"] = "high"

        # Calibrate the moderate threshold on the complete predicted-label group,
        # matching the nested high/moderate semantics used by the completed
        # Scanpy and R workflows.  The moderate selection must contain every
        # high row; query output remains mutually exclusive by assigning high
        # first and only the additional rows to the moderate tier.
        moderate = select_threshold(
            group,
            args.truth_col,
            args.moderate_target_precision,
            args.moderate_min_support,
            required_indices=group.index[high_pass],
        )
        moderate_pass = passes(group, moderate)
        if moderate is not None:
            moderate_only = moderate_pass & ~high_pass
            exclusive_precision = (
                float(
                    group.loc[moderate_only, "predicted_label"]
                    .eq(group.loc[moderate_only, args.truth_col])
                    .mean()
                )
                if moderate_only.any()
                else None
            )
            threshold_rows.append(
                {
                    "predicted_label": label,
                    "tier": "moderate",
                    "target_precision": args.moderate_target_precision,
                    "minimum_support": args.moderate_min_support,
                    **moderate,
                    "exclusive_support": int(moderate_only.sum()),
                    "exclusive_precision": exclusive_precision,
                }
            )
            heldout.loc[group.index[moderate_only], "mapping_tier"] = "moderate_only"

    threshold_columns = [
        "predicted_label",
        "tier",
        "target_precision",
        "minimum_support",
        "confidence_threshold",
        "margin_threshold",
        "calibration_precision",
        "calibration_support",
        "calibration_coverage",
        "exclusive_support",
        "exclusive_precision",
    ]
    thresholds = pd.DataFrame(threshold_rows, columns=threshold_columns)
    thresholds_path = args.out / "mapping_thresholds.tsv"
    thresholds.to_csv(thresholds_path, sep="\t", index=False)

    query["mapping_tier"] = "low_reject"
    for row in threshold_rows:
        label_mask = query["predicted_label"].astype(str).eq(row["predicted_label"])
        tier_mask = (
            label_mask
            & query["confidence"].ge(row["confidence_threshold"])
            & query["margin"].ge(row["margin_threshold"])
        )
        if row["tier"] == "high":
            query.loc[tier_mask, "mapping_tier"] = "high"
        else:
            query.loc[tier_mask & query["mapping_tier"].eq("low_reject"), "mapping_tier"] = (
                "moderate_only"
            )
    query["meets_moderate_or_higher"] = query["mapping_tier"].isin(
        ["high", "moderate_only"]
    )
    status_map = {
        "high": "calibrated_high_broad_candidate",
        "moderate_only": "calibrated_moderate_broad_candidate",
        "low_reject": "rejected_to_qc_or_review",
    }
    query["mapping_status"] = query["mapping_tier"].map(status_map)
    query["fine_anchor_eligible"] = False
    query_path = args.out / "calibrated_query_mapping.tsv.gz"
    query.to_csv(
        query_path,
        sep="\t",
        index=False,
        compression="gzip",
    )
    heldout["mapping_status"] = heldout["mapping_tier"].map(status_map)
    heldout["meets_moderate_or_higher"] = heldout["mapping_tier"].isin(
        ["high", "moderate_only"]
    )
    heldout_path = args.out / "heldout_tier_assignments.tsv.gz"
    heldout.to_csv(
        heldout_path,
        sep="\t",
        index=False,
        compression="gzip",
    )
    validation = (
        heldout.groupby(["predicted_label", "mapping_tier"], dropna=False)
        .agg(
            validation_n=(args.cell_id_col, "size"),
            validation_precision=(
                args.truth_col,
                lambda values: float(
                    np.mean(
                        heldout.loc[values.index, "predicted_label"].to_numpy()
                        == values.to_numpy()
                    )
                ),
            ),
        )
        .reset_index()
    )
    validation_path = args.out / "heldout_tier_validation.tsv"
    validation.to_csv(validation_path, sep="\t", index=False)

    cumulative_rows: list[dict] = []
    for label in sorted(heldout["predicted_label"].astype(str).unique()):
        group = heldout[heldout["predicted_label"].astype(str).eq(label)]
        for tier, mask, target in (
            ("high", group["mapping_tier"].eq("high"), args.high_target_precision),
            (
                "moderate_or_higher",
                group["mapping_tier"].isin(["high", "moderate_only"]),
                args.moderate_target_precision,
            ),
        ):
            selected = group.loc[mask]
            if selected.empty:
                continue
            cumulative_rows.append(
                {
                    "predicted_label": label,
                    "cumulative_tier": tier,
                    "validation_n": int(len(selected)),
                    "validation_precision": float(
                        selected["predicted_label"].eq(selected[args.truth_col]).mean()
                    ),
                    "target_precision": target,
                }
            )
    cumulative_path = args.out / "heldout_cumulative_validation.tsv"
    pd.DataFrame(cumulative_rows).to_csv(cumulative_path, sep="\t", index=False)

    counts = query["mapping_tier"].value_counts().to_dict()
    manifest = {
        "schema_version": "2.0",
        "status": "CALIBRATED_TIERED_EVIDENCE_ONLY",
        "high_target_precision": args.high_target_precision,
        "moderate_target_precision": args.moderate_target_precision,
        "high_min_support": args.high_min_support,
        "moderate_min_support": args.moderate_min_support,
        "labels_with_high_threshold": int(
            thresholds["tier"].eq("high").sum() if len(thresholds) else 0
        ),
        "labels_with_moderate_threshold": int(
            thresholds["tier"].eq("moderate").sum() if len(thresholds) else 0
        ),
        "query_high": int(counts.get("high", 0)),
        "query_moderate_only": int(counts.get("moderate_only", 0)),
        "query_moderate_or_higher": int(
            counts.get("high", 0) + counts.get("moderate_only", 0)
        ),
        "query_low_reject": int(counts.get("low_reject", 0)),
        "query_n": int(len(query)),
        "fine_anchor_eligible": False,
        "heldout_origin": origin["heldout_origin"],
        "calibration_origin_manifest": str(args.calibration_origin_manifest.resolve()),
        "calibration_origin_manifest_sha256": sha256_file(args.calibration_origin_manifest),
        "artifacts": {
            "thresholds": {"path": str(thresholds_path.resolve()), "sha256": sha256_file(thresholds_path)},
            "query_mapping": {"path": str(query_path.resolve()), "sha256": sha256_file(query_path)},
            "heldout_assignments": {"path": str(heldout_path.resolve()), "sha256": sha256_file(heldout_path)},
            "heldout_tier_validation": {"path": str(validation_path.resolve()), "sha256": sha256_file(validation_path)},
            "heldout_cumulative_validation": {"path": str(cumulative_path.resolve()), "sha256": sha256_file(cumulative_path)}
        },
        "calibration_semantics": "Nested cumulative thresholds: high is a subset of moderate-or-higher; mapping_tier uses canonical mutually exclusive high/moderate_only/low_reject values.",
        "warning": "High and moderate_only are mutually exclusive output tiers, but every high row also meets the cumulative moderate-or-higher gate. For an unlabeled frozen-QC observation, state-aware broad-only writeback additionally requires non-OOD, ontology-compatible and profile-scope-safe status; marker/spatial evidence remains an audit and coherent-group challenge.",
    }
    if sum(counts.values()) != len(query):
        raise SystemExit("tier counts do not partition query")
    (args.out / "calibration_summary.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
