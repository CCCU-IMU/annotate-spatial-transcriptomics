#!/usr/bin/env python3
"""Observed-density spatial-anchor evidence with independent held-out calibration."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


def truthy(values: pd.Series) -> pd.Series:
    return values.astype(str).str.lower().isin({"true", "1", "yes", "eligible"})


def predict(
    query: pd.DataFrame,
    reference: pd.DataFrame,
    label_col: str,
    x_col: str,
    y_col: str,
    neighbors: int,
) -> pd.DataFrame:
    k = min(neighbors, len(reference))
    model = NearestNeighbors(n_neighbors=k, metric="euclidean").fit(
        reference[[x_col, y_col]].to_numpy(dtype=float)
    )
    distances, indices = model.kneighbors(query[[x_col, y_col]].to_numpy(dtype=float))
    labels = reference[label_col].astype(str).to_numpy()
    classes = np.unique(labels)
    scores = np.zeros((len(query), len(classes)), dtype=np.float32)
    for index, label in enumerate(classes):
        scores[:, index] = (labels[indices] == label).mean(axis=1)
    order = np.argsort(scores, axis=1)
    top = order[:, -1]
    second = order[:, -2] if len(classes) > 1 else top
    return pd.DataFrame(
        {
            "cell_id": query["cell_id"].astype(str).to_numpy(),
            "predicted_label": classes[top],
            "confidence": scores[np.arange(len(query)), top],
            "margin": scores[np.arange(len(query)), top]
            - scores[np.arange(len(query)), second],
            "reference_neighbors": k,
            "nearest_distance": distances[:, 0],
            "median_neighbor_distance": np.median(distances, axis=1),
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-table", required=True, type=Path)
    parser.add_argument("--query-table", required=True, type=Path)
    parser.add_argument("--label-col", required=True)
    parser.add_argument("--calibration-eligible-col", required=True)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--x-col", default="x")
    parser.add_argument("--y-col", default="y")
    parser.add_argument("--neighbors", type=int, default=50)
    parser.add_argument("--heldout-fraction", type=float, default=0.2)
    parser.add_argument("--minimum-heldout-per-label", type=int, default=20)
    parser.add_argument("--maximum-heldout-per-label", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20260713)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    reference = pd.read_csv(
        args.reference_table, sep="\t", dtype={args.cell_id_col: str}, low_memory=False
    ).rename(columns={args.cell_id_col: "cell_id"})
    query = pd.read_csv(
        args.query_table, sep="\t", dtype={args.cell_id_col: str}, low_memory=False
    ).rename(columns={args.cell_id_col: "cell_id"})
    required_reference = {
        "cell_id",
        args.label_col,
        args.calibration_eligible_col,
        args.x_col,
        args.y_col,
    }
    required_query = {"cell_id", args.x_col, args.y_col}
    if not required_reference.issubset(reference.columns):
        raise SystemExit("reference table lacks required columns")
    if not required_query.issubset(query.columns):
        raise SystemExit("query table lacks required columns")
    if not reference["cell_id"].is_unique or not query["cell_id"].is_unique:
        raise SystemExit("reference or query IDs are duplicated")
    if set(reference["cell_id"]).intersection(set(query["cell_id"])):
        raise SystemExit("spatial reference and query IDs overlap")
    for frame, name in ((reference, "reference"), (query, "query")):
        coordinates = frame[[args.x_col, args.y_col]].to_numpy(dtype=float)
        if not np.isfinite(coordinates).all():
            raise SystemExit(f"{name} contains non-finite coordinates")

    calibration_eligible = truthy(reference[args.calibration_eligible_col])
    rng = np.random.default_rng(args.seed)
    heldout_indices: list[int] = []
    split_rows: list[dict] = []
    for label in sorted(reference.loc[calibration_eligible, args.label_col].astype(str).unique()):
        indices = reference.index[
            calibration_eligible & reference[args.label_col].astype(str).eq(label)
        ].to_numpy()
        rng.shuffle(indices)
        requested = int(round(len(indices) * args.heldout_fraction))
        heldout_n = min(
            args.maximum_heldout_per_label,
            max(args.minimum_heldout_per_label, requested),
            max(0, len(indices) - args.minimum_heldout_per_label),
        )
        if heldout_n < args.minimum_heldout_per_label:
            raise SystemExit(f"insufficient calibration-eligible references for {label}")
        heldout_indices.extend(indices[:heldout_n].tolist())
        split_rows.append(
            {
                "reference_label": label,
                "calibration_eligible_n": len(indices),
                "heldout_n": heldout_n,
                "calibration_train_n": len(indices) - heldout_n,
                "always_train_n": int(
                    ((~calibration_eligible) & reference[args.label_col].astype(str).eq(label)).sum()
                ),
            }
        )
    heldout = reference.loc[heldout_indices].copy()
    train = reference.drop(index=heldout_indices).copy()
    if not set(train[args.label_col].astype(str)) == set(reference[args.label_col].astype(str)):
        raise SystemExit("held-out split removed a reference label from training")

    heldout_predictions = predict(
        heldout, train, args.label_col, args.x_col, args.y_col, args.neighbors
    )
    query_predictions = predict(
        query, train, args.label_col, args.x_col, args.y_col, args.neighbors
    )
    if heldout_predictions["cell_id"].tolist() != heldout["cell_id"].astype(str).tolist():
        raise SystemExit("held-out prediction order contract failed")
    if query_predictions["cell_id"].tolist() != query["cell_id"].astype(str).tolist():
        raise SystemExit("query prediction order contract failed")
    heldout_predictions.to_csv(
        args.out / "heldout_spatial_predictions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    query_predictions.to_csv(
        args.out / "query_spatial_predictions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(
        {
            "cell_id": heldout["cell_id"].astype(str).to_numpy(),
            "true_label": heldout[args.label_col].astype(str).to_numpy(),
        }
    ).to_csv(args.out / "heldout_truth.tsv", sep="\t", index=False)
    pd.DataFrame(split_rows).to_csv(args.out / "reference_split_counts.tsv", sep="\t", index=False)
    train[["cell_id", args.label_col, args.x_col, args.y_col]].to_csv(
        args.out / "spatial_train_reference.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    manifest = {
        "status": "SPATIAL_EVIDENCE_ONLY",
        "n_reference": int(len(reference)),
        "n_train": int(len(train)),
        "n_heldout": int(len(heldout)),
        "n_query": int(len(query)),
        "n_labels": int(reference[args.label_col].nunique()),
        "neighbors": min(args.neighbors, len(train)),
        "query_reference_overlap": 0,
        "heldout_independent": True,
        "query_membership_exact_order": True,
        "fine_anchor_eligible": False,
        "warning": "Observed-density spatial proximity is an independent evidence channel, not a biological label. Calibrate tiers and combine it with expression evidence before writeback.",
    }
    (args.out / "mapping_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
