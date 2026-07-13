#!/usr/bin/env python3
"""Validate a standardized raw-count reference against a frozen membership."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--label-col", required=True)
    parser.add_argument("--exclude-membership", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    membership = pd.read_csv(
        args.membership, sep="\t", dtype={args.cell_id_col: str}, low_memory=False
    )
    errors: list[str] = []
    for column in (args.cell_id_col, args.label_col):
        if column not in membership:
            errors.append(f"membership lacks {column}")
    if errors:
        raise SystemExit("; ".join(errors))
    if not membership[args.cell_id_col].is_unique:
        errors.append("membership IDs are duplicated")

    obj = ad.read_h5ad(args.h5ad)
    obj.obs_names = obj.obs_names.astype(str)
    if not obj.obs_names.is_unique:
        errors.append("H5AD observation IDs are duplicated")
    expected_ids = pd.Index(membership[args.cell_id_col].astype(str))
    if obj.obs_names.tolist() != expected_ids.tolist():
        errors.append("H5AD observation IDs/order differ from frozen membership")
    if args.label_col not in obj.obs:
        errors.append(f"H5AD obs lacks {args.label_col}")
    else:
        observed_labels = obj.obs[args.label_col].astype(str).to_numpy()
        expected_labels = membership[args.label_col].astype(str).to_numpy()
        if not np.array_equal(observed_labels, expected_labels):
            errors.append("H5AD labels/order differ from frozen membership")

    matrix = obj.layers["counts"] if "counts" in obj.layers else obj.X
    matrix = matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    data = matrix.data
    integer_like_fraction = float(np.mean(np.isclose(data, np.rint(data)))) if data.size else 1.0
    if data.size and (not np.isfinite(data).all() or float(data.min()) < 0):
        errors.append("count matrix contains negative or non-finite values")
    if integer_like_fraction < 0.999:
        errors.append("count matrix is not integer-like")
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        errors.append("one or more reference observations have zero counts")

    overlap_n = 0
    if args.exclude_membership:
        excluded = pd.read_csv(
            args.exclude_membership, sep="\t", dtype={args.cell_id_col: str}, low_memory=False
        )
        if args.cell_id_col not in excluded:
            errors.append(f"exclude membership lacks {args.cell_id_col}")
        else:
            overlap_n = len(
                expected_ids.intersection(pd.Index(excluded[args.cell_id_col].astype(str)))
            )
            if overlap_n:
                errors.append(f"reference overlaps excluded/query membership by {overlap_n} IDs")

    counts = (
        membership.groupby(args.label_col, dropna=False)
        .size()
        .rename("n")
        .reset_index()
        .sort_values(args.label_col)
    )
    counts.to_csv(args.out / "reference_label_counts.tsv", sep="\t", index=False)
    manifest = {
        "status": "PASS" if not errors else "FAIL",
        "h5ad": str(args.h5ad.resolve()),
        "membership": str(args.membership.resolve()),
        "n_observations": int(obj.n_obs),
        "n_features": int(obj.n_vars),
        "n_labels": int(membership[args.label_col].nunique()),
        "membership_exact_order": obj.obs_names.tolist() == expected_ids.tolist(),
        "label_exact_order": args.label_col in obj.obs
        and np.array_equal(
            obj.obs[args.label_col].astype(str).to_numpy(),
            membership[args.label_col].astype(str).to_numpy(),
        ),
        "query_reference_overlap": overlap_n,
        "nnz": int(matrix.nnz),
        "integer_like_fraction": integer_like_fraction,
        "count_median": float(np.median(totals)),
        "count_min": float(np.min(totals)),
        "count_max": float(np.max(totals)),
        "fine_anchor_eligible": False,
        "errors": errors,
        "warning": "PASS validates a raw-count reference input only; labels remain calibration anchors and cannot be copied directly to query observations.",
    }
    (args.out / "validation_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
