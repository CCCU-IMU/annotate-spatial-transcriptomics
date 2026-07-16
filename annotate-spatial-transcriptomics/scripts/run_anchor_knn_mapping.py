#!/usr/bin/env python3
"""Diagnostic/small-reference kNN evidence; not the reusable global-Atlas path."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize


def matrix(obj: ad.AnnData, layer: str):
    if layer:
        if layer not in obj.layers:
            raise SystemExit(f"layer {layer!r} is absent")
        return obj.layers[layer]
    return obj.X


def lognorm(x):
    x = x.tocsr().astype(np.float32) if sparse.issparse(x) else sparse.csr_matrix(x, dtype=np.float32)
    totals = np.asarray(x.sum(axis=1)).ravel()
    x = sparse.diags(1e4 / np.maximum(totals, 1)) @ x
    x.data = np.log1p(x.data)
    return x


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--query-h5ad", required=True, type=Path)
    p.add_argument("--reference-h5ad", required=True, type=Path)
    p.add_argument("--reference-label-col", required=True)
    p.add_argument("--exclude-labels", default="")
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--membership", type=Path)
    p.add_argument("--cell-id-col", default="cell_id")
    p.add_argument("--query-layer", default="")
    p.add_argument("--reference-layer", default="")
    p.add_argument("--n-components", type=int, default=50)
    p.add_argument("--neighbors", type=int, default=25)
    p.add_argument("--max-reference-per-label", type=int, default=5000)
    p.add_argument("--seed", type=int, default=20260713)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    query = ad.read_h5ad(args.query_h5ad)
    reference = ad.read_h5ad(args.reference_h5ad)
    query.obs_names = query.obs_names.astype(str)
    reference.obs_names = reference.obs_names.astype(str)
    if not query.obs_names.is_unique or not reference.obs_names.is_unique:
        raise SystemExit("query or reference observation IDs are duplicated")
    if args.reference_label_col not in reference.obs:
        raise SystemExit("reference label column absent")

    membership_n = None
    membership_exact_order = False
    if args.membership is not None:
        membership = pd.read_csv(args.membership, sep="\t", dtype={args.cell_id_col: str})
        if args.cell_id_col not in membership:
            raise SystemExit(f"membership lacks {args.cell_id_col}")
        ids = pd.Index(membership[args.cell_id_col].astype(str), name=args.cell_id_col)
        if not ids.is_unique:
            raise SystemExit("membership contains duplicate query IDs")
        missing = ids.difference(query.obs_names)
        if len(missing):
            raise SystemExit(f"{len(missing)} query IDs absent")
        query = query[ids].copy()
        query.obs_names = ids
        membership_n = len(ids)
        membership_exact_order = query.obs_names.equals(ids)
        if len(query) != membership_n or not membership_exact_order:
            raise SystemExit("membership subset/order contract failed")

    excluded = {x.strip() for x in args.exclude_labels.split(",") if x.strip()}
    raw_label = reference.obs[args.reference_label_col]
    valid = (
        raw_label.notna()
        & raw_label.astype(str).str.len().gt(0)
        & ~raw_label.astype(str).isin(excluded)
    )
    reference = reference[valid].copy()
    genes = query.var_names.intersection(reference.var_names, sort=False)
    if len(genes) < 100:
        raise SystemExit(f"only {len(genes)} shared genes")
    query = query[:, genes]
    reference = reference[:, genes]

    rng = np.random.default_rng(args.seed)
    labels = reference.obs[args.reference_label_col].astype(str)
    keep = []
    for _, idx in labels.groupby(labels).groups.items():
        idx = np.asarray(list(idx))
        keep.extend(rng.choice(idx, min(len(idx), args.max_reference_per_label), replace=False))
    reference = reference[keep].copy()
    labels = reference.obs[args.reference_label_col].astype(str).to_numpy()

    query_matrix = lognorm(matrix(query, args.query_layer))
    reference_matrix = lognorm(matrix(reference, args.reference_layer))
    n_components = min(args.n_components, reference_matrix.shape[0] - 1, reference_matrix.shape[1] - 1)
    svd = TruncatedSVD(n_components=n_components, random_state=args.seed).fit(reference_matrix)
    reference_embedding = normalize(svd.transform(reference_matrix))
    query_embedding = normalize(svd.transform(query_matrix))
    neighbors = min(args.neighbors, len(reference))
    distance, index = NearestNeighbors(n_neighbors=neighbors, metric="cosine").fit(
        reference_embedding
    ).kneighbors(query_embedding)
    weights = np.maximum(1 - distance, 1e-6)
    classes = np.unique(labels)
    scores = np.zeros((len(query), len(classes)), dtype=np.float32)
    for j, label in enumerate(classes):
        scores[:, j] = (weights * (labels[index] == label)).sum(axis=1)
    scores /= np.maximum(scores.sum(axis=1, keepdims=True), 1e-12)
    order = np.argsort(scores, axis=1)
    top = order[:, -1]
    second = order[:, -2] if len(classes) > 1 else top
    output = pd.DataFrame(
        {
            "cell_id": query.obs_names.astype(str),
            "predicted_label": classes[top],
            "confidence": scores[np.arange(len(query)), top],
            "margin": scores[np.arange(len(query)), top]
            - scores[np.arange(len(query)), second],
            "n_shared_genes": len(genes),
            "reference_neighbors": neighbors,
        }
    )
    if not output["cell_id"].is_unique or len(output) != len(query):
        raise SystemExit("prediction boundary contract failed")
    if membership_n is not None and len(output) != membership_n:
        raise SystemExit("prediction rows do not equal membership rows")
    output.to_csv(
        args.out / "anchor_knn_predictions.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    manifest = {
        "status": "EVIDENCE_ONLY",
        "n_query": len(query),
        "membership": str(args.membership.resolve()) if args.membership else None,
        "membership_n": membership_n,
        "membership_exact_order": membership_exact_order if args.membership else None,
        "prediction_rows_equal_membership": membership_n is None or len(output) == membership_n,
        "n_reference": len(reference),
        "shared_genes": len(genes),
        "labels": classes.tolist(),
        "reusable_across_queries": False,
        "query_reference_joint_retraining": True,
        "global_atlas_default_eligible": False,
        "warning": "This script refits SVD per query and is diagnostic/small-reference only. Use a validated fixed Atlas transform/index for all-cell concordance; calibrate thresholds before any QC broad writeback.",
    }
    (args.out / "mapping_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
