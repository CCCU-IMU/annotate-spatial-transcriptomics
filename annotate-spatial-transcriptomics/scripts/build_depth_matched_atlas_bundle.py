#!/usr/bin/env python3
"""Build query-depth-matched atlas train/held-out bundles for multiple frozen routes."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def parse_csv(value: str) -> List[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "route"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--query-h5ad", required=True, type=Path)
    p.add_argument("--reference-h5ad", required=True, type=Path)
    p.add_argument("--reference-label-col", required=True)
    p.add_argument("--reference-kind", choices=["atlas", "internal_anchor"], default="atlas")
    p.add_argument("--membership", required=True, type=Path)
    p.add_argument("--cell-id-col", default="cell_id")
    p.add_argument("--strata-col", required=True)
    p.add_argument("--include-strata", default="")
    p.add_argument("--reference-labels", required=True)
    p.add_argument("--reference-filter-col")
    p.add_argument("--reference-filter-values", default="true,1,passed,pass")
    p.add_argument("--train-per-label", type=int, default=1500)
    p.add_argument("--heldout-per-label", type=int, default=500)
    p.add_argument("--min-query-positive-depth", type=int, default=1)
    p.add_argument("--seed", type=int, default=20260713)
    p.add_argument("--out", required=True, type=Path)
    return p.parse_args()


def integer_counts(obj: ad.AnnData) -> sparse.csr_matrix:
    x = obj.layers["counts"] if "counts" in obj.layers else obj.X
    x = x.tocsr().astype(np.float32) if sparse.issparse(x) else sparse.csr_matrix(x, dtype=np.float32)
    if x.data.size and (np.min(x.data) < 0 or not np.isfinite(x.data).all()):
        raise SystemExit("count matrix contains negative or non-finite values")
    if x.data.size and np.mean(np.isclose(x.data, np.rint(x.data))) < 0.999:
        raise SystemExit("count matrix is not integer-like")
    return x


def thin_matrix(mat: sparse.csr_matrix, target_depths: np.ndarray, rng: np.random.Generator) -> Tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    rows: List[int] = []
    cols: List[int] = []
    vals: List[int] = []
    before = np.zeros(mat.shape[0], dtype=np.int64)
    after = np.zeros(mat.shape[0], dtype=np.int64)
    for j in range(mat.shape[0]):
        row = mat.getrow(j)
        counts = np.rint(row.data).astype(np.int64)
        total = int(counts.sum())
        target = min(total, int(rng.choice(target_depths))) if total > 0 else 0
        before[j] = total
        after[j] = target
        if target <= 0:
            continue
        draw = rng.multinomial(target, counts / total)
        nz = draw > 0
        rows.extend([j] * int(nz.sum()))
        cols.extend(row.indices[nz].tolist())
        vals.extend(draw[nz].tolist())
    out = sparse.csr_matrix((vals, (rows, cols)), shape=mat.shape, dtype=np.float32)
    return out, before, after


def make_adata(
    mat: sparse.csr_matrix,
    source_obs: pd.DataFrame,
    var: pd.DataFrame,
    before: np.ndarray,
    after: np.ndarray,
) -> ad.AnnData:
    obs = source_obs.copy()
    obs["original_shared_counts"] = before
    obs["depth_matched_shared_counts"] = after
    out = ad.AnnData(X=mat, obs=obs, var=var.copy())
    out.layers["counts"] = mat.copy()
    return out


def main() -> None:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    membership = pd.read_csv(args.membership, sep="\t", dtype={args.cell_id_col: str})
    for col in [args.cell_id_col, args.strata_col]:
        if col not in membership:
            raise SystemExit(f"membership lacks {col}")
    if not membership[args.cell_id_col].is_unique:
        raise SystemExit("membership cell IDs are duplicated")

    query = ad.read_h5ad(args.query_h5ad)
    reference = ad.read_h5ad(args.reference_h5ad)
    query.obs_names = query.obs_names.astype(str)
    reference.obs_names = reference.obs_names.astype(str)
    if args.reference_label_col not in reference.obs:
        raise SystemExit("reference label column is absent")

    missing = pd.Index(membership[args.cell_id_col]).difference(query.obs_names)
    if len(missing):
        raise SystemExit(f"{len(missing)} membership IDs are absent from the standardized query")
    wanted_labels = parse_csv(args.reference_labels)
    label = reference.obs[args.reference_label_col].astype(str)
    keep = label.isin(wanted_labels).to_numpy()
    if args.reference_filter_col:
        if args.reference_filter_col not in reference.obs:
            raise SystemExit("reference filter column is absent")
        allowed = {x.lower() for x in parse_csv(args.reference_filter_values)}
        keep &= reference.obs[args.reference_filter_col].astype(str).str.lower().isin(allowed).to_numpy()
    reference = reference[keep].copy()
    labels = reference.obs[args.reference_label_col].astype(str)

    shared = query.var_names.intersection(reference.var_names, sort=False)
    if len(shared) < 500:
        raise SystemExit(f"only {len(shared)} shared genes")
    q = query[:, shared].copy()
    r = reference[:, shared].copy()
    qx = integer_counts(q)
    rx = integer_counts(r)
    q_lookup = pd.Series(np.arange(q.n_obs, dtype=np.int64), index=q.obs_names)

    rng = np.random.default_rng(args.seed)
    train_idx: List[int] = []
    heldout_idx: List[int] = []
    split_rows: List[dict] = []
    label_array = labels.to_numpy()
    for lab in wanted_labels:
        idx = np.flatnonzero(label_array == lab)
        rng.shuffle(idx)
        n_hold = min(args.heldout_per_label, max(0, len(idx) // 3))
        n_train = min(args.train_per_label, max(0, len(idx) - n_hold))
        if n_hold < 20 or n_train < 20:
            raise SystemExit(f"insufficient independent atlas cells for {lab}: train={n_train}, heldout={n_hold}")
        heldout_idx.extend(idx[:n_hold].tolist())
        train_idx.extend(idx[n_hold : n_hold + n_train].tolist())
        split_rows.append({"reference_label": lab, "n_available": len(idx), "n_train": n_train, "n_heldout": n_hold})
    pd.DataFrame(split_rows).to_csv(args.out / "reference_split_counts.tsv", sep="\t", index=False)

    include = set(parse_csv(args.include_strata))
    strata = membership[args.strata_col].astype(str)
    selected_strata = [x for x in strata.drop_duplicates().tolist() if not include or x in include]
    bundle_rows: List[dict] = []
    train_obs_base = r.obs.iloc[train_idx].copy()
    heldout_obs_base = r.obs.iloc[heldout_idx].copy()
    train_mat_base = rx[train_idx]
    heldout_mat_base = rx[heldout_idx]

    for route_index, route in enumerate(selected_strata):
        route_ids = pd.Index(membership.loc[strata.eq(route), args.cell_id_col].astype(str))
        qi = q_lookup.loc[route_ids].to_numpy()
        depths_all = np.asarray(qx[qi].sum(axis=1)).ravel().astype(np.int64)
        depths = depths_all[depths_all >= args.min_query_positive_depth]
        if len(depths) < 20:
            raise SystemExit(f"route {route} has only {len(depths)} positive-depth query observations")
        route_dir = args.out / safe_name(route)
        route_dir.mkdir(parents=True, exist_ok=True)
        route_rng = np.random.default_rng(args.seed + route_index * 1009)
        train_mat, train_before, train_after = thin_matrix(train_mat_base, depths, route_rng)
        heldout_mat, heldout_before, heldout_after = thin_matrix(heldout_mat_base, depths, route_rng)
        train = make_adata(train_mat, train_obs_base, q.var, train_before, train_after)
        heldout = make_adata(heldout_mat, heldout_obs_base, q.var, heldout_before, heldout_after)
        train.write_h5ad(route_dir / "reference_depth_matched_train.h5ad", compression="gzip")
        heldout.write_h5ad(route_dir / "reference_depth_matched_heldout.h5ad", compression="gzip")
        pd.DataFrame(
            {
                "cell_id": heldout.obs_names.astype(str),
                "true_label": heldout.obs[args.reference_label_col].astype(str).to_numpy(),
            }
        ).to_csv(route_dir / "heldout_truth.tsv", sep="\t", index=False)
        pd.DataFrame({args.cell_id_col: route_ids}).to_csv(route_dir / "query_membership.tsv", sep="\t", index=False)
        summary = {
            "status": "DIAGNOSTIC_ONLY",
            "route": route,
            "n_query": int(len(route_ids)),
            "n_query_positive_shared_depth": int(len(depths)),
            "n_query_zero_shared_depth": int(np.sum(depths_all < args.min_query_positive_depth)),
            "query_shared_count_median": float(np.median(depths)),
            "query_shared_count_q10": float(np.quantile(depths, 0.1)),
            "query_shared_count_q90": float(np.quantile(depths, 0.9)),
            "n_shared_genes": int(len(shared)),
            "n_train": int(train.n_obs),
            "n_heldout": int(heldout.n_obs),
            "train_depth_median_after": float(np.median(train_after)),
            "heldout_depth_median_after": float(np.median(heldout_after)),
            "fine_anchor_eligible": False,
            "reference_kind": args.reference_kind,
            "heldout_origin": "reference_self_classification",
            "eligible_for_query_rescue_calibration": False,
            "warning": f"{args.reference_kind} train/held-out cells are independent of the query but both originate from the same reference. This validates reference separability only and cannot calibrate query rescue.",
        }
        (route_dir / "depth_matching_manifest.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
        bundle_rows.append(summary)

    pd.DataFrame(bundle_rows).to_csv(args.out / "atlas_bundle_summary.tsv", sep="\t", index=False)
    manifest = {
        "status": "PASS",
        "query_h5ad": str(args.query_h5ad.resolve()),
        "reference_h5ad": str(args.reference_h5ad.resolve()),
        "reference_label_col": args.reference_label_col,
        "reference_kind": args.reference_kind,
        "reference_labels": wanted_labels,
        "n_shared_genes": int(len(shared)),
        "routes": selected_strata,
        "n_routes": len(selected_strata),
        "query_reference_overlap": int(len(query.obs_names.intersection(reference.obs_names))),
        "heldout_origin": "reference_self_classification",
        "eligible_for_query_rescue_calibration": False,
        "warning": f"All held-out cells originate from the same {args.reference_kind}. Use query-like held-out current-query anchors for any rescue threshold; these files are diagnostic only.",
    }
    if manifest["query_reference_overlap"]:
        raise SystemExit("query/reference observation IDs overlap")
    (args.out / "atlas_bundle_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
