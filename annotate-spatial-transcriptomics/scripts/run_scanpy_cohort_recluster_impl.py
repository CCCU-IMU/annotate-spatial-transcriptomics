#!/usr/bin/env python3
"""Query-only Scanpy cohort reclustering in a joint query/anchor PCA space."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--role-col", default="query_or_anchor")
    parser.add_argument("--anchor-label-col", default="anchor_label")
    parser.add_argument("--query-value", default="query")
    parser.add_argument("--anchor-value", default="anchor")
    parser.add_argument("--resolutions", default="0.2,0.4,0.6,0.8")
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=20)
    parser.add_argument("--n-hvg", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=20260713)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "tables").mkdir(exist_ok=True)
    (args.out / "figures").mkdir(exist_ok=True)

    source = ad.read_h5ad(args.h5ad)
    source.obs_names = source.obs_names.astype(str)
    membership = pd.read_csv(args.membership, sep="\t", dtype={args.cell_id_col: str})
    ids = pd.Index(membership[args.cell_id_col].astype(str))
    if ids.duplicated().any():
        raise RuntimeError("duplicate membership IDs")
    missing = ids.difference(source.obs_names)
    if len(missing):
        raise RuntimeError(f"{len(missing)} membership IDs absent")
    anchor_mode = args.role_col in membership
    if anchor_mode:
        if not membership[args.role_col].isin([args.query_value, args.anchor_value]).all():
            raise RuntimeError("membership role values must be query or anchor")
        if args.anchor_label_col not in membership:
            raise RuntimeError("anchors require anchor_label")
        query_ids = pd.Index(membership.loc[membership[args.role_col] == args.query_value, args.cell_id_col])
        anchor_ids = pd.Index(membership.loc[membership[args.role_col] == args.anchor_value, args.cell_id_col])
        anchor_labels = membership.loc[membership[args.role_col] == args.anchor_value, args.anchor_label_col].astype(str)
        if not len(query_ids) or not len(anchor_ids) or anchor_labels.eq("").any() or anchor_labels.nunique() < 2:
            raise RuntimeError("anchor-assisted mode requires query cells and at least two nonempty anchor labels")
    else:
        query_ids, anchor_ids = ids, pd.Index([])

    joint = source[ids].copy()
    membership = membership.set_index(args.cell_id_col).loc[joint.obs_names]
    for column in membership:
        joint.obs[column] = membership[column].to_numpy()
    if "counts" not in joint.layers:
        joint.layers["counts"] = joint.X.copy()
    total = np.asarray(joint.layers["counts"].sum(axis=1)).ravel()
    zero = total <= 0
    zero_query = joint.obs_names[zero & joint.obs_names.isin(query_ids)]
    zero_anchor = joint.obs_names[zero & joint.obs_names.isin(anchor_ids)]
    pd.DataFrame({"cell_id": zero_query, "route": "qc_holdout", "reason": "zero_count_query_in_selected_layer"}).to_csv(args.out / "tables/zero_count_observations.tsv", sep="\t", index=False)
    joint = joint[~zero].copy()
    query_ids = query_ids.intersection(joint.obs_names)
    anchor_ids = anchor_ids.intersection(joint.obs_names)
    if not len(query_ids) or (anchor_mode and not len(anchor_ids)):
        raise RuntimeError("zero-count filtering removed all query or anchor cells")

    sc.pp.normalize_total(joint, target_sum=1e4)
    sc.pp.log1p(joint)
    sc.pp.highly_variable_genes(joint, n_top_genes=min(args.n_hvg, joint.n_vars), flavor="seurat")
    sc.pp.scale(joint, max_value=10)
    pca_features = int(joint.var["highly_variable"].sum()) if "highly_variable" in joint.var else joint.n_vars
    ncomp = min(max(args.n_pcs, 30), joint.n_obs - 1, pca_features - 1)
    if ncomp < 2:
        raise RuntimeError("cohort too small for PCA")
    sc.tl.pca(joint, n_comps=ncomp, random_state=args.seed, use_highly_variable=True)

    if anchor_mode:
        anchor_frame = pd.DataFrame(joint.obsm["X_pca"][joint.obs_names.isin(anchor_ids), :args.n_pcs], index=joint.obs_names[joint.obs_names.isin(anchor_ids)])
        anchor_label_series = joint.obs.loc[anchor_frame.index, args.anchor_label_col].astype(str)
        centroids = anchor_frame.groupby(anchor_label_series).mean()
        query_frame = pd.DataFrame(joint.obsm["X_pca"][joint.obs_names.isin(query_ids), :args.n_pcs], index=joint.obs_names[joint.obs_names.isin(query_ids)])
        distances = np.stack([np.square(query_frame.to_numpy() - row).sum(axis=1) for row in centroids.to_numpy()], axis=1)
        order = np.argsort(distances, axis=1)
        evidence = pd.DataFrame({
            "cell_id": query_frame.index,
            "nearest_anchor_label": centroids.index.to_numpy()[order[:, 0]],
            "nearest_anchor_distance": distances[np.arange(len(query_frame)), order[:, 0]],
            "anchor_distance_margin": distances[np.arange(len(query_frame)), order[:, 1]] - distances[np.arange(len(query_frame)), order[:, 0]],
        }).set_index("cell_id")
        evidence.to_csv(args.out / "tables/query_anchor_distance_evidence.tsv", sep="\t")
    else:
        evidence = pd.DataFrame(index=query_ids)

    query = joint[query_ids].copy()
    sc.pp.neighbors(query, n_neighbors=min(args.n_neighbors, query.n_obs - 1), use_rep="X_pca", n_pcs=min(args.n_pcs, ncomp), random_state=args.seed)
    sc.tl.umap(query, random_state=args.seed)
    compositions, anchor_summaries, qc_summaries = [], [], []
    source_columns = [x for x in ["source_key", "state_tags", "spatial_tags", "qc_tags", "candidate_lineages"] if x in query.obs]
    qc_columns = [x for x in ["total_counts", "n_genes_by_counts", "nCount_Spatial", "nFeature_Spatial", "pct_counts_mt"] if x in query.obs]
    coordinates = ("sdimx", "sdimy") if {"sdimx", "sdimy"}.issubset(query.obs) else ("x", "y") if {"x", "y"}.issubset(query.obs) else None
    for resolution in map(float, args.resolutions.split(",")):
        tag = str(resolution).replace(".", "p")
        key = f"framework_res{tag}"
        sc.tl.leiden(query, resolution=resolution, key_added=key, random_state=args.seed)
        clusters = query.obs[key].astype(str)
        pd.DataFrame({"cell_id": query.obs_names, "cluster": clusters, "resolution": resolution}).to_csv(args.out / "tables" / f"{key}_clusters.tsv", sep="\t", index=False)
        if clusters.nunique() > 1:
            sc.tl.rank_genes_groups(query, key, method="wilcoxon", pts=True)
            deg = sc.get.rank_genes_groups_df(query, None)
        else:
            deg = pd.DataFrame(columns=["group", "names", "scores", "logfoldchanges", "pvals", "pvals_adj", "pct_nz_group", "pct_nz_reference"])
        deg.to_csv(args.out / "tables" / f"{key}_DEG_all.tsv", sep="\t", index=False)
        deg.groupby("group", group_keys=False).head(100).to_csv(args.out / "tables" / f"{key}_DEG_top100.tsv", sep="\t", index=False)
        for column in source_columns:
            frame = pd.DataFrame({"cluster": clusters, "value": query.obs[column].astype(str)})
            counts = frame.value_counts().rename("n").reset_index()
            counts["fraction"] = counts["n"] / counts.groupby("cluster")["n"].transform("sum")
            counts["resolution"], counts["field"] = resolution, column
            compositions.append(counts[["resolution", "cluster", "field", "value", "n", "fraction"]])
        if qc_columns:
            frame = query.obs[qc_columns].apply(pd.to_numeric, errors="coerce").assign(cluster=clusters, resolution=resolution)
            qc_summaries.append(frame.groupby(["resolution", "cluster"])[qc_columns].median().reset_index())
        if anchor_mode:
            frame = evidence.loc[query.obs_names].assign(cluster=clusters.to_numpy(), resolution=resolution)
            summary = frame.groupby(["resolution", "cluster", "nearest_anchor_label"]).agg(n=("nearest_anchor_distance", "size"), mean_distance=("nearest_anchor_distance", "mean"), mean_margin=("anchor_distance_margin", "mean")).reset_index()
            summary["fraction"] = summary["n"] / summary.groupby(["resolution", "cluster"])["n"].transform("sum")
            anchor_summaries.append(summary)
        sc.pl.umap(query, color=key, show=False)
        plt.savefig(args.out / "figures" / f"{key}_UMAP.png", dpi=320, bbox_inches="tight"); plt.savefig(args.out / "figures" / f"{key}_UMAP.pdf", bbox_inches="tight"); plt.close()
        if coordinates:
            xx = query.obs[coordinates[0]].astype(float).to_numpy(); yy = query.obs[coordinates[1]].astype(float).to_numpy(); categories = sorted(clusters.unique()); lut = {value: i for i, value in enumerate(categories)}
            fig, ax = plt.subplots(figsize=(9, 7)); ax.scatter(xx, yy, c=np.array([lut[x] for x in clusters]), s=.25, cmap="turbo", rasterized=True); ax.invert_yaxis(); ax.set_aspect("equal"); ax.set_axis_off(); ax.set_title(f"Query-only cohort spatial {key}"); fig.savefig(args.out / "figures" / f"{key}_spatial.png", dpi=360, bbox_inches="tight"); fig.savefig(args.out / "figures" / f"{key}_spatial.pdf", bbox_inches="tight"); plt.close(fig)
            for category in categories:
                selected = clusters.to_numpy() == category; fig, ax = plt.subplots(figsize=(9, 7)); ax.scatter(xx[~selected], yy[~selected], s=.12, c="#dddddd", rasterized=True); ax.scatter(xx[selected], yy[selected], s=.5, c="#d62728", rasterized=True); ax.invert_yaxis(); ax.set_aspect("equal"); ax.set_axis_off(); ax.set_title(f"{key} cluster {category} (n={selected.sum()})"); stem = args.out / "figures" / f"{key}_cluster_{category}_highlight"; fig.savefig(str(stem) + ".png", dpi=360, bbox_inches="tight"); fig.savefig(str(stem) + ".pdf", bbox_inches="tight"); plt.close(fig)
    if compositions:
        pd.concat(compositions, ignore_index=True).to_csv(args.out / "tables/cluster_source_state_composition.tsv", sep="\t", index=False)
    if qc_summaries:
        pd.concat(qc_summaries, ignore_index=True).to_csv(args.out / "tables/cluster_QC_summary.tsv", sep="\t", index=False)
    if anchor_summaries:
        pd.concat(anchor_summaries, ignore_index=True).to_csv(args.out / "tables/cluster_anchor_distance_summary.tsv", sep="\t", index=False)
    for column in query.obs.columns:
        if query.obs[column].dtype == object:
            query.obs[column] = query.obs[column].fillna("").astype(str)
    query.write_h5ad(args.out / "cohort_reclustered_query.h5ad", compression="gzip")
    manifest = {
        **vars(args), "anchor_assisted": anchor_mode, "query_only_graph_umap_deg": True,
        "n_query_input": int(len(query_ids) + len(zero_query)), "n_query_analyzed": int(query.n_obs),
        "n_anchors_analyzed": int(len(anchor_ids)), "n_zero_query_qc": int(len(zero_query)), "n_zero_anchor_excluded": int(len(zero_anchor)),
    }
    (args.out / "run_manifest.json").write_text(json.dumps(manifest, default=str, indent=2) + "\n")


if __name__ == "__main__":
    main()
