#!/usr/bin/env python3
"""Build the immutable split-wall GSE233801 reference-only Atlas assets."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import shutil
import zipfile
from pathlib import Path

import anndata as ad
import h5py
import joblib
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def counts(obj: ad.AnnData) -> sparse.csr_matrix:
    value = obj.layers["counts"] if "counts" in obj.layers else obj.X
    matrix = (
        value.tocsr().astype(np.float32)
        if sparse.issparse(value)
        else sparse.csr_matrix(value, dtype=np.float32)
    )
    if matrix.data.size and (
        np.min(matrix.data) < 0 or not np.isfinite(matrix.data).all()
    ):
        raise SystemExit("reference counts contain negative or non-finite values")
    return matrix


def lognorm(matrix: sparse.csr_matrix) -> sparse.csr_matrix:
    total = np.asarray(matrix.sum(axis=1)).ravel()
    result = sparse.diags(1e4 / np.maximum(total, 1)) @ matrix
    result.data = np.log1p(result.data)
    return result.tocsr()


def read_h5ad_csr_counts_subset(
    path: Path, row_index: np.ndarray, column_index: np.ndarray
) -> sparse.csr_matrix:
    """Sequentially read the count CSR and then subset in memory.

    AnnData backed column slicing performs expensive random HDF5 access on this
    100-million-nonzero reference. A single sequential CSR read is faster and
    keeps the fixed-feature boundary explicit.
    """
    with h5py.File(path, "r") as handle:
        key = "layers/counts" if "counts" in handle.get("layers", {}) else "X"
        group = handle[key]
        if str(group.attrs.get("encoding-type", "")) != "csr_matrix":
            raise SystemExit("GSE233801 count matrix is not CSR encoded")
        shape = tuple(int(value) for value in group.attrs["shape"])
        matrix = sparse.csr_matrix(
            (
                group["data"][:],
                group["indices"][:],
                group["indptr"][:],
            ),
            shape=shape,
        )
    return matrix[row_index][:, column_index].astype(np.float32).tocsr()


def predict(
    embedding: np.ndarray, prototypes: np.ndarray, labels: np.ndarray
) -> pd.DataFrame:
    similarity = embedding @ prototypes.T
    order = np.argsort(similarity, axis=1)
    top = order[:, -1]
    second = order[:, -2]
    return pd.DataFrame(
        {
            "predicted_label": labels[top],
            "confidence": (similarity[np.arange(len(top)), top] + 1.0) / 2.0,
            "margin": similarity[np.arange(len(top)), top]
            - similarity[np.arange(len(top)), second],
            "top_similarity": similarity[np.arange(len(top)), top],
        }
    )


def write_deterministic_npz(path: Path, **arrays: np.ndarray) -> None:
    """Write an npz without wall-clock timestamps in the zip members."""
    with zipfile.ZipFile(path, "w") as archive:
        for name in sorted(arrays):
            buffer = io.BytesIO()
            np.lib.format.write_array(buffer, np.asarray(arrays[name]), allow_pickle=False)
            info = zipfile.ZipInfo(f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o600 << 16
            archive.writestr(info, buffer.getvalue())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--cluster-crosswalk", required=True, type=Path)
    ap.add_argument("--fixed-features", required=True, type=Path)
    ap.add_argument("--existing-transform", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--cluster-col", default="reference_cluster")
    ap.add_argument("--dimensions", type=int, default=50)
    ap.add_argument("--train-per-label", type=int, default=5000)
    ap.add_argument("--heldout-per-label", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20260730)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    print("loading GSE233801 metadata", flush=True)

    crosswalk = pd.read_csv(args.cluster_crosswalk, sep="\t", dtype=str).fillna("")
    required = {
        "reference_cluster",
        "framework_broad_label",
        "atlas_capability",
        "include_in_prototype",
    }
    if not required.issubset(crosswalk.columns):
        raise SystemExit("cluster crosswalk lacks required fields")
    if not crosswalk.reference_cluster.is_unique:
        raise SystemExit("cluster crosswalk has duplicated reference_cluster")
    include = crosswalk.include_in_prototype.str.upper().eq("TRUE")
    included = crosswalk.loc[include].copy()
    if (included.framework_broad_label == "").any():
        raise SystemExit("included reference cluster lacks a framework broad label")

    reference_source = ad.read_h5ad(args.reference, backed="r")
    reference_source.obs_names = reference_source.obs_names.astype(str)
    if args.cluster_col not in reference_source.obs:
        raise SystemExit(f"reference lacks {args.cluster_col}")
    feature_table = pd.read_csv(args.fixed_features, sep="\t", dtype=str)
    if "gene" not in feature_table or not feature_table.gene.is_unique:
        raise SystemExit("fixed feature table lacks a unique gene column")
    var_names = pd.Index(reference_source.var_names.astype(str))
    features = [gene for gene in feature_table.gene if gene in var_names]
    if len(features) < 500:
        raise SystemExit("insufficient fixed features in GSE233801")
    if features != feature_table.gene.tolist():
        missing = sorted(set(feature_table.gene) - set(features))
        raise SystemExit(f"frozen feature boundary changed; missing {len(missing)} genes")

    cluster = reference_source.obs[args.cluster_col].astype(str)
    label_by_cluster = dict(
        zip(included.reference_cluster, included.framework_broad_label)
    )
    keep = cluster.isin(label_by_cluster)
    row_index = np.flatnonzero(keep.to_numpy())
    column_index = var_names.get_indexer(features)
    reference_cell_ids = reference_source.obs_names[row_index].astype(str).to_numpy()
    kept_cluster = cluster.iloc[row_index].reset_index(drop=True)
    reference_source.file.close()
    labels_all = kept_cluster.map(label_by_cluster).astype(str).to_numpy()
    print(
        f"reading fixed-feature counts for {len(labels_all)} reference cells",
        flush=True,
    )

    label_names = np.array(sorted(np.unique(labels_all)), dtype=str)
    rng = np.random.default_rng(args.seed)
    train_idx: list[int] = []
    heldout_idx: list[int] = []
    split_rows: list[dict[str, object]] = []
    for label in label_names:
        index = np.flatnonzero(labels_all == label)
        rng.shuffle(index)
        n_hold = min(args.heldout_per_label, max(20, len(index) // 4))
        n_train = min(args.train_per_label, len(index) - n_hold)
        if n_train < 20 or n_hold < 20:
            raise SystemExit(f"underpowered reference label {label}")
        heldout_idx.extend(index[:n_hold])
        train_idx.extend(index[n_hold : n_hold + n_train])
        split_rows.append(
            {
                "framework_broad_label": label,
                "available": len(index),
                "train": n_train,
                "heldout": n_hold,
            }
        )

    matrix = lognorm(
        read_h5ad_csr_counts_subset(args.reference, row_index, column_index)
    )
    print("preparing deterministic reference-only projection/prototypes", flush=True)
    train_index = np.asarray(train_idx, dtype=int)
    heldout_index = np.asarray(heldout_idx, dtype=int)
    if args.existing_transform:
        frozen_model = joblib.load(args.existing_transform)
        if list(frozen_model.get("features", [])) != features:
            raise SystemExit("existing transform has a different fixed-feature boundary")
        svd = frozen_model.get("svd")
        if svd is None or not hasattr(svd, "transform"):
            raise SystemExit("existing transform lacks a reusable reference-only SVD")
        n_dimensions = int(getattr(svd, "n_components", args.dimensions))
        print("applying frozen query-independent GSE233801 transform", flush=True)
    else:
        n_dimensions = min(
            args.dimensions, len(train_index) - 1, len(features) - 1
        )
        svd = TruncatedSVD(
            n_components=n_dimensions, random_state=args.seed
        ).fit(matrix[train_index])
    train_embedding = normalize(svd.transform(matrix[train_index]))
    heldout_embedding = normalize(svd.transform(matrix[heldout_index]))
    train_labels = labels_all[train_index]
    prototypes = normalize(
        np.vstack(
            [
                train_embedding[train_labels == label].mean(axis=0)
                for label in label_names
            ]
        )
    )

    heldout = predict(heldout_embedding, prototypes, label_names)
    heldout.insert(0, "cell_id", reference_cell_ids[heldout_index])
    heldout["true_label"] = labels_all[heldout_index]
    heldout["correct"] = heldout.predicted_label.eq(heldout.true_label)

    feature_path = args.out / "fixed_features.tsv"
    transform_path = args.out / "feature_transform.joblib"
    prototype_path = args.out / "reference_prototypes.npz"
    heldout_path = args.out / "reference_heldout_predictions.tsv.gz"
    crosswalk_path = args.out / "reference_cluster_crosswalk.tsv"
    split_path = args.out / "reference_split.tsv"
    shutil.copy2(args.fixed_features, feature_path)
    shutil.copy2(args.cluster_crosswalk, crosswalk_path)
    if args.existing_transform:
        shutil.copy2(args.existing_transform, transform_path)
    else:
        joblib.dump(
            {"features": features, "svd": svd, "seed": args.seed},
            transform_path,
            compress=3,
        )
    write_deterministic_npz(
        prototype_path, labels=label_names, prototypes=prototypes
    )
    heldout.to_csv(
        heldout_path,
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.DataFrame(split_rows).to_csv(split_path, sep="\t", index=False)

    assets = {
        path.name: sha256(path)
        for path in (
            feature_path,
            transform_path,
            prototype_path,
            heldout_path,
            crosswalk_path,
            split_path,
        )
    }
    manifest = {
        "schema_version": "2.0",
        "bundle_id": "sheep_ovary_GSE233801_split_wall_v2",
        "reference_id": "GSE233801_independent_R_res0p4_split_wall_v002",
        "mapping_engine": "fixed_projection_prototype",
        "n_reference_available": int(len(labels_all)),
        "n_reference_train": int(len(train_index)),
        "n_reference_heldout": int(len(heldout_index)),
        "n_dimensions": int(n_dimensions),
        "n_features": int(len(features)),
        "labels": label_names.tolist(),
        "asset_hashes": assets,
        "reference_self_classification_role": "diagnostic_only",
        "query_classwise_calibration_required_for_direct_rescue": True,
        "query_reference_joint_retraining": False,
        "feature_transform_reused_query_independent_v1": bool(
            args.existing_transform
        ),
    }
    (args.out / "atlas_index_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print("split-wall Atlas assets complete", flush=True)
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
