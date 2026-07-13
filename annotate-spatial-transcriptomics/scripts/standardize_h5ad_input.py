#!/usr/bin/env python3
"""Standardize Scanpy- or stereopy-flavoured H5AD counts for frozen memberships."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List, Optional, Tuple, Union

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--membership", type=Path)
    p.add_argument("--cell-id-col", default="cell_id")
    p.add_argument("--preserve-membership-columns", action="store_true")
    p.add_argument("--input-flavor", choices=["auto", "scanpy", "stereopy"], default="auto")
    p.add_argument("--scanpy-layer")
    p.add_argument("--stereopy-matrix", choices=["raw", "processed"], default="raw")
    p.add_argument("--output-name", default="standardized_query_counts.h5ad")
    p.add_argument("--compression", choices=["gzip", "lzf", "none"], default="gzip")
    return p.parse_args()


def read_membership(path: Optional[Path], cell_id_col: str) -> Optional[pd.Index]:
    if path is None:
        return None
    frame = pd.read_csv(path, sep="\t", dtype={cell_id_col: str})
    if cell_id_col not in frame:
        raise SystemExit(f"membership lacks {cell_id_col}")
    ids = pd.Index(frame[cell_id_col].astype(str), name=cell_id_col)
    if not ids.is_unique:
        raise SystemExit("membership contains duplicate cell IDs")
    return ids


def matrix_audit(x: Union[sparse.spmatrix, np.ndarray]) -> dict:
    m = x.tocsr() if sparse.issparse(x) else sparse.csr_matrix(x)
    if m.data.size and (not np.isfinite(m.data).all() or np.min(m.data) < 0):
        raise SystemExit("expression matrix contains non-finite or negative values")
    rounded = np.rint(m.data)
    integer_fraction = float(np.mean(np.isclose(m.data, rounded))) if m.data.size else 1.0
    totals = np.asarray(m.sum(axis=1)).ravel()
    features = np.diff(m.indptr)
    return {
        "nnz": int(m.nnz),
        "integer_like_fraction": integer_fraction,
        "count_median": float(np.median(totals)),
        "count_min": float(np.min(totals)) if len(totals) else 0.0,
        "count_max": float(np.max(totals)) if len(totals) else 0.0,
        "detected_feature_median": float(np.median(features)),
        "zero_count_observations": int(np.sum(totals <= 0)),
    }


def read_scanpy(path: Path, layer: Optional[str]) -> Tuple[ad.AnnData, dict]:
    obj = ad.read_h5ad(path)
    obj.obs_names = obj.obs_names.astype(str)
    if layer:
        if layer not in obj.layers:
            raise SystemExit(f"scanpy layer {layer!r} is absent")
        obj.X = obj.layers[layer].copy()
    return obj, {"detected_flavor": "scanpy", "source_layer": layer or "X"}


def decode_strings(values: np.ndarray) -> np.ndarray:
    return np.asarray([v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values])


def read_stereopy(path: Path, matrix_choice: str) -> Tuple[ad.AnnData, dict]:
    matrix_key = "exp_matrix@raw" if matrix_choice == "raw" else "exp_matrix"
    cell_key = "cells@raw/obs/_index" if matrix_choice == "raw" else "cells/obs/_index"
    gene_key = "genes@raw/var/_index" if matrix_choice == "raw" else "genes/var/_index"
    position_key = "position@raw" if matrix_choice == "raw" else "position"
    with h5py.File(path, "r") as handle:
        required = [matrix_key, cell_key, gene_key]
        missing = [key for key in required if key not in handle]
        if missing:
            raise RuntimeError("stereopy keys absent: " + ",".join(missing))
        group = handle[matrix_key]
        if group.attrs.get("encoding-type", "") != "csr_matrix":
            raise RuntimeError(f"{matrix_key} is not a stereopy CSR matrix")
        shape = tuple(int(v) for v in group.attrs["shape"])
        data = group["data"][:]
        indices = group["indices"][:]
        indptr = group["indptr"][:]
        x = sparse.csr_matrix((data, indices, indptr), shape=shape)
        obs_names = pd.Index(decode_strings(handle[cell_key][:]), name="cell_id")
        var_names = pd.Index(decode_strings(handle[gene_key][:]), name="gene")
        position = handle[position_key][:] if position_key in handle else None
    if len(obs_names) != x.shape[0] or len(var_names) != x.shape[1]:
        raise RuntimeError("stereopy matrix/name dimensions disagree")
    obs = pd.DataFrame(index=obs_names)
    if position is not None and position.ndim == 2 and position.shape[1] >= 2:
        obs["x"] = position[:, 0]
        obs["y"] = position[:, 1]
    obj = ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=var_names))
    return obj, {
        "detected_flavor": "stereopy",
        "source_layer": matrix_key,
        "stereopy_direct_hdf5_reader": True,
        "source_position_available": position is not None,
    }


def main() -> None:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    membership = read_membership(args.membership, args.cell_id_col)

    failures: List[str] = []
    obj: Optional[ad.AnnData] = None
    source_info: dict = {}
    if args.input_flavor in {"auto", "scanpy"}:
        try:
            obj, source_info = read_scanpy(args.input, args.scanpy_layer)
        except Exception as exc:  # preserve the adapter decision in provenance
            failures.append(f"scanpy:{type(exc).__name__}:{exc}")
            if args.input_flavor == "scanpy":
                raise
    if obj is None and args.input_flavor in {"auto", "stereopy"}:
        try:
            obj, source_info = read_stereopy(args.input, args.stereopy_matrix)
        except Exception as exc:
            failures.append(f"stereopy:{type(exc).__name__}:{exc}")
            raise
    if obj is None:
        raise SystemExit("no input adapter succeeded")

    obj.obs_names = obj.obs_names.astype(str)
    obj.var_names = obj.var_names.astype(str)
    if not obj.obs_names.is_unique:
        raise SystemExit("source observation IDs are not unique")
    if not obj.var_names.is_unique:
        raise SystemExit("source gene IDs are not unique")

    source_n_obs, source_n_vars = obj.shape
    if membership is not None:
        missing = membership.difference(obj.obs_names)
        if len(missing):
            raise SystemExit(f"{len(missing)} membership IDs are absent from the source")
        lookup = pd.Series(np.arange(obj.n_obs, dtype=np.int64), index=obj.obs_names)
        order = lookup.loc[membership].to_numpy()
        obj = obj[order].copy()
        obj.obs_names = membership
        if args.preserve_membership_columns:
            membership_frame = pd.read_csv(
                args.membership, sep="\t", dtype={args.cell_id_col: str}
            ).set_index(args.cell_id_col)
            membership_frame.index = membership_frame.index.astype(str)
            membership_frame = membership_frame.loc[membership]
            for column in membership_frame.columns:
                values = membership_frame[column]
                if pd.api.types.is_object_dtype(values.dtype):
                    obj.obs[column] = values.fillna("").astype(str).to_numpy()
                else:
                    obj.obs[column] = values.to_numpy()
    elif not obj.is_view:
        obj = obj.copy()

    obj.layers.clear()
    obj.raw = None
    obj.X = obj.X.tocsr() if sparse.issparse(obj.X) else sparse.csr_matrix(obj.X)
    audit = matrix_audit(obj.X)
    if audit["integer_like_fraction"] < 0.999:
        raise SystemExit("selected expression matrix is not integer-like raw counts")
    obj.layers["counts"] = obj.X.copy()

    output = args.out / args.output_name
    compression = None if args.compression == "none" else args.compression
    obj.write_h5ad(output, compression=compression)
    manifest = {
        "status": "PASS",
        "input": str(args.input.resolve()),
        "output": str(output.resolve()),
        "source_n_observations": int(source_n_obs),
        "source_n_features": int(source_n_vars),
        "selected_n_observations": int(obj.n_obs),
        "selected_n_features": int(obj.n_vars),
        "membership": str(args.membership.resolve()) if args.membership else None,
        "membership_exact_order": membership is not None and obj.obs_names.equals(membership),
        "membership_columns_preserved": bool(args.preserve_membership_columns),
        "adapter_failures": failures,
        **source_info,
        **audit,
        "warning": "This adapter preserves raw counts and frozen membership only; it does not transfer annotations or clustering results.",
    }
    (args.out / "standardization_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
