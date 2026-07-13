#!/usr/bin/env python3
"""Extract raw positive/anti-marker evidence for candidate-program memberships."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def split_markers(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--program-manifest", required=True, type=Path)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--program-col", default="program")
    parser.add_argument("--layer", default="counts")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(
        args.candidates, sep="\t", dtype={args.cell_id_col: str}, low_memory=False
    )
    programs = pd.read_csv(args.program_manifest, sep="\t", low_memory=False)
    for column in (args.cell_id_col, args.program_col):
        if column not in candidates:
            raise SystemExit(f"candidates lack {column}")
    required_program = {"program", "positive_markers", "anti_markers"}
    if not required_program.issubset(programs.columns):
        raise SystemExit("program manifest lacks required columns")
    candidates = candidates.drop_duplicates([args.cell_id_col, args.program_col]).copy()
    if not set(candidates[args.program_col].astype(str)).issubset(set(programs["program"].astype(str))):
        raise SystemExit("candidate program is absent from program manifest")

    obj = ad.read_h5ad(args.h5ad)
    obj.obs_names = obj.obs_names.astype(str)
    if not obj.obs_names.is_unique:
        raise SystemExit("H5AD observation IDs are duplicated")
    missing = pd.Index(candidates[args.cell_id_col].astype(str)).difference(obj.obs_names)
    if len(missing):
        raise SystemExit(f"{len(missing)} candidate IDs are absent from H5AD")
    matrix = obj.layers[args.layer] if args.layer in obj.layers else obj.X
    matrix = matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    obs_lookup = pd.Series(np.arange(obj.n_obs), index=obj.obs_names)
    var_lookup = pd.Series(np.arange(obj.n_vars), index=obj.var_names.astype(str))

    candidate_rows = []
    gene_rows = []
    nonzero_rows = []
    program_lookup = programs.set_index("program")
    for program, group in candidates.groupby(args.program_col, sort=True):
        spec = program_lookup.loc[str(program)]
        positive_requested = split_markers(spec["positive_markers"])
        anti_requested = split_markers(spec["anti_markers"])
        positive = [gene for gene in positive_requested if gene in var_lookup]
        anti = [gene for gene in anti_requested if gene in var_lookup]
        genes = positive + [gene for gene in anti if gene not in positive]
        ids = group[args.cell_id_col].astype(str).tolist()
        row_indices = obs_lookup.loc[ids].to_numpy()
        gene_indices = var_lookup.loc[genes].to_numpy()
        block = matrix[row_indices][:, gene_indices].tocsr()
        role = np.array(["positive" if gene in positive else "anti" for gene in genes])
        detected = block > 0
        positive_mask = role == "positive"
        anti_mask = role == "anti"
        positive_detected = np.asarray(detected[:, positive_mask].sum(axis=1)).ravel().astype(int)
        anti_detected = np.asarray(detected[:, anti_mask].sum(axis=1)).ravel().astype(int)
        positive_counts = np.asarray(block[:, positive_mask].sum(axis=1)).ravel()
        anti_counts = np.asarray(block[:, anti_mask].sum(axis=1)).ravel()
        for index, cell_id in enumerate(ids):
            row = block.getrow(index)
            entries = [
                f"{genes[column]}:{float(value):g}"
                for column, value in zip(row.indices, row.data)
            ]
            positive_entries = [
                entry for entry in entries if entry.split(":", 1)[0] in positive
            ]
            anti_entries = [entry for entry in entries if entry.split(":", 1)[0] in anti]
            candidate_rows.append(
                {
                    "cell_id": cell_id,
                    "program": program,
                    "positive_detected": int(positive_detected[index]),
                    "anti_detected": int(anti_detected[index]),
                    "positive_raw_counts": float(positive_counts[index]),
                    "anti_raw_counts": float(anti_counts[index]),
                    "detected_positive_gene_counts": ";".join(positive_entries),
                    "detected_anti_gene_counts": ";".join(anti_entries),
                }
            )
            for column, value in zip(row.indices, row.data):
                nonzero_rows.append(
                    {
                        "cell_id": cell_id,
                        "program": program,
                        "gene": genes[column],
                        "marker_role": role[column],
                        "raw_count": float(value),
                    }
                )
        for column, gene in enumerate(genes):
            values = block[:, column]
            detected_n = int(values.getnnz())
            total = float(values.sum())
            gene_rows.append(
                {
                    "program": program,
                    "gene": gene,
                    "marker_role": role[column],
                    "candidate_n": len(ids),
                    "detected_n": detected_n,
                    "detected_fraction": detected_n / len(ids) if ids else 0,
                    "raw_count_sum": total,
                }
            )

    pd.DataFrame(candidate_rows).to_csv(
        args.out / "candidate_program_marker_summary.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(gene_rows).to_csv(
        args.out / "program_gene_detection_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(nonzero_rows).to_csv(
        args.out / "candidate_program_marker_nonzero.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    manifest = {
        "status": "PASS",
        "candidate_program_rows": len(candidate_rows),
        "unique_candidate_observations": int(candidates[args.cell_id_col].nunique()),
        "programs": int(candidates[args.program_col].nunique()),
        "nonzero_marker_rows": len(nonzero_rows),
        "source_layer": args.layer if args.layer in obj.layers else "X",
        "warning": "Raw marker counts are evidence only. Context-gated lineages still require their independent spatial and anti-program gates.",
    }
    (args.out / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
