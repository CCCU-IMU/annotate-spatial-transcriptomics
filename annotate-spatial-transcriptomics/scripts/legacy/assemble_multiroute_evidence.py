#!/usr/bin/env python3
"""Assemble legacy route evidence for migration-only adjudication."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "value"


def read_table(path: Path, cell_id_col: str = "cell_id") -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={cell_id_col: str}, low_memory=False)


def load_tiered_root(root: Path, routes: list[str], prefix: str) -> pd.DataFrame:
    frames = []
    keep = [
        "cell_id",
        "predicted_label",
        "confidence",
        "margin",
        "mapping_tier",
        "mapping_status",
        "fine_anchor_eligible",
    ]
    for route in routes:
        path = root / safe_name(route) / "calibrated_query_mapping.tsv.gz"
        frame = read_table(path)
        missing = sorted(set(keep).difference(frame.columns))
        if missing:
            raise SystemExit(f"{path} lacks {missing}")
        frame = frame[keep].copy()
        frame["evidence_route"] = route
        frames.append(frame)
    out = pd.concat(frames, ignore_index=True)
    if not out["cell_id"].is_unique:
        raise SystemExit(f"{prefix} tiered mapping has duplicate/overlapping IDs")
    return out.rename(
        columns={column: f"{prefix}_{column}" for column in keep if column != "cell_id"}
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--triage", required=True, type=Path)
    parser.add_argument("--external-tiered-root", required=True, type=Path)
    parser.add_argument("--internal-tiered-root", required=True, type=Path)
    parser.add_argument("--spatial-tiered-mapping", type=Path)
    parser.add_argument("--program-metrics", required=True, type=Path)
    parser.add_argument("--query-h5ad", required=True, type=Path)
    parser.add_argument("--anchor-distance", type=Path)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--route-col", default="triage_route")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    triage = read_table(args.triage, args.cell_id_col)
    if args.cell_id_col != "cell_id":
        triage = triage.rename(columns={args.cell_id_col: "cell_id"})
    if not triage["cell_id"].is_unique:
        raise SystemExit("triage IDs are duplicated")
    if args.route_col not in triage:
        raise SystemExit(f"triage lacks {args.route_col}")
    routes = (
        triage.loc[triage["provisional_state"].astype(str).eq("pending_review"), args.route_col]
        .astype(str)
        .drop_duplicates()
        .tolist()
    )
    pending_ids = triage.loc[
        triage["provisional_state"].astype(str).eq("pending_review"), "cell_id"
    ].astype(str)

    external = load_tiered_root(args.external_tiered_root, routes, "atlas")
    internal = load_tiered_root(args.internal_tiered_root, routes, "internal")
    for label, frame in (("external", external), ("internal", internal)):
        if set(frame["cell_id"]) != set(pending_ids):
            raise SystemExit(f"{label} mapping union differs from pending triage")

    program_long = read_table(args.program_metrics)
    required_program = {
        "cell_id",
        "program",
        "positive_hits",
        "positive_fraction",
        "anti_hits",
        "required_groups_pass",
        "strong_program",
    }
    if not required_program.issubset(program_long.columns):
        raise SystemExit("program metrics lack required columns")
    if program_long.duplicated(["cell_id", "program"]).any():
        raise SystemExit("program metrics contain duplicate cell/program rows")
    if set(program_long["cell_id"].astype(str)) != set(triage["cell_id"]):
        raise SystemExit("program metric membership differs from triage")
    programs = program_long["program"].astype(str).drop_duplicates().tolist()
    column_rows = []
    wide_parts = []
    for value in (
        "positive_hits",
        "positive_fraction",
        "anti_hits",
        "required_groups_pass",
        "strong_program",
    ):
        pivot = program_long.pivot(index="cell_id", columns="program", values=value)
        renamed = {}
        for program in programs:
            column = f"program__{safe_name(program)}__{value}"
            renamed[program] = column
            column_rows.append(
                {"program": program, "metric": value, "evidence_matrix_column": column}
            )
        wide_parts.append(pivot.rename(columns=renamed))
    program_wide = pd.concat(wide_parts, axis=1).reset_index()

    query = ad.read_h5ad(args.query_h5ad)
    query.obs_names = query.obs_names.astype(str)
    if not query.obs_names.is_unique or set(query.obs_names) != set(triage["cell_id"]):
        raise SystemExit("query H5AD membership differs from triage")
    counts = query.layers["counts"] if "counts" in query.layers else query.X
    counts = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(counts)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    features = np.diff(counts.indptr)
    obs = query.obs.copy()
    obs.index = obs.index.astype(str)
    coords = pd.DataFrame(
        {
            "cell_id": query.obs_names,
            "x": obs["x"].to_numpy() if "x" in obs else np.nan,
            "y": obs["y"].to_numpy() if "y" in obs else np.nan,
            "total_counts": totals,
            "detected_features": features,
        }
    )
    del query, counts

    evidence = triage.merge(external.drop(columns=["evidence_route"]), on="cell_id", how="left", validate="one_to_one")
    evidence = evidence.merge(
        internal.drop(columns=["evidence_route"]), on="cell_id", how="left", validate="one_to_one"
    )
    if args.spatial_tiered_mapping:
        spatial = read_table(args.spatial_tiered_mapping)
        spatial_keep = [
            "cell_id",
            "predicted_label",
            "confidence",
            "margin",
            "mapping_tier",
            "mapping_status",
            "fine_anchor_eligible",
        ]
        missing_spatial = sorted(set(spatial_keep).difference(spatial.columns))
        if missing_spatial:
            raise SystemExit(f"spatial tiered mapping lacks {missing_spatial}")
        if not spatial["cell_id"].is_unique or set(spatial["cell_id"]) != set(pending_ids):
            raise SystemExit("spatial tiered mapping membership differs from pending triage")
        spatial = spatial[spatial_keep].rename(
            columns={column: f"spatial_{column}" for column in spatial_keep if column != "cell_id"}
        )
        evidence = evidence.merge(spatial, on="cell_id", how="left", validate="one_to_one")
    evidence = evidence.merge(program_wide, on="cell_id", how="left", validate="one_to_one")
    evidence = evidence.merge(coords, on="cell_id", how="left", validate="one_to_one")
    if args.anchor_distance:
        distance = read_table(args.anchor_distance)
        if not distance["cell_id"].is_unique:
            raise SystemExit("anchor-distance IDs are duplicated")
        evidence = evidence.merge(distance, on="cell_id", how="left", validate="one_to_one")
    if len(evidence) != len(triage) or not evidence["cell_id"].is_unique:
        raise SystemExit("assembled evidence row boundary failed")
    pending_mask = evidence["provisional_state"].astype(str).eq("pending_review")
    if evidence.loc[pending_mask, "atlas_mapping_tier"].isna().any() or evidence.loc[
        pending_mask, "internal_mapping_tier"
    ].isna().any():
        raise SystemExit("pending observations lack a tiered mapping row")
    if args.spatial_tiered_mapping and evidence.loc[pending_mask, "spatial_mapping_tier"].isna().any():
        raise SystemExit("pending observations lack a spatial tiered mapping row")

    evidence.to_csv(
        args.out / "multiroute_evidence_matrix.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(column_rows).drop_duplicates().to_csv(
        args.out / "program_column_map.tsv", sep="\t", index=False
    )
    route_mapping = (
        evidence.loc[pending_mask]
        .groupby(
            [
                args.route_col,
                "atlas_predicted_label",
                "atlas_mapping_tier",
                "internal_predicted_label",
                "internal_mapping_tier",
            ],
            dropna=False,
        )
        .size()
        .rename("n")
        .reset_index()
        .sort_values([args.route_col, "n"], ascending=[True, False])
    )
    route_mapping.to_csv(args.out / "route_mapping_cross_counts.tsv", sep="\t", index=False)

    cluster_mapping = (
        evidence.loc[pending_mask]
        .groupby(
            [
                "cluster_res0p1",
                "cluster_res0p4",
                args.route_col,
                "atlas_predicted_label",
                "atlas_mapping_tier",
                "internal_predicted_label",
                "internal_mapping_tier",
            ],
            dropna=False,
        )
        .size()
        .rename("n")
        .reset_index()
        .sort_values(["cluster_res0p1", "cluster_res0p4", "n"], ascending=[True, True, False])
    )
    cluster_mapping.to_csv(args.out / "cluster_mapping_cross_counts.tsv", sep="\t", index=False)

    program_summary = program_long.merge(
        triage[["cell_id", "cluster_res0p1", "cluster_res0p4", args.route_col, "provisional_state"]],
        on="cell_id",
        how="left",
        validate="many_to_one",
    )
    program_summary["anti_dominant"] = program_summary["anti_hits"].to_numpy() > program_summary[
        "positive_hits"
    ].to_numpy()
    program_summary = (
        program_summary[program_summary["provisional_state"].astype(str).eq("pending_review")]
        .groupby(["cluster_res0p1", "cluster_res0p4", args.route_col, "program"], dropna=False)
        .agg(
            n_observations=("cell_id", "size"),
            mean_positive_hits=("positive_hits", "mean"),
            any_positive_fraction=("positive_hits", lambda values: float(values.ge(1).mean())),
            required_groups_pass_fraction=("required_groups_pass", "mean"),
            strong_program_fraction=("strong_program", "mean"),
            mean_anti_hits=("anti_hits", "mean"),
            anti_dominant_fraction=("anti_dominant", "mean"),
        )
        .reset_index()
    )
    program_summary.to_csv(args.out / "cluster_route_program_summary.tsv", sep="\t", index=False)

    qc_summary = (
        evidence.groupby(["cluster_res0p1", "cluster_res0p4", args.route_col, "provisional_state"], dropna=False)
        .agg(
            n_observations=("cell_id", "size"),
            median_total_counts=("total_counts", "median"),
            median_detected_features=("detected_features", "median"),
            q10_total_counts=("total_counts", lambda values: float(values.quantile(0.1))),
            q90_total_counts=("total_counts", lambda values: float(values.quantile(0.9))),
        )
        .reset_index()
    )
    qc_summary.to_csv(args.out / "cluster_route_qc_summary.tsv", sep="\t", index=False)

    manifest = {
        "status": "PASS",
        "n_observations": int(len(evidence)),
        "n_pending": int(pending_mask.sum()),
        "n_direct": int((~pending_mask).sum()),
        "n_programs": len(programs),
        "program_metric_rows": int(len(program_long)),
        "pending_external_mapping_complete": bool(
            evidence.loc[pending_mask, "atlas_mapping_tier"].notna().all()
        ),
        "pending_internal_mapping_complete": bool(
            evidence.loc[pending_mask, "internal_mapping_tier"].notna().all()
        ),
        "pending_spatial_mapping_complete": bool(
            not args.spatial_tiered_mapping
            or evidence.loc[pending_mask, "spatial_mapping_tier"].notna().all()
        ),
        "fine_anchor_eligible": False,
        "warning": "This is an evidence matrix, not a biological writeback. Mapping tiers require independent marker/anti-marker, cluster and spatial adjudication.",
    }
    (args.out / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
