#!/usr/bin/env python3
"""Compare clustering candidates without erasing rare/microcluster observations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score


def safe_scores(left: np.ndarray, right: np.ndarray) -> tuple[float, float]:
    if len(left) < 2:
        return float("nan"), float("nan")
    return adjusted_rand_score(left, right), adjusted_mutual_info_score(left, right)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("grid_summary", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--ranking", type=Path, help="Optional rank table used to shortlist")
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--max-cells", type=int, default=50000)
    ap.add_argument("--seed", type=int, default=20260713)
    ap.add_argument("--microcluster-threshold", type=int, default=100,
                    help="Clusters below this size are excluded only from macro-restricted scores")
    args = ap.parse_args()
    if args.microcluster_threshold < 1:
        raise ValueError("--microcluster-threshold must be positive")
    args.out.mkdir(parents=True, exist_ok=True)

    grid = pd.read_csv(args.grid_summary)
    if args.ranking:
        rank = pd.read_csv(args.ranking, sep="\t").sort_values("quantitative_rank").head(args.top)
        grid = grid[grid.run_id.isin(rank.run_id)]
    elif len(grid) > args.top:
        grid = grid.head(args.top)

    labels: dict[str, pd.Series] = {}
    common_ids: pd.Index | None = None
    for row in grid.itertuples(index=False):
        data = pd.read_csv(row.csv_path, usecols=["cell_id", "cluster"],
                           dtype={"cell_id": str, "cluster": str})
        if data.cell_id.duplicated().any():
            raise RuntimeError(f"duplicate IDs in {row.run_id}")
        series = data.set_index("cell_id").cluster
        common_ids = series.index if common_ids is None else common_ids.intersection(series.index)
        labels[row.run_id] = series
    if common_ids is None or not len(common_ids):
        raise RuntimeError("no common cells")

    runs = list(labels)
    full_labels = {run: labels[run].reindex(common_ids) for run in runs}
    cluster_sizes = {run: full_labels[run].value_counts() for run in runs}
    audit_rows = []
    for run in runs:
        input_sizes = labels[run].value_counts()
        for cluster, n_common in cluster_sizes[run].items():
            audit_rows.append({
                "run_id": run,
                "cluster": cluster,
                "n_common": int(n_common),
                "n_total_in_input": int(input_sizes.get(cluster, 0)),
                "microcluster_threshold": args.microcluster_threshold,
                "small_cluster_review": bool(n_common < args.microcluster_threshold),
                "included_in_all_cell_scores": True,
                "included_in_macro_restricted_scores": bool(n_common >= args.microcluster_threshold),
            })
    pd.DataFrame(audit_rows).to_csv(
        args.out / "clustering_microcluster_audit.tsv", sep="\t", index=False
    )

    sampled_ids = common_ids
    if len(sampled_ids) > args.max_cells:
        sampled_ids = pd.Index(np.random.default_rng(args.seed).choice(
            sampled_ids.to_numpy(), args.max_cells, replace=False
        ))

    pair_rows = []
    migration_rows = []
    for index, run_a in enumerate(runs):
        labels_a = full_labels[run_a].reindex(sampled_ids)
        for run_b in runs[index + 1:]:
            labels_b = full_labels[run_b].reindex(sampled_ids)
            all_ari, all_ami = safe_scores(labels_a.to_numpy(), labels_b.to_numpy())
            a_macro = labels_a.map(cluster_sizes[run_a]).ge(args.microcluster_threshold)
            b_macro = labels_b.map(cluster_sizes[run_b]).ge(args.microcluster_threshold)
            macro_mask = (a_macro & b_macro).to_numpy()
            macro_ari, macro_ami = safe_scores(
                labels_a.to_numpy()[macro_mask], labels_b.to_numpy()[macro_mask]
            )
            pair_rows.append({
                "run_a": run_a,
                "run_b": run_b,
                "n_common_sampled": len(sampled_ids),
                "ARI": all_ari,
                "AMI": all_ami,
                "n_common_sampled_all": len(sampled_ids),
                "ARI_all_observations": all_ari,
                "AMI_all_observations": all_ami,
                "n_common_sampled_macro_restricted": int(macro_mask.sum()),
                "ARI_macro_restricted": macro_ari,
                "AMI_macro_restricted": macro_ami,
                "microcluster_threshold": args.microcluster_threshold,
            })
            contingency = pd.crosstab(labels_a, labels_b, dropna=False)
            for cluster_a, counts in contingency.iterrows():
                for cluster_b, n_sampled in counts[counts > 0].items():
                    migration_rows.append({
                        "run_a": run_a,
                        "cluster_a": cluster_a,
                        "run_b": run_b,
                        "cluster_b": cluster_b,
                        "n_sampled": int(n_sampled),
                        "cluster_a_n_common": int(cluster_sizes[run_a][cluster_a]),
                        "cluster_b_n_common": int(cluster_sizes[run_b][cluster_b]),
                        "cluster_a_small_review": bool(cluster_sizes[run_a][cluster_a] < args.microcluster_threshold),
                        "cluster_b_small_review": bool(cluster_sizes[run_b][cluster_b] < args.microcluster_threshold),
                    })

    pair_columns = [
        "run_a", "run_b", "n_common_sampled", "ARI", "AMI",
        "n_common_sampled_all", "ARI_all_observations",
        "AMI_all_observations", "n_common_sampled_macro_restricted",
        "ARI_macro_restricted", "AMI_macro_restricted", "microcluster_threshold",
    ]
    pairs = pd.DataFrame(pair_rows, columns=pair_columns)
    pairs.to_csv(args.out / "clustering_pairwise_stability.tsv", sep="\t", index=False)
    pd.DataFrame(migration_rows).to_csv(
        args.out / "clustering_pairwise_migration.tsv", sep="\t", index=False
    )

    summary_rows = []
    for run in runs:
        subset = pairs[(pairs.run_a == run) | (pairs.run_b == run)]
        summary_rows.append({
            "run_id": run,
            "mean_ARI_shortlist": subset.ARI_all_observations.mean(),
            "min_ARI_shortlist": subset.ARI_all_observations.min(),
            "mean_AMI_shortlist": subset.AMI_all_observations.mean(),
            "mean_ARI_all_observations": subset.ARI_all_observations.mean(),
            "min_ARI_all_observations": subset.ARI_all_observations.min(),
            "mean_AMI_all_observations": subset.AMI_all_observations.mean(),
            "mean_ARI_macro_restricted": subset.ARI_macro_restricted.mean(),
            "min_ARI_macro_restricted": subset.ARI_macro_restricted.min(),
            "mean_AMI_macro_restricted": subset.AMI_macro_restricted.mean(),
            "comparisons": len(subset),
            "n_small_clusters_for_review": int((cluster_sizes[run] < args.microcluster_threshold).sum()),
        })
    summary = pd.DataFrame(summary_rows)
    if len(summary):
        summary = summary.sort_values(
            ["mean_ARI_macro_restricted", "mean_ARI_all_observations"],
            ascending=False, na_position="last"
        )
    summary.to_csv(args.out / "clustering_stability_summary.tsv", sep="\t", index=False)

    result = {
        "n_runs": len(runs),
        "n_common_observations": len(common_ids),
        "n_cells_sampled": len(sampled_ids),
        "microcluster_threshold": args.microcluster_threshold,
        "runs": runs,
        "all_observation_scores": "retain every sampled observation, including microclusters",
        "macro_restricted_scores": "resolution-ranking aid only; excludes cells in <threshold clusters in either partition",
        "small_cluster_policy": "flag for rare-lineage/technical review; never delete, relabel, or omit from DEG/spatial evidence by size alone",
        "warning": "Stability measures agreement, not biological correctness; inspect full-feature markers, anti-markers, rare-lineage and spatial evidence.",
    }
    (args.out / "clustering_stability.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
