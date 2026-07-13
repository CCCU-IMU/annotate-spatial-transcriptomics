#!/usr/bin/env python3
"""Compare clustering candidates on a deterministic cell sample using ARI/AMI."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("grid_summary", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--ranking", type=Path, help="Optional rank table used to shortlist")
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--max-cells", type=int, default=50000)
    ap.add_argument("--seed", type=int, default=20260713)
    args = ap.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    grid = pd.read_csv(args.grid_summary)
    if args.ranking:
        rank = pd.read_csv(args.ranking, sep="\t").sort_values("quantitative_rank").head(args.top)
        grid = grid[grid.run_id.isin(rank.run_id)]
    elif len(grid) > args.top:
        grid = grid.head(args.top)
    labels = {}; ids = None
    for r in grid.itertuples(index=False):
        d = pd.read_csv(r.csv_path, usecols=["cell_id", "cluster"], dtype={"cell_id": str, "cluster": str})
        if d.cell_id.duplicated().any(): raise RuntimeError(f"duplicate IDs in {r.run_id}")
        s = d.set_index("cell_id").cluster
        ids = s.index if ids is None else ids.intersection(s.index)
        labels[r.run_id] = s
    if ids is None or not len(ids): raise RuntimeError("no common cells")
    if len(ids) > args.max_cells:
        ids = pd.Index(np.random.default_rng(args.seed).choice(ids.to_numpy(), args.max_cells, replace=False))
    runs = list(labels); rows = []
    for i, a in enumerate(runs):
        la = labels[a].reindex(ids).to_numpy()
        for b in runs[i + 1:]:
            lb = labels[b].reindex(ids).to_numpy()
            rows.append({"run_a": a, "run_b": b, "n_common_sampled": len(ids),
                         "ARI": adjusted_rand_score(la, lb), "AMI": adjusted_mutual_info_score(la, lb)})
    pairs = pd.DataFrame(rows); pairs.to_csv(args.out / "clustering_pairwise_stability.tsv", sep="\t", index=False)
    summary_rows = []
    for run in runs:
        x = pairs[(pairs.run_a == run) | (pairs.run_b == run)]
        summary_rows.append({"run_id": run, "mean_ARI_shortlist": x.ARI.mean(), "min_ARI_shortlist": x.ARI.min(),
                             "mean_AMI_shortlist": x.AMI.mean(), "comparisons": len(x)})
    summary = pd.DataFrame(summary_rows).sort_values("mean_ARI_shortlist", ascending=False)
    summary.to_csv(args.out / "clustering_stability_summary.tsv", sep="\t", index=False)
    result = {"n_runs": len(runs), "n_cells_sampled": len(ids), "runs": runs,
              "warning": "Stability measures agreement, not biological correctness; inspect marker and spatial evidence."}
    (args.out / "clustering_stability.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2)); return 0


if __name__ == "__main__": raise SystemExit(main())

