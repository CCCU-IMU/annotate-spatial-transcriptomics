#!/usr/bin/env python3
"""Rank, but never automatically freeze, pool-clustering resolutions."""
from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score


def markers(block: dict) -> set[str]:
    out: set[str] = set()
    for key, value in block.items():
        if key.startswith("anti") or key in {"support_only", "safety", "spatial_expectation"}: continue
        if isinstance(value, list): out.update(map(str, value))
    return out


def anti_markers(block: dict) -> set[str]:
    out: set[str] = set()
    for key, value in block.items():
        if key.startswith("anti") and isinstance(value, list): out.update(map(str, value))
    return out


def tag(path: Path) -> str:
    hit = re.search(r"res(?:olution)?([0-9]+p?[0-9]*)", path.stem, re.I)
    return hit.group(1).replace("p", ".") if hit else path.stem


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cluster-glob", required=True)
    ap.add_argument("--deg-glob", required=True)
    ap.add_argument("--profile", required=True, type=Path)
    ap.add_argument("--lineages", required=True, help="comma-separated profile lineage keys")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--expected-min-compartments", type=int, default=2)
    ap.add_argument("--composition", type=Path, help="cluster_source_state_composition.tsv")
    ap.add_argument("--anchor-summary", type=Path, help="cluster_anchor_distance_summary.tsv")
    ap.add_argument("--qc-summary", type=Path, help="cluster_QC_summary.tsv")
    args = ap.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    profile = json.loads(args.profile.read_text()); lineage_keys = [x.strip() for x in args.lineages.split(",") if x.strip()]
    lineage_markers = {k: markers(profile.get("lineages", {}).get(k, {})) for k in lineage_keys}
    lineage_anti = {k: anti_markers(profile.get("lineages", {}).get(k, {})) for k in lineage_keys}
    if not any(lineage_markers.values()): raise SystemExit("No profile markers found for requested lineages")
    cluster_paths = sorted(map(Path, glob.glob(args.cluster_glob)))
    deg_paths = {tag(Path(p)): Path(p) for p in glob.glob(args.deg_glob)}
    tables: dict[str, pd.DataFrame] = {}; rows = []
    for cp in cluster_paths:
        key = tag(cp); c = pd.read_csv(cp, sep=None, engine="python", dtype={"cell_id": str, "cluster": str})
        if not {"cell_id", "cluster"}.issubset(c): continue
        tables[key] = c[["cell_id", "cluster"]].drop_duplicates("cell_id")
        sizes = c.cluster.value_counts(); n = len(c); tiny_cut = max(20, int(np.ceil(n * .001)))
        tiny_fraction = float(sizes[sizes < tiny_cut].sum() / n); singleton_clusters = float((sizes < tiny_cut).mean())
        coverage = {}; explainable = 0; anti_rates=[]
        dp = deg_paths.get(key)
        if dp and dp.exists():
            d = pd.read_csv(dp, sep="\t"); gene_col = "gene" if "gene" in d else "names" if "names" in d else None
            cl_col = "cluster" if "cluster" in d else "group" if "group" in d else None
            if gene_col and cl_col:
                d[gene_col] = d[gene_col].astype(str); d[cl_col] = d[cl_col].astype(str)
                top = d.groupby(cl_col, group_keys=False).head(100)
                for lineage, genes in lineage_markers.items():
                    denom = max(1, min(5, len(genes)))
                    by = top.groupby(cl_col,group_keys=False)[gene_col].apply(lambda x: len(set(x) & genes) / denom)
                    coverage[lineage] = float(min(1, by.max())) if len(by) else 0.0
                for _, g in top.groupby(cl_col):
                    genes=set(g[gene_col]);hits={k:len(genes & v) for k,v in lineage_markers.items()};best=max(hits,key=hits.get) if hits else None
                    if best is not None and hits[best] >= 2:
                        anti=len(genes & lineage_anti.get(best,set()));anti_rates.append(anti/max(1,hits[best]+anti));
                        if anti <= hits[best]: explainable += 1
        n_clusters = int(sizes.size); explainable_fraction = explainable / max(1, n_clusters)
        program_coverage = float(np.mean(list(coverage.values()))) if coverage else 0.0
        rows.append({"resolution": key, "cluster_file": str(cp.resolve()), "deg_file": str(dp.resolve()) if dp else "",
                     "n_observations": n, "n_clusters": n_clusters, "min_cluster": int(sizes.min()),
                     "median_cluster": float(sizes.median()), "tiny_observation_fraction": tiny_fraction,
                     "tiny_cluster_fraction": singleton_clusters, "program_coverage": program_coverage,
                     "explainable_cluster_fraction": explainable_fraction,"anti_marker_rate":float(np.mean(anti_rates)) if anti_rates else 1.0,
                     **{f"coverage_{k}": v for k, v in coverage.items()}})
    out = pd.DataFrame(rows)
    def add_external(path, kind: str):
        nonlocal out
        if path is None or not path.exists() or out.empty:return
        frame=pd.read_csv(path,sep="\t");frame["resolution"]=frame["resolution"].astype(str).str.replace("p",".",regex=False)
        if kind=="composition" and {"field","fraction","cluster"}.issubset(frame):
            x=frame[frame.field=="source_key"].groupby(["resolution","cluster"]).fraction.max().groupby("resolution").mean().rename("source_dominance").reset_index();out=out.merge(x,on="resolution",how="left")
        elif kind=="anchor" and {"fraction","cluster"}.issubset(frame):
            x=frame.groupby(["resolution","cluster"]).fraction.max().groupby("resolution").mean().rename("anchor_purity").reset_index();out=out.merge(x,on="resolution",how="left")
        elif kind=="qc" and "cluster" in frame:
            metrics=[x for x in frame if x not in {"resolution","cluster"}]
            if metrics:
                values=[]
                for res,g in frame.groupby("resolution"):
                    cvs=[]
                    for metric in metrics:
                        v=pd.to_numeric(g[metric],errors="coerce");m=v.mean();cvs.append(float(v.std()/m) if np.isfinite(m) and m>0 else 0)
                    values.append({"resolution":res,"qc_depth_cv":float(np.nanmean(cvs))})
                out=out.merge(pd.DataFrame(values),on="resolution",how="left")
    add_external(args.composition,"composition");add_external(args.anchor_summary,"anchor");add_external(args.qc_summary,"qc")
    stability = {}
    order = sorted(tables, key=lambda x: float(x) if re.fullmatch(r"[0-9.]+", x) else 1e9)
    for left, right in zip(order, order[1:]):
        m = tables[left].merge(tables[right], on="cell_id", suffixes=("_l", "_r"))
        stability[right] = adjusted_rand_score(m.cluster_l, m.cluster_r) if len(m) else np.nan
    if len(out):
        out["adjacent_ari"] = out.resolution.map(stability)
        out["under_split_penalty"] = (out.n_clusters < args.expected_min_compartments).astype(float)
        for col,default in [("source_dominance",0.5),("anchor_purity",0.5),("qc_depth_cv",0.0)]:
            if col not in out:out[col]=default
            out[col]=out[col].fillna(default)
        out["screening_score"] = (0.30*out.program_coverage + 0.18*out.explainable_cluster_fraction +
                                  0.20*out.adjacent_ari.fillna(out.adjacent_ari.median() if out.adjacent_ari.notna().any() else .5) -
                                  0.12*out.anti_marker_rate + 0.08*out.anchor_purity - 0.08*out.source_dominance -
                                  0.08*out.qc_depth_cv.clip(upper=2) - 0.12*out.tiny_observation_fraction - 0.08*out.tiny_cluster_fraction -
                                  0.25*out.under_split_penalty)
        out["quantitative_rank"] = out.screening_score.rank(ascending=False, method="min").astype(int)
        out = out.sort_values(["quantitative_rank", "n_clusters"])
    out.to_csv(args.out/"pool_resolution_ranking.tsv", sep="\t", index=False)
    summary = {"status": "SHORTLIST_ONLY", "lineages": lineage_keys, "n_candidates": len(out),
               "top_candidates": out.head(3).resolution.tolist() if len(out) else [],
               "mandatory_manual_channels": ["spatial morphology", "anti-marker rate", "source/state/QC composition", "anchor-distance evidence", "cluster migration", "new lineage versus state split"],
               "warning": "The numeric rank cannot freeze a resolution; write an evidence-backed decision to the clustering ledger."}
    (args.out/"pool_resolution_ranking.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2)); return 0 if len(out) else 2


if __name__ == "__main__": raise SystemExit(main())
