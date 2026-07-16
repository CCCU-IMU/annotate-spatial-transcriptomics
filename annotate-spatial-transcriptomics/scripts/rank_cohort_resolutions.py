#!/usr/bin/env python3
"""Shortlist cohort resolutions without forcing a split or freezing a decision."""
from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score


QUESTION_DEFAULTS = {"broad_purity_audit": 1, "targeted_mixture": 2}


def marker_set(block: dict, anti: bool = False) -> set[str]:
    genes: set[str] = set()
    for key, value in block.items():
        is_anti = key.startswith("anti")
        if is_anti != anti or not isinstance(value, list):
            continue
        genes.update(map(str, value))
    return genes


def resolution_tag(path: Path) -> str:
    hit = re.search(r"res(?:olution)?([0-9]+p?[0-9]*)", path.stem, re.I)
    return hit.group(1).replace("p", ".") if hit else path.stem


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster-glob", required=True)
    parser.add_argument("--deg-glob", required=True)
    parser.add_argument("--profile", required=True, type=Path)
    parser.add_argument("--lineages", required=True, help="comma-separated profile lineage keys")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--question-mode", required=True,
        choices=sorted(QUESTION_DEFAULTS),
        help="broad purity audits permit one homogeneous compartment; targeted mixtures may expect competing compartments",
    )
    parser.add_argument(
        "--expected-min-compartments", type=int,
        help="defaults to 1 for broad_purity_audit and 2 for targeted_mixture",
    )
    parser.add_argument("--expected-grid", default="", help="optional comma-separated contract grid")
    parser.add_argument("--composition", type=Path)
    parser.add_argument("--anchor-summary", type=Path)
    parser.add_argument("--qc-summary", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    expected = args.expected_min_compartments
    if expected is None:
        expected = QUESTION_DEFAULTS[args.question_mode]
    if expected < 1:
        raise SystemExit("expected-min-compartments must be >=1")
    if args.question_mode == "broad_purity_audit" and expected != 1:
        raise SystemExit("broad_purity_audit must use expected_min_compartments=1; a broad cohort cannot be forced to split")

    profile = json.loads(args.profile.read_text(encoding="utf-8"))
    lineage_keys = [item.strip() for item in args.lineages.split(",") if item.strip()]
    lineages = profile.get("lineages", {})
    positive = {key: marker_set(lineages.get(key, {})) for key in lineage_keys}
    negative = {key: marker_set(lineages.get(key, {}), anti=True) for key in lineage_keys}
    if not any(positive.values()):
        raise SystemExit("No profile markers found for requested lineages")

    cluster_paths = sorted(map(Path, glob.glob(args.cluster_glob)))
    deg_paths = {resolution_tag(Path(path)): Path(path) for path in glob.glob(args.deg_glob)}
    memberships: dict[str, pd.DataFrame] = {}
    rows: list[dict] = []
    for cluster_path in cluster_paths:
        key = resolution_tag(cluster_path)
        clusters = read_table(cluster_path)
        if not {"cell_id", "cluster"}.issubset(clusters.columns):
            continue
        clusters = clusters[["cell_id", "cluster"]].drop_duplicates("cell_id")
        memberships[key] = clusters
        sizes = clusters["cluster"].value_counts()
        n_observations = len(clusters)
        tiny_cut = max(20, int(np.ceil(n_observations * 0.001)))
        coverage: dict[str, float] = {}
        explainable = 0
        anti_rates: list[float] = []
        deg_path = deg_paths.get(key)
        if deg_path and deg_path.exists():
            deg = pd.read_csv(deg_path, sep="\t")
            gene_col = "gene" if "gene" in deg else "names" if "names" in deg else None
            cluster_col = "cluster" if "cluster" in deg else "group" if "group" in deg else None
            if gene_col and cluster_col:
                deg[gene_col] = deg[gene_col].astype(str)
                deg[cluster_col] = deg[cluster_col].astype(str)
                top = deg.groupby(cluster_col, group_keys=False).head(100)
                for lineage, genes in positive.items():
                    denom = max(1, min(5, len(genes)))
                    values = top.groupby(cluster_col)[gene_col].apply(lambda series: len(set(series) & genes) / denom)
                    coverage[lineage] = float(min(1, values.max())) if len(values) else 0.0
                for _, group in top.groupby(cluster_col):
                    genes = set(group[gene_col])
                    hits = {lineage: len(genes & markers) for lineage, markers in positive.items()}
                    best = max(hits, key=hits.get) if hits else None
                    if best is not None and hits[best] >= 2:
                        anti = len(genes & negative.get(best, set()))
                        anti_rates.append(anti / max(1, hits[best] + anti))
                        if anti <= hits[best]:
                            explainable += 1
        n_clusters = int(len(sizes))
        rows.append({
            "resolution": key,
            "cluster_file": str(cluster_path.resolve()),
            "deg_file": str(deg_path.resolve()) if deg_path else "",
            "n_observations": n_observations,
            "n_clusters": n_clusters,
            "min_cluster": int(sizes.min()),
            "median_cluster": float(sizes.median()),
            "tiny_observation_fraction": float(sizes[sizes < tiny_cut].sum() / n_observations),
            "tiny_cluster_fraction": float((sizes < tiny_cut).mean()),
            "program_coverage": float(np.mean(list(coverage.values()))) if coverage else 0.0,
            "explainable_cluster_fraction": explainable / max(1, n_clusters),
            "anti_marker_rate": float(np.mean(anti_rates)) if anti_rates else 1.0,
            **{f"coverage_{lineage}": value for lineage, value in coverage.items()},
        })

    ranking = pd.DataFrame(rows)

    def merge_external(path: Path | None, kind: str) -> None:
        nonlocal ranking
        if path is None or not path.exists() or ranking.empty:
            return
        frame = pd.read_csv(path, sep="\t")
        frame["resolution"] = frame["resolution"].astype(str).str.replace("p", ".", regex=False)
        if kind == "composition" and {"field", "fraction", "cluster"}.issubset(frame):
            metric = frame[frame.field == "source_key"].groupby(["resolution", "cluster"]).fraction.max().groupby("resolution").mean().rename("source_dominance").reset_index()
        elif kind == "anchor" and {"fraction", "cluster"}.issubset(frame):
            metric = frame.groupby(["resolution", "cluster"]).fraction.max().groupby("resolution").mean().rename("anchor_purity").reset_index()
        elif kind == "qc" and "cluster" in frame:
            values = []
            metrics = [column for column in frame if column not in {"resolution", "cluster"}]
            for resolution, group in frame.groupby("resolution"):
                cvs = []
                for metric_name in metrics:
                    value = pd.to_numeric(group[metric_name], errors="coerce")
                    mean = value.mean()
                    cvs.append(float(value.std() / mean) if np.isfinite(mean) and mean > 0 else 0.0)
                values.append({"resolution": resolution, "qc_depth_cv": float(np.nanmean(cvs))})
            metric = pd.DataFrame(values)
        else:
            return
        ranking = ranking.merge(metric, on="resolution", how="left")

    merge_external(args.composition, "composition")
    merge_external(args.anchor_summary, "anchor")
    merge_external(args.qc_summary, "qc")
    order = sorted(memberships, key=lambda value: float(value) if re.fullmatch(r"[0-9.]+", value) else 1e9)
    adjacent_ari: dict[str, float] = {}
    for left, right in zip(order, order[1:]):
        merged = memberships[left].merge(memberships[right], on="cell_id", suffixes=("_left", "_right"))
        adjacent_ari[right] = adjusted_rand_score(merged.cluster_left, merged.cluster_right) if len(merged) else np.nan

    if not ranking.empty:
        ranking["adjacent_ari"] = ranking.resolution.map(adjacent_ari)
        ranking["under_split_penalty"] = (ranking.n_clusters < expected).astype(float)
        if args.question_mode == "broad_purity_audit":
            ranking["under_split_penalty"] = 0.0
        for column, default in (("source_dominance", 0.5), ("anchor_purity", 0.5), ("qc_depth_cv", 0.0)):
            if column not in ranking:
                ranking[column] = default
            ranking[column] = ranking[column].fillna(default)
        stability = ranking.adjacent_ari.fillna(ranking.adjacent_ari.median() if ranking.adjacent_ari.notna().any() else 0.5)
        ranking["screening_score"] = (
            0.30 * ranking.program_coverage + 0.18 * ranking.explainable_cluster_fraction
            + 0.20 * stability - 0.12 * ranking.anti_marker_rate + 0.08 * ranking.anchor_purity
            - 0.08 * ranking.source_dominance - 0.08 * ranking.qc_depth_cv.clip(upper=2)
            - 0.12 * ranking.tiny_observation_fraction - 0.08 * ranking.tiny_cluster_fraction
            - 0.25 * ranking.under_split_penalty
        )
        ranking["quantitative_rank"] = ranking.screening_score.rank(ascending=False, method="min").astype(int)
        ranking = ranking.sort_values(["quantitative_rank", "n_clusters"])

    ranking.to_csv(args.out / "cohort_resolution_ranking.tsv", sep="\t", index=False)
    observed_grid = sorted(float(value) for value in order if re.fullmatch(r"[0-9.]+", value))
    expected_grid = sorted(float(value) for value in args.expected_grid.split(",") if value.strip())
    complete_grid = not expected_grid or observed_grid == expected_grid
    summary = {
        "status": "SHORTLIST_ONLY" if len(ranking) and complete_grid else "INCOMPLETE_GRID",
        "question_mode": args.question_mode,
        "expected_min_compartments": expected,
        "under_split_penalty_applied": args.question_mode == "targeted_mixture",
        "observed_grid": observed_grid,
        "expected_grid": expected_grid,
        "complete_grid": complete_grid,
        "lineages": lineage_keys,
        "n_candidates": len(ranking),
        "top_candidates": ranking.head(3).resolution.tolist() if len(ranking) else [],
        "allowed_terminal_outcomes": ["homogeneous_parent_confirmed", "subclusters_adjudicated"] if args.question_mode == "broad_purity_audit" else ["subclusters_adjudicated"],
        "homogeneous_parent_confirmed_is_success": args.question_mode == "broad_purity_audit",
        "selection_rule": "Select the resolution with the best integrated evidence for the current biological question. Preserve stable lineages or true subtypes first, avoid state/technical fragmentation second, and use lower complexity only as a tie-breaker.",
        "warning": "The numeric rank only shortlists candidates. A one-cluster broad cohort is not penalized and may close successfully as homogeneous_parent_confirmed after evidence review.",
    }
    (args.out / "cohort_resolution_ranking.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["status"] == "SHORTLIST_ONLY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
