#!/usr/bin/env python3
"""Rank BANKSY grid candidates for review; never auto-freeze a selection."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def entropy(counts: list[int]) -> float:
    total = sum(counts)
    if total <= 0 or len(counts) <= 1:
        return 0.0
    h = -sum((n / total) * math.log(n / total) for n in counts if n)
    return h / math.log(len(counts))


def load_counts(path: Path) -> dict[str, list[int]]:
    result: dict[str, list[int]] = defaultdict(list)
    if not path.exists():
        return result
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            run = row.get("run_id", "")
            n = row.get("n_cells") or row.get("n_observations") or "0"
            result[run].append(int(float(n)))
    return result


def load_broad_evidence(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None or not path.is_file():
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        return {row.get("run_id", ""): row for row in csv.DictReader(handle, delimiter="\t")}


def value(row: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, "") or default)
    except ValueError:
        return default


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("grid_summary", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--target-clusters", type=int, default=20, help="Review target, not a forced selection")
    ap.add_argument("--tiny-threshold", type=int, default=100)
    ap.add_argument("--counts-table", type=Path)
    ap.add_argument("--broad-evidence-table", type=Path, help="Per-run query-derived broad recall/purity evidence TSV")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    with args.grid_summary.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    counts_path = args.counts_table
    if counts_path is None:
        candidates = list(args.grid_summary.parent.parent.glob("plots*/*cluster_counts_by_run.csv"))
        counts_path = candidates[0] if candidates else Path("__missing__")
    counts = load_counts(counts_path)
    broad_evidence = load_broad_evidence(args.broad_evidence_table)
    out = []
    for row in rows:
        run = row["run_id"]
        ncl = int(row["n_clusters"])
        sizes = counts.get(run, [])
        total = sum(sizes)
        tiny_fraction = sum(n for n in sizes if n < args.tiny_threshold) / total if total else float("nan")
        balance = entropy(sizes) if sizes else float("nan")
        count_fit = math.exp(-abs(math.log(max(ncl, 1) / args.target_clusters)))
        tiny_score = 1 - min(tiny_fraction, 1) if math.isfinite(tiny_fraction) else 0.5
        balance_score = balance if math.isfinite(balance) else 0.5
        evidence = broad_evidence.get(run, {})
        evidence_complete = all(key in evidence and str(evidence.get(key, "")).strip() for key in (
            "default_broad_recall", "marker_coherence", "spatial_coherence", "adjacent_stability",
            "large_cluster_purity", "technical_fragmentation", "zero_census_default",
        ))
        if evidence_complete:
            biology_score = (
                0.25 * value(evidence, "default_broad_recall")
                + 0.18 * value(evidence, "marker_coherence")
                + 0.16 * value(evidence, "spatial_coherence")
                + 0.13 * value(evidence, "adjacent_stability")
                + 0.18 * value(evidence, "large_cluster_purity")
                + 0.10 * (1 - min(max(value(evidence, "technical_fragmentation"), 0), 1))
            )
            zero_penalty = min(max(value(evidence, "zero_census_default"), 0), 1)
            score = 0.85 * biology_score + 0.10 * tiny_score + 0.05 * balance_score - 0.20 * zero_penalty
        else:
            biology_score = float("nan")
            score = 0.10 * count_fit + 0.45 * tiny_score + 0.45 * balance_score
        out.append({**row, "total_observations": total or "", "tiny_fraction": tiny_fraction,
                    "normalized_size_entropy": balance, "cluster_count_fit": count_fit,
                    "broad_biology_score": biology_score, "broad_evidence_complete": str(evidence_complete).lower(),
                    "selection_eligible": str(evidence_complete).lower(),
                    "screening_score": score, "selection_status": "shortlist_only",
                    "required_manual_review": "full_catalog_recall+large_cluster_purity+zero_census+markers+spatial+adjacent_resolution_migration"})
    out.sort(key=lambda r: float(r["screening_score"]), reverse=True)
    for i, row in enumerate(out, 1): row["quantitative_rank"] = i
    dest = args.out / "banksy_candidate_ranking.tsv"
    fields = list(out[0]) if out else []
    with dest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(out)
    summary = {
        "grid_summary": str(args.grid_summary.resolve()), "counts_table": str(counts_path),
        "target_clusters_review_value": args.target_clusters, "n_candidates": len(out),
        "broad_evidence_table": str(args.broad_evidence_table.resolve()) if args.broad_evidence_table else "",
        "selection_eligible_candidates": sum(row["selection_eligible"] == "true" for row in out),
        "top_shortlist": [{k: r[k] for k in ["run_id", "resolution", "k_neighbors", "n_clusters", "quantitative_rank"]} for r in out[:6]],
        "warning": "This score is a shortlist only. BANKSY broad selection requires complete query-derived lineage recall, large-label purity, zero-census, marker, spatial and adjacent-resolution review; cluster count never freezes a resolution.",
    }
    (args.out / "banksy_candidate_ranking.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
