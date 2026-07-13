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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("grid_summary", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--target-clusters", type=int, default=20, help="Review target, not a forced selection")
    ap.add_argument("--tiny-threshold", type=int, default=100)
    ap.add_argument("--counts-table", type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    with args.grid_summary.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    counts_path = args.counts_table
    if counts_path is None:
        candidates = list(args.grid_summary.parent.parent.glob("plots*/*cluster_counts_by_run.csv"))
        counts_path = candidates[0] if candidates else Path("__missing__")
    counts = load_counts(counts_path)
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
        score = 0.55 * count_fit + 0.25 * tiny_score + 0.20 * balance_score
        out.append({**row, "total_observations": total or "", "tiny_fraction": tiny_fraction,
                    "normalized_size_entropy": balance, "cluster_count_fit": count_fit,
                    "screening_score": score, "selection_status": "shortlist_only",
                    "required_manual_review": "markers+spatial+adjacent_parameter_stability"})
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
        "top_shortlist": [{k: r[k] for k in ["run_id", "resolution", "k_neighbors", "n_clusters", "quantitative_rank"]} for r in out[:6]],
        "warning": "This score is a shortlist only. An agent must review markers, space, stability and morphology before selection.",
    }
    (args.out / "banksy_candidate_ranking.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

