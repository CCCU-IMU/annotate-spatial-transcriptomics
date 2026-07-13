#!/usr/bin/env python3
"""Select balanced positive markers per label from a standard DEG table."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("deg", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--label-col", default="cluster")
    ap.add_argument("--gene-col", default="gene")
    ap.add_argument("--lfc-col", default="avg_log2FC")
    ap.add_argument("--padj-col", default="p_val_adj")
    ap.add_argument("--min-lfc", type=float, default=0.25)
    ap.add_argument("--max-padj", type=float, default=0.05)
    ap.add_argument("--label-map", type=Path, help="TSV mapping source_cluster to a biological marker_group")
    ap.add_argument("--map-value-col", default="fine_label")
    args = ap.parse_args()
    with args.deg.open(newline="", encoding="utf-8") as handle:
        data = list(csv.DictReader(handle, delimiter="\t"))
    by = defaultdict(list)
    for r in data:
        try: lfc, padj = float(r[args.lfc_col]), float(r[args.padj_col])
        except (ValueError, KeyError): continue
        gene = r.get(args.gene_col, "")
        if gene and lfc >= args.min_lfc and padj <= args.max_padj:
            by[r[args.label_col]].append((padj, -lfc, gene))
    label_map = {}
    if args.label_map:
        with args.label_map.open(newline="", encoding="utf-8") as handle:
            for r in csv.DictReader(handle, delimiter="\t"):
                label_map[str(r["source_cluster"])] = r[args.map_value_col]
    seen = set(); rows = []
    for label in sorted(by):
        taken = 0
        for _, _, gene in sorted(by[label]):
            if gene in seen: continue
            seen.add(gene); rows.append({"gene": gene, "marker_group": label_map.get(str(label), str(label)), "panel": "data_specific"}); taken += 1
            if taken >= args.top: break
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "marker_group", "panel"], delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    print(f"wrote {len(rows)} markers to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
