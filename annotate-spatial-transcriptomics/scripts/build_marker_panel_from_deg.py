#!/usr/bin/env python3
"""Build a current-label canonical and data-specific marker panel from strict DEG tables."""
from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: str, default: float) -> float:
    try:
        result = float(value)
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--broad-deg", required=True, type=Path)
    parser.add_argument("--subtype-deg", required=True, type=Path)
    parser.add_argument("--canonical", required=True, type=Path)
    parser.add_argument("--group-alias", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--top-per-label", type=int, default=4)
    parser.add_argument("--min-pct", type=float, default=0.05)
    parser.add_argument("--max-padj", type=float, default=0.05)
    parser.add_argument(
        "--exclude-regex", default=r"^(LOC|ENSOARG|MT-|RPS|RPL|HBA|HBB|EEF1|GAPDH$|ACTB$|MALAT1$|GNAS$)",
    )
    args = parser.parse_args()

    deg = {"broad": read(args.broad_deg), "subtype": read(args.subtype_deg)}
    for level, rows in deg.items():
        if not rows:
            raise SystemExit(f"{level} DEG is empty")
        required = {"label", "gene", "analysis_view"}
        if not required.issubset(rows[0]):
            raise SystemExit(f"{level} DEG lacks columns: {sorted(required-set(rows[0]))}")
        if any(row.get("analysis_view") != "strict" for row in rows):
            raise SystemExit(f"{level} DEG is not entirely strict evidence")

    label_sets = {
        level: {row.get("label", "") for row in rows if row.get("label", "")}
        for level, rows in deg.items()
    }
    aliases: dict[tuple[str, str], list[str]] = defaultdict(list)
    if args.group_alias:
        for row in read(args.group_alias):
            level, source, target = row.get("level", ""), row.get("source_group", ""), row.get("target_group", "")
            if level not in {"broad", "subtype"} or not source or not target:
                raise SystemExit("Alias rows require level, source_group and target_group")
            aliases[(level, source)].append(target)

    output: list[dict[str, str]] = []
    canonical = read(args.canonical)
    for row in canonical:
        level = row.get("level", "")
        if level not in {"broad", "subtype"} or row.get("panel") != "canonical":
            continue
        source = row.get("marker_group", "")
        targets = aliases.get((level, source), [source])
        for target in targets:
            if target == "__DROP__" or target not in label_sets[level]:
                continue
            output.append({
                "gene": row.get("gene", ""), "marker_group": target,
                "panel": "canonical", "level": level,
            })

    excluded = re.compile(args.exclude_regex, re.IGNORECASE)
    for level, rows in deg.items():
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            gene, label = row.get("gene", ""), row.get("label", "")
            if not gene or not label or excluded.search(gene):
                continue
            padj = number(row.get("p_val_adj", row.get("padj", "")), 1.0)
            logfc = number(row.get("avg_log2FC", row.get("avg_logFC", "")), -math.inf)
            pct1 = number(row.get("pct.1", row.get("pct_in", "")), 0.0)
            pct2 = number(row.get("pct.2", row.get("pct_out", "")), 0.0)
            if padj > args.max_padj or logfc <= 0 or pct1 < args.min_pct:
                continue
            # Balance effect size, contrast and detection. Pure P-value ranking
            # promotes ubiquitous genes in large groups; pure logFC ranking
            # promotes very sparse idiosyncratic genes.
            contrast = max(pct1 - pct2, 0.0)
            marker_score = logfc * max(contrast, 0.01) * math.sqrt(max(pct1, 0.0))
            item = dict(row); item["__rank"] = (-marker_score, -logfc, -contrast, padj, -pct1, gene)
            grouped[label].append(item)
        for label in sorted(label_sets[level]):
            candidates = sorted(grouped.get(label, []), key=lambda row: row["__rank"])
            selected = []
            seen = set()
            for row in candidates:
                gene = row["gene"]
                if gene in seen:
                    continue
                selected.append(gene); seen.add(gene)
                if len(selected) >= args.top_per_label:
                    break
            if not selected:
                raise SystemExit(f"No data-specific marker passed filters for {level} label {label}")
            output.extend(
                {"gene": gene, "marker_group": label, "panel": "data_specific", "level": level}
                for gene in selected
            )

    deduplicated = []
    seen_rows = set()
    for row in output:
        key = (row["gene"], row["marker_group"], row["panel"], row["level"])
        if not row["gene"] or key in seen_rows:
            continue
        seen_rows.add(key); deduplicated.append(row)
    for level in ["broad", "subtype"]:
        for panel in ["canonical", "data_specific"]:
            covered = {row["marker_group"] for row in deduplicated if row["level"] == level and row["panel"] == panel}
            missing = label_sets[level] - covered
            if missing:
                raise SystemExit(f"{level}/{panel} panel lacks current labels: {sorted(missing)}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "marker_group", "panel", "level"], delimiter="\t")
        writer.writeheader(); writer.writerows(deduplicated)
    print(
        f"PASS rows={len(deduplicated)} broad_labels={len(label_sets['broad'])} "
        f"subtype_labels={len(label_sets['subtype'])} out={args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
