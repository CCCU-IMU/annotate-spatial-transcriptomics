#!/usr/bin/env python3
"""Summarize full-transcriptome pseudobulk for broad-level review questions."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from evidence_schema_lib import sha256


TECHNICAL = re.compile(r"^(?:RPL|RPS|MRPL|MRPS|MT-|MTRNR)")


def artifact(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "sha256": sha256(path), "n_bytes": path.stat().st_size}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pseudobulk", required=True, type=Path)
    ap.add_argument("--marker-manifest", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    table = pd.read_csv(args.pseudobulk, sep="\t", low_memory=False)
    markers = pd.read_csv(args.marker_manifest, sep="\t", dtype=str).fillna("")
    required = {"gene", "group_id", "n_observations", "cpm", "detection_fraction"}
    if not required.issubset(table):
        raise SystemExit("pseudobulk table lacks required columns")
    current_groups = sorted(group for group in table.group_id.unique() if group.startswith("current::"))
    recall_groups = sorted(group for group in table.group_id.unique() if group.startswith("recall_question::"))
    if not current_groups or not recall_groups:
        raise SystemExit("pseudobulk requires current and recall-question groups")
    logcpm = table.pivot(index="gene", columns="group_id", values="cpm").fillna(0)
    logcpm = np.log1p(logcpm)
    current_matrix = logcpm[current_groups]
    eligible_genes = np.array([
        gene for gene in current_matrix.index
        if not TECHNICAL.match(str(gene)) and current_matrix.loc[gene].max() > np.log1p(1)
    ], dtype=object)
    variance = current_matrix.loc[eligible_genes].var(axis=1).sort_values(ascending=False)
    informative = variance.head(min(2000, len(variance))).index
    summary_rows: list[dict[str, object]] = []
    similarity_rows: list[dict[str, object]] = []
    differential_rows: list[dict[str, object]] = []
    table_index = table.set_index(["group_id", "gene"])
    for group in recall_groups:
        pieces = group.split("::")
        target = pieces[1]
        origin = pieces[3] if len(pieces) >= 4 and pieces[2] == "from" else ""
        current_target = f"current::{target}"
        if current_target not in logcpm:
            continue
        recall_vector = logcpm.loc[informative, group]
        similarities: list[tuple[str, float]] = []
        for current_group in current_groups:
            value = spearmanr(
                recall_vector.to_numpy(),
                logcpm.loc[informative, current_group].to_numpy(),
            ).correlation
            if not np.isfinite(value):
                value = 0.0
            label = current_group[len("current::"):]
            similarities.append((label, float(value)))
            similarity_rows.append({
                "recall_group": group, "target_broad_label": target,
                "origin_broad_label": origin, "current_reference_broad": label,
                "spearman_top_variable_genes": float(value),
            })
        similarities.sort(key=lambda item: (-item[1], item[0]))
        target_marker_genes = sorted(set(markers.loc[
            markers.broad_label.eq(target)
            & markers.evidence_role.eq("positive_family"), "gene"
        ]) & set(logcpm.index))
        target_marker_cpm = table.loc[
            table.group_id.eq(group) & table.gene.isin(target_marker_genes)
        ]
        n_observations = int(table.loc[table.group_id.eq(group), "n_observations"].iloc[0])
        summary_rows.append({
            "recall_group": group, "target_broad_label": target,
            "origin_broad_label": origin, "n_observations": n_observations,
            "nearest_current_broad": similarities[0][0],
            "nearest_current_spearman": similarities[0][1],
            "target_current_spearman": dict(similarities).get(target, 0.0),
            "target_positive_marker_n_available": len(target_marker_genes),
            "target_positive_marker_n_detected": int((target_marker_cpm.cpm > 0).sum()),
            "target_positive_marker_median_log1p_cpm": float(np.log1p(target_marker_cpm.cpm).median()) if len(target_marker_cpm) else 0.0,
            "target_positive_marker_median_detection_fraction": float(target_marker_cpm.detection_fraction.median()) if len(target_marker_cpm) else 0.0,
        })
        genes = logcpm.index
        recall_raw = table_index.loc[group].reindex(genes).fillna(0)
        target_raw = table_index.loc[current_target].reindex(genes).fillna(0)
        log2fc_target = np.log2((recall_raw.cpm.to_numpy() + 0.1) / (target_raw.cpm.to_numpy() + 0.1))
        if origin and f"current::{origin}" in logcpm:
            origin_raw = table_index.loc[f"current::{origin}"].reindex(genes).fillna(0)
            log2fc_origin = np.log2((recall_raw.cpm.to_numpy() + 0.1) / (origin_raw.cpm.to_numpy() + 0.1))
        else:
            log2fc_origin = np.full(len(genes), np.nan)
        order = np.argsort(-np.abs(log2fc_target))[:100]
        for index in order:
            gene = str(genes[index])
            differential_rows.append({
                "recall_group": group, "target_broad_label": target,
                "origin_broad_label": origin, "gene": gene,
                "recall_cpm": float(recall_raw.cpm.iloc[index]),
                "recall_detection_fraction": float(recall_raw.detection_fraction.iloc[index]),
                "log2fc_vs_current_target": float(log2fc_target[index]),
                "log2fc_vs_current_origin": float(log2fc_origin[index]),
            })
    args.out.mkdir(parents=True, exist_ok=True)
    summary_path = args.out / "broad_cell_type_recall_pseudobulk_summary.tsv"
    similarity_path = args.out / "broad_cell_type_recall_similarity.tsv"
    differential_path = args.out / "broad_cell_type_recall_top_differential.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    pd.DataFrame(similarity_rows).to_csv(similarity_path, sep="\t", index=False)
    pd.DataFrame(differential_rows).to_csv(differential_path, sep="\t", index=False)
    manifest = {
        "schema_version": "2.2", "status": "PASS_EVIDENCE_ONLY",
        "artifact_role": "broad_cell_type_full_transcriptome_pseudobulk_summary",
        "formal_membership_written": False,
        "pseudobulk": artifact(args.pseudobulk),
        "marker_manifest": artifact(args.marker_manifest),
        "informative_gene_n": len(informative),
        "review_group_n": len(summary_rows),
        "artifacts": {
            "summary": artifact(summary_path), "similarity": artifact(similarity_path),
            "top_differential": artifact(differential_path),
        },
        "warning": "Similarity and differential evidence inform one broad-level biological review; they do not assign cell labels.",
    }
    manifest_path = args.out / "broad_cell_type_review_pseudobulk_summary_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
