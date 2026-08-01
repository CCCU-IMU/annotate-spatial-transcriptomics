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
    ap.add_argument("--pseudobulk-manifest", required=True, type=Path)
    ap.add_argument("--broad-evidence-manifest", required=True, type=Path)
    ap.add_argument("--marker-manifest", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    pseudobulk_manifest = json.loads(
        args.pseudobulk_manifest.read_text(encoding="utf-8")
    )
    if (
        pseudobulk_manifest.get("artifact_role")
        != "broad_cell_type_full_transcriptome_pseudobulk"
        or pseudobulk_manifest.get("assay_ancestry")
        != "project_local_non_SCT_raw_counts"
        or Path(str(pseudobulk_manifest.get("pseudobulk", ""))).resolve()
        != args.pseudobulk.resolve()
    ):
        raise SystemExit("pseudobulk summary input lacks project-local raw-count ancestry")
    broad_evidence = json.loads(
        args.broad_evidence_manifest.read_text(encoding="utf-8")
    )
    recall_record = broad_evidence.get("artifacts", {}).get(
        "recall_membership", {}
    )
    if (
        broad_evidence.get("artifact_role")
        != "broad_cell_type_targeted_review_evidence"
        or Path(str(pseudobulk_manifest.get("source_membership", ""))).resolve()
        != Path(str(broad_evidence.get("membership", {}).get("path", ""))).resolve()
        or Path(str(pseudobulk_manifest.get("recall_membership", ""))).resolve()
        != Path(str(recall_record.get("path", ""))).resolve()
    ):
        raise SystemExit("pseudobulk inputs differ from the broad cell-type evidence round")
    table = pd.read_csv(args.pseudobulk, sep="\t", low_memory=False)
    markers = pd.read_csv(args.marker_manifest, sep="\t", dtype=str).fillna("")
    required = {"gene", "group_id", "n_observations", "cpm", "detection_fraction"}
    if not required.issubset(table):
        raise SystemExit("pseudobulk table lacks required columns")
    current_groups = sorted(group for group in table.group_id.unique() if group.startswith("current::"))
    recall_groups = sorted(group for group in table.group_id.unique() if group.startswith("recall_question::"))
    outside_groups = {
        group[len("outside_current::"):]: group
        for group in table.group_id.unique()
        if group.startswith("outside_current::")
    }
    if not current_groups:
        raise SystemExit("pseudobulk requires current broad groups")
    logcpm = table.pivot(index="gene", columns="group_id", values="cpm").fillna(0)
    logcpm = np.log1p(logcpm)
    current_matrix = logcpm[current_groups]
    eligible_genes = np.array([
        gene for gene in current_matrix.index
        if not TECHNICAL.match(str(gene)) and current_matrix.loc[gene].max() > np.log1p(1)
    ], dtype=object)
    variance = (
        current_matrix.loc[eligible_genes]
        .var(axis=1)
        .fillna(0.0)
        .sort_values(ascending=False)
    )
    informative = variance.head(min(2000, len(variance))).index
    summary_rows: list[dict[str, object]] = []
    similarity_rows: list[dict[str, object]] = []
    differential_rows: list[dict[str, object]] = []
    current_summary_rows: list[dict[str, object]] = []
    current_differential_rows: list[dict[str, object]] = []
    table_index = table.set_index(["group_id", "gene"])
    for group in current_groups:
        target = group[len("current::"):]
        outside_group = outside_groups.get(target, "")
        if not outside_group:
            raise SystemExit(f"current pseudobulk lacks target-versus-outside group: {target}")
        target_marker_genes = sorted(set(markers.loc[
            markers.broad_label.eq(target)
            & markers.evidence_role.eq("positive_family"), "gene"
        ]) & set(logcpm.index))
        target_marker_cpm = table.loc[
            table.group_id.eq(group) & table.gene.isin(target_marker_genes)
        ]
        other_groups = [value for value in current_groups if value != group]
        similarities: list[tuple[str, float]] = []
        for other in other_groups:
            value = spearmanr(
                logcpm.loc[informative, group].to_numpy(),
                logcpm.loc[informative, other].to_numpy(),
            ).correlation
            if not np.isfinite(value):
                value = 0.0
            similarities.append((other[len("current::"):], float(value)))
        similarities.sort(key=lambda item: (-item[1], item[0]))
        current_summary_rows.append({
            "broad_label": target,
            "n_observations": int(table.loc[table.group_id.eq(group), "n_observations"].iloc[0]),
            "nearest_other_current_broad": similarities[0][0] if similarities else "",
            "nearest_other_current_spearman": similarities[0][1] if similarities else 0.0,
            "target_positive_marker_n_available": len(target_marker_genes),
            "target_positive_marker_n_detected": int((target_marker_cpm.cpm > 0).sum()),
            "target_positive_marker_median_log1p_cpm": float(np.log1p(target_marker_cpm.cpm).median()) if len(target_marker_cpm) else 0.0,
            "target_positive_marker_median_detection_fraction": float(target_marker_cpm.detection_fraction.median()) if len(target_marker_cpm) else 0.0,
        })
        genes = logcpm.index
        target_raw = table_index.loc[group].reindex(genes).fillna(0)
        outside_raw = table_index.loc[outside_group].reindex(genes).fillna(0)
        log2fc = np.log2(
            (target_raw.cpm.to_numpy() + 0.1)
            / (outside_raw.cpm.to_numpy() + 0.1)
        )
        order = np.argsort(-np.abs(log2fc))[:100]
        for index in order:
            current_differential_rows.append({
                "broad_label": target,
                "gene": str(genes[index]),
                "current_cpm": float(target_raw.cpm.iloc[index]),
                "current_detection_fraction": float(target_raw.detection_fraction.iloc[index]),
                "outside_cpm": float(outside_raw.cpm.iloc[index]),
                "log2fc_current_vs_outside": float(log2fc[index]),
            })
    for group in recall_groups:
        pieces = group.split("::")
        target = pieces[1]
        origin = pieces[3] if len(pieces) >= 4 and pieces[2] == "from" else ""
        current_target = f"current::{target}"
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
            "nearest_current_broad": similarities[0][0] if similarities else "",
            "nearest_current_spearman": similarities[0][1] if similarities else 0.0,
            "target_current_spearman": dict(similarities).get(target, 0.0),
            "target_positive_marker_n_available": len(target_marker_genes),
            "target_positive_marker_n_detected": int((target_marker_cpm.cpm > 0).sum()),
            "target_positive_marker_median_log1p_cpm": float(np.log1p(target_marker_cpm.cpm).median()) if len(target_marker_cpm) else 0.0,
            "target_positive_marker_median_detection_fraction": float(target_marker_cpm.detection_fraction.median()) if len(target_marker_cpm) else 0.0,
        })
        genes = logcpm.index
        recall_raw = table_index.loc[group].reindex(genes).fillna(0)
        if current_target in logcpm:
            target_raw = table_index.loc[current_target].reindex(genes).fillna(0)
            log2fc_target = np.log2(
                (recall_raw.cpm.to_numpy() + 0.1)
                / (target_raw.cpm.to_numpy() + 0.1)
            )
        else:
            log2fc_target = np.full(len(genes), np.nan)
        if origin and f"current::{origin}" in logcpm:
            origin_raw = table_index.loc[f"current::{origin}"].reindex(genes).fillna(0)
            log2fc_origin = np.log2((recall_raw.cpm.to_numpy() + 0.1) / (origin_raw.cpm.to_numpy() + 0.1))
        else:
            log2fc_origin = np.full(len(genes), np.nan)
        if np.isfinite(log2fc_target).any():
            ranking = np.nan_to_num(np.abs(log2fc_target), nan=-1.0)
        elif np.isfinite(log2fc_origin).any():
            ranking = np.nan_to_num(np.abs(log2fc_origin), nan=-1.0)
        else:
            ranking = recall_raw.cpm.to_numpy()
        order = np.argsort(-ranking)[:100]
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
    current_summary_path = args.out / "broad_cell_type_current_pseudobulk_summary.tsv"
    current_differential_path = args.out / "broad_cell_type_current_vs_outside_top_differential.tsv"
    pd.DataFrame(summary_rows, columns=[
        "recall_group", "target_broad_label", "origin_broad_label",
        "n_observations", "nearest_current_broad", "nearest_current_spearman",
        "target_current_spearman", "target_positive_marker_n_available",
        "target_positive_marker_n_detected",
        "target_positive_marker_median_log1p_cpm",
        "target_positive_marker_median_detection_fraction",
    ]).to_csv(summary_path, sep="\t", index=False)
    pd.DataFrame(similarity_rows, columns=[
        "recall_group", "target_broad_label", "origin_broad_label",
        "current_reference_broad", "spearman_top_variable_genes",
    ]).to_csv(similarity_path, sep="\t", index=False)
    pd.DataFrame(differential_rows, columns=[
        "recall_group", "target_broad_label", "origin_broad_label", "gene",
        "recall_cpm", "recall_detection_fraction", "log2fc_vs_current_target",
        "log2fc_vs_current_origin",
    ]).to_csv(differential_path, sep="\t", index=False)
    pd.DataFrame(current_summary_rows).to_csv(current_summary_path, sep="\t", index=False)
    pd.DataFrame(current_differential_rows).to_csv(current_differential_path, sep="\t", index=False)
    manifest = {
        "schema_version": "2.2", "status": "PASS_EVIDENCE_ONLY",
        "artifact_role": "broad_cell_type_full_transcriptome_pseudobulk_summary",
        "formal_membership_written": False,
        "pseudobulk": artifact(args.pseudobulk),
        "pseudobulk_manifest": artifact(args.pseudobulk_manifest),
        "broad_evidence_manifest": artifact(args.broad_evidence_manifest),
        "marker_manifest": artifact(args.marker_manifest),
        "informative_gene_n": len(informative),
        "current_group_n": len(current_summary_rows),
        "review_group_n": len(summary_rows),
        "artifacts": {
            "summary": artifact(summary_path), "similarity": artifact(similarity_path),
            "top_differential": artifact(differential_path),
            "current_summary": artifact(current_summary_path),
            "current_top_differential": artifact(current_differential_path),
        },
        "warning": "Similarity and differential evidence inform one broad-level biological review; they do not assign cell labels.",
    }
    manifest_path = args.out / "broad_cell_type_review_pseudobulk_summary_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
