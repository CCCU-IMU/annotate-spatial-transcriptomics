#!/usr/bin/env python3
"""Build raw-count multi-family challengers for broad labels with zero census.

This is a review trigger, never a label writer.  It deliberately ignores local
signal and spatial proximity: a candidate must be visible through direct genes
from at least two independent positive families.  Fragmented observations are
retained so a vascular branch or another sparse lineage cannot disappear merely
because no five-member spatial component was formed by an earlier scorer.
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csc_matrix

from controller_thresholds import load_controller_thresholds
from evidence_schema_lib import sha256
from lineage_controller_lib import (
    apply_candidate_context,
    candidate_can_release,
    candidate_can_support_broad_review,
    catalog_candidates,
    deterministic_cell_id_set_hash,
    read_tsv,
)


def artifact(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise SystemExit(f"zero-census challenger artifact is missing: {path}")
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "n_bytes": path.stat().st_size,
    }


def load_count_export(root: Path) -> tuple[list[str], list[str], csc_matrix, Path]:
    manifest_path = root / "cell_type_review_count_export_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        manifest.get("artifact_role") != "query_raw_count_cell_type_review_export"
        or manifest.get("assay_ancestry") != "project_local_non_SCT_raw_counts"
        or str(manifest.get("raw_count_assay", "")).upper() == "SCT"
    ):
        raise SystemExit("zero-census challenger requires project-local raw counts")
    gene_map = pd.read_csv(
        root / "cell_type_review_gene_map.tsv", sep="\t", dtype=str
    ).fillna("")
    genes = gene_map.loc[
        gene_map.status.eq("matched"), "requested_gene"
    ].astype(str).tolist()
    cells = pd.read_csv(
        root / "cell_type_review_cells.tsv", sep="\t", dtype=str
    ).fillna("")
    cell_ids = cells.cell_id.astype(str).tolist()
    with gzip.open(root / "cell_type_review_marker_counts.mtx.gz", "rb") as handle:
        counts = csc_matrix(mmread(handle))
    if counts.shape != (len(genes), len(cell_ids)):
        raise SystemExit("zero-census raw-count matrix differs from gene/cell ledgers")
    return genes, cell_ids, counts, manifest_path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--count-export", required=True, type=Path)
    ap.add_argument("--marker-manifest", required=True, type=Path)
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--threshold-registry", required=True, type=Path)
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    policy = load_controller_thresholds(args.threshold_registry)[
        "catalog_wide_lineage_review_policy"
    ]
    min_families = int(policy["minimum_raw_review_direct_families"])
    min_genes = int(policy["minimum_raw_review_direct_genes"])
    min_seed_members = int(policy["minimum_zero_census_direct_seed_observations"])
    min_group_members = int(
        policy["minimum_zero_census_direct_source_group_observations"]
    )

    genes, cell_ids, counts, count_manifest_path = load_count_export(
        args.count_export
    )
    gene_to_row = {gene: index for index, gene in enumerate(genes)}
    membership = pd.read_csv(
        args.membership, sep="\t", dtype=str, low_memory=False
    ).fillna("")
    required_columns = {
        "cell_id", "source_boundary", "source_cluster", "final_broad_label"
    }
    if not required_columns.issubset(membership):
        raise SystemExit("zero-census challenger membership lacks source/broad columns")
    if (
        membership.empty
        or membership.cell_id.eq("").any()
        or membership.cell_id.duplicated().any()
        or set(membership.cell_id.astype(str)) != set(cell_ids)
    ):
        raise SystemExit("zero-census challenger cell universe differs from raw counts")
    membership = membership.set_index("cell_id").loc[cell_ids].reset_index()

    marker = pd.read_csv(
        args.marker_manifest, sep="\t", dtype=str
    ).fillna("")
    marker = marker.loc[
        marker.evidence_role.eq("positive_family")
        & marker.gene.isin(gene_to_row)
    ].copy()
    catalog = catalog_candidates(
        json.loads(args.catalog.read_text(encoding="utf-8"))
    )
    context = apply_candidate_context(
        catalog, read_tsv(args.context_evidence) if args.context_evidence else []
    )
    candidates = {
        candidate_id: candidate
        for candidate_id, candidate in catalog.items()
        if candidate_can_support_broad_review(candidate)
        and candidate_can_release(candidate)
        and str(candidate.get("release_broad_label", ""))
    }
    broad_labels = sorted({
        str(candidate.get("release_broad_label", ""))
        for candidate in candidates.values()
    })
    final_census = Counter(membership.final_broad_label.astype(str))
    atlas_census = Counter(
        membership.atlas_broad.astype(str)
        if "atlas_broad" in membership else []
    )

    candidate_masks: dict[str, dict[str, object]] = {}
    unevaluable: dict[str, list[str]] = defaultdict(list)
    for candidate_id, candidate in sorted(candidates.items()):
        frame = marker.loc[marker.candidate_id.eq(candidate_id)]
        available_families = {
            family: sorted(set(group.gene.astype(str)))
            for family, group in frame.groupby("family_id", sort=True)
        }
        required = [
            str(value) for value in candidate.get("required_positive_families", [])
        ]
        if required and any(family not in available_families for family in required):
            unevaluable[str(candidate["release_broad_label"])].append(candidate_id)
            continue
        tested = required or sorted(available_families)
        if len(tested) < min_families:
            unevaluable[str(candidate["release_broad_label"])].append(candidate_id)
            continue
        family_gene_counts: dict[str, np.ndarray] = {}
        for family in tested:
            rows = [gene_to_row[gene] for gene in available_families[family]]
            family_gene_counts[family] = np.asarray(
                (counts[rows, :] > 0).sum(axis=0)
            ).ravel().astype(np.int16)
        family_n = np.vstack([
            values >= 1 for values in family_gene_counts.values()
        ]).sum(axis=0)
        direct_gene_n = np.vstack(list(family_gene_counts.values())).sum(axis=0)
        seed_families = [
            str(value) for value in candidate.get(
                "seed_required_positive_families", []
            )
            if str(value) in family_gene_counts
        ]
        seed_pool = seed_families or tested
        seed_core = np.vstack([
            family_gene_counts[family] >= 2 for family in seed_pool
        ]).any(axis=0)
        supported = (
            (family_n >= min_families)
            & (direct_gene_n >= min_genes)
            & seed_core
        )
        candidate_masks[candidate_id] = {
            "broad": str(candidate["release_broad_label"]),
            "supported": supported,
            "family_n": family_n,
            "direct_gene_n": direct_gene_n,
        }

    member_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for broad in broad_labels:
        ids = [
            candidate_id for candidate_id, evidence in candidate_masks.items()
            if evidence["broad"] == broad
        ]
        if ids:
            direct_gene = np.vstack([
                candidate_masks[candidate_id]["direct_gene_n"] for candidate_id in ids
            ])
            family_n = np.vstack([
                candidate_masks[candidate_id]["family_n"] for candidate_id in ids
            ])
            support = np.vstack([
                candidate_masks[candidate_id]["supported"] for candidate_id in ids
            ])
            ranking = np.where(support, direct_gene * 100 + family_n, -1)
            winner = ranking.argmax(axis=0)
            any_support = support.any(axis=0)
        else:
            winner = np.zeros(len(cell_ids), dtype=int)
            any_support = np.zeros(len(cell_ids), dtype=bool)
        direct_indices = np.flatnonzero(any_support)
        group_counts: Counter[tuple[str, str]] = Counter()
        for index in direct_indices:
            candidate_id = ids[int(winner[index])]
            row = membership.iloc[int(index)]
            group = (str(row.source_boundary), str(row.source_cluster))
            group_counts[group] += 1
            member_rows.append({
                "broad_label": broad,
                "candidate_id": candidate_id,
                "cell_id": str(row.cell_id),
                "source_boundary": group[0],
                "source_cluster": group[1],
                "current_broad_label": str(row.final_broad_label),
                "direct_family_n": int(family_n[int(winner[index]), index]),
                "direct_gene_n": int(direct_gene[int(winner[index]), index]),
            })
        direct_n = len(direct_indices)
        max_group_n = max(group_counts.values(), default=0)
        final_n = int(final_census[broad])
        if final_n:
            status = "present_broad_not_zero_census"
            rationale = "current_broad_membership_exists"
        elif not ids:
            status = "not_evaluable"
            rationale = "required_direct_marker_families_unavailable"
        elif direct_n >= min_seed_members and max_group_n >= min_group_members:
            status = "review_required"
            rationale = "direct_multifamily_program_requires_bounded_source_review"
        else:
            status = "refuted_no_material_direct_multifamily_program"
            rationale = "direct_multifamily_support_below_review_trigger"
        summary_rows.append({
            "broad_label": broad,
            "candidate_ids": ";".join(ids),
            "context_status": ";".join(sorted({
                str(context.get(candidate_id, {}).get("status", "evaluable"))
                for candidate_id in ids
            })) if ids else "not_evaluable",
            "final_n_observations": final_n,
            "direct_multifamily_n": direct_n,
            "direct_source_group_n": len(group_counts),
            "max_source_group_direct_n": max_group_n,
            "atlas_prediction_n": int(atlas_census[broad]),
            "status": status,
            "rationale": rationale,
        })

    args.out.mkdir(parents=True, exist_ok=True)
    summary_path = args.out / "zero_census_direct_multifamily_challengers.tsv"
    members_path = args.out / "zero_census_direct_multifamily_membership.tsv.gz"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    pd.DataFrame(member_rows, columns=[
        "broad_label", "candidate_id", "cell_id", "source_boundary",
        "source_cluster", "current_broad_label", "direct_family_n",
        "direct_gene_n",
    ]).to_csv(members_path, sep="\t", index=False, compression="gzip")
    manifest = {
        "schema_version": "2.5",
        "status": "PASS",
        "stage": "zero_census_direct_multifamily_challenger_audit",
        "artifact_role": "query_derived_zero_census_challenger",
        "formal_membership_written": False,
        "direct_signal_only": True,
        "spatial_geometry_used_for_identity": False,
        "membership": artifact(args.membership),
        "membership_cell_id_set_sha256": deterministic_cell_id_set_hash(
            membership.to_dict("records")
        ),
        "count_export_manifest": artifact(count_manifest_path),
        "marker_manifest": artifact(args.marker_manifest),
        "candidate_catalog": artifact(args.catalog),
        "threshold_registry": artifact(args.threshold_registry),
        "context_evidence": artifact(args.context_evidence) if args.context_evidence else None,
        "thresholds": {
            "minimum_direct_families": min_families,
            "minimum_direct_genes": min_genes,
            "minimum_direct_seed_observations": min_seed_members,
            "minimum_source_group_observations": min_group_members,
        },
        "artifacts": {
            "summary": artifact(summary_path),
            "membership": artifact(members_path),
        },
        "review_required_broad_n": sum(
            row["status"] == "review_required" for row in summary_rows
        ),
    }
    manifest_path = args.out / "zero_census_direct_challenger_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
