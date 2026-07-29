#!/usr/bin/env python3
"""Build one raw-count precision/recall/spatial evidence packet per broad type.

The script never writes annotation membership.  Broad cell type is the review
unit; source groups and spatial components are internal evidence and bounded
candidate regions.  Generic Stromal is evaluated for current-member precision
but is never recalled from the whole query by a marker-only rule.  Oocyte is
routed to the canonical targeted-cohort workflow.
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree

from evidence_schema_lib import sha256
from lineage_controller_lib import (
    GENERIC_REMAINDER_IDS,
    apply_candidate_context,
    candidate_can_release,
    candidate_can_support_broad_review,
    catalog_candidates,
    deterministic_cell_id_set_hash,
    read_tsv,
)


PALETTE = {
    "Epithelial/mesothelial": "#FF2D9A", "Granulosa": "#FFD60A",
    "Glial/Schwann-like": "#64D2FF", "Immune": "#30D158",
    "Luteal": "#FF3B30", "Oocyte": "#FFFFFF",
    "Smooth muscle": "#BF5AF2", "Stromal/mesenchymal": "#58A6FF",
    "Theca": "#FF9F0A", "Endothelial": "#00E5D4",
    "Pericyte/mural": "#B46CFF",
}


def artifact(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "sha256": sha256(path), "n_bytes": path.stat().st_size}


def top_two_mean(values: list[np.ndarray], n: int) -> np.ndarray:
    if not values:
        return np.zeros(n, dtype=np.float32)
    matrix = np.vstack(values)
    if matrix.shape[0] == 1:
        return matrix[0].astype(np.float32)
    partitioned = np.partition(matrix, matrix.shape[0] - 2, axis=0)
    return partitioned[-2:].mean(axis=0).astype(np.float32)


def connected_components(neighbors: np.ndarray, keep: np.ndarray) -> list[np.ndarray]:
    nodes = set(np.flatnonzero(keep).tolist())
    seen: set[int] = set()
    result: list[np.ndarray] = []
    for start in sorted(nodes):
        if start in seen:
            continue
        seen.add(start)
        stack = [start]
        component: list[int] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for other in np.atleast_1d(neighbors[current]):
                other = int(other)
                if other in nodes and other not in seen:
                    seen.add(other)
                    stack.append(other)
        result.append(np.asarray(sorted(component), dtype=np.int64))
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--count-export", required=True, type=Path)
    ap.add_argument("--marker-manifest", required=True, type=Path)
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--review-manifest", required=True, type=Path)
    ap.add_argument("--coordinates", type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--threshold-registry", required=True, type=Path)
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    threshold_document = json.loads(
        args.threshold_registry.read_text(encoding="utf-8")
    )
    scoring_policy = threshold_document["scoring_policy"]
    review_policy = threshold_document["catalog_wide_lineage_review_policy"]
    export_manifest_path = args.count_export / "cell_type_review_count_export_manifest.json"
    export_manifest = json.loads(export_manifest_path.read_text(encoding="utf-8"))
    if (
        export_manifest.get("artifact_role") != "query_raw_count_cell_type_review_export"
        or export_manifest.get("assay_ancestry") != "project_local_non_SCT_raw_counts"
        or export_manifest.get("raw_count_assay") == "SCT"
    ):
        raise SystemExit("cell-type review export is not project-local raw counts")

    review = json.loads(args.review_manifest.read_text(encoding="utf-8"))
    if (
        review.get("stage") != "post_atlas_catalog_wide_lineage_review"
        or review.get("status") != "ITERATION_REQUIRED"
        or review.get("catalog_wide_double_sided_review") is not True
    ):
        raise SystemExit("broad cell-type evidence requires one open canonical review round")
    review_membership = review.get("membership", {})
    if (
        Path(str(review_membership.get("path", ""))).resolve()
        != args.membership.resolve()
        or review_membership.get("sha256") != sha256(args.membership)
    ):
        raise SystemExit("cell-type evidence membership differs from the open review")
    queue_record = review.get("artifacts", {}).get("review_queue", {})
    queue_path = Path(str(queue_record.get("path", "")))
    if not queue_path.is_file() or queue_record.get("sha256") != sha256(queue_path):
        raise SystemExit("cell-type review queue is missing or stale")
    queue = pd.read_csv(queue_path, sep="\t", dtype=str).fillna("")
    if queue.empty or queue.review_id.eq("").any() or queue.review_id.duplicated().any():
        raise SystemExit("cell-type review queue must contain unique nonempty review_id")
    if queue.target_broad_label.eq("").any() or queue.target_broad_label.duplicated().any():
        raise SystemExit("one open review round must contain exactly one task per broad type")

    marker_manifest = pd.read_csv(args.marker_manifest, sep="\t", dtype=str).fillna("")
    gene_map = pd.read_csv(args.count_export / "cell_type_review_gene_map.tsv", sep="\t", dtype=str).fillna("")
    present_genes = gene_map.loc[gene_map.status.eq("matched"), "requested_gene"].tolist()
    cells = pd.read_csv(args.count_export / "cell_type_review_cells.tsv", sep="\t", dtype={"cell_id": str})
    library = pd.read_csv(args.count_export / "cell_type_review_library_size.tsv.gz", sep="\t", dtype={"cell_id": str})
    with gzip.open(args.count_export / "cell_type_review_marker_counts.mtx.gz", "rb") as handle:
        counts = csc_matrix(mmread(handle)).astype(np.float32)
    if counts.shape != (len(present_genes), len(cells)):
        raise SystemExit("raw-count marker matrix dimensions differ from its gene/cell ledgers")
    membership = pd.read_csv(args.membership, sep="\t", dtype=str, low_memory=False).fillna("")
    coordinate_path = args.coordinates or Path(str(export_manifest.get("coordinates", "")))
    if not coordinate_path.is_file():
        raise SystemExit("cell-type review coordinate ledger is missing")
    coordinates = pd.read_csv(coordinate_path, sep="\t", dtype={"cell_id": str})
    required_membership = {"cell_id", "final_broad_label", "source_boundary", "source_cluster"}
    if not required_membership.issubset(membership):
        raise SystemExit("cell-type review membership lacks broad/source columns")
    data = (
        cells[["cell_id"]]
        .merge(library, on="cell_id", validate="one_to_one")
        .merge(coordinates[["cell_id", "x", "y"]], on="cell_id", validate="one_to_one")
        .merge(membership[list(required_membership)], on="cell_id", validate="one_to_one")
    )
    if data.cell_id.tolist() != cells.cell_id.tolist():
        raise SystemExit("cell order changed while joining cell-type review inputs")
    n = len(data)
    xy = data[["x", "y"]].to_numpy(dtype=float)
    k = min(int(review_policy["spatial_knn_k"]), max(1, n - 1))
    tree = cKDTree(xy)
    neighbors = tree.query(xy, k=k + 1, workers=max(1, args.workers))[1][:, 1:]

    gene_to_row = {gene: index for index, gene in enumerate(present_genes)}
    lib = pd.to_numeric(data.total_raw_counts, errors="coerce").fillna(0).to_numpy(np.float32)
    lib[lib <= 0] = 1.0
    scaled: dict[str, np.ndarray] = {}
    detected: dict[str, np.ndarray] = {}
    local_scaled: dict[str, np.ndarray] = {}
    local_detected: dict[str, np.ndarray] = {}
    for gene in present_genes:
        row = counts.getrow(gene_to_row[gene]).tocoo()
        values = np.log1p(1.0e4 * row.data / lib[row.col]).astype(np.float32)
        p95 = float(np.quantile(values, 0.95)) if values.size else 1.0
        if not np.isfinite(p95) or p95 <= 0:
            p95 = 1.0
        score = np.zeros(n, dtype=np.float32)
        score[row.col] = np.minimum(values / p95, 1.5)
        detection = np.zeros(n, dtype=np.float32)
        detection[row.col] = 1.0
        scaled[gene] = score
        detected[gene] = detection
        local_scaled[gene] = score[neighbors].mean(axis=1).astype(np.float32)
        local_detected[gene] = detection[neighbors].mean(axis=1).astype(np.float32)

    family_rows = marker_manifest.loc[
        marker_manifest.evidence_role.eq("positive_family")
        & marker_manifest.gene.isin(present_genes)
    ]
    family_evidence: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    for (candidate_id, family_id), frame in family_rows.groupby(["candidate_id", "family_id"], sort=True):
        genes = sorted(set(frame.gene))
        direct_gene_n = np.vstack([detected[gene] for gene in genes]).sum(axis=0)
        local_gene_n = np.vstack([
            local_detected[gene]
            >= float(scoring_policy["local_gene_detection_fraction"])
            for gene in genes
        ]).sum(axis=0)
        direct_score = top_two_mean([scaled[gene] for gene in genes], n)
        local_score = top_two_mean([local_scaled[gene] for gene in genes], n)
        coherent = (
            (direct_gene_n >= 2)
            | ((direct_gene_n >= 1) & (local_gene_n >= 2))
            | (local_gene_n >= 3)
        )
        family_evidence[(candidate_id, family_id)] = {
            "score": (
                float(scoring_policy["direct_weight"]) * direct_score
                + float(scoring_policy["local_weight"]) * local_score
            ).astype(np.float32),
            "coherent": coherent,
            "direct_gene_n": direct_gene_n.astype(np.int16),
        }

    candidates = catalog_candidates(json.loads(args.catalog.read_text(encoding="utf-8")))
    context = apply_candidate_context(
        candidates, read_tsv(args.context_evidence) if args.context_evidence else []
    )
    candidates = {
        candidate_id: candidate for candidate_id, candidate in candidates.items()
        if candidate_can_support_broad_review(candidate) and candidate_can_release(candidate)
    }
    candidate_evidence: dict[str, dict[str, np.ndarray]] = {}
    candidate_to_broad: dict[str, str] = {}
    for candidate_id, candidate in sorted(candidates.items()):
        families = sorted({
            family for cid, family in family_evidence if cid == candidate_id
        })
        if not families:
            continue
        required = [
            str(family) for family in candidate.get("required_positive_families", [])
            if (candidate_id, str(family)) in family_evidence
        ]
        tested = required or families
        coherent_matrix = np.vstack([
            family_evidence[(candidate_id, family)]["coherent"] for family in tested
        ])
        required_n = len(required) if required else min(2, len(tested))
        family_n = coherent_matrix.sum(axis=0)
        score = top_two_mean([
            family_evidence[(candidate_id, family)]["score"] for family in tested
        ], n)
        direct_family_n = np.vstack([
            family_evidence[(candidate_id, family)]["direct_gene_n"] >= 1
            for family in tested
        ]).sum(axis=0)
        total_direct_gene_n = np.vstack([
            family_evidence[(candidate_id, family)]["direct_gene_n"] for family in tested
        ]).sum(axis=0)
        candidate_evidence[candidate_id] = {
            "score": score,
            "support": (family_n >= required_n) & (
                score > float(review_policy["minimum_raw_review_family_score"])
            ),
            "direct_seed": (
                direct_family_n >= min(
                    int(review_policy["minimum_raw_review_direct_families"]),
                    len(tested),
                )
            ) & (
                total_direct_gene_n
                >= int(review_policy["minimum_raw_review_direct_genes"])
            ),
            "family_n": family_n,
        }
        candidate_to_broad[candidate_id] = str(candidate.get("release_broad_label", ""))

    broad_candidates: dict[str, list[str]] = defaultdict(list)
    for candidate_id, broad in candidate_to_broad.items():
        broad_candidates[broad].append(candidate_id)
    broad_evidence: dict[str, dict[str, np.ndarray]] = {}
    for broad, ids in broad_candidates.items():
        matrix = np.vstack([candidate_evidence[candidate_id]["score"] for candidate_id in ids])
        support = np.vstack([candidate_evidence[candidate_id]["support"] for candidate_id in ids])
        seed = np.vstack([candidate_evidence[candidate_id]["direct_seed"] for candidate_id in ids])
        winner = np.argmax(np.where(support, matrix, -np.inf), axis=0)
        any_support = support.any(axis=0)
        broad_evidence[broad] = {
            "score": np.where(any_support, matrix[winner, np.arange(n)], 0).astype(np.float32),
            "support": any_support,
            "direct_seed": seed.any(axis=0),
        }

    review_targets = queue[[
        "review_id", "review_mode", "target_broad_label", "unit_signature",
    ]].sort_values("target_broad_label", kind="mergesort")
    missing_target_evidence = sorted(
        set(review_targets.target_broad_label) - set(broad_evidence)
    )
    if missing_target_evidence:
        raise SystemExit(
            "open broad review lacks catalog-derived raw-count evidence: "
            + ", ".join(missing_target_evidence)
        )
    target_metadata = {
        str(row.target_broad_label): row
        for row in review_targets.itertuples(index=False)
    }
    target_marker_inventory: dict[str, dict[str, int]] = {}
    for broad, frame in marker_manifest.loc[
        marker_manifest.evidence_role.eq("positive_family")
    ].groupby("broad_label", sort=True):
        requested = set(frame.gene)
        available = requested & set(present_genes)
        requested_families = set(frame.family_id)
        available_families = set(frame.loc[frame.gene.isin(available), "family_id"])
        target_marker_inventory[str(broad)] = {
            "positive_marker_n_requested": len(requested),
            "positive_marker_n_available": len(available),
            "positive_family_n_requested": len(requested_families),
            "positive_family_n_available": len(available_families),
        }
    specific_broads = [broad for broad in broad_evidence if broad != "Stromal/mesenchymal"]
    summary_rows: list[dict[str, object]] = []
    precision_rows: list[dict[str, object]] = []
    component_rows: list[dict[str, object]] = []
    component_member_rows: list[dict[str, object]] = []
    current_question_by_key: dict[tuple[str, str], dict[str, object]] = {}
    plot_rows: list[dict[str, object]] = []
    args.out.mkdir(parents=True, exist_ok=True)
    plot_dir = args.out / "plots"
    plot_dir.mkdir(exist_ok=True)
    for broad in sorted(target_metadata):
        review_row = target_metadata[broad]
        current = data.final_broad_label.eq(broad).to_numpy()
        target = broad_evidence[broad]
        other_broads = [name for name in specific_broads if name != broad]
        if other_broads:
            other_scores = np.vstack([
                np.where(broad_evidence[name]["support"], broad_evidence[name]["score"], 0)
                for name in other_broads
            ])
            best_other_index = np.argmax(other_scores, axis=0)
            best_other_score = other_scores[best_other_index, np.arange(n)]
            best_other_label = np.asarray(other_broads, dtype=object)[best_other_index]
            best_other_label[best_other_score <= 0] = ""
        else:
            best_other_score = np.zeros(n, dtype=np.float32)
            best_other_label = np.full(n, "", dtype=object)
        margin = target["score"] - best_other_score
        inherited = np.zeros(n, dtype=bool)
        for (boundary, cluster), frame in data.loc[current].groupby(
            ["source_boundary", "source_cluster"], sort=True
        ):
            ids = frame.index.to_numpy()
            source_ids = data.index[
                data.source_boundary.eq(boundary) & data.source_cluster.eq(cluster)
            ].to_numpy()
            group_support_fraction = float(target["support"][source_ids].mean())
            current_fraction = len(ids) / max(1, len(source_ids))
            group_supported = (
                group_support_fraction
                >= float(review_policy["minimum_recall_direct_seed_fraction"])
                and float(np.median(target["score"][source_ids]))
                > float(review_policy["minimum_raw_review_group_median_score"])
            )
            use_inheritance = (
                group_supported
                and current_fraction
                >= float(review_policy["minimum_present_label_supported_fraction"])
            )
            if use_inheritance:
                inherited[ids] = ~target["support"][ids]
            competitor = (
                (best_other_label[ids] != "")
                & (
                    best_other_score[ids]
                    >= target["score"][ids]
                    + float(review_policy["minimum_recall_seed_evidence_margin"])
                )
                & ~target["support"][ids]
            )
            precision_rows.append({
                "broad_label": broad, "source_boundary": boundary,
                "source_cluster": cluster, "n_current": len(ids),
                "target_supported_fraction": float(target["support"][ids].mean()),
                "sparse_inherited_fraction": float(inherited[ids].mean()),
                "competitor_question_fraction": float(competitor.mean()),
                "strongest_competitor": Counter(best_other_label[ids][competitor]).most_common(1)[0][0]
                if competitor.any() else "",
                "group_target_program_supported": str(group_supported).lower(),
            })
        effective_current = current & (target["support"] | inherited)
        competitor_question = (
            current & ~effective_current & (best_other_label != "")
            & (
                best_other_score
                >= target["score"]
                + float(review_policy["minimum_recall_seed_evidence_margin"])
            )
        )
        unsupported_current = current & ~effective_current
        for index in np.flatnonzero(unsupported_current):
            current_question_by_key[(broad, str(data.cell_id.iloc[index]))] = {
                "broad_label": broad,
                "cell_id": str(data.cell_id.iloc[index]),
                "question_reason": (
                    "competing_specific_lineage"
                    if competitor_question[index]
                    else "insufficient_current_identity_support"
                ),
                "challenger_broad_label": str(best_other_label[index]),
                "target_score": float(target["score"][index]),
                "challenger_score": float(best_other_score[index]),
            }
        recall = np.zeros(n, dtype=bool)
        accepted_components: list[np.ndarray] = []
        if broad != "Stromal/mesenchymal" and broad != "Oocyte":
            recall = (
                (~current) & target["support"]
                & (
                    margin
                    >= float(review_policy["minimum_raw_review_recall_margin"])
                )
            )
            for component in connected_components(neighbors, recall):
                seed_n = int(target["direct_seed"][component].sum())
                if (
                    len(component)
                    < int(review_policy["minimum_recall_component_members"])
                    or seed_n
                    < int(review_policy["minimum_recall_direct_seed_members"])
                    or seed_n / len(component)
                    < float(review_policy["minimum_recall_direct_seed_fraction"])
                ):
                    continue
                accepted_components.append(component)
                component_id = f"{broad.replace('/', '_').replace(' ', '_')}__{len(accepted_components):05d}"
                component_rows.append({
                    "broad_label": broad, "component_id": component_id,
                    "n_observations": len(component), "n_direct_seeds": seed_n,
                    "direct_seed_fraction": seed_n / len(component),
                    "median_target_score": float(np.median(target["score"][component])),
                    "median_specific_competitor_margin": float(np.median(margin[component])),
                    "current_broad_census": ";".join(
                        f"{key or 'unresolved'}:{value}" for key, value in sorted(
                            Counter(data.final_broad_label.iloc[component]).items()
                        )
                    ),
                })
                for index in component:
                    component_member_rows.append({
                        "broad_label": broad, "component_id": component_id,
                        "cell_id": data.cell_id.iloc[index],
                        "current_broad_label": data.final_broad_label.iloc[index],
                    })
        recall_n = sum(len(component) for component in accepted_components)
        if broad == "Oocyte":
            status = "canonical_oocyte_targeted_review_required"
        elif competitor_question.any() or recall_n:
            status = "cell_type_targeted_review_required"
        else:
            status = "no_raw_count_marker_spatial_challenger"
        summary_rows.append({
            "review_id": review_row.review_id,
            "review_mode": review_row.review_mode,
            "unit_signature": review_row.unit_signature,
            "broad_label": broad,
            "current_n": int(current.sum()),
            "current_marker_supported_n": int((current & target["support"]).sum()),
            "current_sparse_inherited_n": int((current & inherited).sum()),
            "current_effective_support_fraction": float(effective_current.sum() / max(1, current.sum())),
            "current_competitor_question_n": int(competitor_question.sum()),
            "outside_marker_supported_n": int(recall.sum()),
            "accepted_recall_component_n": len(accepted_components),
            "accepted_recall_observation_n": recall_n,
            "review_status": status,
            **target_marker_inventory.get(broad, {
                "positive_marker_n_requested": 0,
                "positive_marker_n_available": 0,
                "positive_family_n_requested": 0,
                "positive_family_n_available": 0,
            }),
        })
        accepted = np.zeros(n, dtype=bool)
        for component in accepted_components:
            accepted[component] = True
        fig, ax = plt.subplots(figsize=(12, 9), facecolor="#050505")
        ax.set_facecolor("#050505")
        ax.scatter(data.x, data.y, s=0.08, c="#707070", alpha=0.25, linewidths=0, rasterized=True)
        ax.scatter(data.loc[current, "x"], data.loc[current, "y"], s=0.34,
                   c=PALETTE.get(broad, "#FF375F"), alpha=0.98, linewidths=0, rasterized=True)
        if accepted.any():
            ax.scatter(data.loc[accepted, "x"], data.loc[accepted, "y"], s=0.54,
                       facecolors="none", edgecolors="#FFFFFF", linewidths=0.22, rasterized=True)
        ax.set_aspect("equal", adjustable="box")
        ax.invert_yaxis(); ax.axis("off")
        ax.set_title(f"{broad} | current {current.sum():,} | recall questions {accepted.sum():,}",
                     color="white", loc="left", fontsize=13)
        fig.tight_layout()
        plot_path = plot_dir / f"{broad.replace('/', '_').replace(' ', '_')}.png"
        fig.savefig(plot_path, dpi=320, facecolor="#050505")
        plt.close(fig)

    # Couple the two sides of the review.  A recall challenger for lineage B
    # is simultaneously an over-recall question for the cell's current label
    # A.  This is especially important for generic Stromal: its ECM program
    # may remain detectable in a misassigned specific lineage, so Stromal
    # precision cannot be closed from its own marker support alone.
    challenger_cells_by_current: dict[str, set[str]] = defaultdict(set)
    recall_cells_by_target: dict[str, set[str]] = defaultdict(set)
    for row in component_member_rows:
        target = str(row["broad_label"])
        current_label = str(row["current_broad_label"])
        cell = str(row["cell_id"])
        recall_cells_by_target[target].add(cell)
        if current_label and current_label != target:
            challenger_cells_by_current[current_label].add(cell)
    summary_by_broad = {str(row["broad_label"]): row for row in summary_rows}
    for broad, row in summary_by_broad.items():
        row["cross_type_over_recall_question_n"] = len(
            challenger_cells_by_current.get(broad, set())
        )
        if row["cross_type_over_recall_question_n"]:
            row["review_status"] = (
                "canonical_oocyte_targeted_review_required"
                if broad == "Oocyte" else "cell_type_targeted_review_required"
            )
        current = data.final_broad_label.eq(broad).to_numpy()
        recall_ids = recall_cells_by_target.get(broad, set())
        over_ids = challenger_cells_by_current.get(broad, set())
        for cell in sorted(over_ids):
            current_question_by_key.setdefault((broad, cell), {
                "broad_label": broad,
                "cell_id": cell,
                "question_reason": "cross_type_recall_challenger",
                "challenger_broad_label": "",
                "target_score": "",
                "challenger_score": "",
            })
        recall_mask = data.cell_id.isin(recall_ids).to_numpy()
        over_mask = data.cell_id.isin(over_ids).to_numpy()
        fig, ax = plt.subplots(figsize=(12, 9), facecolor="#050505")
        ax.set_facecolor("#050505")
        ax.scatter(data.x, data.y, s=0.08, c="#707070", alpha=0.25,
                   linewidths=0, rasterized=True)
        ax.scatter(data.loc[current, "x"], data.loc[current, "y"], s=0.34,
                   c=PALETTE.get(broad, "#FF375F"), alpha=0.98,
                   linewidths=0, rasterized=True)
        if recall_mask.any():
            ax.scatter(data.loc[recall_mask, "x"], data.loc[recall_mask, "y"],
                       s=0.54, facecolors="none", edgecolors="#FFFFFF",
                       linewidths=0.22, rasterized=True)
        if over_mask.any():
            ax.scatter(data.loc[over_mask, "x"], data.loc[over_mask, "y"],
                       s=0.60, facecolors="none", edgecolors="#00E5FF",
                       linewidths=0.24, rasterized=True)
        legend_handles = [
            Line2D([0], [0], marker="o", linestyle="none", markersize=5,
                   markerfacecolor=PALETTE.get(broad, "#FF375F"),
                   markeredgewidth=0, label="current broad membership"),
        ]
        if recall_mask.any():
            legend_handles.append(
                Line2D([0], [0], marker="o", linestyle="none", markersize=5,
                       markerfacecolor="none", markeredgecolor="#FFFFFF",
                       markeredgewidth=0.8,
                       label="outside recall question (not assigned)")
            )
        if over_mask.any():
            legend_handles.append(
                Line2D([0], [0], marker="o", linestyle="none", markersize=5,
                       markerfacecolor="none", markeredgecolor="#00E5FF",
                       markeredgewidth=0.8,
                       label="current member challenged (not removed)")
            )
        ax.legend(
            handles=legend_handles, loc="lower left", frameon=False,
            labelcolor="white", fontsize=8, handletextpad=0.5,
        )
        ax.set_aspect("equal", adjustable="box"); ax.invert_yaxis(); ax.axis("off")
        ax.set_title(
            f"{broad} | current {current.sum():,} | outside recall questions {recall_mask.sum():,} "
            f"| challenged current members {over_mask.sum():,}",
            color="white", loc="left", fontsize=12,
        )
        fig.tight_layout()
        plot_path = plot_dir / f"{broad.replace('/', '_').replace(' ', '_')}.png"
        fig.savefig(plot_path, dpi=320, facecolor="#050505")
        plt.close(fig)
        plot_rows.append({
            "review_id": target_metadata[broad].review_id,
            "broad_label": broad,
            "path": str(plot_path.resolve()),
            "sha256": sha256(plot_path),
        })

    summary_path = args.out / "broad_cell_type_review_summary.tsv"
    precision_path = args.out / "broad_cell_type_current_member_precision.tsv"
    component_path = args.out / "broad_cell_type_outside_recall_components.tsv"
    component_membership_path = args.out / "broad_cell_type_outside_recall_membership.tsv.gz"
    current_question_path = args.out / "broad_cell_type_current_member_questions.tsv.gz"
    plot_index_path = args.out / "broad_cell_type_spatial_plot_index.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    pd.DataFrame(precision_rows).to_csv(precision_path, sep="\t", index=False)
    pd.DataFrame(component_rows).to_csv(component_path, sep="\t", index=False)
    pd.DataFrame(component_member_rows, columns=[
        "broad_label", "component_id", "cell_id", "current_broad_label",
    ]).to_csv(
        component_membership_path, sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame(sorted(
        current_question_by_key.values(),
        key=lambda row: (str(row["broad_label"]), str(row["cell_id"])),
    ), columns=[
        "broad_label", "cell_id", "question_reason",
        "challenger_broad_label", "target_score", "challenger_score",
    ]).to_csv(
        current_question_path, sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame(plot_rows, columns=[
        "review_id", "broad_label", "path", "sha256",
    ]).to_csv(plot_index_path, sep="\t", index=False)
    manifest = {
        "schema_version": "2.2", "status": "PASS_REVIEW_EVIDENCE_ONLY",
        "artifact_role": "broad_cell_type_targeted_review_evidence",
        "primary_review_unit": "one_annotated_broad_cell_type",
        "source_subclusters_are_internal_evidence_only": True,
        "whole_query_recall_scanned_per_type": True,
        "formal_membership_written": False,
        "generic_stromal_whole_query_recall_forbidden": True,
        "oocyte_route": "canonical_query_only_targeted_cohort",
        "count_export_manifest": artifact(export_manifest_path),
        "marker_manifest": artifact(args.marker_manifest),
        "membership": artifact(args.membership),
        "membership_cell_id_semantic_sha256": deterministic_cell_id_set_hash(
            membership.to_dict("records")
        ),
        "review_manifest": artifact(args.review_manifest),
        "review_queue": artifact(queue_path),
        "coordinates": artifact(coordinate_path),
        "catalog": artifact(args.catalog),
        "threshold_registry": artifact(args.threshold_registry),
        "context_release_eligibility": context,
        "reviewed_broad_n": len(summary_rows),
        "review_queue_n": len(queue),
        "artifacts": {
            "summary": artifact(summary_path), "precision": artifact(precision_path),
            "recall_components": artifact(component_path),
            "recall_membership": artifact(component_membership_path),
            "current_member_questions": artifact(current_question_path),
            "spatial_plot_index": artifact(plot_index_path),
        },
        "required_next_step": "for_each_broad_resolve_precision_recall_molecular_spatial_conclusions;full_transcriptome_DEG_pseudobulk_required_for_any_membership_patch",
    }
    manifest_path = args.out / "broad_cell_type_review_evidence_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
