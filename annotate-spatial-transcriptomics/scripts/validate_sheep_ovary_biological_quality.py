#!/usr/bin/env python3
"""Validate sheep-ovary spatial annotation on three biological endpoints."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors


TARGET_CANDIDATES = {
    "granulosa",
    "oocyte",
    "theca_steroidogenic",
    "theca_structural_perifollicular",
    "vascular_endothelial",
    "pericyte_mural",
    "lymphatic_endothelial",
    "smooth_muscle",
    "stromal_mesenchymal",
    "epithelial_mesothelial",
    "immune",
    "neural_schwann_glia",
    "luteal_steroidogenic",
}
FOLLICLE_CANDIDATES = {
    "theca_steroidogenic": "Theca",
    "vascular_endothelial": "Endothelial",
    "pericyte_mural": "Pericyte/mural",
    "smooth_muscle": "Smooth muscle",
    "stromal_mesenchymal": "Stromal/mesenchymal",
}
FOLLICLE_LAYERS = {
    "granulosa_boundary": {
        "candidates": ("granulosa",),
        "target_labels": ("Granulosa",),
        "minimum_angular_sectors": 8,
    },
    "theca_interna": {
        # Theca interna is an identity-bearing steroidogenic layer.  The
        # structural/perifollicular program remains exploratory and cannot be
        # promoted to formal Theca merely because it forms a ring.
        "candidates": ("theca_steroidogenic",),
        "target_labels": ("Theca",),
        "minimum_angular_sectors": 4,
    },
    "endothelial_interna": {
        "candidates": ("vascular_endothelial", "lymphatic_endothelial"),
        "target_labels": ("Endothelial",),
        "minimum_angular_sectors": 2,
    },
    "pericyte_mural_interna": {
        "candidates": ("pericyte_mural",),
        "target_labels": ("Pericyte/mural",),
        "minimum_angular_sectors": 2,
    },
    "outer_nonvascular_contractile": {
        "candidates": ("smooth_muscle",),
        "target_labels": ("Smooth muscle",),
        "minimum_angular_sectors": 4,
    },
    "outer_stromal_background": {
        "candidates": ("stromal_mesenchymal",),
        "target_labels": ("Stromal/mesenchymal",),
        "minimum_angular_sectors": 6,
    },
}
RESTRICTED_BROADS = {
    # Molecularly supported Theca can recur as many small follicular foci.
    # Compactness is therefore reviewed in follicle ROIs and must not act as
    # a whole-section admission/exclusion rule for the lineage.
    "Granulosa", "Oocyte", "Smooth muscle",
    "Epithelial/mesothelial",
}
BOOLEAN_COLUMNS = {
    "family_coherent", "identity_core_coherent", "identity_core_direct",
    "release_family_coherent", "hard_contradiction", "candidate_seed",
    "technical_flag",
}
SCORE_COLUMNS = [
    "cell_id", "source_boundary", "source_cluster", "candidate_id",
    "release_broad_label", "normalized_evidence", "program_score",
    "family_coherent", "identity_core_coherent", "identity_core_direct",
    "release_family_coherent", "hard_contradiction", "candidate_seed",
    "technical_flag", "x", "y",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "n_bytes": path.stat().st_size,
    }


def as_bool(values: pd.Series) -> pd.Series:
    return values.astype(str).str.lower().isin({"true", "1", "yes", "t"})


def read_membership(path: Path) -> tuple[pd.DataFrame, str]:
    frame = pd.read_csv(path, sep="\t", dtype={"cell_id": str})
    if "cell_id" not in frame or frame.cell_id.empty or frame.cell_id.duplicated().any():
        raise ValueError("membership must contain unique nonempty cell_id")
    label = next(
        (name for name in ("final_broad_label", "broad_label") if name in frame),
        "",
    )
    if not label:
        raise ValueError("membership lacks final_broad_label/broad_label")
    frame[label] = frame[label].fillna("").astype(str)
    return frame, label


def validate_canonical_oocyte_review(
    path: Path, membership: pd.DataFrame, label_col: str,
) -> dict[str, object]:
    """Validate an exact, label-blind canonical Oocyte adjudication.

    A canonical targeted cohort can supersede stale ordinary second-round
    Oocyte scores only for the exact released Oocyte member set.  It cannot
    add observations, use spatial location for admission or rely on zona-only
    evidence.
    """
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("status") != "FROZEN_OOCYTE_MEMBERSHIP":
        raise ValueError("canonical Oocyte review is not frozen")
    reviewed_membership = Path(str(document.get("membership_path", "")))
    if (
        not reviewed_membership.is_file()
        or document.get("membership_sha256") != sha256(reviewed_membership)
    ):
        raise ValueError("canonical Oocyte membership is missing or stale")
    reviewed, reviewed_label_col = read_membership(reviewed_membership)
    reviewed_ids = set(
        reviewed.loc[
            reviewed[reviewed_label_col].eq("Oocyte"), "cell_id"
        ].astype(str)
    )
    current_ids = set(
        membership.loc[membership[label_col].eq("Oocyte"), "cell_id"].astype(str)
    )
    if not current_ids or reviewed_ids != current_ids:
        raise ValueError(
            "canonical Oocyte review does not bind the exact released member set"
        )
    final_n = int(document.get("n_final_oocyte_cellbins", -1))
    cluster_n = int(document.get("n_canonical_cluster_cellbins", -1))
    excluded_n = int(
        document.get(
            "n_direct_hard_somatic_contradiction_retained_in_resident_broad", -1,
        )
    )
    if (
        final_n != len(current_ids)
        or cluster_n < final_n
        or excluded_n != cluster_n - final_n
        or document.get("spatial_location_used_for_admission") is not False
        or document.get("zona_only_admission_forbidden") is not True
        or int(document.get("independent_non_zona_deg_gene_n", 0)) < 2
        or float(document.get("cross_resolution_jaccard", 0) or 0) < 0.80
    ):
        raise ValueError("canonical Oocyte review lacks release-grade evidence")
    return document


def score_path_from_manifest(path: Path) -> Path:
    doc = json.loads(path.read_text(encoding="utf-8"))
    record = doc.get("outputs", {}).get("observation_scores")
    if isinstance(record, dict):
        score_path = Path(str(record.get("path", "")))
        expected = str(record.get("sha256", ""))
        if expected and (not score_path.is_file() or sha256(score_path) != expected):
            raise ValueError(f"stale observation scores in {path}")
    else:
        score_path = Path(str(record or ""))
    if not score_path.is_file():
        raise ValueError(f"missing observation scores in {path}")
    return score_path


def legacy_source_boundary(path: Path) -> str:
    return next(
        (part for part in reversed(path.parts) if part.startswith("initial_cluster_")),
        path.parent.name,
    )


def read_scores(paths: list[Path], allow_diagnostic_legacy: bool = False) -> pd.DataFrame:
    selected: list[pd.DataFrame] = []
    for path in paths:
        header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        missing = set(SCORE_COLUMNS).difference(header)
        legacy_columns = {
            "cell_id", "subcluster_id", "candidate_id", "release_broad_label",
            "support_score", "candidate_supported", "identity_core",
            "direct_identity_core", "hard_contradiction", "x", "y",
        }
        legacy = bool(missing and allow_diagnostic_legacy and legacy_columns <= set(header))
        if missing and not legacy:
            raise ValueError(f"{path} lacks score columns: {sorted(missing)}")
        usecols = sorted(legacy_columns) if legacy else SCORE_COLUMNS
        for chunk in pd.read_csv(
            path, sep="\t", usecols=usecols,
            dtype={"cell_id": str, "candidate_id": str}, chunksize=250_000,
        ):
            if legacy:
                chunk = chunk.rename(columns={
                    "subcluster_id": "source_cluster",
                    "support_score": "normalized_evidence",
                    "candidate_supported": "family_coherent",
                    "identity_core": "identity_core_coherent",
                    "direct_identity_core": "identity_core_direct",
                })
                chunk["source_boundary"] = legacy_source_boundary(path)
                chunk["program_score"] = chunk["normalized_evidence"]
                chunk["release_family_coherent"] = chunk["family_coherent"]
                chunk["candidate_seed"] = chunk["identity_core_coherent"]
                chunk["technical_flag"] = False
                chunk = chunk[SCORE_COLUMNS]
            chunk = chunk.loc[chunk.candidate_id.isin(TARGET_CANDIDATES)].copy()
            if not chunk.empty:
                selected.append(chunk)
    if not selected:
        raise ValueError("no sheep-ovary candidate scores were found")
    scores = pd.concat(selected, ignore_index=True)
    if scores.duplicated(["cell_id", "candidate_id"]).any():
        raise ValueError("candidate scores duplicate cell_id x candidate_id")
    for column in BOOLEAN_COLUMNS:
        scores[column] = as_bool(scores[column])
    for column in ("normalized_evidence", "program_score", "x", "y"):
        scores[column] = pd.to_numeric(scores[column], errors="coerce")
    if scores[["x", "y"]].isna().any().any():
        raise ValueError("candidate scores contain invalid coordinates")
    return scores


def coordinate_frame(scores: pd.DataFrame) -> pd.DataFrame:
    coords = scores[["cell_id", "x", "y"]].drop_duplicates()
    if coords.cell_id.duplicated().any():
        spread = coords.groupby("cell_id")[["x", "y"]].nunique().max(axis=1)
        if (spread > 1).any():
            raise ValueError("one observation has inconsistent coordinates")
        coords = coords.drop_duplicates("cell_id")
    return coords.set_index("cell_id", drop=False)


def membership_coordinate_frame(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    if not {"cell_id", "x", "y"}.issubset(header):
        raise ValueError("coordinate membership lacks cell_id/x/y")
    coords = pd.read_csv(
        path, sep="\t", usecols=["cell_id", "x", "y"],
        dtype={"cell_id": str},
    )
    if coords.cell_id.empty or coords.cell_id.duplicated().any():
        raise ValueError("coordinate membership must contain unique cell IDs")
    coords[["x", "y"]] = coords[["x", "y"]].apply(pd.to_numeric, errors="coerce")
    if coords[["x", "y"]].isna().any().any():
        raise ValueError("coordinate membership contains invalid coordinates")
    return coords.set_index("cell_id", drop=False)


def spatial_scale(xy: np.ndarray) -> float:
    if len(xy) < 2:
        return 1.0
    query_index = np.linspace(
        0, len(xy) - 1, min(len(xy), 50_000), dtype=int,
    )
    model = NearestNeighbors(n_neighbors=2, n_jobs=-1).fit(xy)
    distances, _ = model.kneighbors(xy[query_index])
    positive = distances[:, 1][distances[:, 1] > 0]
    if not len(positive):
        return 1.0
    return float(np.median(positive))


def component_summary(xy: np.ndarray, eps: float, min_size: int) -> tuple[np.ndarray, dict[int, int]]:
    if not len(xy):
        return np.array([], dtype=int), {}
    if len(xy) == 1:
        return np.array([0], dtype=int), {0: 1}
    labels = DBSCAN(eps=eps, min_samples=2, n_jobs=-1).fit_predict(xy)
    counts = Counter(int(value) for value in labels if value >= 0)
    keep = {label: n for label, n in counts.items() if n >= min_size}
    return labels, keep


def candidate_view(scores: pd.DataFrame, candidate_id: str) -> pd.DataFrame:
    return scores.loc[scores.candidate_id == candidate_id].set_index("cell_id", drop=False)


def candidate_hits(frame: pd.DataFrame) -> pd.Series:
    evidence = frame.normalized_evidence.fillna(-np.inf)
    return (
        ~frame.hard_contradiction
        & (
            frame.identity_core_coherent
            | (frame.family_coherent & (evidence >= 0.35))
        )
    )


def layer_candidate_view(
    score_views: dict[str, pd.DataFrame],
    candidates: tuple[str, ...],
    member_ids: set[str],
) -> pd.DataFrame:
    frames = [
        score_views[candidate].loc[
            score_views[candidate].index.intersection(member_ids)
        ].copy()
        for candidate in candidates
        if candidate in score_views and not score_views[candidate].empty
    ]
    if not frames:
        return pd.DataFrame(columns=SCORE_COLUMNS + ["hit"])
    combined = pd.concat(frames, ignore_index=True)
    combined["hit"] = candidate_hits(combined)
    combined = combined.sort_values(
        ["cell_id", "hit", "identity_core_direct", "normalized_evidence"],
        ascending=[True, False, False, False],
    ).drop_duplicates("cell_id")
    return combined


def discriminated_direct_layer_ids(
    score_views: dict[str, pd.DataFrame], member_ids: set[str],
    minimum_margin: float = 0.05,
) -> dict[str, set[str]]:
    """Resolve only directly supported identities inside one follicle ROI.

    The four specific wall identities compete before the generic Stromal
    remainder.  Local-only signal is useful for detecting a coherent layer,
    but it is deliberately excluded from the denominator used to demand a
    formal label.  This prevents one mixed cellbin from being required to
    carry several mutually exclusive broad labels.
    """
    layer_views: dict[str, pd.DataFrame] = {}
    for layer_name in (
        "theca_interna", "endothelial_interna", "pericyte_mural_interna",
        "outer_nonvascular_contractile", "outer_stromal_background",
    ):
        spec = FOLLICLE_LAYERS[layer_name]
        frame = layer_candidate_view(
            score_views, tuple(spec["candidates"]), member_ids,
        ).set_index("cell_id", drop=False)
        if len(frame):
            frame["eligible_direct"] = (
                frame.hit & frame.identity_core_direct
                & frame.release_family_coherent
            )
        layer_views[layer_name] = frame
    result = {name: set() for name in layer_views}
    specific = (
        "theca_interna", "endothelial_interna", "pericyte_mural_interna",
        "outer_nonvascular_contractile",
    )
    all_ids = sorted(set().union(*(
        set(frame.index) for frame in layer_views.values() if len(frame)
    )))
    for cell_id in all_ids:
        proposals: list[tuple[float, str]] = []
        for layer_name in specific:
            frame = layer_views[layer_name]
            if cell_id in frame.index and bool(frame.at[cell_id, "eligible_direct"]):
                proposals.append((
                    float(frame.at[cell_id, "normalized_evidence"]), layer_name,
                ))
        proposals.sort(key=lambda value: (-value[0], value[1]))
        if proposals:
            margin = math.inf if len(proposals) == 1 else proposals[0][0] - proposals[1][0]
            if margin >= minimum_margin:
                result[proposals[0][1]].add(cell_id)
            # An unresolved specific-specific tie must not fall through to the
            # generic Stromal remainder.
            continue
        stromal = layer_views["outer_stromal_background"]
        if cell_id in stromal.index and bool(stromal.at[cell_id, "eligible_direct"]):
            result["outer_stromal_background"].add(cell_id)
    return result


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=fields).to_csv(path, sep="\t", index=False)


def broad_spatial_review(
    membership: pd.DataFrame,
    label_col: str,
    scores: pd.DataFrame,
    coords: pd.DataFrame,
    catalog: dict,
    eps: float,
    diagnostic_legacy: bool = False,
) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    candidate_by_broad: dict[str, set[str]] = {}
    for candidate in catalog.get("candidate_boundaries", []):
        broad = str(candidate.get("release_broad_label", ""))
        candidate_id = str(candidate.get("candidate_id", ""))
        if broad and candidate_id:
            candidate_by_broad.setdefault(broad, set()).add(candidate_id)
    rows: list[dict[str, object]] = []
    issues: list[dict[str, str]] = []
    for broad, member_rows in membership.loc[membership[label_col] != ""].groupby(label_col):
        ids = set(member_rows.cell_id)
        candidate_ids = candidate_by_broad.get(str(broad), set())
        evidence = scores.loc[
            scores.cell_id.isin(ids) & scores.candidate_id.isin(candidate_ids)
        ].copy()
        if evidence.empty:
            support_fraction = 0.0
            contradiction_fraction = 0.0
        else:
            evidence["supported"] = candidate_hits(evidence)
            evidence = evidence.sort_values(
                ["cell_id", "supported", "normalized_evidence"],
                ascending=[True, False, False],
            ).drop_duplicates("cell_id")
            support_fraction = float(evidence.supported.mean())
            contradiction_fraction = float(evidence.hard_contradiction.mean())
        present_coords = coords.loc[coords.index.intersection(ids)]
        component_fraction = math.nan
        component_n = 0
        if str(broad) != "Stromal/mesenchymal" and not present_coords.empty:
            labels, keep = component_summary(
                present_coords[["x", "y"]].to_numpy(float),
                eps * 2.5, max(3, min(20, len(present_coords) // 20)),
            )
            component_n = len(keep)
            component_fraction = float(np.isin(labels, list(keep)).mean()) if keep else 0.0
        status = "PASS"
        rationale = "multigene_support_and_nonextreme_spatial_morphology"
        if diagnostic_legacy and str(broad) == "Oocyte":
            status = "NOT_EVALUABLE"
            rationale = "legacy_recovery_scores_do_not_encode_the_canonical_oocyte_route"
        if str(broad) in RESTRICTED_BROADS and len(ids) >= 20:
            if status == "NOT_EVALUABLE":
                pass
            elif support_fraction < 0.05 or (
                contradiction_fraction > 0.50 and support_fraction < 0.20
            ):
                status = "ITERATION_REQUIRED"
                rationale = "restricted_lineage_lacks_member_level_identity_support"
            elif not math.isnan(component_fraction) and component_fraction < 0.30:
                status = "ITERATION_REQUIRED"
                rationale = "restricted_lineage_is_spatially_diffuse_or_fragmented"
        if status == "ITERATION_REQUIRED":
            issues.append({
                "endpoint": "spatial_celltype_localization",
                "scope_id": str(broad),
                "issue_code": rationale,
                "detail": (
                    f"support={support_fraction:.3f}; contradiction="
                    f"{contradiction_fraction:.3f}; component_fraction="
                    f"{component_fraction:.3f}"
                ),
                "recommended_action": "reopen the complete contributing second-round subcluster/cohort",
            })
        rows.append({
            "broad_label": broad,
            "n_observations": len(ids),
            "candidate_ids": ";".join(sorted(candidate_ids)),
            "identity_supported_fraction": support_fraction,
            "hard_contradiction_fraction": contradiction_fraction,
            "spatial_component_n": component_n,
            "component_supported_fraction": component_fraction,
            "status": status,
            "rationale": rationale,
        })
    return rows, issues


def oocyte_review(
    membership: pd.DataFrame,
    label_col: str,
    scores: pd.DataFrame,
    coords: pd.DataFrame,
    eps: float,
    diagnostic_legacy: bool = False,
) -> tuple[list[dict[str, object]], list[dict[str, str]], dict[str, object]]:
    evidence = candidate_view(scores, "oocyte")
    evidence = evidence.copy()
    evidence["supported"] = candidate_hits(evidence)
    labeled_ids = set(membership.loc[membership[label_col] == "Oocyte", "cell_id"])
    group_rows: list[dict[str, object]] = []
    strong_groups = 0
    for (boundary, cluster), group in evidence.groupby(["source_boundary", "source_cluster"]):
        supported_fraction = float(group.supported.mean())
        contradiction_fraction = float(group.hard_contradiction.mean())
        strong = (
            len(group) >= 8
            and supported_fraction >= 0.25
            and contradiction_fraction <= 0.05
        )
        strong_groups += int(strong)
        group_rows.append({
            "source_boundary": boundary,
            "source_cluster": cluster,
            "n_observations": len(group),
            "identity_supported_fraction": supported_fraction,
            "hard_contradiction_fraction": contradiction_fraction,
            "final_oocyte_n": int(group.cell_id.isin(labeled_ids).sum()),
            "canonical_group_supported": strong,
        })
    issues: list[dict[str, str]] = []
    object_n = 0
    labeled_support = math.nan
    labeled_contradiction = math.nan
    status = "PASS"
    rationale = "no_coherent_oocyte_candidate_and_no_oocyte_label"
    if diagnostic_legacy:
        status = "NOT_EVALUABLE"
        rationale = "legacy_recovery_scores_do_not_encode_the_canonical_oocyte_route"
    if labeled_ids:
        labeled = evidence.loc[evidence.index.intersection(labeled_ids)]
        labeled_support = float(labeled.supported.mean()) if len(labeled) else 0.0
        labeled_contradiction = float(labeled.hard_contradiction.mean()) if len(labeled) else 1.0
        xy = coords.loc[coords.index.intersection(labeled_ids), ["x", "y"]].to_numpy(float)
        labels, keep = component_summary(xy, eps * 3.0, 2)
        object_n = len(keep) + int((labels < 0).sum())
        if diagnostic_legacy:
            pass
        elif strong_groups < 1 or labeled_support < 0.20 or labeled_contradiction > 0.05:
            status = "ITERATION_REQUIRED"
            rationale = "released_oocyte_lacks_canonical_group_support"
        else:
            rationale = "canonical_multimodule_group_and_object_morphology_supported"
    elif strong_groups and not diagnostic_legacy:
        status = "ITERATION_REQUIRED"
        rationale = "coherent_oocyte_group_has_zero_final_census"
    if status == "ITERATION_REQUIRED":
        issues.append({
            "endpoint": "oocyte_annotation_quality",
            "scope_id": "Oocyte",
            "issue_code": rationale,
            "detail": (
                f"final_n={len(labeled_ids)}; canonical_group_n={strong_groups}; "
                f"support={labeled_support}; contradiction={labeled_contradiction}"
            ),
            "recommended_action": "reopen the complete canonical Oocyte cohort; do not expand from zona or location",
        })
    summary = {
        "status": status,
        "rationale": rationale,
        "final_oocyte_n": len(labeled_ids),
        "canonical_supported_group_n": strong_groups,
        "spatial_object_n": object_n,
        "identity_supported_fraction": labeled_support,
        "hard_contradiction_fraction": labeled_contradiction,
        "edge_location_used_as_negative_evidence": False,
    }
    return group_rows, issues, summary


def sector_count(xy: np.ndarray, center: np.ndarray, sectors: int = 12) -> int:
    if not len(xy):
        return 0
    angles = np.mod(np.arctan2(xy[:, 1] - center[1], xy[:, 0] - center[0]), 2 * np.pi)
    bins = np.floor(angles / (2 * np.pi / sectors)).astype(int)
    return int(len(np.unique(np.clip(bins, 0, sectors - 1))))


def follicle_review(
    membership: pd.DataFrame,
    label_col: str,
    scores: pd.DataFrame,
    coords: pd.DataFrame,
    eps: float,
    diagnostic_legacy: bool = False,
) -> tuple[
    list[dict[str, object]], list[dict[str, object]],
    list[dict[str, object]], list[dict[str, str]],
    dict[str, object], pd.DataFrame,
]:
    granulosa_ids = set(membership.loc[membership[label_col] == "Granulosa", "cell_id"])
    gran_coords = coords.loc[coords.index.intersection(granulosa_ids)].copy()
    roi_rows: list[dict[str, object]] = []
    recall_rows: list[dict[str, object]] = []
    layer_rows: list[dict[str, object]] = []
    issues: list[dict[str, str]] = []
    if gran_coords.empty:
        return roi_rows, recall_rows, layer_rows, issues, {
            "status": "NOT_EVALUABLE", "rationale": "no_released_granulosa",
            "follicle_roi_n": 0, "antral_roi_n": 0,
            "antral_cavity_status": "NOT_EVALUABLE",
        }, pd.DataFrame(columns=["cell_id", "follicle_roi_id", "distance_to_granulosa", "distance_shell"])

    min_component = max(12, int(math.ceil(len(coords) * 0.00002)))
    labels, keep = component_summary(
        gran_coords[["x", "y"]].to_numpy(float), eps * 2.5, min_component,
    )
    gran_coords["raw_component"] = labels
    gran_coords = gran_coords.loc[gran_coords.raw_component.isin(keep)].copy()
    if gran_coords.empty:
        issue = {
            "endpoint": "follicle_roi_histology",
            "scope_id": "whole_section",
            "issue_code": "granulosa_has_no_coherent_spatial_component",
            "detail": f"released_granulosa_n={len(granulosa_ids)}",
            "recommended_action": "reopen diffuse Granulosa source subclusters before follicle review",
        }
        return roi_rows, recall_rows, layer_rows, [issue], {
            "status": "ITERATION_REQUIRED",
            "rationale": issue["issue_code"], "follicle_roi_n": 0,
            "antral_roi_n": 0, "antral_cavity_status": "NOT_EVALUABLE",
        }, pd.DataFrame(columns=["cell_id", "follicle_roi_id", "distance_to_granulosa", "distance_shell"])

    ordered = sorted(keep, key=lambda value: (-keep[value], value))
    component_map = {value: f"F{index:03d}" for index, value in enumerate(ordered, 1)}
    gran_coords["follicle_roi_id"] = gran_coords.raw_component.map(component_map)
    model = NearestNeighbors(n_neighbors=1, n_jobs=-1).fit(gran_coords[["x", "y"]].to_numpy(float))
    distances, neighbors = model.kneighbors(coords[["x", "y"]].to_numpy(float))
    nearest = gran_coords.iloc[neighbors[:, 0]].follicle_roi_id.to_numpy()
    roi_membership = coords[["cell_id", "x", "y"]].copy()
    roi_membership["follicle_roi_id"] = nearest
    roi_membership["distance_to_granulosa"] = distances[:, 0]
    roi_membership = roi_membership.loc[roi_membership.distance_to_granulosa <= eps * 12].copy()
    roi_membership["distance_shell"] = pd.cut(
        roi_membership.distance_to_granulosa,
        bins=[-np.inf, eps * 2, eps * 5, eps * 8, eps * 12],
        labels=["granulosa_boundary", "inner_wall", "outer_wall", "outer_stroma"],
    ).astype(str)
    roi_membership["radial_distance_from_center"] = np.nan
    roi_membership["granulosa_reference_radius"] = np.nan
    roi_membership["signed_distance_to_granulosa"] = np.nan
    roi_membership["histological_shell"] = ""
    membership_lookup = membership.set_index("cell_id")[label_col]
    fine_lookup = (
        membership.set_index("cell_id")["final_fine_label"].fillna("").astype(str)
        if "final_fine_label" in membership else pd.Series(dtype=str)
    )
    score_views = {
        candidate: candidate_view(scores, candidate)
        for candidate in TARGET_CANDIDATES
    }
    antral_n = 0
    large_unclear_n = 0
    for raw_component in ordered:
        roi_id = component_map[raw_component]
        gran = gran_coords.loc[gran_coords.raw_component == raw_component]
        xy = gran[["x", "y"]].to_numpy(float)
        center = xy.mean(axis=0)
        radial = np.sqrt(((xy - center) ** 2).sum(axis=1))
        q10, q50 = np.quantile(radial, [0.10, 0.50]) if len(radial) > 1 else (0.0, 0.0)
        hole_ratio = float(q10 / q50) if q50 > 0 else 0.0
        sectors = sector_count(xy, center)
        extent = float(np.ptp(xy, axis=0).max()) if len(xy) > 1 else 0.0
        boundary_radius = float(np.median(radial)) if len(radial) else 0.0
        all_xy = coords[["x", "y"]].to_numpy(float)
        all_radius = np.sqrt(((all_xy - center) ** 2).sum(axis=1))
        inner_n = int((all_radius <= boundary_radius * 0.55).sum())
        wall_n = int((
            (all_radius >= boundary_radius * 0.85)
            & (all_radius <= boundary_radius * 1.15)
        ).sum())
        inner_area = math.pi * max((boundary_radius * 0.55) ** 2, eps ** 2)
        wall_area = math.pi * max(
            (boundary_radius * 1.15) ** 2 - (boundary_radius * 0.85) ** 2,
            eps ** 2,
        )
        cavity_density_ratio = float(
            (inner_n / inner_area) / max(wall_n / wall_area, 1e-12)
        )
        large_geometry = bool(
            len(gran) >= max(80, min_component * 3)
            and boundary_radius >= eps * 6
            and extent >= eps * 12
            and sectors >= 8
            and hole_ratio >= 0.25
        )
        antral = bool(large_geometry and cavity_density_ratio <= 0.35)
        gran_fine = fine_lookup.reindex(gran.cell_id).fillna("")
        annotated_antral_anchor_n = int(
            gran_fine.str.contains("large/antral", case=False, regex=False).sum()
        )
        annotated_antral = bool(
            annotated_antral_anchor_n >= max(50, int(math.ceil(len(gran) * 0.10)))
        )
        # A stage-like fine label is supporting context, not anatomical truth.
        # Only a cavity-bearing spatial object receives the mature wall audit.
        antral_expected = antral
        antral_n += int(antral)
        large_unclear = False
        large_unclear_n += int(large_unclear)
        roi = roi_membership.loc[roi_membership.follicle_roi_id == roi_id].copy()
        roi_ids = set(roi.cell_id)
        gran_angles = np.mod(
            np.arctan2(xy[:, 1] - center[1], xy[:, 0] - center[0]),
            2 * np.pi,
        )
        gran_sectors = np.floor(gran_angles / (2 * np.pi / 24)).astype(int)
        sector_reference = {
            int(sector): float(np.median(radial[gran_sectors == sector]))
            for sector in np.unique(gran_sectors)
        }
        roi_xy = roi[["x", "y"]].to_numpy(float)
        roi_radial = np.sqrt(((roi_xy - center) ** 2).sum(axis=1))
        roi_angles = np.mod(
            np.arctan2(roi_xy[:, 1] - center[1], roi_xy[:, 0] - center[0]),
            2 * np.pi,
        )
        roi_sectors = np.floor(roi_angles / (2 * np.pi / 24)).astype(int)
        roi_reference = np.array([
            sector_reference.get(int(sector), boundary_radius)
            for sector in roi_sectors
        ])
        roi["radial_distance_from_center"] = roi_radial
        roi["granulosa_reference_radius"] = roi_reference
        roi["signed_distance_to_granulosa"] = roi_radial - roi_reference
        roi["histological_shell"] = pd.cut(
            roi.signed_distance_to_granulosa,
            bins=[-np.inf, -eps * 2, eps * 2, eps * 5, eps * 8, np.inf],
            labels=[
                "cavity_side", "granulosa_boundary", "follicle_wall_interna",
                "outer_contractile", "outer_stroma",
            ],
        ).astype(str)
        roi_membership.loc[roi.index, [
            "radial_distance_from_center", "granulosa_reference_radius",
            "signed_distance_to_granulosa", "histological_shell",
        ]] = roi[[
            "radial_distance_from_center", "granulosa_reference_radius",
            "signed_distance_to_granulosa", "histological_shell",
        ]]
        discriminated_direct = discriminated_direct_layer_ids(
            score_views, roi_ids,
        )
        roi_status = "PASS"
        roi_issue_codes: list[str] = []
        medians: dict[str, float] = {}
        for candidate_id, target_broad in FOLLICLE_CANDIDATES.items():
            evidence = score_views[candidate_id]
            local = evidence.loc[evidence.index.intersection(roi_ids)].copy()
            if local.empty:
                continue
            local["hit"] = (
                (~local.hard_contradiction & local.identity_core_direct)
                if diagnostic_legacy else candidate_hits(local)
            )
            hit_ids = set(local.loc[local.hit, "cell_id"])
            hit_roi = roi.loc[roi.cell_id.isin(hit_ids)]
            direct_hit_ids = set(
                local.loc[
                    local.hit & local.identity_core_direct
                    & local.release_family_coherent,
                    "cell_id",
                ]
            )
            direct_hit_roi = roi.loc[roi.cell_id.isin(direct_hit_ids)]
            min_hits = max(10, int(math.ceil(len(roi) * 0.005)))
            recurrence = sector_count(hit_roi[["x", "y"]].to_numpy(float), center)
            coherent = len(hit_roi) >= min_hits and recurrence >= 4
            labels_for_hits = membership_lookup.reindex(hit_roi.cell_id).fillna("")
            assigned_fraction = float((labels_for_hits == target_broad).mean()) if len(hit_roi) else 0.0
            generic_fraction = float(
                labels_for_hits.isin({"", "Stromal/mesenchymal"}).mean()
            ) if len(hit_roi) else 0.0
            labels_for_direct = membership_lookup.reindex(
                direct_hit_roi.cell_id
            ).fillna("")
            direct_assigned_fraction = float(
                (labels_for_direct == target_broad).mean()
            ) if len(direct_hit_roi) else math.nan
            direct_generic_fraction = float(
                labels_for_direct.isin({"", "Stromal/mesenchymal"}).mean()
            ) if len(direct_hit_roi) else math.nan
            median_distance = float(hit_roi.distance_to_granulosa.median()) if len(hit_roi) else math.nan
            medians[candidate_id] = median_distance
            under_recalled = bool(
                antral_expected and coherent and len(direct_hit_roi) >= min_hits
                and target_broad != "Stromal/mesenchymal"
                and direct_assigned_fraction < 0.50
                and direct_generic_fraction >= 0.30
            )
            if under_recalled:
                roi_status = "ITERATION_REQUIRED"
                code = f"{candidate_id}_coherent_program_under_recalled"
                roi_issue_codes.append(code)
                issues.append({
                    "endpoint": "follicle_roi_histology",
                    "scope_id": roi_id,
                    "issue_code": code,
                    "detail": (
                        f"hit_n={len(hit_roi)}; direct_hit_n={len(direct_hit_roi)}; "
                        f"sectors={recurrence}; direct_assigned_fraction="
                        f"{direct_assigned_fraction:.3f}; direct_generic_remainder_fraction="
                        f"{direct_generic_fraction:.3f}"
                    ),
                    "recommended_action": "run one raw-count anatomy-conditioned targeted cohort; geometry may trigger review but cannot write labels",
                })
            recall_rows.append({
                "follicle_roi_id": roi_id,
                "candidate_id": candidate_id,
                "target_broad_label": target_broad,
                "program_hit_n": len(hit_roi),
                "direct_identity_hit_n": len(direct_hit_roi),
                "angular_sector_n": recurrence,
                "coherent_multisector_program": coherent,
                "target_assignment_fraction": assigned_fraction,
                "generic_stromal_or_unresolved_fraction": generic_fraction,
                "direct_identity_assignment_fraction": direct_assigned_fraction,
                "direct_identity_generic_remainder_fraction": direct_generic_fraction,
                "median_distance_to_granulosa": median_distance,
                "under_recalled": under_recalled,
            })
        layer_medians: dict[str, float] = {}
        layer_program_medians: dict[str, float] = {}
        layer_released_outer_quantiles: dict[str, float] = {}
        layer_sectors: dict[str, set[int]] = {}
        layer_coherent: dict[str, bool] = {}
        layer_direct_ids: dict[str, set[str]] = {}
        for layer_name, layer_spec in FOLLICLE_LAYERS.items():
            layer = layer_candidate_view(
                score_views, tuple(layer_spec["candidates"]), roi_ids,
            )
            if diagnostic_legacy and len(layer):
                layer["hit"] = ~layer.hard_contradiction & layer.identity_core_direct
            if layer_name == "granulosa_boundary":
                hit_ids = set(gran.cell_id)
            else:
                hit_ids = set(layer.loc[layer.hit, "cell_id"]) if len(layer) else set()
            hit_roi = roi.loc[roi.cell_id.isin(hit_ids)].copy()
            angular_bins = set(
                np.floor(
                    np.mod(
                        np.arctan2(
                            hit_roi.y.to_numpy(float) - center[1],
                            hit_roi.x.to_numpy(float) - center[0],
                        ),
                        2 * np.pi,
                    ) / (2 * np.pi / 12)
                ).astype(int).tolist()
            ) if len(hit_roi) else set()
            minimum_hits = max(
                5 if layer_name in {
                    "endothelial_interna", "pericyte_mural_interna"
                } else 10,
                int(math.ceil(len(roi) * 0.002)),
            )
            if layer_name == "granulosa_boundary":
                minimum_hits = max(20, int(math.ceil(len(gran) * 0.10)))
            coherent = bool(
                len(hit_roi) >= minimum_hits
                and len(angular_bins) >= int(layer_spec["minimum_angular_sectors"])
            )
            target_labels = set(layer_spec["target_labels"])
            assigned = membership_lookup.reindex(hit_roi.cell_id).fillna("")
            assigned_fraction = float(assigned.isin(target_labels).mean()) if len(hit_roi) else 0.0
            generic_fraction = float(
                assigned.isin({"", "Stromal/mesenchymal"}).mean()
            ) if len(hit_roi) else 0.0
            program_median_signed = float(
                hit_roi.signed_distance_to_granulosa.median()
            ) if len(hit_roi) else math.nan
            assigned_hit_roi = hit_roi.loc[assigned.isin(target_labels).to_numpy()].copy()
            released_median_signed = float(
                assigned_hit_roi.signed_distance_to_granulosa.median()
            ) if len(assigned_hit_roi) else math.nan
            direct_ids = set(
                layer.loc[layer.hit & layer.identity_core_direct, "cell_id"]
            ) if len(layer) else set()
            identity_ids = (
                set(gran.cell_id) if layer_name == "granulosa_boundary"
                else discriminated_direct.get(layer_name, set())
            )
            identity_roi = roi.loc[roi.cell_id.isin(identity_ids)].copy()
            identity_labels = membership_lookup.reindex(identity_roi.cell_id).fillna("")
            if layer_name != "granulosa_boundary" and len(identity_roi):
                wall_mutable = identity_labels.isin({
                    "", "Theca", "Endothelial", "Pericyte/mural", "Smooth muscle",
                    "Stromal/mesenchymal",
                })
                identity_roi = identity_roi.loc[wall_mutable.to_numpy()].copy()
                identity_labels = identity_labels.loc[wall_mutable]
            identity_assignment_fraction = float(
                identity_labels.isin(target_labels).mean()
            ) if len(identity_roi) else math.nan
            identity_generic_fraction = float(
                identity_labels.isin({"", "Stromal/mesenchymal"}).mean()
            ) if len(identity_roi) else math.nan
            identity_program_median = float(
                identity_roi.signed_distance_to_granulosa.median()
            ) if len(identity_roi) else math.nan
            released_identity = identity_roi.loc[
                identity_labels.isin(target_labels).to_numpy()
            ].copy()
            released_identity_median = float(
                released_identity.signed_distance_to_granulosa.median()
            ) if len(released_identity) else math.nan
            layer_medians[layer_name] = released_identity_median
            layer_program_medians[layer_name] = identity_program_median
            layer_released_outer_quantiles[layer_name] = float(
                released_identity.signed_distance_to_granulosa.quantile(0.90)
            ) if len(released_identity) else math.nan
            layer_sectors[layer_name] = angular_bins
            layer_coherent[layer_name] = coherent
            layer_direct_ids[layer_name] = direct_ids
            shell_supported = bool(
                np.isfinite(program_median_signed)
                and (
                    abs(program_median_signed) <= eps * 2
                    if layer_name == "granulosa_boundary"
                    else program_median_signed >= -eps * 2
                )
            )
            layer_issue_codes: list[str] = []
            published_target_n = int(
                membership_lookup.reindex(roi.cell_id).fillna("").isin(target_labels).sum()
            )
            # A missing layer is a valid negative/NOT_EVALUABLE audit.  It is
            # actionable only when a coherent direct identity is under-recalled,
            # or when a published label lacks its corresponding molecular core.
            if published_target_n >= minimum_hits and not coherent:
                layer_issue_codes.append(
                    f"{layer_name}_published_label_lacks_corresponding_program"
                )
            if antral_expected and coherent and assigned_fraction >= 0.30 and not shell_supported:
                layer_issue_codes.append(f"{layer_name}_spatial_shell_mismatch")
            if (
                antral_expected and coherent and layer_name != "outer_stromal_background"
                and len(identity_roi) >= minimum_hits
                and identity_assignment_fraction < 0.70
                and identity_generic_fraction >= 0.30
            ):
                layer_issue_codes.append(f"{layer_name}_label_under_recall")
            for code in layer_issue_codes:
                roi_status = "ITERATION_REQUIRED"
                roi_issue_codes.append(code)
                issues.append({
                    "endpoint": "follicle_roi_histology",
                    "scope_id": roi_id,
                    "issue_code": code,
                    "detail": (
                        f"layer={layer_name}; hit_n={len(hit_roi)}; "
                        f"sectors={len(angular_bins)}; program_median_signed="
                        f"{program_median_signed}; released_median_signed="
                        f"{released_median_signed}; "
                        f"direct_discriminated_n={len(identity_roi)}; "
                        f"identity_assignment_fraction={identity_assignment_fraction:.3f}"
                    ),
                    "recommended_action": (
                        "reopen this exact follicle ROI from raw counts and resolve the "
                        "full radial hierarchy; spatial position cannot assign labels"
                    ),
                })
            layer_rows.append({
                "follicle_roi_id": roi_id,
                "layer_name": layer_name,
                "candidate_ids": ";".join(layer_spec["candidates"]),
                "expected_broad_labels": ";".join(layer_spec["target_labels"]),
                "program_hit_n": len(hit_roi),
                "direct_identity_hit_n": len(direct_ids),
                "discriminated_direct_identity_n": len(identity_roi),
                "angular_sector_n": len(angular_bins),
                "coherent_multisector_program": coherent,
                "target_assignment_fraction": assigned_fraction,
                "direct_identity_assignment_fraction": identity_assignment_fraction,
                "median_signed_distance_to_granulosa": program_median_signed,
                "released_median_signed_distance_to_granulosa": released_median_signed,
                "discriminated_direct_median_signed_distance_to_granulosa": identity_program_median,
                "released_discriminated_direct_median_signed_distance_to_granulosa": released_identity_median,
                "expected_shell_supported": shell_supported,
                "status": "ITERATION_REQUIRED" if layer_issue_codes else (
                    "PASS" if coherent else "NOT_EVALUABLE"
                ),
                "issue_codes": ";".join(layer_issue_codes),
            })

        hierarchy_codes: list[str] = []
        if large_unclear:
            hierarchy_codes.append("large_follicle_cavity_not_clear")
        theca_distance = layer_medians.get("theca_interna", math.nan)
        endothelial_distance = layer_medians.get("endothelial_interna", math.nan)
        pericyte_distance = layer_medians.get("pericyte_mural_interna", math.nan)
        theca_program_distance = layer_program_medians.get("theca_interna", math.nan)
        endothelial_program_distance = layer_program_medians.get(
            "endothelial_interna", math.nan,
        )
        pericyte_program_distance = layer_program_medians.get(
            "pericyte_mural_interna", math.nan,
        )
        smooth_distance = layer_medians.get("outer_nonvascular_contractile", math.nan)
        stromal_distance = layer_medians.get("outer_stromal_background", math.nan)
        stromal_outer_extent = layer_released_outer_quantiles.get(
            "outer_stromal_background", math.nan,
        )
        if antral_expected and (
            layer_coherent.get("theca_interna")
            and layer_coherent.get("endothelial_interna")
            and (
                abs(theca_program_distance - endothelial_program_distance) > eps * 4
                or not (
                    layer_sectors["theca_interna"]
                    & layer_sectors["endothelial_interna"]
                )
            )
        ):
            hierarchy_codes.append("theca_endothelial_interna_not_interleaved")
        if antral_expected and (
            layer_coherent.get("endothelial_interna")
            and layer_coherent.get("pericyte_mural_interna")
            and (
                abs(endothelial_program_distance - pericyte_program_distance) > eps * 4
                or not (
                    layer_sectors["endothelial_interna"]
                    & layer_sectors["pericyte_mural_interna"]
                )
            )
        ):
            hierarchy_codes.append("pericyte_mural_not_adjoining_endothelial_branches")
        inner_distances = [
            value for name, value in (
                ("theca_interna", theca_distance),
                ("endothelial_interna", endothelial_distance),
                ("pericyte_mural_interna", pericyte_distance),
            )
            if layer_coherent.get(name) and np.isfinite(value)
        ]
        if antral_expected and (
            inner_distances
            and layer_coherent.get("outer_nonvascular_contractile")
            and smooth_distance <= float(np.median(inner_distances)) + eps * 0.5
        ):
            hierarchy_codes.append("inner_wall_and_outer_contractile_order_inverted")
        if antral_expected and (
            layer_coherent.get("outer_nonvascular_contractile")
            and layer_coherent.get("outer_stromal_background")
            and np.isfinite(stromal_outer_extent)
            and stromal_outer_extent <= smooth_distance + eps * 0.5
        ):
            hierarchy_codes.append("outer_stromal_background_not_resolved_beyond_contractile_layer")
        gran_direct = layer_direct_ids.get("granulosa_boundary", set())
        theca_direct = layer_direct_ids.get("theca_interna", set())
        direct_overlap_fraction = (
            len(gran_direct & theca_direct) / max(1, min(len(gran_direct), len(theca_direct)))
            if gran_direct and theca_direct and not diagnostic_legacy else 0.0
        )
        basement_proxy_status = "NOT_EVALUABLE"
        if (
            antral_expected
            and layer_coherent.get("granulosa_boundary")
            and layer_coherent.get("theca_interna")
        ):
            basement_proxy_status = "PASS"
            if theca_distance < layer_medians["granulosa_boundary"] - eps * 0.5:
                hierarchy_codes.append("granulosa_theca_boundary_order_inverted")
                basement_proxy_status = "ITERATION_REQUIRED"
            if not diagnostic_legacy and direct_overlap_fraction > 0.25:
                hierarchy_codes.append("granulosa_theca_identity_boundary_blurred")
                basement_proxy_status = "ITERATION_REQUIRED"
        for code in hierarchy_codes:
            roi_status = "ITERATION_REQUIRED"
            roi_issue_codes.append(code)
            issues.append({
                "endpoint": "follicle_roi_histology",
                "scope_id": roi_id,
                "issue_code": code,
                "detail": (
                    f"cavity_density_ratio={cavity_density_ratio:.3f}; "
                    f"granulosa={layer_medians.get('granulosa_boundary')}; "
                    f"theca={theca_distance}; endothelial={endothelial_distance}; "
                    f"pericyte_mural={pericyte_distance}; "
                    f"smooth={smooth_distance}; stromal={stromal_distance}; "
                    f"stromal_outer_extent_q90={stromal_outer_extent}; "
                    f"granulosa_theca_direct_overlap={direct_overlap_fraction:.3f}"
                ),
                "recommended_action": (
                    "review the complete follicle wall as one raw-count ROI; repair only "
                    "expression-supported members and preserve unresolved background"
                ),
            })
        roi_rows.append({
            "follicle_roi_id": roi_id,
            "granulosa_component_n": len(gran),
            "roi_observation_n": len(roi),
            "granulosa_extent": extent,
            "granulosa_angular_sector_n": sectors,
            "radial_hole_ratio": hole_ratio,
            "granulosa_boundary_radius": boundary_radius,
            "cavity_inner_observation_n": inner_n,
            "cavity_wall_observation_n": wall_n,
            "cavity_density_ratio": cavity_density_ratio,
            "annotated_antral_anchor_n": annotated_antral_anchor_n,
            "follicle_stage_geometry": (
                "large_antral_candidate" if antral_expected and not large_unclear else
                "large_cavity_unclear" if large_unclear else "small_or_nonantral"
            ),
            "cavity_structure_status": (
                "PASS" if antral_expected and not large_unclear else
                "ITERATION_REQUIRED" if large_unclear else "NOT_APPLICABLE"
            ),
            "basement_membrane_boundary_proxy_status": basement_proxy_status,
            "layer_sequence_status": (
                "ITERATION_REQUIRED" if hierarchy_codes else
                "PASS" if antral_expected else "NOT_EVALUABLE"
            ),
            "status": roi_status,
            "issue_codes": ";".join(roi_issue_codes),
        })
    status = "ITERATION_REQUIRED" if issues else "PASS"
    return roi_rows, recall_rows, layer_rows, issues, {
        "status": status,
        "rationale": (
            "coherent_follicle_program_is_under_recalled_or_layer_order_is_inverted"
            if issues else "detected_follicle_rois_have_concordant_identity_and_structure"
        ),
        "follicle_roi_n": len(roi_rows),
        "antral_roi_n": antral_n,
        "antral_cavity_status": (
            "ITERATION_REQUIRED" if large_unclear_n else
            "PASS" if antral_n else "NOT_EVALUABLE"
        ),
        "large_cavity_unclear_n": large_unclear_n,
        "histological_sequence": (
            "Granulosa boundary -> steroidogenic Theca interleaved with Endothelial "
            "branches -> Pericyte/mural adjoining endothelial branches -> optional "
            "outer nonvascular mature contractile layer -> Stromal background"
        ),
    }, roi_membership


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument(
        "--coordinate-membership", type=Path,
        help=(
            "Optional complete-section membership used only for x/y geometry "
            "during a bounded ROI re-review. It never supplies labels."
        ),
    )
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--scores", action="append", type=Path, default=[])
    parser.add_argument("--scoring-manifest", action="append", type=Path, default=[])
    parser.add_argument(
        "--expected-roi-review", type=Path,
        help=(
            "Optional frozen pre-iteration ROI review. A targeted repair must "
            "not make a previously confirmed large/antral structure disappear."
        ),
    )
    parser.add_argument(
        "--canonical-oocyte-review", type=Path,
        help=(
            "Optional frozen label-blind canonical Oocyte adjudication. It can "
            "supersede stale ordinary second-round Oocyte scores only when its "
            "exact released Oocyte member set matches the reviewed membership."
        ),
    )
    parser.add_argument(
        "--diagnostic-legacy-scores", action="store_true",
        help=(
            "Read the pre-controller recovery score schema for a diagnostic-only "
            "review. The canonical controller never enables this adapter."
        ),
    )
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    score_paths = [path.resolve() for path in args.scores]
    score_paths.extend(score_path_from_manifest(path) for path in args.scoring_manifest)
    if not score_paths:
        raise SystemExit("at least one --scores/--scoring-manifest is required")
    membership, label_col = read_membership(args.membership)
    scores = read_scores(
        score_paths, allow_diagnostic_legacy=args.diagnostic_legacy_scores,
    )
    if not set(scores.cell_id).issubset(set(membership.cell_id)):
        raise SystemExit("scores contain observations outside membership")
    if args.coordinate_membership:
        coords = membership_coordinate_frame(args.coordinate_membership)
    elif {"x", "y"}.issubset(membership.columns):
        coords = membership_coordinate_frame(args.membership)
    else:
        coords = coordinate_frame(scores)
    if not set(scores.cell_id).issubset(set(coords.cell_id)):
        raise SystemExit("score observations are absent from the coordinate frame")
    eps = spatial_scale(coords[["x", "y"]].to_numpy(float))
    catalog = json.loads(args.catalog.read_text(encoding="utf-8"))

    broad_rows, broad_issues = broad_spatial_review(
        membership, label_col, scores, coords, catalog, eps,
        diagnostic_legacy=args.diagnostic_legacy_scores,
    )
    oocyte_rows, oocyte_issues, oocyte_summary = oocyte_review(
        membership, label_col, scores, coords, eps,
        diagnostic_legacy=args.diagnostic_legacy_scores,
    )
    canonical_oocyte = None
    if args.canonical_oocyte_review:
        try:
            canonical_oocyte = validate_canonical_oocyte_review(
                args.canonical_oocyte_review, membership, label_col,
            )
        except ValueError as error:
            raise SystemExit(str(error)) from error
        final_n = int(canonical_oocyte["n_final_oocyte_cellbins"])
        cluster_n = int(canonical_oocyte["n_canonical_cluster_cellbins"])
        excluded_n = int(
            canonical_oocyte[
                "n_direct_hard_somatic_contradiction_retained_in_resident_broad"
            ]
        )
        retention_fraction = final_n / max(cluster_n, 1)
        broad_issues = [
            issue for issue in broad_issues if issue.get("scope_id") != "Oocyte"
        ]
        for row in broad_rows:
            if row.get("broad_label") == "Oocyte":
                row.update({
                    "identity_supported_fraction": retention_fraction,
                    "hard_contradiction_fraction": 0.0,
                    "spatial_component_n": int(
                        canonical_oocyte.get("n_putative_oocyte_objects", 0)
                    ),
                    "component_supported_fraction": 1.0,
                    "status": "PASS",
                    "rationale": (
                        "exact_label_blind_canonical_oocyte_review_supersedes_"
                        "stale_ordinary_second_round_scores"
                    ),
                })
        oocyte_issues = []
        oocyte_rows = [{
            "source_boundary": "canonical_targeted_cohort",
            "source_cluster": str(
                canonical_oocyte.get("selected_resolution", "canonical_cluster")
            ),
            "n_observations": cluster_n,
            "identity_supported_fraction": retention_fraction,
            "hard_contradiction_fraction": excluded_n / max(cluster_n, 1),
            "final_oocyte_n": final_n,
            "canonical_group_supported": True,
        }]
        oocyte_summary = {
            "status": "PASS",
            "rationale": "exact_label_blind_canonical_group_supported",
            "final_oocyte_n": final_n,
            "canonical_supported_group_n": 1,
            "spatial_object_n": int(
                canonical_oocyte.get("n_putative_oocyte_objects", 0)
            ),
            "identity_supported_fraction": retention_fraction,
            "hard_contradiction_fraction": 0.0,
            "edge_location_used_as_negative_evidence": False,
        }
    (
        roi_rows, recall_rows, layer_rows, follicle_issues,
        follicle_summary, roi_membership,
    ) = follicle_review(
        membership, label_col, scores, coords, eps,
        diagnostic_legacy=args.diagnostic_legacy_scores,
    )
    expected_roi_summary = None
    if args.expected_roi_review:
        expected = pd.read_csv(args.expected_roi_review, sep="\t", dtype=str).fillna("")
        expected_antral_n = int((
            (expected.follicle_stage_geometry == "large_antral_candidate")
            & (expected.cavity_structure_status == "PASS")
        ).sum())
        observed_antral_n = int(follicle_summary.get("antral_roi_n", 0))
        expected_roi_summary = {
            "expected_large_antral_n": expected_antral_n,
            "observed_large_antral_n": observed_antral_n,
            "status": "PASS" if observed_antral_n >= expected_antral_n else "ITERATION_REQUIRED",
        }
        if observed_antral_n < expected_antral_n:
            follicle_issues.append({
                "endpoint": "follicle_roi_histology",
                "scope_id": "targeted_repair",
                "issue_code": "previously_detected_large_antral_roi_lost_after_repair",
                "detail": (
                    f"expected_large_antral_n={expected_antral_n}; "
                    f"observed_large_antral_n={observed_antral_n}"
                ),
                "recommended_action": (
                    "reject the repair candidate; freeze the previously passing "
                    "Granulosa/cavity anchor and reopen only typed failing layers"
                ),
            })
            follicle_summary["status"] = "ITERATION_REQUIRED"
    issues = broad_issues + oocyte_issues + follicle_issues
    broad_status = "ITERATION_REQUIRED" if broad_issues else "PASS"
    status = "ITERATION_REQUIRED" if issues else "PASS"
    args.out.mkdir(parents=True, exist_ok=True)
    broad_path = args.out / "broad_spatial_localization_review.tsv"
    oocyte_path = args.out / "oocyte_canonical_group_review.tsv"
    roi_path = args.out / "follicle_roi_histology_review.tsv"
    recall_path = args.out / "follicle_roi_candidate_recall.tsv"
    layer_path = args.out / "follicle_roi_layer_hierarchy.tsv"
    issue_path = args.out / "biological_quality_next_actions.tsv"
    roi_membership_path = args.out / "follicle_roi_membership.tsv.gz"
    write_tsv(broad_path, broad_rows, [
        "broad_label", "n_observations", "candidate_ids",
        "identity_supported_fraction", "hard_contradiction_fraction",
        "spatial_component_n", "component_supported_fraction", "status", "rationale",
    ])
    write_tsv(oocyte_path, oocyte_rows, [
        "source_boundary", "source_cluster", "n_observations",
        "identity_supported_fraction", "hard_contradiction_fraction",
        "final_oocyte_n", "canonical_group_supported",
    ])
    write_tsv(roi_path, roi_rows, [
        "follicle_roi_id", "granulosa_component_n", "roi_observation_n",
        "granulosa_extent", "granulosa_angular_sector_n", "radial_hole_ratio",
        "granulosa_boundary_radius", "cavity_inner_observation_n",
        "cavity_wall_observation_n", "cavity_density_ratio",
        "annotated_antral_anchor_n",
        "follicle_stage_geometry", "cavity_structure_status",
        "basement_membrane_boundary_proxy_status", "layer_sequence_status",
        "status", "issue_codes",
    ])
    write_tsv(recall_path, recall_rows, [
        "follicle_roi_id", "candidate_id", "target_broad_label", "program_hit_n",
        "direct_identity_hit_n",
        "angular_sector_n", "coherent_multisector_program",
        "target_assignment_fraction", "generic_stromal_or_unresolved_fraction",
        "direct_identity_assignment_fraction",
        "direct_identity_generic_remainder_fraction",
        "median_distance_to_granulosa", "under_recalled",
    ])
    write_tsv(layer_path, layer_rows, [
        "follicle_roi_id", "layer_name", "candidate_ids",
        "expected_broad_labels", "program_hit_n", "direct_identity_hit_n",
        "discriminated_direct_identity_n",
        "angular_sector_n", "coherent_multisector_program",
        "target_assignment_fraction", "direct_identity_assignment_fraction",
        "median_signed_distance_to_granulosa",
        "released_median_signed_distance_to_granulosa",
        "discriminated_direct_median_signed_distance_to_granulosa",
        "released_discriminated_direct_median_signed_distance_to_granulosa",
        "expected_shell_supported", "status", "issue_codes",
    ])
    write_tsv(issue_path, issues, [
        "endpoint", "scope_id", "issue_code", "detail", "recommended_action",
    ])
    roi_membership.to_csv(
        roi_membership_path, sep="\t", index=False, compression="gzip",
    )
    manifest = {
        "schema_version": "2.2",
        "status": status,
        "artifact_role": "biological_quality_review",
        "formal_membership_written": False,
        "diagnostic_legacy_score_adapter_used": args.diagnostic_legacy_scores,
        "geometry_can_trigger_review_but_cannot_assign_labels": True,
        "membership": artifact(args.membership),
        "coordinate_membership": (
            artifact(args.coordinate_membership)
            if args.coordinate_membership else None
        ),
        "candidate_catalog": artifact(args.catalog),
        "observation_score_files": [artifact(path) for path in score_paths],
        "expected_roi_review": (
            artifact(args.expected_roi_review) if args.expected_roi_review else None
        ),
        "canonical_oocyte_review": (
            artifact(args.canonical_oocyte_review)
            if args.canonical_oocyte_review else None
        ),
        "spatial_scale": eps,
        "quality_endpoints": {
            "spatial_celltype_localization": {
                "status": broad_status,
                "present_broad_n": len(broad_rows),
                "review": artifact(broad_path),
            },
            "oocyte_annotation_quality": {
                **oocyte_summary,
                "review": artifact(oocyte_path),
            },
            "follicle_roi_histology": {
                **follicle_summary,
                "targeted_repair_structure_preservation": expected_roi_summary,
                "roi_review": artifact(roi_path),
                "candidate_recall": artifact(recall_path),
                "layer_hierarchy": artifact(layer_path),
                "roi_membership": artifact(roi_membership_path),
            },
        },
        "required_next_action_n": len(issues),
        "required_next_actions": artifact(issue_path),
    }
    manifest_path = args.out / "sheep_ovary_biological_quality_review.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
