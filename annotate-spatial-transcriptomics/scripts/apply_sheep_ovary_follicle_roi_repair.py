#!/usr/bin/env python3
"""Apply one authority-bound sheep-ovary follicle ROI identity repair.

Geometry selects the bounded review cohort. Formal writeback is driven only by
direct, coherent query-derived identity evidence from a raw-count ROI recluster.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path

import pandas as pd

from evidence_schema_lib import sha256
from lineage_controller_lib import (
    apply_candidate_context, candidate_can_release,
    candidate_can_support_broad_review, catalog_candidates,
    deterministic_membership_hash, read_tsv as read_record_tsv,
)


CANDIDATE_TO_BROAD = {
    "theca_steroidogenic": "Theca",
    "vascular_endothelial": "Endothelial",
    "pericyte_mural": "Pericyte/mural",
    "lymphatic_endothelial": "Endothelial",
    "smooth_muscle": "Smooth muscle",
    "stromal_mesenchymal": "Stromal/mesenchymal",
}
SPECIFIC_BROADS = ("Theca", "Endothelial", "Pericyte/mural", "Smooth muscle")
WALL_BROADS = set(SPECIFIC_BROADS) | {"Stromal/mesenchymal"}
BOOL_COLUMNS = (
    "family_coherent", "identity_core_direct", "release_family_coherent",
    "hard_contradiction", "technical_flag",
)
REQUIRED_SCORE_COLUMNS = {
    "cell_id", "candidate_id", "normalized_evidence", "family_coherent",
    "identity_core_direct", "release_family_coherent", "hard_contradiction",
    "technical_flag", "positive_families", "x", "y",
}


def artifact(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "sha256": sha256(path), "n_bytes": path.stat().st_size}


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"cell_id": str}, **kwargs)


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes", "t"})


def parse_bound(values: list[str], name: str) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise SystemExit(f"{name} must use ROI_ID=/absolute/path")
        roi, raw_path = value.split("=", 1)
        path = Path(raw_path).resolve()
        if not roi or roi in result or not path.is_file():
            raise SystemExit(f"invalid or duplicate {name}: {value}")
        result[roi] = path
    return result


def authorized_broads(issue_rows: pd.DataFrame, layer_rows: pd.DataFrame, roi: str) -> set[str]:
    issues = issue_rows.loc[issue_rows.scope_id.astype(str) == roi]
    codes = set(issues.issue_code.astype(str))
    hierarchy_tokens = (
        "order_inverted", "not_interleaved", "background_not_resolved",
        "identity_boundary_blurred", "complete_follicle_wall",
    )
    if any(any(token in code for token in hierarchy_tokens) for code in codes):
        return set(WALL_BROADS)
    allowed: set[str] = set()
    mapping = {
        "theca": "Theca", "endothelial": "Endothelial",
        "vascular": "Endothelial", "pericyte": "Pericyte/mural",
        "mural": "Pericyte/mural", "lymphatic": "Endothelial",
        "smooth": "Smooth muscle", "contractile": "Smooth muscle",
        "stromal": "Stromal/mesenchymal",
    }
    for code in codes:
        for token, broad in mapping.items():
            if token in code:
                allowed.add(broad)
    if not layer_rows.empty:
        failing = layer_rows.loc[
            (layer_rows.follicle_roi_id.astype(str) == roi)
            & (layer_rows.status.astype(str) == "ITERATION_REQUIRED")
        ]
        layer_to_broad = {
            "theca_interna": "Theca",
            "endothelial_interna": "Endothelial",
            "pericyte_mural_interna": "Pericyte/mural",
            "outer_nonvascular_contractile": "Smooth muscle",
            "outer_stromal_background": "Stromal/mesenchymal",
        }
        allowed.update(layer_to_broad[name] for name in failing.layer_name if name in layer_to_broad)
    return allowed


def read_score(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    missing = REQUIRED_SCORE_COLUMNS.difference(header)
    if missing:
        raise SystemExit(f"repair score lacks canonical columns {sorted(missing)}: {path}")
    frame = read_tsv(path)
    for column in BOOL_COLUMNS:
        frame[column] = as_bool(frame[column])
    frame["normalized_evidence"] = pd.to_numeric(frame.normalized_evidence, errors="coerce")
    frame[["x", "y"]] = frame[["x", "y"]].apply(
        pd.to_numeric, errors="coerce"
    )
    if frame[["normalized_evidence", "x", "y"]].isna().any().any():
        raise SystemExit(f"repair score has invalid evidence: {path}")
    if frame.duplicated(["cell_id", "candidate_id"]).any():
        raise SystemExit(f"repair score duplicates cell_id x candidate_id: {path}")
    return frame


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", required=True, type=Path)
    ap.add_argument("--stage-authority", required=True, type=Path)
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--quality-review", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument("--base-scores", action="append", required=True, type=Path)
    ap.add_argument("--repair-score", action="append", default=[])
    ap.add_argument("--repair-ancestry", action="append", default=[])
    ap.add_argument("--specific-margin", type=float, default=0.05)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    if args.specific_margin < 0:
        raise SystemExit("specific margin must be nonnegative")

    authority = json.loads(args.stage_authority.read_text(encoding="utf-8"))
    if (
        authority.get("mode") != "stage_authority"
        or authority.get("phase") != "atlas_and_completeness_review"
        or authority.get("annotation_contract_sha256") != sha256(args.contract)
    ):
        raise SystemExit("stage authority does not permit follicle ROI repair")
    bound_records = [
        ("pre_repair_membership", args.membership),
        ("pre_repair_biological_quality", args.quality_review),
        ("candidate_catalog", args.catalog),
    ]
    if args.context_evidence:
        bound_records.append(("context_evidence", args.context_evidence))
    for key, path in bound_records:
        record = authority.get(key, {})
        if Path(str(record.get("path", ""))).resolve() != path.resolve() or record.get("sha256") != sha256(path):
            raise SystemExit(f"{key} differs from stage authority")

    candidates = catalog_candidates(json.loads(args.catalog.read_text(encoding="utf-8")))
    context_summary = apply_candidate_context(
        candidates,
        read_record_tsv(args.context_evidence) if args.context_evidence else [],
    )
    eligible_candidate_to_broad = {
        candidate_id: broad
        for candidate_id, broad in CANDIDATE_TO_BROAD.items()
        if candidate_id in candidates
        and candidate_can_release(candidates[candidate_id])
        and candidate_can_support_broad_review(candidates[candidate_id])
        and str(candidates[candidate_id].get("release_broad_label", "")) == broad
    }
    if not eligible_candidate_to_broad:
        raise SystemExit("follicle ROI repair has no context-eligible wall candidate")

    review = json.loads(args.quality_review.read_text(encoding="utf-8"))
    if review.get("status") != "ITERATION_REQUIRED":
        raise SystemExit("follicle repair requires an ITERATION_REQUIRED biological review")
    endpoints = review.get("quality_endpoints", {}).get("follicle_roi_histology", {})
    review_paths = {
        key: Path(str(endpoints.get(key, {}).get("path", ""))).resolve()
        for key in ("roi_membership", "layer_hierarchy", "roi_review")
    }
    issue_path = Path(str(review.get("required_next_actions", {}).get("path", ""))).resolve()
    for path in [*review_paths.values(), issue_path]:
        if not path.is_file():
            raise SystemExit("biological review is missing a follicle artifact")
    for key, path in review_paths.items():
        record = endpoints.get(key, {})
        if record.get("sha256") and record.get("sha256") != sha256(path):
            raise SystemExit(f"biological review {key} is stale")
    issues = read_tsv(issue_path).fillna("")
    issues = issues.loc[(issues.endpoint == "follicle_roi_histology") & issues.scope_id.astype(str).str.match(r"^F\d+$")]
    all_issue_rois = sorted(set(issues.scope_id.astype(str)))
    if not all_issue_rois:
        raise SystemExit("no bounded follicle ROI is authorized for repair")
    roi_membership = read_tsv(review_paths["roi_membership"]).fillna("")
    layer_rows = read_tsv(review_paths["layer_hierarchy"]).fillna("")
    roi_review = read_tsv(review_paths["roi_review"]).fillna("")
    repair_rois = [
        roi for roi in all_issue_rois
        if authorized_broads(issues, layer_rows, roi)
    ]
    skipped_non_wall_rois = sorted(set(all_issue_rois) - set(repair_rois))
    if not repair_rois:
        raise SystemExit("no bounded follicle ROI has typed wall-layer repair authority")
    repair_scores = parse_bound(args.repair_score, "repair-score")
    repair_ancestry = parse_bound(args.repair_ancestry, "repair-ancestry")
    if set(repair_scores) != set(repair_rois) or set(repair_ancestry) != set(repair_rois):
        raise SystemExit("every authorized ROI requires exactly one repair score and raw-count ancestry")
    supplied_base = {(path.resolve(), sha256(path)) for path in args.base_scores}
    bound_base = {
        (Path(str(record.get("path", ""))).resolve(), str(record.get("sha256", "")))
        for record in authority.get("base_scores", [])
    }
    supplied_repair = {(path.resolve(), sha256(path)) for path in repair_scores.values()}
    bound_repair = {
        (Path(str(record.get("path", ""))).resolve(), str(record.get("sha256", "")))
        for record in authority.get("repair_scores", [])
    }
    supplied_ancestry = {(path.resolve(), sha256(path)) for path in repair_ancestry.values()}
    bound_ancestry = {
        (Path(str(record.get("path", ""))).resolve(), str(record.get("sha256", "")))
        for record in authority.get("repair_ancestry", [])
    }
    if supplied_base != bound_base or supplied_repair != bound_repair or supplied_ancestry != bound_ancestry:
        raise SystemExit("repair evidence differs from stage authority")

    membership = read_tsv(args.membership).fillna("")
    if membership.cell_id.empty or membership.cell_id.duplicated().any() or "final_broad_label" not in membership:
        raise SystemExit("membership must contain unique cells and final_broad_label")
    original = membership.copy(deep=True).set_index("cell_id", drop=False)
    revised = membership.copy(deep=True).set_index("cell_id", drop=False)
    all_repair_frames: list[pd.DataFrame] = []
    changes: list[dict[str, object]] = []
    typed_authority: dict[str, list[str]] = {}

    for roi in repair_rois:
        ancestry = json.loads(repair_ancestry[roi].read_text(encoding="utf-8"))
        raw_assay = str(ancestry.get("raw_count_assay", ""))
        if ancestry.get("status") != "PASS" or not raw_assay or raw_assay == "SCT" or "raw_counts_SCT" not in str(ancestry.get("clustering_path", "")):
            raise SystemExit(f"{roi} repair does not bind a non-SCT raw-count recluster")
        roi_ids = set(roi_membership.loc[roi_membership.follicle_roi_id.astype(str) == roi, "cell_id"].astype(str))
        if not roi_ids or not roi_ids <= set(revised.index):
            raise SystemExit(f"{roi} membership is empty or outside the full membership")
        score = read_score(repair_scores[roi])
        if set(score.cell_id) != roi_ids:
            raise SystemExit(f"{roi} repair score must exactly cover its bounded ROI")
        all_repair_frames.append(score)
        allowed = authorized_broads(issues, layer_rows, roi)
        if not allowed:
            raise SystemExit(f"{roi} has no typed failing wall layer")
        typed_authority[roi] = sorted(allowed)
        wall_ids = sorted(cell for cell in roi_ids if str(revised.at[cell, "final_broad_label"]) in WALL_BROADS | {""})
        scored = score.loc[
            score.cell_id.isin(wall_ids)
            & score.candidate_id.isin(eligible_candidate_to_broad)
        ].copy()
        scored["target_broad"] = scored.candidate_id.map(eligible_candidate_to_broad)
        scored = scored.loc[scored.target_broad.isin(allowed)]
        scored["eligible"] = (
            scored.family_coherent & scored.identity_core_direct
            & scored.release_family_coherent & ~scored.hard_contradiction
            & ~scored.technical_flag
        )
        scored = scored.loc[scored.eligible].sort_values(
            ["cell_id", "target_broad", "normalized_evidence", "candidate_id"],
            ascending=[True, True, False, True],
        ).drop_duplicates(["cell_id", "target_broad"])
        by_cell = {cell: group for cell, group in scored.groupby("cell_id", sort=True)}
        roi_shell = roi_membership.set_index("cell_id").reindex(wall_ids)
        for cell in wall_ids:
            group = by_cell.get(cell)
            if group is None or group.empty:
                continue
            proposals = []
            for broad in SPECIFIC_BROADS:
                row = group.loc[group.target_broad == broad]
                if not row.empty:
                    best = row.iloc[0]
                    proposals.append((float(best.normalized_evidence), broad, str(best.candidate_id), str(best.positive_families)))
            proposals.sort(key=lambda value: (-value[0], value[1], value[2]))
            winner = None
            if proposals:
                margin = math.inf if len(proposals) == 1 else proposals[0][0] - proposals[1][0]
                if margin >= args.specific_margin:
                    winner = (*proposals[0], margin, proposals[1][1] if len(proposals) > 1 else "")
            if winner is None and not proposals and "Stromal/mesenchymal" in allowed:
                row = group.loc[group.target_broad == "Stromal/mesenchymal"]
                if not row.empty:
                    best = row.iloc[0]
                    winner = (float(best.normalized_evidence), "Stromal/mesenchymal", str(best.candidate_id), str(best.positive_families), math.inf, "")
            if winner is None:
                continue
            evidence_value, new_broad, candidate_id, positive_families, margin, competitor = winner
            old_broad = str(revised.at[cell, "final_broad_label"])
            if old_broad == new_broad:
                continue
            revised.at[cell, "pre_follicle_repair_broad_label"] = old_broad
            revised.at[cell, "final_broad_label"] = new_broad
            # The frozen membership must carry the identity that justified the
            # repaired broad label.  Keeping the pre-repair candidate here makes
            # the downstream completeness audit report a false missing-source
            # failure even though the bounded ROI evidence is valid.
            revised.at[cell, "candidate_id"] = candidate_id
            if "confidence" in revised:
                revised.at[cell, "confidence"] = "high"
            if "final_broad_confidence" in revised:
                revised.at[cell, "final_broad_confidence"] = "high"
            revised.at[cell, "final_state"] = "defined_broad_only"
            revised.at[cell, "qc_reason"] = ""
            revised.at[cell, "assignment_origin"] = "follicle_roi_raw_count_direct_identity_repair"
            revised.at[cell, "broad_freeze_source"] = "post_atlas_follicle_roi_direct_identity"
            revised.at[cell, "follicle_roi_review_id"] = roi
            revised.at[cell, "follicle_roi_repair_status"] = "reassigned"
            revised.at[cell, "follicle_roi_repair_candidate_id"] = candidate_id
            revised.at[cell, "follicle_roi_repair_source"] = "raw_counts_SCT_PCA_SNN_Leiden_full_catalog"
            if "final_fine_label" in revised and str(revised.at[cell, "final_fine_label"]) and old_broad != new_broad:
                revised.at[cell, "final_fine_label"] = ""
                for column in ("final_fine_confidence", "final_fine_candidate_id", "final_fine_assignment_source"):
                    if column in revised:
                        revised.at[cell, column] = ""
            changes.append({
                "cell_id": cell, "follicle_roi_id": roi,
                "histological_shell": str(roi_shell.at[cell, "histological_shell"]) if "histological_shell" in roi_shell else "",
                "old_broad_label": old_broad, "new_broad_label": new_broad,
                "candidate_id": candidate_id, "normalized_evidence": evidence_value,
                "positive_families": positive_families,
                "strongest_competitor": competitor,
                "evidence_margin": margin if math.isfinite(margin) else "",
                "assignment_origin": "follicle_roi_raw_count_direct_identity_discriminator",
            })

    revised = revised.loc[original.index]
    outside = set(original.index) - set(roi_membership.loc[roi_membership.follicle_roi_id.isin(repair_rois), "cell_id"])
    if not original.loc[sorted(outside), "final_broad_label"].equals(revised.loc[sorted(outside), "final_broad_label"]):
        raise SystemExit("repair changed broad labels outside bounded ROIs")
    protected = original.final_broad_label.isin({
        "Granulosa", "Oocyte", "Immune", "Epithelial/mesothelial",
        "Luteal", "Luteal/steroidogenic", "Glial/Schwann-like",
    })
    if not original.loc[protected, "final_broad_label"].equals(revised.loc[protected, "final_broad_label"]):
        raise SystemExit("repair changed a protected non-wall identity")

    base_frames = [read_score(path) for path in args.base_scores]
    base = pd.concat(base_frames, ignore_index=True)
    if base.duplicated(["cell_id", "candidate_id"]).any():
        raise SystemExit("base scores duplicate cell_id x candidate_id")
    repaired_ids = set().union(*(set(frame.cell_id) for frame in all_repair_frames))
    combined = pd.concat(
        [base.loc[~base.cell_id.isin(repaired_ids)], *all_repair_frames],
        ignore_index=True,
    ).sort_values(["cell_id", "candidate_id"], kind="mergesort")
    if combined.duplicated(["cell_id", "candidate_id"]).any():
        raise SystemExit("combined post-repair scores duplicate cell_id x candidate_id")
    if set(combined.cell_id) != set(revised.index):
        raise SystemExit("combined post-repair scores do not cover full membership")
    coordinate_rows = combined[["cell_id", "x", "y"]].drop_duplicates()
    if coordinate_rows.cell_id.duplicated().any():
        raise SystemExit("combined post-repair scores contain inconsistent coordinates")
    if set(coordinate_rows.cell_id) != set(revised.index):
        raise SystemExit("post-repair coordinate ledger does not cover full membership")

    args.out.mkdir(parents=True, exist_ok=True)
    membership_path = args.out / "post_follicle_roi_repair_membership.tsv.gz"
    score_path = args.out / "post_follicle_roi_repair_observation_scores.tsv.gz"
    coordinate_path = args.out / "post_follicle_roi_repair_coordinates.tsv.gz"
    change_path = args.out / "follicle_roi_repair_changes.tsv"
    revised.reset_index(drop=True).to_csv(membership_path, sep="\t", index=False, compression="gzip")
    combined.to_csv(score_path, sep="\t", index=False, compression="gzip")
    coordinate_rows.sort_values("cell_id", kind="mergesort").to_csv(
        coordinate_path, sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame(changes).to_csv(change_path, sep="\t", index=False)
    before = Counter(original.final_broad_label)
    after = Counter(revised.final_broad_label)
    transitions = Counter((row["old_broad_label"], row["new_broad_label"]) for row in changes)
    manifest = {
        "status": "PENDING_POST_REPAIR_BIOLOGICAL_REVIEW",
        "schema_version": "2.2", "stage": "follicle_roi_repair_apply",
        "artifact_role": "review_candidate",
        "stage_authority": artifact(args.stage_authority),
        "pre_repair_membership": artifact(args.membership),
        "repaired_membership": artifact(membership_path),
        "combined_observation_scores": artifact(score_path),
        "coordinate_membership": artifact(coordinate_path),
        "changes": artifact(change_path),
        "repair_ancestry": [
            artifact(path) for _, path in sorted(repair_ancestry.items())
        ],
        "candidate_catalog": artifact(args.catalog),
        "context_evidence": artifact(args.context_evidence) if args.context_evidence else None,
        "context_release_eligibility": context_summary,
        "repair_rois": repair_rois,
        "skipped_non_wall_review_rois": skipped_non_wall_rois,
        "typed_layer_writeback_authority": typed_authority,
        "n_changed_observations": len(changes),
        "before_broad_census": dict(sorted(before.items())),
        "after_broad_census": dict(sorted(after.items())),
        "changed_transitions": {f"{a} -> {b}": n for (a, b), n in sorted(transitions.items())},
        "assignment_rule": "direct identity core; release-family coherence; pairwise normalized-evidence discriminator; optional Stromal exact remainder",
        "minimum_specific_lineage_margin": args.specific_margin,
        "structural_theca_used_as_formal_theca": False,
        "geometry_used_as_assignment_authority": False,
        "outside_roi_label_membership_unchanged": True,
        "protected_non_wall_membership_unchanged": True,
        "semantic_sha256": deterministic_membership_hash(revised.reset_index(drop=True).to_dict("records")),
    }
    manifest_path = args.out / "follicle_roi_repair_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
