#!/usr/bin/env python3
"""Calibrate strict rare-cell seeds without shrinking the full candidate pool."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin(["true", "1"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--screen", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--positive-quantile", type=float, default=0.75)
    parser.add_argument("--contradictory-quantile", type=float, default=0.25)
    parser.add_argument("--modules-required", type=int, default=3)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(args.screen, sep="\t", dtype={args.cell_id_col: str})
    gate = as_bool(data["starting_marker_gate"])
    candidates = data.loc[gate].copy()
    required = {
        "spatial_focus_supported",
        "spatial_focus_id",
        "total_oocyte_program_hits",
        "identity_core_hits",
        "modules_detected",
        "contradictory_somatic_hits",
    }
    missing = required - set(candidates.columns)
    if missing:
        raise SystemExit(f"screen lacks columns: {sorted(missing)}")
    if len(candidates) < 10:
        raise SystemExit("too few starting-gate observations for query-calibrated quantiles")

    total_threshold = float(
        candidates["total_oocyte_program_hits"].quantile(args.positive_quantile)
    )
    identity_threshold = float(
        candidates["identity_core_hits"].quantile(args.positive_quantile)
    )
    contradiction_threshold = float(
        candidates["contradictory_somatic_hits"].quantile(
            args.contradictory_quantile
        )
    )
    spatial_supported = as_bool(candidates["spatial_focus_supported"])
    strict_seed = (
        spatial_supported
        & candidates["total_oocyte_program_hits"].ge(total_threshold)
        & candidates["identity_core_hits"].ge(identity_threshold)
        & candidates["modules_detected"].ge(args.modules_required)
        & candidates["contradictory_somatic_hits"].le(
            contradiction_threshold
        )
    )
    seed_groups = set(candidates.loc[strict_seed, "spatial_focus_id"])
    seed_group_expanded = candidates["spatial_focus_id"].isin(seed_groups)

    candidates["strict_seed_calibrated"] = strict_seed.to_numpy()
    candidates["strict_seed_group_expanded"] = seed_group_expanded.to_numpy()
    candidates["full_targeted_cohort_member"] = True
    candidates["targeted_cohort_role"] = "full_starting_gate_targeted_recluster_cohort"

    # Canonical R-first membership: the complete starting gate is reclustered.
    candidates.to_csv(
        args.out / "oocyte_targeted_recluster_cohort.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    # Backward-compatible seed-focus artifact.  This is supporting evidence,
    # not the canonical recluster membership or a final rare-cell call.
    focus_support = candidates.loc[seed_group_expanded].copy()
    focus_support.to_csv(
        args.out / "calibrated_rare_focus_support.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    summary = {
        "status": "CALIBRATED_TWO_TIER_CANDIDATE_ROUTE",
        "starting_gate": len(candidates),
        "full_targeted_recluster_cohort": len(candidates),
        "spatial_supported": int(spatial_supported.sum()),
        "positive_quantile": args.positive_quantile,
        "contradictory_quantile": args.contradictory_quantile,
        "thresholds": {
            "total_program_hits": total_threshold,
            "identity_core_hits": identity_threshold,
            "contradictory_hits": contradiction_threshold,
            "modules_required": args.modules_required,
        },
        "strict_seeds": int(strict_seed.sum()),
        "strict_seed_groups": len(seed_groups),
        "strict_seed_group_expanded": len(focus_support),
        "canonical_recluster_membership": "oocyte_targeted_recluster_cohort.tsv.gz",
        "warning": (
            "The complete multi-module starting gate is the canonical query-only "
            "targeted recluster cohort. Strict seeds and expanded spatial objects are support, "
            "not the final census or a membership filter. Observations are not "
            "biological-cell counts."
        ),
    }
    (args.out / "calibrated_rare_focus_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
