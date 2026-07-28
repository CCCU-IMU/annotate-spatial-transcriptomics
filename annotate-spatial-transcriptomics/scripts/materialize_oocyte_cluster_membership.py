#!/usr/bin/env python3
"""Materialize Oocyte membership from passing targeted-recluster clusters.

Strict seeds and spatial-object IDs are intentionally ignored here.  They are
cluster-adjudication evidence, not observation-level inclusion filters.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from lineage_controller_lib import (
    apply_candidate_context, candidate_can_release, catalog_candidates,
    read_tsv, sha256, write_tsv,
)


ALLOWED_EXCLUSION_CLASSES = {
    "direct_multifamily_somatic_hard_contradiction",
    "objective_input_qc",
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--canonical-membership", required=True, type=Path)
    ap.add_argument("--passing-clusters", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument("--explicit-exclusions", type=Path)
    ap.add_argument("--cell-id-column", default="cell_id")
    ap.add_argument("--cluster-column", default="recluster_cluster")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    candidates = catalog_candidates(json.loads(args.catalog.read_text(encoding="utf-8")))
    context_summary = apply_candidate_context(
        candidates, read_tsv(args.context_evidence) if args.context_evidence else []
    )
    oocyte_candidates = sorted(
        candidate_id for candidate_id, candidate in candidates.items()
        if str(candidate.get("candidate_role", "")) == "broad"
        and str(candidate.get("release_broad_label", "")) == "Oocyte"
        and candidate_can_release(candidate)
    )
    if not oocyte_candidates:
        raise SystemExit("canonical Oocyte materialization has no context-eligible broad candidate")

    canonical = read_tsv(args.canonical_membership)
    passing_rows = read_tsv(args.passing_clusters)
    if not canonical:
        raise SystemExit("canonical Oocyte membership is empty")
    required = {args.cell_id_column, args.cluster_column}
    missing = required.difference(canonical[0])
    if missing:
        raise SystemExit(f"canonical membership lacks columns: {sorted(missing)}")
    passing = {
        str(row.get(args.cluster_column, "")).strip()
        for row in passing_rows
        if str(row.get(args.cluster_column, "")).strip()
        and str(row.get("adjudication_status", "pass")).strip().lower()
        in {"pass", "passed", "supported", "oocyte"}
    }
    if not passing:
        raise SystemExit("no passing Oocyte recluster cluster was supplied")

    exclusions: dict[str, str] = {}
    if args.explicit_exclusions:
        for row in read_tsv(args.explicit_exclusions):
            cell = str(row.get(args.cell_id_column, "")).strip()
            exclusion_class = str(row.get("exclusion_class", "")).strip()
            if not cell or exclusion_class not in ALLOWED_EXCLUSION_CLASSES:
                raise SystemExit(
                    "explicit Oocyte exclusions require cell_id and an allowed "
                    "exclusion_class"
                )
            if cell in exclusions:
                raise SystemExit("explicit Oocyte exclusions contain duplicate cell IDs")
            exclusions[cell] = exclusion_class

    seen: set[str] = set()
    materialized: list[dict[str, str]] = []
    eligible_n = 0
    for row in canonical:
        cell = str(row[args.cell_id_column]).strip()
        cluster = str(row[args.cluster_column]).strip()
        if not cell or cell in seen:
            raise SystemExit("canonical membership cell IDs must be nonempty and unique")
        seen.add(cell)
        if cluster not in passing:
            continue
        eligible_n += 1
        if cell in exclusions:
            continue
        materialized.append(
            {
                args.cell_id_column: cell,
                args.cluster_column: cluster,
                "final_broad_label": "Oocyte",
                "decision_basis": "canonical_member_of_passing_oocyte_recluster_cluster",
                "fine_anchor_eligible": "false",
            }
        )

    args.out.mkdir(parents=True, exist_ok=True)
    table = args.out / "materialized_oocyte_membership.tsv"
    fields = [
        args.cell_id_column,
        args.cluster_column,
        "final_broad_label",
        "decision_basis",
        "fine_anchor_eligible",
    ]
    unknown_exclusions = sorted(set(exclusions).difference(seen))
    if unknown_exclusions:
        raise SystemExit("explicit Oocyte exclusions are outside the canonical cohort")
    passing_eligible = {
        str(row[args.cell_id_column]).strip()
        for row in canonical
        if str(row[args.cluster_column]).strip() in passing
    }
    outside_passing = sorted(set(exclusions).difference(passing_eligible))
    if outside_passing:
        raise SystemExit("explicit Oocyte exclusions are outside passing clusters")
    write_tsv(table, materialized, fields)
    manifest = {
        "status": "PASS",
        "schema_version": "2.2",
        "canonical_membership_sha256": sha256(args.canonical_membership),
        "passing_clusters_sha256": sha256(args.passing_clusters),
        "candidate_catalog_sha256": sha256(args.catalog),
        "context_evidence_sha256": sha256(args.context_evidence) if args.context_evidence else "",
        "context_release_eligibility": context_summary,
        "eligible_oocyte_candidate_ids": oocyte_candidates,
        "canonical_n": len(canonical),
        "passing_clusters": sorted(passing),
        "eligible_canonical_members_n": eligible_n,
        "explicit_exclusion_n": sum(cell in exclusions for cell in seen),
        "final_oocyte_n": len(materialized),
        "strict_seed_used_as_membership_filter": False,
        "spatial_object_used_as_membership_filter": False,
        "only_allowed_hard_exclusions_applied": True,
        "exclusion_classes": sorted(set(exclusions.values())),
    }
    (args.out / "materialized_oocyte_membership_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
