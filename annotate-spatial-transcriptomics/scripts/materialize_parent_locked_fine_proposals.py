#!/usr/bin/env python3
"""Materialize only high-confidence second-round fine proposals under frozen parents."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from controller_thresholds import load_controller_thresholds
from evidence_schema_lib import sha256
from lineage_controller_lib import (
    apply_candidate_context, candidate_can_release, catalog_candidates, number,
    read_tsv, write_tsv,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument("--threshold-registry", required=True, type=Path)
    ap.add_argument("--fine-audit", action="append", type=Path, default=[])
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    thresholds = load_controller_thresholds(args.threshold_registry)
    fine_policy = thresholds["fine_release_policy"]
    contradiction_ceiling = thresholds[
        "observation_writeback_policy"
    ]["maximum_contradiction_fraction"]

    membership = read_tsv(args.membership)
    by_group: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in membership:
        key = (str(row.get("source_boundary", "")), str(row.get("source_cluster", "")))
        by_group.setdefault(key, []).append(row)

    candidates = catalog_candidates(
        json.loads(args.catalog.read_text(encoding="utf-8"))
    )
    context_summary = apply_candidate_context(
        candidates,
        read_tsv(args.context_evidence) if args.context_evidence else [],
    )
    assignments: dict[str, dict[str, object]] = {}
    source_records: list[dict[str, str]] = []
    context_skipped_records: list[dict[str, str]] = []
    for path in args.fine_audit:
        for row in read_tsv(path):
            if row.get("status") != "supported" or row.get("release_candidate") != "true":
                continue
            candidate_id = str(row.get("candidate_id", ""))
            candidate = candidates.get(candidate_id, {})
            if not candidate_can_release(candidate):
                context_skipped_records.append(row)
                continue
            parent = str(candidate.get("release_broad_label", ""))
            fine = str(candidate.get("release_fine_label", ""))
            if (
                candidate.get("candidate_role") != "fine" or not parent or not fine
                or parent != row.get("parent_broad_label")
                or number(row.get("lineage_supported_fraction"))
                < fine_policy["minimum_lineage_supported_fraction"]
                or number(row.get("lineage_supported_fraction"))
                - number(row.get("strongest_competing_fraction"))
                < fine_policy["minimum_purity_margin"]
                or number(row.get("contradiction_fraction"))
                > contradiction_ceiling
            ):
                raise SystemExit("supported fine proposal violates parent/purity contract")
            key = (str(row.get("cohort_id", "")), str(row.get("subcluster_id", "")))
            members = by_group.get(key, [])
            if not members:
                raise SystemExit(f"fine proposal source subcluster is absent: {key}")
            for member in members:
                if member.get("final_broad_label") != parent:
                    continue
                cell = str(member.get("cell_id", ""))
                proposal = {
                    "cell_id": cell,
                    "parent_broad_label": parent,
                    "final_fine_label": fine,
                    "confidence": "high",
                    "candidate_id": candidate_id,
                    "source_partition": f"{key[0]}::{key[1]}",
                    "assignment_mode": "second_round_parent_locked_whole_subcluster",
                }
                previous = assignments.get(cell)
                if previous and previous["final_fine_label"] != fine:
                    raise SystemExit("one observation received competing supported fine labels")
                assignments[cell] = proposal
            source_records.append(row)

    rows = [assignments[cell] for cell in sorted(assignments)]
    args.out.mkdir(parents=True, exist_ok=True)
    assignment_path = args.out / "parent_locked_fine_assignments.tsv.gz"
    write_tsv(
        assignment_path, rows,
        fields=[
            "cell_id", "parent_broad_label", "final_fine_label", "confidence",
            "candidate_id", "source_partition", "assignment_mode",
        ],
    )
    manifest = {
        "status": "PASS",
        "schema_version": "2.2",
        "stage": "post_broad_freeze_fine_materialization",
        "threshold_registry": {
            "path": str(args.threshold_registry.resolve()),
            "sha256": sha256(args.threshold_registry),
        },
        "broad_labels_modified": False,
        "n_fine_assignments": len(rows),
        "n_supported_sources": len(source_records),
        "n_context_ineligible_supported_sources_skipped": len(
            context_skipped_records
        ),
        "context_release_eligibility": context_summary,
        "context_evidence": (
            {
                "path": str(args.context_evidence.resolve()),
                "sha256": sha256(args.context_evidence),
            }
            if args.context_evidence else None
        ),
        "membership": {
            "path": str(args.membership.resolve()),
            "sha256": sha256(args.membership),
        },
        "assignments": {
            "path": str(assignment_path.resolve()),
            "sha256": sha256(assignment_path),
        },
    }
    manifest_path = args.out / "parent_locked_fine_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
