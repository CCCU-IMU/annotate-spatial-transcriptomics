#!/usr/bin/env python3
"""Apply only calibrated broad rescue to unlabeled post-merge observations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from evidence_schema_lib import sha256
from lineage_controller_lib import (
    apply_candidate_context, candidate_can_release, catalog_candidates,
    deterministic_membership_hash, read_tsv, write_tsv,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--frozen-broad", required=True, type=Path)
    ap.add_argument("--routing", required=True, type=Path)
    ap.add_argument("--atlas-validation", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    validation = json.loads(args.atlas_validation.read_text(encoding="utf-8"))
    if validation.get("status") != "PASS":
        raise SystemExit("Atlas review queue is not fully closed")
    candidates = catalog_candidates(
        json.loads(args.catalog.read_text(encoding="utf-8"))
    )
    context_summary = apply_candidate_context(
        candidates,
        read_tsv(args.context_evidence) if args.context_evidence else [],
    )
    context_eligible_broad = {
        str(candidate.get("release_broad_label", ""))
        for candidate in candidates.values()
        if str(candidate.get("candidate_role", "")) == "broad"
        and candidate_can_release(candidate)
    }
    frozen_rows = read_tsv(args.frozen_broad)
    route_rows = read_tsv(args.routing)
    frozen = {str(row.get("cell_id", "")): row for row in frozen_rows}
    routes = {str(row.get("cell_id", "")): row for row in route_rows}
    if (
        not frozen or "" in frozen or len(frozen) != len(frozen_rows)
        or set(frozen) != set(routes) or len(routes) != len(route_rows)
    ):
        raise SystemExit("Atlas routing must cover frozen broad membership exactly once")

    output: list[dict[str, object]] = []
    rescued = 0
    for cell in sorted(frozen):
        row = dict(frozen[cell])
        route = routes[cell]
        current_broad = str(row.get("final_broad_label", ""))
        proposed = str(route.get("proposed_broad_label", ""))
        route_name = str(route.get("atlas_state_route", ""))
        if current_broad:
            if proposed and proposed != current_broad:
                raise SystemExit("Atlas attempted to overwrite frozen broad membership")
        elif route_name in {"direct_unlabeled_broad_return", "direct_qc_broad_return"}:
            if not proposed:
                raise SystemExit("direct Atlas rescue lacks a proposed broad label")
            if proposed not in context_eligible_broad:
                raise SystemExit(
                    "direct Atlas rescue targets a context-ineligible broad label"
                )
            row["final_state"] = "defined_broad_only"
            row["final_broad_label"] = proposed
            row["confidence"] = "moderate"
            row["assignment_origin"] = "post_merge_atlas_unlabeled_broad_rescue"
            row["unresolved_reason"] = ""
            rescued += 1
        row["atlas_state_route"] = route_name
        row["atlas_broad"] = route.get("atlas_broad", "")
        row["atlas_tier"] = route.get("atlas_tier", "")
        row["atlas_review_id"] = route.get("review_id", "")
        output.append(row)

    args.out.mkdir(parents=True, exist_ok=True)
    membership_path = args.out / "post_atlas_broad_membership.tsv.gz"
    write_tsv(membership_path, output)
    manifest = {
        "status": "PASS",
        "schema_version": "2.2",
        "stage": "atlas_and_completeness_review",
        "formal_broad_overwrite_performed": False,
        "unlabeled_broad_rescue_n": rescued,
        "n_unresolved_biological": sum(
            row.get("final_state") == "unresolved_biological" for row in output
        ),
        "membership": {
            "path": str(membership_path.resolve()),
            "sha256": sha256(membership_path),
            "semantic_sha256": deterministic_membership_hash(output),
        },
        "frozen_broad": {
            "path": str(args.frozen_broad.resolve()),
            "sha256": sha256(args.frozen_broad),
        },
        "routing": {"path": str(args.routing.resolve()), "sha256": sha256(args.routing)},
        "atlas_validation": {
            "path": str(args.atlas_validation.resolve()),
            "sha256": sha256(args.atlas_validation),
        },
        "candidate_catalog": {
            "path": str(args.catalog.resolve()), "sha256": sha256(args.catalog),
        },
        "context_evidence": (
            {"path": str(args.context_evidence.resolve()), "sha256": sha256(args.context_evidence)}
            if args.context_evidence else None
        ),
        "context_release_eligibility": context_summary,
    }
    manifest_path = args.out / "post_atlas_membership_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
