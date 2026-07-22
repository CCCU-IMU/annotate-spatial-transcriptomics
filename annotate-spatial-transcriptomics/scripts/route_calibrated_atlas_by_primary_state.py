#!/usr/bin/env python3
"""Route one calibrated all-cell Atlas result by the frozen primary state.

This router deliberately separates two questions:

* an unlabeled frozen-QC observation may receive a broad-only Atlas return;
* an observation with an existing biological label is only compared and may
  trigger review.  Atlas never silently overwrites that label.

Marker, anti-marker and spatial columns may be carried through as audit or
challenge evidence, but they are not a per-observation prerequisite for the
default moderate/high QC return.  OOD, ontology conflict and profile-excluded
Atlas classes remain hard safety gates.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path


TRUE = {"1", "true", "yes", "y", "pass", "passed"}


def truth(value: object) -> bool:
    return str(value).strip().lower() in TRUE


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--concordance", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    with args.concordance.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        input_fields = list(reader.fieldnames or [])
        frame = list(reader)
    columns = {
        "id": config.get("cell_id_column", "cell_id"),
        "state": config.get("primary_state_column", "primary_state"),
        "primary": config.get("primary_label_column", "primary_broad_label"),
        "atlas": config.get("atlas_label_column", "atlas_broad_label"),
        "tier": config.get("atlas_tier_column", "atlas_tier"),
        "ood": config.get("ood_column", "atlas_ood"),
        "conflict": config.get("ontology_conflict_column", "ontology_conflict"),
    }
    missing = sorted(set(columns.values()).difference(input_fields))
    if missing:
        raise SystemExit(f"concordance table lacks required columns: {missing}")
    cell_ids = [str(row.get(columns["id"], "")).strip() for row in frame]
    if any(not cell for cell in cell_ids) or len(cell_ids) != len(set(cell_ids)):
        raise SystemExit("cell IDs must be nonempty and unique")

    accepted_tiers = {
        str(value).strip().lower()
        for value in config.get("accepted_atlas_tiers", ["high", "moderate_only"])
    }
    qc_states = {
        str(value).strip().lower()
        for value in config.get(
            "qc_states", ["qc_holdout", "low_information_qc_holdout", "pending_qc"]
        )
    }
    excluded_labels = {
        str(value).strip().lower()
        for value in config.get("atlas_direct_return_excluded_labels", [])
    }

    output_rows: list[dict[str, object]] = []
    for row in frame:
        state = str(row[columns["state"]]).strip().lower()
        primary = str(row[columns["primary"]]).strip()
        atlas = str(row[columns["atlas"]]).strip()
        tier = str(row[columns["tier"]]).strip().lower()
        is_qc_unlabeled = state in qc_states and not primary
        calibrated = tier in accepted_tiers and bool(atlas)
        ood = truth(row[columns["ood"]])
        conflict = truth(row[columns["conflict"]])
        excluded = atlas.lower() in excluded_labels

        if is_qc_unlabeled:
            blockers = []
            if not calibrated:
                blockers.append("atlas_below_moderate_or_missing")
            if ood:
                blockers.append("atlas_ood")
            if conflict:
                blockers.append("ontology_conflict")
            if excluded:
                blockers.append("profile_excluded_atlas_class")
            if blockers:
                route = "retain_qc"
                proposed_state = state or "qc_holdout"
                proposed_label = ""
                reason = ";".join(blockers)
                requires_review = False
            else:
                route = "direct_qc_broad_return"
                proposed_state = "defined_broad_only"
                proposed_label = atlas
                reason = "calibrated_moderate_or_higher_unlabeled_qc_non_ood_no_conflict"
                requires_review = False
        elif primary:
            if calibrated and not ood and not conflict and atlas == primary:
                route = "defined_label_atlas_agreement"
                reason = "same_broad_label"
                requires_review = False
            elif calibrated and atlas and atlas != primary:
                route = "defined_label_disagreement_review"
                reason = "materiality_must_be_evaluated_at_cluster_or_cohort_level"
                requires_review = True
            elif ood or conflict:
                route = "defined_label_ood_or_ontology_review"
                reason = "atlas_ood" if ood else "ontology_conflict"
                requires_review = True
            else:
                route = "defined_label_low_atlas_logged"
                reason = "atlas_not_strong_enough_to_challenge"
                requires_review = False
            proposed_state = state
            proposed_label = primary
        else:
            route = "unlabeled_non_qc_state_review"
            proposed_state = state
            proposed_label = ""
            reason = "unlabeled_observation_is_outside_frozen_qc_scope"
            requires_review = True
        output = dict(row)
        output.update({
            "atlas_state_route": route,
            "proposed_state": proposed_state,
            "proposed_broad_label": proposed_label,
            "route_reason": reason,
            "review_required": str(requires_review).lower(),
            "proposed_fine_label": "",
            "fine_anchor_eligible": "false",
            "writeback_status": "proposal_only_requires_atomic_commit",
        })
        output_rows.append(output)

    args.out.mkdir(parents=True, exist_ok=True)
    routed = args.out / "atlas_state_routing.tsv.gz"
    extra_fields = ["atlas_state_route", "proposed_state", "proposed_broad_label", "route_reason",
                    "review_required", "proposed_fine_label", "fine_anchor_eligible", "writeback_status"]
    with gzip.open(routed, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, input_fields + extra_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    census = Counter(str(row["atlas_state_route"]) for row in output_rows)
    manifest = {
        "status": "PASS",
        "n_observations": len(output_rows),
        "route_census": {key: int(value) for key, value in census.items()},
        "direct_qc_broad_return_n": census.get("direct_qc_broad_return", 0),
        "defined_label_review_n": sum(row["review_required"] == "true" for row in output_rows),
        "broad_only": True,
        "fine_anchor_eligible": False,
        "marker_spatial_evidence_role": "audit_or_group_level_challenge_not_per_cell_prerequisite",
        "ledger_writeback_performed": False,
    }
    (args.out / "atlas_state_routing_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
