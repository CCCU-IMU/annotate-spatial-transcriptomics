#!/usr/bin/env python3
"""Register one immutable, evidence-validated broad-class or targeted cohort."""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from evidence_schema_lib import read_tsv, sha256
from validate_cohort_outcome import validate as validate_outcome


FIELDS = [
    "cohort_id", "sample_id", "cohort_type", "question_mode", "source_broad_label",
    "source_run_ids", "source_cluster_ids", "membership_path", "membership_sha256",
    "n_observations", "purpose", "competing_hypotheses", "candidate_resolutions",
    "selected_resolution", "terminal_outcome", "homogeneous_parent_confirmed",
    "applicability", "applicability_rationale", "underpowered_skip_artifact",
    "underpowered_skip_artifact_sha256", "outcome_artifact", "outcome_artifact_sha256",
    "status", "supersedes", "created_at", "closed_at",
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--cohort-type", required=True, choices=["broad_class_recluster", "targeted_recluster", "oocyte_targeted_recluster"])
    parser.add_argument("--question-mode", required=True, choices=["broad_purity_audit", "targeted_mixture"])
    parser.add_argument("--source-broad-label", default="")
    parser.add_argument("--source-run-ids", required=True)
    parser.add_argument("--source-cluster-ids", required=True)
    parser.add_argument("--membership", type=Path)
    parser.add_argument("--purpose", required=True)
    parser.add_argument("--competing-hypotheses", default="")
    parser.add_argument("--candidate-resolutions", default="")
    parser.add_argument("--selected-resolution", default="")
    parser.add_argument("--terminal-outcome", choices=["homogeneous_parent_confirmed", "subclusters_adjudicated"], default="")
    parser.add_argument("--applicability", choices=["applicable", "not_applicable"], required=True)
    parser.add_argument("--applicability-rationale", required=True)
    parser.add_argument("--underpowered-skip-artifact", type=Path)
    parser.add_argument("--outcome-artifact", type=Path)
    parser.add_argument("--status", choices=["validated_done", "not_applicable_reviewed"], required=True)
    parser.add_argument("--supersedes", default="", help="one prior cohort_id replaced by this immutable successor")
    args = parser.parse_args()
    root = args.project_root.resolve()
    project = json.loads((root / "config/project.json").read_text(encoding="utf-8"))
    if args.cohort_type == "broad_class_recluster" and args.question_mode != "broad_purity_audit":
        raise SystemExit("broad-class cohorts must use broad_purity_audit")
    if args.cohort_type != "broad_class_recluster" and args.question_mode != "targeted_mixture":
        raise SystemExit("targeted cohorts must use targeted_mixture")

    membership_path = ""
    membership_hash = ""
    n_observations = 0
    outcome_hash = ""
    skip_hash = ""
    if args.applicability == "applicable":
        if args.status != "validated_done" or not args.membership or not args.membership.is_file():
            raise SystemExit("applicable cohort requires validated_done and a membership artifact")
        members = read_tsv(args.membership)
        ids = [row.get("cell_id", "") for row in members]
        if not members or "cell_id" not in members[0] or "" in ids or len(ids) != len(set(ids)):
            raise SystemExit("membership must contain unique nonempty cell_id")
        if not args.outcome_artifact or not args.outcome_artifact.is_file():
            raise SystemExit("applicable cohort requires a cohort-outcome artifact")
        membership_path = str(args.membership.resolve())
        membership_hash = sha256(args.membership)
        n_observations = len(ids)
        validation = validate_outcome(root, args.outcome_artifact.resolve())
        if validation.get("status") != "PASS":
            raise SystemExit("invalid cohort outcome: " + "; ".join(validation.get("errors", [])))
        outcome = json.loads(args.outcome_artifact.read_text(encoding="utf-8"))
        if outcome.get("cohort_id") != args.cohort_id or outcome.get("query_membership", {}).get("sha256") != membership_hash:
            raise SystemExit("cohort outcome does not bind the registered cohort/query")
        if outcome.get("cohort_type") != args.cohort_type or outcome.get("question_mode") != args.question_mode:
            raise SystemExit("cohort outcome type/question mode differs from registration")
        try:
            registered_grid = [float(value) for value in args.candidate_resolutions.replace(";", ",").split(",") if value.strip()]
            selected = float(args.selected_resolution)
        except ValueError:
            raise SystemExit("candidate resolutions and selected resolution must be numeric")
        if registered_grid != [float(value) for value in outcome.get("candidate_grid", [])]:
            raise SystemExit("registered candidate grid differs from the evidence artifact")
        if selected != float(outcome.get("selected_resolution")):
            raise SystemExit("registered selected resolution differs from the evidence artifact")
        if args.terminal_outcome != outcome.get("terminal_outcome"):
            raise SystemExit("terminal outcome differs from the evidence artifact")
        outcome_hash = sha256(args.outcome_artifact)
    else:
        if args.status != "not_applicable_reviewed" or len(args.applicability_rationale.strip()) < 20:
            raise SystemExit("underpowered skip requires reviewed terminal status and substantive rationale")
        if not args.underpowered_skip_artifact or not args.underpowered_skip_artifact.is_file() or args.underpowered_skip_artifact.stat().st_size == 0:
            raise SystemExit("underpowered skip requires a nonempty formal artifact")
        skip_hash = sha256(args.underpowered_skip_artifact)

    registry = root / "state/recluster_cohort_registry.tsv"
    old = read_tsv(registry)
    if any(row.get("cohort_id") == args.cohort_id for row in old):
        raise SystemExit("cohort_id already exists; create a versioned successor")
    if args.supersedes:
        matches = [row for row in old if row.get("cohort_id") == args.supersedes]
        already = {value for row in old for value in row.get("supersedes", "").replace("|", ",").split(",") if value}
        if len(matches) != 1 or args.supersedes in already:
            raise SystemExit("--supersedes must name one active prior cohort")
    now = datetime.now(timezone.utc).isoformat()
    row = {
        "cohort_id": args.cohort_id, "sample_id": project["sample_id"], "cohort_type": args.cohort_type,
        "question_mode": args.question_mode, "source_broad_label": args.source_broad_label,
        "source_run_ids": args.source_run_ids, "source_cluster_ids": args.source_cluster_ids,
        "membership_path": membership_path, "membership_sha256": membership_hash, "n_observations": n_observations,
        "purpose": args.purpose, "competing_hypotheses": args.competing_hypotheses,
        "candidate_resolutions": args.candidate_resolutions, "selected_resolution": args.selected_resolution,
        "terminal_outcome": args.terminal_outcome,
        "homogeneous_parent_confirmed": str(args.terminal_outcome == "homogeneous_parent_confirmed").lower(),
        "applicability": args.applicability, "applicability_rationale": args.applicability_rationale,
        "underpowered_skip_artifact": str(args.underpowered_skip_artifact.resolve()) if args.underpowered_skip_artifact else "",
        "underpowered_skip_artifact_sha256": skip_hash,
        "outcome_artifact": str(args.outcome_artifact.resolve()) if args.outcome_artifact else "",
        "outcome_artifact_sha256": outcome_hash, "status": args.status, "supersedes": args.supersedes,
        "created_at": now, "closed_at": now,
    }
    with registry.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(old + [row])
    print(json.dumps(row, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
