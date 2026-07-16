#!/usr/bin/env python3
"""Register a direct parent, fine or cross-lineage return with validated evidence."""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from evidence_schema_lib import read_tsv, sha256
from validate_direct_return_evidence import validate as validate_evidence


FIELDS = [
    "return_id", "sample_id", "source_cohort_id", "source_run_id", "source_cluster",
    "membership_path", "membership_sha256", "n_observations", "target_broad_label",
    "target_fine_label", "confidence", "assignment_mode", "rctd_tier", "independent_evidence",
    "evidence_artifact", "evidence_artifact_sha256", "fine_anchor_eligible", "status", "supersedes", "created_at",
]


def truth(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--return-id", required=True)
    parser.add_argument("--source-cohort-id", required=True)
    parser.add_argument("--source-run-id", required=True)
    parser.add_argument("--source-cluster", required=True)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--target-broad-label", required=True)
    parser.add_argument("--target-fine-label", default="")
    parser.add_argument("--confidence", required=True, choices=["moderate", "high"])
    parser.add_argument("--assignment-mode", required=True, choices=["parent_broad_direct", "fine_direct", "cross_lineage_direct", "rctd_assisted"])
    parser.add_argument("--rctd-tier", default="")
    parser.add_argument("--independent-evidence", default="false")
    parser.add_argument("--evidence-artifact", required=True, type=Path)
    parser.add_argument("--fine-anchor-eligible", default="false")
    parser.add_argument("--supersedes", default="", help="one prior return_id replaced by this immutable successor")
    args = parser.parse_args()
    root = args.project_root.resolve()
    if not args.membership.is_file() or not args.evidence_artifact.is_file():
        raise SystemExit("membership and evidence artifacts must exist")
    members = read_tsv(args.membership)
    ids = [row.get("cell_id", "") for row in members]
    if not members or "cell_id" not in members[0] or "" in ids or len(ids) != len(set(ids)):
        raise SystemExit("membership must contain unique nonempty cell_id")
    if args.target_fine_label and (args.confidence != "high" or not truth(args.fine_anchor_eligible)):
        raise SystemExit("fine return requires high confidence and fine-anchor eligibility")
    if args.assignment_mode == "rctd_assisted":
        if args.rctd_tier not in {"low", "moderate", "high"}:
            raise SystemExit("new RCTD-assisted returns require canonical low/moderate/high tier evidence")
        if args.target_fine_label and (args.rctd_tier != "high" or not truth(args.independent_evidence)):
            raise SystemExit("RCTD-assisted fine return requires high tier plus independent evidence")
        if not args.target_fine_label and args.rctd_tier == "low":
            raise SystemExit("low-tier RCTD evidence must return to terminal residual QC")
    validation = validate_evidence(root, args.evidence_artifact.resolve())
    if validation.get("status") != "PASS":
        raise SystemExit("invalid direct-return evidence: " + "; ".join(validation.get("errors", [])))
    evidence = json.loads(args.evidence_artifact.read_text(encoding="utf-8"))
    membership_hash = sha256(args.membership)
    expected = {
        "return_id": args.return_id, "source_cohort_id": args.source_cohort_id,
        "source_subcluster": args.source_cluster, "target_broad_label": args.target_broad_label,
        "target_fine_label": args.target_fine_label, "confidence": args.confidence,
    }
    if any(str(evidence.get(key, "")) != value for key, value in expected.items()) or evidence.get("membership", {}).get("sha256") != membership_hash:
        raise SystemExit("direct-return evidence does not bind the registered return/membership")
    project = json.loads((root / "config/project.json").read_text(encoding="utf-8"))
    registry = root / "state/direct_return_registry.tsv"
    old = read_tsv(registry)
    if any(row.get("return_id") == args.return_id for row in old):
        raise SystemExit("return_id already exists; create a versioned successor")
    if args.supersedes:
        matches = [row for row in old if row.get("return_id") == args.supersedes]
        already = {value for row in old for value in row.get("supersedes", "").replace("|", ",").split(",") if value}
        if len(matches) != 1 or args.supersedes in already:
            raise SystemExit("--supersedes must name one active prior direct return")
    row = {
        "return_id": args.return_id, "sample_id": project["sample_id"], "source_cohort_id": args.source_cohort_id,
        "source_run_id": args.source_run_id, "source_cluster": args.source_cluster,
        "membership_path": str(args.membership.resolve()), "membership_sha256": membership_hash,
        "n_observations": len(ids), "target_broad_label": args.target_broad_label,
        "target_fine_label": args.target_fine_label, "confidence": args.confidence,
        "assignment_mode": args.assignment_mode, "rctd_tier": args.rctd_tier,
        "independent_evidence": args.independent_evidence,
        "evidence_artifact": str(args.evidence_artifact.resolve()),
        "evidence_artifact_sha256": sha256(args.evidence_artifact),
        "fine_anchor_eligible": args.fine_anchor_eligible, "status": "validated_done", "supersedes": args.supersedes,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    with registry.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(old + [row])
    print(json.dumps(row, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
