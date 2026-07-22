#!/usr/bin/env python3
"""Require an upstream recall audit when residual QC is unusually large."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

from evidence_schema_lib import validate_evidence_artifact, validate_json_against_schema


QC_STATES = {"qc_holdout", "low_information_qc_holdout", "pending_qc", "unknown_candidate"}
REQUIRED_REVIEWS = {
    "initial_broad_resolution_recall_review",
    "selected_plus_two_higher_catalog_scan",
    "large_label_embedded_program_review",
    "atlas_tier_census_review",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(newline="", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", required=True, type=Path)
    ap.add_argument("--audit", type=Path)
    ap.add_argument("--state-column", default="annotation_status")
    ap.add_argument("--qc-reason-column", default="qc_reason")
    ap.add_argument("--require-v2", action="store_true")
    ap.add_argument("--fraction-trigger", type=float, default=0.10)
    ap.add_argument("--count-trigger", type=int, default=50000)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    with open_text(args.ledger) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or args.state_column not in rows[0]:
        raise SystemExit(f"ledger is empty or lacks {args.state_column}")
    qc_n = sum(str(row.get(args.state_column, "")).strip().lower() in QC_STATES for row in rows)
    fraction = qc_n / len(rows)
    triggered = qc_n >= args.count_trigger or fraction >= args.fraction_trigger
    errors: list[str] = []
    audit_payload: dict = {}
    if triggered:
        if not args.audit or not args.audit.is_file():
            errors.append("large residual QC requires a hash-bound upstream broad-recall audit")
        else:
            audit_payload = json.loads(args.audit.read_text(encoding="utf-8"))
            if args.require_v2:
                schema = Path(__file__).resolve().parents[1] / "schemas/residual_qc_audit.schema.json"
                audit_payload, schema_errors = validate_json_against_schema(args.audit, schema)
                errors.extend(schema_errors)
            if audit_payload.get("status") not in {"PASS", "CLOSED"}:
                errors.append("residual-QC audit is not PASS/CLOSED")
            if audit_payload.get("cell_ledger_sha256") != sha256(args.ledger):
                errors.append("residual-QC audit is stale for the current cell ledger")
            if int(audit_payload.get("residual_qc_n", -1)) != qc_n:
                errors.append("residual-QC audit count does not match the ledger")
            reviews = set(audit_payload.get("completed_reviews", []))
            missing = sorted(REQUIRED_REVIEWS.difference(reviews))
            if missing:
                errors.append("residual-QC audit lacks required reviews: " + ", ".join(missing))
            if audit_payload.get("unresolved_query_lineage_signals", 1) != 0:
                errors.append("residual-QC audit still has unresolved query-derived lineage signals")
            if args.require_v2:
                reason_census: dict[str, int] = {}
                if args.qc_reason_column not in rows[0]:
                    errors.append(f"v2 residual-QC ledger lacks {args.qc_reason_column}")
                else:
                    for row in rows:
                        if str(row.get(args.state_column, "")).strip().lower() not in QC_STATES:
                            continue
                        reason = str(row.get(args.qc_reason_column, "")).strip()
                        if not reason:
                            errors.append("v2 residual-QC observation lacks a typed qc_reason")
                            continue
                        reason_census[reason] = reason_census.get(reason, 0) + 1
                    declared = {str(k): int(v) for k, v in audit_payload.get("residual_qc_reason_census", {}).items()}
                    if reason_census != declared:
                        errors.append("residual-QC reason census differs from the current ledger")
                for index, artifact in enumerate(audit_payload.get("evidence_artifacts", []), 1):
                    _, artifact_errors = validate_evidence_artifact(
                        args.audit.parent, artifact, f"residual-QC evidence {index}"
                    )
                    errors.extend(artifact_errors)

    result = {
        "status": "PASS" if not errors else "BLOCKED",
        "analysis_n": len(rows),
        "cell_ledger": str(args.ledger.resolve()),
        "cell_ledger_sha256": sha256(args.ledger),
        "residual_qc_n": qc_n,
        "residual_qc_fraction": fraction,
        "audit_triggered": triggered,
        "v2_evidence_required": args.require_v2,
        "count_trigger": args.count_trigger,
        "fraction_trigger": args.fraction_trigger,
        "errors": errors,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
