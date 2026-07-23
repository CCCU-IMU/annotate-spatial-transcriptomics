#!/usr/bin/env python3
"""Apply validated fine assignments without changing locked broad labels.

This writer is deliberately narrow: broad labels are immutable, every new fine
label must be high-confidence and catalog-bound, and all rows carrying a fine
label are normalized to a self-consistent defined_fine state.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


TRUE = {"1", "true", "yes", "high"}


def opener(path: Path, mode: str):
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(mode, encoding="utf-8", newline="")


def fine_catalog(path: Path) -> dict[str, str]:
    doc = json.loads(path.read_text(encoding="utf-8"))
    result: dict[str, str] = {}
    for parent_id, items in doc.get("machine_actionable_fine_candidate_catalog", {}).items():
        if not isinstance(items, list):
            continue
        for item in items:
            if not isinstance(item, dict):
                continue
            label = str(item.get("release_label", "")).strip()
            parent = str(item.get("parent_release_label", "")).strip()
            if not label or not parent:
                raise ValueError(f"fine candidate under {parent_id} lacks release_label or parent_release_label")
            if label in result and result[label] != parent:
                raise ValueError(f"fine release label maps to multiple parents: {label}")
            result[label] = parent
    return result


def read_assignments(path: Path, catalog: dict[str, str]) -> dict[str, dict[str, str]]:
    with opener(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"cell_id", "parent_broad_label", "final_fine_label", "confidence"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError("assignment table lacks required columns: " + ", ".join(sorted(required)))
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            cell = str(row["cell_id"]).strip()
            label = str(row["final_fine_label"]).strip()
            parent = str(row["parent_broad_label"]).strip()
            if not cell or cell in rows:
                raise ValueError(f"empty or duplicate assignment cell_id: {cell!r}")
            if str(row["confidence"]).strip().lower() != "high":
                raise ValueError(f"fine assignment is not high-confidence: {cell}")
            if catalog.get(label) != parent:
                raise ValueError(f"fine label/parent is not catalog-bound: {parent} -> {label}")
            rows[cell] = row
    if not rows:
        raise ValueError("assignment table is empty")
    return rows


def decision_id(sample: str, cell: str, version: str) -> str:
    token = hashlib.sha256(f"{sample}\0{cell}\0{version}".encode()).hexdigest()[:20]
    return f"{sample}_fine_{version.replace('.', '_').replace('-', '_')}_{token}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--assignments", required=True, type=Path)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--decision-version", default="2.1.0")
    parser.add_argument("--source-run-id", required=True)
    parser.add_argument("--validation-artifact", required=True, type=Path)
    args = parser.parse_args()

    catalog = fine_catalog(args.catalog)
    assignments = read_assignments(args.assignments, catalog)
    if not args.validation_artifact.is_file():
        raise FileNotFoundError(args.validation_artifact)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat()
    seen: set[str] = set()
    broad_before: Counter[str] = Counter()
    broad_after: Counter[str] = Counter()
    fine_after: Counter[str] = Counter()
    changed = 0
    repaired_existing_state = 0

    with opener(args.ledger, "rt") as src, opener(args.out, "wt") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("ledger header is missing")
        required = {
            "cell_id", "sample_id", "decision_id", "generation", "state", "broad_label", "fine_label",
            "fine_anchor_eligible", "decision_version", "supersedes", "final_state", "final_broad_label",
            "final_fine_label", "final_fine_confidence", "final_fine_eligible",
        }
        if not required <= set(reader.fieldnames):
            raise ValueError("ledger lacks required fine-state columns: " + ", ".join(sorted(required - set(reader.fieldnames))))
        writer = csv.DictWriter(dst, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            cell = str(row["cell_id"]).strip()
            broad = str(row.get("final_broad_label", "")).strip()
            broad_before[broad] += 1
            assignment = assignments.get(cell)
            existing_fine = str(row.get("final_fine_label", "")).strip()
            needs_existing_repair = bool(existing_fine) and (
                row.get("final_state") != "defined_fine"
                or row.get("state") != "defined_fine"
                or str(row.get("final_fine_confidence", "")).lower() != "high"
                or str(row.get("final_fine_eligible", "")).lower() not in TRUE
                or str(row.get("fine_anchor_eligible", "")).lower() not in TRUE
            )
            if assignment:
                seen.add(cell)
                if broad != str(assignment["parent_broad_label"]).strip():
                    raise ValueError(f"assignment crosses locked broad parent for {cell}: {broad}")
                new_fine = str(assignment["final_fine_label"]).strip()
                row["fine_label"] = new_fine
                row["final_fine_label"] = new_fine
                row["assignment_mode"] = assignment.get("assignment_mode", "stable_parent_cohort_fine_direct") or "stable_parent_cohort_fine_direct"
                row["recluster_cohort_id"] = assignment.get("source_partition", row.get("recluster_cohort_id", ""))
                changed += 1
            elif needs_existing_repair:
                repaired_existing_state += 1
                changed += 1
            if assignment or needs_existing_repair:
                old_decision = row["decision_id"]
                row["decision_id"] = decision_id(row["sample_id"], cell, args.decision_version)
                row["supersedes"] = old_decision
                try:
                    row["generation"] = str(int(float(row.get("generation", "0") or 0)) + 1)
                except ValueError:
                    row["generation"] = "1"
                row["state"] = "defined_fine"
                row["final_state"] = "defined_fine"
                row["fine_confidence"] = "high"
                row["final_fine_confidence"] = "high"
                row["fine_anchor_eligible"] = "true"
                row["final_fine_eligible"] = "true"
                row["evidence_status"] = "validated"
                row["validation_status"] = "passed"
                row["validation_artifact"] = str(args.validation_artifact.resolve())
                row["validation_feature_scope"] = "full_feature"
                row["source_run_id"] = args.source_run_id
                row["route"] = "complete_parent_fine_candidate_audit"
                row["decision_version"] = args.decision_version
                row["closed"] = "true"
                row["next_action"] = "none"
                row["created_at"] = now
            if row.get("final_fine_label"):
                fine_after[row["final_fine_label"]] += 1
            broad_after[str(row.get("final_broad_label", "")).strip()] += 1
            writer.writerow(row)

    missing = sorted(set(assignments) - seen)
    if missing:
        args.out.unlink(missing_ok=True)
        raise ValueError(f"assignment cells absent from ledger: {len(missing)}")
    if broad_before != broad_after:
        args.out.unlink(missing_ok=True)
        raise ValueError("locked broad-label census changed during fine writeback")

    manifest = {
        "status": "PASS",
        "decision_version": args.decision_version,
        "source_run_id": args.source_run_id,
        "n_new_assignments": len(assignments),
        "n_existing_fine_state_repairs": repaired_existing_state,
        "n_changed_rows": changed,
        "broad_labels_modified": False,
        "broad_census": dict(sorted(broad_after.items())),
        "fine_census": dict(sorted(fine_after.items())),
        "validation_artifact": str(args.validation_artifact.resolve()),
        "output_ledger": str(args.out.resolve()),
        "created_at": now,
    }
    args.manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
