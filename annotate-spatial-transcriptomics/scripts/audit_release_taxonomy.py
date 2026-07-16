#!/usr/bin/env python3
"""Fail-closed audit of release taxonomy versus computational cohorts and retained states."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path


BIOLOGICAL_STATUSES = {"defined_fine", "defined_broad_only", "broad_defined", "fine_defined", "broad_only"}
RETAINED_STATUSES = {
    "anatomical_interface",
    "interface_review",
    "low_information_qc_holdout",
    "qc_holdout",
    "technical_state",
    "pending_review",
    "unknown_candidate",
    "excluded_initial_qc",
}
NON_BIOLOGICAL_SUFFIXES = ("_review", "_candidate", "_unresolved", "_holdout", "_cohort")


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def dialect(path: Path) -> str:
    name = path.name[:-3] if path.name.endswith(".gz") else path.name
    return "," if name.endswith(".csv") else "\t"


def norm(value: str) -> str:
    return "_".join(value.strip().lower().replace("/", "_").replace("-", "_").split())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("metadata", type=Path)
    ap.add_argument("--profile", required=True, type=Path)
    ap.add_argument("--broad-column", default="broad_label")
    ap.add_argument("--cohort-column", default="recluster_cohort_id")
    ap.add_argument("--pool-column", help="legacy alias for projects being migrated")
    ap.add_argument("--status-column", default="annotation_status")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    profile = json.loads(args.profile.read_text(encoding="utf-8"))
    taxonomy = profile.get("release_taxonomy", {})
    forbidden = {norm(x) for x in taxonomy.get("forbidden_release_broad_labels", [])}
    retained_names = {norm(x) for x in taxonomy.get("non_biological_retained_states", [])}
    allowed = {norm(x) for x in taxonomy.get("default_biological_broad_classes", [])}
    allowed.update(norm(x) for x in taxonomy.get("evidence_dependent_standalone_classes", {}).keys())

    errors = []
    warnings = []
    biological = Counter()
    retained = Counter()
    statuses = Counter()
    copied_provenance = Counter()
    total = 0

    with open_text(args.metadata) as handle:
        reader = csv.DictReader(handle, delimiter=dialect(args.metadata))
        fields = set(reader.fieldnames or [])
        required = {args.broad_column, args.status_column}
        missing = sorted(required - fields)
        if missing:
            errors.append("missing required columns: " + ", ".join(missing))
        provenance_column = args.pool_column or args.cohort_column
        has_provenance = provenance_column in fields
        if not has_provenance:
            warnings.append(f"cohort provenance column absent: {provenance_column}; identifier-copy check skipped")

        for row in reader:
            total += 1
            label = (row.get(args.broad_column) or "").strip()
            status = (row.get(args.status_column) or "").strip()
            provenance_id = (row.get(provenance_column) or "").strip() if has_provenance else ""
            statuses[status or "<empty>"] += 1
            nlabel = norm(label)
            nprovenance = norm(provenance_id)

            if status in BIOLOGICAL_STATUSES:
                biological[label or "<empty>"] += 1
                if not label:
                    errors.append("biological row has empty broad label")
                if nlabel in forbidden:
                    errors.append(f"forbidden biological broad label: {label}")
                if nlabel in retained_names:
                    errors.append(f"retained-state name used as biological broad label: {label}")
                if allowed and nlabel not in allowed:
                    errors.append(f"biological broad label is outside the profile release vocabulary: {label}")
                if any(nlabel.endswith(s) for s in NON_BIOLOGICAL_SUFFIXES):
                    errors.append(f"workflow-state suffix used in biological broad label: {label}")
                if provenance_id and nlabel == nprovenance:
                    copied_provenance[label] += 1
            elif status in RETAINED_STATUSES:
                retained[label or status or "<empty>"] += 1
            else:
                retained[label or status or "<empty>"] += 1
                warnings.append(f"unrecognized annotation status: {status or '<empty>'}")

    for label, count in copied_provenance.items():
        errors.append(f"cohort/provenance identifier copied directly into biological label ({count} rows): {label}")

    result = {
        "schema_version": "1.0",
        "profile_id": profile.get("profile_id"),
        "metadata": str(args.metadata.resolve()),
        "metadata_sha256": sha256(args.metadata),
        "profile": str(args.profile.resolve()),
        "profile_sha256": sha256(args.profile),
        "total_observations": total,
        "biological_broad_census": dict(sorted(biological.items())),
        "retained_state_census": dict(sorted(retained.items())),
        "annotation_status_census": dict(sorted(statuses.items())),
        "errors": sorted(set(errors)),
        "warnings": sorted(set(warnings)),
        "pass": not errors,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["pass"] else 2


if __name__ == "__main__":
    sys.exit(main())
