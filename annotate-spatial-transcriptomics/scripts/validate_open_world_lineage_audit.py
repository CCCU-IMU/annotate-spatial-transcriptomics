#!/usr/bin/env python3
"""Validate a current-query, literature-informed open-world lineage review."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


OUTCOMES = {"supported_candidate", "not_supported", "not_evaluable"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dotted_lookup(data: dict, dotted: str) -> object:
    value: object = data
    for key in dotted.split("."):
        if not isinstance(value, dict) or key not in value:
            raise KeyError(dotted)
        value = value[key]
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--catalog", type=Path)
    parser.add_argument("--biological-profile", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    ledger = root / "state/cluster_decision_ledger.tsv"
    catalog_path = (
        args.catalog
        or Path(__file__).resolve().parents[1]
        / "references/profiles/sheep_ovary_candidate_lineage_catalog.json"
    ).resolve()
    audit = json.loads(args.audit.read_text(encoding="utf-8"))
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    biological_profile_path = (
        args.biological_profile
        or Path(__file__).resolve().parents[1] / "references/profiles/sheep_ovary.json"
    ).resolve()
    biological_profile = json.loads(biological_profile_path.read_text(encoding="utf-8"))
    errors: list[str] = []
    audit_catalog = Path(str(audit.get("candidate_catalog", ""))) if audit.get("candidate_catalog") else None
    if (
        audit_catalog is None
        or not audit_catalog.is_file()
        or audit_catalog.resolve() != catalog_path
        or audit.get("candidate_catalog_sha256") != sha256(catalog_path)
    ):
        errors.append("open-world audit is missing the current candidate catalog binding")
    if catalog.get("taxonomy_policy", {}).get("presence_is_never_required") is not True:
        errors.append("candidate catalog must explicitly state that biological presence is never required")
    boundaries = catalog.get("candidate_boundaries", [])
    if not isinstance(boundaries, list):
        boundaries = []
        errors.append("candidate catalog boundaries must be a list")
    required_candidates = {
        str(item.get("candidate_id", ""))
        for item in boundaries
        if isinstance(item, dict) and item.get("review_required") is True
    }
    if not required_candidates:
        errors.append("candidate catalog has no required boundary reviews")
    if len(required_candidates) != len(
        [item for item in boundaries if isinstance(item, dict) and item.get("review_required") is True]
    ):
        errors.append("candidate catalog contains duplicate or empty required candidate IDs")
    for item in boundaries:
        if not isinstance(item, dict):
            continue
        pointer = str(item.get("profile_program", ""))
        try:
            dotted_lookup(biological_profile, pointer)
        except KeyError:
            errors.append(f"candidate profile_program is unresolved: {item.get('candidate_id')} -> {pointer}")
    if not ledger.is_file():
        errors.append("cluster decision ledger is missing")
    elif audit.get("cluster_ledger_sha256") != sha256(ledger):
        errors.append("open-world lineage audit is stale for the cluster decision ledger")
    if audit.get("taxonomy_mode") != "open_world_literature_informed":
        errors.append("taxonomy_mode must be open_world_literature_informed")
    reviews = audit.get("candidate_reviews", [])
    if not isinstance(reviews, list):
        reviews = []
    by_id = {str(item.get("candidate_id", "")): item for item in reviews if isinstance(item, dict)}
    missing = sorted(required_candidates - set(by_id))
    if missing:
        errors.append("mandatory candidate-boundary reviews are missing: " + ", ".join(missing))
    unknown = sorted(set(by_id) - required_candidates)
    if unknown:
        errors.append(
            "non-catalog candidates belong in additional_candidates, not candidate_reviews: "
            + ", ".join(unknown)
        )
    for candidate_id, item in by_id.items():
        if item.get("outcome") not in OUTCOMES:
            errors.append(f"invalid outcome for {candidate_id}")
        if len(str(item.get("evidence_summary", "")).strip()) < 5:
            errors.append(f"candidate review lacks evidence summary: {candidate_id}")
        if len(str(item.get("action", "")).strip()) < 3:
            errors.append(f"candidate review lacks action: {candidate_id}")
        if item.get("outcome") == "not_evaluable" and not str(item.get("limitation", "")).strip():
            errors.append(f"not_evaluable candidate lacks limitation: {candidate_id}")
    for section in ("additional_candidates", "unexplained_programs"):
        entries = audit.get(section, [])
        if not isinstance(entries, list):
            errors.append(f"{section} must be a list")
            continue
        for index, item in enumerate(entries):
            if not isinstance(item, dict) or not str(item.get("action", "")).strip():
                errors.append(f"{section}[{index}] lacks an explicit action")
    result = {
        "status": "PASS" if not errors else "FAIL",
        "taxonomy_mode": audit.get("taxonomy_mode"),
        "cluster_ledger": str(ledger.resolve()),
        "cluster_ledger_sha256": sha256(ledger) if ledger.is_file() else None,
        "audit": str(args.audit.resolve()),
        "audit_sha256": sha256(args.audit),
        "candidate_catalog": str(catalog_path),
        "candidate_catalog_sha256": sha256(catalog_path),
        "biological_profile": str(biological_profile_path),
        "biological_profile_sha256": sha256(biological_profile_path),
        "catalog_id": catalog.get("catalog_id"),
        "mandatory_candidates": sorted(required_candidates),
        "reviewed_candidates": sorted(by_id),
        "reviewed_families": sorted(
            {
                str(item.get("family", ""))
                for item in boundaries
                if isinstance(item, dict) and item.get("candidate_id") in by_id
            }
        ),
        "additional_candidate_n": len(audit.get("additional_candidates", [])) if isinstance(audit.get("additional_candidates", []), list) else 0,
        "unexplained_program_n": len(audit.get("unexplained_programs", [])) if isinstance(audit.get("unexplained_programs", []), list) else 0,
        "errors": errors,
    }
    out = args.out or root / "provenance/open_world_lineage_audit_validation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
