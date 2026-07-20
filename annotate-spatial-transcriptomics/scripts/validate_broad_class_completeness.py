#!/usr/bin/env python3
"""Fail-closed query-derived audit of present and zero-census broad lineages."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


TRUE = {"1", "true", "yes", "pass", "passed", "reviewed"}
PRESENT_REQUIRED = ("query_marker_program_review", "deg_review", "spatial_morphology_review", "observation_level_review", "large_label_embedding_review")
ABSENT_REQUIRED = ("query_marker_program_review", "spatial_morphology_review", "observation_level_review", "selected_plus_two_review", "large_label_embedding_review", "qc_ood_review", "technical_missingness_review")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split(value: str) -> list[str]:
    return [x for x in re.split(r"[;,\s]+", value or "") if x]


def active(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    superseded = {x for row in rows for x in split(row.get("supersedes", ""))}
    return [row for row in rows if row.get("review_id", "") not in superseded]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(root: Path, catalog_path: Path) -> dict[str, object]:
    root = root.resolve()
    errors: list[str] = []
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    candidates = {str(x.get("candidate_id", "")): x for x in catalog.get("candidate_boundaries", []) if x.get("review_required")}
    required_zero = {cid for cid, row in candidates.items() if row.get("release_level") == "default_broad_candidate"}
    rows = active(read_tsv(root / "state/broad_class_completeness_registry.tsv"))
    by_candidate: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_candidate.setdefault(row.get("candidate_lineage", ""), []).append(row)
    duplicates = [cid for cid, items in by_candidate.items() if len(items) != 1]
    if duplicates:
        errors.append("duplicate active completeness reviews: " + ", ".join(sorted(duplicates)))
    reviewed: list[str] = []
    zero_census: list[str] = []
    for cid, items in by_candidate.items():
        if len(items) != 1:
            continue
        row = items[0]
        try:
            n = int(float(row.get("final_n_observations", "0") or 0))
        except ValueError:
            errors.append(f"{cid}: invalid final_n_observations")
            continue
        if cid not in candidates:
            errors.append(f"{cid}: candidate is not present in the bound catalog")
        if row.get("status") != "audited":
            errors.append(f"{cid}: completeness review is not audited")
        if not row.get("closure_rationale"):
            errors.append(f"{cid}: closure rationale is missing")
        artifact = Path(row.get("evidence_artifact", ""))
        artifact = artifact if artifact.is_absolute() else root / artifact
        if not artifact.is_file() or sha256(artifact) != row.get("evidence_artifact_sha256"):
            errors.append(f"{cid}: evidence artifact is missing or stale")
        if n > 0:
            reviewed.append(cid)
            if row.get("census_status") != "present" or row.get("decision") not in {"supported", "revised", "downgraded"}:
                errors.append(f"{cid}: present lineage lacks a supported/revised/downgraded decision")
            for field in PRESENT_REQUIRED:
                if row.get(field, "").lower() not in TRUE:
                    errors.append(f"{cid}: present-lineage audit lacks {field}")
        else:
            zero_census.append(cid)
            if cid in required_zero:
                if row.get("census_status") != "absent" or row.get("decision") != "refuted_by_multichannel_query_evidence":
                    errors.append(f"{cid}: zero-census default lineage lacks multichannel query-derived refutation")
                for field in ABSENT_REQUIRED:
                    if row.get(field, "").lower() not in TRUE:
                        errors.append(f"{cid}: zero-census audit lacks {field}")
                if row.get("atlas_role") not in {"challenger_only", "not_used", "discordant_low_confidence"}:
                    errors.append(f"{cid}: Atlas cannot establish absence")
    missing_zero = sorted(cid for cid in required_zero if cid not in by_candidate)
    if missing_zero:
        errors.append("default broad candidates lack a completeness review: " + ", ".join(missing_zero))
    return {"status": "PASS" if not errors else "BLOCKED", "catalog_id": catalog.get("catalog_id"),
            "active_reviews": len(rows), "present_reviewed": sorted(reviewed), "zero_census_reviewed": sorted(zero_census),
            "required_default_candidates": sorted(required_zero), "errors": errors}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = validate(args.project_root, args.catalog.resolve())
    out = args.out or args.project_root / "provenance/broad_class_completeness_validation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
