#!/usr/bin/env python3
"""Ensure every default broad candidate can satisfy the broad-family gate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def resolve(payload: dict, dotted: str):
    value = payload
    for key in dotted.split("."):
        value = value[key]
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    profile = json.loads(args.profile.read_text(encoding="utf-8"))
    catalog = json.loads(args.catalog.read_text(encoding="utf-8"))
    errors: list[str] = []
    rows: list[dict] = []
    for candidate in catalog.get("candidate_boundaries", []):
        if candidate.get("review_required") is not True:
            continue
        candidate_id = candidate.get("candidate_id", "")
        path = candidate.get("profile_program", "")
        try:
            program = resolve(profile, path)
        except (KeyError, TypeError):
            errors.append(f"{candidate_id}: profile program cannot be resolved: {path}")
            continue
        families = program.get("positive_families", {}) if isinstance(program, dict) else {}
        valid = {
            name: sorted({str(gene).strip() for gene in genes if str(gene).strip()})
            for name, genes in families.items()
            if isinstance(genes, list)
        }
        valid = {name: genes for name, genes in valid.items() if genes}
        overlaps = []
        names = sorted(valid)
        for i, left in enumerate(names):
            for right in names[i + 1 :]:
                shared = sorted(set(valid[left]).intersection(valid[right]))
                if shared:
                    overlaps.append({"left": left, "right": right, "shared": shared})
        if len(valid) < 2:
            errors.append(
                f"{candidate_id}: open-world broad scan requires >=2 explicit positive_families, found {len(valid)}"
            )
        if overlaps:
            errors.append(f"{candidate_id}: positive families overlap and are not independent: {overlaps}")
        rows.append(
            {
                "candidate_id": candidate_id,
                "profile_program": path,
                "positive_family_n": len(valid),
                "positive_families": valid,
                "overlaps": overlaps,
            }
        )
    result = {
        "status": "PASS" if not errors else "BLOCKED",
        "schema_version": "2.0",
        "profile": str(args.profile.resolve()),
        "profile_sha256": sha256(args.profile),
        "catalog": str(args.catalog.resolve()),
        "catalog_sha256": sha256(args.catalog),
        "review_required_candidates_checked": len(rows),
        "default_candidates_checked": sum(
            candidate.get("release_level") in {"default_broad_candidate", "context_specific_broad_candidate"}
            for candidate in catalog.get("candidate_boundaries", [])
            if candidate.get("review_required") is True
        ),
        "candidates": rows,
        "errors": errors,
        "broad_presence_rule": "absolute full-feature detection/prevalence and pseudobulk; centered scores are comparative only",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
