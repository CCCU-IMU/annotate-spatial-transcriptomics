#!/usr/bin/env python3
"""Build a complete, catalog-bound marker manifest for broad cell-type review."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from evidence_schema_lib import sha256
from lineage_controller_lib import (
    candidate_can_support_broad_review,
    catalog_candidates,
    write_tsv,
)


GENE = re.compile(r"^(?:[A-Z][A-Z0-9.-]{1,30}|LOC[0-9]+)$")
ANTI_WORDS = ("anti", "contradict", "exclusion", "exclude")
SKIP_WORDS = ("subtype", "state", "spatial", "safety", "policy", "rule", "action")


def resolve(root: dict, dotted: str) -> object:
    value: object = root
    for key in dotted.split("."):
        if not isinstance(value, dict) or key not in value:
            raise SystemExit(f"profile program path does not exist: {dotted}")
        value = value[key]
    return value


def genes(value: object) -> list[str]:
    if not isinstance(value, list):
        return []
    return sorted({str(item) for item in value if GENE.fullmatch(str(item))})


def collect(program: object) -> list[tuple[str, str, str]]:
    if not isinstance(program, dict):
        return []
    rows: list[tuple[str, str, str]] = []
    positive = program.get("positive_families", {})
    if isinstance(positive, dict):
        for family, values in positive.items():
            for gene in genes(values):
                rows.append(("positive_family", str(family), gene))
    for key, value in program.items():
        lowered = str(key).lower()
        if key == "positive_families" or any(word in lowered for word in SKIP_WORDS):
            continue
        role = "anti_family" if any(word in lowered for word in ANTI_WORDS) else "identity_support"
        for gene in genes(value):
            rows.append((role, str(key), gene))
        # Some profiles group broad identity programs one level below the
        # lineage object.  Preserve those lists, but never import subtype/state
        # modules into broad reconstruction evidence.
        if isinstance(value, dict):
            for child, child_value in value.items():
                child_lower = str(child).lower()
                if any(word in child_lower for word in SKIP_WORDS):
                    continue
                child_role = (
                    "anti_family"
                    if any(word in f"{lowered}.{child_lower}" for word in ANTI_WORDS)
                    else "identity_support"
                )
                for gene in genes(child_value):
                    rows.append((child_role, f"{key}.{child}", gene))
    return sorted(set(rows))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    profile = json.loads(args.profile.read_text(encoding="utf-8"))
    candidates = catalog_candidates(json.loads(args.catalog.read_text(encoding="utf-8")))
    rows: list[dict[str, object]] = []
    for candidate_id, candidate in sorted(candidates.items()):
        if not candidate_can_support_broad_review(candidate):
            continue
        path = str(candidate.get("profile_program", "")).strip()
        if not path:
            raise SystemExit(f"broad-review candidate lacks profile_program: {candidate_id}")
        extracted = collect(resolve(profile, path))
        if not any(role == "positive_family" for role, _, _ in extracted):
            raise SystemExit(f"broad-review candidate lacks positive marker families: {candidate_id}")
        for role, family, gene in extracted:
            rows.append({
                "candidate_id": candidate_id,
                "broad_label": str(candidate.get("release_broad_label", "")),
                "candidate_role": str(candidate.get("candidate_role", "")),
                "profile_program": path,
                "evidence_role": role,
                "family_id": family,
                "gene": gene,
            })
    rows = sorted(
        {tuple(row.items()): row for row in rows}.values(),
        key=lambda row: (
            str(row["broad_label"]), str(row["candidate_id"]),
            str(row["evidence_role"]), str(row["family_id"]), str(row["gene"]),
        ),
    )
    args.out.mkdir(parents=True, exist_ok=True)
    manifest_path = args.out / "cell_type_review_marker_manifest.tsv"
    write_tsv(manifest_path, rows)
    summary = {
        "schema_version": "2.2",
        "artifact_role": "cell_type_review_marker_manifest",
        "profile": {"path": str(args.profile.resolve()), "sha256": sha256(args.profile)},
        "catalog": {"path": str(args.catalog.resolve()), "sha256": sha256(args.catalog)},
        "manifest": {"path": str(manifest_path.resolve()), "sha256": sha256(manifest_path)},
        "candidate_n": len({str(row["candidate_id"]) for row in rows}),
        "broad_label_n": len({str(row["broad_label"]) for row in rows}),
        "gene_n": len({str(row["gene"]) for row in rows}),
        "row_n": len(rows),
    }
    summary_path = args.out / "cell_type_review_marker_manifest.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
