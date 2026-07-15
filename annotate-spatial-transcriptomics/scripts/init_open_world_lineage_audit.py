#!/usr/bin/env python3
"""Create a complete, ledger-bound sheep-ovary candidate-boundary audit scaffold."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--catalog", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    root = args.project_root.resolve()
    ledger = root / "state/cluster_decision_ledger.tsv"
    if not ledger.is_file():
        raise SystemExit("cluster decision ledger is required before open-world lineage review")
    catalog_path = (
        args.catalog
        or Path(__file__).resolve().parents[1]
        / "references/profiles/sheep_ovary_candidate_lineage_catalog.json"
    ).resolve()
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    out = args.out or root / "provenance/open_world_lineage_audit.json"
    if out.exists() and not args.force:
        raise SystemExit(f"refusing to overwrite existing audit without --force: {out}")
    reviews = [
        {
            "candidate_id": item["candidate_id"],
            "family": item["family"],
            "release_level": item["release_level"],
            "outcome": "pending",
            "evidence_summary": "",
            "action": "",
            "limitation": "",
        }
        for item in catalog.get("candidate_boundaries", [])
        if isinstance(item, dict) and item.get("review_required") is True
    ]
    audit = {
        "taxonomy_mode": "open_world_literature_informed",
        "catalog_id": catalog.get("catalog_id"),
        "candidate_catalog": str(catalog_path),
        "candidate_catalog_sha256": sha256(catalog_path),
        "cluster_ledger": str(ledger.resolve()),
        "cluster_ledger_sha256": sha256(ledger),
        "candidate_reviews": reviews,
        "additional_candidates": [],
        "unexplained_programs": [],
        "instructions": "Replace every pending outcome using current-query full-feature DEG, anti-marker, stability and spatial evidence. A negative audit is valid; omission and forced presence are not.",
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "SCAFFOLD_CREATED", "audit": str(out), "candidate_n": len(reviews)}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
