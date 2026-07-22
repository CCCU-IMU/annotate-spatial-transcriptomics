#!/usr/bin/env python3
"""Materialize canonical v2 broad/fine labels without mutating the source ledger."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

from evidence_schema_lib import sha256


def open_text(path: Path, mode: str):
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(mode, encoding="utf-8", newline="")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", required=True, type=Path)
    ap.add_argument("--profile", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--broad-column", default="final_broad_label")
    ap.add_argument("--fine-column", default="final_fine_label")
    ap.add_argument("--state-column", default="final_state")
    args = ap.parse_args()
    profile = json.loads(args.profile.read_text(encoding="utf-8"))
    taxonomy = profile.get("release_taxonomy", {})
    aliases = taxonomy.get("broad_aliases", {})
    with open_text(args.ledger, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    missing = {args.broad_column, args.fine_column, args.state_column} - set(fields)
    if missing:
        raise SystemExit("ledger lacks columns: " + ", ".join(sorted(missing)))
    changes = Counter()
    for row in rows:
        broad = row.get(args.broad_column, "")
        if broad not in aliases:
            continue
        target = aliases[broad]
        row[args.broad_column] = target
        changes[f"{broad} -> {target}"] += 1
        if broad == "Pericyte/mural" and not row.get(args.fine_column, ""):
            row[args.fine_column] = "Pericyte/mural"
            row[args.state_column] = "defined_fine"
            changes["preserved Pericyte/mural as fine identity"] += 1
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.out, "wt") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    manifest = {
        "schema_version": "2.0",
        "status": "MIGRATED_REQUIRES_FULL_V2_REVALIDATION",
        "source": {"path": str(args.ledger.resolve()), "sha256": sha256(args.ledger)},
        "output": {"path": str(args.out.resolve()), "sha256": sha256(args.out)},
        "profile": {"path": str(args.profile.resolve()), "sha256": sha256(args.profile)},
        "changes": dict(changes),
        "biological_evidence_reinterpreted": False,
    }
    manifest_path = args.out.with_name(args.out.name + ".taxonomy_migration.json")
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
