#!/usr/bin/env python3
"""Convert exogenous sample context into evaluation-only candidate evidence."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


LUTEAL_COMPATIBLE = {
    "interestrus",
    "inter-estrus",
    "diestrus",
    "dioestrus",
    "luteal",
    "luteal phase",
    "mid luteal",
    "late luteal",
    "pregnant",
    "pregnancy",
    "发情间期",
    "黄体期",
    "妊娠",
}


def normalize(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip().lower().replace("_", " "))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--reproductive-stage", default="")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    species = normalize(args.species)
    tissue = normalize(args.tissue)
    stage = normalize(args.reproductive_stage)
    sheep_ovary = (
        any(token in species for token in ("sheep", "ovis", "ovine", "羊"))
        and any(token in tissue for token in ("ovary", "ovarian", "卵巢"))
    )
    luteal_stage = stage in {normalize(value) for value in LUTEAL_COMPATIBLE}
    rows = [{
        "candidate_id": "luteal_steroidogenic",
        "status": "supported" if sheep_ovary and luteal_stage else "not_evaluable",
        "context_dimension": "reproductive_stage_compatibility",
        "observed_value": args.reproductive_stage,
        "evidence_scope": "evaluation_permission_only",
        "identity_writeback_authority": "false",
        "reason": (
            "stage permits corpus-luteum candidate evaluation; molecular and "
            "spatial identity evidence remains mandatory"
            if sheep_ovary and luteal_stage
            else "stage does not establish luteal compatibility"
        ),
    }]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "status": "PASS",
        "schema_version": "2.2",
        "species": args.species,
        "tissue": args.tissue,
        "reproductive_stage": args.reproductive_stage,
        "semantics": (
            "exogenous context only; permits evaluation of a gated candidate "
            "but cannot score or assign observations"
        ),
        "historical_labels_read": False,
        "membership_read": False,
        "context_evidence": str(args.out.resolve()),
    }
    manifest_path = args.out.with_suffix(args.out.suffix + ".manifest.json")
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
