#!/usr/bin/env python3
"""Apply a bounded membership repair by cell ID, never by row position."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from lineage_controller_lib import (
    deterministic_membership_hash,
    read_tsv,
    sha256,
    write_tsv,
)


def unique_index(rows: list[dict[str, str]], label: str) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        cell_id = str(row.get("cell_id", "")).strip()
        if not cell_id:
            raise SystemExit(f"{label} contains an empty cell_id")
        if cell_id in index:
            raise SystemExit(f"{label} contains duplicate cell_id: {cell_id}")
        index[cell_id] = row
    if not index:
        raise SystemExit(f"{label} is empty")
    return index


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-membership", required=True, type=Path)
    parser.add_argument("--proposal", required=True, type=Path)
    parser.add_argument("--update-column", required=True, action="append")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    base_rows = read_tsv(args.base_membership)
    proposal_rows = read_tsv(args.proposal)
    base = unique_index(base_rows, "base membership")
    proposal = unique_index(proposal_rows, "proposal")
    unknown = sorted(set(proposal).difference(base))
    if unknown:
        raise SystemExit(
            "proposal contains cell IDs outside the base membership: "
            + ", ".join(unknown[:5])
        )
    columns = list(dict.fromkeys(str(value) for value in args.update_column))
    if any(not value or value == "cell_id" for value in columns):
        raise SystemExit("update columns must be nonempty and cannot include cell_id")
    for column in columns:
        if any(column not in row for row in proposal_rows):
            raise SystemExit(f"proposal lacks update column: {column}")

    output_rows: list[dict[str, object]] = []
    for base_row in base_rows:
        cell_id = str(base_row["cell_id"])
        output = dict(base_row)
        patch = proposal.get(cell_id)
        if patch is not None:
            for column in columns:
                output[column] = patch[column]
        for column in columns:
            output.setdefault(column, "")
        output_rows.append(output)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out, output_rows)
    reread = read_tsv(args.out)
    if [row["cell_id"] for row in reread] != [row["cell_id"] for row in base_rows]:
        raise SystemExit("cell-ID patch did not preserve base membership order")
    for row in reread:
        cell_id = row["cell_id"]
        if cell_id in proposal:
            for column in columns:
                if row.get(column, "") != proposal[cell_id].get(column, ""):
                    raise SystemExit(
                        f"cell-ID patch verification failed for {cell_id}:{column}"
                    )
        else:
            original = base[cell_id]
            for column in columns:
                if row.get(column, "") != original.get(column, ""):
                    raise SystemExit(
                        f"non-proposal member changed for {cell_id}:{column}"
                    )

    manifest = {
        "status": "PASS",
        "schema_version": "2.2",
        "operation": "cell_id_membership_patch",
        "base_membership": {
            "path": str(args.base_membership.resolve()),
            "sha256": sha256(args.base_membership),
            "n_observations": len(base_rows),
        },
        "proposal": {
            "path": str(args.proposal.resolve()),
            "sha256": sha256(args.proposal),
            "n_observations": len(proposal_rows),
        },
        "output": {
            "path": str(args.out.resolve()),
            "sha256": sha256(args.out),
            "semantic_sha256": deterministic_membership_hash(output_rows),
        },
        "update_columns": columns,
        "join_key": "cell_id",
        "positional_assignment_used": False,
        "base_order_preserved": True,
        "proposal_ids_subset_of_base": True,
        "nonproposal_values_unchanged": True,
    }
    manifest_path = args.out.parent / f"{args.out.name}.cell_id_patch_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
