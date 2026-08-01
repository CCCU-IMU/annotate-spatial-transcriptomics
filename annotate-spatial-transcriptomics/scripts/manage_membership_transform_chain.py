#!/usr/bin/env python3
"""Create, append to or validate the canonical membership-transform chain."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from lineage_controller_lib import write_tsv
from evidence_schema_lib import validate_json_against_schema
from membership_transform_lib import (
    IDENTITY_NEUTRAL_OPERATIONS,
    TRANSFORM_OPERATIONS,
    artifact,
    changed_rows,
    load_and_validate_chain,
    membership_record,
    rows_by_cell,
    rows_hash,
    validate_chain_document,
    validate_evidence_manifest,
    validate_operation,
)


DELTA_FIELDS = [
    "cell_id", "old_state", "new_state", "old_broad_label",
    "new_broad_label", "old_fine_label", "new_fine_label",
    "old_final_cell_type", "new_final_cell_type", "old_candidate_id",
    "candidate_id", "old_assignment_origin", "new_assignment_origin",
]


def write_chain(document: dict, out: Path) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    path = out / "membership_transform_chain.json"
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    schema = (
        Path(__file__).resolve().parents[1]
        / "schemas/membership_transform_chain.schema.json"
    )
    _, errors = validate_json_against_schema(path, schema)
    if errors:
        raise SystemExit(
            "membership transform chain violates schema: "
            + "; ".join(errors)
        )
    return path


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="command", required=True)
    init = sub.add_parser("init")
    init.add_argument("--membership", required=True, type=Path)
    init.add_argument("--out", required=True, type=Path)
    append = sub.add_parser("append")
    append.add_argument("--chain", required=True, type=Path)
    append.add_argument("--operation", required=True, choices=sorted(TRANSFORM_OPERATIONS))
    append.add_argument("--source", required=True, type=Path)
    append.add_argument("--result", required=True, type=Path)
    append.add_argument("--evidence-manifest", required=True, type=Path)
    append.add_argument("--target-cell-type", default="")
    append.add_argument("--out", required=True, type=Path)
    validate = sub.add_parser("validate")
    validate.add_argument("--chain", required=True, type=Path)
    validate.add_argument("--current-membership", type=Path)
    validate.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    if args.command == "init":
        record = membership_record(args.membership.resolve())
        document = {
            "schema_version": "1.0",
            "status": "PASS",
            "artifact_role": "membership_transform_chain",
            "initial_membership": record,
            "current_membership": record,
            "transform_n": 0,
            "transforms": [],
        }
        path = write_chain(document, args.out)
        print(json.dumps({"status": "PASS", "chain": artifact(path)}, indent=2))
        return 0

    if args.command == "append":
        document = load_and_validate_chain(args.chain)
        current = Path(str(document["current_membership"]["path"]))
        if current.resolve() != args.source.resolve():
            raise SystemExit("new membership transform does not continue the active chain")
        source_rows, source = rows_by_cell(args.source.resolve())
        result_rows, result = rows_by_cell(args.result.resolve())
        if set(source) != set(result):
            raise SystemExit("new membership transform changed the cell universe")
        changes = changed_rows(source, result)
        args.out.mkdir(parents=True, exist_ok=True)
        delta_path = args.out / f"membership_transform_T{len(document['transforms']) + 1:04d}.tsv"
        write_tsv(delta_path, changes, fields=DELTA_FIELDS)
        entry = {
            "transform_id": f"T{len(document['transforms']) + 1:04d}",
            "operation": args.operation,
            "identity_neutral": args.operation in IDENTITY_NEUTRAL_OPERATIONS,
            "target_cell_type": args.target_cell_type,
            "source_membership": membership_record(args.source.resolve()),
            "result_membership": membership_record(args.result.resolve()),
            "delta": artifact(delta_path),
            "delta_semantic_sha256": rows_hash(changes, tuple(DELTA_FIELDS)),
            "changed_observation_n": len(changes),
            "evidence_manifest": artifact(args.evidence_manifest.resolve()),
        }
        updated = dict(document)
        updated["transforms"] = [*document["transforms"], entry]
        updated["transform_n"] = len(updated["transforms"])
        updated["current_membership"] = entry["result_membership"]
        # load_and_validate_chain() above has already validated every existing
        # entry.  Revalidating the full chain after appending duplicates all
        # 468k-row membership reads and can exhaust a scheduler process during
        # final materialization.  Validate the one new contiguous entry here;
        # the resulting chain remains subject to full validation on its next
        # load and at release audit.
        # Reuse the source/result dictionaries already loaded above.  Calling
        # validate_entry() here would allocate a second pair of 468k-row
        # dictionaries while the first pair is still live.
        validate_operation(
            args.operation, source, result, changes,
            args.operation in IDENTITY_NEUTRAL_OPERATIONS,
            args.target_cell_type,
        )
        validate_evidence_manifest(
            args.operation, args.evidence_manifest.resolve(),
            args.result.resolve(), args.target_cell_type,
        )
        path = write_chain(updated, args.out)
        print(json.dumps({"status": "PASS", "chain": artifact(path), "transform": entry}, indent=2))
        return 0

    document = load_and_validate_chain(args.chain, args.current_membership)
    args.out.mkdir(parents=True, exist_ok=True)
    result_path = args.out / "membership_transform_chain_validation.json"
    result = {
        "status": "PASS",
        "chain": artifact(args.chain),
        "transform_n": document["transform_n"],
        "current_membership": document["current_membership"],
    }
    result_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
