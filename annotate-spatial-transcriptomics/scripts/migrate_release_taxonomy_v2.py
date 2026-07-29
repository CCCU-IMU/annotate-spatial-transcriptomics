#!/usr/bin/env python3
"""Non-destructively migrate legacy releases to the v2.2 vascular taxonomy.

Legacy vascular labels are hypotheses, never truth.  The migration preserves
their provenance and makes them unresolved unless an independently validated,
current-query re-adjudication membership is supplied.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

from evidence_schema_lib import sha256
from lineage_controller_lib import catalog_candidates


LEGACY_VASCULAR_LABELS = {
    "Vascular-associated", "Vascular/endothelial", "Vascular/perivascular",
}


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", required=True, type=Path)
    ap.add_argument("--profile", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--broad-column", default="final_broad_label")
    ap.add_argument("--fine-column", default="final_fine_label")
    ap.add_argument("--state-column", default="final_state")
    ap.add_argument("--validated-membership", type=Path)
    ap.add_argument("--validated-membership-manifest", type=Path)
    args = ap.parse_args()

    profile = json.loads(args.profile.read_text(encoding="utf-8"))
    catalog_payload = json.loads(args.catalog.read_text(encoding="utf-8"))
    candidates = catalog_candidates(catalog_payload)
    allowed_broad = {
        str(value.get("release_broad_label", ""))
        for value in candidates.values()
        if value.get("candidate_role") == "broad"
        and value.get("release_broad_label")
    }
    fine_parent = {
        str(value.get("release_fine_label", "")): str(
            value.get("parent_broad_label", "")
            or value.get("release_broad_label", "")
        )
        for value in candidates.values()
        if value.get("candidate_role") == "fine"
        and value.get("release_fine_label")
    }
    forbidden = set(
        catalog_payload.get("taxonomy_policy", {}).get(
            "forbidden_runtime_release_labels", []
        )
    )
    if "Vascular-associated" not in forbidden:
        raise SystemExit("catalog does not forbid legacy Vascular-associated")

    fields, rows = read_rows(args.ledger)
    missing = {
        "cell_id", args.broad_column, args.fine_column, args.state_column,
    } - set(fields)
    if missing:
        raise SystemExit("ledger lacks columns: " + ", ".join(sorted(missing)))
    if any(not row.get("cell_id") for row in rows) or len({row["cell_id"] for row in rows}) != len(rows):
        raise SystemExit("ledger must contain unique nonempty cell_id")

    validated: dict[str, dict[str, str]] = {}
    if bool(args.validated_membership) != bool(args.validated_membership_manifest):
        raise SystemExit("validated membership and manifest must be supplied together")
    if args.validated_membership:
        manifest = json.loads(
            args.validated_membership_manifest.read_text(encoding="utf-8")
        )
        record = manifest.get("membership", {})
        if (
            manifest.get("status") != "PASS"
            or manifest.get("artifact_role")
            != "current_query_taxonomy_readjudication_membership"
            or Path(str(record.get("path", ""))).resolve()
            != args.validated_membership.resolve()
            or record.get("sha256") != sha256(args.validated_membership)
        ):
            raise SystemExit("validated membership is not bound to a PASS re-adjudication")
        validated_fields, validated_rows = read_rows(args.validated_membership)
        required = {"cell_id", "final_broad_label", "final_fine_label"}
        if not required <= set(validated_fields):
            raise SystemExit("validated membership lacks final broad/fine identity")
        for row in validated_rows:
            cell = row.get("cell_id", "")
            broad = row.get("final_broad_label", "")
            fine = row.get("final_fine_label", "")
            if (
                not cell or cell in validated or broad not in allowed_broad
                or broad in forbidden or (fine and fine_parent.get(fine) != broad)
            ):
                raise SystemExit("validated membership contains an invalid identity")
            validated[cell] = row

    extra_fields = [
        "legacy_broad_label", "legacy_fine_label", "taxonomy_migration_status",
        "unresolved_reason", "final_cell_type",
    ]
    for field in extra_fields:
        if field not in fields:
            fields.append(field)
    aliases = profile.get("release_taxonomy", {}).get("broad_aliases", {})
    changes = Counter()
    legacy_ids: set[str] = set()
    for row in rows:
        cell = row["cell_id"]
        broad = row.get(args.broad_column, "")
        fine = row.get(args.fine_column, "")
        if broad in LEGACY_VASCULAR_LABELS:
            legacy_ids.add(cell)
            row["legacy_broad_label"] = broad
            row["legacy_fine_label"] = fine
            row[args.broad_column] = ""
            row[args.fine_column] = ""
            row[args.state_column] = "unresolved_biological"
            row["unresolved_reason"] = (
                "legacy_vascular_taxonomy_requires_current_query_source_subcluster_readjudication"
            )
            row["taxonomy_migration_status"] = "requires_source_subcluster_readjudication"
            row["final_cell_type"] = "QC/Unknown"
            changes["legacy vascular identity -> unresolved_biological"] += 1
            continue
        if broad in aliases:
            target = str(aliases[broad])
            if target in forbidden or target not in allowed_broad:
                raise SystemExit(f"profile alias targets forbidden/unknown label: {target}")
            row[args.broad_column] = target
            broad = target
            changes[f"{row.get('legacy_broad_label') or broad} -> {target}"] += 1
        row["taxonomy_migration_status"] = "unchanged_nonlegacy_identity"
        row["final_cell_type"] = fine or broad or "QC/Unknown"

    unknown_validated = set(validated) - legacy_ids
    if unknown_validated:
        raise SystemExit("validated membership contains cells outside legacy vascular input")
    by_cell = {row["cell_id"]: row for row in rows}
    for cell, patch in validated.items():
        row = by_cell[cell]
        broad = patch["final_broad_label"]
        fine = patch.get("final_fine_label", "")
        row[args.broad_column] = broad
        row[args.fine_column] = fine
        row[args.state_column] = "defined_fine" if fine else "defined_broad_only"
        row["unresolved_reason"] = ""
        row["taxonomy_migration_status"] = "current_query_readjudication_applied"
        row["final_cell_type"] = fine or broad
        changes[f"validated legacy vascular -> {fine or broad}"] += 1

    if any(row.get(args.broad_column, "") in forbidden for row in rows):
        raise SystemExit("migration retained a forbidden runtime label")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.out, "wt") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "schema_version": "2.2",
        "status": "MIGRATED_REQUIRES_FULL_V2_REVALIDATION",
        "source": {"path": str(args.ledger.resolve()), "sha256": sha256(args.ledger)},
        "output": {"path": str(args.out.resolve()), "sha256": sha256(args.out)},
        "profile": {"path": str(args.profile.resolve()), "sha256": sha256(args.profile)},
        "catalog": {"path": str(args.catalog.resolve()), "sha256": sha256(args.catalog)},
        "validated_membership": (
            {"path": str(args.validated_membership.resolve()), "sha256": sha256(args.validated_membership)}
            if args.validated_membership else None
        ),
        "changes": dict(changes),
        "legacy_vascular_n": len(legacy_ids),
        "legacy_vascular_readjudicated_n": len(validated),
        "biological_evidence_reinterpreted": False,
        "legacy_labels_used_as_truth": False,
        "full_catalog_and_biological_review_required": True,
    }
    manifest_path = args.out.with_name(args.out.name + ".taxonomy_migration.json")
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
