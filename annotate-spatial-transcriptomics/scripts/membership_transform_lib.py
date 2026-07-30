#!/usr/bin/env python3
"""Deterministic validation for ordered annotation-membership transforms."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable

from evidence_schema_lib import sha256
from lineage_controller_lib import (
    deterministic_cell_id_set_hash,
    deterministic_membership_hash,
    read_tsv,
)


TRANSFORM_OPERATIONS = {
    "atlas_unlabeled_broad_rescue",
    "post_merge_unresolved_return",
    "follicle_roi_assignment",
    "follicle_roi_depublication",
    "follicle_roi_reconciliation",
    "source_unit_sync",
    "cell_type_review_patch",
    "final_release_materialization",
}

IDENTITY_NEUTRAL_OPERATIONS = {"source_unit_sync"}

IDENTITY_FIELDS = (
    "final_state",
    "final_broad_label",
    "final_fine_label",
    "final_cell_type",
    "candidate_id",
    "state_annotations",
)

SEMANTIC_FIELDS = (
    "source_boundary",
    "source_cluster",
    *IDENTITY_FIELDS,
    "confidence",
    "assignment_origin",
    "qc_reason",
    "unresolved_reason",
)


def artifact(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise ValueError(f"membership transform artifact is missing: {path}")
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "n_bytes": path.stat().st_size,
    }


def rows_by_cell(path: Path) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    rows = read_tsv(path)
    by_cell = {str(row.get("cell_id", "")): row for row in rows}
    if not rows or "" in by_cell or len(by_cell) != len(rows):
        raise ValueError(f"membership must contain unique nonempty cell_id: {path}")
    return rows, by_cell


def semantic_row(row: dict[str, object]) -> tuple[str, ...]:
    return tuple(str(row.get(field, "")) for field in SEMANTIC_FIELDS)


def identity_row(row: dict[str, object]) -> tuple[str, ...]:
    return tuple(str(row.get(field, "")) for field in IDENTITY_FIELDS)


def changed_rows(
    source: dict[str, dict[str, str]], result: dict[str, dict[str, str]]
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for cell_id in sorted(source):
        old = source[cell_id]
        new = result[cell_id]
        if semantic_row(old) == semantic_row(new):
            continue
        rows.append({
            "cell_id": cell_id,
            "old_state": str(old.get("final_state", "")),
            "new_state": str(new.get("final_state", "")),
            "old_broad_label": str(old.get("final_broad_label", "")),
            "new_broad_label": str(new.get("final_broad_label", "")),
            "old_fine_label": str(old.get("final_fine_label", "")),
            "new_fine_label": str(new.get("final_fine_label", "")),
            "old_final_cell_type": str(old.get("final_cell_type", "")),
            "new_final_cell_type": str(new.get("final_cell_type", "")),
            "old_candidate_id": str(old.get("candidate_id", "")),
            "candidate_id": str(new.get("candidate_id", "")),
            "old_assignment_origin": str(old.get("assignment_origin", "")),
            "new_assignment_origin": str(new.get("assignment_origin", "")),
        })
    return rows


def rows_hash(rows: Iterable[dict[str, object]], fields: tuple[str, ...]) -> str:
    digest = hashlib.sha256()
    for row in sorted(
        rows, key=lambda value: tuple(str(value.get(field, "")) for field in fields)
    ):
        digest.update(
            ("\t".join(str(row.get(field, "")) for field in fields) + "\n").encode(
                "utf-8"
            )
        )
    return digest.hexdigest()


def membership_record(path: Path) -> dict[str, object]:
    rows, _ = rows_by_cell(path)
    return {
        **artifact(path),
        "semantic_sha256": deterministic_membership_hash(rows),
        "cell_id_set_sha256": deterministic_cell_id_set_hash(rows),
        "n_observations": len(rows),
    }


def verify_artifact(record: dict, label: str) -> Path:
    path = Path(str(record.get("path", "")))
    if not path.is_file() or record.get("sha256") != sha256(path):
        raise ValueError(f"{label} is missing or stale")
    return path


def verify_membership_record(record: dict, label: str) -> Path:
    path = verify_artifact(record, label)
    observed = membership_record(path)
    for key in (
        "semantic_sha256", "cell_id_set_sha256", "n_observations",
    ):
        if record.get(key) != observed.get(key):
            raise ValueError(f"{label} {key} is stale")
    return path


def validate_operation(
    operation: str,
    source: dict[str, dict[str, str]],
    result: dict[str, dict[str, str]],
    changes: list[dict[str, str]],
    identity_neutral: bool,
    target_cell_type: str,
) -> None:
    if operation not in TRANSFORM_OPERATIONS:
        raise ValueError(f"unsupported membership transform operation: {operation}")
    if identity_neutral != (operation in IDENTITY_NEUTRAL_OPERATIONS):
        raise ValueError("identity-neutral flag disagrees with transform operation")
    if identity_neutral:
        for cell_id in source:
            if identity_row(source[cell_id]) != identity_row(result[cell_id]):
                raise ValueError("identity-neutral transform changed a biological identity")
    if operation == "atlas_unlabeled_broad_rescue":
        for row in changes:
            if row["old_broad_label"] or not row["new_broad_label"]:
                raise ValueError("Atlas transform may only rescue previously unlabeled cells")
    elif operation == "post_merge_unresolved_return":
        for row in changes:
            if row["old_broad_label"] or not row["new_broad_label"]:
                raise ValueError(
                    "post-merge unresolved review may only return previously unlabeled cells"
                )
    elif operation == "follicle_roi_depublication":
        for row in changes:
            if not row["old_broad_label"] or row["new_broad_label"]:
                raise ValueError("depublication must remove, not replace, a broad label")
            if row["new_state"] != "unresolved_biological":
                raise ValueError("depublication must return the observation to unresolved_biological")
    elif operation == "follicle_roi_reconciliation":
        for row in changes:
            old_broad = row["old_broad_label"]
            new_broad = row["new_broad_label"]
            if old_broad == new_broad:
                raise ValueError(
                    "follicle ROI reconciliation cannot record an identity-neutral change"
                )
            if not old_broad and not new_broad:
                raise ValueError(
                    "follicle ROI reconciliation must assign, withdraw or revise a broad identity"
                )
            if old_broad and not new_broad and row["new_state"] != "unresolved_biological":
                raise ValueError(
                    "follicle ROI withdrawal must return the observation to unresolved_biological"
                )
    elif operation == "cell_type_review_patch":
        if not target_cell_type:
            raise ValueError("cell-type review patch lacks its active target")
        for row in changes:
            if target_cell_type not in {
                row["old_broad_label"], row["new_broad_label"]
            }:
                raise ValueError("cell-type review patch changed an unrelated lineage")
    elif operation == "final_release_materialization":
        for row in changes:
            if row["old_broad_label"] != row["new_broad_label"]:
                raise ValueError("final release materialization changed a frozen broad identity")
            if not row["new_broad_label"]:
                if (
                    "qc" not in row["new_state"].lower()
                    or row["new_final_cell_type"] != "QC/Unknown"
                ):
                    raise ValueError(
                        "final release did not convert an unlabeled member to typed QC/Unknown"
                    )
            elif not row["new_final_cell_type"]:
                raise ValueError("final release omitted final_cell_type for a labeled member")


def validate_evidence_manifest(
    operation: str,
    evidence_path: Path,
    result_path: Path,
    target_cell_type: str,
) -> None:
    try:
        document = json.loads(evidence_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("membership transform evidence is not a JSON manifest") from exc
    stage = str(document.get("stage", ""))
    status = str(document.get("status", ""))
    expected: dict[str, tuple[set[str], set[str]]] = {
        "atlas_unlabeled_broad_rescue": (
            {"atlas_and_completeness_review"}, {"PASS"}
        ),
        "post_merge_unresolved_return": (
            {"post_merge_unresolved_component_review"}, {"PASS"}
        ),
        "follicle_roi_assignment": (
            {"follicle_roi_repair_apply"},
            {"PENDING_POST_REPAIR_BIOLOGICAL_REVIEW"},
        ),
        "follicle_roi_depublication": (
            {"follicle_roi_repair_apply"},
            {"PENDING_POST_REPAIR_BIOLOGICAL_REVIEW"},
        ),
        "follicle_roi_reconciliation": (
            {"follicle_roi_repair_apply"},
            {"PENDING_POST_REPAIR_BIOLOGICAL_REVIEW"},
        ),
        "cell_type_review_patch": (
            {"catalog_wide_lineage_review_apply"},
            {"PASS_REQUIRES_NEXT_REVIEW_ROUND"},
        ),
        "final_release_materialization": (
            {"materialize_final_release"},
            {"PASS", "PENDING_USER_REVIEW_HIGH_UNRESOLVED"},
        ),
    }
    if operation in expected:
        allowed_stages, allowed_statuses = expected[operation]
        if stage not in allowed_stages or status not in allowed_statuses:
            raise ValueError(
                f"{operation} is bound to a non-canonical evidence manifest"
            )
    if operation == "cell_type_review_patch":
        if (
            document.get("formal_batch_closure_performed") is not False
            or str(document.get("active_cell_type", "")) != target_cell_type
        ):
            raise ValueError("cell-type review evidence does not bind one active target")
    result_records = [
        document.get("membership", {}),
        document.get("repaired_membership", {}),
    ]
    bound_results = {
        Path(str(record.get("path", ""))).resolve()
        for record in result_records
        if isinstance(record, dict) and str(record.get("path", ""))
    }
    if operation != "source_unit_sync":
        if result_path.resolve() not in bound_results:
            raise ValueError("membership transform evidence does not bind its result")


def validate_entry(entry: dict, expected_source: Path | None = None) -> Path:
    source_path = verify_membership_record(entry.get("source_membership", {}), "transform source")
    result_path = verify_membership_record(entry.get("result_membership", {}), "transform result")
    if expected_source is not None and source_path.resolve() != expected_source.resolve():
        raise ValueError("membership transform chain is not contiguous")
    source_rows, source = rows_by_cell(source_path)
    result_rows, result = rows_by_cell(result_path)
    if set(source) != set(result):
        raise ValueError("membership transform changed the observation universe")
    delta_path = verify_artifact(entry.get("delta", {}), "transform delta")
    delta = read_tsv(delta_path)
    expected = changed_rows(source, result)
    expected_ids = {row["cell_id"] for row in expected}
    delta_ids = [str(row.get("cell_id", "")) for row in delta]
    if "" in delta_ids or len(delta_ids) != len(set(delta_ids)):
        raise ValueError("membership transform delta has empty or duplicate cell_id")
    if set(delta_ids) != expected_ids:
        raise ValueError("membership transform delta does not equal the semantic change set")
    expected_by_cell = {row["cell_id"]: row for row in expected}
    for row in delta:
        canonical = expected_by_cell[str(row["cell_id"])]
        for field in (
            "old_state", "new_state", "old_broad_label", "new_broad_label",
            "old_fine_label", "new_fine_label", "old_candidate_id",
            "old_final_cell_type", "new_final_cell_type", "candidate_id",
            "old_assignment_origin", "new_assignment_origin",
        ):
            if field in row and str(row.get(field, "")) != canonical[field]:
                raise ValueError(f"membership transform delta disagrees on {field}")
    delta_fields = (
        "cell_id", "old_state", "new_state", "old_broad_label",
        "new_broad_label", "old_fine_label", "new_fine_label",
        "old_final_cell_type", "new_final_cell_type", "old_candidate_id",
        "candidate_id", "old_assignment_origin", "new_assignment_origin",
    )
    if entry.get("delta_semantic_sha256") != rows_hash(expected, delta_fields):
        raise ValueError("membership transform delta semantic hash is stale")
    if int(entry.get("changed_observation_n", -1)) != len(expected):
        raise ValueError("membership transform changed-observation count is stale")
    evidence = entry.get("evidence_manifest", {})
    evidence_path = verify_artifact(evidence, "transform evidence manifest")
    validate_operation(
        str(entry.get("operation", "")), source, result, expected,
        bool(entry.get("identity_neutral", False)),
        str(entry.get("target_cell_type", "")),
    )
    validate_evidence_manifest(
        str(entry.get("operation", "")), evidence_path, result_path,
        str(entry.get("target_cell_type", "")),
    )
    return result_path


def validate_chain_document(document: dict, expected_current: Path | None = None) -> Path:
    if (
        document.get("schema_version") != "1.0"
        or document.get("artifact_role") != "membership_transform_chain"
    ):
        raise ValueError("membership transform chain header is invalid")
    initial = verify_membership_record(document.get("initial_membership", {}), "initial membership")
    current = initial
    entries = document.get("transforms", [])
    if not isinstance(entries, list):
        raise ValueError("membership transform chain entries are malformed")
    for index, entry in enumerate(entries, 1):
        if str(entry.get("transform_id", "")) != f"T{index:04d}":
            raise ValueError("membership transform IDs are not deterministic and contiguous")
        current = validate_entry(entry, current)
    bound_current = verify_membership_record(document.get("current_membership", {}), "current membership")
    if current.resolve() != bound_current.resolve():
        raise ValueError("membership transform chain current membership is stale")
    if int(document.get("transform_n", -1)) != len(entries):
        raise ValueError("membership transform count is stale")
    if expected_current is not None and current.resolve() != expected_current.resolve():
        raise ValueError("membership transform chain does not end at the audited membership")
    return current


def load_and_validate_chain(path: Path, expected_current: Path | None = None) -> dict:
    document = json.loads(path.read_text(encoding="utf-8"))
    validate_chain_document(document, expected_current)
    return document
