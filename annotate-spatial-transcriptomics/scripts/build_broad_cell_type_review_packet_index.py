#!/usr/bin/env python3
"""Bind one complete, machine-readable evidence packet to every broad review."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from evidence_schema_lib import sha256
from lineage_controller_lib import read_tsv, write_tsv


def artifact(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "n_bytes": path.stat().st_size,
    }


def bound_artifact(document: dict, key: str) -> Path:
    record = document.get("artifacts", {}).get(key, {})
    path = Path(str(record.get("path", "")))
    if not path.is_file() or record.get("sha256") != sha256(path):
        raise SystemExit(f"broad review evidence artifact is missing or stale: {key}")
    return path


def packet_hash(payload: dict[str, object]) -> str:
    encoded = json.dumps(
        payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def rows_semantic_hash(
    rows: list[dict[str, str]], fields: tuple[str, ...]
) -> str:
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--review-manifest", required=True, type=Path)
    ap.add_argument("--broad-evidence-manifest", required=True, type=Path)
    ap.add_argument("--pseudobulk-summary-manifest", required=True, type=Path)
    ap.add_argument("--threshold-registry", required=True, type=Path)
    ap.add_argument("--biological-quality-review", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    review = json.loads(args.review_manifest.read_text(encoding="utf-8"))
    broad = json.loads(args.broad_evidence_manifest.read_text(encoding="utf-8"))
    pseudobulk = json.loads(
        args.pseudobulk_summary_manifest.read_text(encoding="utf-8")
    )
    thresholds = json.loads(args.threshold_registry.read_text(encoding="utf-8"))[
        "catalog_wide_lineage_review_policy"
    ]
    minimum_markers = int(thresholds["minimum_evaluable_positive_markers"])
    minimum_families = int(thresholds["minimum_evaluable_positive_families"])
    if (
        review.get("stage") != "post_atlas_catalog_wide_lineage_review"
        or review.get("status") != "ITERATION_REQUIRED"
        or broad.get("artifact_role")
        != "broad_cell_type_targeted_review_evidence"
        or pseudobulk.get("artifact_role")
        != "broad_cell_type_full_transcriptome_pseudobulk_summary"
    ):
        raise SystemExit("broad review packet inputs are not canonical evidence artifacts")
    for manifest, key, expected in (
        (broad, "review_manifest", args.review_manifest),
        (pseudobulk, "broad_evidence_manifest", args.broad_evidence_manifest),
        (broad, "threshold_registry", args.threshold_registry),
    ):
        record = manifest.get(key, {})
        if (
            Path(str(record.get("path", ""))).resolve() != expected.resolve()
            or record.get("sha256") != sha256(expected)
        ):
            raise SystemExit(f"broad review packet evidence chain differs for {key}")

    queue_path = bound_artifact(review, "review_queue")
    controller_artifact_paths = {
        key: bound_artifact(review, key)
        for key in (
            "lineage_review_matrix", "present_label_precision_audit",
            "outside_label_recall_components",
            "outside_label_recall_component_membership",
            "outside_label_group_watch", "broad_lineage_review_summary",
            "broad_lineage_review_scope_membership",
        )
    }
    raw_summary_path = bound_artifact(broad, "summary")
    precision_path = bound_artifact(broad, "precision")
    recall_component_path = bound_artifact(broad, "recall_components")
    recall_membership_path = bound_artifact(broad, "recall_membership")
    current_question_path = bound_artifact(broad, "current_member_questions")
    spatial_index_path = bound_artifact(broad, "spatial_plot_index")
    current_summary_path = bound_artifact(pseudobulk, "current_summary")
    recall_summary_path = bound_artifact(pseudobulk, "summary")
    recall_similarity_path = bound_artifact(pseudobulk, "similarity")
    recall_differential_path = bound_artifact(pseudobulk, "top_differential")
    current_differential_path = bound_artifact(
        pseudobulk, "current_top_differential"
    )
    coordinates_record = broad.get("coordinates", {})
    coordinates_path = Path(str(coordinates_record.get("path", "")))
    if (
        not coordinates_path.is_file()
        or coordinates_record.get("sha256") != sha256(coordinates_path)
    ):
        raise SystemExit("broad review spatial evidence is missing or stale")

    queue = read_tsv(queue_path)
    controller_review_rows = read_tsv(
        controller_artifact_paths["broad_lineage_review_summary"]
    )
    controller_precision_rows = read_tsv(
        controller_artifact_paths["present_label_precision_audit"]
    )
    controller_component_rows = read_tsv(
        controller_artifact_paths["outside_label_recall_components"]
    )
    controller_watch_rows = read_tsv(
        controller_artifact_paths["outside_label_group_watch"]
    )
    raw_rows = read_tsv(raw_summary_path)
    precision_rows = read_tsv(precision_path)
    recall_component_rows = read_tsv(recall_component_path)
    read_tsv(recall_membership_path)
    current_question_rows = read_tsv(current_question_path)
    spatial_rows = read_tsv(spatial_index_path)
    current_rows = read_tsv(current_summary_path)
    recall_rows = read_tsv(recall_summary_path)
    recall_similarity_rows = read_tsv(recall_similarity_path)
    recall_differential_rows = read_tsv(recall_differential_path)
    current_differential_rows = read_tsv(current_differential_path)
    raw_by_review = {str(row.get("review_id", "")): row for row in raw_rows}
    current_by_broad = {
        str(row.get("broad_label", "")): row for row in current_rows
    }
    recall_by_broad: dict[str, list[dict[str, str]]] = {}
    for row in recall_rows:
        recall_by_broad.setdefault(
            str(row.get("target_broad_label", "")), []
        ).append(row)
    controller_review_by_broad = {
        str(row.get("broad_label", "")): row
        for row in controller_review_rows
    }

    def by_broad(rows: list[dict[str, str]], column: str) -> dict[str, list[dict[str, str]]]:
        result: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            result.setdefault(str(row.get(column, "")), []).append(row)
        return result

    controller_precision_by_broad = by_broad(
        controller_precision_rows, "broad_label"
    )
    controller_components_by_broad = by_broad(
        controller_component_rows, "target_broad_label"
    )
    controller_watches_by_broad = by_broad(
        controller_watch_rows, "target_broad_label"
    )
    precision_by_broad: dict[str, list[dict[str, str]]] = {}
    for row in precision_rows:
        precision_by_broad.setdefault(
            str(row.get("broad_label", "")), []
        ).append(row)
    components_by_broad: dict[str, list[dict[str, str]]] = {}
    for row in recall_component_rows:
        components_by_broad.setdefault(
            str(row.get("broad_label", "")), []
        ).append(row)
    current_questions_by_broad: dict[str, list[dict[str, str]]] = {}
    for row in current_question_rows:
        current_questions_by_broad.setdefault(
            str(row.get("broad_label", "")), []
        ).append(row)
    spatial_by_review = {
        str(row.get("review_id", "")): row for row in spatial_rows
    }
    if (
        "" in spatial_by_review
        or len(spatial_by_review) != len(spatial_rows)
        or len(spatial_rows) != len(queue)
    ):
        raise SystemExit("broad review spatial plots do not cover the queue exactly")
    for row in spatial_rows:
        path = Path(str(row.get("path", "")))
        if not path.is_file() or row.get("sha256") != sha256(path):
            raise SystemExit("broad review spatial plot is missing or stale")

    quality = None
    quality_artifact = None
    quality_broad: dict[str, dict[str, str]] = {}
    if args.biological_quality_review:
        quality = json.loads(
            args.biological_quality_review.read_text(encoding="utf-8")
        )
        quality_membership = quality.get("membership", {})
        broad_membership = broad.get("membership", {})
        if (
            quality.get("artifact_role") != "biological_quality_review"
            or Path(str(quality_membership.get("path", ""))).resolve()
            != Path(str(broad_membership.get("path", ""))).resolve()
            or quality_membership.get("sha256") != broad_membership.get("sha256")
        ):
            raise SystemExit(
                "sheep-ovary biological quality review differs from the broad review membership"
            )
        endpoint = quality.get("quality_endpoints", {}).get(
            "spatial_celltype_localization", {}
        )
        broad_review_record = endpoint.get("review", {})
        broad_review_path = Path(str(broad_review_record.get("path", "")))
        if (
            not broad_review_path.is_file()
            or broad_review_record.get("sha256") != sha256(broad_review_path)
        ):
            raise SystemExit("sheep-ovary broad spatial review is missing or stale")
        quality_broad = {
            str(row.get("broad_label", "")): row
            for row in read_tsv(broad_review_path)
        }
        quality_artifact = artifact(args.biological_quality_review)

    errors: list[str] = []
    packet_rows: list[dict[str, object]] = []
    for queued in queue:
        review_id = str(queued.get("review_id", ""))
        broad_label = str(queued.get("target_broad_label", ""))
        raw = raw_by_review.get(review_id)
        current = current_by_broad.get(broad_label)
        recalls = sorted(
            recall_by_broad.get(broad_label, []),
            key=lambda row: str(row.get("recall_group", "")),
        )
        if not raw or str(raw.get("broad_label", "")) != broad_label:
            errors.append(f"{review_id}: raw-count evidence row is absent")
            continue
        if broad_label not in controller_review_by_broad:
            errors.append(
                f"{review_id}: selected-resolution controller review row is absent"
            )
            continue
        if str(raw.get("unit_signature", "")) != str(
            queued.get("unit_signature", "")
        ):
            errors.append(f"{review_id}: raw-count evidence scope differs from queue")
            continue
        current_n = int(float(raw.get("current_n", 0) or 0))
        recall_n = int(float(raw.get("accepted_recall_observation_n", 0) or 0))
        competitor_n = int(float(raw.get("current_competitor_question_n", 0) or 0))
        cross_type_n = int(float(raw.get("cross_type_over_recall_question_n", 0) or 0))
        positive_marker_n_available = int(float(
            raw.get("positive_marker_n_available", 0) or 0
        ))
        positive_family_n_available = int(float(
            raw.get("positive_family_n_available", 0) or 0
        ))
        current_differential_n = sum(
            str(row.get("broad_label", "")) == broad_label
            for row in current_differential_rows
        )
        recall_differential_n = sum(
            str(row.get("target_broad_label", "")) == broad_label
            for row in recall_differential_rows
        )
        recall_similarity_n = sum(
            str(row.get("target_broad_label", "")) == broad_label
            for row in recall_similarity_rows
        )
        precision_evaluable = bool(
            current_n > 0 and current
            and precision_by_broad.get(broad_label)
            and current_differential_n > 0
        )
        recall_evaluable = bool(
            positive_marker_n_available >= minimum_markers
            and positive_family_n_available >= minimum_families
            and (
                recall_n == 0
                or (
                    recalls and components_by_broad.get(broad_label)
                    and recall_differential_n > 0
                    and recall_similarity_n > 0
                )
            )
        )
        molecular_evaluable = bool(
            current_n > 0 and current and current_differential_n > 0
            and int(float((current or {}).get(
                "target_positive_marker_n_available", 0
            ) or 0)) >= minimum_markers
            and int(float((current or {}).get(
                "target_positive_marker_n_detected", 0
            ) or 0)) >= minimum_markers
            and (
                recall_n == 0
                or (
                    recalls and recall_differential_n > 0
                    and recall_similarity_n > 0
                )
            )
        )
        spatial = spatial_by_review.get(review_id, {})
        current_precision_question_n = len(
            current_questions_by_broad.get(broad_label, [])
        )
        spatial_evaluable = bool(
            broad.get("membership_cell_id_semantic_sha256")
            and broad.get("review_queue_n") == len(queue)
            and spatial
        )
        ovary_spatial_status = "not_applicable"
        oocyte_review_status = "not_applicable"
        oocyte_canonical_bound = False
        follicle_histology_status = "not_applicable"
        if quality is not None:
            ovary_spatial_status = str(
                quality_broad.get(broad_label, {}).get(
                    "status", "not_evaluable"
                )
            )
            oocyte_review_status = str(
                quality.get("quality_endpoints", {})
                .get("oocyte_annotation_quality", {})
                .get("status", "not_evaluable")
            )
            canonical_record = quality.get("canonical_oocyte_review") or {}
            canonical_path = Path(str(canonical_record.get("path", "")))
            oocyte_canonical_bound = bool(
                canonical_path.is_file()
                and canonical_record.get("sha256") == sha256(canonical_path)
            )
            follicle_histology_status = str(
                quality.get("quality_endpoints", {})
                .get("follicle_roi_histology", {})
                .get("status", "not_evaluable")
            )
        if (
            broad_label == "Oocyte"
            and oocyte_review_status == "PASS"
            and oocyte_canonical_bound
        ):
            # The canonical query-only Oocyte adjudication is already bound to
            # the exact released member set by the ovary quality validator. It
            # supersedes sparse ordinary marker-score questions for those same
            # members; it does not authorize any outside recall.
            current_precision_question_n = 0
        payload = {
            "review_id": review_id,
            "review_mode": queued.get("review_mode", ""),
            "target_broad_label": broad_label,
            "unit_signature": queued.get("unit_signature", ""),
            "raw_count_summary": raw,
            "current_pseudobulk_summary": current or {},
            "recall_pseudobulk_summaries": recalls,
            "controller_review_summary": controller_review_by_broad.get(
                broad_label, {}
            ),
            "controller_precision_source_groups": (
                controller_precision_by_broad.get(broad_label, [])
            ),
            "controller_recall_components": (
                controller_components_by_broad.get(broad_label, [])
            ),
            "controller_group_watches": (
                controller_watches_by_broad.get(broad_label, [])
            ),
            "precision_source_groups": precision_by_broad.get(broad_label, []),
            "current_member_questions": {
                "raw_n": len(current_questions_by_broad.get(broad_label, [])),
                "effective_unresolved_n": current_precision_question_n,
                "semantic_sha256": rows_semantic_hash(
                    current_questions_by_broad.get(broad_label, []),
                    ("broad_label", "cell_id", "question_reason",
                     "challenger_broad_label"),
                ),
            },
            "recall_components": components_by_broad.get(broad_label, []),
            "spatial_plot": spatial,
            "current_differential_n": current_differential_n,
            "recall_differential_n": recall_differential_n,
            "recall_similarity_n": recall_similarity_n,
            "ovary_spatial_status": ovary_spatial_status,
            "oocyte_review_status": oocyte_review_status,
            "oocyte_canonical_bound": oocyte_canonical_bound,
            "follicle_histology_status": follicle_histology_status,
            "broad_evidence_manifest_sha256": sha256(
                args.broad_evidence_manifest
            ),
            "pseudobulk_summary_manifest_sha256": sha256(
                args.pseudobulk_summary_manifest
            ),
            "controller_review_artifact_sha256s": {
                key: sha256(path)
                for key, path in sorted(controller_artifact_paths.items())
            },
        }
        packet_rows.append({
            "review_id": review_id,
            "review_mode": queued.get("review_mode", ""),
            "target_broad_label": broad_label,
            "unit_signature": queued.get("unit_signature", ""),
            "evidence_packet_sha256": packet_hash(payload),
            "current_n": current_n,
            "current_competitor_question_n": competitor_n,
            "cross_type_over_recall_question_n": cross_type_n,
            "current_precision_question_n": current_precision_question_n,
            "outside_recall_question_n": recall_n,
            "positive_marker_n_available": positive_marker_n_available,
            "positive_family_n_available": positive_family_n_available,
            "current_differential_n": current_differential_n,
            "recall_differential_n": recall_differential_n,
            "recall_similarity_n": recall_similarity_n,
            "precision_evaluable": str(precision_evaluable).lower(),
            "recall_evaluable": str(recall_evaluable).lower(),
            "molecular_evaluable": str(molecular_evaluable).lower(),
            "spatial_evaluable": str(spatial_evaluable).lower(),
            "ovary_spatial_status": ovary_spatial_status,
            "oocyte_review_status": oocyte_review_status,
            "oocyte_canonical_bound": str(oocyte_canonical_bound).lower(),
            "follicle_histology_status": follicle_histology_status,
            "broad_evidence_manifest_sha256": sha256(
                args.broad_evidence_manifest
            ),
            "pseudobulk_summary_manifest_sha256": sha256(
                args.pseudobulk_summary_manifest
            ),
        })

    if len(packet_rows) != len(queue):
        errors.append("evidence packets do not cover the review queue exactly once")
    args.out.mkdir(parents=True, exist_ok=True)
    index_path = args.out / "broad_cell_type_review_packet_index.tsv"
    write_tsv(index_path, packet_rows, fields=[
        "review_id", "review_mode", "target_broad_label", "unit_signature",
        "evidence_packet_sha256", "current_n",
        "current_competitor_question_n", "cross_type_over_recall_question_n",
        "current_precision_question_n",
        "outside_recall_question_n", "positive_marker_n_available",
        "positive_family_n_available", "current_differential_n",
        "recall_differential_n", "recall_similarity_n",
        "precision_evaluable", "recall_evaluable", "molecular_evaluable",
        "spatial_evaluable", "ovary_spatial_status", "oocyte_review_status",
        "oocyte_canonical_bound", "follicle_histology_status",
        "broad_evidence_manifest_sha256",
        "pseudobulk_summary_manifest_sha256",
    ])
    manifest = {
        "schema_version": "2.2",
        "status": "PASS" if not errors else "BLOCKED",
        "artifact_role": "broad_cell_type_review_evidence_packet_index",
        "formal_membership_written": False,
        "review_manifest": artifact(args.review_manifest),
        "broad_evidence_manifest": artifact(args.broad_evidence_manifest),
        "pseudobulk_summary_manifest": artifact(
            args.pseudobulk_summary_manifest
        ),
        "threshold_registry": artifact(args.threshold_registry),
        "biological_quality_review": quality_artifact,
        "review_queue": artifact(queue_path),
        "packet_index": artifact(index_path),
        "packet_n": len(packet_rows),
        "errors": errors,
    }
    manifest_path = args.out / "broad_cell_type_review_packet_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
