#!/usr/bin/env python3
"""Validate the staged v2.2 controller chain from provisional plan to release."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from lineage_controller_lib import (
    deterministic_membership_hash, read_tsv, sha256, write_manifest,
)


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def bound(manifest: dict, contract: Path) -> bool:
    return (
        manifest.get("controller_version") == "2.2.0"
        and manifest.get("annotation_contract", {}).get("sha256") == sha256(contract)
    )


def artifact_path(record: object, label: str, errors: list[str]) -> Path | None:
    if not isinstance(record, dict):
        errors.append(f"{label} artifact record is missing")
        return None
    path = Path(str(record.get("path", "")))
    if (
        not path.is_file()
        or record.get("sha256") != sha256(path)
    ):
        errors.append(f"{label} artifact is missing or stale")
        return None
    return path


def artifact_index(records: object, label: str, errors: list[str]) -> set[tuple[str, str]]:
    if not isinstance(records, list):
        errors.append(f"{label} artifact registry is missing")
        return set()
    result: set[tuple[str, str]] = set()
    for index, record in enumerate(records, 1):
        path = artifact_path(record, f"{label} {index}", errors)
        if path:
            result.add((str(path.resolve()), sha256(path)))
    if len(result) != len(records):
        errors.append(f"{label} artifact registry contains duplicates")
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", required=True, type=Path)
    ap.add_argument("--whole-manifest", required=True, type=Path)
    ap.add_argument("--cohort-manifest", required=True, action="append", type=Path)
    ap.add_argument("--local-manifest", action="append", type=Path, default=[])
    ap.add_argument("--merge-manifest", required=True, type=Path)
    ap.add_argument("--atlas-manifest", required=True, type=Path)
    ap.add_argument("--final-manifest", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    errors: list[str] = []

    contract = load(args.contract)
    if contract.get("canonical_lineage_controller", {}).get("phase_order") != [
        "whole_tissue_partition", "cluster_cohort_recluster",
        "local_mixed_subcluster_split", "merge_and_freeze_broad",
        "atlas_and_completeness_review", "materialize_final_release",
    ]:
        errors.append("contract does not bind the staged v2.2 architecture")
    analysis_path = artifact_path(
        contract.get("input_scope", {}).get("analysis_set", {}),
        "contract analysis_set", errors,
    )
    analysis_ids = (
        {str(row.get("cell_id", "")) for row in read_tsv(analysis_path)}
        if analysis_path else set()
    )
    if (
        not analysis_ids or "" in analysis_ids
        or len(analysis_ids) != len(read_tsv(analysis_path))
    ):
        errors.append("contract analysis_set is empty or nonexclusive")

    whole = load(args.whole_manifest)
    if not bound(whole, args.contract) or whole.get("phase") != "whole_tissue_partition":
        errors.append("whole-tissue controller manifest is stale or wrong-stage")
    if whole.get("formal_membership_written") is not False or whole.get("release_authority_written") is not False:
        errors.append("whole-tissue phase acquired forbidden release authority")
    plan_path = artifact_path(
        whole.get("cohort_plan", {}), "whole-tissue cohort plan", errors
    )
    if not plan_path:
        plan_rows: list[dict[str, str]] = []
    else:
        plan_rows = read_tsv(plan_path)
    plan_cohorts = {str(row.get("cohort_id", "")): row for row in plan_rows}
    if not plan_cohorts or "" in plan_cohorts or len(plan_cohorts) != len(plan_rows):
        errors.append("whole-tissue plan has empty or duplicate cohort IDs")
    if any(row.get("formal_label_written") != "false" for row in plan_rows):
        errors.append("whole-tissue plan contains a formal label write")
    cluster_map_path = artifact_path(
        whole.get("cluster_membership", {}),
        "whole-tissue cluster membership", errors,
    )
    cluster_map_rows = read_tsv(cluster_map_path) if cluster_map_path else []
    cluster_map_ids = [str(row.get("cell_id", "")) for row in cluster_map_rows]
    if (
        set(cluster_map_ids) != analysis_ids
        or len(cluster_map_ids) != len(set(cluster_map_ids))
    ):
        errors.append("whole-tissue cohort map does not exactly partition analysis_set")
    for cohort_id, row in plan_cohorts.items():
        membership_path = Path(str(row.get("membership_path", "")))
        if (
            not membership_path.is_file()
            or row.get("membership_sha256") != sha256(membership_path)
        ):
            errors.append(f"{cohort_id}: planned cohort membership is missing or stale")
            continue
        planned_ids = {str(item.get("cell_id", "")) for item in read_tsv(membership_path)}
        mapped_ids = {
            str(item.get("cell_id", "")) for item in cluster_map_rows
            if item.get("cohort_id") == cohort_id
            and item.get("source_initial_cluster") == row.get("source_initial_cluster")
        }
        if (
            not planned_ids or planned_ids != mapped_ids
            or len(planned_ids) != int(row.get("n_observations", -1))
        ):
            errors.append(f"{cohort_id}: planned membership differs from cohort map")

    cohort_manifests: dict[str, dict] = {}
    pending_groups: set[tuple[str, str]] = set()
    for path in args.cohort_manifest:
        manifest = load(path)
        cohort_id = str(manifest.get("cohort_id", ""))
        if (
            not bound(manifest, args.contract)
            or manifest.get("phase") != "cluster_cohort_recluster"
            or manifest.get("formal_membership_written") is not False
            or not cohort_id or cohort_id in cohort_manifests
        ):
            errors.append(f"invalid or duplicate second-round cohort manifest: {path}")
            continue
        cohort_manifests[cohort_id] = manifest
        plan_row = plan_cohorts.get(cohort_id, {})
        if manifest.get("source_initial_cluster") != plan_row.get("source_initial_cluster"):
            errors.append(f"{cohort_id}: source initial cluster differs from plan")
        whole_record = manifest.get("whole_tissue_manifest", {})
        if (
            Path(str(whole_record.get("path", ""))).resolve()
            != args.whole_manifest.resolve()
            or whole_record.get("sha256") != sha256(args.whole_manifest)
        ):
            errors.append(f"{cohort_id}: cohort is not bound to this whole-tissue plan")
        query = manifest.get("query_membership", {})
        if (
            query.get("sha256") != plan_row.get("membership_sha256")
            or int(query.get("n_observations", -1))
            != int(plan_row.get("n_observations", -2))
        ):
            errors.append(f"{cohort_id}: cohort query differs from planned membership")
        ancestry_path = Path(str(manifest.get("raw_count_ancestry", {}).get("path", "")))
        if (
            manifest.get("raw_count_assay") == "SCT"
            or not ancestry_path.is_file()
            or manifest.get("raw_count_ancestry", {}).get("sha256") != sha256(ancestry_path)
        ):
            errors.append(f"{cohort_id}: raw-count ancestry is missing or invalid")
        else:
            ancestry = load(ancestry_path)
            underpowered = (
                manifest.get("cohort_status") == "UNDERPOWERED_NOT_EVALUABLE"
            )
            if (
                ancestry.get("status")
                != ("UNDERPOWERED_NOT_EVALUABLE" if underpowered else "PASS")
                or ancestry.get("clustering_path")
                != (
                    "not_run_fewer_than_three_observations"
                    if underpowered else "raw_counts_SCTv2_PCA_SNN_Leiden"
                )
            ):
                errors.append(f"{cohort_id}: raw-count ancestry state is inconsistent")
        adjud_path = Path(str(manifest.get("adjudication", {}).get("path", "")))
        if not adjud_path.is_file() or manifest.get("adjudication", {}).get("sha256") != sha256(adjud_path):
            errors.append(f"{cohort_id}: second-round adjudication is missing or stale")
            continue
        adjud = load(adjud_path)
        if adjud.get("full_catalog_scan") is not True or adjud.get("provisional_broad_visible_during_scoring") is not False:
            errors.append(f"{cohort_id}: scorer was not full-catalog and provisional-blind")
        for key in (
            "selected_cluster_evidence", "fine_candidate_proposals",
            "state_annotation_proposals", "unmodeled", "cohort_outcome",
        ):
            artifact_path(manifest.get(key, {}), f"{cohort_id} {key}", errors)
        pending_path = Path(str(adjud.get("pending_local_split_membership", {}).get("path", "")))
        if pending_path.is_file():
            for row in read_tsv(pending_path):
                pending_groups.add((str(row.get("source_boundary", "")), str(row.get("source_cluster", ""))))
    if set(cohort_manifests) != set(plan_cohorts):
        errors.append("every initial cluster must enter exactly one second-round cohort")

    local_groups: set[tuple[str, str]] = set()
    for path in args.local_manifest:
        manifest = load(path)
        key = (str(manifest.get("source_boundary", "")), str(manifest.get("source_cluster", "")))
        if (
            not bound(manifest, args.contract)
            or manifest.get("phase") != "local_mixed_subcluster_split"
            or manifest.get("formal_membership_written") is not False
            or not all(key) or key in local_groups
        ):
            errors.append(f"invalid or duplicate local mixed-subcluster manifest: {path}")
        trigger = manifest.get("trigger_manifest", {})
        expected_trigger = next(
            (
                cohort_path for cohort_path in args.cohort_manifest
                if load(cohort_path).get("cohort_id") == key[0]
            ),
            None,
        )
        if (
            expected_trigger is None
            or Path(str(trigger.get("path", ""))).resolve()
            != expected_trigger.resolve()
            or trigger.get("sha256") != sha256(expected_trigger)
        ):
            errors.append(f"{key}: local split is not bound to its cohort trigger")
        artifact_path(
            manifest.get("trigger_membership", {}),
            f"{key} local trigger membership", errors,
        )
        artifact_path(
            manifest.get("candidate_membership", {}),
            f"{key} local candidate membership", errors,
        )
        local_groups.add(key)
    if local_groups != pending_groups:
        errors.append("local observation-level splits do not exactly equal mixed subcluster triggers")

    merge = load(args.merge_manifest)
    if not bound(merge, args.contract) or merge.get("phase") != "merge_and_freeze_broad":
        errors.append("broad-freeze controller manifest is stale or wrong-stage")
    if merge.get("formal_broad_membership_written") is not True:
        errors.append("formal broad membership was not frozen at merge")
    broad_freeze_path = Path(
        str(merge.get("broad_freeze", {}).get("path", ""))
    )
    if (
        not broad_freeze_path.is_file()
        or merge.get("broad_freeze", {}).get("sha256") != sha256(broad_freeze_path)
    ):
        errors.append("broad-freeze manifest is missing or stale")
        broad_freeze = {}
    else:
        broad_freeze = load(broad_freeze_path)
    merge_source_paths = {
        (
            str(Path(str(row.get("path", ""))).resolve()),
            str(row.get("sha256", "")),
        )
        for row in broad_freeze.get("candidate_source_manifests", [])
    }
    expected_source_paths = {
        (str(path.resolve()), sha256(path))
        for path in args.cohort_manifest + args.local_manifest
    }
    if merge_source_paths != expected_source_paths:
        errors.append("broad freeze did not use exactly the canonical cohort/local manifests")
    merge_membership = Path(str(merge.get("membership", {}).get("path", "")))
    if not merge_membership.is_file() or merge.get("membership", {}).get("sha256") != sha256(merge_membership):
        errors.append("frozen broad membership is missing or stale")
        merge_rows: list[dict[str, str]] = []
    else:
        merge_rows = read_tsv(merge_membership)
    if {str(row.get("cell_id", "")) for row in merge_rows} != analysis_ids:
        errors.append("frozen broad membership does not exactly cover analysis_set")
    if (
        broad_freeze.get("analysis_membership", {}).get("sha256")
        != (sha256(analysis_path) if analysis_path else "")
    ):
        errors.append("broad freeze does not bind the contract analysis_set")

    carry_map = {
        "cluster_evidence": "selected_cluster_evidence",
        "fine_candidate_proposals": "fine_candidate_proposals",
        "state_annotation_proposals": "state_annotation_proposals",
        "unmodeled_manifests": "unmodeled",
    }
    for merge_key, cohort_key in carry_map.items():
        expected = {
            (
                str(Path(str(manifest.get(cohort_key, {}).get("path", ""))).resolve()),
                str(manifest.get(cohort_key, {}).get("sha256", "")),
            )
            for manifest in cohort_manifests.values()
        }
        observed = artifact_index(
            merge.get(merge_key), f"merge {merge_key}", errors
        )
        if observed != expected:
            errors.append(f"merge {merge_key} differs from canonical cohort artifacts")

    atlas = load(args.atlas_manifest)
    if not bound(atlas, args.contract) or atlas.get("phase") != "atlas_and_completeness_review":
        errors.append("Atlas/completeness manifest is stale or wrong-stage")
    atlas_prerequisite = atlas.get("prerequisite", {})
    if (
        Path(str(atlas_prerequisite.get("path", ""))).resolve()
        != args.merge_manifest.resolve()
        or atlas_prerequisite.get("sha256") != sha256(args.merge_manifest)
    ):
        errors.append("Atlas phase is not bound to the supplied broad-freeze controller manifest")
    atlas_authority_path = artifact_path(
        atlas.get("stage_authority", {}), "Atlas stage authority", errors
    )
    atlas_authority = load(atlas_authority_path) if atlas_authority_path else {}
    if (
        atlas_authority.get("mode") != "stage_authority"
        or atlas_authority.get("phase") != "atlas_and_completeness_review"
        or atlas_authority.get("annotation_contract_sha256") != sha256(args.contract)
        or atlas_authority.get("frozen_broad", {}).get("sha256")
        != merge.get("membership", {}).get("sha256")
    ):
        errors.append("Atlas stage authority does not bind the frozen broad membership")
    for authority_key, merge_key in (
        ("cluster_evidence", "cluster_evidence"),
        ("fine_candidate_proposals", "fine_candidate_proposals"),
        ("state_annotation_proposals", "state_annotation_proposals"),
        ("unmodeled_manifests", "unmodeled_manifests"),
    ):
        if artifact_index(
            atlas_authority.get(authority_key),
            f"Atlas authority {authority_key}", errors,
        ) != artifact_index(
            merge.get(merge_key), f"merge {merge_key} authority carry", errors
        ):
            errors.append(
                f"Atlas authority {authority_key} differs from broad-freeze carry artifacts"
            )
    for key in ("fine_candidate_proposals", "state_annotation_proposals"):
        if artifact_index(atlas.get(key), f"Atlas {key}", errors) != artifact_index(
            merge.get(key), f"merge {key} carry", errors
        ):
            errors.append(f"Atlas {key} differs from broad-freeze carry artifacts")
    atlas_membership = artifact_path(
        atlas.get("membership", {}), "post-Atlas membership", errors
    )
    atlas_rows = read_tsv(atlas_membership) if atlas_membership else []
    if {str(row.get("cell_id", "")) for row in atlas_rows} != analysis_ids:
        errors.append("post-Atlas membership does not exactly cover analysis_set")
    completeness_path = Path(str(atlas.get("completeness", {}).get("path", "")))
    if not completeness_path.is_file() or atlas.get("completeness", {}).get("sha256") != sha256(completeness_path):
        errors.append("post-merge completeness audit is missing or stale")
    elif load(completeness_path).get("status") != "PASS":
        errors.append("post-merge completeness audit did not pass")

    final_controller = load(args.final_manifest)
    if not bound(final_controller, args.contract) or final_controller.get("phase") != "materialize_final_release":
        errors.append("final controller manifest is stale or wrong-stage")
    if final_controller.get("formal_membership_written") is not True:
        errors.append("final phase did not materialize release membership")
    final_prerequisite = final_controller.get("prerequisite", {})
    if (
        Path(str(final_prerequisite.get("path", ""))).resolve()
        != args.atlas_manifest.resolve()
        or final_prerequisite.get("sha256") != sha256(args.atlas_manifest)
    ):
        errors.append("final phase is not bound to the supplied Atlas controller manifest")
    final_authority_path = artifact_path(
        final_controller.get("stage_authority", {}),
        "final stage authority", errors,
    )
    final_authority = load(final_authority_path) if final_authority_path else {}
    if (
        final_authority.get("mode") != "stage_authority"
        or final_authority.get("phase") != "materialize_final_release"
        or final_authority.get("annotation_contract_sha256") != sha256(args.contract)
        or Path(str(final_authority.get("post_atlas_membership", {}).get("path", ""))).resolve()
        != (atlas_membership.resolve() if atlas_membership else Path("/__missing__"))
        or final_authority.get("post_atlas_membership", {}).get("sha256")
        != atlas.get("membership", {}).get("sha256")
        or final_authority.get("prerequisite_manifest", {}).get("sha256")
        != sha256(args.atlas_manifest)
    ):
        errors.append("final stage authority does not bind the post-Atlas prerequisite")
    if artifact_index(
        final_authority.get("state_annotation_proposals"),
        "final authority state proposals", errors,
    ) != artifact_index(
        atlas.get("state_annotation_proposals"),
        "Atlas state proposal carry", errors,
    ):
        errors.append("final authority state proposals differ from Atlas carry artifacts")
    artifact_path(
        final_authority.get("fine_assignments", {}),
        "final parent-locked fine assignments", errors,
    )
    final_path = Path(str(final_controller.get("membership", {}).get("path", "")))
    if not final_path.is_file() or final_controller.get("membership", {}).get("sha256") != sha256(final_path):
        errors.append("final release membership is missing or stale")
        final_rows: list[dict[str, str]] = []
    else:
        final_rows = read_tsv(final_path)
        if final_controller.get("membership", {}).get("semantic_sha256") != deterministic_membership_hash(final_rows):
            errors.append("final release semantic hash is stale")
    final_ids = [row.get("cell_id", "") for row in final_rows]
    if not final_ids or "" in final_ids or len(final_ids) != len(set(final_ids)):
        errors.append("final membership is empty or nonexclusive")
    if set(final_ids) != analysis_ids:
        errors.append("final membership does not exactly cover analysis_set")
    qc_n = sum(row.get("final_state") == "qc_holdout" for row in final_rows)
    qc_fraction = qc_n / len(final_rows) if final_rows else 1.0
    if qc_n >= 50000 or qc_fraction >= 0.10:
        errors.append("final residual QC violates the v2.2 completion threshold")
    if any(row.get("final_state") == "unresolved_biological" for row in final_rows):
        errors.append("unresolved biological members were not typed at final materialization")

    result = {
        "status": "PASS" if not errors else "BLOCKED",
        "schema_version": "2.2",
        "controller_version": "2.2.0",
        "annotation_contract": str(args.contract.resolve()),
        "annotation_contract_sha256": sha256(args.contract),
        "n_initial_clusters": len(plan_cohorts),
        "n_second_round_cohorts": len(cohort_manifests),
        "n_local_mixed_splits": len(local_groups),
        "n_final": len(final_rows),
        "final_membership": (
            {
                "path": str(final_path.resolve()),
                "sha256": sha256(final_path),
                "semantic_sha256": deterministic_membership_hash(final_rows),
            }
            if final_path.is_file() else {}
        ),
        "residual_qc_n": qc_n,
        "residual_qc_fraction": qc_fraction,
        "errors": errors,
    }
    write_manifest(args.out, result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
