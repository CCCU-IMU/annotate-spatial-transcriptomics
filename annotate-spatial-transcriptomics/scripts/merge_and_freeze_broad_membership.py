#!/usr/bin/env python3
"""Freeze formal broad membership only after every second-round cohort closes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from evidence_schema_lib import sha256
from lineage_controller_lib import (
    apply_candidate_context, candidate_can_release,
    candidate_can_support_broad_review, catalog_candidates,
    deterministic_membership_hash, read_tsv, write_tsv,
)


def artifact_ok(record: dict, path: Path) -> bool:
    return (
        path.is_file()
        and Path(str(record.get("path", ""))).resolve() == path.resolve()
        and record.get("sha256") == sha256(path)
    )


def require_artifact(record: dict, label: str) -> Path:
    path = Path(str(record.get("path", "")))
    if not artifact_ok(record, path):
        raise SystemExit(f"{label} is missing or stale")
    return path


def validate_cohort_source(source: dict, contract: dict, contract_path: Path) -> None:
    adjudication_path = require_artifact(
        source.get("adjudication", {}), "cohort adjudication"
    )
    adjudication = json.loads(adjudication_path.read_text(encoding="utf-8"))
    if (
        adjudication.get("stage") != "cluster_cohort_recluster"
        or adjudication.get("formal_membership_written") is not False
        or adjudication.get("full_catalog_scan") is not True
    ):
        raise SystemExit("cohort adjudication is not a full-catalog proposal stage")
    scoring_path = require_artifact(
        source.get("selected_scoring", {}), "cohort selected scoring"
    )
    scoring = json.loads(scoring_path.read_text(encoding="utf-8"))
    installed_scorer = contract.get("canonical_lineage_controller", {}).get(
        "scripts", {}
    ).get("run_observation_lineage_scoring.R", {})
    if (
        scoring.get("controller_version") != "2.2.0"
        or scoring.get("historical_labels_read") is not False
        or scoring.get("scorer", {}).get("sha256") != installed_scorer.get("sha256")
    ):
        raise SystemExit("cohort selected scoring is not canonical and label blind")
    ancestry_path = require_artifact(
        source.get("raw_count_ancestry", {}), "cohort raw-count ancestry"
    )
    ancestry = json.loads(ancestry_path.read_text(encoding="utf-8"))
    snapshot = contract.get("selected_input_snapshot", {})
    underpowered = source.get("cohort_status") == "UNDERPOWERED_NOT_EVALUABLE"
    ancestry_valid = (
        ancestry.get("status")
        == ("UNDERPOWERED_NOT_EVALUABLE" if underpowered else "PASS")
        and ancestry.get("raw_count_assay") != "SCT"
        and ancestry.get("source_runtime_snapshot", {}).get("sha256")
        == snapshot.get("sha256")
        and (
            ancestry.get("clustering_path") == "raw_counts_SCTv2_PCA_SNN_Leiden"
            if not underpowered
            else ancestry.get("clustering_path")
            == "not_run_fewer_than_three_observations"
        )
    )
    if not ancestry_valid:
        raise SystemExit("cohort did not prove project-local raw-count reclustering")


def validate_local_source(source: dict, contract_path: Path) -> None:
    trigger_path = require_artifact(
        source.get("trigger_manifest", {}), "local split trigger manifest"
    )
    trigger = json.loads(trigger_path.read_text(encoding="utf-8"))
    if (
        trigger.get("status") != "PASS"
        or trigger.get("controller_version") != "2.2.0"
        or trigger.get("phase") != "cluster_cohort_recluster"
        or trigger.get("annotation_contract", {}).get("sha256")
        != sha256(contract_path)
    ):
        raise SystemExit("local split trigger is not a bound second-round cohort")
    require_artifact(source.get("trigger_membership", {}), "local split trigger membership")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", required=True, type=Path)
    ap.add_argument("--stage-authority", required=True, type=Path)
    ap.add_argument("--analysis-membership", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument(
        "--candidate-membership", required=True, type=Path, action="append",
        help="repeat for each terminal cohort/local replacement candidate partition",
    )
    ap.add_argument(
        "--candidate-source-manifest", required=True, type=Path, action="append",
        help="canonical cohort/local controller manifest for each candidate membership",
    )
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    authority = json.loads(args.stage_authority.read_text(encoding="utf-8"))
    if (
        authority.get("mode") != "stage_authority"
        or authority.get("phase") != "merge_and_freeze_broad"
        or authority.get("annotation_contract_sha256") != sha256(args.contract)
    ):
        raise SystemExit("stage authority does not permit broad freeze")
    if contract.get("canonical_lineage_controller", {}).get("controller_version") != "2.2.0":
        raise SystemExit("contract does not bind v2.2.0 controller")
    authority_records = {"candidate_catalog": args.catalog}
    if args.context_evidence:
        authority_records["context_evidence"] = args.context_evidence
    for key, path in authority_records.items():
        record = authority.get(key, {})
        if not artifact_ok(record, path):
            raise SystemExit(f"broad-freeze authority differs for {key}")
    candidates = catalog_candidates(
        json.loads(args.catalog.read_text(encoding="utf-8"))
    )
    context_summary = apply_candidate_context(
        candidates,
        read_tsv(args.context_evidence) if args.context_evidence else [],
    )
    if len(args.candidate_membership) != len(args.candidate_source_manifest):
        raise SystemExit("each candidate membership requires one canonical source manifest")
    authority_sources = authority.get("candidate_source_manifests", [])
    if len(authority_sources) != len(args.candidate_source_manifest):
        raise SystemExit("stage authority does not bind every candidate source manifest")
    authority_source_index = {
        (Path(str(record.get("path", ""))).resolve(), str(record.get("sha256", "")))
        for record in authority_sources
    }
    supplied_source_index = {
        (path.resolve(), sha256(path)) for path in args.candidate_source_manifest
    }
    if authority_source_index != supplied_source_index:
        raise SystemExit("candidate source manifests differ from merge stage authority")

    analysis_rows = read_tsv(args.analysis_membership)
    analysis_ids = [str(row.get("cell_id", "")) for row in analysis_rows]
    if not analysis_ids or "" in analysis_ids or len(analysis_ids) != len(set(analysis_ids)):
        raise SystemExit("analysis membership must contain unique nonempty cell_id")
    analysis_id_set = set(analysis_ids)

    proposals: dict[str, dict[str, str]] = {}
    source_records: list[dict[str, str]] = []
    for path, source_path in zip(
        args.candidate_membership, args.candidate_source_manifest
    ):
        source = json.loads(source_path.read_text(encoding="utf-8"))
        phase = str(source.get("phase", ""))
        if (
            source.get("status") != "PASS"
            or source.get("controller_version") != "2.2.0"
            or source.get("annotation_contract", {}).get("sha256")
            != sha256(args.contract)
            or phase not in {
                "cluster_cohort_recluster", "local_mixed_subcluster_split"
            }
            or source.get("formal_membership_written") is not False
        ):
            raise SystemExit("candidate membership source is not a canonical proposal-stage manifest")
        if phase == "cluster_cohort_recluster":
            validate_cohort_source(source, contract, args.contract)
        else:
            validate_local_source(source, args.contract)
        membership_key = (
            "base_candidate_membership"
            if phase == "cluster_cohort_recluster"
            else "candidate_membership"
        )
        record = source.get(membership_key, {})
        if not artifact_ok(record, path):
            raise SystemExit("candidate membership differs from its canonical source manifest")
        source_records.append({
            "phase": phase,
            "path": str(source_path.resolve()),
            "sha256": sha256(source_path),
            "membership_path": str(path.resolve()),
            "membership_sha256": sha256(path),
        })
        rows = read_tsv(path)
        for row in rows:
            cell = str(row.get("cell_id", ""))
            if not cell or cell not in analysis_id_set:
                raise SystemExit(f"candidate membership contains invalid observation: {cell}")
            if cell in proposals:
                raise SystemExit(
                    "candidate memberships overlap; local replacement must replace, not duplicate, base membership"
                )
            if row.get("proposed_state") not in {"broad_candidate", "unresolved_biological"}:
                raise SystemExit("candidate membership has an invalid pre-freeze state")
            if row.get("proposed_state") == "broad_candidate" and not row.get("proposed_broad_label"):
                raise SystemExit("broad candidate lacks proposed broad label")
            if row.get("proposed_state") == "broad_candidate":
                candidate_id = str(row.get("candidate_id", ""))
                broad = str(row.get("proposed_broad_label", ""))
                candidate = candidates.get(candidate_id, {})
                if (
                    not candidate_can_release(candidate)
                    or not candidate_can_support_broad_review(candidate)
                    or str(candidate.get("release_broad_label", "")) != broad
                ):
                    raise SystemExit(
                        f"broad freeze rejected a context-ineligible or mismatched candidate: {candidate_id}"
                    )
            proposals[cell] = row
    if set(proposals) != analysis_id_set:
        missing = len(analysis_id_set - set(proposals))
        extra = len(set(proposals) - analysis_id_set)
        raise SystemExit(
            f"terminal second-round candidates do not exactly cover analysis set: missing={missing} extra={extra}"
        )

    frozen: list[dict[str, object]] = []
    for cell in sorted(analysis_ids):
        row = proposals[cell]
        broad = str(row.get("proposed_broad_label", ""))
        unresolved = str(row.get("unresolved_reason", ""))
        frozen.append({
            "cell_id": cell,
            "analysis_scope": "analysis_set",
            "source_boundary": row.get("source_boundary", ""),
            "source_cluster": row.get("source_cluster", ""),
            "candidate_id": row.get("candidate_id", ""),
            "final_state": "defined_broad_only" if broad else "unresolved_biological",
            "final_broad_label": broad,
            "final_fine_label": "",
            "confidence": row.get("confidence", ""),
            "assignment_origin": row.get("assignment_origin", ""),
            "qc_reason": "",
            "unresolved_reason": unresolved,
            "broad_frozen": "true",
            "fine_anchor_eligible": "false",
        })
    args.out.mkdir(parents=True, exist_ok=True)
    membership_path = args.out / "frozen_broad_membership.tsv.gz"
    write_tsv(membership_path, frozen)
    manifest = {
        "status": "PASS",
        "schema_version": "2.2",
        "stage": "merge_and_freeze_broad",
        "formal_broad_membership_written": True,
        "fine_membership_written": False,
        "n_observations": len(frozen),
        "n_defined_broad": sum(bool(row["final_broad_label"]) for row in frozen),
        "n_unresolved_biological": sum(not bool(row["final_broad_label"]) for row in frozen),
        "membership": {
            "path": str(membership_path.resolve()),
            "sha256": sha256(membership_path),
            "semantic_sha256": deterministic_membership_hash(frozen),
        },
        "analysis_membership": {
            "path": str(args.analysis_membership.resolve()),
            "sha256": sha256(args.analysis_membership),
        },
        "candidate_memberships": [
            {"path": str(path.resolve()), "sha256": sha256(path)}
            for path in args.candidate_membership
        ],
        "candidate_source_manifests": source_records,
        "candidate_catalog": {
            "path": str(args.catalog.resolve()),
            "sha256": sha256(args.catalog),
        },
        "context_evidence": (
            {
                "path": str(args.context_evidence.resolve()),
                "sha256": sha256(args.context_evidence),
            }
            if args.context_evidence else None
        ),
        "context_release_eligibility": context_summary,
    }
    (args.out / "broad_freeze_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
