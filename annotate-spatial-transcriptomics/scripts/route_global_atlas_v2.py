#!/usr/bin/env python3
"""Authoritative v2 all-cell Atlas comparison and state-aware broad routing."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

from evidence_schema_lib import sha256, validate_artifact_ref


ACCEPTED_TIERS = {"high", "moderate_only"}
QC_STATES = {"qc_holdout", "low_information_qc_holdout", "pending_qc"}
TRUE = {"1", "true", "yes", "pass", "passed"}


def truth(value: object) -> bool:
    return str(value or "").strip().lower() in TRUE


def open_text(path: Path, mode: str):
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(mode, encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def unique(rows: list[dict[str, str]], column: str, label: str) -> dict[str, dict[str, str]]:
    ids = [row.get(column, "").strip() for row in rows]
    if not rows or "" in ids or len(ids) != len(set(ids)):
        raise SystemExit(f"{label} must contain unique nonempty {column}")
    return dict(zip(ids, rows))


def validate_calibration(manifest_path: Path, mapping_path: Path, policy: dict) -> dict[str, dict]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "2.0" or manifest.get("status") != "CALIBRATED_TIERED_EVIDENCE_ONLY":
        raise SystemExit("Atlas calibration manifest is not a v2 tiered calibration")
    if manifest.get("heldout_origin") != policy.get("calibration_origin"):
        raise SystemExit("Atlas calibration origin differs from the active workflow profile")
    origin_path = Path(manifest["calibration_origin_manifest"])
    if not origin_path.is_file() or manifest.get("calibration_origin_manifest_sha256") != sha256(origin_path):
        raise SystemExit("Atlas calibration origin manifest is stale")
    origin = json.loads(origin_path.read_text(encoding="utf-8"))
    if (origin.get("status") != "PASS" or origin.get("heldout_origin") != policy.get("calibration_origin")
            or origin.get("reference_self_classification") is not False
            or int(origin.get("anchor_target_overlap", -1)) != 0):
        raise SystemExit("Atlas calibration origin does not satisfy the current-query held-out contract")
    artifacts = manifest.get("artifacts", {})
    for name, record in artifacts.items():
        _, errors = validate_artifact_ref(manifest_path.parent, record, f"calibration {name}")
        if errors:
            raise SystemExit("; ".join(errors))
    query_record = artifacts.get("query_mapping", {})
    if Path(query_record.get("path", "")).resolve() != mapping_path.resolve() or query_record.get("sha256") != sha256(mapping_path):
        raise SystemExit("Atlas mapping is not the hash-bound calibrated query mapping")
    cumulative_path = Path(artifacts["heldout_cumulative_validation"]["path"])
    rows = read_tsv(cumulative_path)
    result: dict[str, dict] = {}
    min_n = int(policy.get("minimum_heldout_observations_per_return_class", 30))
    max_ece = float(policy.get("maximum_expected_calibration_error", 0.05))
    min_precision = max(float(policy.get("moderate_or_higher_target_precision", 0.90)), 1.0 - max_ece)
    for row in rows:
        if row.get("cumulative_tier") != "moderate_or_higher":
            continue
        label = row.get("predicted_label", "")
        n = int(float(row.get("validation_n", 0) or 0))
        precision = float(row.get("validation_precision", 0) or 0)
        result[label] = {"n": n, "precision": precision, "eligible": n >= min_n and precision >= min_precision}
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell-ledger", required=True, type=Path)
    ap.add_argument("--atlas-mapping", required=True, type=Path)
    ap.add_argument("--calibration-manifest", required=True, type=Path)
    ap.add_argument("--workflow-profile", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--cell-id-col", default="cell_id")
    ap.add_argument("--analysis-scope-col", default="analysis_scope")
    ap.add_argument("--state-col", default="final_state")
    ap.add_argument("--broad-col", default="final_broad_label")
    ap.add_argument("--fine-col", default="final_fine_label")
    ap.add_argument("--cluster-col", default="source_cluster")
    ap.add_argument("--atlas-label-col", default="predicted_label")
    ap.add_argument("--atlas-tier-col", default="mapping_tier")
    ap.add_argument("--ood-col", default="out_of_distribution")
    ap.add_argument("--ontology-conflict-col", default="ontology_conflict")
    ap.add_argument("--min-discordant-n", type=int, default=30)
    ap.add_argument("--min-discordant-fraction", type=float, default=0.20)
    args = ap.parse_args()

    profile = json.loads(args.workflow_profile.read_text(encoding="utf-8"))
    policy = profile.get("external_reference_policy", {})
    crosswalk = {str(k): str(v) for k, v in policy.get("primary_public_atlas_label_crosswalk", {}).items()}
    allowed_returns = set(policy.get("primary_public_atlas_scope", []))
    excluded_returns = set(policy.get("primary_public_atlas_exclusions", []))
    if not crosswalk or allowed_returns & excluded_returns:
        raise SystemExit("active workflow profile has an invalid Atlas crosswalk/scope")
    calibration = validate_calibration(args.calibration_manifest, args.atlas_mapping, policy)

    ledger_rows = read_tsv(args.cell_ledger)
    mapping_rows = read_tsv(args.atlas_mapping)
    ledger = unique(ledger_rows, args.cell_id_col, "cell ledger")
    mapping = unique(mapping_rows, args.cell_id_col, "Atlas mapping")
    analysis_ids = [row[args.cell_id_col] for row in ledger_rows if row.get(args.analysis_scope_col) == "analysis_set"]
    if set(mapping) != set(analysis_ids) or len(mapping) != len(analysis_ids):
        raise SystemExit("Atlas mapping must cover the analysis_set exactly once")

    cluster_sizes = Counter(ledger[cell].get(args.cluster_col, "") for cell in analysis_ids)
    primary_group_sizes = Counter(
        (ledger[cell].get(args.cluster_col, ""), ledger[cell].get(args.broad_col, "").strip())
        for cell in analysis_ids if ledger[cell].get(args.broad_col, "").strip()
    )
    rows: list[dict] = []
    review_candidates: dict[tuple[str, str, str, str], list[str]] = defaultdict(list)
    cell_review_key: dict[str, tuple[str, str, str, str]] = {}
    for cell in analysis_ids:
        current = ledger[cell]
        atlas = mapping[cell]
        state = current.get(args.state_col, "").strip()
        primary = current.get(args.broad_col, "").strip()
        source_label = atlas.get(args.atlas_label_col, "").strip()
        mapped = crosswalk.get(source_label, "")
        tier = atlas.get(args.atlas_tier_col, "").strip()
        ood = truth(atlas.get(args.ood_col, ""))
        ontology = truth(atlas.get(args.ontology_conflict_col, ""))
        class_calibrated = calibration.get(source_label, {}).get("eligible", False)
        in_scope = mapped in allowed_returns and mapped not in excluded_returns
        accepted = tier in ACCEPTED_TIERS and class_calibrated and bool(mapped) and not ood and not ontology and in_scope
        is_qc = not primary and state in QC_STATES
        if is_qc and accepted:
            route = "direct_qc_broad_return"
            proposed_state, proposed_broad = "defined_broad_only", mapped
            review = False
        elif is_qc:
            route = "retain_qc"
            proposed_state, proposed_broad = state, ""
            review = False
        elif primary:
            proposed_state, proposed_broad = state, primary
            if accepted and mapped == primary:
                route, review = "defined_label_atlas_agreement", False
            elif accepted and mapped != primary:
                route, review = "defined_label_disagreement_candidate", False
                key = ("broad_disagreement", current.get(args.cluster_col, ""), primary, mapped)
                review_candidates[key].append(cell)
                cell_review_key[cell] = key
            elif ood or ontology or (source_label and not mapped):
                reason = "ontology_conflict" if ontology else "out_of_distribution" if ood else "unmapped_atlas_label"
                route, review = "defined_label_ood_ontology_or_crosswalk_candidate", False
                key = (reason, current.get(args.cluster_col, ""), primary, mapped or source_label)
                review_candidates[key].append(cell)
                cell_review_key[cell] = key
            else:
                route, review = "defined_label_weak_atlas_logged", False
        else:
            route = "unlabeled_non_qc_state_review"
            proposed_state, proposed_broad, review = state, "", True
        rows.append({
            "cell_id": cell,
            "source_cluster": current.get(args.cluster_col, ""),
            "primary_state": state,
            "primary_broad": primary,
            "primary_fine": current.get(args.fine_col, ""),
            "atlas_source_label": source_label,
            "atlas_broad": mapped,
            "atlas_tier": tier,
            "atlas_class_calibrated": str(class_calibrated).lower(),
            "atlas_scope_pass": str(in_scope).lower(),
            "out_of_distribution": str(ood).lower(),
            "ontology_conflict": str(ontology).lower(),
            "atlas_state_route": route,
            "proposed_state": proposed_state,
            "proposed_broad_label": proposed_broad,
            "proposed_fine_label": "",
            "fine_anchor_eligible": "false",
            "review_required": str(review).lower(),
            "review_id": "",
            "writeback_status": "proposal_only_requires_atomic_commit",
        })

    review_rows = []
    triggered: dict[tuple[str, str, str, str], str] = {}
    for key, cells in sorted(review_candidates.items()):
        reason, cluster, primary, mapped = key
        denominator = primary_group_sizes[(cluster, primary)]
        fraction = len(cells) / denominator
        if len(cells) >= args.min_discordant_n and fraction >= args.min_discordant_fraction:
            review_id = f"atlas_discrepancy__{len(review_rows)+1:04d}"
            triggered[key] = review_id
            review_rows.append({
                "review_id": review_id, "trigger_reason": reason, "source_cluster": cluster,
                "primary_broad": primary, "atlas_broad": mapped,
                "n_trigger": len(cells), "primary_group_n": denominator,
                "cluster_n": cluster_sizes[cluster],
                "trigger_fraction": f"{fraction:.8f}",
                "review_scope": "complete_cluster_or_frozen_cohort",
                "required_action": "one_pass_orthogonal_query_evidence_review",
            })
    for row in rows:
        key = cell_review_key.get(row["cell_id"])
        if key in triggered:
            row["review_required"] = "true"
            row["review_id"] = triggered[key]

    args.out.mkdir(parents=True, exist_ok=True)
    route_path = args.out / "atlas_state_routing.tsv.gz"
    review_path = args.out / "atlas_discrepancy_review_queue.tsv"
    write_tsv(route_path, rows, list(rows[0]))
    write_tsv(review_path, review_rows, ["review_id", "trigger_reason", "source_cluster", "primary_broad", "atlas_broad", "n_trigger", "primary_group_n", "cluster_n", "trigger_fraction", "review_scope", "required_action"])
    census = Counter(row["atlas_state_route"] for row in rows)
    manifest = {
        "schema_version": "2.0",
        "status": "REVIEW_REQUIRED" if review_rows or any(row["review_required"] == "true" for row in rows) else "PASS_NO_REVIEW",
        "authoritative_router": "route_global_atlas_v2.py",
        "query_scope": "complete_analysis_set",
        "n_analysis_set": len(rows),
        "route_census": dict(sorted(census.items())),
        "material_review_thresholds": {"minimum_n": args.min_discordant_n, "minimum_primary_group_fraction": args.min_discordant_fraction},
        "cell_ledger": {"path": str(args.cell_ledger.resolve()), "sha256": sha256(args.cell_ledger)},
        "atlas_mapping": {"path": str(args.atlas_mapping.resolve()), "sha256": sha256(args.atlas_mapping)},
        "calibration_manifest": {"path": str(args.calibration_manifest.resolve()), "sha256": sha256(args.calibration_manifest)},
        "workflow_profile": {"path": str(args.workflow_profile.resolve()), "sha256": sha256(args.workflow_profile)},
        "artifacts": {
            "routing": {"path": str(route_path.resolve()), "sha256": sha256(route_path)},
            "review_queue": {"path": str(review_path.resolve()), "sha256": sha256(review_path)},
        },
        "ledger_writeback_performed": False,
        "fine_anchor_eligible": False,
    }
    (args.out / "atlas_state_routing_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 2 if manifest["status"] == "REVIEW_REQUIRED" else 0


if __name__ == "__main__":
    raise SystemExit(main())
