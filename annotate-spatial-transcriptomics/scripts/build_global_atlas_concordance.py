#!/usr/bin/env python3
"""Compare one calibrated all-cell Atlas mapping with the frozen primary state."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


ACCEPTED_TIERS = {"high", "moderate_only"}
ALL_TIERS = ACCEPTED_TIERS | {"low_reject"}
QC_STATES = {"qc_holdout", "low_information_qc_holdout", "qc_reject", "pending_review", "pending_qc"}


def truth(value: object) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "pass", "passed"}


def open_text(path: Path, mode: str):
    return gzip.open(path, mode, newline="", encoding="utf-8") if path.suffix == ".gz" else path.open(mode, newline="", encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_unique(rows: list[dict[str, str]], column: str, label: str) -> dict[str, dict[str, str]]:
    values = [row.get(column, "").strip() for row in rows]
    if not rows or any(not value for value in values) or len(values) != len(set(values)):
        raise SystemExit(f"{label} must contain unique nonempty {column} values")
    return dict(zip(values, rows))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-ledger", required=True, type=Path)
    parser.add_argument("--atlas-mapping", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--analysis-scope-col", default="analysis_scope")
    parser.add_argument("--primary-state-col", default="final_state")
    parser.add_argument("--primary-broad-col", default="final_broad_label")
    parser.add_argument("--primary-fine-col", default="final_fine_label")
    parser.add_argument("--cluster-col", default="source_cluster")
    parser.add_argument("--atlas-label-col", default="predicted_label")
    parser.add_argument("--atlas-tier-col", default="mapping_tier")
    parser.add_argument("--consensus-col", default="consensus_pass")
    parser.add_argument("--ood-col", default="out_of_distribution")
    parser.add_argument("--ontology-conflict-col", default="ontology_conflict")
    parser.add_argument("--min-discordant-n", type=int, default=30)
    parser.add_argument("--min-discordant-fraction", type=float, default=0.20)
    parser.add_argument("--ood-min-n", type=int, default=5)
    parser.add_argument("--ood-min-fraction", type=float, default=0.50)
    parser.add_argument("--qc-routing-policy", choices=["state_aware_calibrated", "legacy_multichannel"],
                        default="state_aware_calibrated")
    parser.add_argument("--atlas-direct-return-excluded-label", action="append", default=[])
    args = parser.parse_args()
    if args.min_discordant_n < 1 or not 0 < args.min_discordant_fraction <= 1:
        raise SystemExit("invalid material-disagreement threshold")
    if args.ood_min_n < 1 or not 0 < args.ood_min_fraction <= 1:
        raise SystemExit("invalid OOD threshold")

    ledger_rows = read_tsv(args.cell_ledger)
    mapping_rows = read_tsv(args.atlas_mapping)
    ledger_by_id = require_unique(ledger_rows, args.cell_id_col, "cell ledger")
    mapping_by_id = require_unique(mapping_rows, args.cell_id_col, "Atlas mapping")
    required_ledger = {args.analysis_scope_col, args.primary_state_col, args.primary_broad_col, args.cluster_col}
    required_mapping = {args.atlas_label_col, args.atlas_tier_col, args.ood_col}
    if args.qc_routing_policy == "legacy_multichannel":
        required_mapping.add(args.consensus_col)
    if not required_ledger.issubset(ledger_rows[0]):
        raise SystemExit(f"cell ledger lacks columns: {sorted(required_ledger - set(ledger_rows[0]))}")
    if not required_mapping.issubset(mapping_rows[0]):
        raise SystemExit(f"Atlas mapping lacks columns: {sorted(required_mapping - set(mapping_rows[0]))}")
    analysis_ids = [row[args.cell_id_col].strip() for row in ledger_rows if row.get(args.analysis_scope_col) == "analysis_set"]
    if set(mapping_by_id) != set(analysis_ids) or len(mapping_by_id) != len(analysis_ids):
        raise SystemExit("Atlas mapping must cover the analysis_set exactly once")
    unknown_tiers = sorted({row.get(args.atlas_tier_col, "") for row in mapping_rows} - ALL_TIERS)
    if unknown_tiers:
        raise SystemExit(f"Atlas mapping contains unknown tiers: {unknown_tiers}")

    cluster_sizes = Counter(ledger_by_id[cell_id].get(args.cluster_col, "") for cell_id in analysis_ids)
    output_rows: list[dict] = []
    discordance_groups: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    ontology_groups: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    ood_groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for cell_id in analysis_ids:
        primary = ledger_by_id[cell_id]
        mapping = mapping_by_id[cell_id]
        primary_broad = primary.get(args.primary_broad_col, "").strip()
        primary_fine = primary.get(args.primary_fine_col, "").strip() if args.primary_fine_col in primary else ""
        primary_state = primary.get(args.primary_state_col, "").strip()
        cluster = primary.get(args.cluster_col, "").strip()
        atlas_broad = mapping.get(args.atlas_label_col, "").strip()
        tier = mapping.get(args.atlas_tier_col, "").strip()
        accepted_tier = tier in ACCEPTED_TIERS
        consensus = truth(mapping.get(args.consensus_col, ""))
        ood = truth(mapping.get(args.ood_col, ""))
        ontology = truth(mapping.get(args.ontology_conflict_col, "")) if args.ontology_conflict_col in mapping else False
        atlas_scope_pass = atlas_broad.lower() not in {
            label.strip().lower() for label in args.atlas_direct_return_excluded_label
        }
        is_qc = not primary_broad and primary_state in QC_STATES
        if is_qc:
            confidence_gate = accepted_tier and (
                consensus if args.qc_routing_policy == "legacy_multichannel" else True
            )
            if confidence_gate and atlas_broad and not ood and not ontology and atlas_scope_pass:
                status = "qc_writeback_candidate"
            else:
                status = "qc_reject"
        elif not primary_broad:
            status = "unresolved_non_qc_error"
        elif ontology:
            status = "defined_ontology_conflict"
            ontology_groups[(cluster, primary_broad, atlas_broad)].append(cell_id)
        elif ood:
            status = "defined_ood_candidate"
            ood_groups[(cluster, primary_broad)].append(cell_id)
        elif accepted_tier and atlas_broad == primary_broad:
            status = "defined_concordant"
        elif accepted_tier and atlas_broad and atlas_broad != primary_broad:
            status = "defined_discordant_candidate"
            discordance_groups[(cluster, primary_broad, atlas_broad)].append(cell_id)
        else:
            status = "defined_weak_challenge"
        output_rows.append({
            "cell_id": cell_id,
            "source_cluster": cluster,
            "primary_state": primary_state,
            "primary_broad": primary_broad,
            "primary_fine": primary_fine,
            "atlas_broad": atlas_broad,
            "atlas_tier": tier,
            "consensus_pass": str(consensus).lower(),
            "atlas_scope_pass": str(atlas_scope_pass).lower(),
            "out_of_distribution": str(ood).lower(),
            "ontology_conflict": str(ontology).lower(),
            "comparison_status": status,
            "review_required": "false",
            "review_id": "",
        })
    if any(row["comparison_status"] == "unresolved_non_qc_error" for row in output_rows):
        raise SystemExit("one or more analysis cells lack both a primary broad label and a recognized QC state")

    review_rows: list[dict] = []
    triggered: dict[tuple[str, str, str, str], str] = {}
    for (cluster, primary_broad, atlas_broad), ids in sorted(discordance_groups.items()):
        fraction = len(ids) / cluster_sizes[cluster]
        if len(ids) >= args.min_discordant_n and fraction >= args.min_discordant_fraction:
            review_id = f"atlas_discordance__{len(review_rows) + 1:04d}"
            triggered[("defined_discordant_candidate", cluster, primary_broad, atlas_broad)] = review_id
            review_rows.append({
                "review_id": review_id,
                "trigger_type": "material_broad_disagreement",
                "source_cluster": cluster,
                "primary_broad": primary_broad,
                "atlas_broad": atlas_broad,
                "n_trigger": len(ids),
                "cluster_n": cluster_sizes[cluster],
                "trigger_fraction": f"{fraction:.8f}",
                "review_scope": "complete_cluster_or_frozen_cohort",
                "required_action": "one_pass_orthogonal_query_evidence_review",
            })
    for (cluster, primary_broad, atlas_broad), ids in sorted(ontology_groups.items()):
        fraction = len(ids) / cluster_sizes[cluster]
        if len(ids) >= args.min_discordant_n and fraction >= args.min_discordant_fraction:
            review_id = f"atlas_ontology__{len(review_rows) + 1:04d}"
            triggered[("defined_ontology_conflict", cluster, primary_broad, atlas_broad)] = review_id
            review_rows.append({
                "review_id": review_id,
                "trigger_type": "material_ontology_conflict",
                "source_cluster": cluster,
                "primary_broad": primary_broad,
                "atlas_broad": atlas_broad,
                "n_trigger": len(ids),
                "cluster_n": cluster_sizes[cluster],
                "trigger_fraction": f"{fraction:.8f}",
                "review_scope": "complete_cluster_or_frozen_cohort",
                "required_action": "resolve_crosswalk_then_orthogonal_query_evidence_review",
            })
    for (cluster, primary_broad), ids in sorted(ood_groups.items()):
        fraction = len(ids) / cluster_sizes[cluster]
        if len(ids) >= args.ood_min_n and fraction >= args.ood_min_fraction:
            review_id = f"atlas_ood__{len(review_rows) + 1:04d}"
            triggered[("defined_ood_candidate", cluster, primary_broad, "__OOD__")] = review_id
            review_rows.append({
                "review_id": review_id,
                "trigger_type": "coherent_out_of_distribution",
                "source_cluster": cluster,
                "primary_broad": primary_broad,
                "atlas_broad": "",
                "n_trigger": len(ids),
                "cluster_n": cluster_sizes[cluster],
                "trigger_fraction": f"{fraction:.8f}",
                "review_scope": "complete_cluster_or_frozen_cohort",
                "required_action": "open_world_unknown_lineage_review",
            })
    for row in output_rows:
        if row["comparison_status"] == "defined_discordant_candidate":
            key = (row["comparison_status"], row["source_cluster"], row["primary_broad"], row["atlas_broad"])
        elif row["comparison_status"] == "defined_ontology_conflict":
            key = (row["comparison_status"], row["source_cluster"], row["primary_broad"], row["atlas_broad"])
        elif row["comparison_status"] == "defined_ood_candidate":
            key = (row["comparison_status"], row["source_cluster"], row["primary_broad"], "__OOD__")
        else:
            continue
        if key in triggered:
            row["review_required"] = "true"
            row["review_id"] = triggered[key]

    args.out.mkdir(parents=True, exist_ok=True)
    all_cell_path = args.out / "all_cell_atlas_concordance.tsv.gz"
    cluster_path = args.out / "cluster_atlas_concordance.tsv"
    review_path = args.out / "discrepancy_review_queue.tsv"
    qc_writeback_path = args.out / "qc_writeback_membership.tsv"
    qc_reject_path = args.out / "qc_reject_membership.tsv"
    fields = list(output_rows[0])
    write_tsv(all_cell_path, output_rows, fields)
    grouped = Counter((row["source_cluster"], row["primary_broad"], row["atlas_broad"], row["comparison_status"]) for row in output_rows)
    cluster_rows = [{
        "source_cluster": key[0], "primary_broad": key[1], "atlas_broad": key[2],
        "comparison_status": key[3], "n": count, "cluster_n": cluster_sizes[key[0]],
        "fraction_of_cluster": f"{count / cluster_sizes[key[0]]:.8f}",
    } for key, count in sorted(grouped.items())]
    write_tsv(cluster_path, cluster_rows, ["source_cluster", "primary_broad", "atlas_broad", "comparison_status", "n", "cluster_n", "fraction_of_cluster"])
    write_tsv(review_path, review_rows, ["review_id", "trigger_type", "source_cluster", "primary_broad", "atlas_broad", "n_trigger", "cluster_n", "trigger_fraction", "review_scope", "required_action"])
    write_tsv(qc_writeback_path, [{"cell_id": row["cell_id"]} for row in output_rows if row["comparison_status"] == "qc_writeback_candidate"], ["cell_id"])
    write_tsv(qc_reject_path, [{"cell_id": row["cell_id"]} for row in output_rows if row["comparison_status"] == "qc_reject"], ["cell_id"])
    counts = Counter(row["comparison_status"] for row in output_rows)
    manifest = {
        "status": "REVIEW_REQUIRED" if review_rows else "PASS_NO_REVIEW",
        "query_scope": "complete_analysis_set",
        "mapping_mode": "single_all_cell_broad_mapping",
        "qc_writeback_policy": args.qc_routing_policy,
        "atlas_direct_return_excluded_labels": args.atlas_direct_return_excluded_label,
        "cell_ledger": str(args.cell_ledger.resolve()),
        "cell_ledger_sha256": sha256(args.cell_ledger),
        "atlas_mapping": str(args.atlas_mapping.resolve()),
        "atlas_mapping_sha256": sha256(args.atlas_mapping),
        "n_analysis_set": len(analysis_ids),
        "comparison_counts": dict(sorted(counts.items())),
        "review_queue_n": len(review_rows),
        "thresholds": {
            "min_discordant_n": args.min_discordant_n,
            "min_discordant_fraction": args.min_discordant_fraction,
            "ood_min_n": args.ood_min_n,
            "ood_min_fraction": args.ood_min_fraction,
        },
        "artifacts": {
            "all_cell_concordance": {"path": str(all_cell_path.resolve()), "sha256": sha256(all_cell_path)},
            "cluster_concordance": {"path": str(cluster_path.resolve()), "sha256": sha256(cluster_path)},
            "discrepancy_review_queue": {"path": str(review_path.resolve()), "sha256": sha256(review_path)},
            "qc_writeback_membership": {"path": str(qc_writeback_path.resolve()), "sha256": sha256(qc_writeback_path)},
            "qc_reject_membership": {"path": str(qc_reject_path.resolve()), "sha256": sha256(qc_reject_path)},
        },
        "fine_anchor_eligible": False,
    }
    manifest_path = args.out / "global_atlas_concordance_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 2 if review_rows else 0


if __name__ == "__main__":
    raise SystemExit(main())
