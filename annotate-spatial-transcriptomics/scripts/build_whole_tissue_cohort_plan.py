#!/usr/bin/env python3
"""Build provisional one-initial-cluster/one-cohort boundaries without labels."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

from evidence_schema_lib import sha256


FORBIDDEN_RELEASE_FIELDS = {
    "broad_label", "fine_label", "final_broad_label", "final_fine_label",
    "initial_broad_label", "state", "final_state", "qc_reason",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = {str(value).lower() for value in (reader.fieldnames or [])}
        leaked = fields.intersection(FORBIDDEN_RELEASE_FIELDS)
        if leaked:
            raise SystemExit(
                "provisional input contains release fields: " + ", ".join(sorted(leaked))
            )
        return list(reader)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def number(value: object, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def safe(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "cluster"


def evidence_score(row: dict[str, str]) -> float:
    return (
        0.25 * number(row.get("observation_coherent_fraction"))
        + 0.20 * number(row.get("positive_marker_detection_fraction"))
        + 0.15 * number(row.get("spatial_local_support_fraction"))
        + 0.15 * number(row.get("cross_resolution_stable_fraction"))
        + 0.15 * max(0.0, number(row.get("marker_deg_log2fc_mean")))
        + 0.10 * max(0.0, number(row.get("mean_program_score")))
        - 0.25 * number(row.get("hard_contradiction_fraction"))
        - 0.10 * max(0.0, number(row.get("anti_marker_deg_log2fc_mean")))
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--partitions", required=True, type=Path)
    ap.add_argument("--cluster-evidence", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    partitions = read_tsv(args.partitions)
    selected = [row for row in partitions if row.get("resolution_role") == "selected"]
    if not selected or any(not row.get("cell_id") or not row.get("cluster") for row in selected):
        raise SystemExit("selected partition is empty or incomplete")
    cell_ids = [row["cell_id"] for row in selected]
    if len(cell_ids) != len(set(cell_ids)):
        raise SystemExit("selected partition contains duplicate observations")

    catalog = json.loads(args.catalog.read_text(encoding="utf-8"))
    catalog_rows = {
        str(row.get("candidate_id", "")): row
        for row in catalog.get("candidate_boundaries", [])
    }
    evidence = read_tsv(args.cluster_evidence)
    selected_resolution = str(selected[0].get("resolution", ""))
    by_cluster: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in evidence:
        if row.get("resolution_role") != "selected":
            continue
        candidate = catalog_rows.get(str(row.get("candidate_id", "")), {})
        release = str(candidate.get("release_broad_label", ""))
        role = str(candidate.get("candidate_role", ""))
        if not release or role in {"fine", "state", "exploratory"}:
            continue
        by_cluster[str(row.get("source_cluster", ""))].append({
            "candidate_id": str(row.get("candidate_id", "")),
            "release_broad_label": release,
            "score": evidence_score(row),
            "coherent_fraction": number(row.get("observation_coherent_fraction")),
            "stable_fraction": number(row.get("cross_resolution_stable_fraction")),
            "contradiction_fraction": number(row.get("hard_contradiction_fraction")),
        })

    members: dict[str, list[str]] = defaultdict(list)
    for row in selected:
        members[str(row["cluster"])].append(str(row["cell_id"]))
    plan_rows: list[dict[str, object]] = []
    all_membership: list[dict[str, object]] = []
    membership_dir = args.out / "memberships"
    for cluster in sorted(members, key=lambda value: (number(value, float("inf")), value)):
        cohort_id = f"initial_cluster__{safe(cluster)}"
        candidates = sorted(
            by_cluster.get(cluster, []),
            key=lambda row: (-number(row["score"]), str(row["candidate_id"])),
        )
        positive = [
            row for row in candidates
            if number(row["coherent_fraction"]) >= 0.02
            or number(row["stable_fraction"]) >= 0.02
            or number(row["score"]) > 0
        ]
        winner = positive[0] if positive else {}
        runner = positive[1] if len(positive) > 1 else {}
        close_competitor = bool(
            runner and (
                number(winner.get("score")) - number(runner.get("score")) <= 0.10
                or number(runner.get("coherent_fraction")) >= 0.05
            )
        )
        if not winner:
            status = "unknown"
        elif close_competitor:
            status = "mixed"
        else:
            status = "provisional_broad"
        membership_rows = [{"cell_id": cell} for cell in sorted(members[cluster])]
        membership_path = membership_dir / f"{cohort_id}.tsv.gz"
        write_tsv(membership_path, membership_rows, ["cell_id"])
        for cell in sorted(members[cluster]):
            all_membership.append({
                "cell_id": cell,
                "source_initial_cluster": cluster,
                "cohort_id": cohort_id,
            })
        plan_rows.append({
            "cohort_id": cohort_id,
            "source_initial_cluster": cluster,
            "selected_resolution": selected_resolution,
            "n_observations": len(members[cluster]),
            "provisional_status": status,
            "provisional_broad_after_score_freeze": winner.get("release_broad_label", ""),
            "leading_candidate_id": winner.get("candidate_id", ""),
            "leading_score": winner.get("score", ""),
            "runner_up_candidate_id": runner.get("candidate_id", ""),
            "runner_up_score": runner.get("score", ""),
            "watch_candidate_ids": ";".join(str(row["candidate_id"]) for row in positive[1:]),
            "membership_path": str(membership_path.resolve()),
            "membership_sha256": sha256(membership_path),
            "formal_label_written": "false",
        })

    args.out.mkdir(parents=True, exist_ok=True)
    plan_path = args.out / "whole_tissue_cohort_plan.tsv"
    write_tsv(plan_path, plan_rows, list(plan_rows[0]))
    map_path = args.out / "whole_tissue_cluster_membership.tsv.gz"
    write_tsv(map_path, all_membership, ["cell_id", "source_initial_cluster", "cohort_id"])
    if len(all_membership) != len(cell_ids) or {row["cell_id"] for row in all_membership} != set(cell_ids):
        raise SystemExit("cohort plan does not exactly partition the selected analysis set")
    manifest = {
        "status": "PASS",
        "schema_version": "2.2",
        "stage": "whole_tissue_partition",
        "formal_membership_written": False,
        "provisional_only": True,
        "n_observations": len(all_membership),
        "n_initial_clusters": len(plan_rows),
        "cohort_plan": {"path": str(plan_path.resolve()), "sha256": sha256(plan_path)},
        "cluster_membership": {"path": str(map_path.resolve()), "sha256": sha256(map_path)},
    }
    (args.out / "whole_tissue_cohort_plan_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
