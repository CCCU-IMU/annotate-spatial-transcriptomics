#!/usr/bin/env python3
"""Legacy pool partitioner for migration-only projects."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path


def read(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def digest(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""): h.update(block)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--original-query", type=Path, required=True)
    ap.add_argument("--analyzed-membership", type=Path, required=True)
    ap.add_argument("--clusters", type=Path, required=True)
    ap.add_argument("--mapping", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cell-id-col", default="cell_id")
    ap.add_argument("--cluster-col", default="cluster")
    ap.add_argument("--generation", type=int, required=True)
    ap.add_argument("--source-run-id", required=True)
    args = ap.parse_args()
    original = read(args.original_query); analyzed = read(args.analyzed_membership); clusters = read(args.clusters); mapping = read(args.mapping)
    cc = args.cell_id_col
    required_map = {"record_type", "record_value", "target_pool", "provisional_broad_label", "provisional_fine_label", "state", "confidence", "state_tags_add", "next_route", "fine_anchor_eligible"}
    if not mapping or required_map - set(mapping[0]): raise SystemExit(f"mapping missing columns: {sorted(required_map - set(mapping[0] if mapping else []))}")
    original_by_id = {row[cc]: row for row in original}; analyzed_by_id = {row[cc]: row for row in analyzed}; cluster_by_id = {row[cc]: row[args.cluster_col] for row in clusters}
    if len(original_by_id) != len(original) or len(analyzed_by_id) != len(analyzed) or len(cluster_by_id) != len(clusters): raise SystemExit("duplicate cell IDs")
    qmap = {row["record_value"]: row for row in mapping if row["record_type"] == "query_cluster"}
    amap = {row["record_value"]: row for row in mapping if row["record_type"] == "promoted_anchor"}
    if set(cluster_by_id.values()) - set(qmap): raise SystemExit(f"unmapped query clusters: {sorted(set(cluster_by_id.values()) - set(qmap))}")
    output = []
    for cell_id, origin in original_by_id.items():
        if cell_id in cluster_by_id:
            key_type, key_value, rule = "query_cluster", cluster_by_id[cell_id], qmap[cluster_by_id[cell_id]]
        else:
            analyzed_row = analyzed_by_id.get(cell_id)
            if not analyzed_row or analyzed_row.get("query_or_anchor") != "anchor": raise SystemExit(f"original query cell lost from query and promoted-anchor sets: {cell_id}")
            key_type, key_value = "promoted_anchor", analyzed_row.get("anchor_label", "")
            if key_value not in amap: raise SystemExit(f"unmapped promoted anchor label: {key_value}")
            rule = amap[key_value]
        state_tags = ";".join(x for x in [origin.get("state_tags", ""), rule.get("state_tags_add", "")] if x).strip(";")
        output.append({
            "cell_id": cell_id,
            "source_key": origin.get("source_key", ""),
            "parent_decision_id": origin.get("parent_decision_id", origin.get("decision_id", "")),
            "parent_pool_id": origin.get("parent_pool_id", ""),
            "candidate_lineages": origin.get("candidate_lineages", ""),
            "state_tags": state_tags,
            "spatial_tags": origin.get("spatial_tags", ""),
            "qc_tags": origin.get("qc_tags", ""),
            "generation": args.generation,
            "source_run_id": args.source_run_id,
            "partition_record_type": key_type,
            "partition_record_value": key_value,
            "target_pool": rule["target_pool"],
            "provisional_broad_label": rule["provisional_broad_label"],
            "provisional_fine_label": rule["provisional_fine_label"],
            "state": rule["state"],
            "confidence": rule["confidence"],
            "next_route": rule["next_route"],
            "fine_anchor_eligible": rule["fine_anchor_eligible"],
        })
    if len(output) != len(original_by_id) or {x["cell_id"] for x in output} != set(original_by_id): raise SystemExit("original query conservation failure")
    fields = list(output[0]); write(args.out, output, fields)
    counts = Counter((row["target_pool"], row["state"], row["partition_record_type"], row["partition_record_value"]) for row in output)
    summary = [{"target_pool": k[0], "state": k[1], "record_type": k[2], "record_value": k[3], "n": v} for k, v in sorted(counts.items())]
    write(Path(str(args.out) + ".summary.tsv"), summary, list(summary[0]))
    result = {"status": "PASS", "n_original_query": len(original), "n_partitioned": len(output), "n_target_pools": len({x["target_pool"] for x in output}), "sha256": digest(args.out), "output": str(args.out.resolve())}
    Path(str(args.out) + ".manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2)); return 0


if __name__ == "__main__": raise SystemExit(main())
