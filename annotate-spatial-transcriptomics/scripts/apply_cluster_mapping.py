#!/usr/bin/env python3
"""Create a cell ledger by joining a frozen cluster table to reviewed mappings."""

from __future__ import annotations

import argparse
import csv
import gzip
from datetime import datetime, timezone
from pathlib import Path


def open_text(path: Path, mode="rt"):
    return gzip.open(path, mode, newline="", encoding="utf-8") if path.suffix == ".gz" else path.open(mode.replace("t", ""), newline="", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True, type=Path)
    ap.add_argument("--mapping", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--cell-col", default="cell_id")
    ap.add_argument("--cluster-col", default="cluster")
    ap.add_argument("--sample-col", default="sample_id")
    ap.add_argument("--source-method", required=True)
    ap.add_argument("--source-run", required=True)
    ap.add_argument("--decision-version", default="v001")
    args = ap.parse_args()
    with args.mapping.open(newline="", encoding="utf-8") as h:
        maps = list(csv.DictReader(h, delimiter="\t"))
    required = {"source_cluster", "broad_label", "fine_label", "state", "confidence", "evidence_status", "route", "fine_anchor_eligible", "next_action", "closed"}
    if not maps or not required.issubset(maps[0]): raise SystemExit(f"mapping lacks: {sorted(required - set(maps[0] if maps else []))}")
    by = {str(r["source_cluster"]): r for r in maps}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_name(args.out.stem + ".tmp" + args.out.suffix)
    fields = ["sample_id", "cell_id", "decision_id", "spatial_object_id", "analysis_scope", "input_snapshot_id", "source_method", "source_run_id", "source_cluster", "source_key", "parent_decision_id", "parent_pool_id", "pool_snapshot_id", "generation", "route", "route_attempt_id", "route_run_id", "iteration", "state_tags", "spatial_tags", "qc_tags", "candidate_lineages", "state", "broad_label", "fine_label", "confidence", "evidence_status", "validation_status", "validation_artifact", "validation_feature_scope", "fine_anchor_eligible", "decision_version", "supersedes", "closed", "closure_rationale", "next_action", "created_at"]
    now = datetime.now(timezone.utc).isoformat(); n = 0; missing = set(); seen = set()
    inp = gzip.open(args.clusters, "rt", newline="", encoding="utf-8") if args.clusters.suffix == ".gz" else args.clusters.open(newline="", encoding="utf-8")
    out = gzip.open(tmp, "wt", newline="", encoding="utf-8") if args.out.suffix == ".gz" else tmp.open("w", newline="", encoding="utf-8")
    with inp, out:
        reader = csv.DictReader(inp); writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t"); writer.writeheader()
        for row in reader:
            cid = str(row[args.cell_col]); cl = str(row[args.cluster_col]); key = (row.get(args.sample_col, ""), cid)
            if key in seen: raise SystemExit(f"duplicate cell: {key}")
            seen.add(key); m = by.get(cl)
            if m is None: missing.add(cl); continue
            writer.writerow({
                "sample_id": row.get(args.sample_col, ""), "cell_id": cid, "decision_id": m.get("decision_id") or f"{args.decision_version}:{args.source_run}:{cl}", "spatial_object_id": m.get("spatial_object_id", ""), "analysis_scope": m.get("analysis_scope",row.get("analysis_scope","analysis_set")), "input_snapshot_id": m.get("input_snapshot_id", row.get("input_snapshot_id","")),
                "source_method": args.source_method, "source_run_id": args.source_run, "source_cluster": cl,
                "source_key": m.get("source_key",row.get("source_key",f"{args.source_run}|{cl}")), "parent_decision_id": m.get("parent_decision_id",row.get("parent_decision_id","")), "parent_pool_id": m.get("parent_pool_id", row.get("parent_pool_id","whole_tissue")), "pool_snapshot_id": m.get("pool_snapshot_id",row.get("pool_snapshot_id","")), "generation": m.get("generation",row.get("generation",m.get("iteration","1"))), "route": m["route"], "route_attempt_id": m.get("route_attempt_id",row.get("route_attempt_id","")), "route_run_id": m.get("route_run_id", ""), "iteration": m.get("iteration", "1"), "state_tags": m.get("state_tags",row.get("state_tags","")), "spatial_tags": m.get("spatial_tags",row.get("spatial_tags","")), "qc_tags": m.get("qc_tags",row.get("qc_tags","")), "candidate_lineages": m.get("candidate_lineages",row.get("candidate_lineages","")), "state": m["state"],
                "broad_label": m["broad_label"], "fine_label": m["fine_label"], "confidence": m["confidence"],
                "evidence_status": m["evidence_status"], "validation_status": m.get("validation_status", "not_required"), "validation_artifact": m.get("validation_artifact", ""), "validation_feature_scope": m.get("validation_feature_scope", "unknown"), "fine_anchor_eligible": m["fine_anchor_eligible"],
                "decision_version": args.decision_version, "supersedes": m.get("supersedes", ""), "closed": m["closed"],
                "closure_rationale": m.get("closure_rationale", ""), "next_action": m["next_action"], "created_at": now,
            }); n += 1
    if missing:
        tmp.unlink(missing_ok=True); raise SystemExit(f"unmapped clusters: {sorted(missing)}")
    tmp.replace(args.out)
    print(f"wrote {n} observations to {args.out}")
    return 0


if __name__ == "__main__": raise SystemExit(main())
