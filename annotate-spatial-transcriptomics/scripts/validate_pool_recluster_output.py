#!/usr/bin/env python3
"""Fail-closed validation for query-only pool reclustering artifacts."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path


def read_tsv(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def manifest(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    return {row.get("parameter", ""): row.get("value", "") for row in rows}


def truthy(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("output_root", type=Path)
    ap.add_argument("--expected-query", type=int, required=True)
    ap.add_argument("--expected-anchors", type=int, required=True)
    ap.add_argument("--resolutions", required=True)
    ap.add_argument("--require-spatial", action="store_true")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    root = args.output_root
    errors: list[str] = []
    required = [
        root / "RUN_COMPLETE.tsv",
        root / "run_manifest.tsv",
        root / "sessionInfo.txt",
        root / "pool_reclustered_query_seurat.rds",
        root / "joint_query_anchor_pca_seurat.rds",
        root / "tables/analyzed_membership.tsv.gz",
        root / "tables/query_anchor_distance_evidence.tsv",
        root / "tables/cluster_anchor_distance_summary.tsv",
        root / "tables/cluster_source_state_composition.tsv",
        root / "tables/cluster_QC_summary.tsv",
    ]
    for path in required:
        if not path.exists() or path.stat().st_size == 0:
            errors.append(f"missing_or_empty:{path.relative_to(root)}")
    mf = manifest(root / "run_manifest.tsv") if (root / "run_manifest.tsv").exists() else {}
    if not truthy(mf.get("anchor_assisted", "")): errors.append("manifest_anchor_assisted_not_true")
    if not truthy(mf.get("query_only_graph_umap_deg", "")): errors.append("manifest_query_only_not_true")
    if int(float(mf.get("n_query_analyzed", -1))) != args.expected_query: errors.append("manifest_query_count_mismatch")
    if int(float(mf.get("n_anchors_analyzed", -1))) != args.expected_anchors: errors.append("manifest_anchor_count_mismatch")
    expected_ids: set[str] | None = None
    evidence = root / "tables/query_anchor_distance_evidence.tsv"
    if evidence.exists():
        ev = read_tsv(evidence); expected_ids = {row.get("cell_id", "") for row in ev}
        if len(expected_ids) != args.expected_query: errors.append("anchor_evidence_query_count_mismatch")
    summary = []
    for resolution in [x.strip() for x in args.resolutions.split(",") if x.strip()]:
        tag = resolution.replace(".", "p")
        cp = root / f"tables/framework_res{tag}_clusters.tsv"
        dp = root / f"tables/framework_res{tag}_DEG_all.tsv"
        top = root / f"tables/framework_res{tag}_DEG_top100.tsv"
        expected_paths = [cp, dp, top, root / f"figures/framework_res{tag}_UMAP.png", root / f"figures/framework_res{tag}_UMAP.pdf"]
        if args.require_spatial:
            expected_paths.extend([root / f"figures/framework_res{tag}_spatial.png", root / f"figures/framework_res{tag}_spatial.pdf"])
        for path in expected_paths:
            if not path.exists() or path.stat().st_size == 0: errors.append(f"missing_or_empty:{path.relative_to(root)}")
        if cp.exists():
            rows = read_tsv(cp); ids = [row.get("cell_id", "") for row in rows]; clusters = {row.get("cluster", "") for row in rows}
            if len(ids) != args.expected_query or len(set(ids)) != args.expected_query: errors.append(f"cluster_membership_count_or_uniqueness:{resolution}")
            if expected_ids is not None and set(ids) != expected_ids: errors.append(f"cluster_membership_id_mismatch:{resolution}")
            if args.require_spatial:
                for cluster in clusters:
                    for ext in ["png", "pdf"]:
                        hp = root / f"figures/framework_res{tag}_cluster_{cluster}_highlight.{ext}"
                        if not hp.exists() or hp.stat().st_size == 0: errors.append(f"missing_highlight:{resolution}:{cluster}:{ext}")
            summary.append({"resolution": resolution, "n_clusters": len(clusters), "n_query": len(ids)})
    result = {"status": "PASS" if not errors else "FAIL", "output_root": str(root.resolve()), "expected_query": args.expected_query, "expected_anchors": args.expected_anchors, "resolutions": summary, "errors": errors}
    out = args.out or root / "artifact_validation.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
