#!/usr/bin/env python3
"""Inventory transcriptomics inputs and existing progress without mutating inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


KINDS = [
    (re.compile(r"\.h5ad$", re.I), "anndata"),
    (re.compile(r"\.(rds|rda|rdata)$", re.I), "r_object"),
    (re.compile(r"\.(h5|hdf5)$", re.I), "hdf5"),
    (re.compile(r"clusters?\.(csv|tsv)(\.gz)?$", re.I), "cluster_table"),
    (re.compile(r"(umap|pca|tsne).*\.(csv|tsv)(\.gz)?$", re.I), "reduction_table"),
    (re.compile(r"(spatial|coordinate|position).*\.(csv|tsv)(\.gz)?$", re.I), "spatial_table"),
    (re.compile(r"(summary|manifest|registry|ledger).*\.(csv|tsv|json|ya?ml)(\.gz)?$", re.I), "manifest_or_state"),
    (re.compile(r"(done|complete|success).*\.txt$", re.I), "completion_sentinel"),
    (re.compile(r"\.(png|pdf|svg|html)$", re.I), "report_or_figure"),
    (re.compile(r"\.(out|err|log)$", re.I), "log"),
    (re.compile(r"\.(csv|tsv)(\.gz)?$", re.I), "table"),
]


def classify(name: str) -> str:
    for pattern, kind in KINDS:
        if pattern.search(name):
            return kind
    return "other"


def sha256_small(path: Path, limit: int) -> str:
    if path.stat().st_size > limit:
        return "SKIPPED_LARGE"
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input_root", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--hash-max-mb", type=int, default=64)
    ap.add_argument("--max-depth", type=int, default=8)
    args = ap.parse_args()
    root = args.input_root.resolve()
    if not root.is_dir():
        raise SystemExit(f"Input root is not a directory: {root}")
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for base, dirs, files in os.walk(root):
        rel_depth = len(Path(base).relative_to(root).parts)
        if rel_depth >= args.max_depth:
            dirs[:] = []
        dirs[:] = [d for d in dirs if not d.startswith(".git")]
        for name in sorted(files):
            path = Path(base) / name
            stat = path.stat()
            rows.append({
                "path": str(path), "relative_path": str(path.relative_to(root)),
                "kind": classify(name), "size_bytes": stat.st_size,
                "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                "sha256": sha256_small(path, args.hash_max_mb * 1024 * 1024),
            })
    manifest = args.out / "input_manifest.tsv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["path"], delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    counts = Counter(r["kind"] for r in rows)
    bank_grid = [r for r in rows if re.search(r"resolution[_-]?grid[_-]?summary\.csv$", r["relative_path"], re.I)]
    cluster_tables = [r for r in rows if r["kind"] == "cluster_table"]
    summary = {
        "input_root": str(root), "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_files": len(rows), "kind_counts": dict(sorted(counts.items())),
        "banksy_grid_summaries": [r["path"] for r in bank_grid],
        "cluster_tables": len(cluster_tables),
        "warnings": (["No expression object detected"] if not any(r["kind"] in {"anndata", "r_object", "hdf5"} for r in rows) else []),
    }
    (args.out / "input_discovery.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

