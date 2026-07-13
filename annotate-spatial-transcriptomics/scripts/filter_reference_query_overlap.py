#!/usr/bin/env python3
"""Create an immutable reference membership after excluding current-query observations."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

import pandas as pd


def read_table(path: Path, cell_id_col: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={cell_id_col: str}, keep_default_na=False)


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        compression = "gzip" if str(path).endswith(".gz") else None
        frame.to_csv(tmp, sep="\t", index=False, compression=compression)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--query", required=True, type=Path)
    p.add_argument("--reference", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--cell-id-col", default="cell_id")
    p.add_argument("--role-col", default="query_or_anchor")
    p.add_argument("--reference-role", default="anchor")
    p.add_argument("--label-col", default="anchor_label")
    p.add_argument("--min-labels", type=int, default=2)
    p.add_argument("--min-per-label", type=int, default=20)
    a = p.parse_args()

    query = read_table(a.query, a.cell_id_col)
    reference = read_table(a.reference, a.cell_id_col)
    required = {a.cell_id_col, a.role_col, a.label_col}
    if a.cell_id_col not in query or not required.issubset(reference):
        raise SystemExit("query/reference membership columns are incomplete")
    if query[a.cell_id_col].duplicated().any() or reference[a.cell_id_col].duplicated().any():
        raise SystemExit("duplicate observation IDs in query/reference membership")

    anchors = reference[reference[a.role_col].astype(str).eq(a.reference_role)].copy()
    before = anchors[a.label_col].astype(str).value_counts().sort_index()
    overlap = anchors[a.cell_id_col].isin(set(query[a.cell_id_col]))
    clean = anchors.loc[~overlap].copy()
    after = clean[a.label_col].astype(str).value_counts().sort_index()
    if clean[a.cell_id_col].isin(set(query[a.cell_id_col])).any():
        raise SystemExit("reference/query overlap remains after filtering")
    if len(after) < a.min_labels or (after < a.min_per_label).any():
        raise SystemExit(f"reference labels insufficient after overlap removal: {after.to_dict()}")

    atomic_tsv(clean, a.out)
    digest = hashlib.sha256(a.out.read_bytes()).hexdigest()
    labels = sorted(set(before.index) | set(after.index))
    summary = pd.DataFrame({
        "reference_label": labels,
        "n_before": [int(before.get(x, 0)) for x in labels],
        "n_excluded_overlap": [int(before.get(x, 0) - after.get(x, 0)) for x in labels],
        "n_after": [int(after.get(x, 0)) for x in labels],
    })
    summary_path = Path(str(a.out) + ".summary.tsv")
    atomic_tsv(summary, summary_path)
    manifest = {
        "status": "PASS",
        "n_query": int(len(query)),
        "n_reference_candidates": int(len(anchors)),
        "n_excluded_overlap": int(overlap.sum()),
        "n_reference_final": int(len(clean)),
        "n_reference_labels": int(len(after)),
        "query_reference_overlap_final": 0,
        "sha256": digest,
        "output": str(a.out.resolve()),
    }
    Path(str(a.out) + ".manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
