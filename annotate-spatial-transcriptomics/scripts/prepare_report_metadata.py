#!/usr/bin/env python3
"""Build atomic report metadata from the validated single final annotation."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

import pandas as pd


RETAINED_LABELS = {
    "interface_review": "Anatomical interface",
    "qc_holdout": "QC holdout",
    "technical_state": "Technical retained",
    "pending_review": "Pending review",
    "excluded_initial_qc": "Excluded initial QC",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-ledger", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--expected-observations", type=int)
    args = parser.parse_args()
    data = pd.read_csv(args.cell_ledger, sep="\t", dtype={args.cell_id_col: str}, keep_default_na=False)
    required = {
        args.cell_id_col, "analysis_scope", "final_state", "final_broad_label", "final_fine_label",
        "final_confidence", "final_assignment_tier", "final_broad_eligible", "final_fine_eligible",
        "fine_anchor_eligible", "route",
    }
    if not required.issubset(data):
        raise SystemExit(f"cell ledger lacks final report fields: {sorted(required-set(data))}")
    if data[args.cell_id_col].duplicated().any():
        raise SystemExit("duplicate cell IDs in cell ledger")
    if args.expected_observations is not None and len(data) != args.expected_observations:
        raise SystemExit(f"expected {args.expected_observations} observations, found {len(data)}")
    fallback = data.final_state.map(RETAINED_LABELS).fillna("Unresolved")
    # Release figures and evidence tables must use the primary_* columns below.
    # Keep retained/QC states in their own display column so a broad-only return
    # can never be fabricated into a subtype such as "Broad only: Granulosa".
    data["retained_state_display"] = fallback.where(data.final_broad_label.eq(""), "")
    data["broad_display"] = data.final_broad_label.where(data.final_broad_label.ne(""), fallback)
    data["subtype_display"] = data.final_fine_label
    broad_eligible = data.final_broad_eligible.astype(str).str.lower().isin({"true", "1", "yes"})
    fine_eligible = data.final_fine_eligible.astype(str).str.lower().isin({"true", "1", "yes"})
    data["primary_broad_label"] = data.final_broad_label.where(broad_eligible, "")
    data["primary_subtype_label"] = data.final_fine_label.where(fine_eligible, "")
    data["fine_marker_discovery_eligible"] = data.final_fine_eligible
    anchor_confidence = data["final_fine_confidence"] if "final_fine_confidence" in data else data.final_confidence
    data["anchor_reference_eligible"] = (
        anchor_confidence.eq("high")
        & data.fine_anchor_eligible.astype(str).str.lower().isin({"true", "1", "yes"})
    ).map({True: "true", False: "false"})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{args.out.name}.", suffix=".tmp", dir=args.out.parent)
    os.close(fd); temp = Path(name)
    try:
        data.to_csv(temp, sep="\t", index=False, compression="gzip" if str(args.out).endswith(".gz") else None)
        check = pd.read_csv(temp, sep="\t", usecols=[args.cell_id_col], dtype={args.cell_id_col: str}, compression="gzip" if str(args.out).endswith(".gz") else None)
        if len(check) != len(data) or check[args.cell_id_col].nunique() != len(data):
            raise SystemExit("atomic report metadata validation failed")
        os.replace(temp, args.out)
    finally:
        temp.unlink(missing_ok=True)
    result = {
        "status": "PASS", "n_observations": int(len(data)),
        "n_final_broad": int(broad_eligible.sum()),
        "n_final_fine": int(fine_eligible.sum()),
        "n_anchor_reference_eligible": int((data.anchor_reference_eligible == "true").sum()),
        "n_broad_labels": int(data.primary_broad_label[data.primary_broad_label.ne("")].nunique()),
        "n_fine_labels": int(data.primary_subtype_label[data.primary_subtype_label.ne("")].nunique()),
        "release_label_columns": {
            "broad": "primary_broad_label",
            "subtype": "primary_subtype_label",
            "retained_state": "retained_state_display",
        },
        "sha256": hashlib.sha256(args.out.read_bytes()).hexdigest(), "output": str(args.out.resolve()),
    }
    Path(str(args.out) + ".manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
