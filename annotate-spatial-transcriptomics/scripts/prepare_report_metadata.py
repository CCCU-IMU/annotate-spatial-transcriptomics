­r‡^Ñf¥–Ø¦{¿lyÊ'vÃ®¶›­#!/usr/bin/env python3
"""Build atomic report metadata from the validated strict/inclusive/display cell ledger."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

import pandas as pd


RETAINED_LABELS = {
    "interface_review": "Unresolved interface",
    "qc_holdout": "QC holdout",
    "technical_state": "Technical retained",
    "pending_review": "Pending review",
    "excluded_initial_qc": "Excluded initial QC",
}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cell-ledger", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--cell-id-col", default="cell_id")
    p.add_argument("--expected-observations", type=int)
    a = p.parse_args()
    d = pd.read_csv(a.cell_ledger, sep="\t", dtype={a.cell_id_col: str}, keep_default_na=False)
    required = {
        a.cell_id_col, "analysis_scope", "strict_state", "strict_broad_label", "strict_fine_label",
        "inclusive_state", "inclusive_broad_label", "inclusive_fine_label",
        "display_state", "display_broad_label", "display_fine_label",
    }
    if not required.issubset(d):
        raise SystemExit(f"cell ledger lacks report-view fields: {sorted(required-set(d))}")
    if d[a.cell_id_col].duplicated().any():
        raise SystemExit("duplicate cell IDs in cell ledger")
    if a.expected_observations is not None and len(d) != a.expected_observations:
        raise SystemExit(f"expected {a.expected_observations} observations, found {len(d)}")

    fallback = d.display_state.map(RETAINED_LABELS).fillna("Unresolved")
    d["broad_display"] = d.display_broad_label.where(d.display_broad_label.ne(""), fallback)
    d["subtype_display"] = d.display_fine_label
    broad_only = d.display_state.eq("defined_broad_only") & d.display_broad_label.ne("") & d.subtype_display.eq("")
    d.loc[broad_only, "subtype_display"] = "Broad only: " + d.loc[broad_only, "display_broad_label"]
    missing_subtype = d.subtype_display.eq("")
    d.loc[missing_subtype, "subtype_display"] = fallback.loc[missing_subtype]
    d["strict_evidence_eligible"] = (
        d.analysis_scope.eq("analysis_set")
        & d.strict_state.isin({"defined_fine", "defined_broad_only"})
        & d.strict_broad_label.ne("")
    ).map({True: "true", False: "false"})
    d["inclusive_biological_eligible"] = (
        d.analysis_scope.eq("analysis_set")
        & d.inclusive_state.isin({"defined_fine", "defined_broad_only"})
        & d.inclusive_broad_label.ne("")
    ).map({True: "true", False: "false"})
    d["primary_broad_label"] = d.inclusive_broad_label.where(
        d.inclusive_biological_eligible.eq("true"), ""
    )
    d["primary_subtype_label"] = d.strict_fine_label.where(
        d.strict_state.eq("defined_fine") & d.strict_fine_label.ne(""), ""
    )
    d["strict_broad_evidence_label"] = d.strict_broad_label.where(
        d.strict_evidence_eligible.eq("true"), ""
    )

    a.out.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{a.out.name}.", suffix=".tmp", dir=a.out.parent)
    os.close(fd); tmp = Path(name)
    try:
        d.to_csv(tmp, sep="\t", index=False, compression="gzip" if str(a.out).endswith(".gz") else None)
        check = pd.read_csv(
            tmp, sep="\t", usecols=[a.cell_id_col], dtype={a.cell_id_col: str},
            compression="gzip" if str(a.out).endswith(".gz") else None,
        )
        if len(check) != len(d) or check[a.cell_id_col].nunique() != len(d):
            raise SystemExit("atomic report metadata validation failed")
        os.replace(tmp, a.out)
    finally:
        tmp.unlink(missing_ok=True)
    digest = hashlib.sha256(a.out.read_bytes()).hexdigest()
    result = {
        "status": "PASS", "n_observations": int(len(d)),
        "n_strict_evidence_eligible": int((d.strict_evidence_eligible == "true").sum()),
        "n_inclusive_biological_eligible": int((d.inclusive_biological_eligible == "true").sum()),
        "n_display_broad_labels": int(d.broad_display.nunique()),
        "n_display_subtype_labels": int(d.subtype_display.nunique()),
        "sha256": digest, "output": str(a.out.resolve()),
    }
    Path(str(a.out) + ".manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
