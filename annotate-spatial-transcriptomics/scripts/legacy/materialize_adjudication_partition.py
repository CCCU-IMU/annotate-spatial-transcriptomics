#!/usr/bin/env python3
"""Legacy pool-partition writeback retained only for migration audits."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

import pandas as pd


def read_table(path: Path, cell_id_col: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={cell_id_col: str})


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".tsv.gz" if str(path).endswith(".gz") else ".tsv"
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=suffix, dir=path.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        frame.to_csv(tmp, sep="\t", index=False, compression="gzip" if str(path).endswith(".gz") else None)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--original-query", required=True, type=Path)
    p.add_argument("--adjudication", required=True, type=Path)
    p.add_argument("--mapping", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--cell-id-col", default="cell_id")
    p.add_argument("--action-col", default="final_action")
    p.add_argument("--state-col", default="final_state")
    p.add_argument("--broad-label-col", default="final_broad_label")
    p.add_argument("--fine-label-col", default="final_fine_label")
    p.add_argument("--confidence-col", default="final_confidence")
    p.add_argument("--generation", required=True, type=int)
    p.add_argument("--source-run-id", required=True)
    a = p.parse_args()

    q = read_table(a.original_query, a.cell_id_col)
    d = read_table(a.adjudication, a.cell_id_col)
    m = pd.read_csv(a.mapping, sep="\t", dtype=str).fillna("")
    if q[a.cell_id_col].duplicated().any() or d[a.cell_id_col].duplicated().any():
        raise SystemExit("duplicate cell IDs")
    if set(q[a.cell_id_col]) != set(d[a.cell_id_col]):
        raise SystemExit("adjudication does not conserve the original query membership")
    if a.action_col not in d or "action" not in m or "target_pool" not in m:
        raise SystemExit("missing action mapping columns")
    if m.action.duplicated().any():
        raise SystemExit("duplicate action mapping")
    unknown = sorted(set(d[a.action_col].astype(str)) - set(m.action))
    if unknown:
        raise SystemExit(f"unmapped adjudication actions: {unknown}")

    legacy_status = d.get("mapping_status", pd.Series("", index=d.index)).astype(str).str.contains("legacy|medium_high", case=False, regex=True)
    if legacy_status.any():
        raise SystemExit("legacy or combined mapping-tier adjudication is not writeback eligible")
    defined = d.get(a.state_col, pd.Series("", index=d.index)).astype(str).eq("defined_broad_only")
    mapping_like = any(column in d for column in ["mapping_tier", "consensus_tier", "meets_moderate_or_higher"])
    if mapping_like and defined.any():
        tier = d.get("mapping_tier", d.get("consensus_tier", pd.Series("", index=d.index))).astype(str)
        allowed = tier.isin(["high", "moderate_only"])
        moderate = d.get("meets_moderate_or_higher", d.get("writeback_eligible", pd.Series(False, index=d.index))).astype(str).str.lower().isin(["true", "1", "yes"])
        fine = d.get("fine_anchor_eligible", pd.Series(False, index=d.index)).astype(str).str.lower().isin(["true", "1", "yes"])
        validated = d.get("validated_broad_return", d.get("writeback_eligible", pd.Series(False, index=d.index))).astype(str).str.lower().isin(["true", "1", "yes"])
        if (defined & (~allowed | ~moderate | fine | ~validated)).any():
            raise SystemExit("one or more broad rescue rows fail tier, independent validation, or fine-anchor gates")

    keep = [a.cell_id_col, a.action_col]
    for col in [a.state_col, a.broad_label_col, a.fine_label_col, a.confidence_col, "fine_anchor_eligible"]:
        if col in d and col not in keep:
            keep.append(col)
    evidence_cols = [c for c in d.columns if c.endswith("_consistent") or c in {"validated_broad_return", "post_rctd_validated", "strict_oocyte_program", "plasma_review_candidate", "predicted_label_canonical", "proposed_canonical", "cluster"}]
    keep.extend(c for c in evidence_cols if c not in keep)
    z = q.merge(d[keep], on=a.cell_id_col, validate="one_to_one", sort=False)
    z = z.merge(m, left_on=a.action_col, right_on="action", validate="many_to_one", sort=False)
    z["source_run_id"] = a.source_run_id
    z["generation"] = a.generation
    z["partition_action"] = z[a.action_col].astype(str)
    z["state"] = z.get("state_override", "").where(z.get("state_override", "").astype(str).str.len().gt(0), z.get(a.state_col, "")) if "state_override" in z else z.get(a.state_col, "")
    z["provisional_broad_label"] = z.get("broad_label_override", "").where(z.get("broad_label_override", "").astype(str).str.len().gt(0), z.get(a.broad_label_col, "")) if "broad_label_override" in z else z.get(a.broad_label_col, "")
    z["provisional_fine_label"] = z.get("fine_label_override", "").where(z.get("fine_label_override", "").astype(str).str.len().gt(0), z.get(a.fine_label_col, "")) if "fine_label_override" in z else z.get(a.fine_label_col, "")
    z["confidence"] = z.get("confidence_override", "").where(z.get("confidence_override", "").astype(str).str.len().gt(0), z.get(a.confidence_col, "")) if "confidence_override" in z else z.get(a.confidence_col, "")
    if "fine_anchor_eligible_override" in z:
        use = z.fine_anchor_eligible_override.astype(str).str.len().gt(0)
        z.loc[use, "fine_anchor_eligible"] = z.loc[use, "fine_anchor_eligible_override"]
    if "next_route" not in z:
        z["next_route"] = ""

    atomic_tsv(z, a.out)
    digest = hashlib.sha256(a.out.read_bytes()).hexdigest()
    summary = z.groupby(["target_pool", "state", "provisional_broad_label", "partition_action"], dropna=False).size().rename("n").reset_index().sort_values("n", ascending=False)
    summary_path = Path(str(a.out) + ".summary.tsv")
    atomic_tsv(summary, summary_path)
    manifest = {"status": "PASS", "n_original_query": len(q), "n_partitioned": len(z), "n_target_pools": int(z.target_pool.nunique()), "sha256": digest, "source_run_id": a.source_run_id, "generation": a.generation, "output": str(a.out.resolve())}
    manifest_path = Path(str(a.out) + ".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
