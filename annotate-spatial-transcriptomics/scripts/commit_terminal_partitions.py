#!/usr/bin/env python3
"""Atomically supersede current cell decisions from non-overlapping terminal pool partitions."""
from __future__ import annotations

import argparse
import fcntl
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def atomic_frame(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(fd); tmp = Path(name)
    try:
        if str(path).endswith(".gz"):
            with gzip.open(tmp, "wt") as handle:
                frame.to_csv(handle, sep="\t", index=False)
        else:
            frame.to_csv(tmp, sep="\t", index=False)
        os.replace(tmp, path)
    finally:
        if tmp.exists(): tmp.unlink()


def truth(x: pd.Series) -> pd.Series:
    return x.astype(str).str.lower().isin({"true", "1", "yes"})


def choose(frame: pd.DataFrame, *cols: str, default: str = "") -> pd.Series:
    out = pd.Series(default, index=frame.index, dtype=object)
    for col in cols:
        if col in frame:
            value = frame[col].astype(str)
            take = out.astype(str).str.len().eq(0) & value.str.len().gt(0)
            out.loc[take] = value.loc[take]
    return out


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("project_root", type=Path)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--decision-version", required=True)
    p.add_argument("--iteration", required=True, type=int)
    p.add_argument("--sample", required=True)
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args(); root = a.project_root.resolve()
    cell_path = root / "state/cell_ledger.tsv.gz"; cluster_path = root / "state/cluster_decision_ledger.tsv"
    spec = read_tsv(a.manifest)
    required = {"partition_path", "route_id", "validation_artifact", "closure_rationale"}
    if not required.issubset(spec): raise SystemExit(f"manifest missing {sorted(required-set(spec))}")
    pieces = []
    for row in spec.to_dict("records"):
        path = Path(row["partition_path"]); path = path if path.is_absolute() else root / path
        d = read_tsv(path)
        if "cell_id" not in d or "target_pool" not in d: raise SystemExit(f"invalid terminal partition: {path}")
        include = {x for x in row.get("include_target_pools", "").split("|") if x}
        exclude = {x for x in row.get("exclude_target_pools", "").split("|") if x}
        if include: d = d[d.target_pool.isin(include)].copy()
        if exclude: d = d[~d.target_pool.isin(exclude)].copy()
        d["_route_id"] = row["route_id"]; d["_validation_artifact"] = row["validation_artifact"]; d["_closure_rationale"] = row["closure_rationale"]; d["_partition_path"] = str(path.resolve())
        pieces.append(d)
    terminal = pd.concat(pieces, ignore_index=True, sort=False).fillna("")
    if terminal.empty or terminal.cell_id.duplicated().any():
        dup = terminal.loc[terminal.cell_id.duplicated(False), "cell_id"].head().tolist()
        raise SystemExit(f"empty or overlapping terminal partitions: {dup}")
    terminal["_state"] = choose(terminal, "state", "final_state")
    terminal["_broad"] = choose(terminal, "provisional_broad_label", "final_broad_label")
    terminal["_broad"] = terminal["_broad"].replace({"Granulosa": "Follicular somatic", "Vascular-associated": "Stromal/vascular-associated"})
    terminal["_fine"] = choose(terminal, "provisional_fine_label", "final_fine_label")
    terminal["_confidence"] = choose(terminal, "confidence", "final_confidence", default="low")
    terminal["_source_run"] = choose(terminal, "source_run_id", default="terminal_partition")
    terminal["_source_cluster"] = choose(terminal, "partition_record_value", "cluster", default="adjudicated")
    terminal["_fine_anchor"] = choose(terminal, "fine_anchor_eligible", "fine_anchor_eligible_override", default="false")
    terminal["_next_action"] = choose(terminal, "next_route", default="none")
    terminal["_evidence"] = choose(terminal, "state_tags", "partition_action", "final_action", default="validated terminal multiroute partition")
    terminal.loc[terminal._state.isin(["interface_review", "qc_holdout", "pending_review"]), "_next_action"] = "none"
    route_text = terminal["_route_id"].astype(str).str.lower() + " " + terminal["_evidence"].astype(str).str.lower()
    atlas_rows = route_text.str.contains("atlas|calibrated.*mapping", regex=True) & terminal["_state"].eq("defined_broad_only")
    if atlas_rows.any():
        tier_column = "mapping_tier" if "mapping_tier" in terminal else "consensus_tier" if "consensus_tier" in terminal else ""
        if not tier_column:
            raise SystemExit("atlas/mapping broad returns lack observation-level tier proof")
        tier = terminal[tier_column].astype(str)
        accepted = tier.isin(["high", "moderate", "moderate-only"])
        moderate_gate = choose(terminal, "meets_moderate_or_higher", "writeback_eligible", default="false").astype(str).str.lower().isin({"true", "1", "yes"})
        fine_gate = terminal["_fine_anchor"].astype(str).str.lower().isin({"false", "0", "no"})
        if (atlas_rows & (~accepted | ~moderate_gate | ~fine_gate)).any():
            raise SystemExit("atlas/mapping broad returns fail tier or broad-only gates")
    rctd_rows = route_text.str.contains("rctd", regex=False) & terminal["_state"].isin(["defined_fine", "defined_broad_only"])
    if rctd_rows.any():
        if "rctd_tier" not in terminal:
            raise SystemExit("RCTD-assisted returns lack observation-level confidence tier")
        tier = terminal["rctd_tier"].astype(str).str.lower()
        bad_fine = terminal["_state"].eq("defined_fine") & (~tier.eq("extreme") | ~choose(terminal, "independent_fine_evidence", default="false").astype(str).str.lower().isin({"true", "1", "yes"}))
        bad_broad = terminal["_state"].eq("defined_broad_only") & ~tier.isin({"extreme", "high"})
        if (rctd_rows & (bad_fine | bad_broad)).any():
            raise SystemExit("RCTD return violates extreme/fine or high/broad-only gate")
    bad = terminal._state.eq("") | ((terminal._state.isin(["defined_fine", "defined_broad_only"])) & terminal._broad.eq(""))
    if bad.any(): raise SystemExit(f"invalid terminal labels for {int(bad.sum())} observations")

    lock_path = root / "state/terminal_writeback.lock"; lock_path.touch(exist_ok=True)
    with lock_path.open("w") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        cells = read_tsv(cell_path); decisions = read_tsv(cluster_path)
        if cells.cell_id.duplicated().any() or not set(terminal.cell_id).issubset(set(cells.cell_id)):
            raise SystemExit("cell ledger boundary is invalid")
        old_n = len(cells); idx = cells.reset_index().set_index("cell_id").loc[terminal.cell_id, "index"].to_numpy(); old = cells.iloc[idx].copy().reset_index(drop=True)
        terminal = terminal.reset_index(drop=True)
        group_cols = ["_route_id", "_source_run", "_source_cluster", "target_pool", "_state", "_broad", "_fine", "_confidence", "_validation_artifact", "_closure_rationale", "_partition_path"]
        terminal["_decision_key"] = terminal[group_cols].astype(str).agg("|".join, axis=1)
        terminal["_decision_id"] = terminal._decision_key.map(lambda x: f"{a.decision_version}:{hashlib.sha256(x.encode()).hexdigest()[:16]}")
        now = datetime.now(timezone.utc).isoformat()
        replacement = old.copy()
        replacement["source_run_id"] = terminal._source_run.values; replacement["source_cluster"] = terminal._source_cluster.values; replacement["parent_pool_id"] = terminal.target_pool.values; replacement["route"] = terminal._route_id.values
        replacement["state"] = terminal._state.values; replacement["broad_label"] = terminal._broad.values; replacement["fine_label"] = terminal._fine.values; replacement["confidence"] = terminal._confidence.values
        replacement["evidence_status"] = terminal._evidence.values; replacement["fine_anchor_eligible"] = terminal._fine_anchor.values; replacement["decision_version"] = a.decision_version; replacement["supersedes"] = old.decision_id.values
        replacement["closed"] = "True"; replacement["next_action"] = terminal._next_action.values; replacement["iteration"] = str(a.iteration)
        replacement["validation_status"] = terminal._state.map(lambda x: "passed_retained" if x in {"interface_review", "qc_holdout", "pending_review"} else "passed").values
        replacement["closure_rationale"] = terminal._closure_rationale.values; replacement["validation_artifact"] = terminal._validation_artifact.values; replacement["validation_feature_scope"] = "full_feature"
        replacement["decision_id"] = terminal._decision_id.values
        cells.iloc[idx] = replacement[cells.columns].values
        if len(cells) != old_n or cells.cell_id.duplicated().any(): raise SystemExit("cell ledger conservation failed")

        new_rows = []
        for decision_id, g in terminal.groupby("_decision_id", sort=False):
            positions = g.index.to_numpy(); prior = old.iloc[positions]
            first = g.iloc[0]
            new_rows.append({
                "decision_id": decision_id, "decision_version": a.decision_version, "sample_id": a.sample,
                "source_run_id": first["_source_run"], "source_cluster": first["_source_cluster"], "n_observations": str(len(g)),
                "broad_label": first["_broad"], "fine_label": first["_fine"], "state": first["_state"], "confidence": first["_confidence"],
                "evidence_status": first["_evidence"], "route": first["_route_id"], "target_pool": first["target_pool"],
                "fine_anchor_eligible": first["_fine_anchor"], "next_action": first["_next_action"],
                "supersedes": ";".join(sorted(set(prior.decision_id.astype(str)))), "closed": "true",
                "closure_rationale": first["_closure_rationale"],
                "validation_status": "passed_retained" if first["_state"] in {"interface_review", "qc_holdout", "pending_review"} else "passed",
                "validation_artifact": first["_validation_artifact"], "iteration": str(a.iteration), "created_at": now,
                "validation_feature_scope": "full_feature", "spatial_object_count": "", "count_interpretation": "observations_not_inferred_cells"
            })
        add = pd.DataFrame(new_rows).reindex(columns=decisions.columns, fill_value="")
        decisions = pd.concat([decisions[~decisions.decision_id.isin(add.decision_id)], add], ignore_index=True)
        history = root / "state/history"; history.mkdir(parents=True, exist_ok=True); stamp = now.replace(":", "").replace("+00:00", "Z")
        backup_cell = history / f"cell_ledger.pre_{a.decision_version}.{stamp}.tsv.gz"; backup_decisions = history / f"cluster_decision_ledger.pre_{a.decision_version}.{stamp}.tsv"
        if not a.dry_run:
            shutil.copy2(cell_path, backup_cell); shutil.copy2(cluster_path, backup_decisions)
            atomic_frame(cells, cell_path); atomic_frame(decisions, cluster_path)
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
    summary = terminal.groupby(["_state", "_broad", "_fine"], dropna=False).size().rename("n").reset_index()
    outdir = root / "provenance/terminal_writeback"; outdir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(outdir / f"{a.decision_version}_writeback_summary.tsv", sep="\t", index=False)
    report = {"status": "DRY_RUN_PASS" if a.dry_run else "PASS", "decision_version": a.decision_version, "n_updated": len(terminal), "n_ledger_total": old_n, "n_new_decisions": len(new_rows), "partitions": spec.partition_path.tolist(), "cell_backup": str(backup_cell), "decision_backup": str(backup_decisions)}
    (outdir / f"{a.decision_version}_writeback_manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2)); return 0


if __name__ == "__main__": raise SystemExit(main())
