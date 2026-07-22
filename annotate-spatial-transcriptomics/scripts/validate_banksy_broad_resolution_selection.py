#!/usr/bin/env python3
"""Validate that a BANKSY whole-tissue resolution was chosen for broad biology."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(path: Path, annotation_contract: Path | None = None) -> dict[str, object]:
    errors: list[str] = []
    payload = json.loads(path.read_text(encoding="utf-8"))
    candidates = payload.get("candidates") or []
    selected = str(payload.get("selected_resolution", ""))
    if payload.get("method") != "BANKSY" or payload.get("question_mode") != "whole_tissue_broad_annotation":
        errors.append("selection must be a BANKSY whole-tissue broad-annotation decision")
    if not candidates or selected not in {str(x.get("resolution", "")) for x in candidates}:
        errors.append("selected resolution is absent from the complete candidate grid")
    contract = {}
    if annotation_contract:
        contract = json.loads(annotation_contract.read_text(encoding="utf-8"))
        whole = contract.get("whole_tissue_partition", {})
        if whole.get("method") != "BANKSY" or whole.get("candidate_grid_source") != "bound_upstream_input":
            errors.append("v2 BANKSY decision is not bound to an upstream-input resolution grid")
        expected = {float(value) for value in whole.get("candidate_resolutions", [])}
        try:
            observed = {float(row.get("resolution")) for row in candidates}
        except (TypeError, ValueError):
            observed = set()
        if observed != expected:
            errors.append("BANKSY decision candidates differ from the v2 annotation contract")
    required = {
        "full_catalog_lineage_scan", "default_broad_recall_review", "large_cluster_purity_review",
        "zero_census_review", "deg_marker_coherence_review", "spatial_morphology_review",
        "adjacent_resolution_migration_review", "technical_fragmentation_review",
    }
    for row in candidates:
        missing = [key for key in required if row.get(key) not in {True, "true", "pass", "passed", "reviewed"}]
        if missing:
            errors.append(f"resolution {row.get('resolution')}: missing reviews {missing}")
        if row.get("selection_basis") == "cluster_count_only":
            errors.append(f"resolution {row.get('resolution')}: cluster count cannot be the biological selection basis")
    selected_row = next((x for x in candidates if str(x.get("resolution", "")) == selected), {})
    if selected_row and not selected_row.get("selection_rationale"):
        errors.append("selected resolution lacks a biological rationale")
    carry = Path(str(payload.get("carryforward_signal_artifact", "")))
    if not carry.is_file() or sha256(carry) != payload.get("carryforward_signal_sha256"):
        errors.append("weak-signal carry-forward artifact is missing or stale")
    return {"status": "PASS" if not errors else "BLOCKED", "schema_version": "2.0" if annotation_contract else "1.x",
            "decision": str(path.resolve()), "decision_sha256": sha256(path),
            "annotation_contract": str(annotation_contract.resolve()) if annotation_contract else "",
            "annotation_contract_sha256": sha256(annotation_contract) if annotation_contract else "",
            "selected_resolution": selected, "n_candidates": len(candidates), "errors": errors}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("decision", type=Path)
    parser.add_argument("--annotation-contract", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = validate(args.decision.resolve(), args.annotation_contract.resolve() if args.annotation_contract else None)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
