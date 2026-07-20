#!/usr/bin/env python3
"""Validate continuous open-world lineage scans at every annotation boundary."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


TRUE = {"1", "true", "yes", "y", "pass", "passed"}
TERMINAL_COHORT = {"completed", "terminal", "validated_done", "closed", "not_applicable_reviewed"}
SIGNAL_STATES = {"absent", "watch", "candidate", "supported", "refuted"}
RESOLVED_OUTCOMES = {
    "supported_label", "direct_return", "targeted_cohort_opened", "parent_retained_after_refutation",
    "qc_or_unknown", "refuted_by_multichannel_evidence", "merged_as_state_only",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split(value: str) -> list[str]:
    return [x for x in re.split(r"[;,\s]+", str(value or "").strip()) if x]


def active(rows: list[dict[str, str]], id_field: str) -> list[dict[str, str]]:
    superseded = {x for row in rows for x in split(row.get("supersedes", ""))}
    return [row for row in rows if row.get(id_field, "") not in superseded]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def number(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def ordered_resolutions(value: str) -> list[str]:
    values = split(value)
    return [x for _, x in sorted((number(x, float("inf")), x) for x in values)]


def required_resolutions(boundary: dict[str, str]) -> list[str]:
    candidates = ordered_resolutions(boundary.get("candidate_resolutions", ""))
    selected = boundary.get("selected_resolution", "")
    if not selected or selected not in candidates:
        return [selected] if selected else []
    higher = [value for value in candidates if number(value) > number(selected)]
    return [selected] + higher[:2]


def load_candidates(path: Path) -> list[str]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    candidates = [str(row.get("candidate_id", "")) for row in payload.get("candidate_boundaries", [])]
    return sorted(x for x in candidates if x)


def audit(project_root: Path) -> dict[str, object]:
    root = project_root.resolve()
    errors: list[str] = []
    boundaries = active(read_tsv(root / "state/lineage_signal_boundary_registry.tsv"), "boundary_id")
    signals = active(read_tsv(root / "state/lineage_signal_registry.tsv"), "signal_id")
    cohorts = active(read_tsv(root / "state/recluster_cohort_registry.tsv"), "cohort_id")
    boundary_by_cohort = {row.get("source_cohort_id", ""): row for row in boundaries if row.get("source_cohort_id")}

    whole = [row for row in boundaries if row.get("boundary_type") == "whole_tissue"]
    if len(whole) != 1:
        errors.append(f"exactly one active whole_tissue signal boundary is required; found {len(whole)}")
    for cohort in cohorts:
        cohort_type = cohort.get("cohort_type", "")
        if cohort_type not in {"broad_class_recluster", "targeted_recluster"}:
            continue
        if cohort.get("status", "").lower() not in TERMINAL_COHORT:
            continue
        if cohort.get("cohort_id", "") not in boundary_by_cohort:
            errors.append(f"terminal cohort lacks a full-lineage signal boundary: {cohort.get('cohort_id','')}")

    signal_key: dict[tuple[str, str, str, str], list[dict[str, str]]] = {}
    for row in signals:
        key = (row.get("boundary_id", ""), row.get("resolution", ""), row.get("cluster", ""), row.get("candidate_lineage", ""))
        signal_key.setdefault(key, []).append(row)

    for boundary in boundaries:
        boundary_id = boundary.get("boundary_id", "")
        if boundary.get("boundary_type") not in {"whole_tissue", "broad_class_recluster", "targeted_recluster"}:
            errors.append(f"{boundary_id}: invalid boundary_type")
        if boundary.get("status") != "audited":
            errors.append(f"{boundary_id}: boundary is not audited")
        required = required_resolutions(boundary)
        audited = split(boundary.get("audited_resolutions", ""))
        missing_res = [value for value in required if value not in audited]
        if missing_res:
            errors.append(f"{boundary_id}: selected plus next higher resolutions were not audited: {missing_res}")

        catalog_path = resolve(root, boundary.get("candidate_catalog", ""))
        universe_path = resolve(root, boundary.get("cluster_universe_artifact", ""))
        if not catalog_path.is_file() or sha256(catalog_path) != boundary.get("candidate_catalog_sha256"):
            errors.append(f"{boundary_id}: candidate catalog is missing or stale")
            continue
        if not universe_path.is_file() or sha256(universe_path) != boundary.get("cluster_universe_sha256"):
            errors.append(f"{boundary_id}: cluster universe is missing or stale")
            continue
        candidates = load_candidates(catalog_path)
        universe = read_tsv(universe_path)
        if not universe or not {"resolution", "cluster"}.issubset(universe[0]):
            errors.append(f"{boundary_id}: cluster universe must contain resolution and cluster")
            continue
        expected = {
            (str(row.get("resolution", "")), str(row.get("cluster", "")), candidate)
            for row in universe if str(row.get("resolution", "")) in audited for candidate in candidates
        }
        observed = {
            (key[1], key[2], key[3]) for key in signal_key if key[0] == boundary_id
            and key[1] in audited and signal_key[key][0].get("candidate_source") == "catalog"
        }
        missing = expected - observed
        full_universe = {
            (str(row.get("resolution", "")), str(row.get("cluster", "")), candidate)
            for row in universe for candidate in candidates
        }
        all_observed = {
            (key[1], key[2], key[3]) for key in signal_key if key[0] == boundary_id
            and signal_key[key][0].get("candidate_source") == "catalog"
        }
        extra = all_observed - full_universe
        if missing:
            errors.append(f"{boundary_id}: {len(missing)} cluster-by-lineage scans are missing")
        if extra:
            errors.append(f"{boundary_id}: {len(extra)} catalog scans fall outside the declared cluster universe")
        duplicates = [key for key, rows in signal_key.items() if key[0] == boundary_id and len(rows) != 1]
        if duplicates:
            errors.append(f"{boundary_id}: duplicate signal rows exist for {len(duplicates)} keys")
        if boundary.get("unexplained_program_audit", "").lower() not in TRUE:
            errors.append(f"{boundary_id}: unexplained multi-gene programs were not audited")
        unexplained_path = resolve(root, boundary.get("unexplained_program_artifact", ""))
        if not unexplained_path.is_file() or sha256(unexplained_path) != boundary.get("unexplained_program_sha256"):
            errors.append(f"{boundary_id}: unexplained-program audit artifact is missing or stale")
        large = number(boundary.get("n_observations")) >= 50000 or number(boundary.get("analysis_fraction")) >= 0.10
        if large:
            if boundary.get("large_label_triggered", "").lower() not in TRUE:
                errors.append(f"{boundary_id}: large-label purity audit was not triggered")
            large_path = resolve(root, boundary.get("large_label_audit_artifact", ""))
            if not large_path.is_file() or sha256(large_path) != boundary.get("large_label_audit_sha256"):
                errors.append(f"{boundary_id}: large-label audit artifact is missing or stale")

    open_signals = []
    for row in signals:
        signal_id = row.get("signal_id", "")
        status = row.get("signal_status", "")
        if status not in SIGNAL_STATES:
            errors.append(f"{signal_id}: invalid signal_status {status!r}")
            continue
        positive = number(row.get("positive_family_count")) > 0 or bool(split(row.get("positive_genes", "")))
        if positive and status == "absent":
            errors.append(f"{signal_id}: positive lineage evidence cannot be recorded as absent")
        if status in {"watch", "candidate", "supported"}:
            if row.get("review_status") != "resolved":
                open_signals.append(signal_id)
            elif row.get("resolution_outcome") not in RESOLVED_OUTCOMES:
                errors.append(f"{signal_id}: resolved signal lacks an allowed biological outcome")
            if not row.get("required_action"):
                errors.append(f"{signal_id}: non-absent signal lacks a required action")
        if status in {"refuted", "watch", "candidate", "supported"} and row.get("review_status") == "resolved":
            artifact = resolve(root, row.get("evidence_artifact", ""))
            if not artifact.is_file() or sha256(artifact) != row.get("evidence_artifact_sha256"):
                errors.append(f"{signal_id}: resolution evidence is missing or stale")
            if not row.get("closure_rationale"):
                errors.append(f"{signal_id}: resolved signal lacks closure rationale")
    if open_signals:
        errors.append(f"{len(open_signals)} lineage signals remain open: {', '.join(open_signals[:10])}")

    return {
        "status": "PASS" if not errors else "BLOCKED",
        "errors": errors,
        "boundaries": len(boundaries),
        "signal_rows": len(signals),
        "open_signals": open_signals,
        "policy": "continuous_open_world_selected_plus_two_higher_available",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = audit(args.project_root)
    out = args.out or args.project_root / "provenance/lineage_signal_coverage_validation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
