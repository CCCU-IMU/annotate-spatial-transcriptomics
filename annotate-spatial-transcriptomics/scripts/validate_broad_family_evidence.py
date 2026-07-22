#!/usr/bin/env python3
"""Validate the v2 full-feature broad-family evidence Cartesian matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

from evidence_schema_lib import sha256, validate_artifact_ref, validate_json_against_schema


REQUIRED_COLUMNS = {
    "cluster", "n_observations", "candidate_id", "family_name",
    "requested_markers", "available_markers", "missing_markers",
    "available_marker_count", "detected_marker_count", "any_detection_fraction",
    "median_detected_markers_per_observation", "pseudobulk_sum",
    "mean_normalized_expression", "comparative_score_role",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def resolve(payload: dict, dotted: str) -> dict:
    value = payload
    for key in dotted.split("."):
        value = value[key]
    if not isinstance(value, dict):
        raise TypeError(dotted)
    return value


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest", type=Path)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    schema = Path(__file__).resolve().parents[1] / "schemas/broad_family_evidence.schema.json"
    manifest, errors = validate_json_against_schema(args.manifest, schema)
    root = args.manifest.parent
    resolved: dict[str, Path] = {}
    for key in ("source_object", "cluster_membership", "biological_profile", "candidate_catalog", "evidence_table", "validation_manifest"):
        record = manifest.get(key, {}) if isinstance(manifest, dict) else {}
        path, artifact_errors = validate_artifact_ref(root, record, key)
        errors.extend(artifact_errors)
        if path:
            resolved[key] = path

    rows: list[dict[str, str]] = []
    if "validation_manifest" in resolved:
        try:
            validation = json.loads(resolved["validation_manifest"].read_text(encoding="utf-8"))
            if validation.get("status") != "PASS":
                errors.append("full-feature validation manifest did not pass")
            if validation.get("n_observations") is not None and int(validation["n_observations"]) != int(manifest.get("n_observations", -1)):
                errors.append("full-feature validation observation count differs from the evidence manifest")
            if validation.get("normalized_object_sha256") and validation.get("normalized_object_sha256") != manifest.get("source_object", {}).get("sha256"):
                errors.append("full-feature validation manifest is bound to another source object")
        except (OSError, UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError):
            errors.append("full-feature validation manifest is unreadable")
    if "evidence_table" in resolved:
        rows = read_tsv(resolved["evidence_table"])
        fields = set(rows[0]) if rows else set()
        if not rows:
            errors.append("evidence table is empty")
        if missing := sorted(REQUIRED_COLUMNS - fields):
            errors.append("evidence table lacks columns: " + ", ".join(missing))

    expected: set[tuple[str, str, str]] = set()
    if {"biological_profile", "candidate_catalog"} <= set(resolved) and rows:
        profile = json.loads(resolved["biological_profile"].read_text(encoding="utf-8"))
        catalog = json.loads(resolved["candidate_catalog"].read_text(encoding="utf-8"))
        clusters = sorted({row.get("cluster", "") for row in rows})
        if "" in clusters or len(clusters) != int(manifest.get("n_clusters", -1)):
            errors.append("manifest n_clusters differs from the evidence table")
        cluster_counts: dict[str, set[int]] = {}
        for cluster in clusters:
            try:
                cluster_counts[cluster] = {int(row["n_observations"]) for row in rows if row.get("cluster") == cluster}
            except (KeyError, ValueError):
                cluster_counts[cluster] = set()
        if any(len(counts) != 1 or next(iter(counts), 0) <= 0 for counts in cluster_counts.values()):
            errors.append("cluster observation counts are missing, nonpositive or inconsistent across matrix rows")
        elif sum(next(iter(counts)) for counts in cluster_counts.values()) != int(manifest.get("n_observations", -1)):
            errors.append("cluster observation counts do not sum to manifest n_observations")
        candidates = []
        for candidate in catalog.get("candidate_boundaries", []):
            if candidate.get("review_required") is not True:
                continue
            candidate_id = str(candidate.get("candidate_id", ""))
            try:
                program = resolve(profile, str(candidate.get("profile_program", "")))
            except (KeyError, TypeError):
                errors.append(f"cannot resolve profile program for {candidate_id}")
                continue
            families = program.get("positive_families", {})
            if not isinstance(families, dict) or len(families) < 2:
                errors.append(f"{candidate_id} lacks two explicit positive families")
                continue
            candidates.append(candidate_id)
            for cluster in clusters:
                for family in families:
                    expected.add((cluster, candidate_id, family))
        observed = [(row.get("cluster", ""), row.get("candidate_id", ""), row.get("family_name", "")) for row in rows]
        duplicate = sorted(key for key, count in Counter(observed).items() if count != 1)
        if duplicate:
            errors.append(f"duplicate cluster/candidate/family rows: {duplicate[:5]}")
        missing = sorted(expected - set(observed))
        extra = sorted(set(observed) - expected)
        if missing:
            errors.append(f"incomplete Cartesian evidence matrix; missing {len(missing)} rows")
        if extra:
            errors.append(f"evidence matrix contains {len(extra)} unbound rows")
        if set(candidates) != set(manifest.get("candidate_ids", [])):
            errors.append("manifest candidate_ids differ from the bound catalog")
        for line, row in enumerate(rows, 2):
            if row.get("comparative_score_role") != "comparative_not_absence_gate":
                errors.append(f"row {line} permits comparative evidence to act as an absence gate")
            try:
                requested = {x for x in row["requested_markers"].split(",") if x}
                available = {x for x in row["available_markers"].split(",") if x}
                missing_markers = {x for x in row["missing_markers"].split(",") if x}
                if available | missing_markers != requested or available & missing_markers:
                    errors.append(f"row {line} marker availability partition is invalid")
                if int(row["available_marker_count"]) != len(available):
                    errors.append(f"row {line} available marker count is invalid")
                if int(row["detected_marker_count"]) > len(available):
                    errors.append(f"row {line} detected marker count exceeds available markers")
                fraction = float(row["any_detection_fraction"])
                if not 0 <= fraction <= 1:
                    errors.append(f"row {line} detection fraction is outside [0,1]")
                detected = int(row["detected_marker_count"])
                median_detected = float(row["median_detected_markers_per_observation"])
                pseudobulk = float(row["pseudobulk_sum"])
                mean_expression = float(row["mean_normalized_expression"])
                if detected < 0 or not 0 <= median_detected <= len(available):
                    errors.append(f"row {line} detected-marker metrics are invalid")
                if pseudobulk < 0 or mean_expression < 0:
                    errors.append(f"row {line} expression metrics are negative")
            except (KeyError, ValueError):
                errors.append(f"row {line} contains invalid numeric/marker fields")

    result = {
        "status": "PASS" if not errors else "BLOCKED",
        "schema_version": "2.0",
        "manifest": str(args.manifest.resolve()),
        "manifest_sha256": sha256(args.manifest),
        "evidence_rows": len(rows),
        "expected_rows": len(expected),
        "errors": errors,
    }
    out = args.out or args.manifest.with_name("broad_family_evidence_validation.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
