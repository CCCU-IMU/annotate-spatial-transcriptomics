#!/usr/bin/env python3
"""Fail-closed query-derived audit of present and zero-census broad lineages."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

from lineage_decision_lib import observation_writeback_policy


TRUE = {"1", "true", "yes", "pass", "passed", "reviewed"}
PRESENT_REQUIRED = ("query_marker_program_review", "deg_review", "spatial_morphology_review", "observation_level_review", "large_label_embedding_review")
ABSENT_REQUIRED = ("query_marker_program_review", "spatial_morphology_review", "observation_level_review", "selected_plus_two_review", "large_label_embedding_review", "qc_ood_review", "technical_missingness_review")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split(value: str) -> list[str]:
    return [x for x in re.split(r"[;,\s]+", value or "") if x]


def active(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    superseded = {x for row in rows for x in split(row.get("supersedes", ""))}
    return [row for row in rows if row.get("review_id", "") not in superseded]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def framework_tuple(root: Path) -> tuple[int, ...]:
    try:
        value = json.loads((root / "config/project.json").read_text(encoding="utf-8")).get("framework_version", "0")
        return tuple(int(part) for part in str(value).split(".")[:3])
    except (OSError, json.JSONDecodeError, TypeError, ValueError):
        return (0,)


def project_config(root: Path) -> dict:
    try:
        value = json.loads((root / "config/project.json").read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, json.JSONDecodeError, TypeError):
        return {}


def validate_present_evidence(
    path: Path,
    candidate_id: str,
    expected_n: int,
    *,
    policy: dict[str, float],
    require_return_audits: bool = False,
    require_raw_whole_return_audits: bool = False,
    expected_fine_candidates: set[str] | None = None,
) -> list[str]:
    errors: list[str] = []
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return [f"{candidate_id}: present-lineage evidence must be a readable numerical JSON artifact"]
    required = {
        "schema_version", "candidate_lineage", "final_n_observations", "own_program",
        "strongest_competing_program", "observation_level_purity", "spatial_morphology",
    }
    if not isinstance(document, dict) or not required <= set(document):
        return [f"{candidate_id}: present-lineage evidence lacks program, competitor, purity or spatial content"]
    if document.get("schema_version") != "2.0" or document.get("candidate_lineage") != candidate_id:
        errors.append(f"{candidate_id}: present-lineage evidence identity/schema is invalid")
    try:
        if int(document["final_n_observations"]) != expected_n:
            errors.append(f"{candidate_id}: present-lineage evidence count is stale")
        own = document["own_program"]
        competitor = document["strongest_competing_program"]
        purity = document["observation_level_purity"]
        spatial = document["spatial_morphology"]
        family_fractions = own.get("family_detection_fractions", {})
        if int(own.get("positive_family_count", -1)) < 2 or not isinstance(family_fractions, dict) or len(family_fractions) < 2:
            errors.append(f"{candidate_id}: own program lacks two numerical positive families")
        if any(not 0 <= float(value) <= 1 for value in family_fractions.values()):
            errors.append(f"{candidate_id}: own-program family detection fractions are outside [0,1]")
        own_score = float(own["program_score"])
        competitor_score = float(competitor["program_score"])
        if own_score <= competitor_score:
            errors.append(f"{candidate_id}: strongest competing program equals or exceeds the released lineage")
        supported = float(purity["lineage_supported_fraction"])
        competing = float(purity["strongest_competing_fraction"])
        contradiction = float(purity["contradiction_fraction"])
        if any(not 0 <= value <= 1 for value in (supported, competing, contradiction)):
            errors.append(f"{candidate_id}: observation-level purity fractions are outside [0,1]")
        if (
            purity.get("status") != "PASS"
            or supported < policy["present_label_min_lineage_supported_fraction"]
            or supported - competing < policy["present_label_min_purity_margin"]
            or contradiction > policy["maximum_contradiction_fraction"]
        ):
            errors.append(
                f"{candidate_id}: observation-level purity lacks the absolute support floor, margin or contradiction clearance"
            )
        if spatial.get("status") != "PASS" or len(str(spatial.get("rationale", "")).strip()) < 10:
            errors.append(f"{candidate_id}: spatial morphology lacks a substantive PASS review")

        if require_return_audits:
            audits = document.get("source_writeback_audits")
            if not isinstance(audits, list) or not audits:
                errors.append(f"{candidate_id}: present lineage lacks per-source writeback purity audits")
            else:
                seen: set[str] = set()
                accounted = 0
                query_scopes = {"whole_subcluster", "supported_subset", "initial_cluster"}
                for audit in audits:
                    if not isinstance(audit, dict):
                        errors.append(f"{candidate_id}: source writeback audit is not an object")
                        continue
                    audit_id = str(audit.get("return_id", "")).strip()
                    if not audit_id or audit_id in seen:
                        errors.append(f"{candidate_id}: source writeback audit has an empty or duplicate return_id")
                    seen.add(audit_id)
                    try:
                        n = int(audit["n_observations"])
                        accounted += n
                        if n < 1 or len(str(audit["membership_sha256"])) != 64 or audit.get("status") != "PASS":
                            raise ValueError
                        scope = str(audit["return_scope"])
                        if scope in query_scopes:
                            a_supported = float(audit["lineage_supported_fraction"])
                            a_competing = float(audit["strongest_competing_fraction"])
                            a_contradiction = float(audit["contradiction_fraction"])
                            if scope == "supported_subset":
                                minimum = policy["supported_subset_min_lineage_supported_fraction"]
                                margin = policy["supported_subset_min_purity_margin"]
                                purity_ok = (
                                    a_supported >= minimum
                                    and a_supported - a_competing >= margin
                                    and a_contradiction <= policy["maximum_contradiction_fraction"]
                                )
                            else:
                                raw_present = all(
                                    key in audit
                                    for key in (
                                        "raw_two_family_supported_fraction",
                                        "strongest_eligible_competing_raw_fraction",
                                    )
                                )
                                if require_raw_whole_return_audits and not raw_present:
                                    errors.append(f"{candidate_id}: source writeback {audit_id} lacks raw two-family purity")
                                if raw_present:
                                    raw_supported = float(audit["raw_two_family_supported_fraction"])
                                    raw_competing = float(audit["strongest_eligible_competing_raw_fraction"])
                                    purity_ok = (
                                        raw_supported >= policy[
                                            "whole_subcluster_min_raw_two_family_supported_fraction"
                                        ]
                                        and raw_supported - raw_competing >= policy[
                                            "whole_subcluster_min_raw_two_family_margin"
                                        ]
                                        and a_contradiction <= policy["maximum_contradiction_fraction"]
                                    )
                                    if (
                                        raw_competing >= policy[
                                            "whole_subcluster_embedded_competitor_raw_trigger"
                                        ]
                                        and audit.get("embedded_competitor_review_status") != "PASS"
                                    ):
                                        errors.append(
                                            f"{candidate_id}: source writeback {audit_id} lacks closure of a high-coverage embedded competitor"
                                        )
                                else:
                                    purity_ok = (
                                        a_supported >= policy["whole_subcluster_min_lineage_supported_fraction"]
                                        and a_supported - a_competing >= policy["whole_subcluster_min_purity_margin"]
                                        and a_contradiction <= policy["maximum_contradiction_fraction"]
                                    )
                            if not purity_ok:
                                errors.append(f"{candidate_id}: source writeback {audit_id} fails its scope-specific purity gate")
                        elif scope in {"atlas_qc_return", "canonical_oocyte_cluster"}:
                            if audit.get("route_validation_status") != "PASS":
                                errors.append(f"{candidate_id}: source writeback {audit_id} lacks its route-specific validation")
                        else:
                            errors.append(f"{candidate_id}: source writeback {audit_id} has an unknown return_scope")
                    except (KeyError, TypeError, ValueError):
                        errors.append(f"{candidate_id}: source writeback audit {audit_id or '<unknown>'} is malformed")
                if accounted != expected_n:
                    errors.append(f"{candidate_id}: source writeback audits do not exactly account for final membership")

        expected_fine_candidates = expected_fine_candidates or set()
        if expected_fine_candidates:
            audits = document.get("fine_candidate_audit")
            if not isinstance(audits, list):
                errors.append(f"{candidate_id}: optional fine candidates were not systematically audited")
            else:
                by_id = {str(item.get("candidate_id", "")): item for item in audits if isinstance(item, dict)}
                if set(by_id) != expected_fine_candidates:
                    errors.append(f"{candidate_id}: fine-candidate audit does not cover the complete parent-specific catalog")
                for fine_id, audit in by_id.items():
                    status = audit.get("status")
                    if status not in {"supported", "refuted", "not_evaluable"}:
                        errors.append(f"{candidate_id}: fine candidate {fine_id} has an invalid audit status")
                    if status in {"supported", "refuted"}:
                        channels = audit.get("evidence_channels", [])
                        if not isinstance(channels, list) or len(set(channels)) < 2:
                            errors.append(f"{candidate_id}: fine candidate {fine_id} lacks two independent evidence channels")
                    if status == "not_evaluable" and len(str(audit.get("rationale", "")).strip()) < 10:
                        errors.append(f"{candidate_id}: fine candidate {fine_id} lacks a substantive not-evaluable rationale")
    except (KeyError, TypeError, ValueError):
        errors.append(f"{candidate_id}: present-lineage numerical evidence is malformed")
    return errors


def validate(root: Path, catalog_path: Path) -> dict[str, object]:
    root = root.resolve()
    errors: list[str] = []
    config = project_config(root)
    try:
        writeback_policy = observation_writeback_policy(config)
    except (TypeError, ValueError):
        writeback_policy = observation_writeback_policy()
        errors.append("project observation_writeback_policy is invalid")
    require_numerical_present = (
        config.get("numerical_broad_completeness_evidence_required") is True
        or framework_tuple(root) >= (2, 0, 2)
    )
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    fine_catalog = catalog.get("machine_actionable_fine_candidate_catalog", {})
    candidates = {str(x.get("candidate_id", "")): x for x in catalog.get("candidate_boundaries", []) if x.get("review_required")}
    required_zero = {cid for cid, row in candidates.items() if row.get("release_level") == "default_broad_candidate"}
    rows = active(read_tsv(root / "state/broad_class_completeness_registry.tsv"))
    by_candidate: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_candidate.setdefault(row.get("candidate_lineage", ""), []).append(row)
    duplicates = [cid for cid, items in by_candidate.items() if len(items) != 1]
    if duplicates:
        errors.append("duplicate active completeness reviews: " + ", ".join(sorted(duplicates)))
    reviewed: list[str] = []
    zero_census: list[str] = []
    for cid, items in by_candidate.items():
        if len(items) != 1:
            continue
        row = items[0]
        try:
            n = int(float(row.get("final_n_observations", "0") or 0))
        except ValueError:
            errors.append(f"{cid}: invalid final_n_observations")
            continue
        if cid not in candidates:
            errors.append(f"{cid}: candidate is not present in the bound catalog")
        if row.get("status") != "audited":
            errors.append(f"{cid}: completeness review is not audited")
        if not row.get("closure_rationale"):
            errors.append(f"{cid}: closure rationale is missing")
        artifact = Path(row.get("evidence_artifact", ""))
        artifact = artifact if artifact.is_absolute() else root / artifact
        if not artifact.is_file() or sha256(artifact) != row.get("evidence_artifact_sha256"):
            errors.append(f"{cid}: evidence artifact is missing or stale")
        if n > 0:
            reviewed.append(cid)
            if row.get("census_status") != "present" or row.get("decision") not in {"supported", "revised", "downgraded"}:
                errors.append(f"{cid}: present lineage lacks a supported/revised/downgraded decision")
            for field in PRESENT_REQUIRED:
                if row.get(field, "").lower() not in TRUE:
                    errors.append(f"{cid}: present-lineage audit lacks {field}")
            if require_numerical_present and artifact.is_file():
                expected_fine = {
                    str(item.get("candidate_id", ""))
                    for item in fine_catalog.get(cid, [])
                    if isinstance(item, dict) and item.get("candidate_id")
                }
                errors.extend(
                    validate_present_evidence(
                        artifact,
                        cid,
                        n,
                        policy=writeback_policy,
                        require_return_audits=config.get("terminal_return_purity_audit_required") is True,
                        require_raw_whole_return_audits=config.get("raw_two_family_writeback_audit_required") is True,
                        expected_fine_candidates=(
                            expected_fine if config.get("complete_fine_candidate_audit_required") is True else set()
                        ),
                    )
                )
        else:
            zero_census.append(cid)
            if cid in required_zero:
                if row.get("census_status") != "absent" or row.get("decision") != "refuted_by_multichannel_query_evidence":
                    errors.append(f"{cid}: zero-census default lineage lacks multichannel query-derived refutation")
                for field in ABSENT_REQUIRED:
                    if row.get(field, "").lower() not in TRUE:
                        errors.append(f"{cid}: zero-census audit lacks {field}")
                if row.get("atlas_role") not in {"challenger_only", "not_used", "discordant_low_confidence"}:
                    errors.append(f"{cid}: Atlas cannot establish absence")
    missing_zero = sorted(cid for cid in required_zero if cid not in by_candidate)
    if missing_zero:
        errors.append("default broad candidates lack a completeness review: " + ", ".join(missing_zero))
    return {"status": "PASS" if not errors else "BLOCKED", "catalog_id": catalog.get("catalog_id"),
            "numerical_present_lineage_evidence_required": require_numerical_present,
            "active_reviews": len(rows), "present_reviewed": sorted(reviewed), "zero_census_reviewed": sorted(zero_census),
            "required_default_candidates": sorted(required_zero), "errors": errors}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = validate(args.project_root, args.catalog.resolve())
    out = args.out or args.project_root / "provenance/broad_class_completeness_validation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
