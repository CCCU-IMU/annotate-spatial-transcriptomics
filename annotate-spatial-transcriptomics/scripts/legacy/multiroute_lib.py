#!/usr/bin/env python3
"""Legacy v1.5 pool/route validators; migration compatibility only."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path


ROUTE_CLASSES = {
    "biological_anchor_recluster",
    "interface_deconvolution_review",
    "qc_anchor_recluster",
    "qc_atlas_review",
    "context_specific_identity_review",
    "strict_rare_cell_review",
    "diagnostic_supervised_review",
    "explicit_technical_retention",
}
TERMINAL_ROUTE_STATUS = {"validated", "not_applicable_reviewed"}
KNOWN_ROUTE_STATUS = TERMINAL_ROUTE_STATUS | {"prepared", "submitted", "running", "failed_preserved", "cancelled_preserved"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split_ids(value: str) -> set[str]:
    return {x for x in re.split(r"[;,\s]+", (value or "").strip()) if x}


def decision_id(row: dict[str, str]) -> str:
    return row.get("decision_id") or f"{row.get('decision_version','')}:{row.get('source_run_id','')}:{row.get('source_cluster','')}"


def active_decisions(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    superseded: set[str] = set()
    for row in rows:
        superseded.update(split_ids(row.get("supersedes", "")))
    return [row for row in rows if decision_id(row) not in superseded]


def truth(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "pass", "passed"}


def artifact_exists(root: Path, value: str) -> bool:
    if not value:
        return False
    path = Path(value)
    return (path if path.is_absolute() else root / path).exists()


def artifact_path(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def membership_ids(root: Path, value: str) -> set[str]:
    """Read a frozen one-row-per-observation membership artifact."""
    if not value:
        return set()
    rows = read_tsv(artifact_path(root, value))
    if not rows or "cell_id" not in rows[0]:
        return set()
    ids = [str(row.get("cell_id", "")).strip() for row in rows]
    if any(not item for item in ids) or len(ids) != len(set(ids)):
        return set()
    return set(ids)


def validate_membership_artifact(
    root: Path, row: dict[str, str], artifact_field: str, sha_field: str, n_field: str = "n_query",
    allow_empty: bool = False,
) -> tuple[set[str], list[str]]:
    errors: list[str] = []
    value = row.get(artifact_field, "").strip()
    expected_sha = row.get(sha_field, "").strip()
    path = artifact_path(root, value) if value else Path()
    if not value or not path.is_file():
        return set(), [f"missing {artifact_field}"]
    if not expected_sha or sha256(path) != expected_sha:
        errors.append(f"{sha_field} does not match {artifact_field}")
    try:
        expected_n = int(float(row.get(n_field, 0) or 0))
        ids = membership_ids(root, value)
        if not ids:
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt", encoding="utf-8") as handle:
                header = handle.readline().rstrip("\r\n").split("\t")
            if not (allow_empty and expected_n == 0 and "cell_id" in header):
                errors.append(f"{artifact_field} is not a unique nonempty cell_id membership table")
        if expected_n != len(ids):
            errors.append(f"{artifact_field} row count differs from {n_field}")
    except ValueError:
        ids = set()
        errors.append(f"invalid {n_field}")
    return ids, errors


def parse_resolutions(value: str) -> list[float]:
    try:
        return [float(item) for item in re.split(r"[,;\s]+", (value or "").strip()) if item]
    except ValueError:
        return []


def validate_biological_profile(profile: dict) -> list[str]:
    errors = []
    if profile.get("profile_role") != "biological_evidence":
        errors.append(
            "active biological profile must have profile_role=biological_evidence; "
            f"observed {profile.get('profile_role')!r}"
        )
    for key in ["lineages", "evidence_gates", "release_taxonomy", "resolution_policy", "multi_route_policy"]:
        if not isinstance(profile.get(key), dict) or not profile.get(key):
            errors.append(f"biological profile lacks nonempty {key}")
    return errors


def valid_attempt(root: Path, row: dict[str, str], profile: dict | None = None) -> tuple[bool, list[str]]:
    errors: list[str] = []
    route = row.get("route_class", "")
    if route not in ROUTE_CLASSES:
        errors.append(f"unknown route_class={route}")
    status = row.get("status", "")
    if status not in TERMINAL_ROUTE_STATUS:
        errors.append(f"nonterminal status={status}")
    if not artifact_exists(root, row.get("validation_artifact", "")):
        errors.append("missing validation_artifact")
    applicability = row.get("applicability", "")
    if applicability == "not_applicable":
        if status != "not_applicable_reviewed":
            errors.append("not_applicable route lacks not_applicable_reviewed status")
        if len(row.get("applicability_rationale", "").strip()) < 20:
            errors.append("not_applicable rationale is too short")
        return not errors, errors
    if applicability != "applicable":
        errors.append(f"invalid applicability={applicability}")
    if route in {"biological_anchor_recluster", "qc_anchor_recluster"}:
        try:
            if int(float(row.get("n_anchors", 0) or 0)) <= 0:
                errors.append("anchor recluster has no anchors")
        except ValueError:
            errors.append("invalid n_anchors")
        if not truth(row.get("query_only_graph", "")):
            errors.append("anchor recluster lacks query-only graph/DEG boundary")
        _, membership_errors = validate_membership_artifact(
            root, row, "query_membership_artifact", "query_membership_sha256"
        )
        errors.extend("anchor recluster " + error for error in membership_errors)
        if route == "qc_anchor_recluster":
            if row.get("source_state", "") != "qc_holdout":
                errors.append("QC anchor recluster source_state must be qc_holdout")
            if not row.get("pool_snapshot_id", "").strip():
                errors.append("QC anchor recluster lacks frozen QC-holdout pool snapshot")
            if not artifact_exists(root, row.get("outcome_artifact", "")):
                errors.append("QC anchor recluster lacks validated outcome artifact")
            residual_ids, residual_errors = validate_membership_artifact(
                root, row, "residual_qc_membership_artifact", "residual_qc_membership_sha256", "n_qc_retained", allow_empty=True
            )
            errors.extend("QC anchor recluster " + error for error in residual_errors)
            query_ids = membership_ids(root, row.get("query_membership_artifact", ""))
            if residual_ids and not residual_ids.issubset(query_ids):
                errors.append("QC anchor residual membership is not a subset of its full QC input")
        if profile:
            expected = [float(value) for value in profile.get("resolution_policy", {}).get("formal_candidate_resolutions", [])]
            observed = parse_resolutions(row.get("candidate_resolutions", ""))
            if expected and observed != expected:
                errors.append(f"formal candidate_resolutions must equal {expected}; observed {observed}")
            try:
                selected = float(row.get("selected_resolution", ""))
                if expected and selected not in expected:
                    errors.append(f"selected_resolution={selected} is outside the formal candidate grid")
            except ValueError:
                errors.append("anchor recluster lacks numeric selected_resolution")
    if route == "interface_deconvolution_review":
        if not truth(row.get("depth_matched_validation", "")):
            errors.append("applicable interface review lacks depth-matched validation")
        tier_fields = [
            "rctd_extreme_n", "rctd_high_n", "rctd_medium_low_n",
            "rctd_fine_return_n", "rctd_broad_return_n",
        ]
        missing_tiers = [field for field in tier_fields if row.get(field, "") == ""]
        if missing_tiers:
            errors.append("applicable interface review lacks tier fields: " + ",".join(missing_tiers))
        else:
            try:
                tier = {field: int(float(row[field])) for field in tier_fields}
                n_query = int(float(row.get("n_query", 0) or 0))
                n_rerouted = int(float(row.get("n_rerouted", 0) or 0))
                n_defined_fine = int(float(row.get("n_defined_fine", 0) or 0))
                n_defined_broad = int(float(row.get("n_defined_broad_only", 0) or 0))
                if any(value < 0 for value in tier.values()):
                    errors.append("RCTD tier/return counts must be nonnegative")
                if tier["rctd_extreme_n"] + tier["rctd_high_n"] + tier["rctd_medium_low_n"] != n_query:
                    errors.append("RCTD extreme/high/medium-low counts do not partition n_query")
                if tier["rctd_fine_return_n"] > tier["rctd_extreme_n"]:
                    errors.append("RCTD fine returns exceed extreme-confidence evidence")
                if tier["rctd_fine_return_n"] + tier["rctd_broad_return_n"] > tier["rctd_extreme_n"] + tier["rctd_high_n"]:
                    errors.append("RCTD fine+broad returns exceed extreme+high-confidence evidence")
                if tier["rctd_fine_return_n"] != n_defined_fine or tier["rctd_broad_return_n"] != n_defined_broad:
                    errors.append("RCTD return counts disagree with route outcome counts")
                if tier["rctd_fine_return_n"] > 0 and not truth(row.get("independent_fine_evidence", "")):
                    errors.append("RCTD fine return lacks independent marker/recluster/spatial evidence")
                if tier["rctd_medium_low_n"] > n_rerouted:
                    errors.append("medium/low RCTD observations were not all rerouted")
                if tier["rctd_fine_return_n"] + tier["rctd_broad_return_n"] + n_rerouted != n_query:
                    errors.append("RCTD fine, broad and rerouted outcomes do not partition n_query")
                if n_rerouted > 0:
                    if row.get("successor_state", "") != "qc_holdout":
                        errors.append("RCTD rerouted observations must have successor_state=qc_holdout")
                    _, reroute_errors = validate_membership_artifact(
                        root, row, "reroute_membership_artifact", "reroute_membership_sha256", "n_rerouted"
                    )
                    errors.extend("RCTD reroute " + error for error in reroute_errors)
                if tier["rctd_medium_low_n"] > 0 and not row.get("fallback_route_attempt_id", "").strip():
                    errors.append("medium/low RCTD observations lack QC-holdout anchor successor route ID")
                if truth(row.get("fine_anchor_eligible", "")):
                    errors.append("RCTD-assisted route output cannot itself be fine-anchor eligible")
            except (TypeError, ValueError):
                errors.append("invalid numeric RCTD tier/return fields")
    if route == "qc_atlas_review":
        if row.get("source_state", "") != "qc_holdout_residual_after_anchor":
            errors.append("QC atlas source_state must be qc_holdout_residual_after_anchor")
        if not truth(row.get("depth_matched_validation", "")):
            errors.append("QC atlas review lacks depth-matched validation")
        if not truth(row.get("observed_density_spatial_prior", "")):
            errors.append("QC atlas review lacks observed-density spatial prior")
        if row.get("calibration_origin", "") != "query_like_heldout_current_query_anchors":
            errors.append("QC atlas review was not calibrated on query-like held-out current-query anchors")
        _, membership_errors = validate_membership_artifact(
            root, row, "query_membership_artifact", "query_membership_sha256"
        )
        errors.extend("QC atlas review " + error for error in membership_errors)
        if not artifact_exists(root, row.get("calibration_manifest", "")):
            errors.append("QC atlas review lacks calibration origin manifest")
        if not row.get("pool_snapshot_id", "").strip():
            errors.append("QC atlas review lacks its residual QC-holdout pool snapshot")
        if not row.get("parent_pool_snapshot_id", "").strip():
            errors.append("QC atlas review lacks parent complete-QC pool snapshot")
        if not row.get("prerequisite_route_attempt_id", "").strip():
            errors.append("QC atlas review lacks prerequisite QC-anchor route ID")
        if not row.get("prerequisite_outcome_sha256", "").strip():
            errors.append("QC atlas review lacks prerequisite QC-anchor outcome SHA256")
        if not row.get("prerequisite_residual_qc_sha256", "").strip():
            errors.append("QC atlas review lacks prerequisite residual-QC membership SHA256")
    return not errors, errors


def validate_qc_atlas_prerequisite(
    root: Path,
    atlas: dict[str, str],
    attempts: list[dict[str, str]],
    profile: dict,
) -> tuple[bool, list[str]]:
    """Verify that an Atlas query is exactly the residual child of a prior QC-anchor route."""
    errors: list[str] = []
    by_id = {row.get("route_attempt_id", ""): row for row in attempts}
    prerequisite_id = atlas.get("prerequisite_route_attempt_id", "").strip()
    prerequisite = by_id.get(prerequisite_id)
    if not prerequisite or prerequisite.get("route_class") != "qc_anchor_recluster":
        return False, ["QC Atlas is not linked to a preceding QC-holdout anchor recluster"]
    prerequisite_ok, prerequisite_errors = valid_attempt(root, prerequisite, profile)
    if not prerequisite_ok or prerequisite.get("applicability") != "applicable":
        errors.append("QC Atlas prerequisite QC-anchor route is not a valid applicable terminal route")
        errors.extend("prerequisite: " + error for error in prerequisite_errors)
    outcome_value = prerequisite.get("outcome_artifact", "")
    outcome = artifact_path(root, outcome_value) if outcome_value else Path()
    if not outcome_value or not outcome.is_file() or sha256(outcome) != atlas.get("prerequisite_outcome_sha256", ""):
        errors.append("QC Atlas prerequisite outcome SHA256 does not match")
    if prerequisite.get("pool_snapshot_id", "") != atlas.get("parent_pool_snapshot_id", ""):
        errors.append("QC Atlas parent snapshot differs from the complete-QC input snapshot")
    if prerequisite.get("pool_snapshot_id", "") == atlas.get("pool_snapshot_id", ""):
        errors.append("QC Atlas must use a distinct residual-QC child snapshot, not the full QC snapshot")
    residual_sha = prerequisite.get("residual_qc_membership_sha256", "")
    if not residual_sha or residual_sha != atlas.get("prerequisite_residual_qc_sha256", ""):
        errors.append("QC Atlas prerequisite residual-QC membership SHA256 does not match")
    if atlas.get("query_membership_sha256", "") != residual_sha:
        errors.append("QC Atlas query SHA256 is not the preceding QC-anchor residual-QC SHA256")
    residual_ids = membership_ids(root, prerequisite.get("residual_qc_membership_artifact", ""))
    atlas_ids = membership_ids(root, atlas.get("query_membership_artifact", ""))
    if not residual_ids or atlas_ids != residual_ids:
        errors.append("QC Atlas query IDs are not exactly the preceding QC-anchor residual-QC IDs")
    if not prerequisite.get("created_at", "") or not atlas.get("created_at", "") or prerequisite.get("created_at", "") >= atlas.get("created_at", ""):
        errors.append("QC Atlas prerequisite does not precede the Atlas attempt")
    return not errors, errors


def audit_multiroute(root: Path, context: dict, profile: dict) -> dict:
    profile_errors = validate_biological_profile(profile)
    if profile_errors:
        return {
            "status": "BLOCKED", "configuration_errors": profile_errors,
            "large_pool_threshold": None, "large_qc_threshold": None,
            "active_decisions": 0, "route_attempts": 0, "historical_route_attempts": 0,
            "invalid_attempts": [], "gaps": [], "required_views": ["final"],
            "missing_views": ["final"], "phase_status": {},
        }
    decision_history = read_tsv(root / "state/cluster_decision_ledger.tsv")
    decisions = active_decisions(decision_history)
    decision_by_id = {decision_id(row): row for row in decision_history}
    all_attempts = read_tsv(root / "state/route_attempt_registry.tsv")
    superseded_attempts=set()
    for attempt in all_attempts:superseded_attempts.update(split_ids(attempt.get("supersedes","")))
    attempts=[x for x in all_attempts if x.get("route_attempt_id") not in superseded_attempts]
    views = read_tsv(root / "state/annotation_view_registry.tsv")
    policy = profile.get("multi_route_policy", {})
    total = sum(int(float(x.get("n_observations", 0) or 0)) for x in decisions)
    large_min = max(
        int(policy.get("large_pool_min_observations", 1000)),
        int(total * float(policy.get("large_pool_min_fraction", 0.005))),
    )
    large_qc = int(policy.get("large_qc_min_observations", 1000))
    priority = [x.lower() for x in policy.get("priority_lineages_require_anchor_recluster", context.get("priority_lineages", []))]
    by_decision: dict[str, list[dict[str, str]]] = {}
    by_attempt_id: dict[str, dict[str, str]] = {}
    invalid_attempts = []
    for attempt in attempts:
        by_attempt_id[attempt.get("route_attempt_id", "")] = attempt
        for did in split_ids(attempt.get("decision_id", "")):
            by_decision.setdefault(did, []).append(attempt)
    for attempt in attempts:
        ok, errors = valid_attempt(root, attempt, profile)
        if not ok and attempt.get("status") in TERMINAL_ROUTE_STATUS:
            invalid_attempts.append({"route_attempt_id": attempt.get("route_attempt_id", ""), "errors": errors})
        elif attempt.get("status") not in KNOWN_ROUTE_STATUS:
            invalid_attempts.append({"route_attempt_id": attempt.get("route_attempt_id", ""), "errors": [f"unknown status={attempt.get('status','')}"]})

    def attempt_ancestry(decision: dict[str, str]) -> list[dict[str, str]]:
        """Resolve routes attached directly to a decision and transitively through route IDs."""
        did = decision_id(decision)
        queue = list(split_ids(decision.get("route", ""))) + [did]
        seen_tokens: set[str] = set(); found: dict[str, dict[str, str]] = {}
        while queue:
            token = queue.pop()
            if token in seen_tokens: continue
            seen_tokens.add(token)
            # A child route may reference the decision IDs that created its
            # frozen membership rather than the parent route ID directly.
            # Those decisions can later be superseded, but their route field
            # remains the authoritative ancestry edge and must still be
            # traversed to avoid repeating an already completed parent route.
            historical_decision = decision_by_id.get(token)
            if historical_decision:
                queue.extend(split_ids(historical_decision.get("route", "")))
            direct = []
            if token in by_attempt_id: direct.append(by_attempt_id[token])
            direct.extend(by_decision.get(token, []))
            for attempt in direct:
                aid = attempt.get("route_attempt_id", "")
                if aid and aid not in found: found[aid] = attempt
                queue.extend(split_ids(attempt.get("decision_id", "")))
        return list(found.values())

    def has(decision: dict[str, str], route: str, allow_na: bool = False) -> bool:
        for row in attempt_ancestry(decision):
            if row.get("route_class") != route:
                continue
            ok, _ = valid_attempt(root, row, profile)
            if ok and (allow_na or row.get("applicability") == "applicable"):
                return True
        return False

    gaps = []
    for atlas in attempts:
        if atlas.get("route_class") != "qc_atlas_review" or atlas.get("applicability") != "applicable":
            continue
        atlas_ok, _ = valid_attempt(root, atlas, profile)
        if not atlas_ok:
            continue
        prerequisite_ok, prerequisite_errors = validate_qc_atlas_prerequisite(root, atlas, attempts, profile)
        if not prerequisite_ok:
            gaps.append({
                "decision_id": atlas.get("decision_id", "") or "__PROJECT_PHASE__",
                "n_observations": int(float(atlas.get("n_query", 0) or 0)),
                "required_route": "qc_anchor_recluster",
                "reason": "; ".join(prerequisite_errors),
            })
    for row in decisions:
        did = decision_id(row)
        n = int(float(row.get("n_observations", 0) or 0))
        state = row.get("state", "")
        text = " ".join([row.get("broad_label", ""), row.get("fine_label", ""), row.get("target_pool", "")]).lower()
        evidence = " ".join([row.get("route", ""), row.get("evidence_status", ""), row.get("validation_status", "")]).lower()
        zero_count = "zero_count" in evidence
        from_qc = any(x in evidence for x in ["qc_", "atlas", "depth_matched", "mapping_fullfeature"])
        priority_hit = any(x in text for x in priority)
        direct_first_pass = any(x in evidence for x in ["direct_definition", "direct_broad", "freeze preliminary anchor"])
        requires_bio = not from_qc and (
            (state in {"defined_broad_only", "interface_review", "pending_review"} and n >= large_min)
            or (state in {"defined_fine", "defined_broad_only", "interface_review", "pending_review"} and priority_hit and n >= large_min)
            or (policy.get("large_direct_definition_requires_anchor_purity_review", False) and state in {"defined_fine", "defined_broad_only"} and direct_first_pass and n >= large_min)
        )
        if requires_bio and not has(row, "biological_anchor_recluster"):
            gaps.append({"decision_id": did, "n_observations": n, "required_route": "biological_anchor_recluster", "reason": "large/priority biological pool lacks validated balanced-anchor query-only reclustering"})
        if state == "interface_review" and not has(row, "interface_deconvolution_review", allow_na=True):
            gaps.append({"decision_id": did, "n_observations": n, "required_route": "interface_deconvolution_review", "reason": "interface lacks calibrated deconvolution applicability review"})
        if state == "interface_review" and not has(row, "biological_anchor_recluster"):
            gaps.append({"decision_id": did, "n_observations": n, "required_route": "biological_anchor_recluster", "reason": "interface lacks targeted current-query anchor reclustering"})
        for interface_attempt in attempt_ancestry(row):
            if interface_attempt.get("route_class") != "interface_deconvolution_review" or interface_attempt.get("applicability") != "applicable":
                continue
            try:
                medium_low = int(float(interface_attempt.get("rctd_medium_low_n", 0) or 0))
            except ValueError:
                medium_low = 0
            if medium_low <= 0:
                continue
            fallback_id = interface_attempt.get("fallback_route_attempt_id", "").strip()
            fallback = by_attempt_id.get(fallback_id)
            fallback_ok = False
            if fallback and fallback.get("route_class") == "qc_anchor_recluster":
                fallback_ok, _ = valid_attempt(root, fallback, profile)
                fallback_ok = fallback_ok and fallback.get("applicability") == "applicable"
                expected_hash = interface_attempt.get("fallback_expected_membership_sha256", "").strip()
                observed_hash = fallback.get("query_membership_sha256", "").strip()
                fallback_ok = fallback_ok and bool(expected_hash) and expected_hash == observed_hash
                reroute_ids = membership_ids(root, interface_attempt.get("reroute_membership_artifact", ""))
                qc_input_ids = membership_ids(root, fallback.get("query_membership_artifact", ""))
                fallback_ok = fallback_ok and bool(reroute_ids) and reroute_ids.issubset(qc_input_ids)
            if not fallback_ok:
                gaps.append({
                    "decision_id": did,
                    "n_observations": n,
                    "required_route": "qc_anchor_recluster",
                    "reason": f"RCTD medium/low tier ({medium_low}) was not included in a validated complete QC-holdout anchor successor from {interface_attempt.get('route_attempt_id','')}",
                })
        interface_fraction = n / max(total, 1)
        large_interface_policy = policy.get("large_interface_requires_qc_holdout", policy.get("large_interface_requires_atlas_fallback", False))
        if state == "interface_review" and large_interface_policy and (n >= large_min or interface_fraction > float(policy.get("local_interface_max_fraction", 0.02))):
            gaps.append({"decision_id": did, "n_observations": n, "required_route": "reroute_to_biological_pool_or_qc_holdout", "reason": "large/nonlocal reject population cannot close as an anatomical interface and cannot call Atlas directly"})
        if (state == "qc_holdout" or from_qc) and not zero_count:
            if not has(row, "qc_anchor_recluster"):
                gaps.append({"decision_id": did, "n_observations": n, "required_route": "qc_anchor_recluster", "reason": "large post-clustering QC holdout must be anchor-reclustered before atlas"})
            atlas_ancestry_required = state == "qc_holdout" or any(x in evidence for x in ["atlas", "mapping_fullfeature"])
            if atlas_ancestry_required and not has(row, "qc_atlas_review", allow_na=True):
                gaps.append({"decision_id": did, "n_observations": n, "required_route": "qc_atlas_review", "reason": "nonzero QC holdout lacks depth-matched atlas/internal-anchor/observed-density review"})
            ancestry = attempt_ancestry(row)
            anchors = [item for item in ancestry if item.get("route_class") == "qc_anchor_recluster" and item.get("applicability") == "applicable"]
            atlases = [item for item in ancestry if item.get("route_class") == "qc_atlas_review" and item.get("applicability") == "applicable"]
            if anchors and atlases:
                first_anchor = min(item.get("created_at", "") for item in anchors)
                first_atlas = min(item.get("created_at", "") for item in atlases)
                if not first_anchor or not first_atlas or first_atlas <= first_anchor:
                    gaps.append({"decision_id": did, "n_observations": n, "required_route": "qc_anchor_recluster", "reason": "QC atlas route predates the required full QC-pool anchor recluster"})
        if any(x in text for x in ["oocyte", "germ cell", "germ_cell"]) and not (
            has(row, "context_specific_identity_review", allow_na=True)
            or has(row, "strict_rare_cell_review", allow_na=True)
        ):
            gaps.append({"decision_id": did, "n_observations": n, "required_route": "context_specific_identity_review", "reason": "Oocyte/germ identity lacks its contamination-safe context-specific validation route"})

    required_views = {"final"} if policy.get("final_annotation_required", False) else set()
    view_status = {x.get("view"): x for x in views if x.get("status") == "validated" and artifact_exists(root, x.get("artifact", ""))}
    missing_views = sorted(required_views - set(view_status))
    active_interface_n = sum(int(float(row.get("n_observations", 0) or 0)) for row in decisions if row.get("state") == "interface_review")
    active_qc_n = sum(int(float(row.get("n_observations", 0) or 0)) for row in decisions if row.get("state") == "qc_holdout")
    oocyte_required = any(term in " ".join(context.get("priority_lineages", [])).lower() for term in ["oocyte", "germ"]) or any(
        term in (row.get("broad_label", "") + " " + row.get("fine_label", "")).lower()
        for row in decisions for term in ["oocyte", "germ cell", "germ_cell"]
    )
    phase_status = {
        "whole_tissue_broad": "PASS" if decisions else "BLOCKED",
        "route_a_biological_pool": "PASS" if any(item.get("route_class") == "biological_anchor_recluster" and valid_attempt(root, item, profile)[0] for item in attempts) else "BLOCKED",
        "route_b_interface": "PASS" if active_interface_n == 0 or all(has(row, "interface_deconvolution_review", allow_na=True) and has(row, "biological_anchor_recluster") for row in decisions if row.get("state") == "interface_review") else "BLOCKED",
        "route_c_qc_anchor_then_residual_atlas": "PASS" if active_qc_n == 0 or all(has(row, "qc_anchor_recluster") and has(row, "qc_atlas_review", allow_na=True) for row in decisions if row.get("state") == "qc_holdout") else "BLOCKED",
        "route_d_oocyte_specific_if_candidate": "PASS" if not oocyte_required or any(item.get("route_class") in {"context_specific_identity_review", "strict_rare_cell_review"} and valid_attempt(root, item, profile)[0] for item in attempts) else "BLOCKED",
        "final_single_annotation": "PASS" if "final" in view_status else "BLOCKED",
    }
    for phase, status in phase_status.items():
        if status != "PASS":
            gaps.append({"decision_id": "__PROJECT_PHASE__", "n_observations": 0, "required_route": phase, "reason": f"standard sheep-ovary phase {phase} is not complete"})
    return {
        "status": "PASS" if not gaps and not invalid_attempts and not missing_views else "BLOCKED",
        "large_pool_threshold": large_min,
        "large_qc_threshold": large_qc,
        "active_decisions": len(decisions),
        "route_attempts": len(attempts),
        "historical_route_attempts": len(all_attempts),
        "invalid_attempts": invalid_attempts,
        "gaps": gaps,
        "required_views": sorted(required_views),
        "missing_views": missing_views,
        "phase_status": phase_status,
        "configuration_errors": [],
    }
