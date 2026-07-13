#!/usr/bin/env python3
"""Shared fail-closed audit helpers for multi-route annotation state."""

from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path


ROUTE_CLASSES = {
    "biological_anchor_recluster",
    "interface_deconvolution_review",
    "qc_anchor_recluster",
    "qc_atlas_review",
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


def valid_attempt(root: Path, row: dict[str, str]) -> tuple[bool, list[str]]:
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
                if tier["rctd_medium_low_n"] > 0 and not row.get("fallback_route_attempt_id", "").strip():
                    errors.append("medium/low RCTD observations lack atlas/internal-anchor fallback route ID")
                if truth(row.get("fine_anchor_eligible", "")):
                    errors.append("RCTD-assisted route output cannot itself be fine-anchor eligible")
            except (TypeError, ValueError):
                errors.append("invalid numeric RCTD tier/return fields")
    if route == "qc_atlas_review":
        if not truth(row.get("depth_matched_validation", "")):
            errors.append("QC atlas review lacks depth-matched validation")
        if not truth(row.get("observed_density_spatial_prior", "")):
            errors.append("QC atlas review lacks observed-density spatial prior")
    return not errors, errors


def audit_multiroute(root: Path, context: dict, profile: dict) -> dict:
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
        ok, errors = valid_attempt(root, attempt)
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
            ok, _ = valid_attempt(root, row)
            if ok and (allow_na or row.get("applicability") == "applicable"):
                return True
        return False

    gaps = []
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
            if fallback and fallback.get("route_class") == "qc_atlas_review":
                fallback_ok, _ = valid_attempt(root, fallback)
                fallback_ok = fallback_ok and fallback.get("applicability") == "applicable"
            if not fallback_ok:
                gaps.append({
                    "decision_id": did,
                    "n_observations": n,
                    "required_route": "qc_atlas_review",
                    "reason": f"RCTD medium/low tier ({medium_low}) lacks a validated linked atlas/internal-anchor fallback from {interface_attempt.get('route_attempt_id','')}",
                })
        interface_fraction = n / max(total, 1)
        if state == "interface_review" and policy.get("large_interface_requires_atlas_fallback", False) and (n >= large_min or interface_fraction > float(policy.get("local_interface_max_fraction", 0.02))) and not has(row, "qc_atlas_review", allow_na=True):
            gaps.append({"decision_id": did, "n_observations": n, "required_route": "qc_atlas_review", "reason": "large/nonlocal RCTD reject pool requires calibrated atlas/internal-anchor fallback before terminal interface retention"})
        if (state == "qc_holdout" or from_qc) and not zero_count:
            if n >= large_qc and not has(row, "qc_anchor_recluster"):
                gaps.append({"decision_id": did, "n_observations": n, "required_route": "qc_anchor_recluster", "reason": "large post-clustering QC holdout must be anchor-reclustered before atlas"})
            if not has(row, "qc_atlas_review", allow_na=True):
                gaps.append({"decision_id": did, "n_observations": n, "required_route": "qc_atlas_review", "reason": "nonzero QC holdout lacks depth-matched atlas/internal-anchor/observed-density review"})
        if any(x in text for x in ["oocyte", "germ cell", "germ_cell"]) and not has(row, "strict_rare_cell_review", allow_na=True):
            gaps.append({"decision_id": did, "n_observations": n, "required_route": "strict_rare_cell_review", "reason": "Oocyte/germ identity lacks registered strict rare-cell route"})

    required_views = {"strict", "inclusive", "display"} if policy.get("strict_inclusive_display_required", False) else set()
    view_status = {x.get("view"): x for x in views if x.get("status") == "validated" and artifact_exists(root, x.get("artifact", ""))}
    missing_views = sorted(required_views - set(view_status))
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
    }
