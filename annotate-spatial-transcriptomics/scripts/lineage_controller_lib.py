#!/usr/bin/env python3
"""Deterministic v2.2 lineage-controller primitives.

This module owns the release-critical interpretation of observation scores.
R stages compute expression and spatial evidence; this module validates
candidate-local subsets, resolves overlaps without catalog-order assignment,
and closes exact remainders from immutable scores.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from controller_thresholds import (
    local_split_trigger_defaults, observation_writeback_defaults,
)


RELEASE_ROLES = {"broad", "fine"}
GENERIC_REMAINDER_IDS = {"stromal_mesenchymal"}
QC_STATES = {"qc_holdout", "unknown_candidate", "technical_state"}
CONTEXT_SUPPORTED_STATUSES = {"supported", "pass"}
CONTEXT_ALLOWED_STATUSES = CONTEXT_SUPPORTED_STATUSES | {
    "not_evaluable", "refuted",
}
_WRITEBACK_DEFAULTS = observation_writeback_defaults()
_LOCAL_SPLIT_DEFAULTS = local_split_trigger_defaults()


def clamp(value: float) -> float:
    return max(0.0, min(1.0, value))


def truth(value: object) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "pass", "passed"}


def number(value: object, default: float = 0.0) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode.replace("t", ""), encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def iter_tsv(path: Path):
    """Stream a TSV/TSV.GZ without materializing every row."""
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: Iterable[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is not None:
        columns = list(fields)
    else:
        columns = []
        seen: set[str] = set()
        for row in rows:
            for column in row:
                if column not in seen:
                    seen.add(column)
                    columns.append(column)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def catalog_candidates(catalog: dict) -> dict[str, dict]:
    rows = list(catalog.get("candidate_boundaries", []))
    broad_by_label = {
        str(row.get("release_broad_label", "")).strip(): row
        for row in rows
        if str(row.get("candidate_role", "")).lower() == "broad"
        and str(row.get("release_broad_label", "")).strip()
    }
    represented_fine_labels = {
        str(row.get("release_fine_label", "")).strip()
        for row in rows
        if str(row.get("release_fine_label", "")).strip()
    }
    for items in catalog.get("machine_actionable_fine_candidate_catalog", {}).values():
        for item in items if isinstance(items, list) else []:
            if not isinstance(item, dict):
                continue
            fine_label = str(item.get("release_label", "")).strip()
            if fine_label in represented_fine_labels:
                continue
            parent = str(item.get("parent_release_label", "")).strip()
            context = [
                str(item[key]).strip()
                for key in ("context_gate", "required_discriminator")
                if str(item.get(key, "")).strip()
            ]
            parent_candidate = broad_by_label.get(parent, {})
            context_candidate_id = str(
                item.get("context_evidence_candidate_id", "")
                or parent_candidate.get("candidate_id", "")
            ).strip()
            rows.append({
                **item,
                "candidate_role": "fine",
                "release_broad_label": parent,
                "release_fine_label": fine_label,
                "parent_broad_label": parent,
                "writeback_strategy": item.get(
                    "writeback_strategy", "supported_subset_with_parent_lock"
                ),
                "specificity_priority": int(item.get("specificity_priority", 70)),
                "hard_anti_families": item.get("hard_anti_families", []),
                "hard_anti_families_by_observation_unit": item.get(
                    "hard_anti_families_by_observation_unit", {}
                ),
                "soft_anti_families": item.get("soft_anti_families", []),
                "context_requirements": item.get("context_requirements", context),
                # A fine identity must be generated from its own discriminator,
                # not from the parent program alone.  The parent family remains
                # mandatory at group validation so a state-like discriminator
                # cannot create a foreign broad lineage.
                "required_positive_families": item.get(
                    "required_positive_families",
                    ["parent_identity", "fine_discriminator"],
                ),
                "seed_required_positive_families": item.get(
                    "seed_required_positive_families",
                    ["fine_discriminator"],
                ),
                "formal_context_evidence_required": bool(
                    item.get("formal_context_evidence_required")
                    or item.get("context_gate")
                    or parent_candidate.get("formal_context_evidence_required")
                ),
                "context_evidence_candidate_id": context_candidate_id,
                "review_required": True,
                "parent_broad_writeback_strategy": str(
                    parent_candidate.get("writeback_strategy", "")
                ),
            })
    result = {}
    for row in rows:
        candidate_id = str(row.get("candidate_id", "")).strip()
        if not candidate_id or candidate_id in result:
            raise ValueError("candidate catalog has an empty or duplicate candidate_id")
        result[candidate_id] = row
    if not result:
        raise ValueError("candidate catalog is empty")
    broad_strategy = {
        str(candidate.get("release_broad_label", "")): str(
            candidate.get("writeback_strategy", "")
        )
        for candidate in result.values()
        if str(candidate.get("candidate_role", "")) == "broad"
        and not str(candidate.get("release_fine_label", ""))
        and str(candidate.get("release_broad_label", ""))
    }
    for candidate in result.values():
        if str(candidate.get("candidate_role", "")) == "fine":
            candidate.setdefault(
                "parent_broad_writeback_strategy",
                broad_strategy.get(
                    str(candidate.get("release_broad_label", "")), ""
                ),
            )
    return result


def candidate_can_release(candidate: dict, context_ok: bool = True) -> bool:
    if candidate.get("formal_context_evidence_required") and not truth(
        candidate.get("_context_ok")
    ):
        context_ok = False
    return (
        str(candidate.get("candidate_role", "")).lower() in RELEASE_ROLES
        and bool(str(candidate.get("release_broad_label", "")).strip())
        and context_ok
    )


def candidate_can_support_broad_review(candidate: dict) -> bool:
    """Return whether a candidate may reconstruct a broad identity.

    Ordinary fine/state programs are parent-locked and cannot become broad
    recall evidence.  A fine identity participates only when the catalog
    explicitly declares a valid route to the broad parent, such as
    Lymphatic endothelial to Endothelial.
    """
    return (
        str(candidate.get("candidate_role", "")).lower() == "broad"
        or candidate.get("parent_broad_reconstruction_allowed") is True
    )


def apply_candidate_context(
    candidates: dict[str, dict],
    context_rows: Iterable[dict[str, object]] = (),
) -> dict[str, dict[str, str]]:
    """Attach one authoritative evaluation-permission state to every candidate.

    Exogenous context can only make a context-gated candidate evaluable.  It
    never contributes identity evidence.  Missing, conflicting or explicitly
    non-supported context therefore keeps the candidate visible to descriptive
    scans while removing every formal release path.
    """
    by_id: dict[str, dict[str, str]] = {}
    for row in context_rows:
        candidate_id = str(row.get("candidate_id", "")).strip()
        if not candidate_id:
            raise ValueError("context evidence contains an empty candidate_id")
        normalized = {
            "status": str(row.get("status", "")).strip().lower(),
            "reason": str(row.get("reason", "")).strip(),
            "observed_value": str(row.get("observed_value", "")).strip(),
        }
        if normalized["status"] not in CONTEXT_ALLOWED_STATUSES:
            raise ValueError(
                f"context evidence has an invalid status for {candidate_id}: "
                f"{normalized['status'] or '<empty>'}"
            )
        previous = by_id.get(candidate_id)
        if previous is not None and previous != normalized:
            raise ValueError(
                f"context evidence contains conflicting rows for {candidate_id}"
            )
        by_id[candidate_id] = normalized
    accepted_evidence_ids = {
        str(candidate.get("context_evidence_candidate_id", "") or candidate_id).strip()
        for candidate_id, candidate in candidates.items()
    }
    unknown_evidence_ids = sorted(set(by_id) - accepted_evidence_ids)
    if unknown_evidence_ids:
        raise ValueError(
            "context evidence names candidates outside the bound catalog: "
            + ", ".join(unknown_evidence_ids)
        )
    result: dict[str, dict[str, str]] = {}
    for candidate_id, candidate in candidates.items():
        evidence_id = str(
            candidate.get("context_evidence_candidate_id", "") or candidate_id
        ).strip()
        record = by_id.get(evidence_id, {})
        required = bool(candidate.get("formal_context_evidence_required"))
        status = str(record.get("status", ""))
        candidate["_context_ok"] = bool(
            not required or status in CONTEXT_SUPPORTED_STATUSES
        )
        candidate["_context_status"] = (
            "not_required" if not required else status or "not_evaluable"
        )
        candidate["_context_reason"] = str(record.get("reason", ""))
        candidate["_context_evidence_id"] = evidence_id
        result[candidate_id] = {
            "candidate_id": candidate_id,
            "context_evidence_candidate_id": evidence_id,
            "status": str(candidate["_context_status"]),
            "release_eligible": str(candidate_can_release(candidate)).lower(),
            "reason": str(candidate["_context_reason"]),
        }
    return result


def minimum_identity_core_fraction(candidate: dict) -> float:
    """Candidate-specific group prevalence needed for an ordinary identity call."""
    return max(
        0.03,
        min(0.50, number(candidate.get("minimum_identity_core_fraction"), 0.03)),
    )


def group_identity_core_fraction(row: dict[str, str]) -> float:
    return clamp(number(
        row.get("observation_identity_core_fraction"),
        number(row.get("observation_seed_fraction")),
    ))


def group_identity_core_direct_fraction(row: dict[str, str]) -> float:
    return clamp(number(
        row.get("observation_identity_core_direct_fraction"),
        number(row.get("observation_seed_fraction")),
    ))


def group_release_supported_fraction(
    row: dict[str, str], candidate: dict | None = None
) -> float:
    """Support used for whole-subcluster inheritance after an identity core exists.

    Oocyte is the deliberate exception: a canonical cluster can contain a
    sparse non-ZP/maternal identity core while the broader coherent cluster is
    inherited after somatic contradiction review.
    """
    candidate = candidate or {}
    if str(candidate.get("writeback_strategy", "")) == "canonical_cluster_membership":
        return clamp(number(row.get("observation_coherent_fraction")))
    return group_identity_core_fraction(row)


def group_supported_family_count(row: dict[str, str]) -> int:
    explicit = row.get("group_positive_family_supported_count")
    if explicit not in (None, ""):
        return int(number(explicit))
    # Compatibility for pre-v2.2 fixtures. Canonical scorer output always
    # supplies the explicit group-level family count.
    if number(row.get("observation_release_family_coherent_fraction")) >= 0.03:
        return max(2, int(number(row.get("available_positive_family_count"))))
    return 2 if number(row.get("observation_coherent_fraction")) >= 0.03 else 0


def effective_broad_writeback_strategy(candidate: dict) -> str:
    """Return the parent broad strategy that constrains a fine challenger."""
    return str(
        candidate.get("parent_broad_writeback_strategy", "")
        or candidate.get("writeback_strategy", "")
    )


def canonical_cluster_challenger(
    row: dict[str, str], candidate: dict | None = None
) -> bool:
    """Keep a sparse canonical Oocyte-like cluster visible for local recovery.

    This is deliberately not a whole-subcluster release rule.  It only permits
    a canonical-cluster candidate with strong direct/core, DEG and multi-family
    evidence to enter competition when ordinary family-prevalence AND is
    diluted by somatic cellbins.
    """
    candidate = candidate or {}
    if str(candidate.get("writeback_strategy", "")) != "canonical_cluster_membership":
        return False
    marker_deg = number(row.get("marker_deg_log2fc_mean"))
    anti_deg = number(row.get("anti_marker_deg_log2fc_mean"))
    return (
        number(row.get("available_positive_family_count")) >= 2
        and group_supported_family_count(row) >= 2
        and group_identity_core_fraction(row)
        >= max(0.10, minimum_identity_core_fraction(candidate))
        and group_identity_core_direct_fraction(row) >= 0.10
        and number(row.get("observation_release_family_coherent_fraction"))
        >= 0.02
        and number(row.get("positive_marker_detection_fraction")) >= 0.25
        and marker_deg >= 1.50
        and (anti_deg <= 0 or marker_deg - anti_deg >= 1.00)
    )


def group_candidate_detected(row: dict[str, str], candidate: dict | None = None) -> bool:
    """Detect an identity program without treating local-only coherence as identity.

    Ordinary candidates need their declared identity-core prevalence. A rare
    program may remain visible below that prevalence only with positive DEG
    contrast. Both paths require two marker families supported across the
    group, so a shared family such as contractile support cannot create a
    competing lineage by itself.
    """
    candidate = candidate or {}
    # Canonical-cluster identities (currently Oocyte) are intentionally not
    # allowed to fall through to the ordinary low-prevalence detector.  The
    # latter can keep a weak/ambient program visible as watch, but it is not
    # sufficient to claim that a zero census conflicts with a real canonical
    # cluster.  Release, completeness and absence review must all use the same
    # complete canonical-cluster rule.
    if str(candidate.get("writeback_strategy", "")) == "canonical_cluster_membership":
        return canonical_cluster_challenger(row, candidate)
    if number(row.get("available_positive_family_count")) < 2:
        return False
    required_pass = row.get("group_required_positive_families_pass")
    if required_pass not in (None, "") and not truth(required_pass):
        return False
    core = group_identity_core_fraction(row)
    family_count = group_supported_family_count(row)
    detection = number(row.get("positive_marker_detection_fraction"))
    program = number(row.get("mean_program_score"))
    marker_deg = number(row.get("marker_deg_log2fc_mean"))
    anti_deg = max(0.0, number(row.get("anti_marker_deg_log2fc_mean")))
    if family_count < 2 or detection < 0.05 or program < 0.02:
        return False
    ordinary = core >= minimum_identity_core_fraction(candidate)
    rare = core >= 0.005 and marker_deg >= 0.50 and marker_deg - anti_deg >= 0.25
    return ordinary or rare


def group_pseudobulk_contrast(row: dict[str, str]) -> float:
    n = max(1.0, number(row.get("n_observations"), 1.0))
    positive = max(0.0, number(row.get("positive_marker_pseudobulk_sum"))) / n
    anti = max(0.0, number(row.get("anti_marker_pseudobulk_sum"))) / n
    if positive <= 0:
        return 0.0
    return clamp(positive / (positive + anti + 1e-12))


def group_orthogonal_support_count(
    row: dict[str, str], candidate: dict | None = None
) -> int:
    candidate = candidate or {}
    core_min = minimum_identity_core_fraction(candidate)
    marker_deg = number(row.get("marker_deg_log2fc_mean"))
    anti_deg = max(0.0, number(row.get("anti_marker_deg_log2fc_mean")))
    return sum((
        group_supported_family_count(row) >= 2,
        group_identity_core_direct_fraction(row) >= max(0.01, core_min / 4.0),
        marker_deg - anti_deg >= 0.25,
        group_pseudobulk_contrast(row) >= 0.55
        and number(row.get("positive_marker_detection_fraction")) >= 0.05,
        number(row.get("cross_resolution_stable_fraction")) >= 0.50,
    ))


def independent_group_program(
    row: dict[str, str],
    candidate: dict | None = None,
    *,
    maximum_contradiction_fraction: float = _WRITEBACK_DEFAULTS[
        "maximum_contradiction_fraction"
    ],
) -> bool:
    """Distinguish an independent identity from a visible shared/remainder program."""
    candidate = candidate or {}
    candidate_id = str(candidate.get("candidate_id", row.get("candidate_id", "")))
    if candidate_id in GENERIC_REMAINDER_IDS:
        return False
    if not group_candidate_detected(row, candidate):
        return False
    if group_orthogonal_support_count(row, candidate) < 3:
        return False
    if canonical_cluster_challenger(row, candidate):
        return True
    contradiction = number(row.get("hard_contradiction_fraction"))
    core = group_identity_core_fraction(row)
    direct = group_identity_core_direct_fraction(row)
    marker_deg = number(row.get("marker_deg_log2fc_mean"))
    anti_deg = max(0.0, number(row.get("anti_marker_deg_log2fc_mean")))
    strong_specific_program = (
        group_supported_family_count(row) >= 2
        and core >= max(0.05, minimum_identity_core_fraction(candidate))
        and direct >= max(0.02, minimum_identity_core_fraction(candidate) / 4.0)
        and marker_deg - anti_deg >= 0.75
    )
    return (
        contradiction <= maximum_contradiction_fraction
        or strong_specific_program
    )


def contextual_parent_override_labels(candidate: dict | None = None) -> set[str]:
    """Return broad parents that absorb a non-separable lineage-origin program.

    These relationships do not suppress a genuinely separable competing cell
    identity.  They only prevent a shared/contextual program (for example a
    Theca lineage-of-origin trace inside a coherent Luteal compartment) from
    making the entire second-round subcluster look mixed.
    """
    candidate = candidate or {}
    return {
        str(rule.get("writeback_broad_label", "") or rule.get(
            "context_broad_label", ""
        )).strip()
        for rule in candidate.get("contextual_parent_overrides", [])
        if isinstance(rule, dict)
        and str(rule.get("writeback_broad_label", "") or rule.get(
            "context_broad_label", ""
        )).strip()
    }


def observation_direct_identity_seed(
    row: dict[str, str] | None,
    candidate: dict | None = None,
) -> bool:
    """Require a direct, multi-family identity core for mixedness testing.

    Candidate visibility, local smoothing and a shared marker family are not
    sufficient.  This deliberately uses the scorer's frozen direct evidence;
    it is a separability probe and never releases an individual observation.
    """
    row = row or {}
    candidate = candidate or {}
    if str(candidate.get("candidate_id", row.get("candidate_id", ""))) in GENERIC_REMAINDER_IDS:
        return False
    return (
        truth(row.get("identity_core_direct"))
        and truth(row.get("release_family_coherent"))
        and int(number(row.get("positive_family_count"))) >= 2
        and int(number(row.get("positive_gene_count"))) >= 2
        and number(row.get("direct_signal")) >= 0.03
        and not truth(row.get("ambient_suspect"))
        and not hard_contradiction(row)
    )


def direct_identity_members(
    members: Iterable[str],
    candidate_id: str,
    score_index: dict[tuple[str, str], dict[str, str]],
    candidates: dict[str, dict],
) -> set[str]:
    candidate = candidates.get(candidate_id, {})
    return {
        str(cell_id)
        for cell_id in members
        if observation_direct_identity_seed(
            score_index.get((str(cell_id), candidate_id)), candidate
        )
    }


def minimum_exclusive_component_members(n_observations: int) -> int:
    return max(
        5,
        int(math.ceil(
            max(0, n_observations)
            * _LOCAL_SPLIT_DEFAULTS["pairwise_min_exclusive_direct_fraction"]
        )),
    )


def pairwise_separable_identity_components(
    members: Iterable[str],
    left_candidate_id: str,
    right_candidate_id: str,
    score_index: dict[tuple[str, str], dict[str, str]],
    candidates: dict[str, dict],
) -> dict[str, object]:
    """Test whether two broad programs contain mutually exclusive direct cores.

    A pair is separable only when both sides contain a material exclusive
    component.  Near-total nesting therefore records coexpression/background
    rather than routing the complete subcluster into observation-level split.
    """
    member_ids = [str(value) for value in members]
    left = candidates.get(left_candidate_id, {})
    right = candidates.get(right_candidate_id, {})
    left_label = str(left.get("release_broad_label", ""))
    right_label = str(right.get("release_broad_label", ""))
    minimum = minimum_exclusive_component_members(len(member_ids))
    result: dict[str, object] = {
        "left_candidate_id": left_candidate_id,
        "right_candidate_id": right_candidate_id,
        "n_observations": len(member_ids),
        "minimum_exclusive_members": minimum,
        "left_direct_n": 0,
        "right_direct_n": 0,
        "left_only_n": 0,
        "right_only_n": 0,
        "both_n": 0,
        "separable": False,
        "reason": "",
    }
    if (
        not left_label
        or not right_label
        or left_label == right_label
        or left_candidate_id in GENERIC_REMAINDER_IDS
        or right_candidate_id in GENERIC_REMAINDER_IDS
    ):
        result["reason"] = "not_two_independent_broad_identities"
        return result
    left_members = direct_identity_members(
        member_ids, left_candidate_id, score_index, candidates
    )
    right_members = direct_identity_members(
        member_ids, right_candidate_id, score_index, candidates
    )
    left_only = left_members - right_members
    right_only = right_members - left_members
    both = left_members & right_members
    result.update({
        "left_direct_n": len(left_members),
        "right_direct_n": len(right_members),
        "left_only_n": len(left_only),
        "right_only_n": len(right_only),
        "both_n": len(both),
    })
    if len(left_only) >= minimum and len(right_only) >= minimum:
        result["separable"] = True
        result["reason"] = "mutually_exclusive_direct_identity_components"
    elif left_members and right_members:
        result["reason"] = "coexpressed_or_nested_direct_identity"
    else:
        result["reason"] = "one_or_both_direct_identity_components_absent"
    return result


def specific_component_embedded_in_generic_parent(
    members: Iterable[str],
    candidate_id: str,
    score_index: dict[tuple[str, str], dict[str, str]],
    candidates: dict[str, dict],
) -> dict[str, object]:
    """Detect a bounded specific component inside a generic remainder parent."""
    member_ids = [str(value) for value in members]
    direct = direct_identity_members(
        member_ids, candidate_id, score_index, candidates
    )
    minimum = minimum_exclusive_component_members(len(member_ids))
    complement_n = len(member_ids) - len(direct)
    return {
        "candidate_id": candidate_id,
        "direct_identity_n": len(direct),
        "complement_n": complement_n,
        "minimum_component_members": minimum,
        "separable": len(direct) >= minimum and complement_n >= minimum,
    }


def local_split_worthy_group_program(
    row: dict[str, str], candidate: dict | None = None
) -> bool:
    """Detect a specific program eligible for a bounded separability check.

    A mixed second-round subcluster can have a high contradiction fraction by
    construction: mutually exclusive identities are aggregated before the
    observation-level components are separated.  Contradiction therefore
    validates a proposed local subset, but it must not erase the candidate
    program that triggers that check.  Likewise, whole-subcluster DEG can be
    negative when a real minority lineage is diluted by the majority parent.
    This detector therefore requires multi-family, direct identity-core,
    cross-resolution and orthogonal group support, while leaving actual
    separability and contradiction to candidate-local component construction
    and validation.  A positive result is candidate visibility, not a mixed
    subcluster trigger and not a label.  The adjudicator must additionally
    prove either pairwise-exclusive direct identity components or a bounded
    specific component inside a supported generic remainder.
    """
    candidate = candidate or {}
    candidate_id = str(candidate.get("candidate_id", row.get("candidate_id", "")))
    if candidate_id in GENERIC_REMAINDER_IDS:
        return False
    if not group_candidate_detected(row, candidate):
        return False
    if group_orthogonal_support_count(row, candidate) < 3:
        return False
    if canonical_cluster_challenger(row, candidate):
        return True
    core_min = minimum_identity_core_fraction(candidate)
    core = group_identity_core_fraction(row)
    direct = group_identity_core_direct_fraction(row)
    seed = clamp(number(row.get("observation_seed_fraction")))
    release = clamp(number(row.get("observation_release_family_coherent_fraction")))
    stable = number(row.get("cross_resolution_stable_fraction")) >= _LOCAL_SPLIT_DEFAULTS[
        "cross_resolution_stability_minimum"
    ]
    required_families = [
        str(value) for value in candidate.get("required_positive_families", [])
        if str(value)
    ]
    identity_prevalent = (
        core >= max(0.30, 2.0 * core_min)
        and seed >= max(0.05, core_min)
    )
    marker_deg = number(row.get("marker_deg_log2fc_mean"))
    anti_deg = max(0.0, number(row.get("anti_marker_deg_log2fc_mean")))
    ordinary = (
        group_supported_family_count(row) >= 2
        and core >= max(0.05, core_min)
        and direct >= max(0.02, core_min / 4.0)
        and release >= 0.02
        and seed >= max(0.01, core_min / 2.0)
        and (len(required_families) >= 2 or identity_prevalent)
        and stable
    )
    material_minority = (
        group_supported_family_count(row) >= 2
        and core >= _LOCAL_SPLIT_DEFAULTS[
            "minority_min_identity_core_fraction"
        ]
        and direct >= _LOCAL_SPLIT_DEFAULTS[
            "minority_min_direct_core_fraction"
        ]
        and release >= _LOCAL_SPLIT_DEFAULTS[
            "minority_min_release_family_fraction"
        ]
        and seed >= _LOCAL_SPLIT_DEFAULTS["minority_min_seed_fraction"]
        and marker_deg - anti_deg >= _LOCAL_SPLIT_DEFAULTS[
            "minority_min_deg_contrast"
        ]
        and stable
    )
    return ordinary or material_minority


def rare_group_program_watch(
    row: dict[str, str], candidate: dict | None = None,
) -> bool:
    """Record a reproducible sub-material program without opening P41.

    The watch remains available to zero-census and per-broad recall review,
    but a 0.5%-level aggregate trace can no longer route the entire second-
    round subcluster through an expensive observation split.
    """
    candidate = candidate or {}
    if str(candidate.get("candidate_id", "")) in GENERIC_REMAINDER_IDS:
        return False
    core = group_identity_core_fraction(row)
    return (
        group_candidate_detected(row, candidate)
        and group_supported_family_count(row) >= 2
        and core >= _LOCAL_SPLIT_DEFAULTS["rare_watch_min_fraction"]
        and core < _LOCAL_SPLIT_DEFAULTS["minority_min_identity_core_fraction"]
        and group_identity_core_direct_fraction(row)
        >= _LOCAL_SPLIT_DEFAULTS["rare_watch_min_fraction"]
        and number(row.get("cross_resolution_stable_fraction"))
        >= _LOCAL_SPLIT_DEFAULTS["cross_resolution_stability_minimum"]
    )


def dominant_generic_remainder_group(
    row: dict[str, str], candidate: dict | None = None
) -> bool:
    """Recognize a coherent generic parent without making it a challenger.

    Generic Stromal/ECM is allowed to inherit a whole group only after the
    caller has established that no separable specific identity remains. Its
    aggregate contradiction is not reused as a veto: independent specific
    blockers own that decision. The parent itself must instead be nearly
    universal, direct, multigene, stable and spatially connected.
    """
    candidate = candidate or {}
    candidate_id = str(candidate.get("candidate_id", row.get("candidate_id", "")))
    if candidate_id not in GENERIC_REMAINDER_IDS:
        return False
    required_pass = row.get("group_required_positive_families_pass")
    return (
        group_candidate_detected(row, candidate)
        and group_supported_family_count(row) >= 2
        and (required_pass in (None, "") or truth(required_pass))
        and group_identity_core_fraction(row) >= 0.80
        and group_identity_core_direct_fraction(row) >= 0.80
        and number(row.get("group_positive_family_mean_fraction")) >= 0.70
        and number(row.get("positive_marker_detection_fraction")) >= 0.70
        and number(row.get("mean_program_score")) >= 0.15
        and number(row.get("marker_deg_log2fc_mean")) >= -0.75
        and number(row.get("cross_resolution_stable_fraction")) >= 0.50
        and number(row.get("spatial_group_connectivity_fraction")) >= 0.10
        and group_orthogonal_support_count(row, candidate) >= 3
    )


def group_candidate_score(row: dict[str, str], candidate: dict | None = None) -> float:
    """Rank group identities using identity-grade rather than any-family signal."""
    candidate = candidate or {}
    core = group_identity_core_fraction(row)
    direct = group_identity_core_direct_fraction(row)
    release = clamp(number(
        row.get("observation_release_family_coherent_fraction"), core
    ))
    family_mean = clamp(number(
        row.get("group_positive_family_mean_fraction"), release
    ))
    marker_deg = number(row.get("marker_deg_log2fc_mean"))
    anti_deg = max(0.0, number(row.get("anti_marker_deg_log2fc_mean")))
    deg_contrast = clamp((marker_deg - anti_deg + 0.25) / 1.50)
    return (
        0.35 * core
        + 0.15 * direct
        + 0.10 * release
        + 0.10 * family_mean
        + 0.10 * deg_contrast
        + 0.10 * group_pseudobulk_contrast(row)
        + 0.05 * clamp(number(row.get("cross_resolution_stable_fraction")))
        + 0.05 * clamp(number(row.get("mean_program_score")))
        - 0.15 * clamp(number(row.get("hard_contradiction_fraction")))
    )


def positive_family_names(row: dict[str, str]) -> set[str]:
    return {
        value.strip()
        for value in str(row.get("positive_families", "")).split(";")
        if value.strip()
    }


def candidate_core_seed(row: dict[str, str], candidate: dict) -> bool:
    """Return a candidate-identity seed, never a generic propagation bridge.

    Fine identities default to their ``fine_discriminator`` family.  Broad
    candidates may declare one or more identity families.  These seeds create
    spatial components; neighboring support cells may be inherited only in a
    bounded expansion and can never connect two identity cores transitively.
    """
    required = {
        str(value).strip()
        for value in candidate.get("seed_required_positive_families", [])
        if str(value).strip()
    }
    explicit_core = row.get("identity_core_coherent")
    base_core = (
        truth(explicit_core)
        if explicit_core not in (None, "")
        else candidate_anchor(row)
    )
    if not base_core:
        return False
    if required and not positive_family_names(row).intersection(required):
        return False
    return (
        truth(row.get("candidate_seed"))
        or number(row.get("program_score"), -math.inf) >= 0.02
        or number(row.get("direct_signal"), 0.0) > 0.0
    )


def evidence_present(row: dict[str, str]) -> bool:
    return (
        int(number(row.get("positive_family_count"), 0)) > 0
        or int(number(row.get("positive_gene_count"), 0)) > 0
        or number(row.get("direct_signal"), 0) > 0
        or number(row.get("local_signal"), 0) > 0
    )


def hard_contradiction(row: dict[str, str]) -> bool:
    """Only coherent multi-gene direct anti evidence is a hard veto."""
    family_value = row.get("direct_anti_family_count")
    family_count = (
        int(number(family_value, 0))
        if family_value not in (None, "")
        else 2
    )
    return (
        truth(row.get("hard_contradiction"))
        and int(number(row.get("direct_anti_gene_count"), 0)) >= 2
        and family_count >= 2
    )


def supported_seed(row: dict[str, str], minimum_score: float = 0.04) -> bool:
    """Return a release-grade support seed.

    The 0.70 subset support threshold is deliberately absent here, but formal
    support still needs at least two coherent marker families. One-family
    observations remain available to candidate-local expansion and may
    inherit a validated broad identity; they cannot validate a subset alone.
    """
    return (
        evidence_present(row)
        and truth(row.get("family_coherent"))
        and int(number(row.get("positive_family_count"), 0)) >= 2
        and truth(row.get("release_family_coherent", "true"))
        and number(row.get("program_score"), -math.inf) >= minimum_score
        and not hard_contradiction(row)
    )


def candidate_support_seed(
    row: dict[str, str], minimum_score: float = 0.04
) -> bool:
    """Permissive member support; multi-family evidence is checked on the subset."""
    explicit = row.get("candidate_seed")
    seed_ok = (
        truth(explicit)
        if explicit not in (None, "")
        else (
            truth(row.get("family_coherent"))
            and int(number(row.get("positive_family_count"), 0)) >= 1
        )
    )
    return (
        evidence_present(row)
        and seed_ok
        and number(row.get("program_score"), -math.inf) >= minimum_score
        and not hard_contradiction(row)
    )


def candidate_anchor(row: dict[str, str]) -> bool:
    """One coherent positive family may anchor a later group-level proposal."""
    return (
        evidence_present(row)
        and truth(row.get("family_coherent"))
        and int(number(row.get("positive_family_count"), 0)) >= 1
        and not hard_contradiction(row)
    )


def aggregate_score(row: dict[str, str]) -> float:
    return (
        number(row.get("marker_deg_log2fc_mean"))
        + 0.25 * number(row.get("observation_coherent_fraction"))
        + 0.25 * number(row.get("observation_seed_fraction"))
        + 0.10 * number(row.get("mean_program_score"))
        - 0.25 * max(0.0, number(row.get("anti_marker_deg_log2fc_mean")))
    )


def aggregate_program_supported(
    row: dict[str, str],
    candidate: dict,
    family_evidence: dict[str, object],
) -> bool:
    if (
        not row
        or family_evidence.get("status") != "PASS"
        or number(row.get("positive_marker_detection_fraction")) < 0.05
        or number(row.get("mean_program_score")) < 0.02
    ):
        return False
    marker_deg = number(row.get("marker_deg_log2fc_mean"))
    coherent = number(row.get("observation_coherent_fraction"))
    seeded = number(row.get("observation_seed_fraction"))
    anti_deg = number(row.get("anti_marker_deg_log2fc_mean"))
    canonical = (
        str(candidate.get("writeback_strategy", ""))
        == "canonical_cluster_membership"
        and marker_deg >= 1.50
        and coherent >= 0.25
    )
    common = marker_deg >= 0.50 and coherent >= 0.05
    seeded_program = marker_deg >= 0.25 and seeded >= 0.05 and coherent >= 0.10
    anti_compatible = anti_deg <= 0 or marker_deg - anti_deg >= 0.50
    return (canonical or common or seeded_program) and anti_compatible


def group_family_support(
    rows: Iterable[dict[str, str]],
    candidate: dict | None = None,
    minimum_prevalence: float = 0.03,
) -> dict[str, object]:
    """Require two marker families jointly across a proposed biological subset."""
    values = list(rows)
    if not values:
        return {
            "status": "FAIL",
            "supported_families": [],
            "family_prevalence": {},
        }
    counts: dict[str, int] = defaultdict(int)
    for row in values:
        families = {
            value.strip()
            for value in str(row.get("positive_families", "")).split(";")
            if value.strip()
        }
        for family in families:
            counts[family] += 1
    prevalence = {
        family: count / len(values) for family, count in sorted(counts.items())
    }
    supported = sorted(
        family for family, fraction in prevalence.items()
        if fraction >= minimum_prevalence
    )
    required = {
        str(value)
        for value in (candidate or {}).get("required_positive_families", [])
        if str(value)
    }
    passed = len(supported) >= 2 and required.issubset(supported)
    return {
        "status": "PASS" if passed else "FAIL",
        "supported_families": supported,
        "family_prevalence": prevalence,
        "required_families": sorted(required),
    }


def rank_supported(rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    """Rank candidates while making zero-signal rows ineligible to win."""
    eligible = [row for row in rows if candidate_anchor(row)]
    return sorted(
        eligible,
        key=lambda row: (
            -number(row.get("normalized_evidence"), number(row.get("program_score"), -math.inf)),
            -int(number(row.get("specificity_priority"), 0)),
            str(row.get("candidate_id", "")),
        ),
    )


def validate_subset(
    members: list[str],
    target: str,
    score_index: dict[tuple[str, str], dict[str, str]],
    candidate_ids: Iterable[str],
    *,
    catalog: dict[str, dict] | None = None,
    release_level: str = "all",
    aggregate_evidence: dict[str, str] | None = None,
    minimum_supported_fraction: float = _WRITEBACK_DEFAULTS[
        "supported_subset_min_lineage_supported_fraction"
    ],
    minimum_margin: float = _WRITEBACK_DEFAULTS[
        "supported_subset_min_purity_margin"
    ],
    maximum_contradiction_fraction: float = _WRITEBACK_DEFAULTS[
        "maximum_contradiction_fraction"
    ],
) -> dict[str, object]:
    if not members:
        return {"status": "FAIL", "reason": "empty_subset"}
    target_rows = [score_index[(cell, target)] for cell in members]
    family_evidence = group_family_support(
        target_rows, (catalog or {}).get(target, {})
    )
    if aggregate_evidence:
        supported_fraction = group_release_supported_fraction(
            aggregate_evidence, (catalog or {}).get(target, {})
        )
        if (
            number(aggregate_evidence.get("marker_deg_log2fc_mean")) >= 1.50
            and number(aggregate_evidence.get("anti_marker_deg_log2fc_mean")) <= 0
        ):
            contradiction_fraction = 0.0
        else:
            contradiction_fraction = number(
                aggregate_evidence.get("hard_contradiction_fraction")
            )
    else:
        supported_fraction = sum(
            candidate_anchor(row) for row in target_rows
        ) / len(members)
        contradiction_fraction = (
            sum(hard_contradiction(row) for row in target_rows) / len(members)
        )
    competitor_fractions: dict[str, float] = {}
    target_meta = (catalog or {}).get(target, {})
    target_broad = str(target_meta.get("release_broad_label", ""))
    target_role = str(target_meta.get("candidate_role", ""))
    parent_candidate_id = ""
    parent_family_evidence: dict[str, object] = {
        "status": "NOT_REQUIRED", "supported_families": [],
    }
    parent_supported_fraction = 1.0
    parent_contradiction_fraction = 0.0
    parent_identity_status = "NOT_REQUIRED"
    if target_role == "fine" and catalog:
        proposed_parent = str(
            target_meta.get("context_evidence_candidate_id", "")
        ).strip()
        parent_meta = catalog.get(proposed_parent)
        if (
            proposed_parent
            and parent_meta
            and str(parent_meta.get("candidate_role", "")) == "broad"
            and str(parent_meta.get("release_broad_label", "")) == target_broad
        ):
            parent_candidate_id = proposed_parent
            parent_rows = [
                score_index[(cell, parent_candidate_id)] for cell in members
            ]
            parent_family_evidence = group_family_support(
                parent_rows, parent_meta
            )
            parent_supported_fraction = sum(
                candidate_anchor(row) for row in parent_rows
            ) / len(members)
            parent_contradiction_fraction = sum(
                hard_contradiction(row) for row in parent_rows
            ) / len(members)
            parent_required_families = list(
                parent_meta.get("required_positive_families") or []
            )
            parent_family_floor_pass = (
                parent_family_evidence["status"] == "PASS"
                if parent_required_families
                else bool(parent_family_evidence["supported_families"])
            )
            parent_identity_status = (
                "PASS"
                if (
                    parent_family_floor_pass
                    and parent_supported_fraction >= 0.25
                    and parent_contradiction_fraction
                    <= maximum_contradiction_fraction
                )
                else "FAIL"
            )
    for candidate in candidate_ids:
        if candidate == target:
            continue
        if catalog:
            other = catalog[candidate]
            other_broad = str(other.get("release_broad_label", ""))
            other_role = str(other.get("candidate_role", ""))
            if other_broad == target_broad:
                if release_level == "broad":
                    continue
                if target_role != "fine" or other_role != "fine":
                    continue
        values = [score_index.get((cell, candidate)) for cell in members]
        competitor_rows = [row for row in values if row is not None]
        competitor_family = group_family_support(
            competitor_rows, (catalog or {}).get(candidate, {})
        )
        competitor_fractions[candidate] = (
            sum(
                candidate_anchor(row)
                and number(row.get("normalized_evidence"))
                >= number(target_rows[index].get("normalized_evidence")) + 0.05
                for index, row in enumerate(competitor_rows)
            ) / len(members)
            if competitor_family["status"] == "PASS"
            else 0.0
        )
    competitor = max(competitor_fractions.values(), default=0.0)
    competitor_id = max(competitor_fractions, key=competitor_fractions.get, default="")
    margin = supported_fraction - competitor
    passed = (
        supported_fraction >= minimum_supported_fraction
        and margin >= minimum_margin
        and contradiction_fraction <= maximum_contradiction_fraction
        and family_evidence["status"] == "PASS"
        and parent_identity_status != "FAIL"
        and (
            not aggregate_evidence
            or aggregate_program_supported(
                aggregate_evidence, (catalog or {}).get(target, {}),
                family_evidence,
            )
        )
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "lineage_supported_fraction": supported_fraction,
        "strongest_competing_candidate": competitor_id,
        "strongest_competing_fraction": competitor,
        "support_margin": margin,
        "contradiction_fraction": contradiction_fraction,
        "family_support_status": family_evidence["status"],
        "supported_families": ";".join(family_evidence["supported_families"]),
        "parent_candidate_id": parent_candidate_id,
        "parent_identity_status": parent_identity_status,
        "parent_lineage_supported_fraction": parent_supported_fraction,
        "parent_contradiction_fraction": parent_contradiction_fraction,
        "parent_family_support_status": parent_family_evidence["status"],
        "parent_supported_families": ";".join(
            parent_family_evidence.get("supported_families", [])
        ),
    }


def validate_canonical_identity_component(
    members: list[str],
    source_members: list[str],
    target: str,
    score_index: dict[tuple[str, str], dict[str, str]],
    candidate: dict,
    aggregate_evidence: dict[str, str] | None,
) -> dict[str, object]:
    """Recompute the bounded canonical-component exception from raw evidence."""
    if len(members) < 5 or not set(members).issubset(source_members):
        return {"status": "FAIL", "reason": "invalid_canonical_component_membership"}
    if not aggregate_evidence or not canonical_cluster_challenger(
        aggregate_evidence, candidate
    ):
        return {"status": "FAIL", "reason": "canonical_cluster_challenger_absent"}
    target_rows = [score_index[(cell, target)] for cell in members]
    if not all(candidate_core_seed(row, candidate) for row in target_rows):
        return {"status": "FAIL", "reason": "membership_is_not_identity_core_only"}
    family = group_family_support(target_rows, candidate)
    member_set = set(members)
    background = [
        score_index[(cell, target)]
        for cell in source_members
        if cell not in member_set
    ]
    member_program = sum(number(row.get("program_score")) for row in target_rows) / len(target_rows)
    member_direct = sum(number(row.get("direct_signal")) for row in target_rows) / len(target_rows)
    if len(background) >= 20:
        program_delta = member_program - sum(
            number(row.get("program_score")) for row in background
        ) / len(background)
        direct_delta = member_direct - sum(
            number(row.get("direct_signal")) for row in background
        ) / len(background)
    else:
        program_delta = 0.05 if member_program >= 0.04 else 0.0
        direct_delta = 0.05 if member_direct >= 0.08 else 0.0
    spatial_fraction = sum(
        number(row.get("local_seed_fraction")) >= 0.03 for row in target_rows
    ) / len(target_rows)
    passed = (
        family["status"] == "PASS"
        and program_delta >= 0.05
        and direct_delta >= 0.05
        and spatial_fraction >= 0.40
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "reason": (
            "canonical_identity_component"
            if passed
            else "canonical_component_multichannel_support_failed"
        ),
        "lineage_supported_fraction": 1.0,
        "strongest_competing_candidate": "",
        "strongest_competing_fraction": 0.0,
        "support_margin": 1.0,
        "contradiction_fraction": sum(
            hard_contradiction(row) for row in target_rows
        ) / len(target_rows),
        "family_support_status": family["status"],
        "supported_families": ";".join(family["supported_families"]),
        "program_score_delta": program_delta,
        "direct_signal_delta": direct_delta,
        "spatially_supported_fraction": spatial_fraction,
        "validation_mode": "canonical_identity_component",
    }


def resolve_overlap(
    cell_id: str,
    candidate_ids: Iterable[str],
    score_index: dict[tuple[str, str], dict[str, str]],
    catalog: dict[str, dict],
    *,
    discriminator_margin: float = 0.15,
    release_level: str = "all",
) -> tuple[str, str]:
    """Resolve one overlap by normalized evidence, never catalog order."""
    rows = [score_index[(cell_id, candidate)] for candidate in candidate_ids]
    specific = [
        row for row in rows
        if str(row.get("candidate_id")) not in GENERIC_REMAINDER_IDS
    ]
    ranked = rank_supported(specific or rows)
    if not ranked:
        return "", "no_supported_candidate"
    if len(ranked) == 1:
        return str(ranked[0]["candidate_id"]), "single_supported_candidate"
    broad_labels = {
        str(catalog[str(row["candidate_id"])].get("release_broad_label", ""))
        for row in ranked
    }
    first = number(
        ranked[0].get("normalized_evidence"),
        number(ranked[0].get("program_score"), 0),
    )
    second = number(
        ranked[1].get("normalized_evidence"),
        number(ranked[1].get("program_score"), 0),
    )
    if len(broad_labels) == 1 and release_level == "broad":
        return str(ranked[0]["candidate_id"]), "common_broad_identity"
    if first - second < discriminator_margin:
        if len(broad_labels) == 1:
            common_broad = next(iter(broad_labels))
            broad_parent_rows = [
                score_index[(cell_id, candidate_id)]
                for candidate_id, candidate in catalog.items()
                if (cell_id, candidate_id) in score_index
                and candidate.get("release_broad_label") == common_broad
                and not candidate.get("release_fine_label")
                and candidate_can_release(candidate)
                and supported_seed(score_index[(cell_id, candidate_id)])
            ]
            parent_ranked = rank_supported(broad_parent_rows)
            if parent_ranked:
                return str(parent_ranked[0]["candidate_id"]), "common_broad_parent"
        return "", "unresolved_candidate_overlap"
    return str(ranked[0]["candidate_id"]), "pairwise_discriminator"


def group_candidate_summary(
    members: list[str],
    candidate_ids: Iterable[str],
    score_index: dict[tuple[str, str], dict[str, str]],
    catalog: dict[str, dict] | None = None,
) -> list[dict[str, object]]:
    result = []
    for candidate in candidate_ids:
        rows = [score_index[(cell, candidate)] for cell in members]
        family_evidence = group_family_support(
            rows, (catalog or {}).get(candidate, {})
        )
        supported = sum(candidate_anchor(row) for row in rows) / len(rows)
        coherent = sum(truth(row.get("family_coherent")) for row in rows) / len(rows)
        contradictions = sum(hard_contradiction(row) for row in rows) / len(rows)
        positive_scores = [
            number(row.get("program_score")) for row in rows if evidence_present(row)
        ]
        result.append({
            "candidate_id": candidate,
            "supported_fraction": supported,
            "coherent_fraction": coherent,
            "contradiction_fraction": contradictions,
            "mean_supported_score": (
                sum(positive_scores) / len(positive_scores) if positive_scores else 0.0
            ),
            "evidence_n": len(positive_scores),
            "family_support_status": family_evidence["status"],
            "supported_families": ";".join(
                family_evidence["supported_families"]
            ),
        })
    return sorted(
        result,
        key=lambda row: (
            -number(row["supported_fraction"]),
            -number(row["mean_supported_score"]),
            str(row["candidate_id"]),
        ),
    )


def choose_group_parent(
    members: list[str],
    candidate_ids: Iterable[str],
    score_index: dict[tuple[str, str], dict[str, str]],
    catalog: dict[str, dict],
    *,
    preferred_parent: str = "",
    release_level: str = "all",
    blocker_candidate_ids: set[str] | None = None,
    parent_candidate_ids: set[str] | None = None,
    certified_generic_parent_ids: set[str] | None = None,
    minimum_supported_fraction: float = _WRITEBACK_DEFAULTS[
        "supported_subset_min_lineage_supported_fraction"
    ],
    minimum_margin: float = _WRITEBACK_DEFAULTS[
        "supported_subset_min_purity_margin"
    ],
    maximum_contradiction_fraction: float = _WRITEBACK_DEFAULTS[
        "maximum_contradiction_fraction"
    ],
) -> tuple[str, dict[str, object]]:
    """Choose a coarse parent for an exact remainder from immutable scores."""
    certified_generic_parent_ids = set(certified_generic_parent_ids or set())
    if any(
        candidate_id not in GENERIC_REMAINDER_IDS
        for candidate_id in certified_generic_parent_ids
    ):
        raise ValueError(
            "only generic remainder candidates may bypass aggregate contradiction"
        )
    summary = group_candidate_summary(
        members, candidate_ids, score_index, catalog
    )
    releasable = [
        row for row in summary
        if candidate_can_release(catalog[row["candidate_id"]])
        and str(catalog[row["candidate_id"]].get("candidate_role", "")) == "broad"
        and not str(catalog[row["candidate_id"]].get("release_fine_label", ""))
        and (
            parent_candidate_ids is None
            or str(row["candidate_id"]) in parent_candidate_ids
        )
        and row["evidence_n"] > 0
        and row["family_support_status"] == "PASS"
        and (
            row["contradiction_fraction"] <= maximum_contradiction_fraction
            or str(row["candidate_id"]) in certified_generic_parent_ids
        )
    ]
    if not releasable:
        return "", {"status": "UNRESOLVED", "reason": "no_releasable_positive_program"}
    by_id = {str(row["candidate_id"]): row for row in releasable}
    blocker_universe = (
        set(catalog)
        if blocker_candidate_ids is None
        else set(blocker_candidate_ids)
    )
    blockers = [
        row for row in summary
        if str(row["candidate_id"]) in blocker_universe
        and row["evidence_n"] > 0
        and row["family_support_status"] == "PASS"
        and (
            blocker_candidate_ids is not None
            or row["contradiction_fraction"] <= maximum_contradiction_fraction
        )
    ]
    broad_groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in releasable:
        broad_groups[
            str(catalog[str(row["candidate_id"])].get("release_broad_label", ""))
        ].append(row)
    broad_ranked = sorted(
        (
            {
                "broad_label": broad,
                "supported_fraction": max(
                    number(row["supported_fraction"]) for row in rows
                ),
                "representative": sorted(
                    rows,
                    key=lambda row: (
                        -number(row["supported_fraction"]),
                        -number(row["mean_supported_score"]),
                        str(row["candidate_id"]),
                    ),
                )[0],
            }
            for broad, rows in broad_groups.items()
            if broad
        ),
        key=lambda row: (
            -number(row["supported_fraction"]), str(row["broad_label"])
        ),
    )
    if not broad_ranked:
        return "", {"status": "UNRESOLVED", "reason": "no_release_broad_program"}
    winner_broad = broad_ranked[0]
    winner = winner_broad["representative"]
    if preferred_parent in by_id:
        preferred = by_id[preferred_parent]
        preferred_broad = str(
            catalog[preferred_parent].get("release_broad_label", "")
        )
        embedded = [
            row for row in blockers
            if row["candidate_id"] != preferred_parent
            and str(
                catalog[str(row["candidate_id"])].get("release_broad_label", "")
            ) != preferred_broad
            and row["supported_fraction"] >= 0.25
        ]
        if (
            (
                preferred["supported_fraction"] >= minimum_supported_fraction
                or preferred_parent in certified_generic_parent_ids
            )
            and not embedded
        ):
            return preferred_parent, {
                "status": "PASS",
                "reason": (
                    "certified_generic_remainder_parent"
                    if preferred_parent in certified_generic_parent_ids
                    else "coarse_parent_remainder"
                ),
                **preferred,
            }
    competing_broads = [
        row for row in blockers
        if str(catalog[str(row["candidate_id"])].get("release_broad_label", ""))
        != str(winner_broad["broad_label"])
    ]
    embedded_competitors = [
        row for row in competing_broads
        if number(row["supported_fraction"]) >= 0.25
    ]
    if embedded_competitors:
        strongest = sorted(
            embedded_competitors,
            key=lambda row: (
                -number(row["supported_fraction"]),
                -number(row["mean_supported_score"]),
                str(row["candidate_id"]),
            ),
        )[0]
        return "", {
            "status": "UNRESOLVED",
            "reason": "embedded_competing_program",
            "winner": winner["candidate_id"],
            "winner_supported_fraction": winner["supported_fraction"],
            "embedded_candidate": strongest["candidate_id"],
            "embedded_supported_fraction": strongest["supported_fraction"],
        }
    runner_fraction = max(
        (number(row["supported_fraction"]) for row in competing_broads),
        default=0.0,
    )
    if (
        number(winner_broad["supported_fraction"]) >= minimum_supported_fraction
        and number(winner_broad["supported_fraction"]) - runner_fraction >= minimum_margin
    ):
        same_broad = broad_groups[str(winner_broad["broad_label"])]
        coarse = [
            row for row in same_broad
            if not catalog[str(row["candidate_id"])].get("release_fine_label")
            and number(row["supported_fraction"]) >= minimum_supported_fraction
        ]
        selected = sorted(
            coarse or same_broad,
            key=lambda row: (
                -number(row["supported_fraction"]),
                -number(row["mean_supported_score"]),
                str(row["candidate_id"]),
            ),
        )[0]
        return str(selected["candidate_id"]), {
            "status": "PASS",
            "reason": "group_supported_remainder",
            "runner_supported_fraction": runner_fraction,
            **selected,
        }
    return "", {
        "status": "UNRESOLVED",
        "reason": "competing_remainder_programs",
        "winner": winner["candidate_id"],
        "winner_supported_fraction": winner["supported_fraction"],
        "runner_supported_fraction": runner_fraction,
    }


def assignment_row(
    cell_id: str,
    source_boundary: str,
    source_cluster: str,
    candidate_id: str,
    candidate: dict | None,
    origin: str,
    confidence: str,
    *,
    qc_reason: str = "",
) -> dict[str, object]:
    broad = str((candidate or {}).get("release_broad_label", "") or "")
    fine = str((candidate or {}).get("release_fine_label", "") or "")
    if not candidate_id:
        state = "qc_holdout" if qc_reason != "unmodeled_stable_program" else "unknown_candidate"
    elif fine:
        state = "defined_fine"
    else:
        state = "defined_broad_only"
    return {
        "cell_id": cell_id,
        "source_boundary": source_boundary,
        "source_cluster": source_cluster,
        "candidate_id": candidate_id,
        "final_state": state,
        "final_broad_label": broad,
        "final_fine_label": fine,
        "confidence": confidence,
        "assignment_origin": origin,
        "qc_reason": qc_reason,
    }


def deterministic_membership_hash(rows: list[dict[str, object]]) -> str:
    keys = [
        "cell_id", "source_boundary", "source_cluster", "candidate_id",
        "final_state", "final_broad_label", "final_fine_label",
        "final_cell_type",
        "state_annotations", "confidence", "assignment_origin", "qc_reason",
        "unresolved_reason",
    ]
    payload = "\n".join(
        "\t".join(str(row.get(key, "")) for key in keys)
        for row in sorted(rows, key=lambda row: str(row.get("cell_id", "")))
    )
    return hashlib.sha256((payload + "\n").encode("utf-8")).hexdigest()


def deterministic_cell_id_set_hash(rows: Iterable[dict[str, object]]) -> str:
    """Hash membership identity independently of serialization and row order."""
    cell_ids = [str(row.get("cell_id", "")) for row in rows]
    if not cell_ids or "" in cell_ids or len(cell_ids) != len(set(cell_ids)):
        raise ValueError("semantic membership requires unique nonempty cell_id")
    payload = "\n".join(sorted(cell_ids)) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def deterministic_candidate_membership_hash(
    rows: list[dict[str, object]],
) -> str:
    """Hash pre-freeze proposals without implying release-label authority."""
    keys = [
        "cell_id", "source_boundary", "source_cluster", "candidate_id",
        "proposed_state", "proposed_broad_label", "confidence",
        "assignment_origin", "unresolved_reason",
    ]
    payload = "\n".join(
        "\t".join(str(row.get(key, "")) for key in keys)
        for row in sorted(rows, key=lambda row: str(row.get("cell_id", "")))
    )
    return hashlib.sha256((payload + "\n").encode("utf-8")).hexdigest()


def index_scores(rows: list[dict[str, str]]) -> tuple[
    dict[tuple[str, str], dict[str, str]],
    dict[tuple[str, str], list[str]],
    list[str],
]:
    required = {
        "cell_id", "source_boundary", "source_cluster", "candidate_id",
        "program_score", "positive_family_count", "family_coherent",
        "direct_anti_gene_count", "direct_anti_family_count",
        "hard_contradiction",
    }
    if not rows or not required.issubset(rows[0]):
        missing = sorted(required - set(rows[0] if rows else {}))
        raise ValueError("observation score table is empty or lacks: " + ", ".join(missing))
    score_index: dict[tuple[str, str], dict[str, str]] = {}
    groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    candidates: set[str] = set()
    seen_group_cells: set[tuple[str, str, str]] = set()
    for row in rows:
        cell = str(row["cell_id"])
        candidate = str(row["candidate_id"])
        key = (cell, candidate)
        if key in score_index:
            raise ValueError(f"duplicate observation/candidate score: {key}")
        score_index[key] = row
        candidates.add(candidate)
        group_key = (str(row["source_boundary"]), str(row["source_cluster"]))
        group_cell = (*group_key, cell)
        if group_cell not in seen_group_cells:
            groups[group_key].append(cell)
            seen_group_cells.add(group_cell)
    expected = candidates
    per_cell: dict[str, set[str]] = defaultdict(set)
    for cell, candidate in score_index:
        per_cell[cell].add(candidate)
    incomplete = [cell for cell, observed in per_cell.items() if observed != expected]
    if incomplete:
        raise ValueError("score table is not a complete observation × candidate product")
    for key in groups:
        groups[key] = sorted(groups[key])
    return score_index, dict(groups), sorted(candidates)


def validate_unmodeled_programs(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Keep only stable, spatially coherent, nontechnical open-world programs."""
    technical = {
        "ribosomal", "mitochondrial", "stress", "hypoxia", "cell_cycle",
        "cell-cycle", "ecm_intensity", "ecm-only", "ambient", "low_rna",
    }
    accepted = []
    for row in rows:
        resolutions = {
            item.strip() for item in str(row.get("resolutions", "")).replace(",", ";").split(";")
            if item.strip()
        }
        excluded = {
            item.strip().lower()
            for item in str(row.get("excluded_program_classes", "")).replace(",", ";").split(";")
            if item.strip()
        }
        if (
            len(resolutions) >= 2
            and truth(row.get("spatially_coherent"))
            and not excluded.intersection(technical)
            and int(number(row.get("coexpressed_gene_count"), 0)) >= 2
        ):
            accepted.append(row)
    return accepted


def fine_audit_complete(
    catalog_document: dict,
    present_broad_labels: set[str],
    audit_rows: list[dict[str, str]],
) -> tuple[bool, set[tuple[str, str]], set[tuple[str, str]]]:
    expected = {
        (str(item.get("parent_release_label", "")), str(item.get("candidate_id", "")))
        for items in catalog_document.get(
            "machine_actionable_fine_candidate_catalog", {}
        ).values()
        for item in items
        if str(item.get("parent_release_label", "")) in present_broad_labels
    }
    observed = {
        (str(row.get("parent_broad_label", "")), str(row.get("candidate_id", "")))
        for row in audit_rows
        if str(row.get("status", "")) in {"supported", "refuted", "not_evaluable"}
    }
    return expected.issubset(observed), expected, observed


def write_manifest(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
