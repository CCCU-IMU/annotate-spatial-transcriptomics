#!/usr/bin/env python3
"""Load and validate the canonical v2.2 controller-threshold registry."""

from __future__ import annotations

import json
import math
from pathlib import Path


REGISTRY_PATH = (
    Path(__file__).resolve().parents[1]
    / "references/controller_thresholds_v2_2.json"
)

OBSERVATION_POLICY_KEYS = {
    "whole_subcluster_min_lineage_supported_fraction",
    "whole_subcluster_min_purity_margin",
    "whole_subcluster_min_raw_two_family_supported_fraction",
    "whole_subcluster_min_raw_two_family_margin",
    "whole_subcluster_embedded_competitor_raw_trigger",
    "whole_subcluster_dominant_seed_fraction",
    "whole_subcluster_dominant_direct_core_fraction",
    "whole_subcluster_dominant_identity_core_fraction",
    "whole_subcluster_dominant_max_contradiction_fraction",
    "supported_subset_min_lineage_supported_fraction",
    "supported_subset_min_purity_margin",
    "present_label_min_lineage_supported_fraction",
    "present_label_min_purity_margin",
    "maximum_contradiction_fraction",
}


def _fraction_map(section: object, label: str) -> dict[str, float]:
    if not isinstance(section, dict) or not section:
        raise ValueError(f"threshold registry lacks {label}")
    result: dict[str, float] = {}
    for key, value in section.items():
        number = float(value)
        if not math.isfinite(number) or not 0 <= number <= 1:
            raise ValueError(f"{label}.{key} is outside [0,1]")
        result[str(key)] = number
    return result


def load_controller_thresholds(path: Path | None = None) -> dict:
    source = (path or REGISTRY_PATH).resolve()
    document = json.loads(source.read_text(encoding="utf-8"))
    if document.get("schema_version") != "2.2":
        raise ValueError("controller threshold registry is not schema 2.2")

    scoring = _fraction_map(document.get("scoring_policy"), "scoring_policy")
    if abs(scoring["direct_weight"] + scoring["local_weight"] - 1.0) > 1e-12:
        raise ValueError("direct and local scoring weights must sum to one")

    observation = _fraction_map(
        document.get("observation_writeback_policy"),
        "observation_writeback_policy",
    )
    if set(observation) != OBSERVATION_POLICY_KEYS:
        missing = sorted(OBSERVATION_POLICY_KEYS - set(observation))
        extra = sorted(set(observation) - OBSERVATION_POLICY_KEYS)
        raise ValueError(
            "observation policy keys differ from the controller contract: "
            f"missing={missing}, extra={extra}"
        )

    completion = document.get("completion_policy", {})
    if not (
        0 <= float(completion.get("residual_qc_fraction_trigger", -1)) <= 1
        and int(completion.get("residual_qc_count_trigger", -1)) >= 0
    ):
        raise ValueError("invalid residual-QC completion policy")

    local_subset = document.get("local_subset_policy", {})
    if not (
        int(local_subset.get("maximum_second_subset_rounds", -1)) in {0, 1, 2}
        and int(local_subset.get("minimum_component_members", -1)) >= 2
    ):
        raise ValueError("invalid local-subset policy")
    _fraction_map(document.get("fine_release_policy"), "fine_release_policy")
    _fraction_map(document.get("state_release_policy"), "state_release_policy")

    resolution = document.get("resolution_selection", {})
    for purpose in (
        "whole_tissue_cohort_partition", "cohort_identity_resolution"
    ):
        block = resolution.get(purpose, {})
        metric_weights = _fraction_map(
            block.get("metric_weights"), f"resolution_selection.{purpose}.metric_weights"
        )
        _fraction_map(
            block.get("penalty_weights"), f"resolution_selection.{purpose}.penalty_weights"
        )
        if not 0 < sum(metric_weights.values()) <= 1.0 + 1e-12:
            raise ValueError(f"invalid resolution metric-weight total: {purpose}")
    for key in (
        "composite_equivalence_tolerance", "per_metric_equivalence_tolerance"
    ):
        value = float(resolution.get(key, -1))
        if not 0 <= value <= 1:
            raise ValueError(f"invalid resolution-selection tolerance: {key}")
    return document


def observation_writeback_defaults(path: Path | None = None) -> dict[str, float]:
    return dict(load_controller_thresholds(path)["observation_writeback_policy"])


def scoring_defaults(path: Path | None = None) -> dict[str, float]:
    return dict(load_controller_thresholds(path)["scoring_policy"])
