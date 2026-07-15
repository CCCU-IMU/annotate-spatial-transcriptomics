#!/usr/bin/env python3
"""Validate a formal whole-tissue or pool grid against a frozen workflow profile."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from validate_profile_role import load_profile


def parse_grid(value: str) -> list[float]:
    values = [float(item.strip()) for item in value.split(",") if item.strip()]
    if not values or len(values) != len(set(values)):
        raise ValueError("resolution grid must be nonempty and unique")
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workflow-profile", required=True, type=Path)
    parser.add_argument("--scope", required=True, choices=["whole_tissue", "pool"])
    parser.add_argument("--resolutions", required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    try:
        profile = load_profile(args.workflow_profile, "workflow_preprocessing")
        if args.scope == "whole_tissue":
            expected = profile["stereopy_cellbin_pped_contract"]["clustering"]["candidate_resolutions"]
        else:
            expected = profile["cohort_reclustering_contract"]["candidate_resolutions"]
        observed = parse_grid(args.resolutions)
        expected = [float(value) for value in expected]
        if observed != expected:
            raise ValueError(f"formal {args.scope} grid must equal {expected}; observed {observed}")
        minimum = min(expected)
        if any(value < minimum for value in observed):
            raise ValueError(f"resolution below profile minimum {minimum}")
        result = {
            "status": "PASS", "scope": args.scope, "candidate_resolutions": observed,
            "minimum_resolution": minimum, "profile_id": profile.get("profile_id"),
        }
        code = 0
    except Exception as exc:
        result = {"status": "FAIL", "scope": args.scope, "error": str(exc)}
        code = 2
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
