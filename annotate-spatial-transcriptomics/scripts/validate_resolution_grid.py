#!/usr/bin/env python3
"""Validate whole-tissue/cohort grids against the active workflow contract."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from validate_profile_role import load_profile
from workflow_contract_lib import active_workflow_contract
from workflow_contract_lib import sha256


def parse_grid(value: str) -> list[float]:
    values = [float(item.strip()) for item in value.split(",") if item.strip()]
    if not values or len(values) != len(set(values)):
        raise ValueError("resolution grid must be nonempty and unique")
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--workflow-profile", required=True, type=Path)
    parser.add_argument("--scope", required=True, choices=["whole_tissue", "cohort"])
    parser.add_argument("--resolutions", required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    try:
        root = args.project_root.resolve()
        profile = load_profile(args.workflow_profile, "workflow_preprocessing")
        contract = active_workflow_contract(root)
        if contract.get("workflow_profile_id") is None:
            raise ValueError("no active workflow profile is bound to the project")
        if not contract.get("workflow_profile_binding_valid"):
            raise ValueError("active workflow-profile binding is missing or stale")
        if args.workflow_profile.resolve() != Path(contract.get("workflow_profile", "")).resolve() or sha256(args.workflow_profile) != contract.get("workflow_profile_sha256"):
            raise ValueError("supplied workflow profile is not the hash-bound active workflow profile")
        if profile.get("profile_id") != contract.get("workflow_profile_id") and not contract["same_batch_stereopy_cellbin_rfirst"]:
            raise ValueError("active workflow profile ID differs from the supplied profile")
        if args.scope == "whole_tissue":
            expected = contract.get("whole_tissue_resolution_grid")
        else:
            expected = contract.get("cohort_resolution_grid")
        if not isinstance(expected, list) or not expected:
            raise ValueError(f"active workflow profile does not declare a {args.scope} candidate grid")
        observed = parse_grid(args.resolutions)
        expected = [float(value) for value in expected]
        if observed != expected:
            raise ValueError(f"active {args.scope} grid must equal {expected}; observed {observed}")
        result = {
            "status": "PASS", "scope": args.scope, "candidate_resolutions": observed,
            "profile_id": profile.get("profile_id"), "active_workflow_contract": contract,
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
