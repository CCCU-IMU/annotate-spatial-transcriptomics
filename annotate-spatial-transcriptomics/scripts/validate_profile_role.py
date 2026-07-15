#!/usr/bin/env python3
"""Fail closed when a workflow profile is used as a biological profile or vice versa."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REQUIRED = {
    "biological_evidence": {"lineages", "evidence_gates", "annotation_workflow_policy", "release_taxonomy"},
    "workflow_preprocessing": {"workflow", "stereopy_cellbin_pped_contract", "external_reference_policy", "release_policy"},
}


def load_profile(path: Path, expected_role: str) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    role = data.get("profile_role")
    if role != expected_role:
        raise ValueError(f"profile_role={role!r}; expected {expected_role!r}: {path}")
    missing = sorted(REQUIRED[expected_role] - set(data))
    if missing:
        raise ValueError(f"{expected_role} profile lacks required sections: {','.join(missing)}")
    return data


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("profile", type=Path)
    parser.add_argument("--expected-role", required=True, choices=sorted(REQUIRED))
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    try:
        profile = load_profile(args.profile, args.expected_role)
        result = {
            "status": "PASS",
            "profile": str(args.profile.resolve()),
            "profile_id": profile.get("profile_id"),
            "profile_role": profile.get("profile_role"),
            "profile_schema_version": profile.get("profile_schema_version"),
        }
        code = 0
    except Exception as exc:
        result = {"status": "FAIL", "profile": str(args.profile.resolve()), "error": str(exc)}
        code = 2
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
