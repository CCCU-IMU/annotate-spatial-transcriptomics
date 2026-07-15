#!/usr/bin/env python3
"""Validate biological context and its optional tissue profile before analysis."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


REQUIRED = {
    "species", "tissue", "developmental_or_reproductive_stage", "condition",
    "platform", "observation_unit", "primary_questions", "priority_lineages",
}


def norm(value: object) -> str:
    return " ".join(str(value).lower().replace("_", " ").split())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("context", type=Path)
    ap.add_argument("--profile", type=Path)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    ctx = json.loads(args.context.read_text(encoding="utf-8"))
    errors = [f"missing or empty context field: {key}" for key in sorted(REQUIRED) if not ctx.get(key)]
    warnings: list[str] = []
    if not isinstance(ctx.get("primary_questions", []), list): errors.append("primary_questions must be a list")
    if not isinstance(ctx.get("priority_lineages", []), list): errors.append("priority_lineages must be a list")
    if norm(ctx.get("developmental_or_reproductive_stage")) in {"unknown", "na", "n/a"}:
        warnings.append("stage is unknown; stage-specific subtype confidence must be downgraded")
    if norm(ctx.get("observation_unit")) in {"cellbin", "spot"}:
        warnings.append("observation counts are not biological-cell counts")
    profile_id = None
    if args.profile:
        profile = json.loads(args.profile.read_text(encoding="utf-8")); profile_id = profile.get("profile_id")
        if norm(profile.get("species")) != norm(ctx.get("species")): errors.append("profile species does not match context")
        if norm(profile.get("tissue")) != norm(ctx.get("tissue")): errors.append("profile tissue does not match context")
        known = set(profile.get("lineages", {})) | set(profile.get("context_specific_identity_rules", profile.get("rare_cell_rules", {})))
        for lineage in ctx.get("priority_lineages", []):
            key = norm(lineage).replace(" ", "_")
            if key not in known and not any(key in k or k in key for k in known):
                warnings.append(f"priority lineage has no exact profile block: {lineage}")
    result = {"status": "PASS" if not errors else "FAIL", "profile_id": profile_id,
              "errors": errors, "warnings": warnings, "context": ctx}
    out = args.out or args.context.with_name("biological_context_validation.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
