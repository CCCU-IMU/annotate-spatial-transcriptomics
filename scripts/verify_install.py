#!/usr/bin/env python3
from __future__ import annotations

import ast
import json
import sys
from pathlib import Path


REQUIRED = (
    "SKILL.md",
    "agents/openai.yaml",
    "references/context-and-biology.md",
    "references/iterative-controller.md",
    "references/multi-route-controller.md",
    "references/matched-single-cell-reference.md",
    "references/multi-sample-agent-orchestration.md",
    "references/profiles/sheep_ovary_rfirst_profile.json",
    "references/profiles/sheep_ovary_literature_2025_2026.md",
    "references/profiles/sheep_ovary_rfirst_case_reference.md",
    "references/r-first-workflow.md",
    "references/taxonomy-and-pool-design.md",
    "references/quality-standard.md",
    "references/report-contract.md",
    "references/state-schema.md",
    "scripts/discover_inputs.py",
    "scripts/init_annotation_project.py",
    "scripts/autopilot_status.py",
    "scripts/plan_next_iteration.py",
    "scripts/check_completion_gate.py",
    "scripts/audit_release_taxonomy.py",
    "scripts/validate_matched_reference_crosswalk.py",
    "scripts/resolve_workflow_profile.py",
    "scripts/init_annotation_cohort.py",
    "scripts/validate_cohort_state.py",
    "scripts/migrate_project_v1_3_to_v1_4.py",
    "scripts/request_cohort_confirmation.py",
    "scripts/record_cohort_confirmation.py",
    "scripts/request_final_annotation_confirmation.py",
    "scripts/build_report.py",
    "scripts/audit_release.py",
)


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "annotate-spatial-transcriptomics").resolve()
    missing = [name for name in REQUIRED if not (root / name).is_file()]
    if missing:
        print("Installation is incomplete; missing:", file=sys.stderr)
        for name in missing:
            print(f"  - {name}", file=sys.stderr)
        return 1

    failures: list[str] = []
    for path in sorted((root / "scripts").glob("*.py")):
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except (SyntaxError, UnicodeDecodeError) as exc:
            failures.append(f"{path.relative_to(root)}: {exc}")
    for path in sorted(root.rglob("*.json")):
        try:
            json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            failures.append(f"{path.relative_to(root)}: {exc}")
    if failures:
        print("Installation validation failed:", file=sys.stderr)
        for item in failures:
            print(f"  - {item}", file=sys.stderr)
        return 1

    print(f"PASS: {root}")
    print("The Skill structure, Python syntax and JSON assets are valid.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
