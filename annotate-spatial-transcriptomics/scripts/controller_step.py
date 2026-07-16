#!/usr/bin/env python3
"""Translate fail-closed validators into Agent business states without hiding program faults."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from validate_direct_lineage_workflow import audit as audit_workflow


ARTIFACTS = {
    "iteration_plan": "provenance/iteration_plan.json",
    "completion_gate": "provenance/completion_gate.json",
    "master_quality": "provenance/master_quality_review_request.json",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--step", required=True, choices=["workflow_audit", *ARTIFACTS])
    parser.add_argument("--strict-exit-code", action="store_true")
    args = parser.parse_args()
    root = args.project_root.resolve()
    try:
        if args.step == "workflow_audit":
            source = audit_workflow(root)
        else:
            path = root / ARTIFACTS[args.step]
            if not path.is_file():
                source = {"status": "MISSING", "artifact": str(path)}
            else:
                source = json.loads(path.read_text(encoding="utf-8"))
                if not isinstance(source, dict):
                    raise ValueError("controller artifact is not a JSON object")
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError, TypeError) as exc:
        result = {"controller_state": "EXECUTION_FAILURE", "step": args.step, "error": str(exc)}
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 1
    status = str(source.get("status", "")).upper()
    if status in {"PASS", "READY_FOR_COMPLETION_AUDIT", "READY", "APPROVED"}:
        state = "CONTINUE"
    elif status in {"ITERATION_REQUIRED", "BLOCKED", "MISSING"}:
        state = "ITERATION_REQUIRED" if args.step in {"workflow_audit", "iteration_plan"} else "EXPECTED_GATE_BLOCKED"
    else:
        result = {"controller_state": "EXECUTION_FAILURE", "step": args.step, "error": f"unrecognized business status: {status!r}"}
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 1
    result = {"controller_state": state, "step": args.step, "source_status": status, "execution_failure": False}
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 2 if args.strict_exit_code and state != "CONTINUE" else 0


if __name__ == "__main__":
    raise SystemExit(main())
