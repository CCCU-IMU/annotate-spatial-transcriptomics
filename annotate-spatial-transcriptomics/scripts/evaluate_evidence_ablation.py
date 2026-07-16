#!/usr/bin/env python3
"""Verify that anti-marker, spatial and broad-cohort evidence improve benchmark F1."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from run_biological_benchmark import metrics


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("fixture", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    fixture = json.loads(args.fixture.read_text(encoding="utf-8"))
    truth = fixture["truth"]
    predictions = fixture["predictions"]
    scores = {name: metrics(truth, values, fixture.get("unresolved_label", "QC"))["macro_f1"] for name, values in predictions.items()}
    full = scores.get("full_evidence", -1)
    required = ["without_anti_marker", "without_spatial", "without_broad_cohort"]
    errors = [f"missing ablation {name}" for name in required if name not in scores]
    errors.extend(f"full evidence does not outperform {name}" for name in required if name in scores and full <= scores[name])
    result = {"status": "PASS" if not errors else "FAIL", "macro_f1": scores, "errors": errors}
    if args.out:
        args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
