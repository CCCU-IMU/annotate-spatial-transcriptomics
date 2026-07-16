#!/usr/bin/env python3
"""Compute leakage-safe broad-label benchmark metrics after annotation freeze."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from validate_benchmark_isolation import validate as validate_isolation


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def metrics(truth: list[str], prediction: list[str], unresolved: str) -> dict:
    labels = sorted(set(truth) | {value for value in prediction if value and value != unresolved})
    per_class = {}
    for label in labels:
        tp = sum(t == label and p == label for t, p in zip(truth, prediction))
        fp = sum(t != label and p == label for t, p in zip(truth, prediction))
        fn = sum(t == label and p != label for t, p in zip(truth, prediction))
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        per_class[label] = {"precision": precision, "recall": recall, "f1": f1, "support": sum(t == label for t in truth)}
    return {
        "macro_f1": sum(row["f1"] for row in per_class.values()) / max(1, len(per_class)),
        "per_class": per_class,
        "unresolved_fraction": sum(not value or value == unresolved for value in prediction) / max(1, len(prediction)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--frozen-annotation", required=True, type=Path)
    parser.add_argument("--author-labels", required=True, type=Path)
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--truth-col", default="author_broad_label")
    parser.add_argument("--prediction-col", default="final_broad_label")
    parser.add_argument("--rare-labels", default="")
    parser.add_argument("--unresolved-label", default="QC")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    isolation = validate_isolation(args.project_root.resolve(), args.config.resolve(), args.frozen_annotation.resolve())
    if isolation["status"] != "PASS":
        raise SystemExit("benchmark isolation failed: " + "; ".join(isolation["errors"]))
    frozen = read(args.frozen_annotation)
    frozen_hash = hashlib.sha256(args.frozen_annotation.read_bytes()).hexdigest()
    if frozen_hash != isolation.get("frozen_annotation_sha256"):
        raise SystemExit("frozen annotation changed before author-label unblinding")
    author = read(args.author_labels)
    predicted = {row[args.cell_id_col]: row.get(args.prediction_col, "") for row in frozen}
    truth = {row[args.cell_id_col]: row.get(args.truth_col, "") for row in author}
    if set(predicted) != set(truth):
        raise SystemExit("frozen annotation and author labels do not have identical cell IDs")
    ids = sorted(truth)
    truth_values = [truth[cell_id] for cell_id in ids]
    prediction_values = [predicted[cell_id] for cell_id in ids]
    result = metrics(truth_values, prediction_values, args.unresolved_label)
    rare = {value for value in args.rare_labels.split(",") if value}
    rare_fp = {label: sum(t != label and p == label for t, p in zip(truth_values, prediction_values)) for label in sorted(rare)}
    by_sample = defaultdict(list)
    sample_truth: dict[str, list[str]] = defaultdict(list)
    sample_prediction: dict[str, list[str]] = defaultdict(list)
    sample_by_id = {row[args.cell_id_col]: row.get("sample_id", "single") for row in author}
    for cell_id in ids:
        by_sample[sample_by_id[cell_id]].append(int(truth[cell_id] == predicted[cell_id]))
        sample_truth[sample_by_id[cell_id]].append(truth[cell_id])
        sample_prediction[sample_by_id[cell_id]].append(predicted[cell_id])
    per_sample_macro_f1 = {
        sample: metrics(sample_truth[sample], sample_prediction[sample], args.unresolved_label)["macro_f1"]
        for sample in sorted(sample_truth)
    }
    stability_values = list(per_sample_macro_f1.values())
    result.update({
        "status": "PASS", "n_observations": len(ids), "rare_type_false_positives": rare_fp,
        "cross_sample_accuracy": {sample: sum(values) / len(values) for sample, values in sorted(by_sample.items())},
        "cross_sample_macro_f1": per_sample_macro_f1,
        "cross_sample_stability": {
            "minimum_macro_f1": min(stability_values) if stability_values else 0.0,
            "maximum_macro_f1": max(stability_values) if stability_values else 0.0,
            "macro_f1_range": (max(stability_values) - min(stability_values)) if stability_values else 0.0,
        },
        "frozen_annotation_sha256": frozen_hash,
        "author_labels_unblinded_after_freeze": True,
    })
    if hashlib.sha256(args.frozen_annotation.read_bytes()).hexdigest() != frozen_hash:
        raise SystemExit("frozen annotation changed during benchmark evaluation")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
