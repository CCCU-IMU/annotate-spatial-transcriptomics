#!/usr/bin/env python3
"""Opt an existing v1.7 project into continuous lineage-signal coverage."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from init_annotation_project import REGISTRIES


NEW_REGISTRIES = ("lineage_signal_boundary_registry.tsv", "lineage_signal_registry.tsv")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    project_path = root / "config/project.json"
    if not project_path.is_file():
        raise SystemExit("missing config/project.json")
    project = json.loads(project_path.read_text(encoding="utf-8"))
    if project.get("framework_version") not in {"1.7.0", "1.8.0"}:
        raise SystemExit(f"expected framework 1.7.0 or 1.8.0, found {project.get('framework_version')!r}")
    state = root / "state"
    state.mkdir(parents=True, exist_ok=True)
    for filename in NEW_REGISTRIES:
        path = state / filename
        if not path.exists():
            with path.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle, delimiter="\t").writerow(REGISTRIES[filename])
    project.update({
        "framework_version": "1.8.0",
        "continuous_open_world_lineage_scan_required": True,
        "lineage_signal_scan_boundaries": ["whole_tissue", "broad_class_recluster", "targeted_recluster"],
        "lineage_signal_resolution_policy": "selected_plus_two_higher_available",
        "all_cell_marker_spatial_panels_required": True,
    })
    project_path.write_text(json.dumps(project, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "MIGRATED", "project_root": str(root), "framework_version": "1.8.0"}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
