#!/usr/bin/env python3
"""Write a deterministic, non-circular release manifest and SHA256 file."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from dependency_manifest import build as build_dependency_manifest


DEFAULT_ROOTS = (
    "config", "state", "tables", "figures", "spatial_nodes",
    "spatial_genes", "review", "report", "provenance",
)
EXCLUDED = {
    "provenance/release_manifest.tsv",
    "provenance/checksums.sha256",
    # The audit is written after the manifest validates its inputs, so hashing
    # it would create a circular dependency.
    "provenance/release_audit.json",
    "provenance/release_audit.json.deps.json",
    "provenance/release_manifest.tsv.deps.json",
    "provenance/checksums.sha256.deps.json",
    # The autonomous controller is intentionally run after the release audit
    # to prove that the project is terminal.  Its status/next-action outputs
    # are operational sentinels, not immutable biological release evidence.
    "provenance/autopilot_status.json",
    "state/autopilot_next_actions.tsv",
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--roots", nargs="*", default=list(DEFAULT_ROOTS))
    args = parser.parse_args()
    root = args.project_root.resolve()
    provenance = root / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    manifest = provenance / "release_manifest.tsv"
    checksums = provenance / "checksums.sha256"

    paths: list[Path] = []
    for name in args.roots:
        base = root / name
        if not base.exists():
            continue
        paths.extend(
            path for path in base.rglob("*")
            if path.is_file()
            and str(path.relative_to(root)) not in EXCLUDED
            and ".ipynb_checkpoints" not in path.parts
        )
    paths = sorted(set(paths), key=lambda path: str(path.relative_to(root)))
    rows = [
        (str(path.relative_to(root)), path.stat().st_size, digest(path))
        for path in paths
    ]

    manifest_tmp = manifest.with_suffix(manifest.suffix + ".tmp")
    checksum_tmp = checksums.with_suffix(checksums.suffix + ".tmp")
    with manifest_tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "size_bytes", "sha256"])
        writer.writerows(rows)
    with checksum_tmp.open("w", encoding="utf-8") as handle:
        for relative, _, value in rows:
            handle.write(f"{value}  {relative}\n")
    manifest_tmp.replace(manifest)
    checksum_tmp.replace(checksums)
    build_dependency_manifest(manifest, paths, {"asset_class": "release_manifest"})
    build_dependency_manifest(checksums, paths, {"asset_class": "release_checksums"})
    result = {
        "status": "PASS",
        "files": len(rows),
        "manifest": str(manifest),
        "checksums": str(checksums),
        "excluded_non_circular": sorted(EXCLUDED),
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
