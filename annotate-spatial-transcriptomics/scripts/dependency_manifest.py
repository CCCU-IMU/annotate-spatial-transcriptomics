#!/usr/bin/env python3
"""Content-hash dependency manifests for expensive derived assets."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest_path(target: Path) -> Path:
    return Path(str(target) + ".deps.json")


def build(target: Path, dependencies: list[Path], metadata: dict | None = None) -> dict:
    if not target.is_file():
        raise FileNotFoundError(target)
    records = []
    for dependency in dependencies:
        if not dependency.is_file():
            raise FileNotFoundError(dependency)
        records.append({"path": str(dependency.resolve()), "sha256": sha256(dependency)})
    result = {
        "status": "PASS", "target": str(target.resolve()), "target_sha256": sha256(target),
        "dependencies": records, "metadata": metadata or {},
    }
    out = manifest_path(target)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return result


def stale(target: Path, dependencies: list[Path]) -> bool:
    if not target.is_file():
        return True
    manifest = manifest_path(target)
    if not manifest.is_file():
        return True
    try:
        record = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return True
    if (
        record.get("status") != "PASS"
        or record.get("target") != str(target.resolve())
        or record.get("target_sha256") != sha256(target)
    ):
        return True
    expected = {str(path.resolve()): sha256(path) for path in dependencies if path.is_file()}
    if len(expected) != len(dependencies):
        return True
    observed = {row.get("path"): row.get("sha256") for row in record.get("dependencies", [])}
    for path_value, expected_hash in observed.items():
        path = Path(str(path_value))
        if not path.is_file() or sha256(path) != expected_hash:
            return True
    return any(observed.get(path_value) != expected_hash for path_value, expected_hash in expected.items())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("target", type=Path)
    parser.add_argument("dependencies", nargs="+", type=Path)
    args = parser.parse_args()
    print(json.dumps(build(args.target, args.dependencies), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
