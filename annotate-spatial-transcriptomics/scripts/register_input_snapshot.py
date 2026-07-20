#!/usr/bin/env python3
"""Register an immutable input file in a project's input snapshot registry."""

from __future__ import annotations

import argparse
import hashlib
from datetime import datetime, timezone
from pathlib import Path

from registry_io import locked_tsv_update


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--snapshot-id", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--path", required=True, type=Path)
    parser.add_argument("--kind", required=True)
    parser.add_argument("--status", default="frozen")
    parser.add_argument("--sha256-file", type=Path, help="scheduler-produced sha256sum record; avoids rereading a very large input")
    args = parser.parse_args()

    path = args.path.resolve()
    if not path.is_file():
        raise SystemExit(f"input does not exist or is not a file: {path}")
    if args.sha256_file:
        parts = args.sha256_file.read_text(encoding="utf-8").strip().split()
        if len(parts) < 2 or len(parts[0]) != 64 or any(ch not in "0123456789abcdefABCDEF" for ch in parts[0]):
            raise SystemExit("invalid sha256sum record")
        recorded = Path(parts[-1].lstrip("*"))
        if recorded.resolve() != path:
            raise SystemExit("sha256sum record path does not match --path")
        digest = parts[0].lower()
    else:
        digest = sha256(path)
    row = {
        "snapshot_id": args.snapshot_id,
        "sample_id": args.sample,
        "path": str(path),
        "kind": args.kind,
        "size_bytes": str(path.stat().st_size),
        "sha256": digest,
        "status": args.status,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    registry = args.project_root / "state" / "input_snapshot_registry.tsv"
    already = {"value": False}

    def mutate(rows, fields):
        if set(row) - set(fields):
            raise SystemExit(f"registry missing required columns: {sorted(set(row) - set(fields))}")
        prior = [item for item in rows if item.get("snapshot_id") == args.snapshot_id]
        if prior:
            old = prior[0]
            if old.get("sha256") != digest or old.get("path") != str(path):
                raise SystemExit("snapshot ID already exists with a different path or digest")
            already["value"] = True
            return rows
        rows.append(row)
        return rows

    locked_tsv_update(registry, mutate)
    print(args.snapshot_id, "already_registered" if already["value"] else digest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
