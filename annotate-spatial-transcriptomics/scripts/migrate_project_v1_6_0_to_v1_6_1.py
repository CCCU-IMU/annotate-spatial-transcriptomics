#!/usr/bin/env python3
"""Migrate v1.6.0 confidence/RCTD registries to canonical v1.6.1 values."""
from __future__ import annotations

import argparse
import csv
import gzip
import shutil
from datetime import datetime, timezone
from pathlib import Path

from confidence_lib import atlas_canonical, canonical


def read(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    temporary = Path(str(path) + ".tmp")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(temporary, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = root / "provenance" / "migration_backups" / f"v1.6.0_to_v1.6.1_{stamp}"
    backup.mkdir(parents=True, exist_ok=False)

    confidence_files = {
        "state/cell_ledger.tsv": ["confidence", "final_confidence"],
        "state/cell_ledger.tsv.gz": ["confidence", "final_confidence"],
        "state/cluster_decision_ledger.tsv": ["confidence"],
        "state/direct_return_registry.tsv": ["confidence", "rctd_tier"],
    }
    for relative, columns in confidence_files.items():
        path = root / relative
        rows = read(path)
        if not rows:
            continue
        shutil.copy2(path, backup / path.name)
        fields = list(rows[0])
        for row in rows:
            for column in columns:
                value = row.get(column, "").strip()
                if value:
                    row[column] = canonical(value, migration_mode=True)
        write(path, rows, fields)

    route_path = root / "state/route_attempt_registry.tsv"
    routes = read(route_path)
    if routes:
        shutil.copy2(route_path, backup / route_path.name)
        old_fields = list(routes[0])
        new_fields = [field for field in old_fields if field not in {"rctd_extreme_n", "rctd_medium_low_n"}]
        for field in ("rctd_high_n", "rctd_moderate_n", "rctd_low_n"):
            if field not in new_fields:
                new_fields.append(field)
        for row in routes:
            old_high = row.get("rctd_high_n", "0")
            row["rctd_high_n"] = row.get("rctd_extreme_n", "0")
            row["rctd_moderate_n"] = old_high
            row["rctd_low_n"] = row.get("rctd_medium_low_n", "0")
            enum = row.get("atlas_confidence_enum", "").strip()
            if enum and "|" not in enum:
                row["atlas_confidence_enum"] = atlas_canonical(enum, migration_mode=True)
        write(route_path, routes, new_fields)

    print(f"PASS: migrated canonical confidence fields; immutable backup: {backup}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
