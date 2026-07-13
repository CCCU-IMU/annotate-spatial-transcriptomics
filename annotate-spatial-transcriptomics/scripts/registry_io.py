#!/usr/bin/env python3
"""Locked and atomic TSV registry updates for concurrently submitted branches."""

from __future__ import annotations

import csv
import fcntl
import os
import tempfile
from pathlib import Path
from typing import Callable, List, Dict, Tuple


def locked_tsv_update(path: Path, mutate: Callable[[List[Dict[str, str]], List[str]], List[Dict[str, str]]]) -> None:
    lock_path = path.with_suffix(path.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("w") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        if not path.exists():
            raise FileNotFoundError(path)
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            rows = list(reader)
        updated = mutate(rows, fields)
        fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
        try:
            with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader(); writer.writerows(updated)
                handle.flush(); os.fsync(handle.fileno())
            Path(tmp_name).replace(path)
        finally:
            Path(tmp_name).unlink(missing_ok=True)
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
