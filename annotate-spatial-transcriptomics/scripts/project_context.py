#!/usr/bin/env python3
"""Resolve the project biological-context file across supported names."""

from __future__ import annotations

from pathlib import Path
from typing import Optional


CONTEXT_NAMES = ("biological_context.json", "context.json")


def resolve_context_path(project_root: Path, explicit: Optional[Path] = None) -> Path:
    if explicit is not None:
        return explicit
    config = project_root / "config"
    for name in CONTEXT_NAMES:
        candidate = config / name
        if candidate.is_file():
            return candidate
    return config / CONTEXT_NAMES[0]
