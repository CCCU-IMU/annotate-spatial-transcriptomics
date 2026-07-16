#!/usr/bin/env python3
"""Canonical confidence enums and explicit legacy migration mapping."""
from __future__ import annotations

CANONICAL = {"low", "moderate", "high"}
ATLAS = {"high", "moderate_only", "low_reject"}
LEGACY_MAP = {
    "low_confidence": "low", "medium_low": "low",
    "medium": "moderate", "medium_high": "moderate", "moderate_high": "moderate",
    "high_confidence": "high", "very_high": "high", "extreme": "high",
}
RANK = {"low": 0, "moderate": 1, "high": 2}


def normalized(value: object) -> str:
    return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")


def canonical(value: object, migration_mode: bool = False) -> str:
    result = normalized(value)
    if result in CANONICAL:
        return result
    if migration_mode and result in LEGACY_MAP:
        return LEGACY_MAP[result]
    raise ValueError(f"noncanonical confidence {result!r}; new projects allow low/moderate/high only")


def atlas_canonical(value: object, migration_mode: bool = False) -> str:
    result = normalized(value)
    if result in ATLAS:
        return result
    if migration_mode:
        mapping = {"moderate": "moderate_only", "medium": "moderate_only", "medium_high": "moderate_only", "moderate_high": "moderate_only", "low": "low_reject"}
        if result in mapping:
            return mapping[result]
    raise ValueError(f"noncanonical Atlas confidence {result!r}; use high/moderate_only/low_reject")
