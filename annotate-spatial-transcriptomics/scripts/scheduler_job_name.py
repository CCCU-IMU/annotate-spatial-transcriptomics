#!/usr/bin/env python3
"""Create or validate compact, stage-readable scheduler job names."""

from __future__ import annotations

import argparse
import re


STAGES = {
    "00": "INPUT",
    "10": "SCT",
    "11": "SCTQC",
    "20": "RESGRID",
    "21": "RESEVID",
    "30": "BROAD",
    "40": "COHORT",
    "41": "COHORTQC",
    "50": "RCTD",
    "51": "ATLAS",
    "60": "RARE",
    "61": "CONTEXT",
    "70": "WRITEBACK",
    "75": "CONFIRM",
    "80": "FINALDEG",
    "81": "DOTPLOT",
    "82": "SPATIAL",
    "90": "REPORT",
    "99": "AUDIT",
}
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9-]*__P([0-9]{2})_([A-Z0-9]+)(?:_([A-Za-z0-9-]+))?__A([0-9]{2})$")


def clean(value: str, limit: int) -> str:
    value = re.sub(r"[^A-Za-z0-9-]+", "-", value.strip()).strip("-")
    if not value:
        raise ValueError("sample/scope becomes empty after scheduler-safe sanitization")
    return value[:limit].rstrip("-")


def validate(name: str, max_length: int) -> None:
    if len(name) > max_length:
        raise ValueError(f"job name has {len(name)} characters; maximum is {max_length}")
    match = NAME_RE.fullmatch(name)
    if not match:
        raise ValueError("job name does not match SAMPLE__Pnn_STAGE[_SCOPE]__Ann")
    code, label, _, attempt = match.groups()
    if code not in STAGES or label != STAGES[code]:
        raise ValueError(f"stage code/label mismatch; P{code} must be {STAGES.get(code, 'a registered stage')}")
    if int(attempt) < 1:
        raise ValueError("attempt must be A01 or greater")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample")
    parser.add_argument("--stage-code", choices=sorted(STAGES))
    parser.add_argument("--scope", default="")
    parser.add_argument("--attempt", type=int, default=1)
    parser.add_argument("--validate")
    parser.add_argument("--max-length", type=int, default=48)
    args = parser.parse_args()
    if args.validate:
        validate(args.validate, args.max_length)
        print(args.validate)
        return 0
    if not args.sample or not args.stage_code:
        parser.error("--sample and --stage-code are required when not using --validate")
    if not 1 <= args.attempt <= 99:
        parser.error("--attempt must be between 1 and 99")
    sample = clean(args.sample, 16)
    scope = clean(args.scope, 12) if args.scope else ""
    stage = f"P{args.stage_code}_{STAGES[args.stage_code]}"
    name = f"{sample}__{stage}{'_' + scope if scope else ''}__A{args.attempt:02d}"
    validate(name, args.max_length)
    print(name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
