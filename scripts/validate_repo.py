#!/usr/bin/env python3
from __future__ import annotations

import ast
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "annotate-spatial-transcriptomics"
TEXT_SUFFIXES = {".md", ".py", ".R", ".r", ".json", ".yaml", ".yml", ".tsv", ".txt", ".sh"}
BANNED = (
    "/" + "share" + "/" + "org" + "/",
    "bgi" + "_" + "zhangsch",
    "bgi" + "_" + "baiyy",
    "bgi" + "_" + "xinq",
    "D055" + "22A3",
)


def fail(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)


def main() -> int:
    problems: list[str] = []
    required = [
        ROOT / "README.md",
        ROOT / "install.sh",
        ROOT / "LICENSE",
        SKILL / "SKILL.md",
        SKILL / "references" / "multi-route-controller.md",
        SKILL / "references" / "report-contract.md",
    ]
    for path in required:
        if not path.is_file():
            problems.append(f"missing required file: {path.relative_to(ROOT)}")

    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or ".git" in path.parts or ".release_extract" in path.parts:
            continue
        if "__pycache__" in path.parts or path.suffix == ".pyc":
            problems.append(f"generated cache file present: {path.relative_to(ROOT)}")
            continue
        if path.suffix in TEXT_SUFFIXES or path.name in {"LICENSE", "VERSION"}:
            try:
                text = path.read_text(encoding="utf-8")
            except UnicodeDecodeError as exc:
                problems.append(f"not UTF-8 text: {path.relative_to(ROOT)}: {exc}")
                continue
            for token in BANNED:
                if token in text:
                    problems.append(f"private/sample token {token!r} in {path.relative_to(ROOT)}")
            if path.suffix == ".py":
                try:
                    ast.parse(text, filename=str(path))
                except SyntaxError as exc:
                    problems.append(f"Python syntax: {path.relative_to(ROOT)}: {exc}")
            elif path.suffix == ".json":
                try:
                    json.loads(text)
                except json.JSONDecodeError as exc:
                    problems.append(f"JSON syntax: {path.relative_to(ROOT)}: {exc}")

    if problems:
        for problem in problems:
            fail(problem)
        return 1
    print("PASS: repository structure, portability scan, Python syntax and JSON validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
