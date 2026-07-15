#!/usr/bin/env python3
"""Parse generated job sources before scheduler submission without executing the workflow."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def command_for(path: Path, python: str, rscript: str) -> list[str]:
    suffix = path.suffix.lower()
    if suffix == ".py":
        return [python, "-m", "py_compile", str(path)]
    if suffix == ".r":
        return [rscript, "-e", f"parse(file={str(path)!r})"]
    if suffix in {".sh", ".aip"}:
        return ["bash", "-n", str(path)]
    raise ValueError(f"unsupported generated source type: {path}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    checks = []
    for path in args.files:
        if not path.is_file():
            checks.append({"file": str(path), "status": "FAIL", "error": "missing file"})
            continue
        try:
            if path.suffix.lower() == ".json":
                json.loads(path.read_text(encoding="utf-8"))
                checks.append({"file": str(path.resolve()), "status": "PASS", "check": "json_parse"})
            elif path.suffix.lower() in {".tsv", ".txt"}:
                with path.open(newline="", encoding="utf-8") as handle:
                    header = next(csv.reader(handle, delimiter="\t"), [])
                if not header or any(not value.strip() for value in header):
                    raise ValueError("tabular parameter file has an empty header")
                checks.append({"file": str(path.resolve()), "status": "PASS", "check": "tabular_header"})
            else:
                command = command_for(path, args.python, args.rscript)
                run = subprocess.run(command, text=True, capture_output=True, timeout=120)
                checks.append({
                    "file": str(path.resolve()), "status": "PASS" if run.returncode == 0 else "FAIL",
                    "command": command, "returncode": run.returncode,
                    "stderr_tail": run.stderr[-2000:], "stdout_tail": run.stdout[-1000:],
                })
        except Exception as exc:
            checks.append({"file": str(path), "status": "FAIL", "error": str(exc)})
    result = {"status": "PASS" if checks and all(row["status"] == "PASS" for row in checks) else "FAIL", "checks": checks}
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
