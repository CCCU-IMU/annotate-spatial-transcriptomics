#!/usr/bin/env python3
"""Write reproducible Python, platform and R session information for release."""

from __future__ import annotations

import argparse
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    out = args.out or args.project_root / "provenance/release_sessionInfo.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    run = subprocess.run(
        [args.rscript, "-e", "sessionInfo()"], text=True, capture_output=True, timeout=120
    )
    if run.returncode != 0:
        raise SystemExit("R sessionInfo() failed: " + run.stderr[-2000:])
    lines = [
        f"generated_at_utc: {datetime.now(timezone.utc).isoformat()}",
        f"python: {sys.version}",
        f"python_executable: {sys.executable}",
        f"platform: {platform.platform()}",
        f"cwd: {Path.cwd()}",
        f"scheduler_job_id: {os.environ.get('CCP_JOBID', os.environ.get('SLURM_JOB_ID', ''))}",
        "",
        "R sessionInfo():",
        run.stdout,
    ]
    out.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
