#!/usr/bin/env python3
"""Current entry point for Scanpy cohort reclustering."""
from pathlib import Path
import runpy

runpy.run_path(str(Path(__file__).with_name("run_scanpy_cohort_recluster_impl.py")), run_name="__main__")
