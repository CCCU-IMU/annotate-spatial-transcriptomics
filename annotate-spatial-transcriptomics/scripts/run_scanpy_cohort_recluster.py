#!/usr/bin/env python3
"""Current entry point for Scanpy cohort reclustering.

The implementation filename is retained for compatibility with pre-v1.6 projects.
"""
from pathlib import Path
import runpy

runpy.run_path(str(Path(__file__).with_name("run_scanpy_pool_recluster.py")), run_name="__main__")
