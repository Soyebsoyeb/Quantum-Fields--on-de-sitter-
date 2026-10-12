#!/usr/bin/env python
"""Reproduce fig01-fig04 and the matching JSON record in results/.

Usage:  python scripts/reproduce_free_modes.py [--outdir figures]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from desitter_qft.cli import main

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"))
    a = ap.parse_args()
    raise SystemExit(main(["free-modes", "--outdir", a.outdir]))
