"""Command-line entry point:  desitter-qft {free-modes,two-point,tree,validate,all} [--outdir DIR]"""
from __future__ import annotations

import argparse
import os
import sys


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="desitter-qft", description=__doc__)
    ap.add_argument("command", choices=["free-modes", "two-point", "tree", "validate", "all"])
    ap.add_argument("--outdir", default="figures", help="figure output directory (results/ is created beside it)")
    a = ap.parse_args(argv)
    from . import figures as F

    if a.command in ("free-modes", "all"):
        F.make_free_mode_figures(a.outdir)
    if a.command in ("two-point", "all"):
        F.make_two_point_figures(a.outdir)
    if a.command in ("tree", "all"):
        F.make_tree_figures(a.outdir)
    if a.command in ("validate", "all"):
        from .validation import run_all, write_validation_record

        rep = run_all()
        rdir = os.path.join(os.path.dirname(os.path.normpath(a.outdir)), "results")
        os.makedirs(rdir, exist_ok=True)
        write_validation_record(rep, os.path.join(rdir, "validation.json"))
        for r in rep:
            print(f"[{'PASS' if r['passed'] else 'FAIL'}] {r['name']}: {r['value']:.3e} (tol {r['tolerance']:.1e})  -- {r['category']}")
        return 0 if all(r["passed"] for r in rep) else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
