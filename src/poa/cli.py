"""poa command line.

    poa check      are the external tools and the oatkDB databases usable?
    poa run        annotate one sample's mito and/or plastid assembly
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys

from . import PIPELINE_DIR, __version__, run

# (command, version flag, what needs it); the version flag is only for display.
# REQUIRED is what `poa run` needs; the rest only serve stages run by hand
# on a whole dataset (QC 01-05, linearize, rebuilding reference profiles),
# and mafft/GraphAligner have no bioconda build on some platforms.
REQUIRED = [
    ("nhmmscan", "-h", "gene calling"),
    ("nhmmer", "-h", "rRNA calling"),
    ("hmmpress", "-h", "preparing the bundled oatkDB"),
    ("tRNAscan-SE", None, "tRNA calling"),
    ("orfedit", "--version", "RNA-editing-aware ORF correction"),
    ("transsplice", "--version", "spliced-gene reconstruction"),
]
OPTIONAL = [
    ("gfatk", "--version", "dataset QC, linearize"),
    ("gfatools", "version", "dataset QC"),
    ("seqkit", "version", "dataset QC, linearize"),
    ("minimap2", "--version", "dataset QC, linearize"),
    ("samtools", "--version", "linearize"),
    ("GraphAligner", "--version", "linearize"),
    ("mafft", "--version", "rebuilding reference profiles"),
]


def tool_version(cmd: str, flag: str | None) -> str:
    if flag is None:
        return ""
    try:
        proc = subprocess.run([cmd, flag], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return ""
    out = (proc.stdout or proc.stderr).strip().splitlines()
    return out[0][:60] if out else ""


def check(_args: argparse.Namespace) -> int:
    ok = True
    print(f"poa {__version__}")
    print("tools for poa run:")
    for cmd, flag, why in REQUIRED:
        if shutil.which(cmd) is None:
            ok = False
            print(f"  MISSING  {cmd:13} ({why})")
        else:
            print(f"  ok       {cmd:13} {tool_version(cmd, flag)}")
    print("tools for the dataset-scale stages (optional):")
    for cmd, flag, why in OPTIONAL:
        if shutil.which(cmd) is None:
            print(f"  -        {cmd:13} not found ({why})")
        else:
            print(f"  ok       {cmd:13} {tool_version(cmd, flag)}")
    print("oatkDB databases:")
    sys.path.insert(0, str(PIPELINE_DIR / "common"))
    import oatkdb
    for organelle in ("mito", "pltd"):
        try:
            print(f"  ok       {organelle:13} {oatkdb.fam_path(organelle)}")
        except (OSError, subprocess.CalledProcessError) as e:
            ok = False
            print(f"  FAILED   {organelle:13} {e}")
    print("all good" if ok else "missing pieces above - see README.md, Install")
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="poa", description="Plant organelle annotator")
    ap.add_argument("--version", action="version", version=f"poa {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)
    p = sub.add_parser("check", help="check external tools and the oatkDB databases")
    p.set_defaults(func=check)
    p = sub.add_parser("run", help="annotate one sample's mito and/or plastid assembly",
                       description=run.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    run.add_arguments(p)
    p.set_defaults(func=run.run)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
