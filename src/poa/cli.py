"""poa command line.

poa check      are the external tools and the oatkDB databases usable?
poa run        annotate one sample's mito and/or plastid assembly
"""

from __future__ import annotations

import shutil
import subprocess
import sys

import click

from . import __version__, run
from .pipeline.common.log import setup_logging

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
    ("gfatk", "--version", "QC graph checks (poa run with a GFA), linearize"),
    ("gfatools", "version", "dataset QC"),
    ("seqkit", "version", "dataset QC, linearize"),
    ("minimap2", "--version", "QC cross-organelle check (poa run with both organelles), linearize"),
    ("samtools", "--version", "linearize"),
    ("GraphAligner", "--version", "linearize"),
    ("mafft", "--version", "rebuilding reference profiles"),
]

INPUT = click.Path(exists=True, dir_okay=False)


def verbatim(doc: str) -> str:
    """Help text shown as written: click rewraps paragraphs unless each
    starts with a \\b line, which would mangle the stage diagram and the
    output table in run's docstring."""
    return "\n\n".join("\b\n" + para for para in doc.strip().split("\n\n"))


def tool_version(cmd: str, flag: str | None) -> str:
    if flag is None:
        return ""
    try:
        proc = subprocess.run([cmd, flag], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return ""
    out = (proc.stdout or proc.stderr).strip().splitlines()
    return out[0][:60] if out else ""


# An option each of these must have: the bundled profiles record editing
# sites in a column only a version with it can read (orfedit 0.2.0,
# transsplice 0.3.0).
NEEDS_OPTION = {"orfedit": "--novel-edit-penalty", "transsplice": "--novel-edit-penalty"}


def has_option(cmd: str, option: str) -> bool:
    try:
        proc = subprocess.run([cmd, "scan-batch", "--help"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return option in proc.stdout


def check_tools() -> int:
    ok = True
    print(f"poa {__version__}")
    print("tools for poa run:")
    for cmd, flag, why in REQUIRED:
        if shutil.which(cmd) is None:
            ok = False
            print(f"  MISSING  {cmd:13} ({why})")
        elif cmd in NEEDS_OPTION and not has_option(cmd, NEEDS_OPTION[cmd]):
            ok = False
            print(f"  TOO OLD  {cmd:13} {tool_version(cmd, flag)} (no {NEEDS_OPTION[cmd]})")
        else:
            print(f"  ok       {cmd:13} {tool_version(cmd, flag)}")
    print("tools for the dataset-scale stages (optional):")
    for cmd, flag, why in OPTIONAL:
        if shutil.which(cmd) is None:
            print(f"  -        {cmd:13} not found ({why})")
        else:
            print(f"  ok       {cmd:13} {tool_version(cmd, flag)}")
    print("oatkDB databases:")
    from .pipeline.common import oatkdb

    for organelle in ("mito", "pltd"):
        try:
            print(f"  ok       {organelle:13} {oatkdb.fam_path(organelle)}")
        except (OSError, subprocess.CalledProcessError) as e:
            ok = False
            print(f"  FAILED   {organelle:13} {e}")
    print("all good" if ok else "missing pieces above - see README.md, Install")
    return 0 if ok else 1


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="poa", message="%(prog)s %(version)s")
def main() -> None:
    """Plant organelle annotator."""
    setup_logging()


@main.command()
def check() -> None:
    """Check external tools and the oatkDB databases."""
    sys.exit(check_tools())


@main.command(name="run", help=verbatim(run.__doc__), short_help="Annotate one sample's mito and/or plastid assembly.")
@click.option("--mito", type=INPUT, metavar="FASTA", help="mitochondrial contigs (.fasta, may be .gz)")
@click.option("--mito-gfa", type=INPUT, metavar="GFA", help="mitochondrial assembly graph, for unitig-level output")
@click.option("--pltd", type=INPUT, metavar="FASTA", help="plastid contigs (.fasta, may be .gz)")
@click.option("--pltd-gfa", type=INPUT, metavar="GFA", help="plastid assembly graph, for unitig-level output")
@click.option("-o", "--outdir", required=True, type=click.Path(file_okay=False), help="output directory")
@click.option(
    "-n",
    "--name",
    help="sample name used in output file names (default: the first input's file name up to its first '.')",
)
@click.option(
    "--lineage",
    type=click.Choice(["auto", "angiosperm", "moss", "liverwort", "fern", "horsetail", "lycophyte"]),
    default="auto",
    show_default=True,
    help="reference lineage for mito editing profiles and spliced-gene templates: "
    "auto picks the one the sample's mito genes fit best. moss and liverwort have "
    "references of their own, from RefSeq genomes. fern (ferns; forced, also "
    "quillworts and hornworts) uses the angiosperm references with fern templates and "
    "profiles, and reads in-frame stop codons through as U-to-C RNA edits. horsetail "
    "(Equisetum) has exon templates and editing profiles from Equisetum RNA-seq. "
    "lycophyte (clubmosses, Lycopodiaceae) has them from RNA-edited clubmoss genes and "
    "reads stops through",
)
@click.option(
    "--genus",
    help="genus for the QC genome-size check (default: the first word of --name, e.g. Arabidopsis_thaliana)",
)
@click.option("-t", "--threads", type=click.IntRange(min=1), default=4, show_default=True, help="threads")
@click.option("--force", is_flag=True, help="replace an existing OUTDIR/work")
def run_command(**kwargs) -> None:
    sys.exit(run.run(**kwargs))


if __name__ == "__main__":
    main()
