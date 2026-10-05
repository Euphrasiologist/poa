"""poa run: annotate one sample's mitochondrial and/or plastid assembly.

Takes the contig FASTA for each organelle (from oatk or any other
assembler), optionally with its assembly graph (GFA), and runs the
annotation pipeline on it:

  denovo_annotation 00-04 -> editing 01 -> trans_splicing 01 (mito)
  -> unitig_coords 00-01 + qc_basic_stats 06 (with a GFA) -> gff_export 01

The stages are the scripts under poa/pipeline/, run against a data root
built in OUTDIR/work (data/{mito,plastid}/<name>/<name>.{mito,pltd}.ctg.fasta
and .gfa). qc_basic_stats 01-05 are not run: they judge an assembly against
the rest of a dataset (genus-level size, a dataset-derived core-gene set,
oatk's own gene calls), so with no qc_summary.tsv every stage annotates
whatever it is given.

Output:
  OUTDIR/<name>.<organelle>.gff          contig-level GFF3
  OUTDIR/<name>.<organelle>.unitig.gff   unitig-level GFF3 (with a GFA)
  OUTDIR/tables/                         gene calls, edits, spliced genes, QC
  OUTDIR/logs/<stage>.log                each stage's stderr
  OUTDIR/work/                           the data root (intermediate files)
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

from . import PIPELINE_DIR

ORGANELLES = {"mito": "mito", "pltd": "plastid"}  # suffix -> data/ subdirectory
# a stage that reports failures this way failed, even with exit status 0
# (same rule as tests/smoke/run_smoke.sh)
FAILURE_RE = re.compile(r"failed=[1-9]|\[err\]| failed for ")
NAME_RE = re.compile(r"^[A-Za-z0-9_-]+$")
# per-stage result tables copied to OUTDIR/tables/
TABLES = [
    "denovo_annotation/results/gene_calls.tsv",
    "editing/results/edited_gene_calls.tsv",
    "editing/results/predicted_edits.tsv",
    "trans_splicing/results/reconstructed_genes.tsv",
    "trans_splicing/results/reconstructed_exons.tsv",
    "trans_splicing/results/reconstructed_junctions.tsv",
    "trans_splicing/results/reconstructed_edits.tsv",
    "unitig_coords/results/linearization_qc.tsv",
    "qc_basic_stats/results/low_depth_paths.tsv",
]


class StageFailed(Exception):
    pass


def add_arguments(p: argparse.ArgumentParser) -> None:
    p.add_argument("--mito", metavar="FASTA", help="mitochondrial contigs (.fasta, may be .gz)")
    p.add_argument("--mito-gfa", metavar="GFA", help="mitochondrial assembly graph, for unitig-level output")
    p.add_argument("--pltd", metavar="FASTA", help="plastid contigs (.fasta, may be .gz)")
    p.add_argument("--pltd-gfa", metavar="GFA", help="plastid assembly graph, for unitig-level output")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("-n", "--name", help="sample name used in output file names "
                   "(default: the first input's file name up to its first '.')")
    p.add_argument("-t", "--threads", type=int, default=4, help="default 4")
    p.add_argument("--force", action="store_true", help="replace an existing OUTDIR/work")


def default_name(path: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]", "_", Path(path).name.split(".")[0]) or "sample"


def check_input(src: str) -> None:
    path = Path(src)
    if not path.is_file():
        raise SystemExit(f"[err] no such file: {src}")
    if path.stat().st_size == 0:
        raise SystemExit(f"[err] empty input: {src}")


def stage_input(src: str, dst: Path) -> None:
    """Copy one input into the data root, decompressing .gz; reject empty files."""
    src_path = Path(src)
    opener = gzip.open if src_path.suffix == ".gz" else open
    with opener(src_path, "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    if dst.stat().st_size == 0:
        raise SystemExit(f"[err] empty input: {src}")


def run_stage(label: str, script: str, args: list[str], work: Path, logs: Path) -> None:
    log_path = logs / f"{label}.log"
    t0 = time.time()
    print(f"[poa] {label} ...", end="", file=sys.stderr, flush=True)
    env = dict(os.environ, PLANT_ORGANELLE_DATA_ROOT=str(work))
    with open(log_path, "w") as log:
        proc = subprocess.run([sys.executable, str(PIPELINE_DIR / script), *args],
                              cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT)
    text = log_path.read_text()
    if proc.returncode != 0 or FAILURE_RE.search(text):
        print(" FAILED", file=sys.stderr)
        tail = "\n".join(text.splitlines()[-20:])
        raise StageFailed(f"{label} failed (exit {proc.returncode}); last lines of {log_path}:\n{tail}")
    print(f" done ({time.time() - t0:.0f}s)", file=sys.stderr)


def run(args: argparse.Namespace) -> int:
    fastas = {o: getattr(args, o) for o in ORGANELLES if getattr(args, o)}
    gfas = {o: getattr(args, f"{o}_gfa") for o in ORGANELLES if getattr(args, f"{o}_gfa")}
    if not fastas:
        sys.exit("[err] give --mito and/or --pltd (contig FASTA)")
    for o in gfas:
        if o not in fastas:
            sys.exit(f"[err] --{o}-gfa needs --{o}: poa annotates contigs, and the graph alone "
                     f"doesn't say which path through it is the genome")
    name = args.name or default_name(next(iter(fastas.values())))
    if not NAME_RE.match(name):
        sys.exit(f"[err] --name may only contain letters, digits, '_' and '-': {name!r}")

    for path in (*fastas.values(), *gfas.values()):
        check_input(path)

    out = Path(args.outdir).resolve()
    work, logs = out / "work", out / "logs"
    if work.exists():
        if not args.force:
            sys.exit(f"[err] {work} exists - pass --force to replace it")
        shutil.rmtree(work)
    logs.mkdir(parents=True, exist_ok=True)
    try:
        for o in fastas:
            d = work / "data" / ORGANELLES[o] / name
            d.mkdir(parents=True)
            stage_input(fastas[o], d / f"{name}.{o}.ctg.fasta")
            if o in gfas:
                stage_input(gfas[o], d / f"{name}.{o}.gfa")
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise
    (work / "analysis" / "qc_basic_stats" / "results").mkdir(parents=True)

    orgs = list(fastas)
    organelle = orgs[0] if len(orgs) == 1 else "both"
    jobs = str(len(orgs))
    per_job = str(max(1, args.threads // len(orgs)))
    threads = str(max(1, args.threads))
    every_qc = ["--qc-status", "pass,flag,fail"]
    gene_calls = ["--gene-calls", str(work / "analysis/denovo_annotation/results/gene_calls.tsv")]
    stages = [
        ("nhmmscan", "denovo_annotation/src/00_run_nhmmscan.py",
         ["--organelle", organelle, "--jobs", jobs, "--cpu", per_job, *every_qc]),
        ("trnascan", "denovo_annotation/src/01_run_trnascan.py",
         ["--organelle", organelle, "--jobs", jobs, "--thread", per_job, *every_qc]),
        ("rrna", "denovo_annotation/src/02_run_barrnap.py",
         ["--organelle", organelle, "--jobs", jobs, "--threads", per_job, *every_qc]),
        ("combined_gff", "denovo_annotation/src/03_build_combined_gff.py", ["--organelle", organelle]),
        ("gene_calls", "denovo_annotation/src/04_build_gene_calls.py", ["--organelle", organelle]),
        ("editing", "editing/src/01_scan_editing.py",
         ["--organelle", organelle, "--threads", threads, *gene_calls]),
    ]
    if "mito" in orgs:
        stages.append(("trans_splicing", "trans_splicing/src/01_reconstruct.py",
                       ["--threads", threads, *gene_calls]))
    if gfas:
        stages += [
            ("unitig_map", "unitig_coords/src/00_build_unitig_map.py",
             ["--organelle", organelle, "--jobs", jobs, *every_qc]),
            ("linearization_qc", "unitig_coords/src/01_linearization_qc.py", ["--organelle", organelle]),
        ]
    if "mito" in gfas:
        stages.append(("low_depth_paths", "qc_basic_stats/src/06_low_depth_paths.py", []))
    stages.append(("gff_export", "gff_export/src/01_build_gff.py", ["--organelle", organelle, *gene_calls]))

    print(f"[poa] {name}: {', '.join(orgs)}{' (+ graph: ' + ', '.join(gfas) + ')' if gfas else ''}; "
          f"output in {out}", file=sys.stderr)
    try:
        for label, script, stage_args in stages:
            run_stage(label, script, stage_args, work, logs)
    except StageFailed as e:
        print(f"[err] {e}", file=sys.stderr)
        return 1

    results = work / "analysis" / "gff_export" / "results"
    written = []
    for o in orgs:
        for src, dst in ((results / f"{name}.{o}.gff", out / f"{name}.{o}.gff"),
                         (results / "unitig" / f"{name}.{o}.unitig.gff", out / f"{name}.{o}.unitig.gff")):
            if src.exists():
                shutil.copyfile(src, dst)
                written.append(dst)
        if not (out / f"{name}.{o}.gff").exists():
            print(f"[err] gff_export wrote no {o} GFF ({results})", file=sys.stderr)
            return 1
    (out / "tables").mkdir(exist_ok=True)
    for rel in TABLES:
        src = work / "analysis" / rel
        if src.exists():
            shutil.copyfile(src, out / "tables" / src.name)
    for path in written:
        n_genes = sum(1 for line in open(path) if "\tgene\t" in line)
        print(f"[poa] {path.name}: {n_genes} genes", file=sys.stderr)
    return 0
