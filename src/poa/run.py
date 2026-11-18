"""poa run: annotate one sample's mitochondrial and/or plastid assembly.

Takes the contig FASTA for each organelle (from oatk or any other
assembler), optionally with its assembly graph (GFA), and runs the
annotation pipeline on it:

  denovo_annotation 00-04 -> editing 01 -> trans_splicing 01 (mito)
  -> unitig_coords 00-01 + qc_basic_stats 06 (with a GFA) -> gff_export 01
  -> qc (pass/flag/fail per organelle)

Each stage is a library call, run in this process against a data root
built in OUTDIR/work (data/{mito,plastid}/<name>/<name>.{mito,pltd}.ctg.fasta
and .gfa), with its output (and that of the tools it runs) going to
OUTDIR/logs/<stage>.log.
The dataset QC gate (qc_basic_stats 01-05) judges an assembly partly
against the rest of a dataset, so it isn't run; every stage annotates
whatever it is given. The qc stage applies the same rules to this sample
alone (qc_basic_stats/sample_qc.py), with the dataset-derived parts - the
core-gene set and genus genome sizes - shipped as reference data, and
writes a pass/flag/fail line per organelle to OUTDIR/<name>.qc.tsv. Graph
checks need the GFA, and the cross-organelle check both organelles; a
check that can't run is listed under not_assessed.

Output:
  OUTDIR/<name>.<organelle>.gff          contig-level GFF3
  OUTDIR/<name>.<organelle>.unitig.gff   unitig-level GFF3 (with a GFA)
  OUTDIR/<name>.qc.tsv                   pass/flag/fail and reasons, per organelle
  OUTDIR/tables/                         gene calls, edits, spliced genes, QC
  OUTDIR/logs/<stage>.log                each stage's stderr
  OUTDIR/work/                           the data root (intermediate files)
"""

from __future__ import annotations

import csv
import gzip
import logging
import os
import re
import shutil
import sys
import time
import traceback
from contextlib import contextmanager
from pathlib import Path
from typing import Callable, Iterator

from .pipeline.common.errors import StageError
from .pipeline.denovo_annotation.combined_gff import build_combined_gff
from .pipeline.denovo_annotation.gene_calls import build_gene_calls
from .pipeline.denovo_annotation.nhmmscan import run_nhmmscan
from .pipeline.denovo_annotation.rrna import run_rrna
from .pipeline.denovo_annotation.trnascan import run_trnascan
from .pipeline.editing.scan import scan as scan_editing
from .pipeline.gff_export.build_gff import build_gff
from .pipeline.qc_basic_stats.low_depth_paths import low_depth_paths
from .pipeline.qc_basic_stats.sample_qc import sample_qc
from .pipeline.trans_splicing.reconstruct import reconstruct
from .pipeline.unitig_coords.linearization_qc import linearization_qc
from .pipeline.unitig_coords.unitig_map import build_unitig_maps

log = logging.getLogger(__name__)

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
    "editing/results/lineage.tsv",
    "trans_splicing/results/reconstructed_genes.tsv",
    "trans_splicing/results/reconstructed_exons.tsv",
    "trans_splicing/results/reconstructed_junctions.tsv",
    "trans_splicing/results/reconstructed_edits.tsv",
    "unitig_coords/results/linearization_qc.tsv",
    "qc_basic_stats/results/low_depth_paths.tsv",
]


class StageFailed(Exception):
    pass


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


@contextmanager
def fds_to(log) -> Iterator[None]:
    """Send this process's stdout and stderr (fds 1 and 2, so the tools a
    stage runs too) to `log` for the duration."""
    sys.stdout.flush()
    sys.stderr.flush()
    saved = os.dup(1), os.dup(2)
    try:
        os.dup2(log.fileno(), 1)
        os.dup2(log.fileno(), 2)
        yield
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        for fd in saved:
            os.close(fd)


def run_stage(label: str, stage: Callable[[], object], logs: Path) -> None:
    """Call `stage()`, with its output going to logs/<label>.log."""
    log_path = logs / f"{label}.log"
    t0 = time.time()
    print(f"[poa] {label} ...", end="", file=sys.stderr, flush=True)
    failed = False
    with open(log_path, "w") as fh:
        try:
            with fds_to(fh):
                stage()
        except StageError as e:
            failed = True
            print(f"[err] {e}", file=fh)
        except Exception:
            failed = True
            traceback.print_exc(file=fh)
    text = log_path.read_text()
    if failed or FAILURE_RE.search(text):
        print(" FAILED", file=sys.stderr)
        tail = "\n".join(text.splitlines()[-20:])
        raise StageFailed(f"{label} failed; last lines of {log_path}:\n{tail}")
    print(f" done ({time.time() - t0:.0f}s)", file=sys.stderr)


def run(
    *,
    outdir: str | os.PathLike,
    mito: str | None = None,
    mito_gfa: str | None = None,
    pltd: str | None = None,
    pltd_gfa: str | None = None,
    name: str | None = None,
    lineage: str = "auto",
    genus: str | None = None,
    threads: int = 4,
    force: bool = False,
) -> int:
    """Annotate one sample (see the module docstring); returns an exit status.
    Bad input exits via SystemExit with an [err] message."""
    inputs = {"mito": (mito, mito_gfa), "pltd": (pltd, pltd_gfa)}
    fastas = {o: fa for o, (fa, _) in inputs.items() if fa}
    gfas = {o: gfa for o, (_, gfa) in inputs.items() if gfa}
    if not fastas:
        sys.exit("[err] give --mito and/or --pltd (contig FASTA)")
    for o in gfas:
        if o not in fastas:
            sys.exit(
                f"[err] --{o}-gfa needs --{o}: poa annotates contigs, and the graph alone "
                f"doesn't say which path through it is the genome"
            )
    name = name or default_name(next(iter(fastas.values())))
    if not NAME_RE.match(name):
        sys.exit(f"[err] --name may only contain letters, digits, '_' and '-': {name!r}")

    for path in (*fastas.values(), *gfas.values()):
        check_input(path)

    out = Path(outdir).resolve()
    work, logs = out / "work", out / "logs"
    if work.exists():
        if not force:
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
    jobs = len(orgs)
    per_job = max(1, threads // len(orgs))
    threads = max(1, threads)
    every_qc = ("pass", "flag", "fail")
    gene_calls = work / "analysis/denovo_annotation/results/gene_calls.tsv"
    scan = dict(jobs=jobs, qc_status=every_qc)
    stages = [
        ("nhmmscan", lambda: run_nhmmscan(work, orgs, cpu=per_job, **scan)),
        ("trnascan", lambda: run_trnascan(work, orgs, threads=per_job, **scan)),
        ("rrna", lambda: run_rrna(work, orgs, threads=per_job, **scan)),
        ("combined_gff", lambda: build_combined_gff(work, orgs)),
        ("gene_calls", lambda: build_gene_calls(work, orgs)),
        (
            "editing",
            lambda: scan_editing(work, orgs, gene_calls=gene_calls, threads=threads, lineage=lineage),
        ),
    ]
    if "mito" in orgs:
        stages.append(
            (
                "trans_splicing",
                lambda: reconstruct(work, gene_calls=gene_calls, threads=threads, lineage=lineage),
            )
        )
    if gfas:
        stages += [
            ("unitig_map", lambda: build_unitig_maps(work, orgs, **scan)),
            ("linearization_qc", lambda: linearization_qc(work, orgs)),
        ]
    if "mito" in gfas:
        stages.append(("low_depth_paths", lambda: low_depth_paths(work)))
    stages.append(("gff_export", lambda: build_gff(work, orgs, gene_calls=gene_calls)))
    stages.append(("qc", lambda: sample_qc(work, orgs, gene_calls=gene_calls, genus=genus)))

    print(
        f"[poa] {name}: {', '.join(orgs)}{' (+ graph: ' + ', '.join(gfas) + ')' if gfas else ''}; " f"output in {out}",
        file=sys.stderr,
    )
    try:
        for label, stage in stages:
            run_stage(label, stage, logs)
    except StageFailed as e:
        log.error(e)
        return 1

    results = work / "analysis" / "gff_export" / "results"
    written = []
    for o in orgs:
        for src, dst in (
            (results / f"{name}.{o}.gff", out / f"{name}.{o}.gff"),
            (results / "unitig" / f"{name}.{o}.unitig.gff", out / f"{name}.{o}.unitig.gff"),
        ):
            if src.exists():
                shutil.copyfile(src, dst)
                written.append(dst)
        if not (out / f"{name}.{o}.gff").exists():
            log.error(f"gff_export wrote no {o} GFF ({results})")
            return 1
    (out / "tables").mkdir(exist_ok=True)
    for rel in TABLES:
        src = work / "analysis" / rel
        if src.exists():
            shutil.copyfile(src, out / "tables" / src.name)
    for path in written:
        n_genes = sum(1 for line in open(path) if "\tgene\t" in line)
        print(f"[poa] {path.name}: {n_genes} genes", file=sys.stderr)
    qc_src = work / "analysis" / "qc_basic_stats" / "results" / "sample_qc.tsv"
    if qc_src.exists():
        qc_out = out / f"{name}.qc.tsv"
        shutil.copyfile(qc_src, qc_out)
        with open(qc_out) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                why = f" ({r['status_reasons']})" if r["status_reasons"] else ""
                skipped = re.sub(r"\([^)]*\)", "", r["not_assessed"]).replace(";", ", ")
                skipped = f"; not assessed: {skipped}" if skipped else ""
                print(f"[poa] QC {r['organelle']}: {r['status']}{why}{skipped} -> {qc_out.name}", file=sys.stderr)
    return 0
