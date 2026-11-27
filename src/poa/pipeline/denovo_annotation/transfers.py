"""DNA one organelle has taken from the other: plastid DNA in a mitogenome
(an "MTPT", common in seed plants) and, rarely, mitochondrial DNA in a
plastome. Two searches, each off unless asked for (poa run --cross-search,
--cross-align):

cross_search: each organelle's contigs are searched with the OTHER
organelle's oatkDB gene models too (nhmmscan, as for gene calling, E-value
below 1e-5). A hit is kept when no model of the contigs' own organelle
scores as well on the same stretch (overlapping it by at least
SAME_LOCUS_OVERLAP of the shorter, same strand) within MIN_MARGIN bits,
from the gene-calling scan's own tblout: plastid psbA in a mitogenome
stays, while the mitogenome's own rrn18 hit by the plastid rrn16 model, or
its rpl2 by the plastid rpl2 model, does not - their mito models fit them
better. So does a plastid trnS-GGA hit as well by the mito model (oatkDB's
mito trnS-GGA is itself plastid-like). Several
of the other organelle's models on one stretch (a tRNA hit by several
anticodon models) leave the best. What is kept is the other organelle's
gene, or a piece of it: `plastid_derived` in a mitogenome, `mito_derived`
in a plastome.

cross_align (both organelles): the plastid contigs are aligned to the mito
contigs (blastn -task dc-megablast, E-value below 1e-10; a stretch the
mitogenome has twice is reported twice, and one from the plastid's
inverted repeat twice), and every alignment of at least MIN_SHARED_BP at
MIN_IDENTITY or more is a stretch the two genomes share, reported in both
with where it lies in the other. Below MIN_IDENTITY, alignments are the
two organelles' rRNA genes, related from long before any transfer (on
five test samples, 60 of the 89 alignments under 0.8 lay over the
mitogenome's own rRNA genes, and none from 0.8 to 0.95 did). Of the
plastid genes cross_search found in those five mitogenomes (>= 150 bp),
the alignments hold 132 of 143; minimap2 (asm20) held 81, which is why
this is BLAST. An alignment does not say which way the DNA
moved, so these are `shared_with` the other organelle, not marked derived:
plastid-to-mitochondrion transfer is by far the commoner, but carrot's
plastome, for one, has mitochondrial DNA in its inverted repeat. A
cross_search hit inside a shared stretch says which way that part went.

Output: analysis/denovo_annotation/results/organelle_transfers.tsv, a row
per stretch per organelle it lies in (0-based half-open, like
gene_calls.tsv): species organelle contig start end strand method name
score e_value coverage identity aligned_bp own_best_score other_organelle
other_contig other_start other_end. gff_export writes each row as a
`biological_region` in that organelle's GFF.
        analysis/denovo_annotation/work/nhmmscan/<species>.<organelle>.cross.tblout
src/05_organelle_transfers.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable

import pandas as pd

from poa.pipeline.common import oatkdb
from poa.pipeline.common import species_discovery as sd
from poa.pipeline.common.errors import StageError
from poa.pipeline.common.log import setup_logging
from poa.pipeline.common.qc_gate import qc_eligible_species
from poa.pipeline.denovo_annotation import targets
from poa.pipeline.denovo_annotation.nhmmscan import DEFAULT_EVALUE, nhmmscan_exe
from poa.pipeline.denovo_annotation.nhmmscan_gff import ENV_FROM, ENV_TO, EVALUE, QUERY, SCORE, STRAND, TARGET, fam_model_lengths

log = logging.getLogger(__name__)

OTHER = {"mito": "pltd", "pltd": "mito"}
COLUMNS = [
    "species", "organelle", "contig", "start", "end", "strand", "method", "name", "score", "e_value", "coverage",
    "identity", "aligned_bp", "own_best_score", "other_organelle", "other_contig", "other_start", "other_end",
]
RESULTS = Path("analysis") / "denovo_annotation" / "results" / "organelle_transfers.tsv"
WORK = Path("analysis") / "denovo_annotation" / "work" / "nhmmscan"
# an own-organelle model hit this much of the shorter of the two is on the same stretch
SAME_LOCUS_OVERLAP = 0.5
# bits a cross hit must beat every own-organelle model on its stretch by
MIN_MARGIN = 10.0
BLAST_EVALUE = "1e-10"
MIN_SHARED_BP = 100
MIN_IDENTITY = 0.80
TEAM301_BLAST = Path("/software/team301/ncbi-blast-2.16.0+/bin")


def read_tblout(path: Path, max_evalue: float | None = None) -> pd.DataFrame:
    """nhmmscan --tblout -> model, contig, start (0-based), end, strand, e_value, score."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        a, b = int(f[ENV_FROM]), int(f[ENV_TO])
        e = float(f[EVALUE])
        if max_evalue is not None and not e < max_evalue:
            continue
        rows.append(
            dict(model=f[TARGET], contig=f[QUERY], start=min(a, b) - 1, end=max(a, b), strand=f[STRAND][0], e_value=e, score=float(f[SCORE]))
        )
    return pd.DataFrame(rows, columns=["model", "contig", "start", "end", "strand", "e_value", "score"])


def same_stretch(s1: int, e1: int, s2: int, e2: int) -> bool:
    return min(e1, e2) - max(s1, s2) >= SAME_LOCUS_OVERLAP * min(e1 - s1, e2 - s2)


def foreign_hits(cross: pd.DataFrame, own: pd.DataFrame) -> pd.DataFrame:
    """The cross hits no own-organelle model scores as well on, best first,
    one per stretch; with own_best_score, the best own score there (0 if none)."""
    kept: list[dict] = []
    own_by = {k: g for k, g in own.groupby(["contig", "strand"])} if not own.empty else {}
    for h in cross.sort_values("score", ascending=False).to_dict("records"):
        if any(k["contig"] == h["contig"] and k["strand"] == h["strand"] and same_stretch(k["start"], k["end"], h["start"], h["end"]) for k in kept):
            continue
        o = own_by.get((h["contig"], h["strand"]))
        here = [] if o is None else [r.score for r in o.itertuples() if same_stretch(r.start, r.end, h["start"], h["end"])]
        best = max(here, default=0.0)
        if h["score"] - best < MIN_MARGIN:
            continue
        kept.append({**h, "own_best_score": best})
    return pd.DataFrame(kept, columns=[*cross.columns, "own_best_score"])


def cross_search_one(species: str, organelle: str, ctg_fasta: str, cpu: int, work_dir: Path, force: bool) -> list[dict] | None:
    """Rows for one assembly, or None if it can't be searched."""
    own_tblout = work_dir / f"{species}.{organelle}.tblout"
    if not own_tblout.exists():
        log.warning(f"cross_search failed for {species} ({organelle}): no gene-calling scan ({own_tblout}) to compare with")
        return None
    other = OTHER[organelle]
    fam = oatkdb.fam_path(other)
    tblout = work_dir / f"{species}.{organelle}.cross.tblout"
    if force or not tblout.exists() or tblout.stat().st_size == 0:
        proc = subprocess.run(
            [nhmmscan_exe(), "--cpu", str(cpu), "--tblout", str(tblout), str(fam), ctg_fasta],
            capture_output=True, text=True, timeout=3600,
        )
        if proc.returncode != 0:
            log.warning(f"cross_search failed for {species} ({organelle}): {proc.stderr.strip()[:300]}")
            return None
    lengths = fam_model_lengths(fam)
    hits = foreign_hits(read_tblout(tblout, float(DEFAULT_EVALUE)), read_tblout(own_tblout))
    return [
        dict(
            species=species, organelle=organelle, contig=h["contig"], start=h["start"], end=h["end"], strand=h["strand"],
            method="cross_search", name=h["model"], score=h["score"], e_value=h["e_value"],
            coverage=round(min(1.0, (h["end"] - h["start"]) / lengths[h["model"]]), 3) if h["model"] in lengths else None,
            own_best_score=h["own_best_score"], other_organelle=other,
        )
        for h in hits.to_dict("records")
    ]


def blast_exe(name: str) -> str | None:
    found = shutil.which(name)
    return found or (str(TEAM301_BLAST / name) if (TEAM301_BLAST / name).exists() else None)


def cross_align_one(species: str, mito_fasta: str, pltd_fasta: str, threads: int) -> list[dict] | None:
    """Rows (both organelles) for one species' pair of assemblies, or None."""
    blastn, makeblastdb = blast_exe("blastn"), blast_exe("makeblastdb")
    if blastn is None or makeblastdb is None:
        log.warning(f"cross_align failed for {species}: blastn/makeblastdb not found")
        return None
    with tempfile.TemporaryDirectory(prefix="poa_cross_align.") as tmp:
        db = Path(tmp) / "mito"
        made = subprocess.run([makeblastdb, "-in", mito_fasta, "-dbtype", "nucl", "-out", str(db)], capture_output=True, text=True)
        if made.returncode != 0:
            log.warning(f"cross_align failed for {species}: makeblastdb: {made.stderr.strip()[-300:]}")
            return None
        proc = subprocess.run(
            [
                blastn, "-task", "dc-megablast", "-query", pltd_fasta, "-db", str(db), "-evalue", BLAST_EVALUE,
                "-num_threads", str(threads), "-max_target_seqs", "50", "-max_hsps", "1000",
                "-outfmt", "6 qseqid qstart qend sseqid sstart send nident length",
            ],
            capture_output=True, text=True, timeout=3600,
        )
    if proc.returncode != 0:
        log.warning(f"cross_align failed for {species}: blastn: {proc.stderr.strip()[-300:]}")
        return None
    rows = []
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 8:
            continue
        q, qs, qe, t, ts, te, ident, length = f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]), int(f[6]), int(f[7])
        if length < MIN_SHARED_BP or ident / length < MIN_IDENTITY:
            continue
        # the plastid is the query, always on its + strand; strand: the mito stretch's, relative to it
        strand = "+" if ts <= te else "-"
        ts, te = min(ts, te), max(ts, te)
        common = dict(species=species, method="cross_align", strand=strand, identity=round(ident / length, 3), aligned_bp=length)
        rows.append(dict(common, organelle="mito", contig=t, start=ts - 1, end=te, name="shared with the plastid", other_organelle="pltd", other_contig=q, other_start=qs - 1, other_end=qe))
        rows.append(dict(common, organelle="pltd", contig=q, start=qs - 1, end=qe, name="shared with the mitogenome", other_organelle="mito", other_contig=t, other_start=ts - 1, other_end=te))
    return rows


def find_transfers(
    root: Path,
    organelles: Iterable[str] = targets.ORGANELLES,
    *,
    cross_search: bool = True,
    cross_align: bool = True,
    species: set[str] | None = None,
    cpu: int = 2,
    jobs: int = 1,
    qc_status: Iterable[str] = ("pass",),
    force: bool = False,
) -> Path:
    """Run the searches asked for over every assembly under data root
    `root` and write RESULTS (rows of earlier runs for other species kept
    when `species` is given). Returns its path."""
    root = Path(root)
    organelles = targets.check_organelles(organelles)
    qc_eligible = qc_eligible_species(root, qc_status)
    work_dir = root / WORK
    work_dir.mkdir(parents=True, exist_ok=True)
    found = {o: targets.assemblies(root, o, species, qc_eligible) for o in organelles}
    rows: list[dict] = []

    def collect(out: list[dict] | None) -> str:
        if out is None:
            return "failed"
        rows.extend(out)
        return "scanned"

    if cross_search:
        for organelle in organelles:
            try:
                oatkdb.fam_path(OTHER[organelle])
            except (OSError, subprocess.CalledProcessError) as e:
                raise StageError(f"{OTHER[organelle]}: can't prepare the oatkDB database: {e}") from e
            counts = targets.scan_all(
                found[organelle], lambda r: collect(cross_search_one(r.species, organelle, r.ctg_fasta, cpu, work_dir, force)), jobs
            )
            targets.report("cross_search", organelle, counts, work_dir)
    if cross_align:
        if set(organelles) != set(targets.ORGANELLES):
            raise StageError("cross_align needs both organelles")
        if blast_exe("blastn") is None or blast_exe("makeblastdb") is None:
            raise StageError("cross_align needs BLAST+ (blastn, makeblastdb) on PATH")
        plastids = {r.species: r.ctg_fasta for r in found["pltd"]}
        pairs = [r for r in found["mito"] if r.species in plastids]
        counts = targets.scan_all(pairs, lambda r: collect(cross_align_one(r.species, r.ctg_fasta, plastids[r.species], cpu)), jobs)
        targets.report("cross_align", "mito+pltd", counts, work_dir)

    out = root / RESULTS
    out.parent.mkdir(parents=True, exist_ok=True)
    new = pd.DataFrame(rows, columns=COLUMNS)
    if species is not None and out.exists():
        old = pd.read_csv(out, sep="\t")
        new = pd.concat([old[~old.species.isin(species)], new], ignore_index=True)
    new = new.sort_values(["species", "organelle", "contig", "start", "method"], kind="stable")
    for col in ("start", "end", "other_start", "other_end", "aligned_bp"):
        new[col] = new[col].astype("Int64")
    new.to_csv(out, sep="\t", index=False)
    by = ", ".join(f"{o} {m} {n}" for (o, m), n in new.groupby(["organelle", "method"]).size().items())
    log.info(f"organelle_transfers: {len(new)} rows -> {out}" + (f" ({by})" if by else ""))
    return out


def main(argv: list[str] | None = None) -> int:
    setup_logging()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--no-cross-search", action="store_true", help="skip the search with the other organelle's models")
    ap.add_argument("--no-cross-align", action="store_true", help="skip aligning the plastid to the mitogenome")
    ap.add_argument("--cpu", type=int, default=2, help="threads per process (default 2)")
    ap.add_argument("--jobs", type=int, default=1, help="processes run concurrently (default 1)")
    ap.add_argument("--qc-status", default="pass", help="comma-separated subset of pass,flag,fail (default: pass)")
    ap.add_argument("--force", action="store_true", help="search again even where a cross tblout exists")
    args = ap.parse_args(argv)
    root = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
    orgs = list(targets.ORGANELLES) if args.organelle == "both" else [args.organelle]
    try:
        find_transfers(
            root,
            orgs,
            cross_search=not args.no_cross_search,
            cross_align=not args.no_cross_align and len(orgs) == 2,
            species=sd.load_species_list(Path(args.species_list)) if args.species_list else None,
            cpu=args.cpu,
            jobs=args.jobs,
            qc_status=args.qc_status.split(","),
            force=args.force,
        )
    except StageError as e:
        log.error(e)
        return 1
    return 0
