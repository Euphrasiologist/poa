"""Build per-exon-slot and whole-gene reference profiles for `transsplice`.

Source: `../reference/raw/*.fasta` + `junctions.tsv` - real GenBank RefSeq
mitochondrial genomes (Arabidopsis, Beta vulgaris, maize, rice, Vitis),
fetched and independently re-verified (every exon re-translated from its
own raw coordinates and checked against the curated RefSeq protein - see
`../README.md`), plus, for the cis-spliced genes (reference_genes.py),
six more angiosperm RefSeq mitogenomes fetched and verified the same way
by fetch_references.py. Genes: reference_genes.GENES.

Two kinds of profile per gene, both written under `../reference/profiles/
<gene>/`:
- `exon{N}.pssm` + `exon{N}.threshold`: one per exon slot (N = 1..K),
  built from the loci that have a canonical (modal) exon count - a locus
  missing an intron entirely (e.g. Beta vulgaris `rps3`) or flagged
  UNRESOLVED (rice `rps3` - a likely assembly artifact, confirmed not a
  transcription error on our end) is simply absent from a slot's training
  set rather than distorting it. `.threshold` is FLOOR_FRACTION of the
  minimum score of each slot's own training sequences run back through
  `orfedit scan-batch` against the profile just built from them - the
  same scoring code path `transsplice` uses at runtime.
- `exon{N}.frame`: `lead<TAB>trail` - how many bases at the exon's 5'
  end finish a codon begun in the previous exon, and how many at its 3'
  end start one the next exon finishes (the modal value across training
  loci). Exon profiles are built from the exon's *full codons only*,
  translated in their true frame: translating every exon from its own
  first base put every exon after a phase-1/2 junction out of frame
  (nad1 exons 2/5, nad2 exons 3-5, nad5 exons 2/4, rps3 exon 2), and
  joining exons on codon boundaries then dropped the split codon's bases,
  frame-shifting the rest of the gene - transsplice adds lead/trail back.
- `whole_gene.pssm`: one full-length profile per gene, from every
  included locus's curated RefSeq protein (regardless of its own exon
  count - Beta vulgaris' intron-less `rps3` protein is still a real,
  correct full-length `rps3`, just not useful for the per-exon slots).

`--lineage moss` (or any other non-angiosperm lineage in
reference_genes.LINEAGE_MITOGENOMES) reads reference/lineages/<lineage>/raw/
and writes reference/lineages/<lineage>/profiles/ instead - the same
build, on that lineage's RefSeq loci (fetch_references.py --lineage).
A gene the lineage has lost (moss rps10) gets no profile, and a gene with
no intron there (moss rps3, rpl2) gets a one-slot template.

src/00_build_exon_profiles.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from poa import PIPELINE_DIR
from poa.pipeline.common import lineage as lin
from poa.pipeline.common.pssm import align_and_build, composition, mark_edited, run_mafft, translate, read_fasta, write_profile
from poa.pipeline.common.refseq_mito import n_edited, unedited
from poa.pipeline.common.tools import cargo_tool
from poa.pipeline.trans_splicing.reference_genes import (
    GENES,
    LINEAGE_EDITED_PROTEINS,
    LINEAGE_PLASTID,
    mito_spliced,
    LINEAGE_MITOGENOMES,
    LINEAGE_PLASTOMES,
    LINEAGE_UNEDITED_MITOGENOMES,
    PLASTID_SPLICED,
    PLASTOMES,
    reference_dir,
)
from poa.pipeline.common.log import setup_logging
from poa.pipeline.trans_splicing import column_mask

log = logging.getLogger(__name__)

# Slot floor = this fraction of the lowest training self-score. With exons
# translated in frame, self-scores are high (nad5 exon 2: 1013), and a
# floor at the bare minimum rejected even Arabidopsis's own nad5 exons 1-2
# (a training genome) whenever the candidate window didn't span the whole
# exon. It was 0.5 with the old slot profiles. The SLOT_PSEUDOCOUNT ones
# are sharper: a reference exon's self-score rose by about half, a
# divergent species' exon by a fifth to two-fifths, so at 0.5 ferns and
# lycophytes (on the angiosperm templates) lost most of their exons - 117
# complete mito genes in 29 DToL assemblies fell to 27 - and angiosperms
# lost some too. Lower is better as far down as tried: complete mito genes
# in 955 assemblies 6,617 at 0.5, 7,025 at 0.25, 7,277 at 0.15, 7,552 at
# 0.075 (7,620 at 0.025), genes with a frameshift call 1,159, 1,039, 818,
# 779 (789); plastid and moss likewise, with no complete gene lost; and
# held out on the 24 RefSeq genomes placement is the same throughout (426
# of 477 exons right). The cost is a few more short exons placed that
# resemble no reference exon (plastid clpP exon 1: 111 of 1,108, against
# 64 of 1,061 at 0.25). Negative self-scores (tiny exons) are kept as-is.
FLOOR_FRACTION = 0.075

# Exon-slot profiles score against the whole protein's composition, with
# this many pseudo-observations per column shared out in proportion to it
# (pssm.build_profile). Before, the background was the slot alignment's own
# composition with a flat pseudocount, which inverts the profile of a short
# exon: a residue the exon lacks scored as a match in every column (nad5
# exon 3, 7 codons: +1.6 for 16 of 20 residues, self-score floor -2.4).
# Held out on the 24 RefSeq genomes: 2 placed 74% of exons within 3 bp
# (5: 72%, 10: 71%, the old profiles 65%); every exon of 31-60 bp against
# 2 of 15. None: the old profiles.
SLOT_PSEUDOCOUNT = 2.0

# Alignment columns gappier than this are left out of a lineage's
# profiles (pssm.GAP_COLUMN_THRESHOLD, 0.9, unless given here). Fern: at
# 0.9 a stretch only one reference has stays in - a C-terminal tail past
# where the other ferns' genes stop - and with stop codons read through
# transsplice aligned on into it: 35 of 89 spliced genes in 11 DToL ferns
# ran past a stop codon their RNA shows is never edited.
LINEAGE_GAP_THRESHOLD = {"fern": 0.5, "horsetail": 0.5, "lycophyte": 0.5}
# Fewest loci a template is built from (2 unless given here). Horsetail:
# the lineage is one genus, Equisetum, whose species' proteins differ by a
# few percent, and some genes have RNA-seq deep enough in one species
# only; a profile of that one edited protein still places the others'.
# Lycophyte: one clubmoss (Phylloglossum) has a transcriptome deep enough
# to edit every gene, and its proteins alone, held out by genus on three
# RefSeq clubmosses, placed 88% of exons within 3 bp, against 78% for the
# four unedited RefSeq proteins together.
LINEAGE_MIN_LOCI = {"horsetail": 1, "lycophyte": 1}
# Pseudo-observations per column for a lineage's whole-gene and exon-slot
# profiles, where SLOT_PSEUDOCOUNT and the whole-gene default are too
# many: against one to four proteins (horsetail, lycophyte) 2
# pseudo-observations outweigh the data, a mismatch costs half a bit, and
# an RNA edit is hardly worth its penalty - held out, the other Equisetum
# species' edits were found at 50%.
LINEAGE_PSEUDOCOUNT = {"horsetail": 0.25, "lycophyte": 0.25}


# Loci that are not GenBank records: genes of DToL assemblies with their
# RNA-seq editing applied (reference/lineages/fern/raw_rna/README.md), in
# the same format as raw/, which fetch_references.py rewrites.
RNA_RAW_NAME = "raw_rna"


def raw_dirs(ref_dir: Path) -> list[Path]:
    return [d for d in (ref_dir / "raw", ref_dir / RNA_RAW_NAME) if d.is_dir()]


def load_all_raw_data(*raw_dirs: Path) -> dict[tuple[str, str, str], dict]:
    """{(accession, gene, locus_tag): {exon_num_int: dna_seq, "refseq": protein_seq}}"""
    data: dict[tuple[str, str, str], dict] = defaultdict(dict)
    for fasta_path in (path for raw_dir in raw_dirs for path in sorted(raw_dir.glob("*.fasta"))):
        for header, seq in read_fasta(fasta_path).items():
            parts = header.split("|")
            accession, gene, locus_tag, record_type = parts[0], parts[1], parts[2], parts[3]
            key = (accession, gene, locus_tag)
            if record_type.startswith("exon"):
                n = int(record_type[len("exon") :].split("of")[0])
                data[key][n] = seq
            elif record_type == "translation_REFSEQ_ANNOTATED":
                data[key]["refseq"] = seq
    return data


def load_junctions(*raw_dirs: Path) -> dict[tuple[str, str, str], dict]:
    """{(accession, gene, locus_tag): {"n_exons": int, "excluded": bool, "reason": str}}"""
    out: dict[tuple[str, str, str], dict] = {}
    for raw_dir in raw_dirs:
        with open(raw_dir / "junctions.tsv") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                key = (row["accession"], row["gene"], row["locus_tag"])
                if key in out:
                    continue
                excluded = "UNRESOLVED" in row["notes"]
                # an RNA-derived locus whose reads are not the assembly's plant (rna_refs.py): its
                # protein is used, its differences from the genome are not taken for editing sites
                out[key] = {"n_exons": int(row["n_exons"]), "excluded": excluded, "sites": NO_SITES_NOTE not in row["notes"]}
    return out


# Lineages whose templates take a gene's full intron set, from loci that
# do not all have it (union_exons). Angiosperm cox2 has two intron
# positions: the Nymphaea, Nelumbo and Sedum genes have both, beet's has
# the first only, Arabidopsis' and the cucurbits' the second only, the
# legumes' neither. Slots from the loci of the commonest exon count mixed
# the two two-exon forms, and cox2 was placed in none of 18 RefSeq genomes.
LINEAGE_UNION_INTRONS = ("angiosperm",)
NO_SITES_NOTE = "sites not used"
# junctions of two loci are the same intron when they fall within this
# many columns of the proteins' alignment
JUNCTION_COLUMNS = 4


def union_exons(
    gene: str, included: list, raw: dict, junctions: dict, work_dir: Path, min_loci: int
) -> tuple[dict, dict] | None:
    """A gene's loci recut to its full intron set, or None when every
    locus already has it.

    Each locus' introns are placed on an alignment of the loci's proteins;
    positions shared by at least min_loci loci make the set. A locus
    without one of them is cut at the homologous codon, in the phase most
    loci have there, so every locus gives every exon slot; an intron only
    one locus has is ignored (its two exons stay joined). Returns ({locus:
    {exon number: DNA, "refseq": protein}}, {locus: {"n_exons": K,
    "excluded": False}}) for the loci with exons and a protein."""
    loci = [k for k in included if "refseq" in raw[k] and k in junctions]
    exons = {k: [raw[k][n] for n in sorted(n for n in raw[k] if isinstance(n, int))] for k in loci}
    loci = [k for k in loci if exons[k]]
    if len(loci) < 2:
        return None
    fasta, aligned = work_dir / f"{gene}.union.protein.fasta", work_dir / f"{gene}.union.aligned.fasta"
    fasta.write_text("".join(f">{i}\n{raw[k]['refseq']}\n" for i, k in enumerate(loci)))
    if not run_mafft(fasta, aligned):
        return None
    rows = read_fasta(aligned)
    # per locus: residue index -> alignment column, and back
    columns = {}
    for i, k in enumerate(loci):
        columns[k] = [c for c, ch in enumerate(rows[str(i)]) if ch != "-"]
    cuts = []  # (column, locus, nt position in the joined CDS)
    for k in loci:
        cum = 0
        for seq in exons[k][:-1]:
            cum += len(seq)
            residue = min(cum // 3, len(columns[k]) - 1)
            cuts.append((columns[k][residue], k, cum))
    cuts.sort(key=lambda c: c[0])
    clusters: list[list] = []
    for cut in cuts:
        if clusters and cut[0] - clusters[-1][-1][0] <= JUNCTION_COLUMNS:
            clusters[-1].append(cut)
        else:
            clusters.append([cut])
    shared = [c for c in clusters if len({k for _, k, _ in c}) >= min_loci]
    counts = Counter(len(exons[k]) for k in loci)
    modal = counts.most_common(1)[0][0]
    has = {k: {i for i, c in enumerate(shared) if any(k2 == k for _, k2, _ in c)} for k in loci}
    if all(len(has[k]) == len(shared) for k in loci) and modal == len(shared) + 1:
        return None
    if len(shared) + 1 < modal:
        # the alignment did not bring the loci's introns together: leave the gene as it is
        log.warning(f"{gene}: {modal - 1} intron(s) in most loci but {len(shared)} shared position(s) - not recut")
        return None
    new_raw, new_junctions = {}, {}
    for k in loci:
        cds = "".join(exons[k])
        points = []
        for c in shared:
            own = [nt for _, k2, nt in c if k2 == k]
            if own:
                points.append(own[0])
                continue
            column = sorted(col for col, _, _ in c)[len(c) // 2]
            phase = Counter(nt % 3 for _, _, nt in c).most_common(1)[0][0]
            residue = sum(1 for col in columns[k] if col < column)
            points.append(min(3 * residue + phase, len(cds)))
        points = sorted(set(points))
        if len(points) != len(shared) or points[0] <= 0 or points[-1] >= len(cds):
            continue  # the locus lacks the stretch an intron sits in
        bounds = [0, *points, len(cds)]
        new_raw[k] = {n + 1: cds[a:b] for n, (a, b) in enumerate(zip(bounds, bounds[1:]))}
        new_raw[k]["refseq"] = raw[k]["refseq"]
        new_junctions[k] = {"n_exons": len(shared) + 1, "excluded": False}
    log.info(
        f"{gene}: full intron set of {len(shared)} from loci with {dict(counts)} exons "
        f"({len(new_raw)} loci recut, {len(loci) - len(new_raw)} left out)"
    )
    return new_raw, new_junctions


def build_and_write_profile(
    protein_seqs: dict[str, str],
    out_path: Path,
    label: str,
    work_dir: Path,
    background: dict[str, float] | None = None,
    pseudocount_mass: float | None = None,
    trim_ends: tuple[bool, bool] = (False, False),
    gap_threshold: float | None = None,
    min_sequences: int = 2,
    groups: dict[str, str] | None = None,
) -> bool:
    # groups: the proteins mark their edited residues (lower case), counted by these (pssm.align_and_build)
    if len(protein_seqs) < min_sequences:
        log.warning(f"{label}: only {len(protein_seqs)} sequences, need >= {min_sequences} - skipping")
        return False
    columns = align_and_build(
        protein_seqs,
        work_dir,
        label,
        **({} if gap_threshold is None else {"gap_threshold": gap_threshold}),
        background=background,
        pseudocount_mass=pseudocount_mass,
        trim_ends=trim_ends,
        groups=groups,
        sites=groups is not None,
    )
    if columns is None:
        log.warning(f"{label}: no usable alignment - skipping")
        return False
    write_profile(columns, out_path)
    return True


def compute_min_self_score(
    gene: str,
    slot: int,
    dna_seqs: dict[str, str],
    profile_path: Path,
    work_dir: Path,
    orfedit: str,
    stop_edit_penalty: float | None = None,
) -> float | None:
    """Score every training sequence for this slot back against the profile
    just built from it, via a real `orfedit scan-batch` call (same code
    path `transsplice` uses at runtime) - the floor is the minimum of these,
    not an invented cutoff."""
    manifest_rows = []
    for name, seq in dna_seqs.items():
        fasta_path = work_dir / f"selfscore.{gene}.exon{slot}.{name}.fasta"
        fasta_path.write_text(f">self\n{seq}\n")
        manifest_rows.append(
            {
                "species": name,
                "organelle": "mito",
                "gene": f"{gene}_exon{slot}",
                "fasta": str(fasta_path),
                "profile": str(profile_path),
                "contig": "self",
                "start": 0,
                "end": len(seq),
                "strand": "+",
            }
        )
    manifest_path = work_dir / f"selfscore.{gene}.exon{slot}.manifest.tsv"
    cols = ["species", "organelle", "gene", "fasta", "profile", "contig", "start", "end", "strand"]
    with open(manifest_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in manifest_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    out_path = work_dir / f"selfscore.{gene}.exon{slot}.out.tsv"
    edits_path = work_dir / f"selfscore.{gene}.exon{slot}.edits.tsv"
    proc = subprocess.run(
        [
            orfedit,
            "scan-batch",
            "--manifest",
            str(manifest_path),
            "--gene-calls-out",
            str(out_path),
            "--edits-out",
            str(edits_path),
            "--threads",
            "1",
            "--flank",
            "0",
            # a lineage that edits U to C (fern): score its exons as they
            # will be scored when found, in-frame stops read through
            *(["--stop-edit-penalty", str(stop_edit_penalty)] if stop_edit_penalty is not None else []),
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )
    if proc.returncode != 0:
        log.warning(f"{gene} exon{slot}: self-scoring orfedit call failed: {proc.stderr.strip()[:300]}")
        return None
    scores = []
    with open(out_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            scores.append(float(row["score"]))
    return min(scores) if scores else None


def build_exon_profiles(
    lineage: str = lin.DEFAULT_LINEAGE,
    genes: list[str] | None = None,
    work_root: Path | None = None,
    organelle: str = "mito",
    exclude: tuple[str, ...] = (),
    out_dir: Path | None = None,
    slot_pseudocount: float | None = SLOT_PSEUDOCOUNT,
    slots_only: bool = False,
    use_mask: bool = True,
) -> Path:
    """(Re)build `lineage`'s exon-slot and whole-gene templates in the
    package's reference directory, from its raw/ exon sets. genes: only
    these (default reference_genes.GENES) - rebuilding a gene overwrites
    its whole_gene.pssm, including one bootstrap.py refined. Alignments and
    self-scoring files go under work_root/analysis/trans_splicing/work
    (default: the current directory). slots_only: leave an existing
    whole_gene.pssm alone (the exon slots still score against the RefSeq
    proteins' composition). Returns the profiles directory."""
    if lineage not in LINEAGE_MITOGENOMES:
        raise ValueError(f"lineage must be one of {sorted(LINEAGE_MITOGENOMES)}, not {lineage!r}")
    # Lineages that edit U to C (fern): the genome's own translation of an
    # exon has stop codons where the protein has Gln or Arg, so slot
    # profiles are cut from the annotated (edited) protein instead.
    own_plastid = organelle == "pltd" and lineage in LINEAGE_PLASTID
    stop_edit = lin.stop_edit_penalty(lineage) if organelle == "mito" or own_plastid else None
    genes = list(genes or (PLASTID_SPLICED if organelle == "pltd" else mito_spliced(lineage)))
    # references whose proteins are the RNA-edited ones: slot profiles are
    # cut from the protein
    edited_proteins = (organelle == "mito" or own_plastid) and (stop_edit is not None or lineage in LINEAGE_EDITED_PROTEINS)
    orfedit = cargo_tool("orfedit")
    ref_dir = reference_dir(PIPELINE_DIR, organelle, lineage)
    sources = raw_dirs(ref_dir)
    # references whose proteins are not the edited ones (fern)
    if organelle == "mito":
        exclude = tuple(exclude) + tuple(LINEAGE_UNEDITED_MITOGENOMES.get(lineage, ()))
    # The shared plastid set, once it has RNA-edited loci (raw_rna/, genes
    # of DToL angiosperms): built from those alone, slots cut from their
    # proteins. The RefSeq plastomes' proteins are the genomes' own
    # translations, and against templates of them 2 of 17 edits called in
    # six angiosperms' spliced plastid genes were in their RNA.
    if organelle == "pltd" and not own_plastid and (ref_dir / RNA_RAW_NAME).is_dir():
        exclude = tuple(exclude) + tuple(PLASTOMES)
        edited_proteins = True
    profiles_dir = Path(out_dir) if out_dir else ref_dir / "profiles"
    work_dir = Path(work_root or os.getcwd()) / "analysis" / "trans_splicing" / "work"
    if organelle == "pltd":
        work_dir = work_dir / ("plastid" if not own_plastid else f"plastid_{lineage}")
    elif lineage != lin.DEFAULT_LINEAGE:
        work_dir = work_dir / lineage
    work_dir.mkdir(parents=True, exist_ok=True)
    profiles_dir.mkdir(parents=True, exist_ok=True)

    raw = load_all_raw_data(*sources)
    junctions = load_junctions(*sources)
    organisms = {**PLASTOMES, **LINEAGE_PLASTOMES.get(lineage, {})} if organelle == "pltd" else LINEAGE_MITOGENOMES[lineage]
    mask = column_mask.load_mask(lineage) if use_mask and organelle == "mito" else {}
    gap_threshold = LINEAGE_GAP_THRESHOLD.get(lineage) if organelle == "mito" or own_plastid else None
    min_loci = LINEAGE_MIN_LOCI.get(lineage, 2) if organelle == "mito" or own_plastid else 2
    sharp = LINEAGE_PSEUDOCOUNT.get(lineage) if organelle == "mito" or own_plastid else None
    if sharp is not None and slot_pseudocount is not None:
        slot_pseudocount = sharp

    for gene in genes:
        loci = [key for key in raw if key[1] == gene and key[0] not in exclude]
        included = [key for key in loci if not junctions.get(key, {}).get("excluded", False)]
        excluded = [key for key in loci if key not in included]
        if excluded:
            log.info(f"{gene}: excluding {[k[2] for k in excluded]} (flagged in junctions.tsv)")
        # angiosperm: a record annotated with edited proteins can still give
        # a gene or two unedited (beet ccmFc, rice nad4): those loci out too
        if lineage == lin.DEFAULT_LINEAGE and organelle == "mito":
            counts = {
                key: n_edited("".join(raw[key][n] for n in sorted(k for k in raw[key] if isinstance(k, int))), raw[key]["refseq"])
                for key in included
                if "refseq" in raw[key]
            }
            # (not a locus whose protein is from RNA-seq: few edits there is
            # what its reads show)
            plain = {key for key in unedited(counts) if not key[0].startswith("RNA_")}
            if plain:
                log.info(f"{gene}: excluding {sorted(k[0] for k in plain)} (annotated unedited: {sorted(counts.values())})")
                included = [key for key in included if key not in plain]

        gene_dir = profiles_dir / gene

        # Whole-gene profile: every included locus's curated RefSeq protein,
        # regardless of its own exon count.
        whole_seqs = {}
        for key in included:
            acc, _, locus = key
            if "refseq" in raw[key]:
                whole_seqs[f"{acc}_{locus}"] = raw[key]["refseq"]
        # ... with the residues editing made marked, where the proteins are the
        # edited ones (lin.SITE_LINEAGES): the editing sites of the whole-gene
        # profile and, cut from the same proteins, of the exon slots'
        marked, groups = whole_seqs, None
        if lineage in lin.SITE_LINEAGES and edited_proteins:
            marked, groups = {}, {}
            for key in included:
                acc, _, locus = key
                if "refseq" not in raw[key]:
                    continue
                name, protein = f"{acc}_{locus}", raw[key]["refseq"]
                dna = "".join(raw[key][n] for n in sorted(k for k in raw[key] if isinstance(k, int)))
                marked[name] = mark_edited(dna, protein) if junctions.get(key, {}).get("sites", True) else protein.upper()
                groups[name] = acc.split("_")[1] if acc.startswith("RNA_") else organisms.get(acc, acc).split()[0]
        if not (slots_only and (gene_dir / "whole_gene.pssm").exists()):
            built = build_and_write_profile(
                marked,
                gene_dir / "whole_gene.pssm",
                f"{gene}.whole_gene",
                work_dir,
                gap_threshold=gap_threshold,
                min_sequences=min_loci,
                groups=groups,
                **({} if sharp is None else dict(background=composition(whole_seqs.values()), pseudocount_mass=sharp)),
            )
            # columns the lineage's own species don't have (column_mask.py)
            if built and mask.get(gene):
                n = column_mask.apply_mask(gene_dir / "whole_gene.pssm", mask[gene])
                log.info(f"{gene}: dropped {n} masked whole-gene column(s)")
        # exon slots score against the whole protein's composition, not
        # their own few residues' (see pssm.build_profile)
        slot_background = composition(whole_seqs.values()) if slot_pseudocount is not None else None

        # A gene whose loci differ in which introns they have: slots for the
        # full set, every locus recut to it (union_exons)
        if lineage in LINEAGE_UNION_INTRONS and organelle == "mito":
            recut = union_exons(gene, included, raw, junctions, work_dir, min_loci)
            if recut is not None:
                raw = {**raw, **recut[0]}
                junctions = {**junctions, **recut[1]}
                included = [key for key in included if key in recut[0]]

        # Per-exon-slot profiles: only loci with the canonical (modal) exon
        # count contribute a given slot - naturally excludes Beta vulgaris'
        # intron-less rps3 from exon2 (it simply has no exon2 record).
        n_exons_seen = Counter(junctions[key]["n_exons"] for key in included if key in junctions)
        if not n_exons_seen:
            log.warning(f"{gene}: no junction data for included loci - skipping exon slots")
            continue
        k = n_exons_seen.most_common(1)[0][0]
        log.info(f"{gene}: {len(included)} loci, canonical exon count K={k} " f"(seen: {dict(n_exons_seen)})")

        # lead/trail per (locus, slot) from the locus's own cumulative exon
        # lengths - the junction phase, which the translation must respect
        frames: dict[tuple, dict[int, tuple[int, int]]] = {}
        for key in included:
            if junctions.get(key, {}).get("n_exons") != k or not all(n in raw[key] for n in range(1, k + 1)):
                continue
            cum, per = 0, {}
            for n in range(1, k + 1):
                length = len(raw[key][n])
                lead = (3 - cum % 3) % 3
                trail = (cum + length) % 3
                per[n] = (lead, trail)
                cum += length
            frames[key] = per

        # the annotated protein's residues for each slot's whole codons
        slot_proteins: dict[int, dict[str, str]] = {n: {} for n in range(1, k + 1)}
        if edited_proteins:
            for key, per in frames.items():
                acc, _, locus = key
                # (with its editing sites marked, where it has them)
                protein, cum = marked.get(f"{acc}_{locus}", raw[key].get("refseq", "")), 0
                for n in range(1, k + 1):
                    lead, trail = per[n]
                    first, n_codons = (cum + lead) // 3, (len(raw[key][n]) - lead - trail) // 3
                    slot_proteins[n][f"{acc}_{locus}"] = protein[first : first + n_codons]
                    cum += len(raw[key][n])

        for slot in range(1, k + 1):
            dna_seqs = {}
            for key in included:
                acc, _, locus = key
                # Require this locus's OWN exon count to match the
                # canonical K, not just "has a record at this slot index" -
                # an intron-less locus (e.g. Beta vulgaris rps3, n_exons=1)
                # stores its single whole-CDS record as "exon1of1", which
                # would otherwise corrupt exon1's profile with a ~551aa
                # full-length sequence instead of the real ~24aa leader
                # exon (caught empirically: min_self_score was -483 before
                # this check).
                if junctions.get(key, {}).get("n_exons") != k:
                    continue
                if slot in raw[key] and key in frames:
                    lead, trail = frames[key][slot]
                    seq = raw[key][slot]
                    dna_seqs[f"{acc}_{locus}"] = seq[lead : len(seq) - trail]
            if len(dna_seqs) < min_loci:
                log.warning(f"{gene} exon{slot}: only {len(dna_seqs)} sequences - skipping")
                continue
            if not edited_proteins:
                protein_seqs = {name: translate(seq).rstrip("*") for name, seq in dna_seqs.items()}
            else:
                protein_seqs = slot_proteins[slot]
            profile_path = gene_dir / f"exon{slot}.pssm"
            # intron-side ends only (pssm.build_profile's trim_ends), and
            # only with the new profiles: None rebuilds the old ones as
            # they were
            trim = (slot > 1, slot < k) if slot_pseudocount is not None else (False, False)
            if not build_and_write_profile(
                protein_seqs,
                profile_path,
                f"{gene}.exon{slot}",
                work_dir,
                slot_background,
                slot_pseudocount,
                trim,
                gap_threshold,
                min_loci,
                groups,
            ):
                continue
            min_score = compute_min_self_score(gene, slot, dna_seqs, profile_path, work_dir, orfedit, stop_edit)
            if min_score is None:
                log.warning(f"{gene} exon{slot}: could not compute min_self_score - skipping")
                continue
            floor = min_score * FLOOR_FRACTION if min_score > 0 else min_score
            (gene_dir / f"exon{slot}.threshold").write_text(f"{floor:.4f}\n")
            lead, trail = Counter(frames[key][slot] for key in frames).most_common(1)[0][0]
            (gene_dir / f"exon{slot}.frame").write_text(f"{lead}\t{trail}\n")
            log.info(
                f"{gene} exon{slot}: {len(dna_seqs)} training seqs, "
                f"min_self_score={min_score:.3f} -> {profile_path.name}"
            )

    return profiles_dir


def main(argv: list[str] | None = None) -> int:
    setup_logging()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--genes",
        nargs="+",
        default=None,
        help="build only these (default: all the organelle's). Rebuilding a gene overwrites its whole_gene.pssm, "
        "including one 02_bootstrap_refine.py refined from this dataset",
    )
    ap.add_argument(
        "--lineage",
        default=lin.DEFAULT_LINEAGE,
        choices=sorted(LINEAGE_MITOGENOMES),
        help="which lineage's references to build from (default: %(default)s)",
    )
    ap.add_argument("--organelle", default="mito", choices=["mito", "pltd"], help="pltd: the plastid templates")
    ap.add_argument(
        "--no-mask",
        action="store_true",
        help="keep every whole-gene profile column, ignoring the lineage's whole_gene_mask.tsv "
        "(column_mask.py; do this before deriving a new mask)",
    )
    ap.add_argument("--exclude", nargs="+", default=[], metavar="ACC", help="leave these accessions out (held-out tests)")
    ap.add_argument("--out-dir", type=Path, default=None, help="write the templates here, not into the package")
    ap.add_argument(
        "--slot-pseudocount",
        type=lambda x: None if x.lower() == "none" else float(x),
        default=SLOT_PSEUDOCOUNT,
        help="exon-slot profiles score against the whole protein's composition, with this many "
        "pseudo-observations per column (default: %(default)s). none: the slot alignment's own "
        "composition and flat pseudocounts, as before",
    )
    ap.add_argument(
        "--slots-only",
        action="store_true",
        help="rebuild the exon-slot profiles and floors only, keeping each gene's whole_gene.pssm "
        "(one 02_bootstrap_refine.py refined, say)",
    )
    args = ap.parse_args(argv)
    build_exon_profiles(
        args.lineage,
        args.genes,
        Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd()),
        args.organelle,
        tuple(args.exclude),
        args.out_dir,
        args.slot_pseudocount,
        args.slots_only,
        not args.no_mask,
    )
    return 0
