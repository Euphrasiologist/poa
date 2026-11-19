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
from poa.pipeline.common.pssm import align_and_build, composition, translate, read_fasta, write_profile
from poa.pipeline.common.tools import cargo_tool
from poa.pipeline.trans_splicing.reference_genes import GENES, LINEAGE_MITOGENOMES, PLASTID_SPLICED, reference_dir
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


def load_all_raw_data(raw_dir: Path) -> dict[tuple[str, str, str], dict]:
    """{(accession, gene, locus_tag): {exon_num_int: dna_seq, "refseq": protein_seq}}"""
    data: dict[tuple[str, str, str], dict] = defaultdict(dict)
    for fasta_path in raw_dir.glob("*.fasta"):
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


def load_junctions(raw_dir: Path) -> dict[tuple[str, str, str], dict]:
    """{(accession, gene, locus_tag): {"n_exons": int, "excluded": bool, "reason": str}}"""
    out: dict[tuple[str, str, str], dict] = {}
    with open(raw_dir / "junctions.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["accession"], row["gene"], row["locus_tag"])
            if key in out:
                continue
            excluded = "UNRESOLVED" in row["notes"]
            out[key] = {"n_exons": int(row["n_exons"]), "excluded": excluded}
    return out


def build_and_write_profile(
    protein_seqs: dict[str, str],
    out_path: Path,
    label: str,
    work_dir: Path,
    background: dict[str, float] | None = None,
    pseudocount_mass: float | None = None,
    trim_ends: tuple[bool, bool] = (False, False),
) -> bool:
    if len(protein_seqs) < 2:
        log.warning(f"{label}: only {len(protein_seqs)} sequences, need >= 2 - skipping")
        return False
    columns = align_and_build(
        protein_seqs, work_dir, label, background=background, pseudocount_mass=pseudocount_mass, trim_ends=trim_ends
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
    stop_edit = lin.stop_edit_penalty(lineage) if organelle == "mito" else None
    genes = list(genes or (PLASTID_SPLICED if organelle == "pltd" else GENES))
    orfedit = cargo_tool("orfedit")
    ref_dir = reference_dir(PIPELINE_DIR, organelle, lineage)
    raw_dir = ref_dir / "raw"
    profiles_dir = Path(out_dir) if out_dir else ref_dir / "profiles"
    work_dir = Path(work_root or os.getcwd()) / "analysis" / "trans_splicing" / "work"
    if organelle == "pltd":
        work_dir = work_dir / "plastid"
    elif lineage != lin.DEFAULT_LINEAGE:
        work_dir = work_dir / lineage
    work_dir.mkdir(parents=True, exist_ok=True)
    profiles_dir.mkdir(parents=True, exist_ok=True)

    raw = load_all_raw_data(raw_dir)
    junctions = load_junctions(raw_dir)
    mask = column_mask.load_mask(lineage) if use_mask and organelle == "mito" else {}

    for gene in genes:
        loci = [key for key in raw if key[1] == gene and key[0] not in exclude]
        included = [key for key in loci if not junctions.get(key, {}).get("excluded", False)]
        excluded = [key for key in loci if key not in included]
        if excluded:
            log.info(f"{gene}: excluding {[k[2] for k in excluded]} (flagged in junctions.tsv)")

        gene_dir = profiles_dir / gene

        # Whole-gene profile: every included locus's curated RefSeq protein,
        # regardless of its own exon count.
        whole_seqs = {}
        for key in included:
            acc, _, locus = key
            if "refseq" in raw[key]:
                whole_seqs[f"{acc}_{locus}"] = raw[key]["refseq"]
        if not (slots_only and (gene_dir / "whole_gene.pssm").exists()):
            built = build_and_write_profile(whole_seqs, gene_dir / "whole_gene.pssm", f"{gene}.whole_gene", work_dir)
            # columns the lineage's own species don't have (column_mask.py)
            if built and mask.get(gene):
                n = column_mask.apply_mask(gene_dir / "whole_gene.pssm", mask[gene])
                log.info(f"{gene}: dropped {n} masked whole-gene column(s)")
        # exon slots score against the whole protein's composition, not
        # their own few residues' (see pssm.build_profile)
        slot_background = composition(whole_seqs.values()) if slot_pseudocount is not None else None

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
        if stop_edit is not None:
            for key, per in frames.items():
                acc, _, locus = key
                protein, cum = raw[key].get("refseq", ""), 0
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
            if len(dna_seqs) < 2:
                log.warning(f"{gene} exon{slot}: only {len(dna_seqs)} sequences - skipping")
                continue
            if stop_edit is None:
                protein_seqs = {name: translate(seq).rstrip("*") for name, seq in dna_seqs.items()}
            else:
                protein_seqs = slot_proteins[slot]
            profile_path = gene_dir / f"exon{slot}.pssm"
            # intron-side ends only (pssm.build_profile's trim_ends), and
            # only with the new profiles: None rebuilds the old ones as
            # they were
            trim = (slot > 1, slot < k) if slot_pseudocount is not None else (False, False)
            if not build_and_write_profile(
                protein_seqs, profile_path, f"{gene}.exon{slot}", work_dir, slot_background, slot_pseudocount, trim
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
