"""Whole-gene profile columns to drop because the species a lineage's
templates are used on do not have them.

A lineage's whole-gene profile is built from its reference proteins as
annotated. Where there are few references and some are unvalidated
assemblies, the profile picks up stretches most real genes lack - the
fern nad7 profile began with 22 columns that only the Woodwardia and
Blechnopsis records have, 30 residues upstream of where Haplopteris and
Osmunda start the protein - and every gene aligned to it then leaves
those columns empty.

`write_mask` reads a transsplice genes table made with the *unmasked*
templates (build them with `00_build_exon_profiles.py --no-mask`) and
lists, per gene, the profile columns that at least MIN_UNALIGNED of the
complete reconstructions leave unaligned. `exon_profiles` drops the
listed columns from whole_gene.pssm whenever it builds that lineage. The
mask can only remove columns: it does not move the profile towards the
reconstructions, so it cannot learn a misplaced exon the way refining a
profile from reconstructions can.

    python -m poa.pipeline.trans_splicing.column_mask GENES_TSV --lineage fern [--species REGEX]
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from poa.pipeline.common import lineage as lin

PIPELINE_DIR = Path(__file__).resolve().parents[1]
MASK_NAME = "whole_gene_mask.tsv"
MIN_UNALIGNED = 0.8
MIN_GENES = 5


def mask_path(lineage: str) -> Path:
    return lin.trans_splicing_reference(PIPELINE_DIR, lineage) / MASK_NAME


def load_mask(lineage: str) -> dict[str, set[int]]:
    """{gene: 0-based columns of the unmasked whole-gene profile to drop}"""
    path = mask_path(lineage)
    mask: dict[str, set[int]] = {}
    if path.exists():
        with open(path) as fh:
            for row in csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"):
                mask.setdefault(row["gene"], set()).add(int(row["column"]))
    return mask


def apply_mask(profile_path: Path, columns: set[int]) -> int:
    """Drop `columns` (0-based) from a .pssm file in place; how many went."""
    lines = Path(profile_path).read_text().splitlines(keepends=True)
    kept = [lines[0]] + [line for i, line in enumerate(lines[1:]) if i not in columns]
    Path(profile_path).write_text("".join(kept))
    return len(lines) - len(kept)


def unsupported_columns(proteins: list[str], min_unaligned: float = MIN_UNALIGNED) -> list[int]:
    """Columns ('-' in transsplice's protein column) unaligned in at least
    min_unaligned of `proteins`, which all span the same profile."""
    if not proteins:
        return []
    n = len(proteins)
    return [i for i in range(len(proteins[0])) if sum(p[i] == "-" for p in proteins) >= min_unaligned * n]


def write_mask(genes_tsv: Path, lineage: str, species: str | None = None) -> Path:
    by_gene: dict[str, list[str]] = {}
    with open(genes_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["complete"].lower() != "true" or (species and not re.search(species, row["species"])):
                continue
            by_gene.setdefault(row["gene"], []).append(row["protein"])
    out = mask_path(lineage)
    with open(out, "w") as fh:
        fh.write(
            f"# whole-gene profile columns (0-based, of the unmasked build) that at least {MIN_UNALIGNED:.0%} of\n"
            f"# complete reconstructions leave unaligned; exon_profiles drops them. See column_mask.py.\n"
            "gene\tcolumn\tn_genes\tn_unaligned\n"
        )
        for gene in sorted(by_gene):
            proteins = [p for p in by_gene[gene] if len(p) == len(by_gene[gene][0])]
            if len(proteins) < MIN_GENES:
                continue
            for i in unsupported_columns(proteins):
                fh.write(f"{gene}\t{i}\t{len(proteins)}\t{sum(p[i] == '-' for p in proteins)}\n")
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genes_tsv", help="transsplice genes table from a run on the unmasked templates")
    ap.add_argument("--lineage", required=True)
    ap.add_argument("--species", default=None, help="only species matching this regular expression")
    args = ap.parse_args(argv)
    out = write_mask(Path(args.genes_tsv), args.lineage, args.species)
    mask = load_mask(args.lineage)
    print(f"{out}: " + ", ".join(f"{g} {len(c)}" for g, c in sorted(mask.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
