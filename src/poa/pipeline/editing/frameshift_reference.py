"""Write editing/reference/frameshifts_dtol.tsv from a dataset's frameshifts.

    python -m poa.pipeline.editing.frameshift_reference FRAMESHIFTS_TSV

FRAMESHIFTS_TSV is editing/results/frameshifts.tsv from a run over the
whole dataset (or orfedit's own --frameshifts-out, gene labels with their
copy), with the default lineage's profiles. Writes, per organelle, gene
and profile column, the number of assemblies with a frameshift there;
scan.add_recurrence counts those within a few columns of each new
frameshift.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from poa.pipeline.editing.scan import FRAMESHIFT_REFERENCE


def build_frameshift_reference(frameshifts: Path, out: Path = FRAMESHIFT_REFERENCE) -> None:
    fs = pd.read_csv(frameshifts, sep="\t")
    fs["gene"] = fs["gene"].str.split("|").str[0]
    ref = (
        fs.groupby(["organelle", "gene", "profile_column"])["species"]
        .nunique()
        .rename("n_assemblies")
        .reset_index()
        .sort_values(["organelle", "gene", "profile_column"])
    )
    ref.to_csv(out, sep="\t", index=False)
    n = fs[["species", "organelle"]].drop_duplicates().shape[0]
    print(f"[info] frameshift_reference: {len(fs)} frameshifts in {n} assemblies -> {out}", file=sys.stderr)


if __name__ == "__main__":
    build_frameshift_reference(Path(sys.argv[1]))
