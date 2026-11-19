"""Write editing/reference/frameshifts_dtol.tsv from a dataset's frameshifts.

    python -m poa.pipeline.editing.frameshift_reference EDITING_TSV [RECONSTRUCTED_TSV]

EDITING_TSV is editing/results/frameshifts.tsv from a run over the whole
dataset (or orfedit's own --frameshifts-out, gene labels with their copy),
RECONSTRUCTED_TSV trans_splicing/results/reconstructed_frameshifts.tsv,
both with the default lineage's profiles. Writes, per source (editing,
reconstruction), organelle, gene and profile column, the number of
assemblies with a frameshift there; scan.add_recurrence counts those
within a few columns of each new frameshift.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from poa.pipeline.editing.scan import FRAMESHIFT_REFERENCE


def count(path: Path, source: str) -> pd.DataFrame:
    fs = pd.read_csv(path, sep="\t")
    fs["gene"] = fs["gene"].str.split("|").str[0]
    out = fs.groupby(["organelle", "gene", "profile_column"])["species"].nunique().rename("n_assemblies").reset_index()
    n = fs[["species", "organelle"]].drop_duplicates().shape[0]
    print(f"[info] frameshift_reference: {source}: {len(fs)} frameshifts in {n} assemblies", file=sys.stderr)
    return out.assign(source=source)


def build_frameshift_reference(
    editing: Path, reconstructed: Path | None = None, out: Path = FRAMESHIFT_REFERENCE
) -> None:
    parts = [count(editing, "editing")]
    if reconstructed is not None:
        parts.append(count(reconstructed, "reconstruction"))
    ref = pd.concat(parts)[["source", "organelle", "gene", "profile_column", "n_assemblies"]]
    ref.sort_values(["source", "organelle", "gene", "profile_column"]).to_csv(out, sep="\t", index=False)
    print(f"[info] frameshift_reference: -> {out}", file=sys.stderr)


if __name__ == "__main__":
    build_frameshift_reference(Path(sys.argv[1]), Path(sys.argv[2]) if len(sys.argv) > 2 else None)
