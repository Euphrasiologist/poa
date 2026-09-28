#!/usr/bin/env python3
"""Flag species whose mito and plastid marker-gene trees disagree on where
it sits phylogenetically - a contamination/mislabeling signal that's
independent of the assembly-graph checks in qc_basic_stats (a species can
have one perfectly clean-looking subgraph and still be the wrong species'
mito paired with the wrong species' plastid).

Not merged into qc_basic_stats/results/qc_summary.tsv's pass/flag/fail gate:
this needs both trees to already exist, which themselves depend on a first
QC pass having already selected which species to include - folding it back
into the gate would be circular. Treat this as a second, independent lens,
not a replacement for the main QC gate.

Method: reuses the pairwise ML distance matrices IQ-TREE already wrote
(mito.mldist / pltd.mldist - no new tree-building or alignment needed).
For each species present in both, find its nearest-neighbor SET in each
tree (all species tied at the minimum distance - both trees have many
exact ties, especially mito with only ~9 marker genes vs pltd's ~60, so a
single arbitrary nearest neighbor would be noisy). Compare the GENUS
membership of the two neighbor sets: if they share no genus at all, the
mito and plastid trees disagree about this species' closest relative -
flag it. A species with congeners sampled elsewhere in the dataset that
land in completely different, unrelated parts of the two trees is a real
red flag; note this can also arise from genuinely low taxon density around
a species (few/no close relatives sampled at all, in which case "nearest
neighbor" is a distant one on both trees and less informative) or, rarely,
real biology (hybrid origin, introgression) - treat as "worth a look",
same framing as the other cross-organelle checks in this suite, not an
automatic verdict.

Output: analysis/phylogeny/results/tree_discordance.tsv
Columns: species mito_nearest_genera pltd_nearest_genera
         mito_nearest_dist pltd_nearest_dist genus_overlap discordant
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
TREES_DIR = ANALYSIS_DIR / "phylogeny" / "results" / "trees"


def load_mldist(path: Path) -> tuple[list[str], list[list[float]]]:
    lines = path.read_text().splitlines()
    n = int(lines[0].strip())
    names, matrix = [], []
    for line in lines[1:1 + n]:
        parts = line.split()
        names.append(parts[0])
        matrix.append([float(x) for x in parts[1:1 + n]])
    return names, matrix


def nearest_neighbor_set(names: list[str], matrix: list[list[float]], i: int, eps: float = 1e-6) -> tuple[set, float]:
    """All j != i within eps of the minimum distance from i - handles the
    many exact ties in these matrices rather than picking one arbitrarily."""
    dists = [(j, d) for j, d in enumerate(matrix[i]) if j != i]
    if not dists:
        return set(), float("nan")
    min_d = min(d for _, d in dists)
    neighbors = {names[j] for j, d in dists if d <= min_d + eps}
    return neighbors, min_d


def genus(species: str) -> str:
    return species.split("_")[0]


def main():
    mito_path = TREES_DIR / "mito.mldist"
    pltd_path = TREES_DIR / "pltd.mldist"
    if not mito_path.exists() or not pltd_path.exists():
        print(f"[err] missing {mito_path} or {pltd_path} - run 04_build_tree.py for both organelles first",
              file=sys.stderr)
        sys.exit(1)

    mito_names, mito_matrix = load_mldist(mito_path)
    pltd_names, pltd_matrix = load_mldist(pltd_path)
    mito_idx = {n: i for i, n in enumerate(mito_names)}
    pltd_idx = {n: i for i, n in enumerate(pltd_names)}
    common = sorted(set(mito_names) & set(pltd_names))

    rows = []
    n_discordant = 0
    for sp in common:
        mito_nb, mito_d = nearest_neighbor_set(mito_names, mito_matrix, mito_idx[sp])
        pltd_nb, pltd_d = nearest_neighbor_set(pltd_names, pltd_matrix, pltd_idx[sp])
        mito_genera = {genus(n) for n in mito_nb}
        pltd_genera = {genus(n) for n in pltd_nb}
        overlap = mito_genera & pltd_genera
        discordant = len(overlap) == 0
        if discordant:
            n_discordant += 1
        rows.append({
            "species": sp,
            "mito_nearest_genera": ";".join(sorted(mito_genera)),
            "pltd_nearest_genera": ";".join(sorted(pltd_genera)),
            "mito_nearest_dist": round(mito_d, 6),
            "pltd_nearest_dist": round(pltd_d, 6),
            "genus_overlap": ";".join(sorted(overlap)),
            "discordant": discordant,
        })

    import pandas as pd
    out_path = ANALYSIS_DIR / "phylogeny" / "results" / "tree_discordance.tsv"
    columns = ["species", "mito_nearest_genera", "pltd_nearest_genera",
               "mito_nearest_dist", "pltd_nearest_dist", "genus_overlap", "discordant"]
    pd.DataFrame(rows, columns=columns).to_csv(out_path, sep="\t", index=False)
    print(f"[info] tree_discordance: {len(rows)} species compared, {n_discordant} discordant -> {out_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
