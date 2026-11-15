"""How much of an assembly graph a resolved contig set accounts for.

Used to gate gfatk-resolve promotions (`linearize/src/05_promote_resolve.py`,
`07_restore_pathfinder.py`): a circuit that walks only part of the graph
(Cyperaceae plastids came out 34-55 kb of a ~150-230 kb genome) must not
replace anything. Placement uses `unitig_coords`' sequence matching, so a
number here means the same thing as `pct_unitig_content_placed` in
`unitig_coords/results/linearization_qc.tsv`.
"""

from __future__ import annotations

from pathlib import Path

from poa.pipeline.unitig_coords import unitig_map as _um


def unitig_content_placed(gfa: Path, ctg_fasta: Path) -> float:
    """Fraction (0-1) of the graph's unique unitig bp found in `ctg_fasta`
    (forward or reverse complement; circular joins and trimmed repeat
    variants handled as in unitig_coords/unitig_map.py)."""
    segs, _ = _um.read_gfa(gfa)
    _, ctgs = _um.read_ctg_fasta(ctg_fasta)
    total = sum(len(s) for s in segs.values())
    if not total or not ctgs:
        return 0.0
    doubled = {cid: c + c for cid, c in ctgs.items()}
    placed = 0
    for useq in segs.values():
        max_trim = min(_um.MAX_TRIM_BP, len(useq) // 4)
        rc = _um.revcomp(useq)
        if any(
            _um.find_all_with_tolerance(doubled[cid], len(c), needle, max_trim)[0]
            for cid, c in ctgs.items()
            for needle in (useq, rc)
        ):
            placed += len(useq)
    return placed / total
