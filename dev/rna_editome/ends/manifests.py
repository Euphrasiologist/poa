"""manifests.py : from the poa runs in ends/runs, the orfedit manifest each run's final single-exon calls came from, pooled into
ends/sites.manifest.tsv (species corrected with the two editing-site costs) and ends/flat.manifest.tsv (one cost)."""
import glob
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent
out = {"sites": [], "flat": []}
for t in sorted(glob.glob(str(H / "runs" / "*" / "tables" / "lineage.tsv"))):
    run = Path(t).parents[1]; L = pd.read_csv(t, sep="\t", dtype=str).fillna("").iloc[0]
    assert L.lineage == "angiosperm", (run.name, L.lineage)  # one lineage's options only, below
    kind = "sites" if L.edit_sites == "yes" else "flat"
    m = run / "work/analysis/editing/work" / ("manifest.angiosperm.sites.tsv" if kind == "sites" else "manifest.angiosperm.tsv")
    out[kind].append(pd.read_csv(m, sep="\t", dtype=str)); print(run.name, kind, len(out[kind][-1]))
for k, v in out.items():
    pd.concat(v).to_csv(H / f"{k}.manifest.tsv", sep="\t", index=False); print(k, sum(map(len, v)), "rows")
