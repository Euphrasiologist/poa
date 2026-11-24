"""call_frac.py TAG ORG SPECIES... : every predicted edit of pilot/TAG/<species> (single-exon corrections) in ORG against the
species' pileups: read depth and the fraction of reads with the edit. Writes pilot/TAG.ORG.calls.tsv and prints classes."""
import sys, glob, pandas as pd
from pathlib import Path
H = Path(__file__).resolve().parent.parent; sys.path.insert(0, str(H)); import rna
tag, org, species = sys.argv[1], sys.argv[2], sys.argv[3:]
PILE = "pile" if org == "mito" else "pltd.pile"; rows = []
for sp in species:
    pfs = sorted(glob.glob(str(H / f"tol.{sp}.{PILE}.tsv.gz")) + glob.glob(str(H / f"pub.{sp}.*[0-9].{PILE}.tsv.gz")))
    f = H / "pilot" / tag / sp / "tables" / "predicted_edits.tsv"
    if not pfs or not f.exists(): continue
    P = pd.read_csv(f, sep="\t"); P = P[P.organelle == org]
    piles = [rna.load_pile(p) for p in pfs]
    for r in P.itertuples():
        n = c = t = g = a = 0
        for pile in piles:
            x = pile.get((r.contig, int(r.genomic_pos) + 1))
            if x: n += x[1]; c += x[2]["C"]; t += x[2]["T"]; g += x[2]["G"]; a += x[2]["A"]
        ref = piles[0].get((r.contig, int(r.genomic_pos) + 1), ("?",))[0]
        alt = t if ref == "C" else a if ref == "G" else 0
        rows.append((sp, r.gene.split(".")[0].split("@")[0], r.profile_column, r.resulting_aa, ref, n, alt))
D = pd.DataFrame(rows, columns=["species", "gene", "profile_column", "aa", "ref", "depth", "alt"]); D["frac"] = D.alt / D.depth.clip(lower=1)
D["cls"] = "no reads (<10)"; D.loc[D.depth >= 10, "cls"] = "unedited (<2%)"; D.loc[(D.depth >= 10) & (D.frac >= 0.02), "cls"] = "low (2-20%)"; D.loc[(D.depth >= 10) & (D.frac >= 0.2) & (D.alt >= 3), "cls"] = "edited (>=20%)"
D.to_csv(H / "pilot" / f"{tag}.{org}.calls.tsv", sep="\t", index=False)
print(tag, org, D.cls.value_counts().to_dict(), "ref bases", D.ref.value_counts().to_dict())
u = D[D.cls == "unedited (<2%)"]; print("unedited calls by gene:", u.gene.value_counts().head(25).to_dict())
print("recurrent unedited columns:", u.groupby(["gene", "profile_column", "aa"]).size().sort_values(ascending=False).head(30).to_dict())
