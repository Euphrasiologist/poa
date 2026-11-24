"""ab_calls.py A B : single-exon calls of the same species in two sets of poa runs (pilot/A, pilot/B), call by call (same hit):
how many changed their start, their end, their aligned columns; whether they end at a stop; frame changes and internal stops;
by lineage called. Writes pilot/A.vs.B.calls.tsv."""
import sys, glob
import pandas as pd
A, B = sys.argv[1:3]; key = ["species", "organelle", "gene", "contig", "hit_start", "hit_end"]
def load(tag):
    out = []
    for f in glob.glob(f"pilot/{tag}/*/tables/edited_gene_calls.tsv"):
        c = pd.read_csv(f, sep="\t"); L = pd.read_csv(f.replace("edited_gene_calls", "lineage"), sep="\t").lineage.iloc[0]; c["lineage"] = L; out.append(c)
    d = pd.concat(out).drop_duplicates(key)
    if "n_internal_stops" not in d: d["n_internal_stops"] = 0
    return d.set_index(key)
a, b = load(A), load(B); j = a.join(b, lsuffix="_a", rsuffix="_b", how="inner")
print(f"{len(a)} calls in {A}, {len(b)} in {B}, {len(j)} of the same hit; only in {A}: {len(a.index.difference(b.index))}, only in {B}: {len(b.index.difference(a.index))}")
j["start_moved"] = (j.corrected_start_a != j.corrected_start_b); j["end_moved"] = (j.corrected_end_a != j.corrected_end_b)
j["dcol"] = j.n_columns_aligned_b - j.n_columns_aligned_a; j["dlen"] = (j.corrected_end_b - j.corrected_start_b) - (j.corrected_end_a - j.corrected_start_a)
j["full_a"] = j.n_columns_aligned_a >= 0.99 * j.n_columns_total_a; j["full_b"] = j.n_columns_aligned_b >= 0.99 * j.n_columns_total_b
g = j.groupby(["lineage_b", "organelle"]).agg(calls=("dcol", "size"), span_changed=("dlen", lambda x: int(((j.loc[x.index, "start_moved"]) | (j.loc[x.index, "end_moved"])).sum())),
        stop_a=("includes_stop_a", "sum"), stop_b=("includes_stop_b", "sum"), full_a=("full_a", "sum"), full_b=("full_b", "sum"), cols_up=("dcol", lambda x: int((x > 0).sum())), cols_down=("dcol", lambda x: int((x < 0).sum())),
        fs_a=("n_frameshifts_a", "sum"), fs_b=("n_frameshifts_b", "sum"), istops_b=("n_internal_stops_b", "sum"), genes_istop_b=("n_internal_stops_b", lambda x: int((x > 0).sum())), edits_a=("n_edits_a", "sum"), edits_b=("n_edits_b", "sum"))
pd.set_option("display.width", 250); print(g.to_string())
big = j[j.dcol <= -10].reset_index()
print(f"\ncalls with 10 or more fewer columns aligned: {len(big)}"); print(big.groupby(["lineage_b", "organelle", "gene"]).size().sort_values(ascending=False).head(25).to_string())
up = j[j.dcol >= 10].reset_index(); print(f"\ncalls with 10 or more columns more: {len(up)}"); print(up.groupby(["lineage_b", "organelle", "gene"]).size().sort_values(ascending=False).head(12).to_string())
s = j[j.n_internal_stops_b > 0].reset_index(); print(f"\ninternal stops by gene:"); print(s.groupby(["lineage_b", "organelle", "gene"]).agg(genes=("gene", "size"), stops=("n_internal_stops_b", "sum")).sort_values("genes", ascending=False).head(30).to_string())
j.reset_index().drop(columns=["protein_a", "protein_b"]).to_csv(f"pilot/{A}.vs.{B}.calls.tsv", sep="\t", index=False)
