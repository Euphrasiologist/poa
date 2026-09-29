"""The genes trans_splicing builds templates for and reconstructs - one
list shared by 00a_fetch_references.py, 00_build_exon_profiles.py,
01_reconstruct.py and gff_export, so they can't drift apart.

transsplice itself doesn't care whether a junction is cis or trans (it
classifies each junction per species after assignment), so cis-spliced
genes are just more templates. Mito only for now; plastid cis-spliced
genes and bryophyte-specific intron sets are later stages.
"""
TRANS_SPLICED_MITO = ["nad1", "nad2", "nad5", "rps3"]
# angiosperm mito cis-spliced genes (stage 2a)
CIS_SPLICED_MITO = ["nad4", "nad7", "ccmFc", "cox2", "rpl2", "rps10"]
GENES = TRANS_SPLICED_MITO + CIS_SPLICED_MITO

# GenBank gene-name spellings -> ours (compared lower-cased)
GENE_SYNONYMS = {"ccmfc": "ccmFc", "ccmf-c": "ccmFc", "yejv": "ccmFc", "coxii": "cox2", "cox-2": "cox2"}

# angiosperm RefSeq mitogenomes: the original five (nad1/2/5/rps3 templates)
# plus six more spread across the tree for the cis-spliced genes
REFERENCE_MITOGENOMES = {
    "NC_037304.1": "Arabidopsis thaliana",
    "NC_002511.2": "Beta vulgaris",
    "NC_007982.1": "Zea mays",
    "NC_011033.1": "Oryza sativa",
    "NC_012119.1": "Vitis vinifera",
    "NC_006581.1": "Nicotiana tabacum",
    "NC_020455.1": "Glycine max",
    "NC_017855.1": "Daucus carota",
    "NC_021152.1": "Liriodendron tulipifera",
    "NC_016740.1": "Phoenix dactylifera",
    "NC_014487.1": "Silene latifolia",
}
