"""The genes trans_splicing builds templates for and reconstructs - one
list shared by fetch_references, exon_profiles, reconstruct, gff_export
and editing's lineage_profiles, so they can't drift apart.

transsplice itself doesn't care whether a junction is cis or trans (it
classifies each junction per species after assignment), so cis-spliced
genes are just more templates. Mito only for now; plastid cis-spliced
genes are a later stage.

Templates come per lineage (LINEAGE_MITOGENOMES): "angiosperm" lives in
reference/{raw,profiles}/ as it always has, any other lineage in
reference/lineages/<name>/{raw,profiles}/. Which lineage a species uses is
decided once, by editing.scan (common/lineage.py).
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

# moss (Bryophyta) RefSeq mitogenomes, spread across the phylum: Sphagnum,
# Andreaea, Polytrichopsida, Tetraphis, Buxbaumia and the Bryopsida orders.
# Moss mitogenomes keep one intron set across the phylum - nad1 3 exons,
# nad2 2, nad4 2, nad5 4, nad7 3, ccmFc 2, cox2 4 (3 in a few), rps3 and
# rpl2 1, rps10 lost - where angiosperm templates expect 5/5/4/5/5/2/2/2/2,
# so the angiosperm templates can never fill a moss gene.
MOSS_MITOGENOMES = {
    "NC_024521.1": "Sphagnum palustre",
    "NC_079917.1": "Sphagnum subsecundum",
    "NC_068856.1": "Andreaea regularis",
    "NC_046759.1": "Andreaea wangiana",
    "NC_024520.1": "Atrichum angustatum",
    "NC_039775.1": "Polytrichum commune",
    "NC_024290.1": "Tetraphis pellucida",
    "NC_024518.1": "Buxbaumia aphylla",
    "NC_007945.1": "Physcomitrella patens",
    "NC_024523.1": "Funaria hygrometrica",
    "NC_027515.1": "Syntrichia filaris",
    "NC_029452.1": "Racomitrium lanuginosum",
    "NC_050308.1": "Chorisodontium aciphyllum",
    "NC_036945.1": "Mielichhoferia elongata",
    "NC_024519.1": "Bartramia pomiformis",
    "NC_024522.1": "Orthotrichum stellatum",
    "NC_024516.1": "Hypnum imponens",
    "NC_016121.1": "Anomodon rugelii",
    "NC_024515.1": "Climacium americanum",
    "NC_031212.1": "Brachythecium rivulare",
}

LINEAGE_MITOGENOMES = {"angiosperm": REFERENCE_MITOGENOMES, "moss": MOSS_MITOGENOMES}
