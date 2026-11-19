"""The genes trans_splicing builds templates for and reconstructs - one
list shared by fetch_references, exon_profiles, reconstruct, gff_export
and editing's lineage_profiles, so they can't drift apart.

transsplice itself doesn't care whether a junction is cis or trans (it
classifies each junction per species after assignment), so cis-spliced
genes are just more templates. Mito only for now; plastid cis-spliced
genes have their own template set (PLASTID_SPLICED, reference/plastid/).

Templates come per lineage (LINEAGE_MITOGENOMES): "angiosperm" lives in
reference/{raw,profiles}/ as it always has, any other lineage in
reference/lineages/<name>/{raw,profiles}/. Which lineage a species uses is
decided once, by editing.scan (common/lineage.py).
"""

TRANS_SPLICED_MITO = ["nad1", "nad2", "nad5", "rps3"]
# angiosperm mito cis-spliced genes (stage 2a)
CIS_SPLICED_MITO = ["nad4", "nad7", "ccmFc", "cox2", "rpl2", "rps10"]
GENES = TRANS_SPLICED_MITO + CIS_SPLICED_MITO

# plastid genes with introns (rps12: exon 1 trans-spliced to exons 2-3).
# Not petB, petD, rpl16: their first exons are 2-3 codons, too short to
# template or find, and editing corrects their main exon as a single-exon
# gene.
# One template set for every land plant: plastid introns are far more
# stable than mito ones; a species lacking one (grass clpP/rpoC1, Pinus
# ndh genes) just doesn't fill a slot.
PLASTID_SPLICED = ["atpF", "clpP", "ndhA", "ndhB", "rpl2", "rpoC1", "rps12", "rps16", "ycf3"]
PLASTID_SYNONYMS = {"clpp1": "clpP", "pafi": "ycf3", "ycf3": "ycf3"}
# RefSeq plastomes across land plants: eudicots, monocots, magnoliids,
# Amborella, a conifer, a fern, a moss and a liverwort
PLASTOMES = {
    "NC_000932.1": "Arabidopsis thaliana",
    "NC_001879.2": "Nicotiana tabacum",
    "NC_002202.1": "Spinacia oleracea",
    "NC_007957.1": "Vitis vinifera",
    "NC_007942.1": "Glycine max",
    "NC_001320.1": "Oryza sativa",
    "NC_001666.2": "Zea mays",
    "NC_008326.1": "Liriodendron tulipifera",
    "NC_005086.1": "Amborella trichopoda",
    "NC_001631.1": "Pinus thunbergii",
    "NC_004766.1": "Adiantum capillus-veneris",
    "NC_005087.2": "Physcomitrium patens",
    "NC_001319.1": "Marchantia polymorpha",
}


def reference_dir(code_dir, organelle: str = "mito", lineage: str | None = None):
    """trans_splicing/reference/ for a mito lineage (common/lineage), or
    reference/plastid/ - one set for every plastid."""
    from pathlib import Path

    from poa.pipeline.common import lineage as lin

    if organelle == "pltd":
        return Path(code_dir) / "trans_splicing" / "reference" / "plastid"
    return lin.trans_splicing_reference(code_dir, lineage or lin.DEFAULT_LINEAGE)

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

# Leptosporangiate fern mitogenomes with annotation. None is RefSeq: five
# Polypodiales (the three Woodwardia, Blechnopsis and Dryopteris single
# circles are assembler output nobody has validated), Osmunda, and
# Haplopteris ensiformis, the one fern mitogenome with a tested structure
# (nine circular chromosomes, one accession each) and an RNA-seq editome.
# Fern genes differ from the angiosperm templates in exon count - nad1 4
# exons, nad5 4 or 5, rps3 3, rps10 1, no ccmFc - and none of their group
# II introns is trans-spliced. Their annotated proteins are the edited
# ones in both directions (U-to-C as well as C-to-U), which
# fetch_references allows for when it checks them against the genome.
FERN_MITOGENOMES = {
    "MW732172.1": "Dryopteris crassirhizoma",
    "PV999925.1": "Woodwardia japonica",
    "PV999926.1": "Woodwardia unigemmata",
    "PV999927.1": "Woodwardia prolifera",
    "PQ143019.1": "Blechnopsis orientalis",
    "PQ202832.1": "Osmunda japonica",
    **{f"OM8675{n}.1": "Haplopteris ensiformis" for n in range(45, 54)},
}

LINEAGE_MITOGENOMES = {"angiosperm": REFERENCE_MITOGENOMES, "moss": MOSS_MITOGENOMES, "fern": FERN_MITOGENOMES}
