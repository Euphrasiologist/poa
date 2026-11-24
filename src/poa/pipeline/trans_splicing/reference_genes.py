"""The genes trans_splicing builds templates for and reconstructs - one
list (GENES, plus a lineage's own in LINEAGE_EXTRA_SPLICED) shared by fetch_references, exon_profiles, reconstruct, gff_export
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
# Genes with introns in one lineage that are single-exon in angiosperms.
# For such a lineage they are reconstructed from its exon templates like
# the rest, and editing leaves them alone (a single-exon correction across
# an intron reports the intron as frameshifts).
# Fern: atp1, cox1 and rps14 have one intron each in all eleven DToL ferns
# with RNA-seq. Horsetail: cox1 has three (Equisetum, same evidence).
# Lycophyte (Lycopodiaceae): atp6 and cob have two introns, atp9 and cox1
# three, cox3, nad3 and rps14 one, in all four RefSeq mitogenomes.
# Liverwort: atp9, nad3 and rps14 have one intron, atp1, cox3 and nad4L
# two, cob three and cox1 up to nine, in most of the RefSeq mitogenomes
# (atp1 has none in a third of them, nad4L one in a few).
LINEAGE_EXTRA_SPLICED = {
    "liverwort": ["atp1", "atp9", "cob", "cox1", "cox3", "nad3", "nad4L", "rps14"],
    "fern": ["atp1", "cox1", "rps14"],
    "horsetail": ["cox1"],
    "lycophyte": ["atp6", "atp9", "cob", "cox1", "cox3", "nad3", "rps14"],
}


# Genes of GENES a lineage is not given a template for. Liverwort nad7 is a
# pseudogene in every liverwort but Haplomitrium and Treubia: a template of
# those two drew 187 frameshift calls in 55 of 64 DToL liverworts.
LINEAGE_NOT_RECONSTRUCTED = {"liverwort": ["nad7"]}


def mito_spliced(lineage: str | None = None) -> list[str]:
    """The mito genes reconstructed for `lineage`; with None, for any."""
    if lineage is None:
        extra = sorted({g for genes in LINEAGE_EXTRA_SPLICED.values() for g in genes})
    else:
        extra = LINEAGE_EXTRA_SPLICED.get(lineage, [])
    left_out = LINEAGE_NOT_RECONSTRUCTED.get(lineage, [])
    return [g for g in GENES + [g for g in extra if g not in GENES] if g not in left_out]


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


# Lineages with plastid references of their own, used with --lineage on
# top of the land-plant plastid set: exon templates in
# reference/lineages/<x>/plastid/ and editing profiles in
# editing/reference/lineages/<x>/pltd/. Fern: plastid genes carry 300-900
# C-to-U and 30-80 U-to-C edits; against the land-plant references 59% of
# the edits poa called in 11 DToL ferns' plastid genes were in their RNA.
# Horsetail, lycophyte, moss and liverwort: the land-plant profiles are
# mostly angiosperm, whose proteins poa "corrected" their plastid genes
# towards - 713 edits called in 431 genes of 8 horsetails where the RNA
# shows 50, six of them among the calls. Horsetail and lycophyte
# references are RNA-edited genes of DToL assemblies, as for fern; moss
# and liverwort ones are RefSeq plastomes (LINEAGE_PLASTOMES).
LINEAGE_PLASTID = ("fern", "horsetail", "lycophyte", "moss", "liverwort")
# RefSeq plastomes, one per genus, for the lineages whose plastid
# references come from GenBank annotation: mosses and liverworts edit
# their plastid transcripts little or not at all, so the annotated
# proteins are the ones to align to. Not Takakia or Haplomitrium, which
# edit hundreds of sites and are annotated unedited.
MOSS_PLASTOMES = {
    "NC_060705.1": "Sphagnum multifibrosum",
    "NC_070057.1": "Andreaea regularis",
    "NC_060348.1": "Polytrichum commune",
    "NC_058541.1": "Atrichum angustatum",
    "NC_024291.1": "Tetraphis pellucida",
    "NC_046056.1": "Buxbaumia aphylla",
    "NC_046057.1": "Diphyscium foliosum",
    "NC_058544.1": "Funaria hygrometrica",
    "NC_037465.1": "Physcomitrium patens",
    "NC_081863.1": "Encalypta ciliata",
    "NC_044155.1": "Fissidens nobilis",
    "NC_080897.1": "Dicranum hengduanensis",
    "NC_050352.1": "Syntrichia filaris",
    "NC_069304.1": "Didymodon constrictus",
    "NC_085544.1": "Niphotrichum japonicum",
    "NC_050047.1": "Bartramia pomiformis",
    "NC_058542.1": "Bryum argenteum",
    "NC_056136.1": "Pohlia cruda",
    "NC_054293.1": "Mnium marginatum",
    "NC_058546.1": "Orthotrichum stellatum",
    "NC_058548.1": "Ulota hutchinsiae",
    "NC_058547.1": "Ptychomnion cygnisetum",
    "NC_051864.1": "Climacium dendroides",
    "NC_069303.1": "Anomodon rugelii",
    "NC_025668.1": "Sanionia uncinata",
    "NC_049070.1": "Calliergonella cuspidata",
}
LIVERWORT_PLASTOMES = {
    "NC_042505.1": "Marchantia polymorpha",
    "NC_042887.1": "Riccia fluitans",
    "NC_039590.1": "Dumortiera hirsuta",
    "NC_050394.1": "Wiesnerella denudata",
    "NC_088074.1": "Cyathodium cavernarum",
    "NC_042418.1": "Reboulia hemisphaerica",
    "NC_080368.1": "Pellia neesiana",
    "NC_019628.1": "Apopellia endiviifolia",
    "NC_058632.1": "Pallavicinia lyellii",
    "NC_080366.1": "Moerckia flotoviana",
    "NC_080365.1": "Fossombronia foveolata",
    "NC_082503.1": "Aneura maxima",
    "NC_082060.1": "Pleurozia gigantea",
    "NC_082430.1": "Porella gracillima",
    "NC_015402.1": "Ptilidium pulcherrimum",
    "NC_043781.1": "Radula japonica",
    "NC_079849.1": "Frullania muscicola",
    "NC_043782.1": "Jubula hutchinsiae",
    "NC_043777.1": "Cheilolejeunea xanthocarpa",
    "NC_043778.1": "Cololejeunea lanciloba",
    "NC_051002.1": "Scapania ampliata",
    "NC_053890.1": "Diplophyllum taxifolium",
    "NC_051003.1": "Douinia plicata",
    "NC_071973.1": "Calypogeia suecica",
    "NC_063926.1": "Nowellia curvifolia",
    "NC_040133.1": "Gymnomitrion concinnatum",
    "NC_043788.1": "Heteroscyphus argutus",
    "NC_043784.1": "Plagiochila chinensis",
    "NC_043785.1": "Bazzania praerupta",
    "NC_043779.1": "Schistochila macrodonta",
}
LINEAGE_PLASTOMES = {"moss": MOSS_PLASTOMES, "liverwort": LIVERWORT_PLASTOMES}


def reference_dir(code_dir, organelle: str = "mito", lineage: str | None = None):
    """trans_splicing/reference/ for a mito lineage (common/lineage), or
    reference/plastid/ - one set for every plastid, unless `lineage` has
    plastid references of its own (LINEAGE_PLASTID)."""
    from pathlib import Path

    from poa.pipeline.common import lineage as lin

    if organelle == "pltd":
        if lineage in LINEAGE_PLASTID:
            return lin.trans_splicing_reference(code_dir, lineage) / "plastid"
        return Path(code_dir) / "trans_splicing" / "reference" / "plastid"
    return lin.trans_splicing_reference(code_dir, lineage or lin.DEFAULT_LINEAGE)

# GenBank gene-name spellings -> ours (compared lower-cased)
GENE_SYNONYMS = {"ccmfc": "ccmFc", "ccmf-c": "ccmFc", "yejv": "ccmFc", "coxii": "cox2", "cox-2": "cox2"}

# RefSeq angiosperm mitogenomes annotated with the RNA-edited proteins:
# every record of the 451 RefSeq angiosperm mitogenomes (2026-11) in which
# at least 25 coding sequences carry an RNA-editing exception, 240 to 650
# edited residues each. One Beta of three; Asclepias is the one asterid
# and rice the one monocot. A few of their genes are annotated unedited
# all the same (beet ccmB, rice nad4 and cox2): those loci are left out
# gene by gene (common/refseq_mito.unedited).
ANGIOSPERM_EDITED_MITOGENOMES = {
    "NC_037304.1": "Arabidopsis thaliana",
    "NC_008285.1": "Brassica napus",
    "NC_002511.2": "Beta vulgaris",
    "NC_015994.1": "Beta macrocarpa",
    "NC_014487.1": "Silene latifolia",
    "NC_022796.1": "Asclepias syriaca",
    "NC_014043.1": "Citrullus lanatus",
    "NC_014050.1": "Cucurbita pepo",
    "NC_016743.2": "Lotus japonicus",
    "NC_016742.1": "Pongamia pinnata",
    "NC_069572.1": "Sedum plumbizincicola",
    "NC_030753.1": "Nelumbo nucifera",
    "NC_037468.1": "Nymphaea colorata",
    "NC_011033.1": "Oryza sativa",
}
# angiosperm RefSeq mitogenomes: the original five (nad1/2/5/rps3
# templates), six more spread across the tree for the cis-spliced genes,
# and the rest of the edited ones above. The records whose proteins are
# the genome's own translation (maize, grape, tobacco, soybean, carrot,
# Liriodendron, date palm) stay as a record of the exon structures; the
# templates and editing profiles are built from the edited ones alone
# (LINEAGE_UNEDITED_MITOGENOMES, LINEAGE_EDITED_MITOGENOMES): from the
# eleven together poa had the edited residue at 32 of 547 edited positions
# of five RefSeq genomes.
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
    **{acc: name for acc, name in ANGIOSPERM_EDITED_MITOGENOMES.items()},
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

# Liverwort (Marchantiophyta) RefSeq mitogenomes, one per genus across the
# three classes: Haplomitriopsida (Haplomitrium, Treubia), Marchantiopsida
# (complex thalloids) and Jungermanniopsida (simple thalloids and leafy
# liverworts). Liverwort genes keep an intron set of their own
# (LINEAGE_EXTRA_SPLICED; cox2 3 exons, nad2, nad4, nad5 and rpl2 2, nad1,
# rps3, rps10 and ccmFc 1), none of it trans-spliced; nad7 is a pseudogene
# in all but Haplomitrium and Treubia.
LIVERWORT_MITOGENOMES = {
    "NC_043889.1": "Haplomitrium hookeri",
    "NC_016122.1": "Treubia lacunosa",
    "NC_037508.1": "Marchantia polymorpha subsp. ruderalis",
    "NC_001660.1": "Marchantia paleacea",
    "NC_043906.1": "Riccia fluitans",
    "NC_042873.1": "Dumortiera hirsuta",
    "NC_053538.1": "Wiesnerella denudata",
    "NC_088073.1": "Cyathodium cavernarum",
    "NC_035978.1": "Asperifolia arguta",
    "NC_080363.1": "Pellia epiphylla",
    "NC_080360.1": "Apopellia endiviifolia",
    "NC_080362.1": "Pallavicinia lyellii",
    "NC_080361.1": "Moerckia flotoviana",
    "NC_043887.1": "Fossombronia foveolata",
    "NC_026901.1": "Aneura pinguis",
    "NC_043888.1": "Metzgeria furcata",
    "NC_013444.1": "Pleurozia purpurea",
    "NC_043908.1": "Porella platyphylla",
    "NC_043907.1": "Ptilidium ciliare",
    "NC_082157.1": "Scapania undulata",
    "NC_035979.1": "Calypogeia suecica",
    "NC_063925.1": "Nowellia curvifolia",
    "NC_054215.1": "Diplophyllum taxifolium",
    "NC_054214.1": "Douinia plicata",
    "NC_043909.1": "Solenostoma sphaerocarpum",
    "NC_040132.1": "Gymnomitrion concinnatum",
    "NC_037041.1": "Trilophozia quinquedentata",
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

# Ferns whose annotated proteins are the fully edited ones, from RNA-seq:
# Haplopteris (above) and four Salviniales whose mitochondrial genes are
# deposited one record per gene, 37 each (Kwok van der Giezen et al. 2025,
# Mol Biol Evol 42:msaf241). Each carries about 860 C-to-U and 440-780
# U-to-C edits; the other fern records' proteins have the stop-removing
# edits and nothing else. Editing profiles are built from these alone
# (editing/lineage_profiles.py): a profile from unedited proteins gives
# orfedit nothing to correct towards.
SALVINIALES_GENE_RECORDS = {
    "Azolla rubra": range(554530, 554567),
    "Azolla pinnata": range(554567, 554604),
    "Marsilea mutica": range(554604, 554641),
    "Salvinia molesta": range(554641, 554678),
}
FERN_EDITED_MITOGENOMES = {
    **{acc: name for acc, name in FERN_MITOGENOMES.items() if name == "Haplopteris ensiformis"},
    **{f"PQ{n}.1": name for name, numbers in SALVINIALES_GENE_RECORDS.items() for n in numbers},
}

# The fern records whose annotated proteins are not the edited ones (only
# their stop codons are edited out). They stay fetched, as a record of the
# exon structures, but are left out of the fern exon templates: with them
# in, each edited position's profile column is split between the edited
# and the unedited residue, and transsplice found 9% of the C-to-U edits
# RNA-seq shows in 11 DToL ferns' spliced genes, against 83% without.
FERN_UNEDITED_MITOGENOMES = tuple(acc for acc, name in FERN_MITOGENOMES.items() if name != "Haplopteris ensiformis")
LINEAGE_UNEDITED_MITOGENOMES = {
    "fern": FERN_UNEDITED_MITOGENOMES,
    "angiosperm": tuple(acc for acc in REFERENCE_MITOGENOMES if acc not in ANGIOSPERM_EDITED_MITOGENOMES),
}
# Horsetails (Equisetum, the one living genus) have no annotated
# mitogenome in GenBank: their references are all RNA-derived
# (reference/lineages/horsetail/raw_rna/README.md).
HORSETAIL_MITOGENOMES: dict[str, str] = {}
# Clubmosses (Lycopodiaceae): the four RefSeq mitogenomes, two subfamilies.
# Their genes keep an intron set of their own (LINEAGE_EXTRA_SPLICED; nad2
# 3 exons, nad4 2, cox2 4, rps3 3) and have lost the ccm genes and nad7.
# Their annotated proteins are unedited, so the templates and profiles
# are built from RNA-edited loci instead (lineages/lycophyte/README.md)
# and these stay as the record of the exon structures. Not Selaginella or
# Isoetes, whose mitogenomes are nothing like these.
LYCOPHYTE_MITOGENOMES = {
    "NC_017755.1": "Huperzia squarrosa",
    "NC_071971.1": "Huperzia crispata",
    "NC_086557.1": "Phylloglossum drummondii",
    "NC_080981.1": "Lycopodium japonicum",
}
# The clubmoss RefSeq proteins are unedited too (bar a few stop codons).
LINEAGE_UNEDITED_MITOGENOMES["lycophyte"] = tuple(LYCOPHYTE_MITOGENOMES)
LINEAGE_MITOGENOMES = {
    "angiosperm": REFERENCE_MITOGENOMES,
    "moss": MOSS_MITOGENOMES,
    "liverwort": LIVERWORT_MITOGENOMES,
    "fern": FERN_MITOGENOMES,
    "horsetail": HORSETAIL_MITOGENOMES,
    "lycophyte": LYCOPHYTE_MITOGENOMES,
}
# Lineages whose reference proteins are the RNA-edited ones throughout, so
# exon-slot profiles are cut from the protein, not translated from the
# genome (exon_profiles.py).
LINEAGE_EDITED_PROTEINS = ("angiosperm", "fern", "horsetail", "lycophyte")
# the genomes a lineage's editing profiles come from, where they are not
# all of LINEAGE_MITOGENOMES
LINEAGE_EDITED_MITOGENOMES = {
    "angiosperm": ANGIOSPERM_EDITED_MITOGENOMES,
    "fern": FERN_EDITED_MITOGENOMES,
    "lycophyte": {},
}
