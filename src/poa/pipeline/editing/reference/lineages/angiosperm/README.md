# Angiosperm editing references from RNA-seq

`rna_proteins.fasta` and `rna_proteins.pltd.fasta` (headers `species|gene`) are single-exon mitochondrial and plastid
protein genes of Darwin Tree of Life angiosperm assemblies, each the gene poa called with the C-to-U editing the
species' own RNA-seq shows applied: 1349 mitochondrial proteins of 31 genes from 104 species, and 4198 plastid proteins of
56 genes from 106. The spliced genes of the same species (403 mitochondrial and 597 plastid loci) are in
`trans_splicing/reference/raw_rna/` and `trans_splicing/reference/plastid/raw_rna/`. With the RefSeq mitogenomes
annotated with edited proteins (`reference_genes.ANGIOSPERM_EDITED_MITOGENOMES`) they are what
`00b_build_lineage_profiles.py --lineage angiosperm` builds `profiles_refseq/mito/` from, and alone what
`--organelle pltd` builds `pltd/` from.

Residues in lower case are those a C-to-U edit made. From them each profile takes its editing sites (the `edit`
column; `editing/README.md`).

## Method

`dev/rna_editome/rna_refs.py`, run with `RNA_REVERSE=0`. The RNA-seq is public (ENA; the runs are below), one paired
run a species, mapped with bwa-mem2 to the species' own mitochondrial and plastid contigs
(`dev/rna_editome/map4.sh`). A site is edited where at least two reads and a fifth of the reads show U for the
genome's C. A gene is used when at least 90% of it has 10 or more reads and, with those edits applied, it
translates without an internal stop: for a single-exon gene the reading frame poa corrected, for a spliced gene the
exons transsplice placed, complete and without a frameshift call.

The reads are often not from the plant that was assembled, and where the two differ by a C/T the reads look edited.
No angiosperm edits U to C, so the sites where the reads show C for the genome's T count those differences. A species
with more than max(3, 5% of its C-to-U sites) of them is "noisy": its proteins are used, in upper case throughout,
and its loci are noted `sites not used` in `junctions.tsv`, so that they shape the profile but are no evidence of
where editing happens. One with more than max(50, half) is left out: the reads are another species'.
Mitochondrion: 71 ok, 40 noisy. Plastid: 75 noisy, 33 ok, 3 left out.

## Species and runs

| species | run | mito C-to-U / U-to-C sites | mito | plastid C-to-U / U-to-C sites | plastid |
|---|---|---|---|---|---|
| Acorus calamus | SRR16293892 | 113 / 4 | ok | 84 / 32 | noisy |
| Ailanthus altissima | SRR22105849 | 108 / 9 | noisy | 99 / 28 | noisy |
| Alliaria petiolata | SRR13765005 | 437 / 0 | ok | 29 / 1 | ok |
| Anchusa arvensis | ERR15799809 | 411 / 16 | ok | 54 / 20 | noisy |
| Arabidopsis thaliana | DRR095136 | 430 / 3 | ok | 26 / 1 | ok |
| Arctium lappa | SRR17633632 | 311 / 25 | noisy | 51 / 15 | noisy |
| Asparagus officinalis | SRR27137566 | 626 / 7 | ok | 166 / 52 | noisy |
| Bellardia viscosa | ERR13317821 | 48 / 31 | noisy | 93 / 36 | noisy |
| Berberis vulgaris | ERR18094426 | 467 / 15 | ok | 94 / 14 | noisy |
| Betula pendula | SRR15533024 | 59 / 10 | noisy | 141 / 50 | noisy |
| Brassica nigra | SRR11906213 | 163 / 0 | ok | 40 / 3 | ok |
| Bryonia dioica | SRR11573048 | 81 / 6 | noisy | 77 / 36 | noisy |
| Buxus sempervirens | ERR10378015 | 480 / 10 | ok | 102 / 1 | ok |
| Carex elata | ERR15551511 | 88 / 6 | noisy | 28 / 2 | ok |
| Centaurium intermedium | ERR17338211 | 188 / 16 | noisy | 40 / 6 | noisy |
| Ceratophyllum demersum | SRR35525968 | 107 / 8 | noisy | 20 / 4 | noisy |
| Chelidonium majus | SRR30505276 | 353 / 5 | ok | 50 / 15 | noisy |
| Claytonia sibirica | ERR15985871 | 151 / 17 | noisy | 31 / 1 | ok |
| Convolvulus arvensis | SRR37662152 | 513 / 0 | ok | 56 / 28 | noisy |
| Cornus suecica | ERR18094435 | 461 / 22 | ok | 47 / 6 | noisy |
| Crocosmia x crocosmiiflora | SRR5961547 | 249 / 11 | ok | 51 / 16 | noisy |
| Cruciata laevipes | ERR18094433 | 230 / 0 | ok | 23 / 5 | noisy |
| Daucus carota | SRR9217752 | 115 / 0 | ok | 23 / 5 | noisy |
| Diapensia lapponica | ERR17004075 | 506 / 3 | ok | 32 / 0 | ok |
| Dioscorea communis | SRR8819992 | 27 / 0 | ok | 59 / 17 | noisy |
| Drosera rotundifolia | SRR8948654 | 33 / 4 | noisy | 10 / 4 | noisy |
| Ecballium elaterium | SRR37300033 | 32 / 8 | noisy | 38 / 11 | noisy |
| Eryngium campestre | SRR37662193 | 629 / 22 | ok | 156 / 57 | noisy |
| Erythranthe guttata | SRR1556037 | 12 / 0 | ok | 38 / 4 | noisy |
| Euonymus europaeus | ERR11814143 | 265 / 20 | noisy | 46 / 20 | noisy |
| Euphorbia maculata | SRR34068846 | 53 / 6 | noisy | 68 / 19 | noisy |
| Filago lutescens | ERR12828821 | 230 / 16 | noisy | 48 / 3 | ok |
| Frangula alnus | ERR18094438 | 549 / 3 | ok | 58 / 0 | ok |
| Fraxinus excelsior | SRR36502356 | 156 / 29 | noisy | 85 / 22 | noisy |
| Fuchsia magellanica | ERR17792440 | 434 / 15 | ok | 8 / 0 | ok |
| Galium verum | SRR37662158 | 425 / 29 | noisy | 269 / 258 | left out |
| Geranium pyrenaicum | ERR1753274 | 9 / 0 | ok | 18 / 5 | noisy |
| Geranium robertianum | ERR1753296 | 5 / 1 | ok | 36 / 6 | noisy |
| Hedera helix | ERR9435019 | 784 / 94 | noisy | 640 / 629 | left out |
| Helianthemum oelandicum | ERR15951534 | 465 / 1 | ok | 14 / 0 | ok |
| Hibiscus verdcourtii | SRR32995397 | 299 / 4 | ok | 56 / 4 | noisy |
| Hordeum vulgare | SRR30479339 | 283 / 20 | noisy | 28 / 10 | noisy |
| Houttuynia cordata | SRR34656552 | 140 / 9 | noisy | 49 / 10 | noisy |
| Humulus lupulus | SRR12329022 | 218 / 1 | ok | 59 / 11 | noisy |
| Hypericum perforatum | ERR15005011 | 322 / 0 | ok | 31 / 14 | noisy |
| Ilex aquifolium | ERR12642428 | 309 / 48 | noisy | 154 / 39 | noisy |
| Impatiens parviflora | ERR15551475 | 429 / 5 | ok | 40 / 0 | ok |
| Juglans regia | SRR1068175 | 198 / 7 | ok | 15 / 5 | noisy |
| Juncus conglomeratus | ERR17532079 | 41 / 1 | ok | 34 / 3 | ok |
| Koelreuteria paniculata | SRR37056986 | 23 / 2 | ok | 20 / 8 | noisy |
| Lactuca saligna | SRR30764674 | 557 / 0 | ok | 41 / 8 | noisy |
| Lactuca serriola | SRR30764666 | 546 / 1 | ok | 37 / 0 | ok |
| Levisticum officinale | ERR15286565 | 367 / 8 | ok | 44 / 12 | noisy |
| Linum bienne | ERR18064660 | 201 / 1 | ok | 21 / 1 | ok |
| Lobelia urens | ERR14656571 | 404 / 11 | ok | 63 / 25 | noisy |
| Lolium perenne | SRR14655083 | 194 / 21 | noisy | 27 / 9 | noisy |
| Lythrum salicaria | ERR17338149 | 522 / 1 | ok | 40 / 0 | ok |
| Malus domestica | ERR12138253 | 311 / 13 | ok | 88 / 13 | noisy |
| Medicago sativa | SRR10742442 | 167 / 16 | noisy | 38 / 6 | noisy |
| Myrica gale | SRR22815379 | 233 / 19 | noisy | 66 / 24 | noisy |
| Myriophyllum alterniflorum | ERR15696450 | 565 / 0 | ok | 49 / 0 | ok |
| Narcissus pseudonarcissus | SRR5788585 | 97 / 10 | noisy | 106 / 32 | noisy |
| Nardus stricta | ERR15183986 | 45 / 6 | noisy | 12 / 6 | noisy |
| Nicandra physalodes | SRR31078684 | 215 / 35 | noisy | 101 / 30 | noisy |
| Nuphar pumila | ERR17004077 | 803 / 0 | ok | 92 / 0 | ok |
| Nymphaea alba | ERR17890722 | 794 / 0 | ok | 79 / 0 | ok |
| Nymphoides peltata | SRR1994392 | 150 / 16 | noisy | 53 / 19 | noisy |
| Oxyria digyna | SRR22045469 | 447 / 2 | ok | 33 / 0 | ok |
| Papaver somniferum | SRR8305583 | 86 / 11 | noisy | 57 / 14 | noisy |
| Paulownia tomentosa | SRR11787961 | 521 / 9 | ok | 163 / 55 | noisy |
| Pentaglottis sempervirens | ERR15551453 | 394 / 4 | ok | 34 / 3 | ok |
| Phragmites australis | SRR3233387 | 376 / 8 | ok | 41 / 10 | noisy |
| Platanthera chlorantha | SRR27827629 | 197 / 7 | ok | 42 / 16 | noisy |
| Platanus x hispanica | ERR15286576 | 844 / 0 | ok | 91 / 0 | ok |
| Polygala serpyllifolia | ERR13731911 | 75 / 7 | noisy | 43 / 4 | noisy |
| Pontederia cordata | ERR17229639 | 648 / 0 | ok | 62 / 0 | ok |
| Populus alba | SRR6411249 | 279 / 20 | noisy | 62 / 8 | noisy |
| Populus nigra | SRR23436602 | 332 / 4 | ok | 48 / 6 | noisy |
| Portulaca oleracea | SRR10241278 | 166 / 4 | ok | 44 / 3 | ok |
| Potamogeton perfoliatus | ERR15820896 | 102 / 10 | noisy | 82 / 38 | noisy |
| Quercus cerris | SRR22815343 | 413 / 3 | ok | 57 / 1 | ok |
| Ranunculus repens | ERR17468095 | 307 / 5 | ok | 45 / 5 | noisy |
| Reseda luteola | ERR17338141 | 536 / 1 | ok | 35 / 0 | ok |
| Rhodiola rosea | SRR5936537 | 175 / 1 | ok | 5 / 4 | noisy |
| Ribes rubrum | ERR15696453 | 605 / 0 | ok | 69 / 5 | noisy |
| Salix alba | SRR40750569 | 61 / 17 | noisy | 45 / 19 | noisy |
| Sambucus nigra | ERR15140864 | 618 / 25 | ok | 103 / 28 | noisy |
| Saponaria officinalis | SRR25745560 | 170 / 4 | ok | 34 / 16 | noisy |
| Scheuchzeria palustris | ERR16121729 | 620 / 11 | ok | 70 / 5 | noisy |
| Scrophularia umbrosa | ERR17229716 | 280 / 9 | ok | 62 / 24 | noisy |
| Solanum crispum | SRR18909654 | 69 / 22 | noisy | 76 / 32 | noisy |
| Spirodela polyrhiza | SRR14319991 | 566 / 19 | ok | 36 / 25 | noisy |
| Suaeda vera | SRR12121641 | 362 / 4 | ok | 33 / 6 | noisy |
| Succisa pratensis | ERR18094450 | 266 / 6 | ok | 67 / 19 | noisy |
| Tellima grandiflora | ERR17532086 | 667 / 2 | ok | 65 / 0 | ok |
| Teucrium chamaedrys | SRR29428738 | 378 / 9 | ok | 88 / 43 | noisy |
| Trifolium subterraneum | DRR144433 | 101 / 6 | noisy | 29 / 3 | ok |
| Typha angustifolia | SRR15541138 | 139 / 10 | noisy | 63 / 20 | noisy |
| Ulmus procera | ERR17762815 | 371 / 1 | ok | 47 / 0 | ok |
| Urtica dioica | ERR17204400 | 198 / 3 | ok | 42 / 18 | noisy |
| Vaccinium myrtillus | SRR7523932 | 342 / 58 | noisy | 59 / 56 | left out |
| Valeriana pyrenaica | ERR17854984 | 292 / 22 | noisy | 38 / 1 | ok |
| Valerianella locusta | ERR16727270 | 265 / 7 | ok | 35 / 7 | noisy |
| Veronica verna | ERR18094440 | 78 / 3 | ok | 29 / 3 | ok |
| Vicia sativa | SRR11809687 | 86 / 1 | ok | 40 / 3 | ok |
| Vinca minor | ERR13802782 | 83 / 10 | noisy | 65 / 17 | noisy |
| Viola tricolor | SRR15458583 | 165 / 1 | ok | 33 / 14 | noisy |
| Viscum album | ERR10968077 | 40 / 1 | ok | 36 / 8 | noisy |
| Vitis vinifera | SRR30692097 | 55 / 4 | noisy | 112 / 23 | noisy |
| Xerophyllum tenax | ERR17229642 | 630 / 10 | ok | 98 / 7 | noisy |
| Zostera marina | SRR16285855 | 24 / 4 | noisy | 1 / 1 | ok |
