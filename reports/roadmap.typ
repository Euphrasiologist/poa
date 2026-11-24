#import "@local/report:0.1.0": *
#show: report.with(
  title: "poa, orfedit, transsplice and gfatk: roadmap",
  subtitle: "State on 19 November 2026 and what remains",
  abstract: [
    All four tools have unpushed commits from this week. The annotation
    side is in good shape for seed plants and mosses and is known to be
    wrong for ferns, which need reverse RNA editing before anything else
    will help them. `gfatk resolve` is now reproducible and measures depth
    correctly, but its single-circuit model leaves much of a fern
    mitochondrial graph out. The dataset is half-updated and should be
    re-run once, after the decisions listed in section 2.
  ],
  toc: false,
)

= Where things stand

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, left, left),
    table.header([*repository*], [*head (local)*], [*what is unpushed*]),
    [gfatk], [`5fa62d6`], [0.6.3: overlaps written once; reproducible `resolve`; read-guided walks through repeats; links followed on either strand in `linear`, `path`, `trim`, `extract`; depth from well-aligned reads. Tag 0.6.3 on this commit.],
    [orfedit], [`03a704c`], [0.2.0: frameshift windows; faster alignment; `--stop-search`.],
    [transsplice], [`ed0474c`], [0.3.0: frameshifts per exon; staged parallel batch; short exons placed by whole-gene score; gene ends held to start and stop codons.],
    [poa], [`506759c`], [Spliced-gene frameshifts, plastid spliced genes, rebuilt exon templates with floors at a quarter of the reference self-score, rebuilt frameshift reference.],
  ),
  caption: [Local commits. Nothing has been pushed; pushing, tagging and bioconda are the user's to do.],
)

The new exon templates in poa need the new transsplice: with transsplice
before `ed0474c` they drop nad5 exon 3. poa's `environment.yml`,
`install.sh` and recipe still pin orfedit 0.1.0, transsplice 0.2.1 and
gfatk 0.6.2.

= Decisions needed

+ *The 28 dataset assemblies whose re-resolved circuit fails the promotion
  gate.* They still hold sequences made with gfatk before 0.6.2, which walk
  junctions the graphs do not have. Options: return them to oatk's
  original files (what the pipeline's own rule says), promote the valid
  but partial circuits, or leave them until `resolve` can place more of
  those graphs.
+ *The 26 assemblies already rewritten.* They were written from the build
  before the depth fix. With the depth fix 25 of the 54 pass the gate, not
  26 (_Allium vineale_ falls to 70% placed), and _Silene uniflora_ comes
  out as five circuits. They should be regenerated with the tagged 0.6.3
  and the gate applied again before anything downstream is re-run.
+ *Exon-slot floors.* A quarter of the lowest reference self-score is
  committed. Lower is better on every measure tried: at 0.075 the 955 mito
  assemblies give 7,552 complete spliced genes and 779 with a frameshift
  call, against 7,025 and 1,039 at 0.25, and held-out placement is
  unchanged (426 of 477 exons). Proposed: set 0.075, after checking the
  plastid and moss runs at that value.

= Work, in order

== 1. Reverse editing in orfedit #badge("IN PROGRESS", fill: accent)

Ferns edit U to C as well as C to U. In the Polypodiales reference
mitogenomes the genome carries about nine in-frame stop codons per gene
that editing removes (428 in 46 genes in _Woodwardia japonica_). orfedit
models C to U only, so each of those reads as a broken gene, and the
rules added this week that end a gene at the first in-frame stop cut fern
genes short. The change: allow U-to-C edits in the codon hypotheses, so
that TAA and TAG can read as Gln and TGA as Arg, switched on by lineage;
with it on, the stop-finding rules in orfedit and transsplice must not
take an editable stop as a gene's end. Validation: the fern references'
own annotated proteins, and the dataset's 35 fern assemblies before and
after; seed plants and mosses must come out unchanged with it off.

== 2. Fern and lycophyte lineages

A leptosporangiate fern template set from the five Polypodiales
references (_Dryopteris crassirhizoma_, three _Woodwardia_, _Blechnopsis
orientalis_), with their exon counts: nad1 and nad5 have four exons, rps3
three, rps10 one, and ccmFc is absent. A separate lycophyte set from the
four Lycopodiaceae references (nad2 three exons, nad4 two, cox2 four, no
nad7). Each needs a rule in the lineage caller and held-out validation.
Horsetails have no annotated reference and stay on whichever set fits
best. Until this is done the README should say that fern proteins are not
reliable.

== 3. Release chain

+ gfatk: push, tag 0.6.3 on `5fa62d6`, bioconda.
+ orfedit 0.2.0 to crates.io and bioconda (after step 1, or as 0.2.0 now
  and 0.3.0 with reverse editing).
+ transsplice: `cargo update -p orfedit`, release 0.3.0.
+ poa: bump the three pins, then release; the bioconda recipe for poa
  itself is still to be submitted, and needs the user's confirmation
  before a pull request is opened.

== 4. Dataset re-run

After decisions 1 and 2 and the releases: regenerate the gfatk-resolved
assemblies, then run the dataset's downstream tables once (QC,
annotation, the poa suite), then rebuild poa's frameshift reference with
those assemblies included. Two things in the dataset repository to know
about: I regenerated `resolve_summary.tsv` from the current targets file,
and someone else was running `linearize` there the same day.

The dataset's `_resolve_one.sh` recruits any read with a primary hit to
the graph and aligns with GraphAligner's `vg` preset. Two-thirds of the
reads it recruits to a mitochondrial graph are only partly organellar.
gfatk now ignores those for depth, but filtering at recruitment and
trying the `dbg` preset are both untested and worth a look.

== 5. gfatk `resolve` on fern mitochondrial graphs

The valid circuits place 30 to 86% of the 28 failing graphs. The missing
sequence is real (three-quarters of it does not match the circuit) and
reachable; what keeps it out is that each junction may be crossed once
and each segment only as often as its depth allows. Relaxing both brings
most of it in but roughly doubles the circuit. Either the copy numbers
are still too low for these graphs or the genome is not one circle; the
depth fix makes the first testable now. Smaller open items from the code
review: connectivity is rooted in one component of the graph, a segment
named in a link but absent panics, `stats` reports "circular" from a test
that ignores orientation and doubles overlap totals when links are
written both ways. `resolve` also takes about twice as long since the
depth fix, for a reason not yet looked at.

== 6. Annotation loose ends

- ccmFc carries a frameshift call in 49% of mito assemblies and nad7
  calls rose from 37 to 56; neither is explained.
- Genes whose exon count differs from the template (intron loss) cannot
  be reconstructed; one candidate cannot fill two slots.
- orfedit's corrected coordinates can now run a few codons past the last
  residue of its `protein` column.
- The planted-indel benchmarks for spliced genes have not been re-run
  since the exon-placement changes.
- The editing stage can still call frameshifts inside spliced genes.
- The whole-gene profiles of nad1, nad2, nad5 and rps3 are RefSeq-only
  again and could be refined from the new reconstructions.
- A full mito transsplice run is 17 minutes at 16 threads, up from 13.
