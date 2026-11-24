# Pilot results (2026-11-23)

Two annotated RefSeq genomes per taxon group plus two controls. The first seven species were
run with `cpgme` on a laptop (macOS arm64); the second species per group and *Lycorma* gene
bodies on Raven. Inputs and versions for each species are in `results/<accession>/provenance.json`.
Figures: [pilot_cds.png](pilot_cds.png), [pilot_gene_body.png](pilot_gene_body.png),
[pilot_upstream.png](pilot_upstream.png). Full table: [pilot_summary.tsv](pilot_summary.tsv).

## Summary (longest CDS per gene; gene body where computed)

| Group | Species | Genes | CDS median | CDS low class: share, mean | Ashman's D, CDS / gene body | Bimodal, CDS / gene body |
|---|---|---|---|---|---|---|
| coral | *Acropora millepora* | 29,970 | 0.71 | 40 %, 0.49 | 1.96 / 2.42 | no / yes |
| coral | *Stylophora pistillata* | 34,954 | 0.69 | 38 %, 0.46 | 2.03 / 2.73 | yes / yes |
| bivalve | *Magallana gigas* | 25,930 | 0.63 | 45 %, 0.42 | 2.74 / 3.53 | yes / yes |
| bivalve | *Crassostrea virginica* | 28,340 | 0.65 | 42 %, 0.42 | 2.70 / 3.76 | yes / yes |
| hemiptera | *Lycorma delicatula* | 13,877 | 0.51 | 72 %, 0.46 | 1.81 / 2.43 | no / yes |
| hemiptera | *Nilaparvata lugens* | 17,648 | 0.75 | 35 %, 0.58 | 1.31 / 0.42 | no / no |
| tunicate | *Ciona robusta* | 14,101 | 0.79 | 49 %, 0.59 | 2.54 / 3.13 | yes / yes |
| tunicate | *Styela clava* | 19,222 | 0.82 | 43 %, 0.60 | 2.52 / 3.23 | yes / yes |
| sea star | *Acanthaster planci* | 18,389 | 0.65 | 53 %, 0.52 | 2.21 / 3.88 | yes / yes |
| sea star | *Asterias rubens* | 15,938 | 0.67 | 50 %, 0.55 | 1.54 / 2.09 | no / yes |
| control (methylated) | *Apis mellifera* | 9,726 | 0.84 | 48 %, 0.53 | 2.96 / 3.30 | yes / yes |
| control (unmethylated) | *Drosophila melanogaster* | 13,168 | 0.88 | n/a | 0.93 / 0.92 | no / no |

Bimodality rule (adopted 2026-11-23): BIC favours two components by more than 10 **and**
Ashman's D > 2. Ashman's D measures how well the two fitted modes separate. The "low class"
columns come from the two-component fit whether or not the distribution is called bimodal.

## What the pilot shows

1. **The controls behave as expected.** *Drosophila* has one mode near 0.9 in CDS and gene
   bodies (D about 0.9). *Apis* is clearly bimodal in both.
2. **Every focal species with a genome is bimodal in gene bodies.** *Magallana*, *Ciona* and
   *Acanthaster* are also bimodal in CDS. *Acropora* CDS falls just under the cut-off
   (D 1.96): the two classes are there but overlap. That is consistent with the bimodal coral
   pattern of Dimond & Roberts 2015, which used transcripts, but the cut-off makes it
   borderline. Only the bimodality was checked; the functional partitioning needs GO
   annotation, which has not been run.
3. **Gene bodies separate the two classes more cleanly than CDS in every species**
   (higher Ashman's D). Intronic CpG depletion carries extra signal. CDS stays the
   cross-species unit (no genome needed, comparable across tiers), but gene bodies should be
   computed wherever a genome exists. Report D as a continuous measure, not only the yes/no
   call.
4. ***Lycorma* is the outlier.** It has the lowest median CpG O/E, and 72 % of genes fall in
   the low class. Its CDS distribution is one large low mode with a long high tail rather than
   two separated modes, so it is not called bimodal (D 1.81). This could reflect heavy, broad
   gene-body methylation, or an AT-biased mutation spectrum (its GpC O/E median is about 1.0,
   which argues against a pure composition effect). With one species this is a lead, not a
   result. Gene bodies (on Raven) and other Hemiptera, including non-invasive ones, are needed
   to tell lineage from invasiveness.
   *Update:* in gene bodies *Lycorma* is bimodal (D 2.43; low class 51 % at 0.41), like the
   other focal species. Its CDS shape is a CDS-level effect, not an absent second class.
5. **Upstream (1 kb) regions are not useful as a methylation proxy here.** No species is
   bimodal under the adopted rule. Under BIC alone, *Drosophila* upstream looked "bimodal"
   (ΔBIC 305). That reflects promoter CpG islands and intergenic sequence, not germline
   methylation history.

## Bimodality rule

With more than 10,000 genes, BIC favours two components for almost any non-Gaussian shape.
The rule now also requires Ashman's D > 2 (`cpgme/fit.py`, with a unit test on a skewed
single-peak distribution). `summary.json` keeps the BIC-only call as `bimodal_bic` for
comparison.

## Operational notes

- Compute is trivial: all seven species re-run in about 1 minute from cached downloads.
  Downloads dominate, at roughly 1 GB per minute from NCBI.
- One NCBI download (*Drosophila*) hung for over 30 minutes without the socket timeout
  firing. `fetch.py` now has a watchdog that aborts and retries after 120 s with no progress
  (unit-tested; not yet tested against a live stall).
- Disk: about 4.2 GB of downloads for the seven species; can be deleted after a run.
  The 2.2 Gb *Lycorma* genome was left to Raven as planned.

## GO Slim test: *Acropora millepora* (2026-11-23)

`cpgme annotate` on Raven (eggNOG-mapper 2.1.13, eggNOG DB 5.0.2, diamond 2.0.11; GO
release 2026-07-26; 40 cores, about 6 minutes). Figures: [goslim_gene_body.png](goslim_gene_body.png),
[goslim_cds.png](goslim_cds.png). Table: [goslim_tests.tsv](goslim_tests.tsv).

| Feature | Housekeeping median (n) | Inducible median (n) | Low class, HK / IND | Rank-biserial | p (HK lower) |
|---|---|---|---|---|---|
| CDS | 0.57 (924) | 0.68 (5,174) | 59 % / 42 % | −0.25 | 5 × 10⁻³⁵ |
| Gene body | 0.61 (942) | 0.77 (5,196) | 59 % / 35 % | −0.36 | 5 × 10⁻⁷⁰ |

1. **The 2015 pattern holds.** Genes in housekeeping GO Slim terms have lower CpG O/E (more
   historical germline methylation) than genes in inducible terms, in CDS and more strongly
   in gene bodies.
2. **The term ranking agrees.** Of 71 terms with ≥ 20 genes (gene bodies), all twelve
   housekeeping terms rank below the median (mean rank 10); the lowest four are
   cytoplasmic translation, mitochondrial gene expression, translation and ribosome
   biogenesis. Inducible terms average rank 55, led by extracellular matrix organization,
   inflammatory response and cell adhesion. Terms outside both sets also behave sensibly:
   DNA repair and autophagy sit low; nervous, circulatory, muscle and endocrine system
   processes sit high.
3. **Caveats.** One species. Only 36 % of genes have GO terms from eggNOG (10,995 of 30,136),
   so the test rests on the better-characterised genes. About 1,700 genes carry both a
   housekeeping and an inducible term and are excluded. The 2015 paper used transcriptomes
   and an older GO Slim, so this matches its pattern, not its exact numbers.

Phase 0 criterion "reproduce the bimodal distribution and GO Slim ordering for *Acropora*":
met for the GO Slim ordering in both features; bimodality is met in gene bodies and
borderline in CDS (above).

## GO Slim test: all pilot species (2026-11-23)

Same run settings on Raven. *Lycorma* is CDS only (gene bodies not yet run). Rank-biserial
< 0 means housekeeping genes rank lower; p is one-sided (housekeeping lower).

| Group | Species | GO coverage | Feature | HK median | IND median | Low class, HK / IND | Rank-biserial | p |
|---|---|---|---|---|---|---|---|---|
| coral | *Acropora millepora* | 36 % | CDS | 0.57 | 0.68 | 59 % / 42 % | −0.25 | 5 × 10⁻³⁵ |
| | | | gene body | 0.61 | 0.77 | 59 % / 35 % | −0.36 | 5 × 10⁻⁷⁰ |
| bivalve | *Magallana gigas* | 46 % | CDS | 0.41 | 0.64 | 83 % / 45 % | −0.47 | 9 × 10⁻¹¹⁰ |
| | | | gene body | 0.48 | 0.70 | 84 % / 46 % | −0.39 | 7 × 10⁻⁷⁹ |
| hemiptera | *Lycorma delicatula* | 66 % | CDS | 0.48 | 0.52 | 86 % / 77 % | −0.10 | 7 × 10⁻⁷ |
| tunicate | *Ciona robusta* | 62 % | CDS | 0.69 | 0.80 | 63 % / 49 % | −0.16 | 1 × 10⁻¹² |
| | | | gene body | 0.68 | 0.79 | 72 % / 52 % | −0.22 | 7 × 10⁻²⁴ |
| sea star | *Acanthaster planci* | 59 % | CDS | 0.57 | 0.64 | 75 % / 55 % | −0.24 | 4 × 10⁻³³ |
| | | | gene body | 0.49 | 0.54 | 89 % / 64 % | −0.23 | 2 × 10⁻³¹ |
| control (methylated) | *Apis mellifera* | 73 % | CDS | 0.60 | 0.88 | 76 % / 46 % | −0.33 | 5 × 10⁻⁴⁴ |
| | | | gene body | 0.81 | 1.25 | 78 % / 44 % | −0.39 | 7 × 10⁻⁶³ |
| control (unmethylated) | *Drosophila melanogaster* | 74 % | CDS | 0.90 | 0.87 | 49 % / 61 % | +0.18 | 1 |
| | | | gene body | 0.95 | 0.93 | 28 % / 32 % | +0.12 | 1 |

"Low class" shares come from each species' two-component fit, whether or not the
distribution is called bimodal.

1. **The pattern holds in every species with gene-body methylation.** All five focal species
   and *Apis* show lower CpG O/E in housekeeping genes, in both features where both were run.
   Effect sizes are moderate (rank-biserial −0.10 to −0.47), largest in *Magallana* and *Apis*.
2. **The negative control passes.** *Drosophila*, which lacks gene-body methylation, shows no
   housekeeping-lower gap; housekeeping genes are slightly *higher* (rank-biserial +0.12 to
   +0.18). So the split in the other species is not produced by GO annotation, gene length or
   codon usage alone, at least not in a way that also shows up in a fly.
3. ***Lycorma* fits the pattern, weakly.** The gap exists (p = 7 × 10⁻⁷) but is the
   smallest (rank-biserial −0.10), and both classes sit mostly in the low class (86 % / 77 %).
   This is consistent with its unusual distribution above: CpG depletion is broad, reaching
   inducible genes too. Gene bodies are still needed.
4. **GO coverage differs** (36 % in *Acropora* to 74 % in *Drosophila*), so the set of genes
   tested is not the same across species. Orthology-based comparisons (Phase 4) will help.
5. **The term figure** ([goslim_gene_body.png](goslim_gene_body.png)) shows the same ordering
   across species: translation, RNA processing, DNA replication and cell-cycle terms low;
   adhesion, extracellular matrix, signalling and organ-system terms high; *Drosophila* flat.

## Second species per group (2026-11-23)

*Stylophora pistillata*, *Crassostrea virginica*, *Nilaparvata lugens*, *Styela clava* and
*Asterias rubens* (manifest notes give the reason for each), plus *Lycorma* gene bodies. Same
settings on Raven. GO Slim test, housekeeping (HK) vs inducible (IND):

| Group | Species | GO coverage | Feature | HK median | IND median | Low class, HK / IND | Rank-biserial | p |
|---|---|---|---|---|---|---|---|---|
| coral | *Stylophora pistillata* | 35 % | CDS | 0.52 | 0.65 | 63 % / 41 % | −0.28 | 3 × 10⁻⁴⁴ |
| | | | gene body | 0.57 | 0.76 | 69 % / 39 % | −0.39 | 2 × 10⁻⁸⁹ |
| bivalve | *Crassostrea virginica* | 45 % | CDS | 0.42 | 0.65 | 81 % / 43 % | −0.47 | 2 × 10⁻¹¹⁶ |
| | | | gene body | 0.44 | 0.70 | 83 % / 44 % | −0.41 | 2 × 10⁻⁸⁹ |
| hemiptera | *Lycorma delicatula* | 66 % | gene body | 0.42 | 0.58 | 82 % / 49 % | −0.37 | 3 × 10⁻⁶⁸ |
| hemiptera | *Nilaparvata lugens* | 57 % | CDS | 0.79 | 0.79 | 12 % / 14 % | +0.04 | 0.96 |
| | | | gene body | 0.74 | 0.73 | 8 % / 6 % | +0.05 | 1 |
| tunicate | *Styela clava* | 56 % | CDS | 0.66 | 0.82 | 66 % / 42 % | −0.27 | 6 × 10⁻³⁸ |
| | | | gene body | 0.71 | 0.91 | 73 % / 43 % | −0.36 | 1 × 10⁻⁶³ |
| sea star | *Asterias rubens* | 61 % | CDS | 0.59 | 0.67 | 67 % / 50 % | −0.23 | 3 × 10⁻²⁸ |
| | | | gene body | 0.59 | 0.61 | 88 % / 70 % | −0.16 | 5 × 10⁻¹⁵ |

1. **The pattern replicates in four of the five new species.** *Stylophora*, *C. virginica*,
   *Styela* and *Asterias* are bimodal in gene bodies and show lower CpG O/E in housekeeping
   genes, with effect sizes in the range of their first-species counterparts. *C. virginica*
   and *M. gigas* are nearly identical (rank-biserial −0.47 in CDS for both).
2. ***Nilaparvata lugens* looks like the unmethylated control.** One narrow mode (gene body
   D 0.42), no housekeeping/inducible difference, and few genes in a "low" class. Unlike
   *Drosophila*, the mode sits well below 1 (median 0.72 vs 0.92), while GpC O/E is about 1.0,
   so CpGs are depleted fairly evenly across genes rather than in a subset. Possible readings:
   loss of gene-body methylation in this lineage with residual depletion, or methylation
   spread uniformly across genes. Published *N. lugens* WGBS would distinguish these.
3. **The two planthoppers differ sharply.** *Lycorma* (Fulgoridae) has two classes and a
   strong functional split in gene bodies; *Nilaparvata* (Delphacidae) has neither. Within
   Hemiptera, the CpG O/E structure is lineage-dependent, so invasion contrasts here must be
   made between close relatives.
4. ***Asterias* is the weakest of the methylated species** (CDS not bimodal, ΔBIC < 0; gene
   body D 2.09; gap −0.16), against a strong pattern in *Acanthaster* (D 3.88).

## Next steps

1. WGBS validation (running on Raven): *M. gigas* gonad (SRR1045858) and eight female gonads
   (PRJNA806944), *Apis* sperm (GSM1897585). Then *N. lugens* WGBS to explain its flat profile.
2. A covariate check: does the housekeeping/inducible gap survive controlling for GC content,
   gene length and expression breadth?
3. Phase 1: build the full manifest from NCBI (all annotated genomes in the five groups).

## Five more bivalves (2026-11-23)

Added *Ostrea edulis*, *Mytilus edulis*, *Pecten maximus*, *Mercenaria mercenaria* and
*Dreissena polymorpha* (chromosome-level RefSeq reference genomes; CDS + gene bodies, run on a
laptop, about 6 minutes for all five including downloads). Bivalves now span Ostreida,
Mytilida, Pectinida and Venerida (incl. Dreissenidae).

| Species | Genes | CDS median | CDS low class: share, mean | Ashman's D, CDS / gene body | Bimodal, CDS / gene body |
|---|---|---|---|---|---|
| *Ostrea edulis* | 27,586 | 0.66 | 42 %, 0.42 | 2.75 / 3.26 | yes / yes |
| *Mytilus edulis* | 35,671 | 0.60 | 36 %, 0.29 | 2.80 / 3.53 | yes / yes |
| *Pecten maximus* | 25,959 | 0.68 | 42 %, 0.46 | 1.93 / 2.89 | no / yes |
| *Mercenaria mercenaria* | 38,088 | 0.72 | 33 %, 0.40 | 2.84 / 3.62 | yes / yes |
| *Dreissena polymorpha* | 36,148 | 0.80 | 45 %, 0.49 | 3.10 / 3.60 | yes / yes |

1. **All seven bivalves are bimodal in gene bodies** (D 2.9–3.8), and six of seven in CDS.
   Bivalves are the most consistently bimodal group in the set.
2. ***Pecten* CDS is the one miss** (D 1.93, BIC favours one component), like *Acropora* CDS:
   the gene-body call is clear, so this is overlap of the two CDS classes, not a missing class.
3. ***Mytilus* has the deepest low class** (CDS low mean 0.29 vs 0.40–0.49 in the others),
   suggesting the most heavily methylated gene set (strongest historical CpG depletion).
4. ***Dreissena* is shifted up** (CDS median 0.80, high mode at 1.0, GpC O/E 0.96) yet still
   cleanly bimodal (D 3.1). The whole distribution is higher, not flatter; whether that is
   composition or less methylation needs WGBS or the GpC comparison across groups.
5. The oysters (*Magallana*, *Crassostrea*, *Ostrea*) are nearly indistinguishable
   (CDS low mean 0.42 in all three).
6. Upstream regions: no new bivalve is bimodal under the adopted rule.

GO Slim annotation for these five is pending: `proteins.faa.gz` is written for each, but the
local eggNOG database (`data/_eggnog/eggnog.db.gz`) is a truncated download, so emapper must
run on Raven.

## Two non-invasive *Lycorma* relatives (2026-11-23)

Added *Laodelphax striatellus* (small brown planthopper, Delphacidae; native Eurasian rice
pest) and *Macrosteles quadrilineatus* (aster leafhopper, Cicadellidae; native to North
America). These are the only annotated Auchenorrhyncha genomes on NCBI besides *Lycorma*,
*Nilaparvata* and the invasive *Homalodisca vitripennis*. Four were requested; other
Fulgoroidea genomes (e.g. *Catonia carolina*, *Cyphopterum fauveli*) are unannotated.

| Species | Genes | CDS median | Ashman's D, CDS / gene body | Bimodal, CDS / gene body |
|---|---|---|---|---|
| *Lycorma delicatula* (invasive) | 13,877 | 0.51 | 1.81 / 2.43 | no / yes |
| *Nilaparvata lugens* | 17,648 | 0.75 | 1.31 / 0.42 | no / no |
| *Laodelphax striatellus* | 15,487 | 0.80 | 1.09 / 0.11 | no / no |
| *Macrosteles quadrilineatus* | 21,708 | 0.62 | 0.50 / 0.23 | no / no |

1. **None of the three non-invasive relatives is bimodal**, in CDS or gene bodies.
   *Lycorma* is the only auchenorrhynchan with two separated gene-body classes, and it also
   has by far the lowest median CpG O/E.
2. The two delphacids (*Nilaparvata*, *Laodelphax*) are very similar: CDS median 0.75–0.80,
   gene-body D below 0.5. *Macrosteles* has a lower median (0.62) but still one broad mode.
3. This sharpens the *Lycorma* lead but does not separate invasiveness from lineage:
   *Lycorma* is the only fulgorid, and the relatives come from two other families. An
   annotated non-invasive fulgorid would be the decisive comparison.

**Parser change.** *Laodelphax* has a GenBank (submitter) annotation with no NCBI GeneIDs.
`features.py` now falls back to `locus_tag` when a record has no GeneID. All 19 species
were re-run: features are byte-identical for every RefSeq species except 13 *Drosophila*
mitochondrial genes, whose IDs changed from `gene:` to `locus_tag:` with identical values;
mixture summaries differ only by floating-point noise (relative ≤ 1e-13).
