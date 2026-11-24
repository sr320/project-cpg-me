# Pilot results (2026-11-23)

One annotated RefSeq genome per taxon group plus two controls, run with `cpgme` on a laptop
(macOS arm64). Inputs and versions for each species are in `results/<accession>/provenance.json`.
Figures: [pilot_cds.png](pilot_cds.png), [pilot_gene_body.png](pilot_gene_body.png),
[pilot_upstream.png](pilot_upstream.png). Full table: [pilot_summary.tsv](pilot_summary.tsv).

## Summary (longest CDS per gene; gene body where computed)

| Group | Species | Genes | CDS median | CDS low class: share, mean | Ashman's D, CDS / gene body | Bimodal, CDS / gene body |
|---|---|---|---|---|---|---|
| coral | *Acropora millepora* | 29,970 | 0.71 | 40 %, 0.49 | 1.96 / 2.42 | no / yes |
| bivalve | *Magallana gigas* | 25,930 | 0.63 | 45 %, 0.42 | 2.74 / 3.53 | yes / yes |
| hemiptera | *Lycorma delicatula* | 13,877 | 0.51 | 72 %, 0.46 | 1.81 / not run | no / not run |
| tunicate | *Ciona robusta* | 14,101 | 0.79 | 49 %, 0.59 | 2.54 / 3.13 | yes / yes |
| sea star | *Acanthaster planci* | 18,389 | 0.65 | 53 %, 0.52 | 2.21 / 3.88 | yes / yes |
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

## Next steps

1. On Raven: *Lycorma* gene bodies (set `genome_features` to yes in the manifest).
2. A covariate check: does the housekeeping/inducible gap survive controlling for GC content,
   gene length and expression breadth?
3. Validate against measured methylation: *M. gigas* and *Apis* WGBS gene-level data.
4. Phase 1: build the full manifest from NCBI (all annotated genomes in the five groups).
