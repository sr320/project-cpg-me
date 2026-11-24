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

## Next steps

1. On Raven: *Lycorma* gene bodies, and eggNOG-mapper on all seven proteomes, to test the
   GO Slim pattern from the 2015 paper (housekeeping genes low CpG O/E, inducible genes high).
2. Validate against measured methylation: *M. gigas* and *Apis* WGBS gene-level data.
3. Phase 1: build the full manifest from NCBI (all annotated genomes in the five groups).
