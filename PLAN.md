# project-cpg-me: plan

Comparative germline methylation (CpG O/E) across invertebrate taxa, computed from public
sequence data on whatever machine is free, merged into one dataset, and tested against
resilience and invasion traits.

Starting point: Dimond & Roberts 2015 (Mol Ecol, doi:10.1111/mec.13414). Across six coral
transcriptomes, CpG O/E was bimodal; low-CpG O/E (historically methylated) genes were enriched
for housekeeping functions, high-CpG O/E (sparsely methylated) genes for inducible functions,
and environmentally responsive genes sat in the high-CpG O/E class. This project asks whether
that structure, measured the same way in many species, carries signal about how a species
copes with new or changing environments.

## 1. Questions

1. How do the shape of the CpG O/E distribution (bimodality, fraction of low-CpG O/E genes,
   separation of the two modes) and its functional partitioning vary within and among taxa?
2. For shared orthologs, which genes keep the same methylation class across taxa and which
   switch?
3. Do species-level CpG O/E features covary with resilience or invasion traits once phylogeny
   is accounted for?

Question 3 is a hypothesis test, not a given. "Prediction" is the last phase and only
happens if 1 and 2 show a signal worth modelling.

## 2. Taxa

| Group | Seed species (genome-annotated first) | Contrast available |
|---|---|---|
| Bivalves | *Crassostrea gigas*, *C. virginica*, *Ostrea lurida*, *Mytilus galloprovincialis*, *M. edulis*, *Dreissena polymorpha*, *D. rostriformis*, *Corbicula fluminea*, *Mercenaria*, *Pecten*, *Panopea* | global invaders vs. range-restricted natives |
| Corals | *Acropora* spp., *Pocillopora*, *Stylophora pistillata*, *Montipora capitata*, *Porites* spp., *Orbicella faveolata*; the six species from the 2015 paper for replication | bleaching-susceptible vs. tolerant |
| Lanternfly relatives | *Lycorma delicatula* plus other Hemiptera: planthoppers (*Nilaparvata*, *Laodelphax*, *Sogatella*), *Halyomorpha halys*, *Bemisia tabaci*, aphids, *Homalodisca* | invasive pests vs. non-invasive congeners |
| Tunicates | *Ciona intestinalis/robusta*, *Botryllus schlosseri*, *Botrylloides*, *Styela clava*, *Didemnum vexillum*, *Halocynthia*, *Oikopleura* | several of the best-known marine invaders |
| Sea stars | *Acanthaster planci*, *Asterias rubens*, *A. amurensis*, *Pisaster ochraceus*, *Pycnopodia helianthoides*, *Patiria miniata* | outbreak/invasive species; wasting-disease susceptible vs. resistant |

Two things to settle before building the manifest:

- **"Flies."** The spotted lanternfly (*Lycorma delicatula*; invasive in Japan, Korea, and the
  US) is a planthopper (Hemiptera: Fulgoridae), not a dipteran. Its relatives are Hemiptera,
  which have functional DNMT1/DNMT3 and gene body methylation. True flies (Diptera) have lost
  most of that machinery and show a single CpG O/E mode near 1. The plan uses Hemiptera as the
  group and keeps two or three Diptera (*Drosophila*, *Aedes*) as a negative control.
  Confirmed 2026-11-23: Hemiptera is the focal insect group.
- **Outgroup/controls.** One well-characterised species per mode: *Apis mellifera* and
  *Nematostella vectensis* (clear bimodality), *Drosophila* (none).

## 3. Data sources

Priority order per species; record which tier was used, because tiers are not interchangeable.

1. **Annotated genome** (NCBI Datasets, RefSeq then GenBank; Ensembl Metazoa as a second
   source). Gives CDS, mRNA, and full gene bodies including introns.
2. **Transcriptome assembly** (NCBI TSA, or the assemblies behind published papers). This is
   what the 2015 paper used. Coral TSA assemblies need symbiont and bacterial contigs removed.
3. **Raw RNA-seq only** (SRA). Deferred; assembling is the only expensive step in the whole
   project and should be used only to fill important gaps.

Validation data: published WGBS/MBD-BS gene-level methylation where it exists (*C. gigas*,
*C. virginica*, *Acropora*, *Stylophora*, *Montipora*, *Pocillopora*, *Ciona*, *Nilaparvata*,
*Apis*) to check how well CpG O/E tracks measured methylation in each taxon.

Trait data: WRiMS/WoRMS and GISD (introduced/invasive status), NEMESIS, OBIS/GBIF (range
size, latitudinal and thermal breadth), EPPO/CABI for insects, coral trait database
(bleaching response), literature for thermal limits and disease susceptibility.

## 4. Per-species computation

One command per accession, same on every machine, no cross-species dependencies.

```
cpgme run GCF_902806645.1        # fetch -> extract -> compute -> fit -> write
```

Steps:

1. Fetch genome, GFF, CDS, protein with `datasets download genome accession`; verify checksums.
2. Build feature sets: longest-isoform CDS per gene, mRNA, full gene body, introns, and 1 kb
   upstream. CDS is the cross-species comparable unit; gene body is the biologically closer one.
3. Compute per feature: length, counts of C, G, CG, GC, then

   `CpG O/E = (n_CG / (n_C * n_G)) * (L^2 / (L - 1))`

   plus GpC O/E as a control for base-composition effects, and GC content. Drop features
   under 300 bp or with more than 5% N.
4. Fit 1- and 2-component Gaussian mixtures; record BIC for each, component means, SDs,
   mixing proportion, the intersection threshold, and a per-gene posterior for the low class.
5. Functional annotation: eggNOG-mapper on the protein set (GO, KEGG, COG), mapped to GO Slim
   as in the 2015 paper.
6. Write outputs (below) and a provenance file.

Output contract, one directory per accession:

```
results/<accession>/
  features.tsv.gz    gene_id, feature, length, n_C, n_G, n_CG, n_GC, gc, cpg_oe, gpc_oe, p_low
  summary.json       n_genes, mixture fits, BIC delta, frac_low, mode means, threshold
  annotation.tsv.gz  gene_id, GO, GO Slim, KEGG, COG, eggNOG OG
  provenance.json    accession, source tier, tool versions, git commit, host, date, checksums
```

Everything downstream reads only these files. Sequence data is never committed.

## 5. Running on many machines

The per-species job is small (minutes, a few GB of RAM for most genomes), so the design goal
is zero coordination, not throughput.

- **Manifest.** `manifest/species.tsv` in git is the single source of truth: taxon group,
  species, taxid, accession, tier, status of trait data. Adding a species is a pull request.
- **Sharding instead of locking.** Each machine runs `cpgme run-shard --shard k/N`, which
  takes the accessions whose hash falls in its slice. No machine needs to know what the others
  are doing, and re-running is harmless because a job is skipped if `provenance.json` exists
  with matching manifest row and code version.
- **Same environment everywhere.** One locked conda/pixi environment (Linux x86-64 and
  macOS arm64), so laptops and Raven produce identical numbers. A test fixture (one small
  genome with known output) runs on each machine before it contributes.
- **Division of labour.** Participating machines (decided 2026-11-23): hete, chlo, raven and
  klone, one shard each (`manifest/machines.tsv`, `--shard k/4`).
  - hete, chlo, klone: fetch plus CpG O/E plus mixture fit for standard-size genomes. klone
    (UW Hyak) runs through SLURM and needs Duo 2FA, so its jobs are submitted by hand.
  - Raven: large genomes (some bivalves and Hemiptera are multi-Gb), eggNOG-mapper,
    OrthoFinder, any transcriptome assembly, and the final merge.
- **Collecting results.** Each machine rsyncs `results/<accession>/` to Gannet (decided
  2026-11-23; path in `manifest/machines.tsv`). Raven builds the merged tables from there.
  Outputs are small enough that the merged tables can also be versioned as release assets.
- **Status.** `cpgme status` compares manifest against the shared results directory and
  prints what is done, missing, or stale.

## 6. Merge and comparative analysis

1. **Merged tables** (Parquet): all per-gene rows; one species-level row per accession.
2. **Orthology.** OrthoFinder across all proteomes, with BUSCO metazoa single-copy genes as a
   conservative core set. Each orthogroup gets a vector of methylation class across species.
3. **Species-level features.** Bimodality (BIC delta, dip test), fraction low, mode
   separation, median CpG O/E of stress-response/immune/housekeeping GO Slim sets, and the
   CpG O/E of genes in expanded families.
4. **Similarities.** Orthogroups with conserved class across all five groups; test whether
   the housekeeping/inducible split from corals holds everywhere.
5. **Differences.** Orthogroups that switch class between clades or between contrasted
   species within a clade; lineage-specific shifts in the whole distribution.
6. **Expression link.** Where public stress-response RNA-seq exists, test whether
   differentially expressed genes sit in the high-CpG O/E class, as in the 2015 paper.

## 7. Trait association and prediction

- Build `traits/species_traits.tsv` with a source citation for every value.
- Test association with phylogenetic methods (PGLS, phylogenetic logistic regression) on a
  time-calibrated tree (TimeTree plus published phylogenies). Species are not independent
  data points; uncorrected correlations across five phyla would mostly recover phylogeny.
- Prefer within-clade paired contrasts (invasive vs. non-invasive congeners; tolerant vs.
  susceptible corals) over pooled cross-phylum models.
- Only then fit predictive models, evaluated by leaving out whole clades, and report them
  with the effect sizes from the contrasts rather than as a standalone classifier.

## 8. Known limits

- CpG O/E records methylation over evolutionary time in the germline. It is not the current
  methylome and says nothing directly about plastic, environmentally induced methylation.
- Genome and transcriptome tiers give different values for the same species; compare within
  tier, and use species with both to quantify the offset.
- Annotation quality varies a lot. Fragmented gene models and missing genes shift the
  distribution; BUSCO completeness goes into every species row as a covariate and a filter.
- GC content and mutation bias differ across lineages; GpC O/E is carried through for this.
- Invasion status depends partly on opportunity (shipping, aquaculture) and on which species
  get sequenced. Invaders are over-represented among sequenced genomes.
- The number of independent contrasts per group will be small. Some groups may only support
  descriptive results.

## 9. Phases

| Phase | Work | Done when |
|---|---|---|
| 0. Pilot | Implement `cpgme run`; process two species per group plus controls; reproduce the bimodal distribution and GO Slim ordering for *Acropora* from the 2015 paper | Pilot species match published patterns; controls behave as expected |
| 1. Manifest | Query NCBI for all annotated genomes in the five groups; add TSA assemblies for gaps; review list | Manifest merged, with tier and BUSCO threshold per row |
| 2. Distributed run | Lock environment, validate each machine on the fixture, run shards, collect results | `cpgme status` shows all rows done |
| 3. Validation | Compare CpG O/E with published gene-level methylation per group | Correlation and misclassification rate reported per taxon |
| 4. Merge | Orthology, merged tables, species-level features | Parquet tables and a descriptive report |
| 5. Traits | Assemble trait table; phylogeny; contrasts and PGLS | Association results with effect sizes |
| 6. Prediction | Clade-held-out models, if phase 5 supports it | Written assessment of what is and is not predictable |

## 10. Repository layout

```
manifest/     species.tsv, machines.tsv
cpgme/        fetch, extract, compute, fit, annotate, merge, status
envs/         locked environment
tests/        fixture genome and expected output
traits/       species_traits.tsv with sources
analysis/     Quarto notebooks for phases 3-6
results/      gitignored; synced to shared storage
```

## 11. Decisions needed

Settled:
- Hemiptera is the insect group, with Diptera as negative control.
- Merged results live on Gannet (2026-11-23). Exact path still to be filled in
  `manifest/machines.tsv`.
- Participating machines are hete, chlo, raven and klone (2026-11-23): four shards, with
  large genomes, eggNOG-mapper, OrthoFinder and the merge on Raven.

Open:
1. Whether transcriptome-only species are in scope for the first pass or deferred.
