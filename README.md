# project-cpg-me

Comparative germline DNA methylation across invertebrates, estimated as CpG O/E from public
genome annotations. See [PLAN.md](PLAN.md) for the full plan.

## Setup (any machine: laptop, desktop, Raven)

Needs only [uv](https://docs.astral.sh/uv/). No conda or NCBI CLI required: data come from the
NCBI Datasets REST API.

```bash
uv sync --extra report
uv run pytest -q          # run on every new machine before it contributes results
```

## Run

```bash
uv run cpgme run GCF_013753865.1            # CDS only (no genome download)
uv run cpgme run GCF_013753865.1 --genome   # plus gene bodies and 1 kb upstream
uv run cpgme run-shard --shard 2/4          # this machine's share of manifest/species.tsv
uv run cpgme status                         # what is done / pending
uv run python analysis/pilot_report.py      # figures + summary table across species
```

## Functional annotation and the GO Slim test

After `run`, `annotate` writes one protein per gene (the longest CDS's protein), annotates it
with eggNOG-mapper, maps GO terms to GO Slim biological-process terms, and tests the 2015
prediction that housekeeping genes have lower CpG O/E than inducible genes.

```bash
# Where eggNOG-mapper and its database are installed (Raven):
uv run cpgme annotate GCF_013753865.1 --emapper-data /path/to/eggnog_data --cpu 24

# Elsewhere: write the proteins, annotate them on Raven, bring the result back.
uv run cpgme annotate GCF_013753865.1 --proteins-only     # -> results/<acc>/proteins.faa.gz
#   on Raven: emapper.py -i proteins.faa --itype proteins -m diamond --cpu 24 \
#             --data_dir /path/to/eggnog_data --output emapper --output_dir out/
uv run cpgme annotate GCF_013753865.1 --emapper-annotations out/emapper.emapper.annotations

uv run python analysis/goslim_report.py     # cross-species table + term-by-species figure
```

The housekeeping and inducible term sets are fixed in `cpgme/annotate.py` (chosen before any
results). goslim_generic lacks "response to stress", "translation" and "signal transduction",
so those three are added to the slim. GO files are cached in `data/_ontology/`; their release
is recorded in `annotation_provenance.json`.

Sharding is by hash of the accession, so machines need no coordination; finished accessions
are skipped unless the code changes. `data/` (downloads) and `results/` are not in git.

## Outputs, per accession (`results/<accession>/`)

| File | Contents |
|---|---|
| `features.tsv.gz` | one row per gene × feature (`cds`, `gene_body`, `upstream`): counts, GC, CpG O/E, GpC O/E, filter flag, posterior of the low-CpG O/E class |
| `summary.json` | 1- vs 2-component Gaussian mixture fits per feature (BIC, means, weights, threshold, Ashman's D) |
| `provenance.json` | inputs and md5s, code hash, git commit, host, package versions, timestamps |
| `proteins.faa.gz` | one protein per gene, header = gene_id (from `annotate`) |
| `annotation.tsv.gz` | per gene: eggNOG ortholog, OGs, COG, name, description, KEGG, GO, GO Slim |
| `goslim_cpg.tsv` | per GO Slim term × feature: n genes, median CpG O/E, share in low class, Mann-Whitney vs other annotated genes (BH q) |
| `goslim_test.json` | housekeeping vs inducible genes per feature: medians, share low, rank-biserial, one-sided p |
| `annotation_provenance.json` | emapper version/command, GO release, inputs md5 |

CpG O/E = (n_CG / (n_C × n_G)) × L² / (L − 1). Features < 300 bp or > 5 % ambiguous bases are
excluded from fits. Based on Dimond & Roberts 2015, *Mol Ecol* doi:10.1111/mec.13414.
