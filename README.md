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

Sharding is by hash of the accession, so machines need no coordination; finished accessions
are skipped unless the code changes. `data/` (downloads) and `results/` are not in git.

## Outputs, per accession (`results/<accession>/`)

| File | Contents |
|---|---|
| `features.tsv.gz` | one row per gene × feature (`cds`, `gene_body`, `upstream`): counts, GC, CpG O/E, GpC O/E, filter flag, posterior of the low-CpG O/E class |
| `summary.json` | 1- vs 2-component Gaussian mixture fits per feature (BIC, means, weights, threshold, Ashman's D) |
| `provenance.json` | inputs and md5s, code hash, git commit, host, package versions, timestamps |

CpG O/E = (n_CG / (n_C × n_G)) × L² / (L − 1). Features < 300 bp or > 5 % ambiguous bases are
excluded from fits. Based on Dimond & Roberts 2015, *Mol Ecol* doi:10.1111/mec.13414.
