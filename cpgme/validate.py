"""Validate CpG O/E against measured gene-level methylation (WGBS).

Input is a Bismark CpG coverage file (`*.bismark.cov.gz` from
`bismark_methylation_extractor --bedGraph`; columns chrom, start, end, percent,
n_methylated, n_unmethylated; 1-based). Reads are aligned to the same assembly as
`cpgme run`, so gene coordinates come straight from its GFF.

Per gene (protein-coding gene span from the GFF):
  n_cpg       CpG cytosines covered by at least MIN_COV reads (both strands count)
  meth        weighted methylation level: sum(methylated) / sum(total) over those CpGs
  frac_mcpg   share of those CpGs with methylation >= 0.5
A gene with fewer than MIN_CPGS covered CpGs gets no value.

Comparison, per CpG O/E feature (cds, gene_body):
  Spearman rho between CpG O/E and meth (negative is expected)
  AUC of low CpG O/E for predicting a methylated gene (meth >= METH_THRESHOLD)
  Misclassification of the low-CpG O/E class (p_low > 0.5) against methylated genes

Outputs (results/<accession>/validation/<label>/):
  gene_methylation.tsv.gz   gene_id, n_cpg, meth, frac_mcpg
  validation.json           statistics above, inputs and thresholds
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

from .features import protein_coding_genes

MIN_COV = 5
MIN_CPGS = 5
METH_THRESHOLD = 0.5


def read_coverage(path: Path, min_cov: int = MIN_COV) -> dict[str, tuple[np.ndarray, ...]]:
    """{chrom: (positions sorted, n_meth, n_total)} for CpGs with total >= min_cov."""
    pos, meth, tot = defaultdict(list), defaultdict(list), defaultdict(list)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("track"):
                continue
            f = line.split("\t")
            m, u = int(f[4]), int(f[5])
            if m + u < min_cov:
                continue
            pos[f[0]].append(int(f[1]))
            meth[f[0]].append(m)
            tot[f[0]].append(m + u)
    out = {}
    for chrom in pos:
        p = np.array(pos[chrom], dtype=np.int64)
        order = np.argsort(p, kind="stable")
        out[chrom] = (p[order], np.array(meth[chrom], dtype=np.int64)[order],
                      np.array(tot[chrom], dtype=np.int64)[order])
    return out


def gene_methylation(cov: dict[str, tuple[np.ndarray, ...]], genes_by_seq: dict[str, list[dict]],
                     min_cpgs: int = MIN_CPGS) -> dict[str, dict]:
    """Per-gene methylation over the gene span (1-based, closed)."""
    out: dict[str, dict] = {}
    for chrom, genes in genes_by_seq.items():
        if chrom not in cov:
            continue
        p, m, t = cov[chrom]
        cm = np.concatenate([[0], np.cumsum(m)])
        ct = np.concatenate([[0], np.cumsum(t)])
        level = m / t
        cmc = np.concatenate([[0], np.cumsum(level >= 0.5)])
        for g in genes:
            i, j = np.searchsorted(p, g["start"], "left"), np.searchsorted(p, g["end"], "right")
            n = j - i
            if n < min_cpgs:
                continue
            out[g["gene_id"]] = {"n_cpg": int(n),
                                 "meth": float((cm[j] - cm[i]) / (ct[j] - ct[i])),
                                 "frac_mcpg": float((cmc[j] - cmc[i]) / n)}
    return out


def compare(cpg: dict[str, tuple[float, float]], meth: dict[str, dict],
            threshold: float = METH_THRESHOLD) -> dict:
    """cpg: {gene_id: (cpg_oe, p_low)} for one feature."""
    genes = [g for g in cpg if g in meth]
    if len(genes) < 10:
        return {"n_genes": len(genes), "note": "too few genes"}
    oe = np.array([cpg[g][0] for g in genes])
    p_low = np.array([cpg[g][1] for g in genes])
    lvl = np.array([meth[g]["meth"] for g in genes])
    methylated = lvl >= threshold
    rho, p_rho = spearmanr(oe, lvl)
    res = {"n_genes": len(genes), "spearman_rho": float(rho), "spearman_p": float(p_rho),
           "frac_methylated": float(methylated.mean()),
           "median_meth": float(np.median(lvl))}
    if 0 < methylated.sum() < len(genes):
        u = mannwhitneyu(oe[~methylated], oe[methylated], alternative="greater").statistic
        res["auc_low_oe_predicts_methylated"] = float(u / (methylated.sum() * (~methylated).sum()))
    low = p_low > 0.5
    ok = np.isfinite(p_low)
    if ok.any():
        tp = int((low & methylated & ok).sum())
        fp = int((low & ~methylated & ok).sum())
        fn = int((~low & methylated & ok).sum())
        tn = int((~low & ~methylated & ok).sum())
        n = tp + fp + fn + tn
        res["low_class"] = {
            "tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "misclassification": (fp + fn) / n,
            "precision": tp / (tp + fp) if tp + fp else float("nan"),
            "recall": tp / (tp + fn) if tp + fn else float("nan"),
            "median_meth_low": float(np.median(lvl[low & ok])) if (low & ok).any() else None,
            "median_meth_high": float(np.median(lvl[~low & ok])) if (~low & ok).any() else None,
        }
    return res


def validate(accession: str, results_dir: Path, data_dir: Path, cov_path: Path, label: str,
             source: dict | None = None, min_cov: int = MIN_COV, min_cpgs: int = MIN_CPGS,
             threshold: float = METH_THRESHOLD, log=print) -> Path:
    from .annotate import read_features
    from .fetch import _md5

    outdir = Path(results_dir) / accession
    features = outdir / "features.tsv.gz"
    if not features.exists():
        raise FileNotFoundError(f"{features} missing; run `cpgme run {accession}` first")
    gffs = list((Path(data_dir) / accession).glob("*.gff"))
    if len(gffs) != 1:
        raise FileNotFoundError(f"expected one GFF in {Path(data_dir) / accession}, found {gffs}")

    log(f"{accession}: reading {cov_path}")
    cov = read_coverage(cov_path, min_cov)
    genes = protein_coding_genes(gffs[0])
    shared = set(cov) & set(genes)
    if not shared:
        raise ValueError(f"no chromosome names shared between {cov_path} and {gffs[0]}; "
                         "were reads aligned to the same assembly?")
    meth = gene_methylation(cov, genes, min_cpgs)
    log(f"{accession}: {len(meth)} genes with >= {min_cpgs} CpGs at >= {min_cov}x")

    vdir = outdir / "validation" / label
    vdir.mkdir(parents=True, exist_ok=True)
    with gzip.open(vdir / "gene_methylation.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_id", "n_cpg", "meth", "frac_mcpg"])
        for gid, r in sorted(meth.items()):
            w.writerow([gid, r["n_cpg"], f"{r['meth']:.4f}", f"{r['frac_mcpg']:.4f}"])

    cpg = read_features(features)
    stats = {f: compare(cpg[f], meth, threshold) for f in ("cds", "gene_body") if f in cpg}
    n_cpg_cov = int(sum(len(v[0]) for v in cov.values()))
    (vdir / "validation.json").write_text(json.dumps({
        "accession": accession, "label": label, "source": source or {},
        "coverage_file": str(cov_path), "coverage_md5": _md5(Path(cov_path)),
        "min_cov": min_cov, "min_cpgs": min_cpgs, "meth_threshold": threshold,
        "n_cpg_covered": n_cpg_cov, "n_genes_with_meth": len(meth),
        "features": stats,
        "date": dt.datetime.now(dt.timezone.utc).isoformat(),
    }, indent=2))
    for f, s in stats.items():
        if "spearman_rho" in s:
            lc = s.get("low_class", {})
            log(f"{accession} {label} {f}: n={s['n_genes']} rho={s['spearman_rho']:.3f} "
                f"AUC={s.get('auc_low_oe_predicts_methylated', float('nan')):.3f} "
                f"misclass={lc.get('misclassification', float('nan')):.3f}")
    return vdir
