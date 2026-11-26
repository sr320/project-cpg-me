"""Covariate check: does the housekeeping/inducible CpG O/E gap survive GC content and length?

Housekeeping and inducible genes can differ in GC content and length, and both shape CpG O/E
without any methylation. Three checks per species and feature (CDS, gene body):

  1. Residual test. Over all genes passing the filter, regress CpG O/E on GC, GC^2,
     log10 length, its square and GC x log10 length. Then compare the residuals of
     housekeeping and inducible genes (one-sided Mann-Whitney, rank-biserial), as in the
     original test. If the gap is GC or length, the residual gap disappears.
  2. Linear model on classified genes only: CpG O/E ~ housekeeping (0/1), without and with
     the same covariates. Reports how much of the raw difference is left.
  3. GpC O/E composition control: the same housekeeping/inducible test on GpC O/E, which
     reflects composition but not CpG depletion by methylation. A gap here would point to
     composition rather than methylation.

  4. Where WGBS validation exists, the direct version: are housekeeping genes more
     methylated than inducible ones in measured gene-body methylation, before and after the
     same GC and length adjustment? This separates real methylation from composition when
     the adjustment over-corrects (methylation itself lowers GC over time).

Expression breadth is not included: it needs RNA-seq across tissues for every species.

Usage:  uv run python analysis/covariate_check.py [--results results]
                                                  [--out analysis/covariate_check.tsv]
Reads   results/<acc>/features.tsv.gz and annotation.tsv.gz (from `cpgme run` and
        `cpgme annotate`), and validation/<label>/gene_methylation.tsv.gz if present.
Writes  <out> (checks 1-3) and <out stem>_wgbs.tsv (check 4).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cpgme.annotate import HOUSEKEEPING, INDUCIBLE  # noqa: E402

FEATURES = ("cds", "gene_body")


def gene_classes(annotation: Path) -> pd.Series:
    """gene_id -> 'hk' | 'ind' for genes in exactly one class (as in class_test)."""
    a = pd.read_csv(annotation, sep="\t", usecols=["gene_id", "goslim"], dtype=str)
    a = a.dropna(subset=["goslim"])
    terms = a.goslim.str.split(",").map(set)
    hk = terms.map(lambda t: bool(t & HOUSEKEEPING.keys()))
    ind = terms.map(lambda t: bool(t & INDUCIBLE.keys()))
    cls = pd.Series(np.where(hk & ~ind, "hk", np.where(ind & ~hk, "ind", "")), index=a.gene_id)
    return cls[cls != ""]


def design(df: pd.DataFrame) -> np.ndarray:
    gc = df.gc.to_numpy() - 0.4
    ll = np.log10(df.length.to_numpy()) - 3
    return np.column_stack([np.ones(len(df)), gc, gc ** 2, ll, ll ** 2, gc * ll])


def ols(y: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Coefficients and standard errors."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    cov = np.linalg.pinv(X.T @ X) * (resid @ resid) / dof
    return beta, np.sqrt(np.diag(cov))


def mwu_less(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """One-sided p (a lower) and rank-biserial (< 0: a ranks lower)."""
    r = stats.mannwhitneyu(a, b, alternative="less")
    return float(r.pvalue), float(2 * r.statistic / (len(a) * len(b)) - 1)


def check(df: pd.DataFrame, cls: pd.Series) -> dict:
    df = df.copy()
    X = design(df)
    beta, _ = ols(df.cpg_oe.to_numpy(), X)
    df["resid"] = df.cpg_oe.to_numpy() - X @ beta
    r2 = 1 - df.resid.var() / df.cpg_oe.var()

    df["cls"] = cls.reindex(df.gene_id).to_numpy()
    hk, ind = df[df.cls == "hk"], df[df.cls == "ind"]
    p_raw, rb_raw = mwu_less(hk.cpg_oe.to_numpy(), ind.cpg_oe.to_numpy())
    p_res, rb_res = mwu_less(hk.resid.to_numpy(), ind.resid.to_numpy())
    p_gpc, rb_gpc = mwu_less(hk.gpc_oe.to_numpy(), ind.gpc_oe.to_numpy())

    c = pd.concat([hk, ind])
    y, is_hk = c.cpg_oe.to_numpy(), (c.cls == "hk").to_numpy().astype(float)
    b0, _ = ols(y, np.column_stack([np.ones(len(c)), is_hk]))
    b1, se1 = ols(y, np.column_stack([design(c), is_hk]))
    t = b1[-1] / se1[-1]
    return {
        "n_housekeeping": len(hk), "n_inducible": len(ind),
        "median_gc_hk": hk.gc.median(), "median_gc_ind": ind.gc.median(),
        "median_len_hk": hk.length.median(), "median_len_ind": ind.length.median(),
        "covariate_r2": r2,
        "rb_raw": rb_raw, "p_raw": p_raw,
        "rb_residual": rb_res, "p_residual": p_res,
        "diff_raw": b0[-1], "diff_adjusted": b1[-1],
        "p_adjusted": float(stats.t.cdf(t, len(c) - len(b1))),  # one-sided: hk lower
        "frac_retained": b1[-1] / b0[-1] if b0[-1] else np.nan,
        "rb_gpc": rb_gpc, "p_gpc_hk_lower": p_gpc,
    }


def wgbs_check(gb: pd.DataFrame, meth: pd.DataFrame, cls: pd.Series) -> dict:
    """Housekeeping vs inducible in measured gene-body methylation (rank-biserial < 0:
    housekeeping more methylated), raw and adjusted for GC and length."""
    d = meth.merge(gb[["gene_id", "gc", "length"]], on="gene_id")
    X = design(d)
    beta, _ = ols(d.meth.to_numpy(), X)
    d["resid"] = d.meth.to_numpy() - X @ beta
    d["cls"] = cls.reindex(d.gene_id).to_numpy()
    hk, ind = d[d.cls == "hk"], d[d.cls == "ind"]
    p_raw, rb_raw = mwu_less(-hk.meth.to_numpy(), -ind.meth.to_numpy())
    p_adj, rb_adj = mwu_less(-hk.resid.to_numpy(), -ind.resid.to_numpy())
    return {"n_housekeeping": len(hk), "n_inducible": len(ind),
            "median_meth_hk": hk.meth.median(), "median_meth_ind": ind.meth.median(),
            "frac_methylated_hk": (hk.meth >= 0.5).mean(),
            "frac_methylated_ind": (ind.meth >= 0.5).mean(),
            "rb_meth": rb_raw, "p_meth": p_raw, "rb_meth_adjusted": rb_adj, "p_meth_adjusted": p_adj}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--results", type=Path, default=Path("results"))
    ap.add_argument("--out", type=Path, default=Path("analysis/covariate_check.tsv"))
    a = ap.parse_args()

    rows, wrows = [], []
    for d in sorted(a.results.iterdir()):
        ann, feat = d / "annotation.tsv.gz", d / "features.tsv.gz"
        if not (ann.exists() and feat.exists()):
            continue
        organism = json.loads((d / "summary.json").read_text())["organism"]
        cls = gene_classes(ann)
        f = pd.read_csv(feat, sep="\t", usecols=["gene_id", "feature", "length", "gc",
                                                 "cpg_oe", "gpc_oe", "pass_filter"])
        f = f[(f.pass_filter == 1) & f.cpg_oe.notna() & f.gpc_oe.notna()]
        for feature in FEATURES:
            sub = f[f.feature == feature]
            if sub.empty:
                continue
            rows.append({"organism": organism, "accession": d.name, "feature": feature,
                         **check(sub, cls)})
        for gm in sorted((d / "validation").glob("*/gene_methylation.tsv.gz")):
            wrows.append({"organism": organism, "accession": d.name, "dataset": gm.parent.name,
                          **wgbs_check(f[f.feature == "gene_body"],
                                       pd.read_csv(gm, sep="\t"), cls)})
    if not rows:
        raise SystemExit(f"no annotated species in {a.results}")
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False, float_format="%.4g")

    show = out[["organism", "feature", "n_housekeeping", "n_inducible", "median_gc_hk",
                "median_gc_ind", "covariate_r2", "rb_raw", "rb_residual", "p_residual",
                "frac_retained", "rb_gpc"]]
    with pd.option_context("display.width", 200, "display.max_columns", 20,
                           "display.float_format", "{:.3g}".format):
        print(show.to_string(index=False))
    print(f"\nwrote {a.out}")
    if wrows:
        w = pd.DataFrame(wrows)
        wout = a.out.with_name(a.out.stem + "_wgbs.tsv")
        w.to_csv(wout, sep="\t", index=False, float_format="%.4g")
        with pd.option_context("display.width", 200, "display.max_columns", 20,
                               "display.float_format", "{:.3g}".format):
            print(w.drop(columns=["accession", "p_meth", "p_meth_adjusted"]).to_string(index=False))
        print(f"\nwrote {wout}")


if __name__ == "__main__":
    main()
