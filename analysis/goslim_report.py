"""GO Slim report: does the housekeeping/inducible split of Dimond & Roberts 2015 hold?

Usage:  uv run python analysis/goslim_report.py [--results results] [--feature cds]
                                                 [--out analysis/goslim]
Reads   results/<acc>/goslim_cpg.tsv and goslim_test.json (from `cpgme annotate`).
Writes  <out>_tests.tsv   housekeeping vs inducible test, one row per species x feature
        <out>_<feature>.png  per-term shift in median CpG O/E, terms x species
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INK = "#0b0b0b"
INK2 = "#52514e"


def load(results: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    terms, tests = [], []
    for d in sorted(results.iterdir()):
        tj, tt = d / "goslim_test.json", d / "goslim_cpg.tsv"
        if not (tj.exists() and tt.exists()):
            continue
        organism = json.loads((d / "summary.json").read_text())["organism"]
        t = json.loads(tj.read_text())
        for feature, r in t["tests"].items():
            tests.append({"organism": organism, "accession": d.name, "feature": feature,
                          "n_genes": t["n_genes"], "n_with_goslim": t["n_with_goslim"], **r})
        df = pd.read_csv(tt, sep="\t")
        df["organism"] = organism
        terms.append(df)
    if not terms:
        raise SystemExit(f"no annotated species in {results}; run `cpgme annotate` first")
    return pd.concat(terms, ignore_index=True), pd.DataFrame(tests)


def plot(terms: pd.DataFrame, feature: str, path: Path) -> None:
    t = terms[terms.feature == feature].copy()
    if t.empty:
        return
    t["shift"] = t.median_cpg_oe - t.median_other
    wide = t.pivot_table(index="name", columns="organism", values="shift")
    order = wide.mean(axis=1).sort_values().index
    wide = wide.loc[order]
    cls = t.drop_duplicates("name").set_index("name")["class"].fillna("").reindex(order)
    sig = t.pivot_table(index="name", columns="organism", values="q").reindex(order)

    fig, ax = plt.subplots(figsize=(1.1 * wide.shape[1] + 3.2, 0.22 * len(order) + 1.2))
    lim = np.nanmax(np.abs(wide.values)) or 1
    xs, ys = np.meshgrid(np.arange(wide.shape[1]), np.arange(len(order)))
    v = wide.values
    ok = np.isfinite(v)
    sc = ax.scatter(xs[ok], ys[ok], c=v[ok], cmap="RdBu_r", vmin=-lim, vmax=lim,
                    s=np.where(sig.values[ok] < 0.05, 60, 18), edgecolors="none")
    ax.set_xticks(range(wide.shape[1]))
    ax.set_xticklabels(wide.columns, rotation=40, ha="right", fontsize=7, style="italic")
    ax.set_yticks(range(len(order)))
    tag = {"housekeeping": "  [HK]", "inducible": "  [IND]"}
    ax.set_yticklabels([f"{n}{tag.get(c, '')}" for n, c in zip(order, cls)], fontsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=INK2)
    cb = fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.02)
    cb.set_label("median CpG O/E: term − other annotated genes", fontsize=7, color=INK2)
    cb.ax.tick_params(labelsize=6)
    ax.set_title(f"GO Slim terms, {feature.replace('_', ' ')}\nlarge dot: q < 0.05; "
                 "[HK] housekeeping, [IND] inducible (a priori)", fontsize=8, color=INK,
                 loc="right")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=Path("results"))
    ap.add_argument("--out", type=Path, default=Path("analysis/goslim"))
    a = ap.parse_args()
    terms, tests = load(a.results)
    tests.to_csv(f"{a.out}_tests.tsv", sep="\t", index=False, float_format="%.4g")
    for feature in sorted(terms.feature.unique()):
        plot(terms, feature, Path(f"{a.out}_{feature}.png"))
    cols = ["organism", "feature", "n_housekeeping", "n_inducible", "median_housekeeping",
            "median_inducible", "rank_biserial", "p_housekeeping_lower", "supports_2015"]
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(tests[cols].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
