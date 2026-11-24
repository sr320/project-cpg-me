"""Pilot report: CpG O/E distributions and mixture fits for every species in results/.

Usage:  uv run python analysis/pilot_report.py [--results results] [--out analysis/pilot]
Writes  <out>_summary.tsv and <out>_<feature>.png (one small-multiple panel per species).
        --merge  update only the rows of species in --results in an existing <out>_summary.tsv
                 and skip the figures (they need every species' per-gene table)
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LOW = "#eb6834"    # categorical slot 2 - low CpG O/E component
HIGH = "#2a78d6"   # categorical slot 1 - high CpG O/E component
HIST = "#c9c8c2"   # neutral fill for the observed distribution
INK = "#0b0b0b"
INK2 = "#52514e"
GROUP_ORDER = ["coral", "bivalve", "hemiptera", "tunicate", "seastar", "control"]


def load(results: Path, manifest: Path) -> list[dict]:
    meta = {}
    if manifest.exists():
        with open(manifest) as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                meta[r["accession"]] = r
    out = []
    for d in sorted(results.iterdir()):
        if not (d / "summary.json").exists():
            continue
        s = json.loads((d / "summary.json").read_text())
        s["group"] = meta.get(d.name, {}).get("group", "")
        s["role"] = meta.get(d.name, {}).get("role", "")
        s["table"] = pd.read_csv(d / "features.tsv.gz", sep="\t")
        out.append(s)
    out.sort(key=lambda s: (GROUP_ORDER.index(s["group"]) if s["group"] in GROUP_ORDER
                            else 99, s["organism"]))
    return out


def pdf(x, mu, sd, w):
    return w * np.exp(-0.5 * ((x - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))


def plot(species: list[dict], feature: str, path: Path) -> None:
    sp = [s for s in species if feature in s["features"]]
    if not sp:
        return
    ncol = min(4, len(sp))
    nrow = int(np.ceil(len(sp) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 2.7 * nrow), sharex=True,
                             squeeze=False)
    bins = np.linspace(0, 2.0, 101)
    grid = np.linspace(0, 2.0, 500)
    for ax, s in zip(axes.flat, sp):
        t = s["table"]
        x = t.loc[(t.feature == feature) & (t.pass_filter == 1), "cpg_oe"].dropna()
        f = s["features"][feature]
        ax.hist(x, bins=bins, density=True, color=HIST, edgecolor="white", linewidth=0.5)
        if "low_mean" in f:
            style = "-" if f["bimodal"] else ":"
            ax.plot(grid, pdf(grid, f["low_mean"], f["low_sd"], f["low_weight"]),
                    color=LOW, lw=2, ls=style)
            ax.plot(grid, pdf(grid, f["high_mean"], f["high_sd"], f["high_weight"]),
                    color=HIGH, lw=2, ls=style)
            if f["bimodal"] and f.get("threshold") is not None:
                ax.axvline(f["threshold"], color=INK2, lw=1, ls="--")
        ax.set_title(f"{s['organism']}", fontsize=9, style="italic", color=INK, loc="left")
        if f.get("bimodal"):
            verdict = (f"bimodal, D {f['ashman_d']:.2f}\n"
                       f"low {f['low_weight']:.0%} at {f['low_mean']:.2f}")
        elif "ashman_d" in f:
            verdict = (f"not bimodal, D {f['ashman_d']:.2f}\n"
                       f"median {f['median']:.2f}")
        else:
            verdict = "too few genes"
        ax.text(0.98, 0.95, f"{s['group']}\nn={f['n']:,}\n{verdict}", transform=ax.transAxes,
                ha="right", va="top", fontsize=7, color=INK2)
        ax.set_yticks([])
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(INK2)
        ax.tick_params(colors=INK2, labelsize=7)
    for ax in axes.flat[len(sp):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel(f"CpG O/E ({feature.replace('_', ' ')})", fontsize=8, color=INK2)
    handles = [plt.Line2D([], [], color=LOW, lw=2, label="low CpG O/E component"),
               plt.Line2D([], [], color=HIGH, lw=2, label="high CpG O/E component"),
               plt.Line2D([], [], color=INK2, lw=1, ls="--", label="threshold"),
               plt.Line2D([], [], color=INK2, lw=2, ls=":", label="not bimodal (ΔBIC ≤ 10 or D ≤ 2)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def summary_table(species: list[dict]) -> pd.DataFrame:
    rows = []
    for s in species:
        for feature, f in s["features"].items():
            rows.append({
                "group": s["group"], "organism": s["organism"], "accession": s["accession"],
                "feature": feature, "n": f.get("n"), "median": f.get("median"),
                "bimodal": f.get("bimodal"), "bimodal_bic": f.get("bimodal_bic"),
                "delta_bic": f.get("delta_bic"),
                "ashman_d": f.get("ashman_d"), "low_mean": f.get("low_mean"),
                "low_weight": f.get("low_weight"), "high_mean": f.get("high_mean"),
                "threshold": f.get("threshold"), "gpc_oe_median": f.get("gpc_oe_median"),
            })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=Path("results"))
    ap.add_argument("--manifest", type=Path, default=Path("manifest/species.tsv"))
    ap.add_argument("--out", type=Path, default=Path("analysis/pilot"))
    ap.add_argument("--merge", action="store_true",
                    help="update rows in an existing summary; do not redraw figures")
    a = ap.parse_args()
    species = load(a.results, a.manifest)
    tab = summary_table(species)
    path = Path(f"{a.out}_summary.tsv")
    if a.merge and path.exists():
        old = pd.read_csv(path, sep="\t")
        old = old[~old.accession.isin(tab.accession)]
        tab = pd.concat([old, tab], ignore_index=True)
        rank = {g: i for i, g in enumerate(GROUP_ORDER)}
        tab = tab.sort_values(["group", "organism", "accession"],
                              key=lambda c: c.map(rank).fillna(99) if c.name == "group" else c,
                              kind="stable")
    tab.to_csv(path, sep="\t", index=False, float_format="%.4g")
    if not a.merge:
        for feature in ("cds", "gene_body", "upstream"):
            plot(species, feature, Path(f"{a.out}_{feature}.png"))
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(tab.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
