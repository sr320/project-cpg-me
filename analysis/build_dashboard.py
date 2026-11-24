"""Build docs/data/dashboard.json for the GitHub Pages dashboard from results/.

Usage:  uv run python analysis/build_dashboard.py [--results results] [--out docs/data/dashboard.json]
        [--merge]   keep species already in --out that are not in --results (for adding a
                    few species on a machine without the full merged results)

Only summaries go into the JSON (fits, binned histograms, GO Slim statistics, validation
statistics), so the page stays small and no sequence data is published.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np

BINS = np.linspace(0, 2, 81)  # 0.025-wide bins over CpG O/E 0-2
FEATURES = ("cds", "gene_body", "upstream")
GROUP_ORDER = ["coral", "bivalve", "hemiptera", "tunicate", "seastar", "control"]
FIT_KEYS = ["n", "n_total", "median", "mean", "bimodal", "bimodal_bic", "delta_bic",
            "ashman_d", "low_mean", "low_sd", "low_weight", "high_mean", "high_sd",
            "high_weight", "threshold", "gpc_oe_median", "frac_low_posterior"]


def manifest(path: Path) -> dict[str, dict]:
    with open(path) as fh:
        rows = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
        return {r["accession"]: r for r in rows if r.get("accession")}


def histograms(features_tsv: Path) -> dict[str, list[int]]:
    vals = defaultdict(list)
    with gzip.open(features_tsv, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["pass_filter"] == "1" and r["cpg_oe"]:
                vals[r["feature"]].append(float(r["cpg_oe"]))
    out = {}
    for f, v in vals.items():
        counts, _ = np.histogram(np.clip(v, 0, 2 - 1e-9), BINS)
        out[f] = counts.astype(int).tolist()
    return out


def _r(x, nd=4):
    if x is None:
        return None
    if isinstance(x, bool):
        return x
    if isinstance(x, (int, np.integer)):
        return int(x)
    x = float(x)
    return None if not np.isfinite(x) else round(x, nd)


def _sig(x):
    """p- and q-values: keep 3 significant figures instead of fixed decimals."""
    if x is None:
        return None
    x = float(x)
    return None if not np.isfinite(x) else float(f"{x:.3g}")


def _val(k, v):
    return _sig(v) if k.startswith("p_") or k in ("p", "q", "spearman_p") else _r(v)


def species_record(d: Path, meta: dict) -> dict:
    s = json.loads((d / "summary.json").read_text())
    prov = json.loads((d / "provenance.json").read_text())
    rec = {
        "accession": d.name, "organism": s["organism"], "taxid": s.get("taxid"),
        "group": meta.get("group", ""), "role": meta.get("role", ""),
        "tier": meta.get("tier", ""), "notes": meta.get("notes", ""),
        "assembly": prov.get("assembly_name"), "annotation": prov.get("annotation"),
        "fits": {f: {k: _r(v.get(k)) for k in FIT_KEYS if k in v}
                 for f, v in s["features"].items()},
        "hist": histograms(d / "features.tsv.gz"),
    }
    gt = d / "goslim_test.json"
    if gt.exists():
        t = json.loads(gt.read_text())
        rec["annotation_stats"] = {k: t[k] for k in
                                   ("n_genes", "n_emapper_hits", "n_with_go", "n_with_goslim")}
        rec["goslim_test"] = {f: {k: _val(k, v) for k, v in r.items()} for f, r in t["tests"].items()}
    vroot = d / "validation"
    if vroot.exists():
        rec["validation"] = []
        for vj in sorted(vroot.glob("*/validation.json")):
            v = json.loads(vj.read_text())
            rec["validation"].append({
                "label": v["label"], "source": v.get("source", {}).get("description", ""),
                "n_cpg_covered": v.get("n_cpg_covered"),
                "n_genes_with_meth": v.get("n_genes_with_meth"),
                "meth_threshold": v.get("meth_threshold"), "min_cov": v.get("min_cov"),
                "features": {f: {k: (_val(k, x) if not isinstance(x, dict) else
                                     {kk: _r(xx) for kk, xx in x.items()})
                                 for k, x in r.items()} for f, r in v["features"].items()},
            })
    return rec


def goslim_terms(results: Path, species: list[dict]) -> dict:
    """{feature: {"terms": [...], "cells": [[term_idx, species_idx, shift, q, n, median]]}}"""
    out = {}
    for feature in ("cds", "gene_body"):
        rows = []
        for si, sp in enumerate(species):
            p = results / sp["accession"] / "goslim_cpg.tsv"
            if not p.exists():
                continue
            with open(p) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r["feature"] == feature:
                        rows.append((r["name"], r["class"], si, float(r["median_cpg_oe"]),
                                     float(r["median_other"]), float(r["q"]), int(r["n_genes"])))
        if not rows:
            continue
        shift = defaultdict(list)
        cls = {}
        for name, c, si, med, other, q, n in rows:
            shift[name].append(med - other)
            cls[name] = c
        terms = sorted(shift, key=lambda t: np.mean(shift[t]))
        idx = {t: i for i, t in enumerate(terms)}
        out[feature] = {
            "terms": [{"name": t, "class": cls[t]} for t in terms],
            "cells": [[idx[name], si, _r(med - other), _sig(q), n, _r(med)]
                      for name, c, si, med, other, q, n in rows],
        }
    return out


def _merge_goslim(new: dict, old: dict, species: list[dict]) -> dict:
    """Add GO Slim cells from a previous dashboard for species that have none locally."""
    if not old:
        return new
    pos = {s["accession"]: i for i, s in enumerate(species)}
    for feature, block in old["goslim"].items():
        have = {c[1] for c in new.get(feature, {}).get("cells", [])}
        tgt = new.setdefault(feature, {"terms": [], "cells": []})
        tidx = {t["name"]: i for i, t in enumerate(tgt["terms"])}
        for ti, si, *rest in block["cells"]:
            acc = old["species"][si]
            if acc not in pos or pos[acc] in have:
                continue
            term = block["terms"][ti]
            if term["name"] not in tidx:
                tidx[term["name"]] = len(tgt["terms"])
                tgt["terms"].append(term)
            tgt["cells"].append([tidx[term["name"]], pos[acc], *rest])
    # Rows are ordered by mean shift across species, as in goslim_terms().
    for block in new.values():
        shift = defaultdict(list)
        for ti, _, sh, *_ in block["cells"]:
            if sh is not None:
                shift[ti].append(sh)
        order = sorted(range(len(block["terms"])),
                       key=lambda i: np.mean(shift[i]) if shift[i] else np.inf)
        remap = {old_i: new_i for new_i, old_i in enumerate(order)}
        block["terms"] = [block["terms"][i] for i in order]
        block["cells"] = [[remap[c[0]], *c[1:]] for c in block["cells"]]
    return new


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=Path("results"))
    ap.add_argument("--manifest", type=Path, default=Path("manifest/species.tsv"))
    ap.add_argument("--out", type=Path, default=Path("docs/data/dashboard.json"))
    ap.add_argument("--merge", action="store_true",
                    help="keep species in the existing --out that are not in --results")
    a = ap.parse_args()
    meta = manifest(a.manifest)
    species = [species_record(d, meta.get(d.name, {})) for d in sorted(a.results.iterdir())
               if (d / "summary.json").exists()]
    old_goslim = {}
    if a.merge and a.out.exists():
        old = json.loads(a.out.read_text())
        local = {s["accession"] for s in species}
        kept = [s for s in old["species"] if s["accession"] not in local]
        # GO Slim cells index species by position; carry them over re-indexed.
        old_goslim = {"species": [s["accession"] for s in old["species"]],
                      "goslim": old.get("goslim", {})}
        for s in kept:
            m = meta.get(s["accession"], {})
            s.update({k: m.get(k, s.get(k, "")) for k in ("group", "role", "tier", "notes")})
        species += kept
    species.sort(key=lambda s: (GROUP_ORDER.index(s["group"]) if s["group"] in GROUP_ORDER
                                else 99, s["role"] != "focal", s["organism"]))
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                                text=True, check=True).stdout.strip()
    except Exception:
        commit = "unknown"
    data = {
        "generated": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "git_commit": commit,
        "bins": [_r(b, 3) for b in BINS.tolist()],
        "manifest_size": len(meta),
        "species": species,
        "goslim": _merge_goslim(goslim_terms(a.results, species), old_goslim, species),
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(data, separators=(",", ":")))
    n_val = sum(len(s.get("validation", [])) for s in species)
    print(f"{a.out}: {len(species)} species, {n_val} validation datasets, "
          f"{a.out.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
