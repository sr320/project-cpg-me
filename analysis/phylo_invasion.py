"""Phylogenetic test: do invasive Hemiptera differ in CpG O/E profile once relatedness is modelled?

The dashboard's trait-space test shuffles invasion labels within clades, which removes the
big clade differences but treats species inside a clade as independent. Hemiptera are not:
13 of the invasive species are aphids or close relatives. Here the same question is asked
with D-PGLS (Adams 2014, Evolution 68:2675): a multivariate PGLS / phylogenetic PERMANOVA in
which the profile matrix and the invasion indicator are both transformed by the inverse
Cholesky factor of the Brownian-motion covariance C, and significance comes from
randomising residuals of the null model (RRPP; Collyer et al. 2015, Heredity 115:357).

Tree      Open Tree of Life synthetic tree, induced subtree for the 38 hemipteran species
          (api.opentreeoflife.org v3; cached in analysis/hemiptera_otl.tre). Topology only:
          branch lengths are Grafen's (rho = 1), with equal branch lengths as a sensitivity
          check. Polytomies stay polytomies (a star under BM).
Profiles  as in the dashboard: 'summary' (12 distribution statistics, median/MAD scaled,
          clipped at +-3), 'shape' (share of genes in 80 CpG O/E bins), 'goslim' (median
          shift per GO Slim term, missing = 0). Gene body and CDS.
Labels    manifest/traits.tsv. Contrasts: invasive vs native; invasive vs all others;
          invasive + introduced vs native.
Also      Kmult (Adams 2014, Syst Biol 63:685): phylogenetic signal of each profile.

Usage:  uv run python analysis/phylo_invasion.py [--perm 9999] [--refresh-tree]
Writes  analysis/phylo_invasion.tsv, analysis/phylo_invasion.png (with --extra report)
"""

from __future__ import annotations

import argparse, csv, json, re, urllib.request
from pathlib import Path

import numpy as np

TREE = Path("analysis/hemiptera_otl.tre")
CONTRASTS = {
    "invasive_vs_native": lambda s: {"invasive": 1, "native": 0}.get(s),
    "invasive_vs_others": lambda s: 1 if s == "invasive" else 0,
    "nonnative_vs_native": lambda s: {"invasive": 1, "introduced": 1, "native": 0}.get(s),
}


# ---------- tree ----------
def fetch_tree(names: list[str]) -> str:
    def post(path, body):
        r = urllib.request.Request("https://api.opentreeoflife.org/v3/" + path, data=json.dumps(body).encode(),
                                   headers={"content-type": "application/json"})
        return json.load(urllib.request.urlopen(r, timeout=120))
    m = post("tnrs/match_names", {"names": names, "context_name": "Insects"})
    ott = {r["name"]: r["matches"][0]["taxon"]["ott_id"] for r in m["results"] if r["matches"]}
    missing = set(names) - set(ott)
    if missing:
        raise SystemExit(f"Open Tree has no match for: {sorted(missing)}")
    return post("tree_of_life/induced_subtree", {"ott_ids": list(ott.values()), "label_format": "name"})["newick"]


def parse_newick(s: str) -> dict:
    """Nested dicts {'name', 'children'}; quoted labels allowed, branch lengths ignored."""
    pos = 0

    def label():
        nonlocal pos
        if pos < len(s) and s[pos] == "'":
            end = s.index("'", pos + 1)
            name, pos = s[pos + 1:end], end + 1
        else:
            m = re.match(r"[^,():;]*", s[pos:])
            name, pos = m.group(0), pos + m.end()
        if pos < len(s) and s[pos] == ":":
            m = re.match(r":[^,();]*", s[pos:])
            pos += m.end()
        return name

    def node():
        nonlocal pos
        kids = []
        if s[pos] == "(":
            pos += 1
            kids.append(node())
            while s[pos] == ",":
                pos += 1
                kids.append(node())
            pos += 1  # ')'
        return {"name": label(), "children": kids}

    return node()


def collapse(n: dict) -> dict:
    """Drop single-child (taxonomy) nodes."""
    while len(n["children"]) == 1:
        n = n["children"][0]
    return {"name": n["name"], "children": [collapse(c) for c in n["children"]]}


def tips(n):
    return [n["name"].replace("_", " ")] if not n["children"] else [t for c in n["children"] for t in tips(c)]


def bm_cov(root: dict, mode: str) -> tuple[list[str], np.ndarray]:
    """Brownian covariance C (tips x tips), total depth scaled to 1.
    grafen: node height = (tips below - 1), so branch lengths are differences in height.
    equal:  every branch length 1, tips then extended to make the tree ultrametric."""
    names = tips(root)
    idx = {t: i for i, t in enumerate(names)}
    n = len(names)
    C = np.zeros((n, n))
    if mode == "grafen":
        H = len(names) - 1
        def walk(nd, depth_above):
            # depth_above: distance from root to this node (root = 0), in Grafen units
            below = tips(nd)
            h = len(below) - 1
            d = H - h
            ii = [idx[t] for t in below]
            # shared path from root to this node is d: every pair below shares at least d
            C[np.ix_(ii, ii)] = d
            for c in nd["children"]:
                walk(c, d)
        walk(root, 0)
        np.fill_diagonal(C, H)
        return names, C / H
    # equal branch lengths: depth counted in edges, then ultrametricised at the tips
    def walk(nd, d):
        below = tips(nd)
        ii = [idx[t] for t in below]
        C[np.ix_(ii, ii)] = d
        for c in nd["children"]:
            walk(c, d + 1)
    walk(root, 0)
    D = C.max() + 1
    np.fill_diagonal(C, D)
    return names, C / D


# ---------- profiles ----------
def profiles(d: dict, feat: str, kind: str, keep: list[dict]) -> tuple[list[dict], np.ndarray]:
    rows, X = [], []
    if kind == "goslim":
        g = d["goslim"].get(feat)
        by = {}
        for ti, si, shift, *_ in g["cells"]:
            by.setdefault(si, np.zeros(len(g["terms"])))[ti] = shift
    for s in keep:
        fit = s["fits"].get(feat)
        if kind == "summary":
            if not fit:
                continue
            dbic = fit["delta_bic"]
            v = [fit["median"], fit["mean"], fit["low_mean"], fit["high_mean"], fit["low_sd"], fit["high_sd"],
                 fit["low_weight"], fit["frac_low_posterior"], fit["ashman_d"],
                 np.sign(dbic) * np.log10(1 + abs(dbic)), fit["gpc_oe_median"],
                 (s.get("goslim_test") or {}).get(feat, {}).get("rank_biserial")]
            if any(x is None for x in v):
                continue
        elif kind == "shape":
            h = s["hist"].get(feat)
            if not h:
                continue
            v = np.array(h, float) / sum(h)
        else:
            si = d["species"].index(s)
            if si not in by:
                continue
            v = by[si]
        rows.append(s)
        X.append(np.asarray(v, float))
    X = np.array(X)
    if kind == "summary":  # robust scaling, as in the dashboard
        med = np.median(X, 0)
        mad = 1.4826 * np.median(np.abs(X - med), 0)
        X = np.clip((X - med) / np.where(mad > 1e-12, mad, 1), -3, 3)
    return rows, X


# ---------- statistics ----------
def ss_res(X: np.ndarray, Y: np.ndarray) -> float:
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    R = Y - X @ beta
    return float((R * R).sum())


def dpgls(Y: np.ndarray, x: np.ndarray, C: np.ndarray, perm: int, rng) -> dict:
    """D-PGLS with RRPP. Returns F, R2, Z (effect size of log F vs its null distribution), p."""
    n = len(x)
    P = np.linalg.inv(np.linalg.cholesky(C))
    Yt, X0, X1 = P @ Y, P @ np.ones((n, 1)), P @ np.column_stack([np.ones(n), x])
    b0, *_ = np.linalg.lstsq(X0, Yt, rcond=None)
    fit0, res0 = X0 @ b0, Yt - X0 @ b0

    def F_of(Yp):
        s0, s1 = ss_res(X0, Yp), ss_res(X1, Yp)
        return ((s0 - s1) / 1) / (s1 / (n - 2)), (s0 - s1) / s0

    F, R2 = F_of(Yt)
    Fs = np.array([F_of(fit0 + res0[rng.permutation(n)])[0] for _ in range(perm)])
    lF, lFs = np.log(F), np.log(np.r_[F, Fs])
    return {"F": F, "R2": R2, "Z": (lF - lFs.mean()) / lFs.std(ddof=1),
            "p": (1 + (Fs >= F - 1e-12).sum()) / (perm + 1)}


def kmult(Y: np.ndarray, C: np.ndarray, perm: int, rng) -> tuple[float, float]:
    """Multivariate Blomberg's K (Adams 2014): 1 = Brownian expectation, < 1 weaker signal."""
    n = len(Y)
    Ci, one = np.linalg.inv(C), np.ones((n, 1))
    expected = (np.trace(C) - n / (one.T @ Ci @ one).item()) / (n - 1)

    def K(Yp):
        a = np.linalg.solve(one.T @ Ci @ one, one.T @ Ci @ Yp)
        R = Yp - one @ a
        return (np.trace(R.T @ R) / np.trace(R.T @ Ci @ R)) / expected

    k = K(Y)
    ks = np.array([K(Y[rng.permutation(n)]) for _ in range(perm)])
    return k, (1 + (ks >= k).sum()) / (perm + 1)


# ---------- figure ----------
def figure(root, names, d, traits, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sp = {s["organism"]: s for s in d["species"] if s["group"] == "hemiptera" and s["role"] == "focal"}
    bins = np.array(d["bins"])
    col = {"invasive": "#eb6834", "introduced": "#eda100", "native": "#2a78d6", "uncertain": "#898781"}
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 10), gridspec_kw={"width_ratios": [1.1, 1]}, sharey=True)
    y = {t: i for i, t in enumerate(names)}
    xs = {}

    def draw(nd, depth):
        if not nd["children"]:
            t = nd["name"].replace("_", " ")
            xs[t] = depth
            return y[t]
        ys = [draw(c, depth + 1) for c in nd["children"]]
        for c, yc in zip(nd["children"], ys):
            cx = depth + 1
            ax.plot([depth, depth], [min(ys), max(ys)], color="#52514e", lw=1)
            ax.plot([depth, cx], [yc, yc], color="#52514e", lw=1)
        return sum(ys) / len(ys)

    draw(root, 0)
    D = max(xs.values())
    for t, x in xs.items():
        ax.plot([x, D], [y[t], y[t]], color="#c3c2b7", lw=0.6, ls=":")
        st = traits.get(t, "uncertain")
        ax.text(D + 0.2, y[t], t, va="center", fontsize=8.5, style="italic", color=col[st])
    ax.set_xlim(-0.3, D + 9)
    ax.axis("off")
    H = np.full((len(names), 60), np.nan)
    for t, i in y.items():
        h = sp.get(t, {}).get("hist", {}).get("gene_body") or sp.get(t, {}).get("hist", {}).get("cds")
        if h:
            h = np.array(h[:60], float)
            H[i] = h / h.sum()
    ax2.imshow(H, aspect="auto", cmap="Greys", extent=[bins[0], bins[60], len(names) - 0.5, -0.5],
               interpolation="nearest")
    ax2.set_xlabel("CpG O/E (gene body; CDS for transcriptome species)")
    ax2.set_ylim(len(names) - 0.5, -0.5)
    ax2.set_yticks([])
    for st, c in col.items():
        ax.plot([], [], "s", color=c, label=st)
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    fig.suptitle("Hemiptera: Open Tree topology, invasion status (label colour) and CpG O/E distribution",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=150)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("docs/data/dashboard.json"))
    ap.add_argument("--traits", type=Path, default=Path("manifest/traits.tsv"))
    ap.add_argument("--out", type=Path, default=Path("analysis/phylo_invasion.tsv"))
    ap.add_argument("--perm", type=int, default=9999)
    ap.add_argument("--refresh-tree", action="store_true")
    a = ap.parse_args()
    d = json.loads(a.data.read_text())
    traits = {r["species"]: r["invasion"] for r in csv.DictReader(
        (l for l in open(a.traits) if not l.startswith("#")), delimiter="\t")}
    hemi = [s for s in d["species"] if s["group"] == "hemiptera" and s["role"] == "focal"]
    if a.refresh_tree or not TREE.exists():
        TREE.write_text(fetch_tree(sorted(s["organism"] for s in hemi)) + "\n")
    root = collapse(parse_newick(TREE.read_text().strip()))
    rng = np.random.default_rng(20261127)
    covs = {m: bm_cov(root, m) for m in ("grafen", "equal")}
    out = []
    for feat in ("gene_body", "cds"):
        for kind in ("summary", "shape", "goslim"):
            rows, Y = profiles(d, feat, kind, hemi)
            names = [s["organism"] for s in rows]
            for cname, f in CONTRASTS.items():
                lab = [f(traits.get(n)) for n in names]
                keep = [i for i, v in enumerate(lab) if v is not None]
                x = np.array([lab[i] for i in keep], float)
                Yk = Y[keep]
                for model in ("none", "grafen", "equal"):
                    if model == "none":
                        C = np.eye(len(keep))
                    else:
                        tn, Cf = covs[model]
                        ii = [tn.index(names[i]) for i in keep]
                        C = Cf[np.ix_(ii, ii)]
                    r = dpgls(Yk, x, C, a.perm, rng)
                    out.append({"feature": feat, "profile": kind, "contrast": cname, "phylogeny": model,
                                "n_1": int(x.sum()), "n_0": int(len(x) - x.sum()), **r, "K": None})
                    print(f"{feat:9} {kind:7} {cname:20} {model:6} n={int(x.sum())}/{int(len(x)-x.sum())} "
                          f"R2={r['R2']:.3f} Z={r['Z']:.2f} p={r['p']:.4f}", flush=True)
            for model in ("grafen", "equal"):
                tn, Cf = covs[model]
                ii = [tn.index(n) for n in names]
                k, p = kmult(Y, Cf[np.ix_(ii, ii)], a.perm, rng)
                out.append({"feature": feat, "profile": kind, "contrast": "Kmult_signal", "phylogeny": model,
                            "n_1": len(names), "n_0": 0, "F": None, "R2": None, "Z": None, "p": p, "K": k})
                print(f"{feat:9} {kind:7} Kmult {model}: K={k:.2f} p={p:.4f}")
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in out:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    try:
        figure(root, tips(root), d, traits, a.out.with_suffix(".png"))
    except ImportError:
        print("matplotlib not installed (uv sync --extra report): figure skipped")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
