"""Functional annotation and the GO Slim test of Dimond & Roberts 2015.

Steps, per accession (after `cpgme run`):

1. Write one protein per gene: the protein of the longest CDS, so annotation rows join to
   `features.tsv.gz` on gene_id (NCBI GeneID).
2. Annotate with eggNOG-mapper (diamond mode). This needs the eggNOG database (~50 GB) and
   is meant for Raven; an existing `*.emapper.annotations` file can be passed instead.
3. Map each gene's GO terms to GO Slim biological-process terms by walking up the GO graph
   (is_a and part_of).
4. Per GO Slim term and feature type: gene count, median CpG O/E, share of genes in the
   low-CpG O/E class, and a Mann-Whitney test against all other annotated genes.
5. The a priori test: genes in housekeeping terms should have lower CpG O/E (more
   historical germline methylation) than genes in inducible terms.

Outputs (results/<accession>/):
  proteins.faa.gz             one protein per gene, header = gene_id
  annotation.tsv.gz           gene_id, eggNOG fields, GO terms, GO Slim terms
  goslim_cpg.tsv              per GO Slim term x feature statistics
  goslim_test.json            housekeeping vs inducible test per feature
  annotation_provenance.json  inputs, ontology versions, code hash
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import platform
import shutil
import socket
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

from . import __version__
from .features import longest_cds
from .seqio import read_fasta

GO_BASIC_URL = "https://current.geneontology.org/ontology/go-basic.obo"
GO_SLIM_URL = "https://current.geneontology.org/ontology/subsets/goslim_generic.obo"

# goslim_generic has no general stress-response, translation or signal-transduction term;
# the 2015 paper's categories need them, so they are added to the slim.
EXTRA_SLIM_TERMS = {
    "GO:0006950": "response to stress",
    "GO:0006412": "translation",
    "GO:0007165": "signal transduction",
}

# A priori classes for the 2015 test. Housekeeping: constitutive cellular machinery,
# expected low CpG O/E. Inducible: environment-responsive and signalling functions,
# expected high CpG O/E. Chosen before looking at any annotation results. DNA repair is
# left out: in GO it is a descendant of "response to stress", so it would sit in both.
HOUSEKEEPING = {
    "GO:0006260": "DNA replication",
    "GO:0016071": "mRNA metabolic process",
    "GO:0006399": "tRNA metabolic process",
    "GO:0042254": "ribosome biogenesis",
    "GO:0006412": "translation",
    "GO:0002181": "cytoplasmic translation",
    "GO:0006457": "protein folding",
    "GO:0000278": "mitotic cell cycle",
    "GO:0007059": "chromosome segregation",
    "GO:0140053": "mitochondrial gene expression",
    "GO:0006913": "nucleocytoplasmic transport",
    "GO:0006886": "intracellular protein transport",
}
INDUCIBLE = {
    "GO:0007155": "cell adhesion",
    "GO:0023052": "signaling",
    "GO:0007165": "signal transduction",
    "GO:0002376": "immune system process",
    "GO:0098542": "defense response to other organism",
    "GO:0006954": "inflammatory response",
    "GO:0006950": "response to stress",
    "GO:0012501": "programmed cell death",
    "GO:0048856": "anatomical structure development",
    "GO:0030154": "cell differentiation",
    "GO:0030198": "extracellular matrix organization",
}

ANNOT_COLUMNS = ["gene_id", "protein_id", "seed_ortholog", "evalue", "eggnog_ogs",
                 "cog_category", "preferred_name", "description", "kegg_ko", "go",
                 "goslim", "goslim_names"]
FEATURES = ("cds", "gene_body")
MIN_TERM_GENES = 20


# ---- ontology ------------------------------------------------------------------------

class GO:
    """GO graph from an OBO file: parents via is_a and part_of, alt_id and obsolete
    handling, and ancestor sets."""

    def __init__(self, obo: Path):
        self.name: dict[str, str] = {}
        self.namespace: dict[str, str] = {}
        self.parents: dict[str, set[str]] = defaultdict(set)
        self.alt: dict[str, str] = {}
        self.obsolete: set[str] = set()
        self.version = ""
        self._anc: dict[str, frozenset[str]] = {}
        term = None
        with open(obo) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("data-version:"):
                    self.version = line.split(":", 1)[1].strip()
                if line.startswith("["):
                    term = {} if line == "[Term]" else None
                    continue
                if term is None or ":" not in line:
                    continue
                key, _, val = line.partition(": ")
                if key == "id":
                    term = {"id": val}
                elif key == "name":
                    self.name[term["id"]] = val
                elif key == "namespace":
                    self.namespace[term["id"]] = val
                elif key == "alt_id":
                    self.alt[val] = term["id"]
                elif key == "is_obsolete" and val == "true":
                    self.obsolete.add(term["id"])
                elif key == "is_a":
                    self.parents[term["id"]].add(val.split()[0])
                elif key == "relationship" and val.startswith("part_of "):
                    self.parents[term["id"]].add(val.split()[1])

    def canonical(self, go_id: str) -> str | None:
        go_id = self.alt.get(go_id, go_id)
        if go_id not in self.name or go_id in self.obsolete:
            return None
        return go_id

    def ancestors(self, go_id: str) -> frozenset[str]:
        """The term and all its ancestors."""
        if go_id in self._anc:
            return self._anc[go_id]
        out, stack = {go_id}, [go_id]
        while stack:
            for p in self.parents.get(stack.pop(), ()):
                if p not in out:
                    out.add(p)
                    stack.append(p)
        self._anc[go_id] = frozenset(out)
        return self._anc[go_id]


def slim_terms(slim_obo: Path, go: GO) -> dict[str, str]:
    """Biological-process slim terms {id: name}, plus EXTRA_SLIM_TERMS."""
    ids = [t for t in GO(slim_obo).name if go.namespace.get(t) == "biological_process"]
    terms = {t: go.name[t] for t in ids if go.canonical(t) == t}
    terms.update({t: go.name.get(t, n) for t, n in EXTRA_SLIM_TERMS.items()})
    return terms


def map_to_slim(go_ids: list[str], go: GO, slim: dict[str, str]) -> set[str]:
    hits: set[str] = set()
    for g in go_ids:
        c = go.canonical(g)
        if c is not None:
            hits |= go.ancestors(c) & slim.keys()
    return hits


def ensure_ontology(cache: Path) -> tuple[Path, Path]:
    from .fetch import _download_to
    cache.mkdir(parents=True, exist_ok=True)
    paths = []
    for url in (GO_BASIC_URL, GO_SLIM_URL):
        p = cache / url.rsplit("/", 1)[1]
        if not p.exists():
            _download_to(url, p)
        paths.append(p)
    return paths[0], paths[1]


# ---- proteins and eggNOG -------------------------------------------------------------

def write_gene_proteins(cds_fasta: Path, prot_fasta: Path, out: Path,
                        gff: Path | None = None) -> dict[str, str]:
    """Write the longest-CDS protein for each gene with header '>gene_id protein_id'.
    Returns {gene_id: protein_id} for the genes written."""
    want = {rec["protein_id"]: gid for gid, rec in longest_cds(cds_fasta, gff).items()
            if rec["protein_id"]}
    written: dict[str, str] = {}
    with gzip.open(out, "wt") as fh:
        for header, seq in read_fasta(prot_fasta):
            pid = header.split()[0]
            gid = want.get(pid)
            if gid is None or gid in written:
                continue
            fh.write(f">{gid} {pid}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
            written[gid] = pid
    return written


def run_emapper(proteins: Path, outdir: Path, data_dir: Path | None, cpu: int) -> Path:
    exe = shutil.which("emapper.py")
    if exe is None:
        raise RuntimeError(
            "emapper.py not found. Install eggNOG-mapper and its database (Raven), or pass "
            "--emapper-annotations with an existing *.emapper.annotations file.")
    outdir.mkdir(parents=True, exist_ok=True)
    plain = outdir / "proteins.faa"
    with gzip.open(proteins, "rt") as src, open(plain, "w") as dst:
        shutil.copyfileobj(src, dst)
    cmd = [exe, "-i", str(plain), "--itype", "proteins", "-m", "diamond",
           "--output", "emapper", "--output_dir", str(outdir), "--cpu", str(cpu),
           "--override"]
    if data_dir:
        cmd += ["--data_dir", str(data_dir)]
    subprocess.run(cmd, check=True)
    plain.unlink()
    return outdir / "emapper.emapper.annotations"


def parse_emapper(path: Path) -> dict[str, dict]:
    """Read an eggNOG-mapper v2 .annotations file -> {query: fields}. '-' becomes ''."""
    out: dict[str, dict] = {}
    header = None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#query"):
                header = line[1:].rstrip("\n").split("\t")
                continue
            if line.startswith("#") or not line.strip():
                continue
            if header is None:
                raise ValueError(f"{path}: no '#query' header line; not an emapper v2 file?")
            rec = dict(zip(header, line.rstrip("\n").split("\t")))
            rec = {k: ("" if v == "-" else v) for k, v in rec.items()}
            out[rec["query"]] = rec
    return out


# ---- statistics ----------------------------------------------------------------------

def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN passes through)."""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty_like(pv)
    adj[order] = np.minimum(ranked, 1)
    q[ok] = adj
    return q


def _mwu(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Two-sided Mann-Whitney p and rank-biserial effect (positive: a > b)."""
    if len(a) == 0 or len(b) == 0:
        return float("nan"), float("nan")
    res = mannwhitneyu(a, b, alternative="two-sided")
    return float(res.pvalue), float(2 * res.statistic / (len(a) * len(b)) - 1)


def read_features(path: Path) -> dict[str, dict[str, tuple[float, float]]]:
    """{feature: {gene_id: (cpg_oe, p_low)}} for rows passing the filter."""
    out: dict[str, dict[str, tuple[float, float]]] = defaultdict(dict)
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["pass_filter"] != "1" or not r["cpg_oe"]:
                continue
            p_low = float(r["p_low"]) if r["p_low"] else float("nan")
            out[r["feature"]][r["gene_id"]] = (float(r["cpg_oe"]), p_low)
    return out


def goslim_stats(cpg: dict[str, tuple[float, float]], gene_slim: dict[str, set[str]],
                 slim: dict[str, str], min_genes: int = MIN_TERM_GENES) -> list[dict]:
    """Per GO Slim term: compare its genes with all other annotated genes."""
    annotated = [g for g in cpg if gene_slim.get(g)]
    rows = []
    for term, name in sorted(slim.items()):
        inside = [g for g in annotated if term in gene_slim[g]]
        if len(inside) < min_genes:
            continue
        a = np.array([cpg[g][0] for g in inside])
        b = np.array([cpg[g][0] for g in annotated if term not in gene_slim[g]])
        pl = np.array([cpg[g][1] for g in inside])
        p, r = _mwu(a, b)
        rows.append({
            "goslim": term, "name": name, "n_genes": len(inside),
            "median_cpg_oe": float(np.median(a)), "mean_cpg_oe": float(a.mean()),
            "frac_low": float(np.nanmean(pl > 0.5)) if np.isfinite(pl).any() else float("nan"),
            "median_other": float(np.median(b)) if len(b) else float("nan"),
            "rank_biserial": r, "p": p,
            "class": ("housekeeping" if term in HOUSEKEEPING
                      else "inducible" if term in INDUCIBLE else ""),
        })
    for row, q in zip(rows, bh(np.array([r["p"] for r in rows]))):
        row["q"] = float(q)
    return rows


def class_test(cpg: dict[str, tuple[float, float]], gene_slim: dict[str, set[str]]) -> dict:
    """Housekeeping vs inducible genes. Genes in both classes are excluded."""
    hk, ind, n_both = [], [], 0
    for g, terms in gene_slim.items():
        if g not in cpg:
            continue
        in_hk, in_ind = bool(terms & HOUSEKEEPING.keys()), bool(terms & INDUCIBLE.keys())
        n_both += in_hk and in_ind
        if in_hk and not in_ind:
            hk.append(cpg[g])
        elif in_ind and not in_hk:
            ind.append(cpg[g])
    hk_a, ind_a = np.array(hk).reshape(-1, 2), np.array(ind).reshape(-1, 2)
    if len(hk_a) and len(ind_a):
        res = mannwhitneyu(hk_a[:, 0], ind_a[:, 0], alternative="less")
        p_less = float(res.pvalue)
        rb = float(2 * res.statistic / (len(hk_a) * len(ind_a)) - 1)
    else:
        p_less = rb = float("nan")
    med = lambda a: float(np.median(a[:, 0])) if len(a) else float("nan")  # noqa: E731
    low = lambda a: float(np.nanmean(a[:, 1] > 0.5)) if len(a) else float("nan")  # noqa: E731
    return {
        "n_housekeeping": len(hk_a), "n_inducible": len(ind_a), "n_both_excluded": n_both,
        "median_housekeeping": med(hk_a), "median_inducible": med(ind_a),
        "frac_low_housekeeping": low(hk_a), "frac_low_inducible": low(ind_a),
        # rank-biserial < 0: housekeeping genes rank lower, as the 2015 paper predicts
        "rank_biserial": rb,
        "p_housekeeping_lower": p_less,
        "supports_2015": bool(np.isfinite(p_less) and p_less < 0.05 and rb < 0),
    }


# ---- driver --------------------------------------------------------------------------

def _code_hash() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:12]


def annotate(accession: str, results_dir: Path, data_dir: Path,
             emapper_annotations: Path | None = None, emapper_data: Path | None = None,
             cpu: int = 4, proteins_only: bool = False, log=print) -> Path:
    from .fetch import _md5, fetch
    from .tsa import fetch_tsa, is_tsa, tsa_cds, write_proteins

    outdir = Path(results_dir) / accession
    features = outdir / "features.tsv.gz"
    if not features.exists():
        raise FileNotFoundError(f"{features} missing; run `cpgme run {accession}` first")
    started = dt.datetime.now(dt.timezone.utc)

    log(f"{accession}: one protein per gene")
    proteins = outdir / "proteins.faa.gz"
    if is_tsa(accession):
        files = fetch_tsa(accession, Path(data_dir) / accession)
        files["PROT_FASTA"] = files["TSA_FASTA"]  # proteins are translated from the contigs
        gene_protein = write_proteins(tsa_cds(files["TSA_FASTA"]), proteins)
    else:
        files = fetch(accession, Path(data_dir) / accession, genome=False)
        if "PROT_FASTA" not in files:
            raise ValueError(f"{accession}: NCBI package has no protein FASTA")
        gene_protein = write_gene_proteins(files["CDS_FASTA"], files["PROT_FASTA"], proteins,
                                           files.get("GENOME_GFF"))
    log(f"{accession}: {len(gene_protein)} proteins -> {proteins}")
    if proteins_only:
        return proteins

    if emapper_annotations is None:
        log(f"{accession}: running eggNOG-mapper")
        emapper_annotations = run_emapper(proteins, outdir / "emapper", emapper_data, cpu)
    emap = parse_emapper(Path(emapper_annotations))
    unknown = set(emap) - set(gene_protein)
    if emap and len(unknown) > 0.5 * len(emap):
        raise ValueError(f"{emapper_annotations}: most queries are not gene_ids of "
                         f"{accession}; was it run on {proteins.name}?")

    log(f"{accession}: mapping GO -> GO Slim")
    go_obo, slim_obo = ensure_ontology(Path(data_dir) / "_ontology")
    go = GO(go_obo)
    slim = slim_terms(slim_obo, go)

    gene_slim: dict[str, set[str]] = {}
    with gzip.open(outdir / "annotation.tsv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, ANNOT_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for gid, pid in gene_protein.items():
            rec = emap.get(gid, {})
            gos = [g for g in rec.get("GOs", "").split(",") if g]
            terms = map_to_slim(gos, go, slim)
            gene_slim[gid] = terms
            w.writerow({
                "gene_id": gid, "protein_id": pid,
                "seed_ortholog": rec.get("seed_ortholog", ""), "evalue": rec.get("evalue", ""),
                "eggnog_ogs": rec.get("eggNOG_OGs", ""),
                "cog_category": rec.get("COG_category", ""),
                "preferred_name": rec.get("Preferred_name", ""),
                "description": rec.get("Description", ""),
                "kegg_ko": rec.get("KEGG_ko", ""), "go": ",".join(gos),
                "goslim": ",".join(sorted(terms)),
                "goslim_names": ";".join(slim[t] for t in sorted(terms)),
            })

    log(f"{accession}: GO Slim statistics")
    cpg = read_features(features)
    stat_cols = ["feature", "goslim", "name", "class", "n_genes", "median_cpg_oe",
                 "mean_cpg_oe", "frac_low", "median_other", "rank_biserial", "p", "q"]
    tests = {}
    with open(outdir / "goslim_cpg.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, stat_cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for feature in FEATURES:
            if feature not in cpg:
                continue
            for row in goslim_stats(cpg[feature], gene_slim, slim):
                w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v)
                            for k, v in {"feature": feature, **row}.items()})
            tests[feature] = class_test(cpg[feature], gene_slim)

    n_go = sum(bool(emap.get(g, {}).get("GOs")) for g in gene_protein)
    (outdir / "goslim_test.json").write_text(json.dumps({
        "accession": accession,
        "n_genes": len(gene_protein),
        "n_emapper_hits": sum(g in emap for g in gene_protein),
        "n_with_go": n_go,
        "n_with_goslim": sum(bool(t) for t in gene_slim.values()),
        "housekeeping_terms": HOUSEKEEPING,
        "inducible_terms": INDUCIBLE,
        "tests": tests,
    }, indent=2))

    (outdir / "annotation_provenance.json").write_text(json.dumps({
        "accession": accession,
        "cpgme_version": __version__,
        "annotate_code_hash": _code_hash(),
        "emapper_annotations": str(emapper_annotations),
        "emapper_annotations_md5": _md5(Path(emapper_annotations)),
        "emapper_header": _emapper_header(Path(emapper_annotations)),
        "go_version": go.version,
        "goslim_version": GO(slim_obo).version,
        "extra_slim_terms": EXTRA_SLIM_TERMS,
        "protein_fasta_md5": _md5(files["PROT_FASTA"]),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "started": started.isoformat(),
        "finished": dt.datetime.now(dt.timezone.utc).isoformat(),
    }, indent=2))
    for f, t in tests.items():
        log(f"{accession} {f}: housekeeping median {t['median_housekeeping']:.3f} "
            f"(n={t['n_housekeeping']}) vs inducible {t['median_inducible']:.3f} "
            f"(n={t['n_inducible']}); p(hk lower)={t['p_housekeeping_lower']:.2g}")
    log(f"{accession}: annotation done -> {outdir}")
    return outdir


def _emapper_header(path: Path) -> list[str]:
    """The '##' comment lines eggNOG-mapper writes (version, command, database)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    lines = []
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            lines.append(line[2:].strip())
    return lines
