"""Per-accession pipeline: fetch -> features -> counts -> mixture fits -> outputs.

Output contract (results/<accession>/):
  features.tsv.gz   one row per gene x feature
  summary.json      mixture fits per feature type
  provenance.json   what was run, where, on what inputs
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import platform
import socket
import subprocess
from pathlib import Path

import numpy as np

from . import __version__
from .compute import count
from .features import genomic_features, longest_cds, protein_coding_genes
from .fetch import _md5, assembly_report, fetch
from .fit import fit_distribution

MIN_LENGTH = 300
MAX_FRAC_AMBIG = 0.05
COLUMNS = ["gene_id", "gene", "feature", "length", "n_ambig", "n_C", "n_G", "n_CG", "n_GC",
           "gc", "cpg_oe", "gpc_oe", "pass_filter", "p_low"]


def _git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                              text=True, cwd=Path(__file__).parent, check=True).stdout.strip()
    except Exception:
        return "unknown"


def _code_hash() -> str:
    """Hash of the modules that determine output values, so results are invalidated when
    the computation changes (but not when only downloading or the CLI changes)."""
    h = hashlib.sha256()
    for name in ("compute.py", "features.py", "fit.py", "run.py", "seqio.py"):
        h.update((Path(__file__).parent / name).read_bytes())
    return h.hexdigest()[:12]


def is_done(outdir: Path, genome: bool) -> bool:
    prov = outdir / "provenance.json"
    if not prov.exists():
        return False
    p = json.loads(prov.read_text())
    return p.get("code_hash") == _code_hash() and p.get("genome_features") == genome


def run(accession: str, results_dir: Path, data_dir: Path, genome: bool = False,
        force: bool = False, log=print) -> Path:
    outdir = Path(results_dir) / accession
    if not force and is_done(outdir, genome):
        log(f"{accession}: up to date, skipping")
        return outdir
    outdir.mkdir(parents=True, exist_ok=True)
    started = dt.datetime.now(dt.timezone.utc)

    log(f"{accession}: assembly report")
    report = assembly_report(accession)
    log(f"{accession}: fetching ({'with' if genome else 'without'} genome)")
    files = fetch(accession, Path(data_dir) / accession, genome=genome)

    rows: list[dict] = []
    log(f"{accession}: longest CDS per gene")
    cds = longest_cds(files["CDS_FASTA"])
    for gid, rec in cds.items():
        rows.append(_row(gid, rec["gene"], "cds", count(rec["seq"])))

    if genome:
        log(f"{accession}: gene bodies and upstream regions")
        genes = protein_coding_genes(files["GENOME_GFF"])
        names = {gid: rec["gene"] for gid, rec in cds.items()}
        for gid, feature, c in genomic_features(files["GENOME_FASTA"], genes):
            rows.append(_row(gid, names.get(gid, ""), feature, c))

    log(f"{accession}: fitting mixtures")
    summary = {"accession": accession,
               "organism": report["organism"]["organism_name"],
               "taxid": report["organism"]["tax_id"],
               "filters": {"min_length": MIN_LENGTH, "max_frac_ambig": MAX_FRAC_AMBIG},
               "features": {}}
    for feature in sorted({r["feature"] for r in rows}):
        idx = [i for i, r in enumerate(rows) if r["feature"] == feature]
        vals = np.array([rows[i]["cpg_oe"] if rows[i]["pass_filter"] else np.nan for i in idx])
        fit, p_low = fit_distribution(vals)
        gpc = np.array([rows[i]["gpc_oe"] for i in idx if rows[i]["pass_filter"]], dtype=float)
        fit["gpc_oe_median"] = float(np.nanmedian(gpc)) if gpc.size else None
        fit["n_total"] = len(idx)
        summary["features"][feature] = fit
        for i, p in zip(idx, p_low):
            rows[i]["p_low"] = p

    with gzip.open(outdir / "features.tsv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k)) for k in COLUMNS})
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

    import sklearn
    provenance = {
        "accession": accession,
        "organism": summary["organism"],
        "source_tier": "genome_annotation",
        "assembly_name": report.get("assembly_info", {}).get("assembly_name"),
        "annotation": report.get("annotation_info", {}).get("name"),
        "genome_features": genome,
        "cpgme_version": __version__,
        "code_hash": _code_hash(),
        "git_commit": _git_commit(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "sklearn": sklearn.__version__,
        "started": started.isoformat(),
        "finished": dt.datetime.now(dt.timezone.utc).isoformat(),
        "inputs_md5": {k: _md5(p) for k, p in files.items()},
    }
    (outdir / "provenance.json").write_text(json.dumps(provenance, indent=2))
    log(f"{accession}: done -> {outdir}")
    return outdir


def _row(gid: str, gene: str, feature: str, c) -> dict:
    return {
        "gene_id": gid, "gene": gene, "feature": feature,
        "length": c.length, "n_ambig": c.n_ambig, "n_C": c.n_C, "n_G": c.n_G,
        "n_CG": c.n_CG, "n_GC": c.n_GC, "gc": c.gc, "cpg_oe": c.cpg_oe, "gpc_oe": c.gpc_oe,
        "pass_filter": c.length >= MIN_LENGTH and c.frac_ambig <= MAX_FRAC_AMBIG
                       and np.isfinite(c.cpg_oe),
        "p_low": float("nan"),
    }


def _fmt(v):
    if isinstance(v, (bool, np.bool_)):
        return int(v)
    if isinstance(v, float):
        return "" if not np.isfinite(v) else f"{v:.6g}"
    return v
