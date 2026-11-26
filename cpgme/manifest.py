"""Phase 1: list candidate annotated genomes for each focal group from NCBI Datasets.

For every group taxon, all current assemblies with a gene annotation are pulled from the
Datasets v2 API and reduced to one assembly per species:

  1. RefSeq annotation over GenBank (submitter) annotation,
  2. then higher BUSCO completeness (RefSeq annotations report it; GenBank ones do not, and
     are scored here by scripts/busco_candidates.sh into manifest/busco.tsv),
  3. then the more recent annotation, then higher contig N50,

unless one of the species' assemblies is already in the manifest, which is then kept.

Each species gets a proposed `include` call: yes when BUSCO completeness is at least the
threshold, `check` when there is no BUSCO score, no when below the
threshold or the assembly is flagged atypical. Species already in the manifest are marked
so they are not added twice. The output is a review list, not the manifest itself.
"""

from __future__ import annotations

import csv
import json
import time
import urllib.parse
from pathlib import Path

from .fetch import API, _get

# group -> (NCBI taxon id, taxon name). Corals are the stony corals (Scleractinia);
# other anthozoans are listed under "anthozoa_other" for review only.
GROUPS = {
    "coral": (6125, "Scleractinia"),
    "anthozoa_other": (6101, "Anthozoa"),
    "bivalve": (6544, "Bivalvia"),
    "hemiptera": (7524, "Hemiptera"),
    "tunicate": (7712, "Tunicata"),
    "seastar": (7588, "Asteroidea"),
}
BUSCO_MIN = 0.80

COLUMNS = ["group", "species", "taxid", "accession", "source", "annotation_provider",
           "annotation_date", "busco_source", "busco_lineage", "busco_complete", "busco_fragmented",
           "protein_coding_genes", "assembly_level", "contig_n50", "genome_mb",
           "n_annotated_assemblies", "atypical", "in_manifest", "include", "reason"]


def taxon_reports(taxid: int) -> list[dict]:
    reports, token = [], None
    while True:
        q = {"filters.has_annotation": "true", "filters.assembly_version": "current",
             "page_size": 200}
        if token:
            q["page_token"] = token
        data = json.loads(_get(f"{API}/genome/taxon/{taxid}/dataset_report?"
                               + urllib.parse.urlencode(q)))
        reports += data.get("reports", [])
        token = data.get("next_page_token")
        if not token:
            return reports
        time.sleep(0.4)


def _row(r: dict) -> dict:
    ann = r.get("annotation_info", {})
    busco = ann.get("busco", {})
    st = r.get("assembly_stats", {})
    info = r.get("assembly_info", {})
    return {
        "species": r["organism"]["organism_name"],
        "taxid": r["organism"]["tax_id"],
        "accession": r["accession"],
        "source": "refseq" if r["accession"].startswith("GCF_") else "genbank",
        "annotation_provider": ann.get("provider", ""),
        "annotation_date": ann.get("release_date", ""),
        "busco_source": "ncbi" if busco else "",
        "busco_lineage": busco.get("busco_lineage", ""),
        "busco_complete": round(busco["complete"], 3) if "complete" in busco else "",
        "busco_fragmented": round(busco["fragmented"], 3) if "fragmented" in busco else "",
        "protein_coding_genes": ann.get("stats", {}).get("gene_counts", {})
                                   .get("protein_coding", ""),
        "assembly_level": info.get("assembly_level", ""),
        "contig_n50": st.get("contig_n50", ""),
        "genome_mb": round(int(st.get("total_sequence_length", 0)) / 1e6, 1),
        "atypical": "yes" if info.get("atypical", {}).get("is_atypical") else "",
    }


def _rank(row: dict) -> tuple:
    b = row["busco_complete"]
    return (row["source"] == "refseq", b if b != "" else -1, row["annotation_date"],
            row["contig_n50"] or 0)


def own_busco(path: Path | None) -> dict[str, dict]:
    """BUSCO scores computed here (scripts/busco_candidates.sh) for annotations NCBI did not
    score: accession -> row of manifest/busco.tsv."""
    if not path or not path.exists():
        return {}
    with open(path) as fh:
        return {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}


def candidates(groups: dict = GROUPS, manifest: Path | None = None,
               busco: Path | None = None) -> list[dict]:
    extra = own_busco(busco)
    have_acc, have_tax = set(), set()
    if manifest and manifest.exists():
        with open(manifest) as fh:
            for m in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                have_acc.add(m.get("accession", ""))
                have_tax.add(str(m.get("taxid", "")))
    seen_tax = set()
    out = []
    for group, (taxid, _) in groups.items():
        by_species: dict[int, list[dict]] = {}
        for r in taxon_reports(taxid):
            row = _row(r)
            by_species.setdefault(row["taxid"], []).append(row)
        for tid, rows in by_species.items():
            if tid in seen_tax:  # e.g. a scleractinian already listed under coral
                continue
            seen_tax.add(tid)
            # Keep an assembly already in the manifest (results and validation depend on it).
            kept = [r for r in rows if r["accession"] in have_acc]
            best = kept[0] if kept else max(rows, key=_rank)
            best["group"] = group
            best["n_annotated_assemblies"] = len(rows)
            best["in_manifest"] = ("same_accession" if best["accession"] in have_acc else
                                   "other_accession" if str(tid) in have_tax else "")
            if best["busco_complete"] == "" and best["accession"] in extra:
                e = extra[best["accession"]]
                best.update(busco_source="cpgme", busco_lineage=e["lineage"],
                            busco_complete=round(float(e["complete"]), 3),
                            busco_fragmented=round(float(e["fragmented"]), 3))
            b = best["busco_complete"]
            if group == "anthozoa_other":
                best["include"], best["reason"] = "no", "not a stony coral; review only"
            elif best["atypical"]:
                best["include"], best["reason"] = "no", "assembly flagged atypical"
            elif b == "":
                best["include"], best["reason"] = "check", "no BUSCO score yet"
            elif b < BUSCO_MIN:
                best["include"], best["reason"] = "no", f"BUSCO complete < {BUSCO_MIN:.0%}"
            else:
                best["include"], best["reason"] = "yes", ""
            out.append(best)
        time.sleep(0.4)
    out.sort(key=lambda r: (list(groups).index(r["group"]), r["species"]))
    return out


def write(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        fh.write(f"# Candidate annotated genomes from NCBI Datasets, one per species. "
                 f"include=yes: BUSCO complete >= {BUSCO_MIN:.0%}; check: no BUSCO score yet. "
                 f"busco_source=cpgme: scored here, protein mode, one protein per gene.\n")
        w = csv.DictWriter(fh, COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
