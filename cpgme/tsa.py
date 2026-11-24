"""Transcriptome (TSA) tier: CDS from a Transcriptome Shotgun Assembly, no genome needed.

Used for species with no annotated genome. The TSA contigs are downloaded from NCBI's
GenBank TSA FTP area and a coding sequence is called on each contig (longest ORF), which
stands in for the "longest CDS per gene" of the genome tier.

Accessions are TSA project versions written like an NCBI master record of that version,
e.g. GCWQ01000000 (project GCWQ, version 01). Only the `cds` feature is produced: there
are no gene bodies or upstream regions without a genome.

ORF calling (TransDecoder-like, simplified):
  - all six frames; an ORF runs from the first ATG after a stop codon to the next stop;
  - an ORF may be 5'-partial (no stop upstream within the contig: it starts at the first
    complete codon of the frame) and/or 3'-partial (runs off the contig end), as
    TransDecoder allows;
  - the longest ORF per contig is kept if it is at least MIN_ORF_NT long;
  - contigs whose ORF protein is identical to one already kept are dropped, which removes
    exact duplicates and isoforms sharing a full ORF, but not partial overlaps.

Caveats, compared with the genome tier: contigs are not grouped into genes (isoforms and
fragments of one gene can each contribute a row), lowly expressed genes are missing, and
whole-organism libraries can include symbiont or gut transcripts.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

from .fetch import _download_to
from .seqio import read_fasta

TSA_FTP = "https://ftp.ncbi.nlm.nih.gov/genbank/tsa"
ACCESSION = re.compile(r"^([A-Z]{4})(\d{2})000000$")
MIN_ORF_NT = 300
STOPS = {"TAA", "TAG", "TGA"}
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")

CODON_TABLE = dict(zip(
    [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def is_tsa(accession: str) -> bool:
    return ACCESSION.match(accession) is not None


def _parts(accession: str) -> tuple[str, int]:
    m = ACCESSION.match(accession)
    if not m:
        raise ValueError(f"{accession} is not a TSA project accession like GCWQ01000000")
    return m.group(1), int(m.group(2))


def fetch_tsa(accession: str, outdir: Path) -> dict[str, Path]:
    """Download the TSA contigs (FASTA) and master record (GenBank flat file)."""
    prefix, version = _parts(accession)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    base = f"{TSA_FTP}/{prefix[0]}"
    files = {"TSA_FASTA": f"tsa.{prefix}.{version}.fsa_nt.gz",
             "TSA_MASTER": f"tsa.{prefix}.mstr.gbff.gz"}
    paths = {}
    for kind, name in files.items():
        dest = outdir / name
        if not dest.exists():
            _download_to(f"{base}/{name}", dest)
        paths[kind] = dest
    return paths


def master_info(master_gbff: Path) -> dict:
    """Organism, taxid and assembly method from the TSA master record."""
    with gzip.open(master_gbff, "rt") as fh:
        text = fh.read()
    org = re.search(r"^  ORGANISM\s+(.+)$", text, re.M)
    taxid = re.search(r'/db_xref="taxon:(\d+)"', text)
    method = re.search(r"Assembly Method\s+::\s+(.+)", text)
    name = re.search(r"Assembly Name\s+::\s+(.+)", text)
    return {
        "organism": org.group(1).strip() if org else "",
        "taxid": int(taxid.group(1)) if taxid else None,
        "assembly_method": method.group(1).strip() if method else None,
        "assembly_name": name.group(1).strip() if name else None,
    }


def revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def translate(cds: str) -> str:
    return "".join(CODON_TABLE.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))


def longest_orf(seq: str, min_nt: int = MIN_ORF_NT) -> dict | None:
    """Longest ORF over six frames, or None if none reaches min_nt (stop codon included)."""
    seq = seq.upper()
    best = None
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        n = len(s)
        for frame in range(3):
            start = None  # open ORF start in this frame
            after_stop = False  # has a stop been seen upstream in this frame?
            for i in range(frame, n - 2, 3):
                codon = s[i:i + 3]
                if codon in STOPS:
                    if start is not None:
                        best = _better(best, s, start, i + 3, strand, frame,
                                       five_partial=not after_stop and s[start:start + 3] != "ATG",
                                       three_partial=False)
                    start = None
                    after_stop = True
                elif start is None:
                    if codon == "ATG" or not after_stop:
                        start = i
            if start is not None:
                end = start + (n - start) // 3 * 3
                best = _better(best, s, start, end, strand, frame,
                               five_partial=not after_stop and s[start:start + 3] != "ATG",
                               three_partial=True)
    if best is None or len(best["seq"]) < min_nt:
        return None
    return best


def _better(best, s, start, end, strand, frame, five_partial, three_partial):
    if best is not None and end - start <= len(best["seq"]):
        return best
    return {"seq": s[start:end], "strand": strand, "frame": frame, "start": start,
            "five_partial": five_partial, "three_partial": three_partial}


def tsa_cds(tsa_fasta: Path, min_nt: int = MIN_ORF_NT) -> dict[str, dict]:
    """{contig_id: {"seq", "protein", "partial", "contig_length"}}, one ORF per contig,
    identical proteins kept once (first contig wins)."""
    out: dict[str, dict] = {}
    seen: set[str] = set()
    for header, seq in read_fasta(tsa_fasta):
        orf = longest_orf(seq, min_nt)
        if orf is None:
            continue
        protein = translate(orf["seq"]).rstrip("*")
        if protein in seen:
            continue
        seen.add(protein)
        partial = ("5'" if orf["five_partial"] else "") + ("3'" if orf["three_partial"] else "")
        out[header.split()[0]] = {"seq": orf["seq"], "protein": protein,
                                  "partial": partial, "contig_length": len(seq)}
    return out


def write_proteins(cds: dict[str, dict], out: Path) -> dict[str, str]:
    """Write one protein per contig as '>tsa:<contig> <contig>'; return {gene_id: contig}."""
    written = {}
    with gzip.open(out, "wt") as fh:
        for contig, rec in cds.items():
            gid = f"tsa:{contig}"
            fh.write(f">{gid} {contig}\n")
            p = rec["protein"]
            for i in range(0, len(p), 80):
                fh.write(p[i:i + 80] + "\n")
            written[gid] = contig
    return written

