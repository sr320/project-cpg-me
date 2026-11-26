"""Build per-gene feature sequences and count them.

Features
--------
cds        longest CDS isoform per gene (from NCBI cds_from_genomic.fna). This is the
           cross-species comparable unit and needs no genome download.
gene_body  full protein-coding gene span from the GFF, introns included (needs genome).
upstream   1 kb upstream of the gene on its own strand, clipped at sequence ends
           (needs genome).

Gene IDs are NCBI GeneIDs ("GeneID:12345"), which link CDS records to GFF genes. Submitter
(GenBank) annotations often have no GeneID; they fall back to the locus tag
("locus_tag:ABC_00001"), which the GFF genes carry. Some GenBank CDS FASTA headers carry
neither a GeneID nor a locus tag; their gene is then looked up by protein_id in the GFF.
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Iterator

from .compute import Counts, count
from .seqio import read_fasta, read_gff

TAG = re.compile(r"\[(\w+)=([^\]]*)\]")
UPSTREAM_BP = 1000


def _gene_id_from_cds_header(tags: dict[str, str]) -> str | None:
    for xref in tags.get("db_xref", "").split(","):
        if xref.startswith("GeneID:"):
            return xref
    if tags.get("locus_tag"):
        return f"locus_tag:{tags['locus_tag']}"
    if "gene" in tags:
        return f"gene:{tags['gene']}"
    return None


def protein_gene_ids(gff: Path) -> dict[str, str]:
    """{protein_id: gene_id} from the GFF CDS rows (GeneID, else locus tag), for CDS FASTA
    headers that carry no gene identifier."""
    out: dict[str, str] = {}
    for rec in read_gff(gff, {"CDS"}):
        a = rec["attrs"]
        pid = a.get("protein_id")
        if not pid or pid in out:
            continue
        gid = next((x for x in a.get("Dbxref", "").split(",") if x.startswith("GeneID:")), None)
        if gid is None and a.get("locus_tag"):
            gid = f"locus_tag:{a['locus_tag']}"
        if gid:
            out[pid] = gid
    return out


def longest_cds(cds_fasta: Path, gff: Path | None = None) -> dict[str, dict]:
    """Return {gene_id: {"seq", "transcript", "protein_id", "gene", "n_isoforms"}}.
    With `gff`, CDS records whose header has no gene identifier are assigned to a gene via
    their protein_id."""
    best: dict[str, dict] = {}
    n_iso: dict[str, int] = defaultdict(int)
    by_protein: dict[str, str] | None = None
    for header, seq in read_fasta(cds_fasta):
        tags = dict(TAG.findall(header))
        if tags.get("pseudo") == "true":
            continue
        gid = _gene_id_from_cds_header(tags)
        if gid is None and gff is not None and tags.get("protein_id"):
            if by_protein is None:
                by_protein = protein_gene_ids(gff)
            gid = by_protein.get(tags["protein_id"])
        if gid is None:
            continue
        n_iso[gid] += 1
        if gid not in best or len(seq) > len(best[gid]["seq"]):
            best[gid] = {
                "seq": seq,
                "transcript": header.split()[0],
                "protein_id": tags.get("protein_id", ""),
                "gene": tags.get("gene", ""),
                "partial": tags.get("partial", ""),
            }
    for gid, rec in best.items():
        rec["n_isoforms"] = n_iso[gid]
    return best


def protein_coding_genes(gff: Path) -> dict[str, list[dict]]:
    """Return {seqid: [gene records]} for protein-coding genes with a GeneID (or, failing
    that, a locus tag)."""
    by_seq: dict[str, list[dict]] = defaultdict(list)
    for rec in read_gff(gff, {"gene"}):
        attrs = rec["attrs"]
        if attrs.get("gene_biotype", "protein_coding") != "protein_coding":
            continue
        gid = next((x for x in attrs.get("Dbxref", "").split(",") if x.startswith("GeneID:")), None)
        if gid is None and attrs.get("locus_tag"):
            gid = f"locus_tag:{attrs['locus_tag']}"
        if gid is None:
            continue
        by_seq[rec["seqid"]].append({
            "gene_id": gid,
            "start": rec["start"],
            "end": rec["end"],
            "strand": rec["strand"],
        })
    return by_seq


def genomic_features(genome_fasta: Path, genes_by_seq: dict[str, list[dict]]
                     ) -> Iterator[tuple[str, str, Counts]]:
    """Stream the genome one sequence at a time; yield (gene_id, feature, Counts)."""
    for header, seq in read_fasta(genome_fasta):
        seqid = header.split()[0]
        genes = genes_by_seq.get(seqid)
        if not genes:
            continue
        n = len(seq)
        for g in genes:
            s, e = g["start"] - 1, g["end"]  # 0-based half-open
            yield g["gene_id"], "gene_body", count(seq[s:e])
            if g["strand"] == "-":
                us, ue = e, min(n, e + UPSTREAM_BP)
            else:
                us, ue = max(0, s - UPSTREAM_BP), s
            if ue > us:
                yield g["gene_id"], "upstream", count(seq[us:ue])
