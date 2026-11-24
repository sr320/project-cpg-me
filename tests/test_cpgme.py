import gzip
import math

import numpy as np
import pytest

from cpgme.cli import in_shard
from cpgme.compute import count
from cpgme.features import genomic_features, longest_cds, protein_coding_genes
from cpgme.fit import fit_distribution


# ---- compute -------------------------------------------------------------------------

def test_cpg_oe_hand_computed():
    c = count("ACGTCGAT")  # L=8, C=2, G=2, CG=2, GC=0
    assert (c.length, c.n_C, c.n_G, c.n_CG, c.n_GC) == (8, 2, 2, 2, 0)
    assert c.cpg_oe == pytest.approx((2 / (2 * 2)) * (64 / 7))
    assert c.gpc_oe == 0


def test_case_and_ambiguity():
    c = count("acgNNcg")
    assert c.length == 5 and c.n_ambig == 2
    assert c.n_CG == 2  # soft-masked bases count; no CG across the N run


def test_strand_symmetry():
    seq = "ATGCGCGTTACGGATCCAGTCGA" * 5
    rc = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
    assert count(seq).cpg_oe == pytest.approx(count(rc).cpg_oe)


def test_random_sequence_near_one():
    rng = np.random.default_rng(1)
    seq = "".join(rng.choice(list("ACGT"), 200_000))
    c = count(seq)
    assert c.cpg_oe == pytest.approx(1, abs=0.03)
    assert c.gpc_oe == pytest.approx(1, abs=0.03)


def test_no_c_gives_nan():
    assert math.isnan(count("ATATATGG").cpg_oe)


# ---- features ------------------------------------------------------------------------

@pytest.fixture
def fixture_files(tmp_path):
    cds = tmp_path / "cds.fna"
    cds.write_text(
        ">lcl|chr1_cds_XP_1.1_1 [gene=geneA] [db_xref=GeneID:1] [protein_id=XP_1.1]\nATGCGA\n"
        ">lcl|chr1_cds_XP_2.1_2 [gene=geneA] [db_xref=GeneID:1] [protein_id=XP_2.1]\nATGCGACGT\nTAA\n"
        ">lcl|chr1_cds_XP_3.1_3 [gene=geneB] [db_xref=GeneID:2] [protein_id=XP_3.1]\nATGTAA\n"
        ">lcl|chr1_cds_4 [gene=geneC] [db_xref=GeneID:3] [pseudo=true]\nATGCCCCCCCCCCTAA\n"
    )
    genome = tmp_path / "genome.fna.gz"
    chrom = "A" * 10 + "CGCGCGCGCG" + "T" * 10 + "GGGGGCCCCC" + "A" * 10  # 50 bp
    with gzip.open(genome, "wt") as fh:
        fh.write(">chr1 test\n" + chrom[:25] + "\n" + chrom[25:] + "\n")
    gff = tmp_path / "genomic.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tRefSeq\tgene\t11\t20\t.\t+\t.\tID=gene-A;Dbxref=GeneID:1;gene_biotype=protein_coding\n"
        "chr1\tRefSeq\tgene\t31\t40\t.\t-\t.\tID=gene-B;Dbxref=GeneID:2;gene_biotype=protein_coding\n"
        "chr1\tRefSeq\tgene\t1\t5\t.\t+\t.\tID=gene-R;Dbxref=GeneID:9;gene_biotype=lncRNA\n"
        "chr1\tRefSeq\tmRNA\t11\t20\t.\t+\t.\tID=rna-A;Dbxref=GeneID:1\n"
    )
    return cds, genome, gff, chrom


def test_longest_cds(fixture_files):
    cds, *_ = fixture_files
    best = longest_cds(cds)
    assert set(best) == {"GeneID:1", "GeneID:2"}  # pseudogene dropped
    assert best["GeneID:1"]["seq"] == "ATGCGACGTTAA"  # multi-line record joined
    assert best["GeneID:1"]["protein_id"] == "XP_2.1"
    assert best["GeneID:1"]["n_isoforms"] == 2


def test_gene_body_and_upstream(fixture_files):
    _, genome, gff, chrom = fixture_files
    genes = protein_coding_genes(gff)
    assert [g["gene_id"] for g in genes["chr1"]] == ["GeneID:1", "GeneID:2"]
    out = {(gid, f): c for gid, f, c in genomic_features(genome, genes)}
    assert out[("GeneID:1", "gene_body")] == count(chrom[10:20])
    assert out[("GeneID:1", "upstream")] == count(chrom[0:10])  # + strand: before start
    assert out[("GeneID:2", "gene_body")] == count(chrom[30:40])
    assert out[("GeneID:2", "upstream")] == count(chrom[40:50])  # - strand: after end


# ---- fit -----------------------------------------------------------------------------

def test_bimodal_recovered():
    rng = np.random.default_rng(0)
    x = np.concatenate([rng.normal(0.35, 0.08, 3000), rng.normal(0.85, 0.1, 5000)])
    s, p = fit_distribution(x)
    assert s["bimodal"]
    assert s["low_mean"] == pytest.approx(0.35, abs=0.02)
    assert s["high_mean"] == pytest.approx(0.85, abs=0.02)
    assert s["low_weight"] == pytest.approx(3 / 8, abs=0.03)
    assert 0.35 < s["threshold"] < 0.85
    assert p[:3000].mean() > 0.9 and p[3000:].mean() < 0.1


def test_unimodal_not_called_bimodal():
    x = np.random.default_rng(0).normal(0.9, 0.12, 8000)
    s, _ = fit_distribution(x)
    assert not s["bimodal"]


def test_nan_passthrough():
    x = np.r_[np.random.default_rng(0).normal(0.8, 0.1, 200), np.nan]
    _, p = fit_distribution(x)
    assert math.isnan(p[-1]) and np.isfinite(p[:-1]).all()


# ---- sharding ------------------------------------------------------------------------

def test_shards_partition():
    accs = [f"GCF_{i:09d}.1" for i in range(500)]
    for n in (1, 3, 7):
        assigned = [sum(in_shard(a, k, n) for k in range(1, n + 1)) for a in accs]
        assert all(x == 1 for x in assigned)


# ---- download watchdog ---------------------------------------------------------------

def test_watchdog_aborts_stalled_download(monkeypatch):
    import io
    import threading

    from cpgme import fetch

    class Stalled:
        def __init__(self):
            self.sent, self.closed = False, threading.Event()

        def read(self, n):
            if not self.sent:
                self.sent = True
                return b"x" * 10
            self.closed.wait(30)  # blocks until the watchdog closes us
            raise OSError("connection closed")

        def close(self):
            self.closed.set()

    monkeypatch.setattr(fetch, "STALL_SECONDS", 1)
    with pytest.raises(TimeoutError):
        fetch._copy_with_watchdog(Stalled(), io.BytesIO())


def test_skewed_unimodal_large_n_not_bimodal():
    # A single skewed peak: BIC still prefers 2 components, but they overlap (D < 2).
    x = np.random.default_rng(0).gamma(9, 0.1, 15000)
    s, _ = fit_distribution(x)
    assert s["bimodal_bic"] and s["ashman_d"] < 2
    assert not s["bimodal"]


def test_locus_tag_fallback_for_genbank_annotation(tmp_path):
    """Submitter annotations have no GeneID; CDS and GFF genes link through locus_tag."""
    cds = tmp_path / "cds.fna"
    cds.write_text(
        ">lcl|s1_cds_RZF1.1_1 [locus_tag=LS_001] [protein_id=RZF1.1]\nATGCGA\n"
        ">lcl|s1_cds_RZF2.1_2 [locus_tag=LS_001] [protein_id=RZF2.1]\nATGCGACGTTAA\n")
    gff = tmp_path / "g.gff"
    gff.write_text(
        "s1\tGenbank\tgene\t11\t20\t.\t+\t.\tID=gene-LS_001;gene_biotype=protein_coding;"
        "locus_tag=LS_001\n")
    from cpgme.features import longest_cds, protein_coding_genes
    best = longest_cds(cds)
    assert set(best) == {"locus_tag:LS_001"}
    assert best["locus_tag:LS_001"]["protein_id"] == "RZF2.1"
    genes = protein_coding_genes(gff)
    assert [g["gene_id"] for g in genes["s1"]] == ["locus_tag:LS_001"]
