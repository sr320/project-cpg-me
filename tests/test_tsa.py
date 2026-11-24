import gzip
import random

from cpgme.tsa import STOPS, is_tsa, longest_orf, master_info, revcomp, translate, tsa_cds

_CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
           if a + b + c not in STOPS and a + b + c != "ATG"]


def body(n: int, seed: int = 1) -> str:
    """n random sense codons (no stop, no ATG) whose other five frames all contain stops."""
    rng = random.Random(seed)
    while True:
        s = "".join(rng.choice(_CODONS) for _ in range(n))
        r = revcomp(s)
        if all(any(x[i:i + 3] in STOPS for i in range(f, len(x) - 2, 3))
               for x, frames in ((s, (1, 2)), (r, (0, 1, 2))) for f in frames):
            return s


def test_is_tsa():
    assert is_tsa("GCWQ01000000")
    assert not is_tsa("GCF_013753865.1")
    assert not is_tsa("GCWQ00000000.1")


def test_translate():
    assert translate("ATGGCCTAA") == "MA*"


def test_complete_orf_plus_strand():
    cds = "ATG" + body(120) + "TAA"
    seq = "TTTAA" + cds + "GGGTT"
    orf = longest_orf(seq, min_nt=300)
    assert orf["seq"] == cds and orf["strand"] == "+"
    assert not orf["three_partial"]


def test_orf_on_minus_strand():
    cds = "ATG" + body(120) + "TGA"
    orf = longest_orf("CC" + revcomp(cds) + "A", min_nt=300)
    assert orf["seq"] == cds and orf["strand"] == "-"


def test_five_prime_partial_runs_from_contig_start():
    # no stop upstream: ORF starts at the first codon, not the first ATG
    seq = body(50) + "ATG" + body(60, seed=2) + "TAG"
    orf = longest_orf(seq, min_nt=300)
    assert orf["start"] == 0 and orf["five_partial"]
    assert orf["seq"] == seq


def test_after_stop_needs_atg():
    seq = "TAA" + body(50) + "ATG" + body(100, seed=2) + "TAG"
    orf = longest_orf(seq, min_nt=300)
    assert orf["seq"].startswith("ATG") and not orf["five_partial"]
    assert len(orf["seq"]) == 3 + 300 + 3


def test_three_prime_partial():
    seq = "TAA" + "ATG" + body(120)  # runs off the end
    orf = longest_orf(seq, min_nt=300)
    assert orf["three_partial"] and orf["seq"] == "ATG" + body(120)


def test_too_short_returns_none():
    assert longest_orf("TAA" + "ATG" + "GCC" * 20 + "TAA", min_nt=300) is None


def test_tsa_cds_drops_duplicate_proteins(tmp_path):
    cds = "ATG" + body(120) + "TAA"
    fa = tmp_path / "tsa.fsa_nt.gz"
    with gzip.open(fa, "wt") as fh:
        fh.write(f">C1.1 contig one\nTT{cds}AA\n>C2.1 dup on other strand\n{revcomp(cds)}\n"
                 f">C3.1 short\nTAAATGGCCTAA\n")
    out = tsa_cds(fa)
    assert list(out) == ["C1.1"]
    assert out["C1.1"]["protein"] == translate(cds)[:-1]


def test_master_info(tmp_path):
    p = tmp_path / "m.gbff.gz"
    with gzip.open(p, "wt") as fh:
        fh.write('SOURCE      Cyrpoptus belfragei\n  ORGANISM  Cyrpoptus belfragei\n'
                 '            Eukaryota; Metazoa.\n'
                 '            Assembly Method       :: SOAPdenovo-Trans v1.02\n'
                 '            Assembly Name         :: Cybel\n'
                 '                     /db_xref="taxon:1464866"\n')
    info = master_info(p)
    assert info == {"organism": "Cyrpoptus belfragei", "taxid": 1464866,
                    "assembly_method": "SOAPdenovo-Trans v1.02", "assembly_name": "Cybel"}
