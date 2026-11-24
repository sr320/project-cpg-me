import csv
import gzip
import json

import numpy as np
import pytest

from cpgme import annotate as A

OBO = """format-version: 1.2
data-version: releases/test

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0006950
name: response to stress
namespace: biological_process
is_a: GO:0008150

[Term]
id: GO:0006260
name: DNA replication
namespace: biological_process
is_a: GO:0008150

[Term]
id: GO:0006979
name: response to oxidative stress
namespace: biological_process
alt_id: GO:9999999
is_a: GO:0006950

[Term]
id: GO:0006270
name: DNA replication initiation
namespace: biological_process
relationship: part_of GO:0006260

[Term]
id: GO:0000001
name: old term
namespace: biological_process
is_obsolete: true

[Term]
id: GO:0005634
name: nucleus
namespace: cellular_component

[Typedef]
id: part_of
name: part of
"""

SLIM_OBO = """format-version: 1.2
data-version: releases/test/slim

[Term]
id: GO:0006260
name: DNA replication
namespace: biological_process

[Term]
id: GO:0005634
name: nucleus
namespace: cellular_component
"""


@pytest.fixture
def onto(tmp_path):
    go_p, slim_p = tmp_path / "go-basic.obo", tmp_path / "goslim_generic.obo"
    go_p.write_text(OBO)
    slim_p.write_text(SLIM_OBO)
    return go_p, slim_p


def test_class_sets_disjoint_in_ontology():
    """No housekeeping term may be an ancestor or descendant of an inducible term in the
    real GO graph (checked only when the ontology has been downloaded)."""
    from pathlib import Path
    obo = Path("data/_ontology/go-basic.obo")
    if not obo.exists():
        pytest.skip("GO not downloaded")
    go = A.GO(obo)
    for a in A.HOUSEKEEPING:
        for b in A.INDUCIBLE:
            assert b not in go.ancestors(a) and a not in go.ancestors(b), (a, b)


def test_go_graph(onto):
    go = A.GO(onto[0])
    assert go.version == "releases/test"
    assert go.ancestors("GO:0006979") == {"GO:0006979", "GO:0006950", "GO:0008150"}
    assert "GO:0006260" in go.ancestors("GO:0006270")  # via part_of
    assert go.canonical("GO:9999999") == "GO:0006979"  # alt_id
    assert go.canonical("GO:0000001") is None  # obsolete
    assert go.canonical("GO:1234567") is None  # unknown


def test_slim_mapping(onto):
    go = A.GO(onto[0])
    slim = A.slim_terms(onto[1], go)
    # BP terms from the slim file plus the extras that exist; cellular_component dropped
    assert "GO:0006260" in slim and "GO:0005634" not in slim
    assert "GO:0006950" in slim
    assert A.map_to_slim(["GO:9999999"], go, slim) == {"GO:0006950"}
    assert A.map_to_slim(["GO:0006270", "GO:0000001"], go, slim) == {"GO:0006260"}
    assert A.map_to_slim([], go, slim) == set()


def test_write_gene_proteins(tmp_path):
    cds = tmp_path / "cds.fna"
    cds.write_text(
        ">c1 [gene=a] [db_xref=GeneID:1] [protein_id=XP_1.1]\nATGAAA\n"
        ">c2 [gene=a] [db_xref=GeneID:1] [protein_id=XP_2.1]\nATGAAACCCTAA\n"
        ">c3 [gene=b] [db_xref=GeneID:2] [protein_id=XP_3.1]\nATGTAA\n")
    prot = tmp_path / "protein.faa"
    prot.write_text(">XP_1.1 short\nMK\n>XP_2.1 long [Org]\nMKP\n>XP_3.1 b\nM\n>XP_9.1 x\nMM\n")
    out = tmp_path / "p.faa.gz"
    assert A.write_gene_proteins(cds, prot, out) == {"GeneID:1": "XP_2.1", "GeneID:2": "XP_3.1"}
    with gzip.open(out, "rt") as fh:
        assert fh.read() == ">GeneID:1 XP_2.1\nMKP\n>GeneID:2 XP_3.1\nM\n"


def _emapper_file(path, rows):
    cols = ["query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
            "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko"]
    with open(path, "w") as fh:
        fh.write("## emapper-2.1.12\n## command: emapper.py -i x\n##\n")
        fh.write("#" + "\t".join(cols) + "\n")
        for q, gos in rows:
            fh.write("\t".join([q, "1.x", "1e-50", "100", "OG@1", "Metazoa", "L", "desc",
                                "nm", gos or "-", "-", "ko:K1"]) + "\n")
        fh.write("## 3 queries scanned\n")


def test_parse_emapper(tmp_path):
    p = tmp_path / "x.emapper.annotations"
    _emapper_file(p, [("GeneID:1", "GO:0006260,GO:0006270"), ("GeneID:2", None)])
    rec = A.parse_emapper(p)
    assert set(rec) == {"GeneID:1", "GeneID:2"}
    assert rec["GeneID:1"]["GOs"] == "GO:0006260,GO:0006270"
    assert rec["GeneID:2"]["GOs"] == ""  # '-' -> ''
    assert A._emapper_header(p)[0] == "emapper-2.1.12"


def test_bh_matches_reference():
    p = np.array([0.01, 0.04, 0.03, np.nan, 0.5])
    q = A.bh(p)
    # sorted 0.01,0.03,0.04,0.5 with m=4: 0.04, 0.0533, 0.0533, 0.5
    assert q[0] == pytest.approx(0.04)
    assert q[2] == pytest.approx(0.16 / 3) and q[1] == pytest.approx(0.16 / 3)
    assert np.isnan(q[3]) and q[4] == pytest.approx(0.5)


def _synthetic(n=300, seed=0):
    """Housekeeping genes drawn low, inducible high, others mixed."""
    rng = np.random.default_rng(seed)
    cpg, gene_slim = {}, {}
    hk, ind = "GO:0006260", "GO:0007155"
    for i in range(n):
        g = f"GeneID:{i}"
        if i % 3 == 0:
            cpg[g], gene_slim[g] = (rng.normal(0.4, 0.08), 0.95), {hk}
        elif i % 3 == 1:
            cpg[g], gene_slim[g] = (rng.normal(0.9, 0.08), 0.05), {ind}
        else:
            cpg[g], gene_slim[g] = (rng.normal(0.7, 0.2), 0.5), {"GO:0005975"}
    gene_slim["GeneID:0"] = {hk, ind}  # in both classes: excluded from the class test
    return cpg, gene_slim


def test_class_test_direction():
    cpg, gene_slim = _synthetic()
    t = A.class_test(cpg, gene_slim)
    assert t["n_housekeeping"] == 99 and t["n_inducible"] == 100 and t["n_both_excluded"] == 1
    assert t["median_housekeeping"] < 0.5 < t["median_inducible"]
    assert t["frac_low_housekeeping"] == 1 and t["frac_low_inducible"] == 0
    assert t["rank_biserial"] < -0.9 and t["supports_2015"]
    # reversed data must not support the hypothesis
    flipped = {g: (2 - v, p) for g, (v, p) in cpg.items()}
    assert not A.class_test(flipped, gene_slim)["supports_2015"]


def test_goslim_stats():
    cpg, gene_slim = _synthetic()
    slim = {"GO:0006260": "DNA replication", "GO:0007155": "cell adhesion",
            "GO:0005975": "carbohydrate metabolic process", "GO:0006950": "response to stress"}
    rows = {r["goslim"]: r for r in A.goslim_stats(cpg, gene_slim, slim)}
    assert "GO:0006950" not in rows  # below MIN_TERM_GENES
    assert rows["GO:0006260"]["class"] == "housekeeping"
    assert rows["GO:0006260"]["rank_biserial"] < 0 < rows["GO:0007155"]["rank_biserial"]
    assert rows["GO:0006260"]["q"] < 1e-10


def test_annotate_end_to_end(tmp_path, onto, monkeypatch):
    from cpgme import fetch as F

    acc = "GCF_TEST.1"
    data = tmp_path / "data" / acc
    data.mkdir(parents=True)
    cds, prot = data / "cds_from_genomic.fna", data / "protein.faa"
    rng = np.random.default_rng(1)
    feat_rows, emap = [], []
    with open(cds, "w") as c, open(prot, "w") as p:
        for i in range(60):
            gid = f"GeneID:{i}"
            c.write(f">c{i} [gene=g{i}] [db_xref={gid}] [protein_id=XP_{i}.1]\nATGAAATAA\n")
            p.write(f">XP_{i}.1 p\nMK\n")
            hk = i % 2 == 0
            v = rng.normal(0.4 if hk else 0.9, 0.05)
            for feature in ("cds", "gene_body"):
                feat_rows.append({"gene_id": gid, "feature": feature, "cpg_oe": f"{v:.4f}",
                                  "pass_filter": 1, "p_low": 0.9 if hk else 0.1})
            emap.append((gid, "GO:0006270" if hk else "GO:0006979"))
    res = tmp_path / "results" / acc
    res.mkdir(parents=True)
    with gzip.open(res / "features.tsv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, ["gene_id", "feature", "cpg_oe", "pass_filter", "p_low"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(feat_rows)
    ann = tmp_path / "x.emapper.annotations"
    _emapper_file(ann, emap)

    monkeypatch.setattr(F, "fetch", lambda *a, **k: {"CDS_FASTA": cds, "PROT_FASTA": prot})
    monkeypatch.setattr(A, "ensure_ontology", lambda cache: onto)

    A.annotate(acc, tmp_path / "results", tmp_path / "data", emapper_annotations=ann,
               log=lambda *_: None)

    t = json.loads((res / "goslim_test.json").read_text())
    assert t["n_genes"] == 60 and t["n_with_goslim"] == 60
    for feature in ("cds", "gene_body"):
        r = t["tests"][feature]
        assert r["n_housekeeping"] == 30 and r["n_inducible"] == 30
        assert r["supports_2015"]
    with gzip.open(res / "annotation.tsv.gz", "rt") as fh:
        first = next(csv.DictReader(fh, delimiter="\t"))
    assert first["goslim"] == "GO:0006260" and first["goslim_names"] == "DNA replication"
    prov = json.loads((res / "annotation_provenance.json").read_text())
    assert prov["go_version"] == "releases/test"
    assert prov["emapper_header"][0] == "emapper-2.1.12"


def test_annotate_rejects_foreign_emapper(tmp_path, onto, monkeypatch):
    from cpgme import fetch as F

    acc = "GCF_TEST.1"
    res = tmp_path / "results" / acc
    res.mkdir(parents=True)
    with gzip.open(res / "features.tsv.gz", "wt") as fh:
        fh.write("gene_id\tfeature\tcpg_oe\tpass_filter\tp_low\n")
    cds, prot = tmp_path / "cds.fna", tmp_path / "prot.faa"
    cds.write_text(">c [gene=a] [db_xref=GeneID:1] [protein_id=XP_1.1]\nATG\n")
    prot.write_text(">XP_1.1\nM\n")
    ann = tmp_path / "x.emapper.annotations"
    _emapper_file(ann, [("XP_1.1", "GO:0006260"), ("XP_2.1", "")])
    monkeypatch.setattr(F, "fetch", lambda *a, **k: {"CDS_FASTA": cds, "PROT_FASTA": prot})
    monkeypatch.setattr(A, "ensure_ontology", lambda cache: onto)
    with pytest.raises(ValueError, match="not gene_ids"):
        A.annotate(acc, tmp_path / "results", tmp_path / "data", emapper_annotations=ann,
                   log=lambda *_: None)
