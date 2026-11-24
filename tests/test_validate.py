import gzip
import json

import numpy as np
import pytest

from cpgme import validate as V


def _cov(path, rows):
    with gzip.open(path, "wt") as fh:
        for chrom, pos, m, u in rows:
            pct = 100 * m / (m + u) if m + u else 0
            fh.write(f"{chrom}\t{pos}\t{pos}\t{pct:.2f}\t{m}\t{u}\n")


def test_read_coverage_filters_and_sorts(tmp_path):
    p = tmp_path / "x.cov.gz"
    _cov(p, [("c1", 30, 5, 5), ("c1", 10, 1, 1), ("c1", 20, 9, 1), ("c2", 5, 0, 6)])
    cov = V.read_coverage(p, min_cov=5)
    pos, m, t = cov["c1"]
    assert pos.tolist() == [20, 30] and m.tolist() == [9, 5] and t.tolist() == [10, 10]
    assert cov["c2"][0].tolist() == [5]


def test_gene_methylation_bounds_and_weighting():
    cov = {"c1": (np.array([10, 20, 30, 40, 50, 60]), np.array([10, 10, 0, 0, 5, 9]),
                  np.array([10, 10, 10, 10, 10, 10]))}
    genes = {"c1": [{"gene_id": "A", "start": 10, "end": 50, "strand": "+"},   # 5 CpGs
                    {"gene_id": "B", "start": 11, "end": 50, "strand": "+"},   # 4: too few
                    {"gene_id": "C", "start": 1, "end": 60, "strand": "-"}],
             "c9": [{"gene_id": "D", "start": 1, "end": 100, "strand": "+"}]}
    out = V.gene_methylation(cov, genes, min_cpgs=5)
    assert set(out) == {"A", "C"}
    assert out["A"]["n_cpg"] == 5 and out["A"]["meth"] == pytest.approx(25 / 50)
    assert out["A"]["frac_mcpg"] == pytest.approx(3 / 5)  # 1.0, 1.0, 0.5
    assert out["C"]["meth"] == pytest.approx(34 / 60)


def test_compare_detects_negative_relationship():
    rng = np.random.default_rng(0)
    cpg, meth = {}, {}
    for i in range(400):
        methylated = i % 2 == 0
        cpg[f"g{i}"] = (rng.normal(0.45 if methylated else 0.9, 0.08), 0.95 if methylated else 0.05)
        meth[f"g{i}"] = {"meth": rng.uniform(0.6, 0.95) if methylated else rng.uniform(0, 0.1)}
    s = V.compare(cpg, meth)
    assert s["n_genes"] == 400 and s["spearman_rho"] < -0.7
    assert s["auc_low_oe_predicts_methylated"] > 0.95
    lc = s["low_class"]
    assert lc["misclassification"] == 0 and lc["tp"] == 200 and lc["tn"] == 200
    # no relationship: AUC near 0.5
    shuffled = dict(zip(meth, rng.permutation(list(meth.values()))))
    assert abs(V.compare(cpg, shuffled)["auc_low_oe_predicts_methylated"] - 0.5) < 0.1


def test_validate_end_to_end(tmp_path):
    acc = "GCF_T.1"
    data = tmp_path / "data" / acc
    data.mkdir(parents=True)
    gff_lines, feat_rows, cov_rows = ["##gff-version 3"], [], []
    for i in range(40):
        s = 1000 * i + 1
        gff_lines.append(f"chr1\tRefSeq\tgene\t{s}\t{s + 499}\t.\t+\t.\t"
                         f"ID=g{i};Dbxref=GeneID:{i};gene_biotype=protein_coding")
        m = i % 2 == 0
        feat_rows.append(f"GeneID:{i}\tgene_body\t{0.4 if m else 0.9}\t1\t{0.9 if m else 0.1}")
        cov_rows += [("chr1", s + k * 10, 9 if m else 0, 1 if m else 10) for k in range(8)]
    (data / "genomic.gff").write_text("\n".join(gff_lines) + "\n")
    res = tmp_path / "results" / acc
    res.mkdir(parents=True)
    with gzip.open(res / "features.tsv.gz", "wt") as fh:
        fh.write("gene_id\tfeature\tcpg_oe\tpass_filter\tp_low\n" + "\n".join(feat_rows) + "\n")
    cov = tmp_path / "s.bismark.cov.gz"
    _cov(cov, cov_rows)
    vdir = V.validate(acc, tmp_path / "results", tmp_path / "data", cov, "test",
                      log=lambda *_: None)
    v = json.loads((vdir / "validation.json").read_text())
    s = v["features"]["gene_body"]
    assert v["n_genes_with_meth"] == 40 and s["n_genes"] == 40
    assert s["spearman_rho"] < -0.8 and s["low_class"]["misclassification"] == 0


def test_validate_rejects_other_assembly(tmp_path):
    acc = "GCF_T.1"
    (tmp_path / "data" / acc).mkdir(parents=True)
    (tmp_path / "data" / acc / "genomic.gff").write_text(
        "chr1\tR\tgene\t1\t100\t.\t+\t.\tDbxref=GeneID:1;gene_biotype=protein_coding\n")
    (tmp_path / "results" / acc).mkdir(parents=True)
    with gzip.open(tmp_path / "results" / acc / "features.tsv.gz", "wt") as fh:
        fh.write("gene_id\tfeature\tcpg_oe\tpass_filter\tp_low\n")
    cov = tmp_path / "s.cov.gz"
    _cov(cov, [("NC_other", 5, 5, 5)])
    with pytest.raises(ValueError, match="same assembly"):
        V.validate(acc, tmp_path / "results", tmp_path / "data", cov, "x", log=lambda *_: None)
