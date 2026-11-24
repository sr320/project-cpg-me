"""Dinucleotide counts and observed/expected ratios.

CpG O/E = (n_CG / (n_C * n_G)) * (L^2 / (L - 1))   (Gavery & Roberts 2010; Dimond & Roberts 2015)

L is the length after removing ambiguous bases. GpC O/E is computed the same way and
serves as a composition control: it should sit near 1 regardless of methylation history.
Both ratios are strand-symmetric, so sequences never need reverse-complementing.
"""

from __future__ import annotations

from dataclasses import dataclass

ACGT = frozenset("ACGT")


@dataclass(frozen=True)
class Counts:
    length: int  # unambiguous bases (A, C, G, T)
    n_ambig: int
    n_C: int
    n_G: int
    n_CG: int
    n_GC: int

    @property
    def gc(self) -> float:
        return (self.n_C + self.n_G) / self.length if self.length else float("nan")

    @property
    def frac_ambig(self) -> float:
        total = self.length + self.n_ambig
        return self.n_ambig / total if total else float("nan")

    def _oe(self, n_obs: int) -> float:
        L = self.length
        if L < 2 or self.n_C == 0 or self.n_G == 0:
            return float("nan")
        return (n_obs / (self.n_C * self.n_G)) * (L * L / (L - 1))

    @property
    def cpg_oe(self) -> float:
        return self._oe(self.n_CG)

    @property
    def gpc_oe(self) -> float:
        return self._oe(self.n_GC)


def count(seq: str) -> Counts:
    """Count bases and dinucleotides in one sequence (case-insensitive).

    Dinucleotides spanning an ambiguous base are not counted, so N runs never create
    or hide CG/GC steps.
    """
    s = seq.upper()
    n_A, n_C, n_G, n_T = s.count("A"), s.count("C"), s.count("G"), s.count("T")
    length = n_A + n_C + n_G + n_T
    # "CG" and "GC" cannot overlap with themselves, so str.count is exact.
    return Counts(
        length=length,
        n_ambig=len(s) - length,
        n_C=n_C,
        n_G=n_G,
        n_CG=s.count("CG"),
        n_GC=s.count("GC"),
    )
