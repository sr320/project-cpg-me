"""Gaussian mixture fits to a CpG O/E distribution.

Fits 1- and 2-component mixtures (as in Dimond & Roberts 2015, which used normalmixEM),
calls bimodality when BIC favours 2 components, the modes are well separated
(Ashman's D) and each component holds a real share of genes, and reports the low-CpG O/E component (historically heavily
methylated genes), the density intersection between components, and per-gene posterior
probability of belonging to the low component.
"""

from __future__ import annotations

import numpy as np
from sklearn.mixture import GaussianMixture

SEED = 20151010  # fixed so every machine reproduces the same fit
BIC_MIN = 10        # BIC(1) - BIC(2) must exceed this ...
ASHMAN_D_MIN = 2.0  # ... and Ashman's D must exceed this ...
WEIGHT_MIN = 0.05   # ... and each component must hold at least this share of genes
                    # (a handful of extreme outliers can otherwise get a component of their
                    # own with a huge D, e.g. one gene body at CpG O/E 55 in Rhynocoris)


def _fit(x: np.ndarray, k: int) -> GaussianMixture:
    return GaussianMixture(n_components=k, n_init=10, random_state=SEED).fit(x.reshape(-1, 1))


def _normal_pdf(x, mu, sd):
    return np.exp(-0.5 * ((x - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))


def _intersection(w, mu, sd) -> float | None:
    """x between the two means where the weighted component densities are equal."""
    lo, hi = sorted(mu)
    grid = np.linspace(lo, hi, 20001)
    diff = w[0] * _normal_pdf(grid, mu[0], sd[0]) - w[1] * _normal_pdf(grid, mu[1], sd[1])
    sign_change = np.nonzero(np.diff(np.sign(diff)))[0]
    if len(sign_change) == 0:
        return None
    return float(grid[sign_change[0]])


def fit_distribution(values: np.ndarray) -> tuple[dict, np.ndarray]:
    """Return (summary dict, p_low per value; NaN where the input is NaN)."""
    values = np.asarray(values, dtype=float)
    ok = np.isfinite(values)
    x = values[ok]
    p_low = np.full(values.shape, np.nan)
    summary: dict = {"n": int(x.size)}
    if x.size < 50:
        summary["note"] = "too few genes to fit"
        return summary, p_low

    summary.update({
        "mean": float(x.mean()), "median": float(np.median(x)), "sd": float(x.std(ddof=1)),
    })
    g1, g2 = _fit(x, 1), _fit(x, 2)
    bic1, bic2 = float(g1.bic(x.reshape(-1, 1))), float(g2.bic(x.reshape(-1, 1)))

    mu = g2.means_.ravel()
    sd = np.sqrt(g2.covariances_.ravel())
    w = g2.weights_.ravel()
    low = int(np.argmin(mu))
    high = 1 - low

    # Ashman's D: separation of the two modes; D > 2 is usually taken as clean bimodality.
    ashman_d = float(np.sqrt(2) * abs(mu[0] - mu[1]) / np.sqrt(sd[0] ** 2 + sd[1] ** 2))

    p_low[ok] = g2.predict_proba(x.reshape(-1, 1))[:, low]
    summary.update({
        "bic_1": bic1,
        "bic_2": bic2,
        "delta_bic": bic1 - bic2,  # > 0 favours two components
        # BIC alone over-calls with >10k genes (any skew favours 2 components), so a
        # bimodal call also needs the modes to be well separated.
        "bimodal_bic": bool(bic1 - bic2 > BIC_MIN),
        "bimodal": bool(bic1 - bic2 > BIC_MIN and ashman_d > ASHMAN_D_MIN
                        and w.min() >= WEIGHT_MIN),
        "low_mean": float(mu[low]), "low_sd": float(sd[low]), "low_weight": float(w[low]),
        "high_mean": float(mu[high]), "high_sd": float(sd[high]), "high_weight": float(w[high]),
        "ashman_d": ashman_d,
        "threshold": _intersection(w, mu, sd),
        "frac_low_posterior": float(np.mean(p_low[ok] > 0.5)),
    })
    return summary, p_low
