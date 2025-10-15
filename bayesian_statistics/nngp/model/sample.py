"""Utility routines for multinomial NNGP inference.

This module collects the reusable numerical pieces that appear in the
multinomial NNGP Gibbs sampler:

- kernel helpers for building Vecchia/NNGP factors
- vectorised helpers for computing softmax probabilities
- conditional Gaussian updates for the coefficient field under the
  Pólya–Gamma augmentation
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np
from pypolyagamma import PyPolyaGamma
from scipy.special import logsumexp

from .nngp import NNGPFactors


def rbf_kernel(
    X1: np.ndarray,
    X2: np.ndarray,
    lengthscale: float,
    variance: float,
) -> np.ndarray:
    """Squared exponential (RBF) kernel."""
    X1 = np.atleast_2d(X1)
    X2 = np.atleast_2d(X2)
    a2 = np.sum(X1**2, axis=1, keepdims=True)
    b2 = np.sum(X2**2, axis=1, keepdims=True).T
    sqd = a2 + b2 - 2.0 * (X1 @ X2.T)
    return variance * np.exp(-0.5 * sqd / (lengthscale**2))


@dataclass
class LocalNNGPKernel:
    """Isotropic RBF kernel with convenience wrapper."""

    lengthscale: float = 0.15
    variance: float = 1.0

    def K(self, X1: np.ndarray, X2: np.ndarray) -> np.ndarray:
        return rbf_kernel(X1, X2, self.lengthscale, self.variance)


class PolyaGammaSampler:
    """Thin wrapper around :class:`pypolyagamma.PyPolyaGamma`.

    The external library is not vectorised, so we provide a small helper to
    draw many independent samples while keeping the calling code tidy.
    """

    def __init__(self) -> None:
        self._pg = PyPolyaGamma()

    def sample(self, b: np.ndarray, c: np.ndarray) -> np.ndarray:
        b = np.asarray(b, dtype=float)
        c = np.asarray(c, dtype=float)
        if b.shape != c.shape:
            raise ValueError("shape mismatch between b and c in Polya-Gamma draw")
        out = np.empty_like(c, dtype=float)
        flat_b = b.ravel()
        flat_c = c.ravel()
        flat_out = out.ravel()
        for idx in range(flat_out.size):
            flat_out[idx] = self._pg.pgdraw(flat_b[idx], flat_c[idx])
        return out


def compute_eta(beta: np.ndarray, W: np.ndarray) -> np.ndarray:
    """Compute linear predictors for each category and location.

    Parameters
    ----------
    beta
        Array with shape (K-1, p+1, n). The first axis corresponds to the
        non-baseline categories, the second axis to the intercept plus
        covariate coefficients, and the final axis to spatial locations.
    W
        Design matrix of shape (n, p+1). The first column must be ones.
    """
    K_minus_1, _, n_obs = beta.shape
    eta = np.zeros((K_minus_1, n_obs), dtype=float)
    for k in range(K_minus_1):
        # beta[k] has shape (p+1, n); W is (n, p+1)
        eta[k] = np.sum(W * beta[k].T, axis=1)
    return eta


def softmax_with_baseline(eta: np.ndarray) -> np.ndarray:
    """Convert linear predictors to probabilities.

    The softmax is parameterised with the last category acting as the
    baseline. ``eta`` therefore has shape (K-1, n). The returned array has
    shape (K, n) and includes the baseline probabilities in the final row.
    """
    exp_eta = np.exp(eta - eta.max(axis=0, keepdims=True))
    sum_exp = exp_eta.sum(axis=0, keepdims=True)
    baseline = 1.0 + sum_exp
    probs_non_baseline = exp_eta / baseline
    probs_baseline = 1.0 / baseline
    return np.vstack([probs_non_baseline, probs_baseline])


def compute_log_rest_terms(eta: np.ndarray) -> np.ndarray:
    """Compute ``log(1 + sum_{ℓ≠k} exp(η_ℓ))`` for every category/location.

    Parameters
    ----------
    eta
        Linear predictors of shape (K-1, n).

    Returns
    -------
    np.ndarray
        Array of shape (K-1, n) containing the desired log-rest terms.
    """
    if eta.ndim != 2:
        raise ValueError("eta must be a 2D array with shape (K-1, n)")
    # logsumexp over baseline (0) plus all categories
    stacked = np.vstack([np.zeros((1, eta.shape[1])), eta])
    log_sum_total = logsumexp(stacked, axis=0)
    # For each category compute log( sum_total - exp(eta_k) ) in a numerically stable way
    diff = np.clip(np.exp(eta - log_sum_total), 0.0, 1.0 - 1e-12)
    log_rest = log_sum_total + np.log1p(-diff)
    return log_rest


def update_beta_category(
    beta_k: np.ndarray,
    eta_k: np.ndarray,
    W: np.ndarray,
    factors_by_feature: Sequence[NNGPFactors],
    order: np.ndarray,
    omega: np.ndarray,
    kappa_tilde: np.ndarray,
    rng: np.random.Generator,
    min_variance: float = 1e-12,
) -> Tuple[np.ndarray, np.ndarray]:
    """Gibbs update for a single category across all locations.

    Parameters
    ----------
    beta_k
        Coefficient surface for category ``k`` with shape (p+1, n).
    eta_k
        Linear predictor for category ``k`` with shape (n,). Updated in place.
    W
        Design matrix with intercept as the first column.
    factors
        Vecchia/NNGP factors for the observed locations.
    order
        Ordering used when building ``factors``; ensures neighbours precede the
        current location.
    omega
        Pólya–Gamma draws for this category, shape (n,).
    kappa_tilde
        ``κ + ω · C`` terms (see documentation), shape (n,).
    rng
        NumPy random generator used for sampling.
    min_variance
        Lower bound to keep the prior conditional variance positive.
    """
    p_plus_1, n_obs = beta_k.shape
    if W.shape != (n_obs, p_plus_1):
        raise ValueError("shape mismatch between beta and design matrix W")
    if len(factors_by_feature) != p_plus_1:
        raise ValueError("expected one NNGP factor set per feature column")

    beta_new = beta_k.copy()
    eta_new = eta_k.copy()

    for j in range(p_plus_1):
        w_col = W[:, j]
        factors = factors_by_feature[j]
        for idx in order:
            neighbours = factors.neighbor_idx[idx]
            if neighbours.size:
                mu_prior = float(factors.a_rows[idx] @ beta_new[j, neighbours])
            else:
                mu_prior = 0.0
            d_i = float(max(factors.d[idx], min_variance))
            w_ij = float(w_col[idx])

            if w_ij == 0.0:
                mean = mu_prior
                var = d_i
            else:
                eta_without = eta_new[idx] - w_ij * beta_new[j, idx]
                omega_i = float(omega[idx])
                denom = omega_i * w_ij * w_ij + 1.0 / d_i
                var = 1.0 / denom
                linear = (float(kappa_tilde[idx]) - omega_i * eta_without) * w_ij
                mean = var * (linear + mu_prior / d_i)

            beta_draw = rng.normal(mean, np.sqrt(max(var, min_variance)))
            if w_ij != 0.0:
                eta_new[idx] = (
                    eta_new[idx] - w_ij * beta_new[j, idx]
                ) + w_ij * beta_draw
            beta_new[j, idx] = beta_draw

    return beta_new, eta_new
