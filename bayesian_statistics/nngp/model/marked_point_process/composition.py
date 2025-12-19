"""Composition (mark) component for the marked point process model.

This module contains functions for:
- Computing κ̃ = y - N/2 for mark PG augmentation
- Sampling PG(N, η) for mark component
- Computing softmax probabilities
- Computing linear predictors η

Note: Most functions delegate to existing sample.py implementations
to ensure consistency with the existing multinomial model.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from pypolyagamma import PyPolyaGamma

# Import existing implementations for consistency
from bayesian_statistics.nngp.model.sample import (
    compute_eta as _original_compute_eta,
)
from bayesian_statistics.nngp.model.sample import (
    softmax_with_baseline as _original_softmax,
)


def compute_kappa_tilde_mark(
    counts_k: np.ndarray,
    total_counts: np.ndarray,
) -> np.ndarray:
    """Compute κ̃ for mark component.

    κ̃_ik = y_ik - N_i/2

    where y_ik is the count for category k at site i,
    and N_i is the total count at site i.

    This differs from the point process case (κ = y - 1/2).

    Parameters
    ----------
    counts_k : np.ndarray
        Counts for category k, shape (n_sites,).
    total_counts : np.ndarray
        Total counts per site, shape (n_sites,).

    Returns
    -------
    np.ndarray
        κ̃ values, shape (n_sites,).
    """
    return counts_k - total_counts / 2.0


def sample_xi_mark(
    eta_k: np.ndarray,
    total_counts: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample PG(N_i, η_ik) for mark component.

    For marks, we use b = N_i (total counts) instead of b = 1.
    ξ_ik ~ PG(N_i, η_ik) for i = 1,...,n

    Parameters
    ----------
    eta_k : np.ndarray
        Linear predictor for category k, shape (n_sites,).
    total_counts : np.ndarray
        Total counts per site (N_i), shape (n_sites,).
    rng : np.random.Generator
        Random number generator (not used directly, PG has its own RNG).

    Returns
    -------
    np.ndarray
        Sampled PG values, shape (n_sites,).
    """
    pg = PyPolyaGamma()
    n = len(eta_k)
    xi = np.zeros(n)
    for i in range(n):
        xi[i] = pg.pgdraw(float(total_counts[i]), float(eta_k[i]))
    return xi


def softmax_with_baseline(eta: np.ndarray) -> np.ndarray:
    """Convert linear predictors to probabilities.

    Delegates to sample.py implementation for consistency.

    The softmax is parameterised with the last category acting as the
    baseline. eta has shape (K-1, n). The returned array has
    shape (K, n) and includes the baseline probabilities in the final row.

    Parameters
    ----------
    eta : np.ndarray
        Linear predictors, shape (K-1, n).

    Returns
    -------
    np.ndarray
        Probabilities, shape (K, n). Last row is baseline.
    """
    return _original_softmax(eta)


def compute_eta(beta: np.ndarray, W: np.ndarray) -> np.ndarray:
    """Compute linear predictors for each category and location.

    Delegates to sample.py implementation for consistency.

    Parameters
    ----------
    beta : np.ndarray
        Coefficient array, shape (K-1, p+1, n).
    W : np.ndarray
        Design matrix, shape (n, p+1). First column must be ones.

    Returns
    -------
    np.ndarray
        Linear predictors, shape (K-1, n).
    """
    return _original_compute_eta(beta, W)
