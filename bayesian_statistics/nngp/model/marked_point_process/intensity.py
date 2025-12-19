"""Intensity (point process) component for the marked point process model.

This module contains functions for:
- Computing existence probability q(s,t) = sigmoid(η_int)
- Sampling pseudo-absence points via Poisson thinning
- Updating λ* from Gamma posterior
- Sampling PG(1, η) for point process
- Updating (β_int, u_int) via NNGP Gibbs step
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np
from pypolyagamma import PyPolyaGamma

from bayesian_statistics.nngp.model.nngp import NNGPFactors
from bayesian_statistics.nngp.model.sample import update_beta_category


def compute_q(eta: np.ndarray) -> np.ndarray:
    """Compute existence probability q(s,t) = sigmoid(η).

    Parameters
    ----------
    eta : np.ndarray
        Linear predictor values.

    Returns
    -------
    np.ndarray
        Probability values in [0, 1].
    """
    return 1.0 / (1.0 + np.exp(-eta))


def sample_pseudo_absence(
    lambda_star: float,
    eta_grid: np.ndarray,
    valid_mask: np.ndarray,
    region_volume: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample pseudo-absence points via Poisson thinning.

    U ~ IPP(λ*(1-q)) using Poisson thinning algorithm:
    1. N ~ Poisson(λ* × |D|)
    2. Sample N candidate points uniformly from valid grid
    3. Accept each with probability 1-q(s,t)

    Parameters
    ----------
    lambda_star : float
        Intensity upper bound parameter.
    eta_grid : np.ndarray
        Linear predictor at all grid points, shape (n_grid,).
    valid_mask : np.ndarray
        Boolean mask indicating valid grid cells, shape (n_grid,).
    region_volume : float
        Area of the spatial region |D|.
    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    np.ndarray
        Indices of accepted pseudo-absence points (into grid array).
    """
    # Step 1: Sample number of candidates from Poisson(λ* × |D|)
    n_candidates = rng.poisson(lambda_star * region_volume)

    if n_candidates == 0:
        return np.array([], dtype=int)

    # Step 2: Sample candidate indices from valid grid
    valid_indices = np.where(valid_mask)[0]
    if len(valid_indices) == 0:
        return np.array([], dtype=int)

    candidate_indices = rng.choice(valid_indices, size=n_candidates, replace=True)

    # Step 3: Thinning - accept with probability 1-q(s,t)
    q_candidates = compute_q(eta_grid[candidate_indices])
    u = rng.uniform(size=n_candidates)
    accept_mask = u < (1.0 - q_candidates)

    return candidate_indices[accept_mask]


def update_lambda_star(
    n_total: int,
    region_volume: float,
    prior_shape: float,
    prior_rate: float,
    rng: np.random.Generator,
) -> float:
    """Sample λ* from its Gamma posterior.

    λ* | X, U ~ Gamma(m_0 + n, r_0 + |D|)

    where n = n_X + n_U is the total number of points.

    Parameters
    ----------
    n_total : int
        Total number of points (sites + pseudo-absence).
    region_volume : float
        Area of the spatial region |D|.
    prior_shape : float
        Shape parameter m_0 of Gamma prior.
    prior_rate : float
        Rate parameter r_0 of Gamma prior.
    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    float
        Sampled λ* value.
    """
    posterior_shape = prior_shape + n_total
    posterior_rate = prior_rate + region_volume
    return rng.gamma(posterior_shape, 1.0 / posterior_rate)


def sample_omega_intensity(
    eta: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample PG(1, η) for point process.

    For point process, we use b=1 for all points.
    ω_i ~ PG(1, η_i) for i = 1,...,n

    Parameters
    ----------
    eta : np.ndarray
        Linear predictor values, shape (n,).
    rng : np.random.Generator
        Random number generator (not used directly, PG has its own RNG).

    Returns
    -------
    np.ndarray
        Sampled PG values, shape (n,).
    """
    pg = PyPolyaGamma()
    n = len(eta)
    omega = np.zeros(n)
    for i in range(n):
        omega[i] = pg.pgdraw(1.0, eta[i])
    return omega


def compute_kappa_intensity(y: np.ndarray) -> np.ndarray:
    """Compute κ for point process.

    κ_i = y_i - 1/2
    where y_i = 1 for observed sites (X), y_i = 0 for pseudo-absence (U).

    Parameters
    ----------
    y : np.ndarray
        Binary indicator: 1 for sites, 0 for pseudo-absence.

    Returns
    -------
    np.ndarray
        κ values.
    """
    return y - 0.5


def compute_eta_intensity(beta: np.ndarray, W: np.ndarray) -> np.ndarray:
    """Compute linear predictor η_int for the intensity component.

    η_int,i = Σ_j W[i,j] × β[j,i]

    This is analogous to compute_eta in sample.py but for a single
    category (the point process intensity).

    Parameters
    ----------
    beta : np.ndarray
        Coefficient array with shape (p+1, n).
    W : np.ndarray
        Design matrix with shape (n, p+1).

    Returns
    -------
    np.ndarray
        Linear predictor values with shape (n,).
    """
    # eta[i] = sum_j W[i,j] * beta[j,i]
    return np.sum(W * beta.T, axis=1)


def update_beta_intensity(
    beta_int: np.ndarray,
    eta_int: np.ndarray,
    W: np.ndarray,
    factors_by_feature: Sequence[NNGPFactors],
    order: np.ndarray,
    omega: np.ndarray,
    kappa: np.ndarray,
    rng: np.random.Generator,
    min_variance: float = 1e-12,
    prior_mean_by_feature: Optional[Sequence[Optional[np.ndarray]]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Gibbs update for intensity coefficients (step (c) in sec7.tex).

    This implements the NNGP-based Gibbs update for (β_int, u_int).
    The spatial effects are embedded in the site-specific coefficients β_int[j, i].

    The update is:
        θ_int = (β_int, u_int) ~ N(m_int, V_int)

    where m_int and V_int are computed from ω, X, U.

    This function delegates to update_beta_category for the actual computation,
    following the DRY principle.

    Parameters
    ----------
    beta_int : np.ndarray
        Current intensity coefficients, shape (p+1, n).
    eta_int : np.ndarray
        Current linear predictor, shape (n,).
    W : np.ndarray
        Design matrix with intercept, shape (n, p+1).
    factors_by_feature : Sequence[NNGPFactors]
        NNGP factors for each feature (length p+1).
    order : np.ndarray
        Morton ordering for NNGP traversal.
    omega : np.ndarray
        Polya-Gamma samples ω_i ~ PG(1, η_int,i), shape (n,).
    kappa : np.ndarray
        κ_i = y_i - 0.5 where y_i = 1 for X, 0 for U, shape (n,).
    rng : np.random.Generator
        Random number generator.
    min_variance : float
        Minimum variance for numerical stability.
    prior_mean_by_feature : Optional[Sequence[Optional[np.ndarray]]]
        Optional prior means for each feature.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Updated (beta_int, eta_int).
    """
    # Delegate to the existing update_beta_category function
    # This ensures DRY principle and consistent behavior with mark updates
    return update_beta_category(
        beta_k=beta_int,
        eta_k=eta_int,
        W=W,
        factors_by_feature=factors_by_feature,
        order=order,
        omega=omega,
        kappa_tilde=kappa,
        rng=rng,
        min_variance=min_variance,
        prior_mean_by_feature=prior_mean_by_feature,
    )
