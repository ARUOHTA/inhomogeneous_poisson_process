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
from typing import Dict, Optional, Sequence, Tuple

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


@dataclass
class DistanceDependentNNGPKernel:
    """Distance-dependent variance RBF kernel.

    A kernel where the local standard deviation increases with distance from
    an origin point (e.g., obsidian source). Points farther from the origin
    have larger variance, allowing them to exert broader spatial influence.

    This addresses the archaeological intuition that distant observations
    imply transport routes must exist, and should influence a wider spatial
    region to ensure route continuity.

    This class stores pre-computed distance information and presents a standard
    K(X1, X2) interface for compatibility with existing NNGP infrastructure.

    Parameters
    ----------
    lengthscale : float
        RBF kernel length scale (spatial correlation range).
    base_variance : float
        Base variance near the origin.
    distance_scaling : float
        Distance scaling coefficient gamma. When gamma=0, this reduces to
        a standard isotropic kernel.
    scaling_type : str
        Scaling function type: "linear" or "exponential".
    coord_to_distance : dict
        Mapping from (x, y) coordinates (as tuple) to distance from origin.
        Used internally to look up distances for kernel evaluation.
    """

    lengthscale: float = 0.15
    base_variance: float = 1.0
    distance_scaling: float = 0.0
    scaling_type: str = "linear"
    coord_to_distance: Optional[Dict[Tuple[float, float], float]] = None

    def __post_init__(self):
        """Initialize coordinate-to-distance mapping if needed."""
        if self.coord_to_distance is None:
            self.coord_to_distance = {}

    @classmethod
    def from_coordinates_and_distances(
        cls,
        coords: np.ndarray,
        distances: np.ndarray,
        lengthscale: float = 0.15,
        base_variance: float = 1.0,
        distance_scaling: float = 0.0,
        scaling_type: str = "linear",
    ):
        """Factory method to create kernel with distance lookup table.

        Parameters
        ----------
        coords : np.ndarray
            Coordinates of all points, shape (n, 2).
        distances : np.ndarray
            Distance from origin for each point, shape (n,).
        lengthscale, base_variance, distance_scaling, scaling_type :
            Kernel parameters.

        Returns
        -------
        DistanceDependentNNGPKernel
            Kernel instance with pre-built coordinate-to-distance mapping.
        """
        coords = np.atleast_2d(coords)
        distances = np.asarray(distances, dtype=float).ravel()
        if coords.shape[0] != distances.size:
            raise ValueError("coords and distances must have same length")

        # Build lookup dictionary
        coord_to_distance = {
            (float(coords[i, 0]), float(coords[i, 1])): float(distances[i])
            for i in range(coords.shape[0])
        }

        return cls(
            lengthscale=lengthscale,
            base_variance=base_variance,
            distance_scaling=distance_scaling,
            scaling_type=scaling_type,
            coord_to_distance=coord_to_distance,
        )

    def _get_distances(self, X: np.ndarray) -> np.ndarray:
        """Look up distances for given coordinates."""
        X = np.atleast_2d(X)
        distances = np.zeros(X.shape[0], dtype=float)
        for i in range(X.shape[0]):
            key = (float(X[i, 0]), float(X[i, 1]))
            if key in self.coord_to_distance:
                distances[i] = self.coord_to_distance[key]
            else:
                # If coordinate not in lookup, return 0 (fallback)
                distances[i] = 0.0
        return distances

    def local_std(self, distances: np.ndarray) -> np.ndarray:
        """Compute local standard deviation based on distance from origin.

        Parameters
        ----------
        distances : np.ndarray
            Distance from origin (Z-scores or normalized distances), shape (n,).

        Returns
        -------
        np.ndarray
            Local standard deviation, shape (n,).
        """
        sigma_0 = np.sqrt(self.base_variance)
        if self.distance_scaling == 0.0:
            return np.full_like(distances, sigma_0, dtype=float)

        if self.scaling_type == "linear":
            return sigma_0 * (1.0 + self.distance_scaling * distances)
        elif self.scaling_type == "exponential":
            return sigma_0 * np.exp(self.distance_scaling * distances)
        else:
            raise ValueError(f"Unknown scaling_type: {self.scaling_type}")

    def K(self, X1: np.ndarray, X2: np.ndarray) -> np.ndarray:
        """Compute distance-dependent variance kernel matrix.

        Compatible with standard kernel interface K(X1, X2).
        Distances are looked up from internal coordinate-to-distance mapping.

        Parameters
        ----------
        X1 : np.ndarray
            First set of coordinates, shape (n1, 2).
        X2 : np.ndarray
            Second set of coordinates, shape (n2, 2).

        Returns
        -------
        np.ndarray
            Kernel matrix, shape (n1, n2).
        """
        # Look up distances
        distances1 = self._get_distances(X1)
        distances2 = self._get_distances(X2)

        # Base RBF kernel (variance=1 normalization)
        X1 = np.atleast_2d(X1)
        X2 = np.atleast_2d(X2)
        a2 = np.sum(X1**2, axis=1, keepdims=True)
        b2 = np.sum(X2**2, axis=1, keepdims=True).T
        sqd = a2 + b2 - 2.0 * (X1 @ X2.T)
        K_base = np.exp(-0.5 * sqd / (self.lengthscale**2))

        # Local standard deviation
        sigma1 = self.local_std(distances1)  # (n1,)
        sigma2 = self.local_std(distances2)  # (n2,)

        # Apply distance-dependent variance
        sigma_outer = sigma1[:, np.newaxis] * sigma2[np.newaxis, :]  # (n1, n2)
        return K_base * sigma_outer


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

    For numerical stability, we subtract the maximum across all categories
    (including the baseline eta=0) at each location.
    """
    # Find max across all categories (including baseline=0)
    # Shape: (1, n)
    eta_max = np.maximum(eta.max(axis=0, keepdims=True), 0.0)

    # Subtract max for stability
    eta_stable = eta - eta_max  # (K-1, n)
    baseline_stable = 0.0 - eta_max  # (1, n), baseline eta is 0

    # Compute exponentials
    exp_eta = np.exp(eta_stable)  # (K-1, n)
    exp_baseline = np.exp(baseline_stable)  # (1, n)

    # Normalize
    sum_exp = exp_eta.sum(axis=0, keepdims=True) + exp_baseline  # (1, n)
    probs_non_baseline = exp_eta / sum_exp  # (K-1, n)
    probs_baseline = exp_baseline / sum_exp  # (1, n)

    return np.vstack([probs_non_baseline, probs_baseline])


def logratio_to_probs(logratio: np.ndarray) -> np.ndarray:
    """Convert log-ratios to probabilities.

    Given log-ratios g_k = log(p_k) - log(p_K) for k=1,...,K-1 where K is
    the baseline category, this function recovers the original probabilities
    p_1, ..., p_K that sum to 1.

    This is the mathematically correct inverse of the log-ratio transformation
    and does NOT use the max-subtraction trick (which would change the baseline).

    Parameters
    ----------
    logratio
        Log-ratio array with shape (K-1, n). Each column contains
        g_k = log(p_k / p_K) for the non-baseline categories.

    Returns
    -------
    np.ndarray
        Probability array with shape (K, n). Each column sums to 1.
        The last row contains the baseline category probabilities.

    Notes
    -----
    Mathematical derivation:
        Given g_k = log(p_k / p_K), we have exp(g_k) = p_k / p_K.
        Using the constraint sum_k p_k = 1:
            p_K * (1 + sum_{k=1}^{K-1} exp(g_k)) = 1
            => p_K = 1 / (1 + sum_{k=1}^{K-1} exp(g_k))
            => p_k = exp(g_k) * p_K for k=1,...,K-1
    """
    exp_logratio = np.exp(logratio)  # (K-1, n), each entry is p_k / p_K
    sum_exp = exp_logratio.sum(axis=0, keepdims=True)  # (1, n)
    prob_baseline = 1.0 / (1.0 + sum_exp)  # (1, n), this is p_K
    probs_non_baseline = exp_logratio * prob_baseline  # (K-1, n), these are p_1, ..., p_{K-1}
    return np.vstack([probs_non_baseline, prob_baseline])


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


def weighted_inverse_softmax(
    Z: np.ndarray,
    w: np.ndarray,
    tau: float = 1.0,
    alpha: float = 1.0,
) -> np.ndarray:
    """Compute distance-based baseline probabilities via weighted inverse softmax.

    Parameters
    ----------
    Z
        Distance Z-scores with shape (n, K). Smaller values indicate closer proximity.
    w
        Source importance weights with shape (K,). All entries must be positive.
    tau
        Temperature parameter controlling concentration. Smaller values yield
        sharper distributions around the nearest source.
    alpha
        Importance weight exponent. Setting alpha=0 ignores the weights entirely.

    Returns
    -------
    np.ndarray
        Baseline probabilities with shape (n, K). Each row sums to one.
    """
    if Z.ndim != 2:
        raise ValueError("Z must be a 2D array with shape (n, K)")
    if w.ndim != 1:
        raise ValueError("w must be a 1D array with shape (K,)")
    if Z.shape[1] != w.shape[0]:
        raise ValueError("Z.shape[1] must equal w.shape[0]")
    if not np.all(w > 0):
        raise ValueError("all entries in w must be positive")
    if tau <= 0:
        raise ValueError("tau must be positive")

    logw = np.log(w)
    bias = alpha * logw
    logits = (-Z / tau) + bias
    return np.exp(logits - logsumexp(logits, axis=1, keepdims=True))


def compute_distance_features(
    Z: np.ndarray,
    w: np.ndarray,
    tau: float = 1.0,
    alpha: float = 1.0,
) -> np.ndarray:
    """Compute log-ratio distance features from Z-scores.

    Parameters
    ----------
    Z
        Distance Z-scores with shape (n, K).
    w
        Source importance weights with shape (K,).
    tau
        Temperature parameter.
    alpha
        Importance weight exponent.

    Returns
    -------
    np.ndarray
        Log-ratio features with shape (n, K-1). These are defined as
        g_k = log(p_0k) - log(p_0K) where p_0 is the weighted inverse softmax
        of the distance Z-scores.
    """
    p0 = weighted_inverse_softmax(Z, w, tau, alpha)
    log_p0 = np.log(p0 + 1e-12)
    log_baseline = log_p0[:, -1:]
    g = log_p0[:, :-1] - log_baseline
    return g


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
    prior_mean_by_feature: Optional[Sequence[Optional[np.ndarray]]] = None,
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
    factors_by_feature
        Vecchia/NNGP factors for the observed locations, one per feature.
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
    prior_mean_by_feature
        Optional sequence of prior means for each feature (length p+1).
        If provided, prior_mean_by_feature[j] should be None (for zero mean)
        or a 1D array of shape (n,) containing the prior mean for feature j.
        Used for hierarchical GP priors with non-zero mean (e.g., distance prior).
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

        # Get prior mean for this feature (None means zero mean)
        prior_mean_j = None
        if prior_mean_by_feature is not None and j < len(prior_mean_by_feature):
            prior_mean_j = prior_mean_by_feature[j]

        for idx in order:
            neighbours = factors.neighbor_idx[idx]

            # Compute NNGP conditional prior mean
            if neighbours.size:
                # Conditional mean: mu_j(s_i) + a_i^T (beta_j(neighbors) - mu_j(neighbors))
                neighbor_values = beta_new[j, neighbours]
                if prior_mean_j is not None:
                    # Non-zero prior mean case
                    neighbor_mean = prior_mean_j[neighbours]
                    mu_prior = float(prior_mean_j[idx] +
                                   factors.a_rows[idx] @ (neighbor_values - neighbor_mean))
                else:
                    # Zero prior mean case (standard)
                    mu_prior = float(factors.a_rows[idx] @ neighbor_values)
            else:
                # No neighbors: use marginal prior mean
                mu_prior = 0.0 if prior_mean_j is None else float(prior_mean_j[idx])

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
