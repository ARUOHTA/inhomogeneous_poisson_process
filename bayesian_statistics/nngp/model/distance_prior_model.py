"""Distance-prior multinomial NNGP model with hierarchical GP formulation.

This module implements a multinomial NNGP model where distance information to
source locations is incorporated as a **prior mean** for the intercept GP,
rather than as an additive term. This hierarchical formulation ensures that:

1. In data-rich areas: GP posterior adapts from the distance prior based on data
2. In data-sparse areas: GP posterior stays close to the distance prior
3. Interpretation: "distance effect as baseline, GP adjusts locally"

This solves the problem where additive models dilute the distance prior in
unobserved locations.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from tqdm import tqdm

from bayesian_statistics.models.preprocessing.data_preprocessor import (
    ObsidianDataPreprocessor,
)

from .multinomial_model import (
    MultinomialDataset,
    MultinomialNNGPConfig,
    _build_kernel_list,
)
from .nngp import FactorCache
from .sample import (
    LocalNNGPKernel,
    PolyaGammaSampler,
    compute_distance_features,
    compute_eta,
    compute_log_rest_terms,
    logratio_to_probs,
    softmax_with_baseline,
    update_beta_category,
)


@dataclass
class DistancePriorConfig(MultinomialNNGPConfig):
    """Configuration for distance-prior multinomial NNGP model.

    Extends MultinomialNNGPConfig with distance prior hyperparameters.
    """

    # Distance prior hyperparameters
    tau: float = 1.0  # Temperature for weighted inverse softmax
    alpha: float = 1.0  # Weight on source importance
    source_weights: Optional[Sequence[float]] = None  # Importance weights (K,)
    lambda_fixed: Optional[Sequence[float]] = None  # Fixed scaling (K-1,) or None

    # Per-feature prior means for covariates
    prior_mean_by_feature: Optional[Sequence[Optional[float]]] = None
    # Length must be p+1 (intercept + covariates)
    # prior_mean_by_feature[0] must be None (intercept uses distance prior)
    # prior_mean_by_feature[j] for j>=1: prior mean for covariate j-1


@dataclass
class DistancePriorDataset:
    """Dataset with distance features for prior mean.

    Contains all fields from MultinomialDataset plus distance-based prior means.
    """

    # Base fields from MultinomialDataset
    period: int
    origins: list[str]
    coords: np.ndarray
    counts: np.ndarray
    total_counts: np.ndarray
    design_matrix_sites: np.ndarray
    grid_points: np.ndarray
    design_matrix_grid: np.ndarray
    site_ids: np.ndarray
    variable_names: list[str]

    # Distance Z-scores and features
    distance_zscores_sites: np.ndarray  # (n_sites, K) Z-scores
    distance_zscores_grid: np.ndarray  # (n_grid, K) Z-scores
    distance_features_sites: np.ndarray  # (n_sites, K-1) log-ratio g_k
    distance_features_grid: np.ndarray  # (n_grid, K-1) log-ratio g_k

    # Prior means for intercept (computed from distance features)
    prior_mean_intercept_sites: np.ndarray  # (K-1, n_sites) lambda * g_k
    prior_mean_intercept_grid: np.ndarray  # (K-1, n_grid) lambda * g_k

    # Hyperparameters used for distance prior
    tau: float  # Temperature parameter
    alpha: float  # Importance weight exponent
    source_weights_full: np.ndarray  # (K,) Full weights including baseline

    # Prior means for covariate features (constant across space)
    prior_mean_covariate_features: Optional[np.ndarray] = None
    # shape: (p,) where p is number of covariates (excluding intercept)
    # prior_mean_covariate_features[j-1] = prior mean for covariate j (j>=1)
    # None means all covariates have zero prior mean (backward compatible)

    def num_categories(self) -> int:
        return int(self.counts.shape[1])

    def num_sites(self) -> int:
        return int(self.coords.shape[0])


@dataclass
class DistancePriorResults:
    """Results from distance-prior MCMC sampling.

    Contains posterior samples and model configuration.
    """

    dataset: DistancePriorDataset
    config: DistancePriorConfig
    factor_cache: FactorCache
    kernels: list[LocalNNGPKernel]
    beta_samples: np.ndarray  # (n_save, K-1, p+1, n_sites)
    lambda_values: np.ndarray  # (K-1,) fixed scaling factors

    def predict_probabilities(
        self,
        location: str = "sites",
        sample_conditional: bool = False,
        seed: Optional[int] = None,
    ) -> np.ndarray:
        """Predict category probabilities at sites or grid.

        Parameters
        ----------
        location
            Either "sites" or "grid".
        sample_conditional
            If True and location="grid", sample from GP conditional.
            If False, use posterior mean projection.
        seed
            Random seed for conditional sampling.

        Returns
        -------
        np.ndarray
            Probability array with shape (K, n_locations), where K includes baseline.
        """
        K_minus_1 = self.beta_samples.shape[1]
        n_samples = self.beta_samples.shape[0]

        if location == "sites":
            W = self.dataset.design_matrix_sites
            n_loc = self.dataset.num_sites()
            beta_mean = self.beta_samples.mean(axis=0)  # (K-1, p+1, n_sites)
        elif location == "grid":
            W = self.dataset.design_matrix_grid
            n_loc = W.shape[0]
            if sample_conditional:
                # Sample beta_grid from GP conditional
                rng = np.random.default_rng(seed)
                beta_grid_samples = np.zeros((n_samples, K_minus_1, W.shape[1], n_loc))
                A_grid_list = self.factor_cache.A_grid_all()
                D_grid_list = self.factor_cache.D_grid_all()

                for sample_idx in range(n_samples):
                    beta_sites = self.beta_samples[sample_idx]  # (K-1, p+1, n_sites)
                    for k in range(K_minus_1):
                        for j in range(W.shape[1]):
                            A_grid = A_grid_list[j]
                            D_grid = D_grid_list[j]
                            beta_sites_j = beta_sites[k, j]  # (n_sites,)

                            # Conditional mean
                            mu_cond = A_grid @ beta_sites_j  # (n_grid,)

                            # Add prior mean adjustment for intercept
                            if j == 0:
                                prior_mean_grid = (
                                    self.dataset.prior_mean_intercept_grid[k]
                                )
                                prior_mean_sites = (
                                    self.dataset.prior_mean_intercept_sites[k]
                                )
                                mu_cond += prior_mean_grid - (A_grid @ prior_mean_sites)
                            # NOTE: Covariate prior means (j >= 1) are constant across space,
                            # so they cancel out in the conditional distribution:
                            # mu_cond = A_grid @ (beta_sites + c) = A_grid @ beta_sites + A_grid @ c
                            # Adjustment: c - A_grid @ c = c - c = 0 (since A_grid @ 1 = 1)
                            # Therefore, no adjustment needed for constant prior means.

                            # Sample from conditional
                            std_cond = np.sqrt(np.maximum(D_grid, 1e-12))
                            beta_grid_samples[sample_idx, k, j] = rng.normal(
                                mu_cond, std_cond
                            )

                beta_mean = beta_grid_samples.mean(axis=0)  # (K-1, p+1, n_grid)
            else:
                # Use posterior mean projection
                A_grid_list = self.factor_cache.A_grid_all()
                beta_sites_mean = self.beta_samples.mean(axis=0)  # (K-1, p+1, n_sites)
                beta_mean = np.zeros((K_minus_1, W.shape[1], n_loc))
                for k in range(K_minus_1):
                    for j in range(W.shape[1]):
                        A_grid = A_grid_list[j]
                        beta_mean[k, j] = A_grid @ beta_sites_mean[k, j]

                        # Add prior mean adjustment for intercept
                        if j == 0:
                            prior_mean_grid = self.dataset.prior_mean_intercept_grid[k]
                            prior_mean_sites = self.dataset.prior_mean_intercept_sites[
                                k
                            ]
                            beta_mean[k, j] += prior_mean_grid - (
                                A_grid @ prior_mean_sites
                            )
        else:
            raise ValueError(f"location must be 'sites' or 'grid', got {location}")

        # Compute eta and probabilities
        eta = compute_eta(beta_mean, W)  # (K-1, n_loc)
        probs = softmax_with_baseline(eta)  # (K, n_loc)
        return probs

    def decompose_effects(
        self,
        location: str = "sites",
    ) -> dict[str, np.ndarray]:
        """Decompose posterior probabilities into component effects.

        This decomposition separates the distance prior from data-driven adjustments.
        Note: effects are non-additive due to softmax nonlinearity.

        Parameters
        ----------
        location
            Either "sites" or "grid".

        Returns
        -------
        dict
            Dictionary with keys:
            - "distance": distance prior effect only (lambda_k * g_k)
            - "intercept_adjustment": data-driven adjustment to distance prior (delta_k)
            - "intercept": full intercept (distance + adjustment, for reference)
            - "covariate_{name}": individual covariate effects
            - "covariates_only": all covariates excluding intercept
            - "full": complete model (distance + adjustment + covariates)

        Notes
        -----
        The hierarchical model structure is:
            beta_0k(s) = lambda_k * g_k(s) + delta_k(s)
        where:
            - lambda_k * g_k(s) is the distance prior (fixed baseline)
            - delta_k(s) is the data-driven adjustment (learned from observations)

        The decomposition allows us to see:
            - "distance": pure distance effect (no data influence)
            - "intercept_adjustment": how data pulls away from distance baseline
            - "intercept": combined effect (should match "distance" where data is sparse)
        """
        K_minus_1 = self.beta_samples.shape[1]
        beta_mean = self.beta_samples.mean(axis=0)  # (K-1, p+1, n_sites)

        if location == "sites":
            W = self.dataset.design_matrix_sites
            n_loc = self.dataset.num_sites()
            g = self.dataset.distance_features_sites.T  # (K-1, n_sites)
        elif location == "grid":
            W = self.dataset.design_matrix_grid
            n_loc = W.shape[0]
            g = self.dataset.distance_features_grid.T  # (K-1, n_grid)

            # Project beta to grid
            A_grid_list = self.factor_cache.A_grid_all()
            beta_grid = np.zeros((K_minus_1, W.shape[1], n_loc))
            for k in range(K_minus_1):
                for j in range(W.shape[1]):
                    beta_grid[k, j] = A_grid_list[j] @ beta_mean[k, j]
                    if j == 0:
                        prior_mean_grid = self.dataset.prior_mean_intercept_grid[k]
                        prior_mean_sites = self.dataset.prior_mean_intercept_sites[k]
                        beta_grid[k, j] += prior_mean_grid - (
                            A_grid_list[j] @ prior_mean_sites
                        )
            beta_mean = beta_grid
        else:
            raise ValueError(f"location must be 'sites' or 'grid', got {location}")

        effects = {}

        # Helper function: softmax without numerical stabilization
        # Used for decompose_effects to preserve baseline category and ensure
        # consistency with logratio_to_probs
        def softmax_no_stabilization(eta: np.ndarray) -> np.ndarray:
            """Softmax with baseline, without max-subtraction."""
            exp_eta = np.exp(eta)
            sum_exp = exp_eta.sum(axis=0, keepdims=True)
            baseline = 1.0 + sum_exp
            probs_non_baseline = exp_eta / baseline
            probs_baseline = 1.0 / baseline
            return np.vstack([probs_non_baseline, probs_baseline])

        # 1. Distance prior effect only (baseline expectation)
        # g is already in log-ratio form: g_k = log(p_0k) - log(p_0K)
        # Scale by lambda and convert back to probabilities
        logratio_dist = self.lambda_values[:, np.newaxis] * g  # (K-1, n_loc)
        effects["distance"] = logratio_to_probs(logratio_dist)

        # 2. Data-driven adjustment to distance prior
        # beta_0k(s) = lambda_k * g_k(s) + delta_k(s)
        # delta_k(s) represents how data pulls away from distance prior
        eta_intercept = beta_mean[:, 0, :]  # (K-1, n_loc)
        eta_adjustment = eta_intercept - logratio_dist  # delta_k(s)
        effects["intercept_adjustment"] = softmax_no_stabilization(eta_adjustment)

        # 3. Full intercept (distance + adjustment, for reference)
        effects["intercept"] = softmax_no_stabilization(eta_intercept)

        # 4. Individual covariates
        p_plus_1 = beta_mean.shape[1]
        for j in range(1, p_plus_1):
            covar_name = self.dataset.variable_names[j - 1]
            w_j = W[:, j]
            eta_j = beta_mean[:, j, :] * w_j
            effects[f"covariate_{covar_name}"] = softmax_no_stabilization(eta_j)

        # 5. All covariates only (excluding intercept)
        if p_plus_1 > 1:
            eta_covariates_only = compute_eta(beta_mean[:, 1:, :], W[:, 1:])
            effects["covariates_only"] = softmax_no_stabilization(eta_covariates_only)

        # 6. Full model (distance + adjustment + covariates)
        eta_full = compute_eta(beta_mean, W)
        effects["full"] = softmax_no_stabilization(eta_full)

        return effects


def prepare_distance_prior_dataset(
    preprocessor: ObsidianDataPreprocessor,
    period: int,
    origins: Sequence[str],
    variable_names: Sequence[str],
    distance_column_names: Sequence[str],
    grid_subsample_ratio: float = 0.1,
    drop_zero_total_sites: bool = True,
    tau: float = 1.0,
    alpha: float = 1.0,
    source_weights: Optional[Sequence[float]] = None,
    lambda_values: Optional[Sequence[float]] = None,
    prior_mean_by_feature: Optional[Sequence[Optional[float]]] = None,
) -> DistancePriorDataset:
    """Prepare dataset with distance-based prior means.

    Parameters
    ----------
    preprocessor
        Preprocessor containing site and grid data.
    period
        Time period to extract.
    origins
        Category names including baseline (length K).
    variable_names
        Covariate names (excluding intercept).
    distance_column_names
        Column names for distance Z-scores (length K-1, excluding baseline).
    grid_subsample_ratio
        Fraction of grid points to use.
    drop_zero_total_sites
        Whether to drop sites with zero total count.
    tau
        Temperature parameter for weighted inverse softmax.
    alpha
        Importance weight exponent.
    source_weights
        Source importance weights (length K-1). If None, use uniform weights.
    lambda_values
        Scaling factors for prior mean (length K-1). If None, use 1.0 for all.
    prior_mean_by_feature
        Optional prior means for each feature (length p+1).
        prior_mean_by_feature[0] must be None (intercept uses distance prior).
        prior_mean_by_feature[j] for j>=1: prior mean for covariate j-1.
        If None, all covariates have zero prior mean (default).

    Returns
    -------
    DistancePriorDataset
        Dataset with distance features and prior means.
    """
    from .multinomial_model import prepare_multinomial_dataset

    # Prepare base dataset
    base_dataset = prepare_multinomial_dataset(
        preprocessor=preprocessor,
        period=period,
        origins=origins,
        variable_names=variable_names,
        grid_subsample_ratio=grid_subsample_ratio,
        drop_zero_total_sites=drop_zero_total_sites,
    )

    K = len(origins)
    K_minus_1 = K - 1

    # Extract distance data from preprocessor DataFrames
    site_df = preprocessor.df_sites.sort("遺跡ID")
    elevation_df = preprocessor.df_elevation.sort(["y", "x"])

    # Distance at sites (full, before filtering)
    dist_sites_full = site_df.select(distance_column_names).to_numpy().astype(float)

    # Distance at grid (full)
    dist_grid_full = elevation_df.select(distance_column_names).to_numpy().astype(float)

    # Apply same filtering as base_dataset
    # base_dataset.site_ids contains the IDs of sites that were kept
    # We need to select the corresponding rows from dist_sites_full
    # Since site_df is sorted by 遺跡ID and site_ids are the kept IDs,
    # we can use site_ids as indices
    dist_sites = dist_sites_full[base_dataset.site_ids]

    # Subsample grid to match base_dataset
    # Replicate the same subsampling logic as prepare_multinomial_dataset
    if grid_subsample_ratio < 1.0:
        n_grid = dist_grid_full.shape[0]
        keep = max(int(np.floor(n_grid * grid_subsample_ratio)), 1)
        indices = np.linspace(0, n_grid - 1, keep, dtype=int)
        dist_grid = dist_grid_full[indices]
    else:
        dist_grid = dist_grid_full

    if dist_sites.shape[1] != K_minus_1:
        raise ValueError(
            f"Expected {K_minus_1} distance columns, got {dist_sites.shape[1]}"
        )

    # Compute Z-scores
    mean_dist = dist_grid_full.mean(axis=0, keepdims=True)
    std_dist = dist_grid_full.std(axis=0, keepdims=True)
    std_dist = np.where(std_dist < 1e-12, 1.0, std_dist)

    Z_sites = (dist_sites - mean_dist) / std_dist
    Z_grid = (dist_grid - mean_dist) / std_dist

    # Add dummy distance for baseline category (high Z-score = far away)
    baseline_z_sites = np.full((Z_sites.shape[0], 1), 3.0, dtype=float)
    baseline_z_grid = np.full((Z_grid.shape[0], 1), 3.0, dtype=float)
    Z_sites_with_baseline = np.column_stack([Z_sites, baseline_z_sites])
    Z_grid_with_baseline = np.column_stack([Z_grid, baseline_z_grid])

    # Set source weights
    if source_weights is None:
        w = np.ones(K, dtype=float)
    else:
        if len(source_weights) != K_minus_1:
            raise ValueError(
                f"source_weights must have length {K_minus_1}, got {len(source_weights)}"
            )
        w = np.array(list(source_weights) + [1.0], dtype=float)

    # Compute distance features (log-ratio)
    g_sites = compute_distance_features(
        Z_sites_with_baseline, w, tau=tau, alpha=alpha
    )  # (n_sites, K-1)
    g_grid = compute_distance_features(
        Z_grid_with_baseline, w, tau=tau, alpha=alpha
    )  # (n_grid, K-1)

    # Set lambda values
    if lambda_values is None:
        lambda_vec = np.ones(K_minus_1, dtype=float)
    else:
        if len(lambda_values) != K_minus_1:
            raise ValueError(
                f"lambda_values must have length {K_minus_1}, got {len(lambda_values)}"
            )
        lambda_vec = np.array(lambda_values, dtype=float)

    # Compute prior means for intercept: mu_k(s) = lambda_k * g_k(s)
    prior_mean_intercept_sites = lambda_vec[:, np.newaxis] * g_sites.T  # (K-1, n_sites)
    prior_mean_intercept_grid = lambda_vec[:, np.newaxis] * g_grid.T  # (K-1, n_grid)

    # Process prior means for covariate features
    p = len(variable_names)  # Number of covariates
    p_plus_1 = p + 1  # Including intercept

    prior_mean_covariates = None
    if prior_mean_by_feature is not None:
        # Validate length
        if len(prior_mean_by_feature) != p_plus_1:
            raise ValueError(
                f"prior_mean_by_feature must have length {p_plus_1} "
                f"(1 intercept + {p} covariates), got {len(prior_mean_by_feature)}"
            )

        # Validate intercept is None
        if prior_mean_by_feature[0] is not None:
            raise ValueError(
                "prior_mean_by_feature[0] must be None "
                "(intercept uses distance prior, not a constant prior mean)"
            )

        # Extract covariate prior means (skip intercept at index 0)
        covariate_means = [
            0.0 if x is None else float(x) for x in prior_mean_by_feature[1:]
        ]
        prior_mean_covariates = np.array(covariate_means, dtype=float)

    return DistancePriorDataset(
        # Base fields from MultinomialDataset
        period=base_dataset.period,
        origins=base_dataset.origins,
        coords=base_dataset.coords,
        counts=base_dataset.counts,
        total_counts=base_dataset.total_counts,
        design_matrix_sites=base_dataset.design_matrix_sites,
        grid_points=base_dataset.grid_points,
        design_matrix_grid=base_dataset.design_matrix_grid,
        site_ids=base_dataset.site_ids,
        variable_names=base_dataset.variable_names,
        # Distance fields
        distance_zscores_sites=Z_sites_with_baseline,
        distance_zscores_grid=Z_grid_with_baseline,
        distance_features_sites=g_sites,
        distance_features_grid=g_grid,
        prior_mean_intercept_sites=prior_mean_intercept_sites,
        prior_mean_intercept_grid=prior_mean_intercept_grid,
        # Hyperparameters
        tau=tau,
        alpha=alpha,
        source_weights_full=w,  # Save the full weight array (K,) for reuse
        # Covariate prior means
        prior_mean_covariate_features=prior_mean_covariates,
    )


def _build_prior_means_for_category(
    category_idx: int,
    dataset: DistancePriorDataset,
    n_sites: int,
) -> list[Optional[np.ndarray]]:
    """Build prior mean arrays for each feature of a given category.

    This helper function constructs the prior mean specification for all features
    (intercept + covariates) of a specific category, to be passed to update_beta_category.

    Parameters
    ----------
    category_idx
        Category index (0 to K-2).
    dataset
        Dataset containing prior mean information.
    n_sites
        Number of observation sites.

    Returns
    -------
    list[Optional[np.ndarray]]
        List of length p+1, where:
        - Element 0: prior mean for intercept (shape (n_sites,), from distance prior)
        - Element j (j>=1): prior mean for covariate j-1
          - If prior_mean_covariate_features is None or value is 0: returns None (zero mean)
          - Otherwise: returns constant array of shape (n_sites,) with the prior mean value

    Examples
    --------
    >>> # Intercept uses distance prior, covariates use zero mean
    >>> prior_means = _build_prior_means_for_category(0, dataset, 100)
    >>> prior_means[0].shape  # (100,) - distance prior
    >>> prior_means[1]  # None - zero mean for covariate 1

    >>> # With non-zero covariate prior mean
    >>> dataset.prior_mean_covariate_features = np.array([0.5, -0.3])
    >>> prior_means = _build_prior_means_for_category(0, dataset, 100)
    >>> prior_means[1]  # array([0.5, 0.5, ..., 0.5]) - constant prior mean
    >>> prior_means[2]  # array([-0.3, -0.3, ..., -0.3])
    """
    p_plus_1 = dataset.design_matrix_sites.shape[1]
    prior_means = []

    for j in range(p_plus_1):
        if j == 0:
            # Intercept: use distance-based prior mean (spatial)
            prior_means.append(dataset.prior_mean_intercept_sites[category_idx])
        else:
            # Covariate j-1 (since j=0 is intercept)
            covariate_idx = j - 1

            if dataset.prior_mean_covariate_features is None:
                # No prior means specified: use zero mean
                prior_means.append(None)
            else:
                mean_value = dataset.prior_mean_covariate_features[covariate_idx]
                if mean_value == 0.0:
                    # Explicit zero mean
                    prior_means.append(None)
                else:
                    # Non-zero constant prior mean (replicated across all sites)
                    prior_means.append(np.full(n_sites, mean_value, dtype=float))

    return prior_means


def run_mcmc_with_distance_prior(
    dataset: DistancePriorDataset,
    config: DistancePriorConfig,
    kernel: Optional[LocalNNGPKernel | Sequence[LocalNNGPKernel]] = None,
) -> DistancePriorResults:
    """Run MCMC sampling with hierarchical distance prior.

    In this formulation, the distance effect is incorporated as the prior mean
    of the intercept GP, not as an additive term. This ensures that:
    - Data-rich areas: GP posterior adapts from distance prior
    - Data-sparse areas: GP posterior stays near distance prior

    Parameters
    ----------
    dataset
        Prepared dataset with distance prior means.
    config
        MCMC configuration with distance prior hyperparameters.
    kernel
        Optional kernel specification. Can be:
        - None: use config's kernel parameters (default)
        - LocalNNGPKernel: single kernel for all features
        - Sequence[LocalNNGPKernel]: one kernel per feature (length p+1)

    Returns
    -------
    DistancePriorResults
        Posterior samples and model configuration.
    """
    rng = np.random.default_rng(config.seed)
    K_minus_1 = dataset.num_categories() - 1
    n_sites = dataset.num_sites()
    p_plus_1 = dataset.design_matrix_sites.shape[1]

    # Set lambda values
    if config.lambda_fixed is not None:
        lambda_vec = np.array(config.lambda_fixed, dtype=float)
        if lambda_vec.shape[0] != K_minus_1:
            raise ValueError(
                f"lambda_fixed must have length {K_minus_1}, got {lambda_vec.shape[0]}"
            )
    else:
        lambda_vec = np.ones(K_minus_1, dtype=float)

    # Build kernels using the shared kernel construction logic
    # This properly handles kernel_lengthscale_by_feature and kernel_variance_by_feature
    kernels = _build_kernel_list(config, p_plus_1, kernel)

    # Build NNGP factors using FactorCache
    print("Building NNGP factors...")
    factor_cache = FactorCache(
        s_coords=dataset.coords,
        grid_points=dataset.grid_points,
        M=config.neighbor_count,
        kernels=kernels,
    )

    # Get ordering and factors
    order = factor_cache.order_S
    factors_S_list = factor_cache.factors_S_all()

    # Initialize
    print("Initializing parameters...")
    W_sites = dataset.design_matrix_sites
    totals = dataset.total_counts
    counts = dataset.counts.T  # (K-1, n_sites)

    # Initialize beta with prior mean for intercept
    beta = np.zeros((K_minus_1, p_plus_1, n_sites), dtype=float)
    for k in range(K_minus_1):
        beta[k, 0, :] = dataset.prior_mean_intercept_sites[k]  # Intercept
        # Covariates stay at 0

    eta = compute_eta(beta, W_sites)  # (K-1, n_sites)

    # Pólya-Gamma sampler
    pg_sampler = PolyaGammaSampler()

    # Storage
    n_save = (config.n_iter - config.burn_in) // config.thinning
    beta_samples = np.zeros((n_save, K_minus_1, p_plus_1, n_sites), dtype=float)

    # MCMC loop
    print(f"Running MCMC for {config.n_iter} iterations...")
    save_idx = 0
    for iteration in tqdm(range(config.n_iter)):
        # 1. Sample omega (Pólya-Gamma)
        log_rest = compute_log_rest_terms(eta)  # (K-1, n_sites)
        for k in range(K_minus_1):
            C = log_rest[k]
            psi = eta[k] - C
            omega = pg_sampler.sample(totals, psi)
            kappa = counts[k] - 0.5 * totals
            kappa_tilde = kappa + omega * C

            # 2. Sample beta_k with per-feature prior means
            prior_means = _build_prior_means_for_category(k, dataset, n_sites)

            beta_k, eta_k = update_beta_category(
                beta[k],
                eta[k],
                W_sites,
                factors_S_list,
                order,
                omega,
                kappa_tilde,
                rng,
                prior_mean_by_feature=prior_means,
            )
            beta[k] = beta_k
            eta[k] = eta_k

        # Save samples
        if (
            iteration >= config.burn_in
            and (iteration - config.burn_in) % config.thinning == 0
        ):
            beta_samples[save_idx] = beta
            save_idx += 1

    print(f"MCMC complete. Saved {save_idx} samples.")

    return DistancePriorResults(
        dataset=dataset,
        config=config,
        factor_cache=factor_cache,
        kernels=kernels,
        beta_samples=beta_samples,
        lambda_values=lambda_vec,
    )


__all__ = [
    "DistancePriorConfig",
    "DistancePriorDataset",
    "DistancePriorResults",
    "prepare_distance_prior_dataset",
    "run_mcmc_with_distance_prior",
]
