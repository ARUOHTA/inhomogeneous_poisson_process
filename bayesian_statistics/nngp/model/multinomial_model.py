"""Multinomial NNGP model and MCMC driver."""

from __future__ import annotations

from collections.abc import Sequence as SequenceABC
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from tqdm import tqdm

from bayesian_statistics.models.preprocessing.data_preprocessor import (
    ObsidianDataPreprocessor,
)

from .nngp import FactorCache
from .sample import (
    LocalNNGPKernel,
    PolyaGammaSampler,
    compute_eta,
    compute_log_rest_terms,
    softmax_with_baseline,
    update_beta_category,
)


@dataclass
class MultinomialNNGPConfig:
    """Hyper-parameters controlling the Gibbs sampler."""

    n_iter: int = 2000
    burn_in: int = 500
    thinning: int = 1
    neighbor_count: int = 25
    kernel_lengthscale: float = 0.05
    kernel_variance: float = 1.0
    kernel_lengthscale_by_feature: Optional[Sequence[float]] = None
    kernel_variance_by_feature: Optional[Sequence[float]] = None
    seed: Optional[int] = None

    def n_saved(self) -> int:
        usable = max(self.n_iter - self.burn_in, 0)
        if self.thinning <= 0:
            raise ValueError("thinning must be at least 1")
        return usable // self.thinning


@dataclass
class MultinomialDataset:
    """Container for preprocessed inputs."""

    period: int
    origins: List[str]
    coords: np.ndarray
    counts: np.ndarray
    total_counts: np.ndarray
    design_matrix_sites: np.ndarray
    grid_points: np.ndarray
    design_matrix_grid: np.ndarray
    site_ids: np.ndarray
    variable_names: List[str]

    def num_categories(self) -> int:
        return int(self.counts.shape[1])

    def num_sites(self) -> int:
        return int(self.coords.shape[0])


@dataclass
class MultinomialNNGPResults:
    """Posterior samples and helpers for summarising the model."""

    dataset: MultinomialDataset
    config: MultinomialNNGPConfig
    factor_cache: FactorCache
    kernels: List[LocalNNGPKernel]
    beta_samples: np.ndarray  # shape (n_save, K-1, p+1, n_sites)

    def num_draws(self) -> int:
        return int(self.beta_samples.shape[0])

    def num_categories(self) -> int:
        return int(self.dataset.num_categories())

    def posterior_site_mean(self) -> np.ndarray:
        """Posterior mean of site-level composition (shape: K × n_sites)."""
        W = self.dataset.design_matrix_sites
        accumulator = np.zeros(
            (self.num_categories(), self.dataset.num_sites()), dtype=float
        )
        for beta in self.beta_samples:
            eta = compute_eta(beta, W)
            probs = softmax_with_baseline(eta)
            accumulator += probs
        return accumulator / max(self.num_draws(), 1)

    def posterior_grid_mean(
        self,
        sample_conditional: bool = False,
        seed: Optional[int] = None,
    ) -> np.ndarray:
        """Posterior mean of grid-level composition.

        Parameters
        ----------
        sample_conditional
            If ``True`` draw grid coefficients from their conditional
            distribution instead of using only the conditional mean. This is a
            cheap way to include Vecchia uncertainty at prediction time.
        seed
            Optional random seed used when ``sample_conditional`` is ``True``.
        """
        if self.dataset.design_matrix_grid.size == 0:
            raise ValueError("grid design matrix is empty")

        rng = np.random.default_rng(seed) if sample_conditional else None
        W_grid = self.dataset.design_matrix_grid
        n_grid = W_grid.shape[0]
        K_minus_1 = self.beta_samples.shape[1]
        accumulator = np.zeros((self.num_categories(), n_grid), dtype=float)
        A_grid_list = self.factor_cache.A_grid_all()
        p_plus_1 = W_grid.shape[1]
        if len(A_grid_list) != p_plus_1:
            raise ValueError(
                "design matrix feature count does not match cached kernels"
            )
        sqrt_d_grid_list = [
            self.factor_cache.sqrt_d_grid_for_feature(j) for j in range(p_plus_1)
        ]

        for beta in self.beta_samples:
            eta_grid = np.zeros((K_minus_1, n_grid), dtype=float)
            for k in range(K_minus_1):
                component = beta[k]  # (p+1, n_sites)
                eta_vec = np.zeros(n_grid, dtype=float)
                for j in range(p_plus_1):
                    proj = A_grid_list[j] @ component[j]
                    if sample_conditional:
                        noise = rng.normal(size=n_grid) * sqrt_d_grid_list[j]
                        proj = proj + noise
                    eta_vec += W_grid[:, j] * proj
                eta_grid[k] = eta_vec
            probs = softmax_with_baseline(eta_grid)
            accumulator += probs
        return accumulator / max(self.num_draws(), 1)


def _prepare_counts(
    preprocessor: ObsidianDataPreprocessor,
    period: int,
    origins: Sequence[str],
) -> Tuple[np.ndarray, np.ndarray]:
    if len(origins) < 2:
        raise ValueError(
            "origins must contain at least one non-baseline category and a baseline"
        )

    counts_non_baseline: List[np.ndarray] = []
    totals: Optional[np.ndarray] = None
    for origin in origins[:-1]:
        totals, target_counts = preprocessor.preprocess_obsidian_data(period, origin)
        counts_non_baseline.append(target_counts.astype(float))
    if totals is None:
        raise RuntimeError("preprocess_obsidian_data returned no counts")
    counts_mat = np.stack(counts_non_baseline, axis=1)
    baseline = np.maximum(totals.astype(float) - counts_mat.sum(axis=1), 0.0)
    counts = np.column_stack([counts_mat, baseline])
    return counts, counts.sum(axis=1)


def prepare_multinomial_dataset(
    preprocessor: ObsidianDataPreprocessor,
    period: int,
    origins: Sequence[str],
    variable_names: Sequence[str],
    grid_subsample_ratio: float = 1.0,
    drop_zero_total_sites: bool = False,
) -> MultinomialDataset:
    """Collect all arrays required by the multinomial NNGP sampler."""
    site_df = preprocessor.df_sites.sort("遺跡ID")
    coords = site_df.select(["経度", "緯度"]).to_numpy().astype(float)
    site_ids = site_df.select("遺跡ID").to_numpy().astype(int).ravel()

    counts, totals = _prepare_counts(preprocessor, period, origins)
    if counts.shape[0] != site_ids.shape[0]:
        raise ValueError(
            "Mismatch between number of sites in counts and metadata; "
            "ensure site identifiers are contiguous and start at zero."
        )
    if drop_zero_total_sites:
        mask = totals > 0
        if not np.any(mask):
            raise ValueError("drop_zero_total_sites left no observations")
        coords = coords[mask]
        counts = counts[mask]
        totals = totals[mask]
        site_ids = site_ids[mask]

    W_grids, W_sites = preprocessor.create_explanatory_variables(list(variable_names))
    if drop_zero_total_sites:
        W_sites = W_sites[mask]
    design_matrix_sites = np.column_stack([np.ones(W_sites.shape[0]), W_sites]).astype(
        float
    )

    elevation_df = preprocessor.df_elevation.sort(["y", "x"])
    grid_points = elevation_df.select(["x", "y"]).to_numpy().astype(float)
    design_matrix_grid = np.column_stack([np.ones(W_grids.shape[0]), W_grids]).astype(
        float
    )

    if not (0 < grid_subsample_ratio <= 1.0):
        raise ValueError("grid_subsample_ratio must be in (0, 1]")
    if grid_subsample_ratio < 1.0:
        n_grid = grid_points.shape[0]
        keep = max(int(np.floor(n_grid * grid_subsample_ratio)), 1)
        indices = np.linspace(0, n_grid - 1, keep, dtype=int)
        grid_points = grid_points[indices]
        design_matrix_grid = design_matrix_grid[indices]

    return MultinomialDataset(
        period=period,
        origins=list(origins),
        coords=coords,
        counts=counts,
        total_counts=totals.astype(float),
        design_matrix_sites=design_matrix_sites,
        grid_points=grid_points,
        design_matrix_grid=design_matrix_grid,
        site_ids=site_ids,
        variable_names=list(variable_names),
    )


def _expand_feature_parameter(
    values: Optional[Sequence[float]],
    default_value: float,
    size: int,
    name: str,
) -> List[float]:
    if values is None:
        return [float(default_value)] * size
    expanded = list(values)
    if len(expanded) != size:
        raise ValueError(f"{name} must have length {size}, got {len(expanded)}")
    return [float(v) for v in expanded]


def _build_kernel_list(
    config: MultinomialNNGPConfig,
    n_features: int,
    kernel_override: Optional[LocalNNGPKernel | Sequence[LocalNNGPKernel]],
) -> List[LocalNNGPKernel]:
    if kernel_override is not None:
        if isinstance(kernel_override, LocalNNGPKernel):
            return [kernel_override] * n_features
        if isinstance(kernel_override, SequenceABC):
            kernels = list(kernel_override)
            if len(kernels) != n_features:
                raise ValueError(
                    f"kernel override length must match number of features ({n_features})"
                )
            return [
                LocalNNGPKernel(lengthscale=k.lengthscale, variance=k.variance)
                for k in kernels
            ]
        raise TypeError("kernel must be a LocalNNGPKernel or a sequence of them")

    lengthscales = _expand_feature_parameter(
        config.kernel_lengthscale_by_feature,
        config.kernel_lengthscale,
        n_features,
        "kernel_lengthscale_by_feature",
    )
    variances = _expand_feature_parameter(
        config.kernel_variance_by_feature,
        config.kernel_variance,
        n_features,
        "kernel_variance_by_feature",
    )
    if (
        config.kernel_lengthscale_by_feature is None
        and config.kernel_variance_by_feature is None
    ):
        base_kernel = LocalNNGPKernel(
            lengthscale=config.kernel_lengthscale,
            variance=config.kernel_variance,
        )
        return [base_kernel] * n_features
    return [
        LocalNNGPKernel(lengthscale=ls, variance=var)
        for ls, var in zip(lengthscales, variances)
    ]


def _initialise_beta(dataset: MultinomialDataset) -> np.ndarray:
    K_minus_1 = dataset.num_categories() - 1
    p_plus_1 = dataset.design_matrix_sites.shape[1]
    n_sites = dataset.num_sites()
    beta = np.zeros((K_minus_1, p_plus_1, n_sites), dtype=float)

    counts_sum = dataset.counts.sum(axis=0)
    baseline = counts_sum[-1]
    for k in range(K_minus_1):
        numerator = counts_sum[k] + 0.5
        denominator = baseline + 0.5
        global_logit = np.log(numerator / denominator)
        beta[k, 0, :] = global_logit
    return beta


def run_mcmc(
    dataset: MultinomialDataset,
    config: MultinomialNNGPConfig,
    kernel: Optional[LocalNNGPKernel | Sequence[LocalNNGPKernel]] = None,
) -> MultinomialNNGPResults:
    if config.n_saved() <= 0:
        raise ValueError("configuration would save zero posterior samples")

    W_sites = dataset.design_matrix_sites
    W_grid = dataset.design_matrix_grid
    counts = dataset.counts
    totals = dataset.total_counts
    K = dataset.num_categories()

    K_minus_1 = K - 1
    n_sites = dataset.num_sites()
    p_plus_1 = W_sites.shape[1]
    print(f"1. building NNGP factor cache for {n_sites} sites...")
    kernels = _build_kernel_list(config, p_plus_1, kernel)
    print(f"2. using {config.neighbor_count} neighbors for NNGP approximation...")
    factor_cache = FactorCache(
        s_coords=dataset.coords,
        grid_points=dataset.grid_points,
        M=config.neighbor_count,
        kernels=kernels,
    )
    print("3. running MCMC...")
    order = factor_cache.order_S
    factors_S_list = factor_cache.factors_S_all()

    beta = _initialise_beta(dataset)
    eta = compute_eta(beta, W_sites)

    rng = np.random.default_rng(config.seed)
    pg_sampler = PolyaGammaSampler()

    n_save = config.n_saved()
    beta_samples = np.zeros((n_save, K_minus_1, p_plus_1, n_sites), dtype=float)
    print(f"4. saving {n_save} posterior samples...")
    save_idx = 0
    for iteration in tqdm(range(config.n_iter)):
        for k in range(K_minus_1):
            log_rest = compute_log_rest_terms(eta)
            C = log_rest[k]
            psi = eta[k] - C
            omega = pg_sampler.sample(totals, psi)
            kappa = counts[:, k] - 0.5 * totals
            kappa_tilde = kappa + omega * C
            beta_k, eta_k = update_beta_category(
                beta[k],
                eta[k],
                W_sites,
                factors_S_list,
                order,
                omega,
                kappa_tilde,
                rng,
            )
            beta[k] = beta_k
            eta[k] = eta_k

        if (
            iteration >= config.burn_in
            and (iteration - config.burn_in) % config.thinning == 0
        ):
            beta_samples[save_idx] = beta
            save_idx += 1

    return MultinomialNNGPResults(
        dataset=dataset,
        config=config,
        factor_cache=factor_cache,
        kernels=kernels,
        beta_samples=beta_samples,
    )


__all__ = [
    "MultinomialNNGPConfig",
    "MultinomialDataset",
    "MultinomialNNGPResults",
    "prepare_multinomial_dataset",
    "run_mcmc",
]
