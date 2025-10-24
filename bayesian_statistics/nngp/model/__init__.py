# Subpackage for modeling utilities (grid sampling, local NNGP, etc.)

from .distance_prior_model import (
    DistancePriorConfig,
    DistancePriorDataset,
    DistancePriorResults,
    prepare_distance_prior_dataset,
    run_mcmc_with_distance_prior,
)
from .multinomial_model import (
    MultinomialDataset,
    MultinomialNNGPConfig,
    MultinomialNNGPResults,
    prepare_multinomial_dataset,
    run_mcmc,
)
from .nngp import FactorCache, NNGPFactors, build_cross_factors, build_nngp_factors
from .sample import (
    DistanceDependentNNGPKernel,
    LocalNNGPKernel,
    PolyaGammaSampler,
    compute_distance_features,
    compute_eta,
    logratio_to_probs,
    softmax_with_baseline,
    weighted_inverse_softmax,
)

__all__ = [
    "MultinomialDataset",
    "MultinomialNNGPConfig",
    "MultinomialNNGPResults",
    "prepare_multinomial_dataset",
    "run_mcmc",
    "DistancePriorConfig",
    "DistancePriorDataset",
    "DistancePriorResults",
    "prepare_distance_prior_dataset",
    "run_mcmc_with_distance_prior",
    "FactorCache",
    "NNGPFactors",
    "build_cross_factors",
    "build_nngp_factors",
    "DistanceDependentNNGPKernel",
    "LocalNNGPKernel",
    "PolyaGammaSampler",
    "compute_eta",
    "logratio_to_probs",
    "softmax_with_baseline",
    "weighted_inverse_softmax",
    "compute_distance_features",
]
