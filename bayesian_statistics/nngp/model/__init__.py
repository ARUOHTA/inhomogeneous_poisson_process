# Subpackage for modeling utilities (grid sampling, local NNGP, etc.)

from .multinomial_model import (
    MultinomialDataset,
    MultinomialNNGPConfig,
    MultinomialNNGPResults,
    prepare_multinomial_dataset,
    run_mcmc,
)
from .nngp import FactorCache, NNGPFactors, build_cross_factors, build_nngp_factors
from .sample import (
    LocalNNGPKernel,
    PolyaGammaSampler,
    compute_eta,
    softmax_with_baseline,
)

__all__ = [
    "MultinomialDataset",
    "MultinomialNNGPConfig",
    "MultinomialNNGPResults",
    "prepare_multinomial_dataset",
    "run_mcmc",
    "FactorCache",
    "NNGPFactors",
    "build_cross_factors",
    "build_nngp_factors",
    "LocalNNGPKernel",
    "PolyaGammaSampler",
    "compute_eta",
    "softmax_with_baseline",
]
