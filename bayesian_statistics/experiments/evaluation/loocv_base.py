"""Base LOOCV evaluation module.

Provides subsample LOOCV evaluation with Aitchison distance.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np

from bayesian_statistics.models.evaluation.metrics import CompositionalMetrics


@dataclass
class LOOCVResult:
    """Single LOOCV result for one site."""

    period: int
    site_id: int
    observed: np.ndarray
    predicted: np.ndarray
    aitchison_distance: float


class SubsampleLOOCVEvaluator:
    """Subsample LOOCV evaluator using Aitchison distance.

    Attributes
    ----------
    n_samples : int
        Number of sites to subsample for LOOCV.
    seed : int
        Random seed for reproducibility.
    """

    def __init__(self, n_samples: int = 20, seed: int = 42):
        """Initialize evaluator.

        Parameters
        ----------
        n_samples : int
            Number of sites to subsample for LOOCV.
        seed : int
            Random seed for reproducibility.
        """
        self.n_samples = n_samples
        self.seed = seed
        self.metrics = CompositionalMetrics()
        self._random = random.Random(seed)

    def _reset_random_state(self) -> None:
        """Reset random state for reproducibility."""
        self._random = random.Random(self.seed)

    def _get_subsample(self, site_ids: List[int]) -> List[int]:
        """Get subsample of site IDs.

        Parameters
        ----------
        site_ids : list of int
            All available site IDs.

        Returns
        -------
        list of int
            Subsampled site IDs.
        """
        if len(site_ids) <= self.n_samples:
            return list(site_ids)
        return self._random.sample(site_ids, self.n_samples)

    def _compute_mean_distance(self, results: List[LOOCVResult]) -> float:
        """Compute mean Aitchison distance.

        Parameters
        ----------
        results : list of LOOCVResult
            LOOCV results.

        Returns
        -------
        float
            Mean Aitchison distance.
        """
        if not results:
            return np.nan
        distances = [r.aitchison_distance for r in results]
        return np.mean(distances)

    def _compute_std_distance(self, results: List[LOOCVResult]) -> float:
        """Compute standard deviation of Aitchison distance.

        Parameters
        ----------
        results : list of LOOCVResult
            LOOCV results.

        Returns
        -------
        float
            Standard deviation of Aitchison distance.
        """
        if not results:
            return np.nan
        distances = [r.aitchison_distance for r in results]
        return np.std(distances)

    def compute_summary(
        self, results: List[LOOCVResult], period: Optional[int] = None
    ) -> Dict[str, Any]:
        """Compute summary statistics for LOOCV results.

        Parameters
        ----------
        results : list of LOOCVResult
            LOOCV results.
        period : int, optional
            Period index. If None, inferred from results.

        Returns
        -------
        dict
            Summary statistics.
        """
        if not results:
            return {
                "mean_aitchison_distance": np.nan,
                "std_aitchison_distance": np.nan,
                "n_samples": 0,
                "period": period,
            }

        # Infer period from results if not provided
        if period is None and results:
            period = results[0].period

        return {
            "mean_aitchison_distance": self._compute_mean_distance(results),
            "std_aitchison_distance": self._compute_std_distance(results),
            "n_samples": len(results),
            "period": period,
        }

    def evaluate_single_site(
        self,
        site_id: int,
        period: int,
        observed: np.ndarray,
        predicted: np.ndarray,
    ) -> LOOCVResult:
        """Evaluate prediction for a single site.

        Parameters
        ----------
        site_id : int
            Site ID.
        period : int
            Period index.
        observed : np.ndarray
            Observed composition ratio.
        predicted : np.ndarray
            Predicted composition ratio.

        Returns
        -------
        LOOCVResult
            Evaluation result.
        """
        distance = self.metrics.aitchison_distance(observed, predicted)
        return LOOCVResult(
            period=period,
            site_id=site_id,
            observed=observed,
            predicted=predicted,
            aitchison_distance=distance,
        )
