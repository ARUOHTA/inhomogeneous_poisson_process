"""NW (Nadaraya-Watson) model runner."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict, List, Optional

import numpy as np

from ..config import ExperimentConfig
from ..output import ProgressManager, get_progress_manager

if TYPE_CHECKING:
    from bayesian_statistics.models.preprocessing.data_preprocessor import (
        ObsidianDataPreprocessor,
    )


class NWRunner:
    """Runner for Nadaraya-Watson estimator.

    This class wraps the NadarayaWatsonEstimator to provide a consistent
    interface matching MMCPRunner for easy comparison.

    Attributes
    ----------
    config : ExperimentConfig
        Experiment configuration.
    sigma : float
        Kernel bandwidth for grid.
    sigma_for_sites : float
        Kernel bandwidth for sites.
    progress : ProgressManager
        Progress manager for output control.
    """

    def __init__(
        self,
        config: ExperimentConfig,
        sigma: float = 500,
        sigma_for_sites: float = 0.1,
        progress: ProgressManager | None = None,
    ):
        """Initialize NWRunner.

        Parameters
        ----------
        config : ExperimentConfig
            Experiment configuration.
        sigma : float
            Kernel bandwidth for grid estimation.
        sigma_for_sites : float
            Kernel bandwidth for site estimation.
        progress : ProgressManager, optional
            Progress manager for output control. If None, uses global default.
        """
        self.config = config
        self.sigma = sigma
        self.sigma_for_sites = sigma_for_sites
        self.progress = progress or get_progress_manager()
        self._estimator = None
        self._is_fitted = False

    def _create_estimator(self):
        """Create NadarayaWatsonEstimator instance.

        Returns
        -------
        NadarayaWatsonEstimator
            Estimator instance.
        """
        from bayesian_statistics.models.composition.nadaraya_watson import (
            NadarayaWatsonEstimator,
        )

        return NadarayaWatsonEstimator(
            sigma=self.sigma,
            sigma_for_sites=self.sigma_for_sites,
            variable_names=list(self.config.distance_column_names),
        )

    def fit(
        self,
        preprocessor: "ObsidianDataPreprocessor",
    ) -> None:
        """Fit the estimator to compute kernel weights.

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            Data preprocessor.
        """
        self._estimator = self._create_estimator()
        self._estimator.fit(preprocessor)
        self._is_fitted = True

    def _run_nw(
        self,
        preprocessor: "ObsidianDataPreprocessor",
        period: int,
    ) -> Dict[str, Any]:
        """Run NW estimation for a single period.

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            Data preprocessor.
        period : int
            Time period index.

        Returns
        -------
        dict
            Dictionary containing site_probs, grid_probs, estimator.
        """
        if not self._is_fitted:
            self.fit(preprocessor)

        # Get origins excluding "その他"
        origins = self.config.origins[:-1]  # Exclude "その他"
        n_origins = len(self.config.origins)

        # Collect results for each origin
        site_probs_list = []
        grid_probs_list = []

        for origin in origins:
            result = self._estimator.predict_single(preprocessor, period, origin)
            grid_probs_list.append(result["ratio_mesh"].ravel())
            site_probs_list.append(result["ratio_sites"])

        # Stack into arrays
        # Shape: (n_origins-1, n_points)
        site_probs_partial = np.array(site_probs_list)
        grid_probs_partial = np.array(grid_probs_list)

        # Add "その他" as complement (1 - sum of others)
        site_others = np.maximum(0, 1 - site_probs_partial.sum(axis=0))
        grid_others = np.maximum(0, 1 - grid_probs_partial.sum(axis=0))

        site_probs = np.vstack([site_probs_partial, site_others[np.newaxis, :]])
        grid_probs = np.vstack([grid_probs_partial, grid_others[np.newaxis, :]])

        return {
            "site_probs": site_probs,
            "grid_probs": grid_probs,
            "estimator": self._estimator,
        }

    def run_single_period(
        self,
        preprocessor: "ObsidianDataPreprocessor",
        period: int,
    ) -> Dict[str, Any]:
        """Run NW model for a single period.

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            Data preprocessor.
        period : int
            Time period index (0-4).

        Returns
        -------
        dict
            Dictionary containing:
            - site_probs: Predicted probabilities at sites (K, n_sites)
            - grid_probs: Predicted probabilities at grid (K, n_grid)
            - period: Time period index
            - estimator: The fitted estimator
        """
        result = self._run_nw(preprocessor, period)

        return {
            "site_probs": result["site_probs"],
            "grid_probs": result["grid_probs"],
            "period": period,
            "estimator": result["estimator"],
        }

    def run_all_periods(
        self,
        preprocessor: "ObsidianDataPreprocessor",
        periods: Optional[List[int]] = None,
    ) -> Dict[int, Dict[str, Any]]:
        """Run NW model for all specified periods.

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            Data preprocessor.
        periods : list of int, optional
            List of period indices to run. If None, runs all 5 periods.

        Returns
        -------
        dict
            Dictionary mapping period index to result dict.
        """
        if periods is None:
            periods = list(self.config.time_periods.keys())

        # Fit once for all periods
        if not self._is_fitted:
            self.fit(preprocessor)

        all_results = {}
        for period in periods:
            period_name = self.config.time_periods[period]
            self.progress.subsection(f"Period {period}: {period_name}")

            with self.progress.nested():
                result = self.run_single_period(preprocessor, period)
            all_results[period] = result

            self.progress.info(f"Sites: {result['site_probs'].shape[1]}")

        return all_results
