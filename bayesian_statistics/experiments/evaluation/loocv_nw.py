"""NW-specific LOOCV evaluation.

Provides leave-one-out cross validation for Nadaraya-Watson estimator.
Wraps the existing LOOCVEvaluator from models.evaluation.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np

from bayesian_statistics.experiments.config import ExperimentConfig
from bayesian_statistics.experiments.output import ProgressManager, get_progress_manager
from bayesian_statistics.models.composition.nadaraya_watson import (
    NadarayaWatsonEstimator,
)
from bayesian_statistics.models.preprocessing.data_preprocessor import (
    ObsidianDataPreprocessor,
)

from .loocv_base import LOOCVResult, SubsampleLOOCVEvaluator


class NWLOOCVEvaluator(SubsampleLOOCVEvaluator):
    """LOOCV evaluator for Nadaraya-Watson estimator.

    Uses subsample approach for computational efficiency.
    Simpler than MMCP as NW doesn't require MCMC.

    Attributes
    ----------
    config : ExperimentConfig
        Experiment configuration.
    progress : ProgressManager
        Progress manager for output control.
    """

    def __init__(
        self,
        config: ExperimentConfig,
        progress: ProgressManager | None = None,
    ):
        """Initialize NW LOOCV evaluator.

        Parameters
        ----------
        config : ExperimentConfig
            Experiment configuration with LOOCV settings.
        progress : ProgressManager, optional
            Progress manager for output control. If None, uses global default.
        """
        super().__init__(
            n_samples=config.loocv_n_samples,
            seed=config.loocv_seed,
        )
        self.config = config
        self.progress = progress or get_progress_manager()

    def evaluate_period(
        self,
        preprocessor: ObsidianDataPreprocessor,
        period: int,
    ) -> Dict[str, Any]:
        """Evaluate NW model for a single period using subsample LOOCV.

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            Data preprocessor.
        period : int
            Period index.

        Returns
        -------
        dict
            Evaluation results with summary statistics.
        """
        # Get site information
        site_ids = preprocessor.df_sites["遺跡ID"].unique().sort().to_list()
        subsample_ids = self._get_subsample(list(range(len(site_ids))))

        # Compute observed compositions for all sites
        observed_compositions = self._compute_observed_compositions(
            preprocessor, period
        )

        results = []
        iterator = self.progress.progress(
            subsample_ids,
            desc=f"LOOCV Period {period}",
        )

        for site_idx in iterator:
            site_id = site_ids[site_idx]
            if site_id not in observed_compositions:
                continue

            try:
                result = self._evaluate_single_site_loocv(
                    preprocessor=preprocessor,
                    period=period,
                    site_id=site_id,
                    site_idx=site_idx,
                    observed_compositions=observed_compositions,
                )
                if result is not None:
                    results.append(result)
            except Exception as e:
                self.progress.detail(f"Warning: Site {site_id} evaluation failed: {e}")

        summary = self.compute_summary(results, period=period)
        summary["results"] = results

        return summary

    def _evaluate_single_site_loocv(
        self,
        preprocessor: ObsidianDataPreprocessor,
        period: int,
        site_id: int,
        site_idx: int,
        observed_compositions: Dict[int, np.ndarray],
    ) -> LOOCVResult | None:
        """Evaluate prediction for a single site with leave-one-out.

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            Data preprocessor.
        period : int
            Period index.
        site_id : int
            Site ID to leave out.
        site_idx : int
            Site index.
        observed_compositions : dict
            Precomputed observed compositions.

        Returns
        -------
        LOOCVResult or None
            Evaluation result, or None if evaluation failed.
        """
        observed = observed_compositions.get(site_id)
        if observed is None:
            return None

        # Skip if all zeros
        if np.sum(observed) == 0:
            return None

        # Create excluded preprocessor
        excluded_preprocessor = self._create_excluded_preprocessor(
            preprocessor, site_id
        )

        if excluded_preprocessor.df_sites.height < 3:
            return None

        # Fit NW on excluded data
        nw = NadarayaWatsonEstimator(
            sigma=self.config.nw_sigma,
            sigma_for_sites=self.config.nw_sigma_for_sites,
            variable_names=self.config.nw_variable_names,
        )
        nw.fit(excluded_preprocessor)

        # Get target site coordinates
        target_site = preprocessor.df_sites.filter(
            preprocessor.df_sites["遺跡ID"] == site_id
        )
        if target_site.height == 0:
            return None

        target_lon = target_site["経度"].item()
        target_lat = target_site["緯度"].item()

        # Predict at target location
        predicted = self._predict_at_location(
            nw, excluded_preprocessor, period, target_lon, target_lat
        )

        return self.evaluate_single_site(
            site_id=site_id,
            period=period,
            observed=observed,
            predicted=predicted,
        )

    def _create_excluded_preprocessor(
        self,
        preprocessor: ObsidianDataPreprocessor,
        exclude_site_id: int,
    ) -> ObsidianDataPreprocessor:
        """Create preprocessor with one site excluded.

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            Original preprocessor.
        exclude_site_id : int
            Site ID to exclude.

        Returns
        -------
        ObsidianDataPreprocessor
            Preprocessor with site excluded.
        """
        import polars as pl

        # Create new preprocessor
        excluded = ObsidianDataPreprocessor(preprocessor.data_dir)
        excluded.load_data()

        # Filter out the excluded site
        excluded._df_obsidian = excluded.df_obsidian.filter(
            pl.col("遺跡ID") != exclude_site_id
        )
        excluded._df_sites = excluded.df_sites.filter(
            pl.col("遺跡ID") != exclude_site_id
        )

        return excluded

    def _compute_observed_compositions(
        self,
        preprocessor: ObsidianDataPreprocessor,
        period: int,
    ) -> Dict[int, np.ndarray]:
        """Compute observed compositions for all sites.

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            Data preprocessor.
        period : int
            Period index.

        Returns
        -------
        dict
            Dictionary mapping site_id -> composition array.
        """
        import polars as pl

        origins = self.config.origins
        compositions = {}

        period_data = preprocessor.df_obsidian.filter(pl.col("時期") == period)
        site_ids = preprocessor.df_sites["遺跡ID"].unique().sort().to_list()

        for site_id in site_ids:
            site_data = period_data.filter(pl.col("遺跡ID") == site_id)
            total = site_data.height

            if total == 0:
                continue

            composition = np.zeros(len(origins))
            for i, origin in enumerate(origins):
                if origin == "その他":
                    # Count everything not in other categories
                    other_count = site_data.filter(
                        ~pl.col("産地カテゴリ").is_in(origins[:-1])
                    ).height
                    composition[i] = other_count / total
                else:
                    count = site_data.filter(pl.col("産地カテゴリ") == origin).height
                    composition[i] = count / total

            # Normalize
            if np.sum(composition) > 0:
                composition = composition / np.sum(composition)

            compositions[site_id] = composition

        return compositions

    def _predict_at_location(
        self,
        nw: NadarayaWatsonEstimator,
        preprocessor: ObsidianDataPreprocessor,
        period: int,
        target_lon: float,
        target_lat: float,
    ) -> np.ndarray:
        """Predict composition at a specific location.

        Parameters
        ----------
        nw : NadarayaWatsonEstimator
            Fitted NW estimator.
        preprocessor : ObsidianDataPreprocessor
            Data preprocessor.
        period : int
            Period index.
        target_lon : float
            Target longitude.
        target_lat : float
            Target latitude.

        Returns
        -------
        np.ndarray
            Predicted composition ratios.
        """
        origins = self.config.origins
        predicted = np.zeros(len(origins))

        # Get mesh grid
        lon_mesh, lat_mesh = preprocessor.create_meshgrid()

        # Find nearest grid point
        distances = np.sqrt((lon_mesh - target_lon) ** 2 + (lat_mesh - target_lat) ** 2)
        min_idx = np.unravel_index(np.argmin(distances), distances.shape)
        grid_flat_idx = min_idx[0] * lon_mesh.shape[1] + min_idx[1]

        # Predict for each origin
        for i, origin in enumerate(origins[:-1]):  # Exclude "その他"
            try:
                counts, target_counts = preprocessor.preprocess_obsidian_data(
                    period, origin
                )

                # Get weights for this grid point
                weights = nw.weights[grid_flat_idx, :]

                # Compute weighted ratio
                valid_mask = counts > 0
                if np.sum(valid_mask) > 0:
                    weighted_total = np.sum(weights[valid_mask] * counts[valid_mask])
                    weighted_target = np.sum(
                        weights[valid_mask] * target_counts[valid_mask]
                    )

                    if weighted_total > 0:
                        predicted[i] = weighted_target / weighted_total
            except Exception:
                predicted[i] = 0.0

        # Compute "その他" as remainder
        predicted[-1] = max(0.0, 1.0 - np.sum(predicted[:-1]))

        # Normalize
        total = np.sum(predicted)
        if total > 0:
            predicted = predicted / total

        return predicted

    def evaluate_all_periods(
        self,
        preprocessor: ObsidianDataPreprocessor,
    ) -> Dict[int, Dict[str, Any]]:
        """Evaluate NW model for all periods.

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            Data preprocessor.

        Returns
        -------
        dict
            Dictionary mapping period -> evaluation results.
        """
        all_results = {}

        for period in self.config.time_periods.keys():
            period_name = self.config.time_periods[period]
            self.progress.info(f"\nPeriod {period} ({period_name}) LOOCV")

            with self.progress.nested():
                all_results[period] = self.evaluate_period(
                    preprocessor=preprocessor,
                    period=period,
                )

        return all_results
