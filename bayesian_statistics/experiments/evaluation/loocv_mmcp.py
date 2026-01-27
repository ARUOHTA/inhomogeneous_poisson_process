"""MMCP-specific LOOCV evaluation.

Provides leave-one-out cross validation for Marked Multinomial Composition Process.
Uses subsample approach for computational efficiency.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from bayesian_statistics.experiments.config import ExperimentConfig
from bayesian_statistics.experiments.output import ProgressManager, get_progress_manager
from bayesian_statistics.models.preprocessing.data_preprocessor import (
    ObsidianDataPreprocessor,
)
from bayesian_statistics.nngp.model.marked_point_process import (
    prepare_marked_point_process_dataset,
)

from .loocv_base import LOOCVResult, SubsampleLOOCVEvaluator


class MMCPLOOCVEvaluator(SubsampleLOOCVEvaluator):
    """LOOCV evaluator for MMCP model.

    Uses subsample approach: randomly selects a subset of sites
    for leave-one-out evaluation to reduce computation time.

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
        """Initialize MMCP LOOCV evaluator.

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
        mmcp_runner: Any,
    ) -> Dict[str, Any]:
        """Evaluate MMCP model for a single period using subsample LOOCV.

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            Data preprocessor.
        period : int
            Period index.
        mmcp_runner : MMCPRunner
            MMCP runner instance.

        Returns
        -------
        dict
            Evaluation results with summary statistics.
        """
        # Get available site IDs for this period
        dataset = prepare_marked_point_process_dataset(
            preprocessor=preprocessor,
            period=period,
            origins=self.config.origins,
            distance_column_names=self.config.distance_column_names,
            tau=self.config.tau,
            alpha=self.config.alpha,
            source_weights=self.config.source_weights,
        )

        site_ids = list(range(dataset.num_sites()))
        subsample_ids = self._get_subsample(site_ids)

        results = []
        iterator = self.progress.progress(
            subsample_ids,
            desc=f"LOOCV Period {period}",
        )

        for site_idx in iterator:
            try:
                result = self._evaluate_single_site_loocv(
                    preprocessor=preprocessor,
                    period=period,
                    site_idx=site_idx,
                    mmcp_runner=mmcp_runner,
                )
                if result is not None:
                    results.append(result)
            except Exception as e:
                self.progress.detail(f"Warning: Site {site_idx} evaluation failed: {e}")

        summary = self.compute_summary(results, period=period)
        summary["results"] = results

        return summary

    def _evaluate_single_site_loocv(
        self,
        preprocessor: ObsidianDataPreprocessor,
        period: int,
        site_idx: int,
        mmcp_runner: Any,
    ) -> Optional[LOOCVResult]:
        """Evaluate prediction for a single site with leave-one-out.

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            Data preprocessor.
        period : int
            Period index.
        site_idx : int
            Site index to leave out.
        mmcp_runner : MMCPRunner
            MMCP runner instance.

        Returns
        -------
        LOOCVResult or None
            Evaluation result, or None if evaluation failed.
        """
        # Create dataset excluding target site
        dataset_full = prepare_marked_point_process_dataset(
            preprocessor=preprocessor,
            period=period,
            origins=self.config.origins,
            distance_column_names=self.config.distance_column_names,
            tau=self.config.tau,
            alpha=self.config.alpha,
            source_weights=self.config.source_weights,
        )

        # Get observed composition for target site
        counts = dataset_full.counts[site_idx]
        total = counts.sum()
        if total == 0:
            return None
        observed = counts / total

        # Get site coordinates
        site_coords = dataset_full.site_coords[site_idx]

        # Create excluded dataset and run model
        # Note: This is a simplified approach - full implementation would
        # exclude the site and retrain the model
        dataset_excluded = self._create_excluded_dataset(dataset_full, site_idx)

        if dataset_excluded.num_sites() < 3:
            return None

        # Run MCMC on excluded dataset
        # For efficiency, use fewer iterations for LOOCV
        loocv_config = self._create_loocv_config()

        try:
            from bayesian_statistics.nngp.model.marked_point_process import (
                MarkedPointProcessSampler,
            )

            sampler = MarkedPointProcessSampler(
                dataset_excluded,
                loocv_config,
                show_progress=False,
            )
            results = sampler.run(show_progress=False)

            # Predict at target site location
            predicted = self._predict_at_location(
                results, site_coords, dataset_excluded
            )

            return self.evaluate_single_site(
                site_id=site_idx,
                period=period,
                observed=observed,
                predicted=predicted,
            )
        except Exception as e:
            import traceback
            self.progress.detail(f"LOOCV error for site {site_idx}: {e}")
            self.progress.detail(traceback.format_exc())
            return None

    def _create_excluded_dataset(self, dataset: Any, exclude_idx: int) -> Any:
        """Create dataset with one site excluded.

        Parameters
        ----------
        dataset : MarkedPointProcessDataset
            Original dataset.
        exclude_idx : int
            Index of site to exclude.

        Returns
        -------
        MarkedPointProcessDataset
            Dataset with site excluded.
        """
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessDataset,
        )

        n_sites = dataset.num_sites()
        keep_mask = np.ones(n_sites, dtype=bool)
        keep_mask[exclude_idx] = False

        return MarkedPointProcessDataset(
            site_coords=dataset.site_coords[keep_mask],
            counts=dataset.counts[keep_mask],
            origins=dataset.origins,
            site_ids=dataset.site_ids[keep_mask] if dataset.site_ids is not None else None,
            total_counts=dataset.total_counts[keep_mask] if dataset.total_counts is not None else None,
            design_matrix_intensity=dataset.design_matrix_intensity[keep_mask] if dataset.design_matrix_intensity is not None else None,
            design_matrix_marks=dataset.design_matrix_marks[keep_mask] if dataset.design_matrix_marks is not None else None,
            grid_coords=dataset.grid_coords,
            design_matrix_grid_intensity=dataset.design_matrix_grid_intensity,
            design_matrix_grid_marks=dataset.design_matrix_grid_marks,
            valid_grids=dataset.valid_grids,
            region=dataset.region,
            period=dataset.period,
            distance_features_sites=dataset.distance_features_sites[keep_mask] if dataset.distance_features_sites is not None else None,
            distance_features_grid=dataset.distance_features_grid,
            prior_mean_intercept_sites=dataset.prior_mean_intercept_sites[:, keep_mask] if dataset.prior_mean_intercept_sites is not None else None,
            prior_mean_intercept_grid=dataset.prior_mean_intercept_grid,
            # Note: volume is computed in __post_init__, not passed to constructor
        )

    def _create_loocv_config(self) -> Any:
        """Create reduced config for LOOCV runs.

        Returns
        -------
        MarkedPointProcessConfig
            Configuration with reduced iterations for faster LOOCV.
        """
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessConfig,
        )

        # Use fewer iterations for LOOCV (1/5 of normal)
        n_iter = max(200, self.config.n_iter // 5)
        burn_in = n_iter // 4
        thinning = max(1, self.config.thinning)

        return MarkedPointProcessConfig(
            n_iter=n_iter,
            burn_in=burn_in,
            thinning=thinning,
            neighbor_count=self.config.neighbor_count,
        )

    def _predict_at_location(
        self,
        results: Any,
        site_coords: np.ndarray,
        dataset: Any,
    ) -> np.ndarray:
        """Predict composition at a specific location.

        Parameters
        ----------
        results : MarkedPointProcessResults
            MCMC results.
        site_coords : np.ndarray
            Target location coordinates [lon, lat].
        dataset : MarkedPointProcessDataset
            Dataset used for prediction.

        Returns
        -------
        np.ndarray
            Predicted composition ratios.
        """
        # Find nearest grid point
        grid_coords = dataset.grid_coords
        distances = np.sqrt(
            (grid_coords[:, 0] - site_coords[0]) ** 2
            + (grid_coords[:, 1] - site_coords[1]) ** 2
        )
        nearest_idx = np.argmin(distances)

        # Get prediction at nearest grid point
        grid_probs = results.predict_probabilities(
            location="grid", sample_conditional=False
        )

        return grid_probs[:, nearest_idx]

    def evaluate_all_periods(
        self,
        preprocessor: ObsidianDataPreprocessor,
        mmcp_runner: Any,
    ) -> Dict[int, Dict[str, Any]]:
        """Evaluate MMCP model for all periods.

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            Data preprocessor.
        mmcp_runner : MMCPRunner
            MMCP runner instance.

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
                    mmcp_runner=mmcp_runner,
                )

        return all_results
