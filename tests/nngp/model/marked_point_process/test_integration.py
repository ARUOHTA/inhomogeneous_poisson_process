"""Integration tests for the marked point process model."""

from __future__ import annotations

import numpy as np
import pytest


class TestFullPipeline:
    """End-to-end tests for the marked point process model."""

    def test_full_mcmc_runs_without_error(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """Full MCMC should run without error on small data."""
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessConfig,
            MarkedPointProcessDataset,
            MarkedPointProcessSampler,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=origins,
            grid_coords=small_grid_coords,
            valid_grids=np.ones(small_grid_coords.shape[0], dtype=bool),
            region=region,
        )

        config = MarkedPointProcessConfig(
            n_iter=20,
            burn_in=5,
            thinning=1,
            seed=42,
            neighbor_count=5,
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run(show_progress=False)

        assert results is not None
        assert len(results.lambda_star_samples) == config.n_saved()

    def test_posterior_samples_are_finite(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """All posterior samples should be finite."""
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessConfig,
            MarkedPointProcessDataset,
            MarkedPointProcessSampler,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=origins,
            grid_coords=small_grid_coords,
            valid_grids=np.ones(small_grid_coords.shape[0], dtype=bool),
            region=region,
        )

        config = MarkedPointProcessConfig(
            n_iter=20,
            burn_in=5,
            seed=42,
            neighbor_count=5,
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run(show_progress=False)

        assert np.all(np.isfinite(results.lambda_star_samples))
        assert np.all(np.isfinite(results.beta_mark_samples))

    def test_posterior_probabilities_sum_to_one(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """Composition probabilities should sum to 1 at each site."""
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessConfig,
            MarkedPointProcessDataset,
            MarkedPointProcessSampler,
            compute_eta,
            softmax_with_baseline,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=origins,
            grid_coords=small_grid_coords,
            valid_grids=np.ones(small_grid_coords.shape[0], dtype=bool),
            region=region,
        )

        config = MarkedPointProcessConfig(
            n_iter=20,
            burn_in=5,
            seed=42,
            neighbor_count=5,
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run(show_progress=False)

        # Check last sample
        beta_last = results.beta_mark_samples[-1]
        eta = compute_eta(beta_last, dataset.design_matrix_marks)
        probs = softmax_with_baseline(eta)

        # Probabilities should sum to 1
        np.testing.assert_array_almost_equal(probs.sum(axis=0), np.ones(dataset.num_sites()))

    def test_reproducibility_with_same_seed(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """Same seed should produce same results."""
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessConfig,
            MarkedPointProcessDataset,
            MarkedPointProcessSampler,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=origins,
            grid_coords=small_grid_coords,
            valid_grids=np.ones(small_grid_coords.shape[0], dtype=bool),
            region=region,
        )

        config = MarkedPointProcessConfig(
            n_iter=10,
            burn_in=3,
            seed=12345,
            neighbor_count=5,
        )

        # Run twice with same seed
        sampler1 = MarkedPointProcessSampler(dataset, config)
        results1 = sampler1.run(show_progress=False)

        sampler2 = MarkedPointProcessSampler(dataset, config)
        results2 = sampler2.run(show_progress=False)

        # Results should be identical
        np.testing.assert_array_almost_equal(
            results1.lambda_star_samples,
            results2.lambda_star_samples,
        )


class TestMinimalExample:
    """Tests with minimal valid examples."""

    def test_two_sites_two_categories(self):
        """Simplest case: 2 sites, 2 categories."""
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessConfig,
            MarkedPointProcessDataset,
            MarkedPointProcessSampler,
        )

        coords = np.array([[0.0, 0.0], [1.0, 1.0]])
        counts = np.array([[5, 5], [8, 2]], dtype=float)
        origins = ["A", "B"]
        grid_coords = np.array([[0.5, 0.5]])

        dataset = MarkedPointProcessDataset(
            site_coords=coords,
            counts=counts,
            origins=origins,
            grid_coords=grid_coords,
            valid_grids=np.array([True]),
        )

        config = MarkedPointProcessConfig(
            n_iter=10,
            burn_in=2,
            seed=42,
            neighbor_count=1,  # Only 2 sites
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run(show_progress=False)

        assert results is not None
        assert np.all(np.isfinite(results.lambda_star_samples))


class TestIntensitySpatialVariation:
    """Tests for spatially varying intensity predictions."""

    def test_predict_intensity_varies_spatially(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """Intensity predictions should vary spatially after MCMC."""
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessConfig,
            MarkedPointProcessDataset,
            MarkedPointProcessSampler,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=origins,
            grid_coords=small_grid_coords,
            valid_grids=np.ones(small_grid_coords.shape[0], dtype=bool),
            region=region,
        )

        config = MarkedPointProcessConfig(
            n_iter=30,
            burn_in=10,
            seed=42,
            neighbor_count=5,
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run(show_progress=False)

        # Intensity at sites should not be uniform
        intensity_sites = results.predict_intensity(location="sites")
        assert not np.allclose(intensity_sites, intensity_sites[0])

        # Intensity at grid should also vary
        intensity_grid = results.predict_intensity(location="grid")
        assert not np.allclose(intensity_grid, intensity_grid[0])

    def test_beta_int_samples_not_all_zero(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """Beta_int samples should not remain all zero after MCMC."""
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessConfig,
            MarkedPointProcessDataset,
            MarkedPointProcessSampler,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=origins,
            grid_coords=small_grid_coords,
            valid_grids=np.ones(small_grid_coords.shape[0], dtype=bool),
            region=region,
        )

        config = MarkedPointProcessConfig(
            n_iter=30,
            burn_in=10,
            seed=42,
            neighbor_count=5,
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run(show_progress=False)

        # Beta_int samples should have been updated
        assert not np.allclose(results.beta_int_samples, 0)
        assert np.all(np.isfinite(results.beta_int_samples))


class TestGridPrediction:
    """Tests for grid prediction functionality."""

    def test_predict_probabilities_at_sites(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """predict_probabilities should work at sites."""
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessConfig,
            MarkedPointProcessDataset,
            MarkedPointProcessSampler,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=origins,
            grid_coords=small_grid_coords,
            valid_grids=np.ones(small_grid_coords.shape[0], dtype=bool),
            region=region,
        )

        config = MarkedPointProcessConfig(
            n_iter=20,
            burn_in=5,
            seed=42,
            neighbor_count=5,
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run(show_progress=False)

        # Test site predictions
        probs = results.predict_probabilities(location="sites")
        K = dataset.num_categories()
        n_sites = dataset.num_sites()

        assert probs.shape == (K, n_sites)
        assert np.all(np.isfinite(probs))
        assert np.all(probs >= 0)
        assert np.all(probs <= 1)
        np.testing.assert_array_almost_equal(probs.sum(axis=0), np.ones(n_sites))

    def test_predict_probabilities_at_grid(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """predict_probabilities should work at grid points."""
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessConfig,
            MarkedPointProcessDataset,
            MarkedPointProcessSampler,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=origins,
            grid_coords=small_grid_coords,
            valid_grids=np.ones(small_grid_coords.shape[0], dtype=bool),
            region=region,
        )

        config = MarkedPointProcessConfig(
            n_iter=20,
            burn_in=5,
            seed=42,
            neighbor_count=5,
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run(show_progress=False)

        # Test grid predictions
        probs = results.predict_probabilities(location="grid")
        K = dataset.num_categories()
        n_grid = dataset.num_grid()

        assert probs.shape == (K, n_grid)
        assert np.all(np.isfinite(probs))
        assert np.all(probs >= 0)
        assert np.all(probs <= 1)
        np.testing.assert_array_almost_equal(probs.sum(axis=0), np.ones(n_grid))

    def test_predict_intensity_at_grid(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """predict_intensity should return positive values."""
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessConfig,
            MarkedPointProcessDataset,
            MarkedPointProcessSampler,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=origins,
            grid_coords=small_grid_coords,
            valid_grids=np.ones(small_grid_coords.shape[0], dtype=bool),
            region=region,
        )

        config = MarkedPointProcessConfig(
            n_iter=20,
            burn_in=5,
            seed=42,
            neighbor_count=5,
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run(show_progress=False)

        # Test intensity predictions
        intensity = results.predict_intensity(location="grid")
        n_grid = dataset.num_grid()

        assert intensity.shape == (n_grid,)
        assert np.all(np.isfinite(intensity))
        assert np.all(intensity >= 0)
