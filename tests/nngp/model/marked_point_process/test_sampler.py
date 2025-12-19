"""Tests for MarkedPointProcessSampler."""

from __future__ import annotations

import numpy as np
import pytest


class TestSamplerInitialization:
    """Tests for sampler initialization."""

    def test_sampler_initializes_without_error(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        origins: list,
    ):
        """Sampler should initialize without error."""
        from bayesian_statistics.nngp.model.marked_point_process.config import (
            MarkedPointProcessConfig,
        )
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )
        from bayesian_statistics.nngp.model.marked_point_process.sampler import (
            MarkedPointProcessSampler,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=origins,
        )
        config = MarkedPointProcessConfig(n_iter=10, burn_in=2)
        sampler = MarkedPointProcessSampler(dataset, config)

        assert sampler is not None

    def test_sampler_has_rng(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        origins: list,
    ):
        """Sampler should have a random generator."""
        from bayesian_statistics.nngp.model.marked_point_process.config import (
            MarkedPointProcessConfig,
        )
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )
        from bayesian_statistics.nngp.model.marked_point_process.sampler import (
            MarkedPointProcessSampler,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=origins,
        )
        config = MarkedPointProcessConfig(n_iter=10, seed=42)
        sampler = MarkedPointProcessSampler(dataset, config)

        assert sampler.rng is not None


class TestSamplerRun:
    """Tests for sampler run method."""

    def test_sampler_runs_without_error(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """Sampler should run without error on small data."""
        from bayesian_statistics.nngp.model.marked_point_process.config import (
            MarkedPointProcessConfig,
        )
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )
        from bayesian_statistics.nngp.model.marked_point_process.sampler import (
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
            n_iter=5,
            burn_in=2,
            thinning=1,
            seed=42,
            neighbor_count=5,  # Small for test
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run()

        assert results is not None

    def test_sampler_returns_correct_sample_count(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """Results should have correct number of saved samples."""
        from bayesian_statistics.nngp.model.marked_point_process.config import (
            MarkedPointProcessConfig,
        )
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )
        from bayesian_statistics.nngp.model.marked_point_process.sampler import (
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
            burn_in=4,
            thinning=2,
            seed=42,
            neighbor_count=5,
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run()

        # (10 - 4) / 2 = 3 samples
        expected_n_samples = config.n_saved()
        assert results.lambda_star_samples.shape[0] == expected_n_samples


class TestSamplerLambdaStar:
    """Tests for lambda* sampling."""

    def test_lambda_star_samples_positive(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """All lambda* samples should be positive."""
        from bayesian_statistics.nngp.model.marked_point_process.config import (
            MarkedPointProcessConfig,
        )
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )
        from bayesian_statistics.nngp.model.marked_point_process.sampler import (
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
            burn_in=2,
            seed=42,
            neighbor_count=5,
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run()

        assert np.all(results.lambda_star_samples > 0)


class TestSamplerNoNaN:
    """Tests for numerical stability."""

    def test_samples_have_no_nan(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
        origins: list,
        region: list,
    ):
        """All samples should be finite (no NaN or Inf)."""
        from bayesian_statistics.nngp.model.marked_point_process.config import (
            MarkedPointProcessConfig,
        )
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )
        from bayesian_statistics.nngp.model.marked_point_process.sampler import (
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
            burn_in=2,
            seed=42,
            neighbor_count=5,
        )

        sampler = MarkedPointProcessSampler(dataset, config)
        results = sampler.run()

        assert np.all(np.isfinite(results.lambda_star_samples))
