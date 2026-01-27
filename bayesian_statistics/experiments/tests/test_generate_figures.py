"""Tests for figure generation functions.

TDD approach: Write failing tests first, then fix the implementation.
"""

from __future__ import annotations

from typing import TYPE_CHECKING
from unittest.mock import MagicMock

import numpy as np
import pytest

from bayesian_statistics.experiments.config import ExperimentConfig


class TestDistancePriorComputation:
    """Test the distance prior p0 computation for Figure 5.5."""

    def test_compute_p0_from_distance_features(self):
        """p0 should be computed from distance_features_grid using softmax.

        MarkedPointProcessDataset has:
        - distance_features_grid: (n_grid, K-1) - log-ratio features
        - prior_mean_intercept_grid: (K-1, n_grid) - lambda * features

        p0 should be computed by applying softmax to the log-ratios.
        """
        # Arrange: Mock dataset with distance features
        n_grid = 100
        K_minus_1 = 4  # 4 categories excluding baseline

        # distance_features_grid contains log-ratio: log(p_k) - log(p_K)
        distance_features = np.random.randn(n_grid, K_minus_1)

        # Act: Compute p0 from log-ratios using softmax
        # log(p_k / p_K) = g_k => p_k = exp(g_k) * p_K
        # Sum constraint: sum(p) = 1 => p_K + sum(exp(g_k) * p_K) = 1
        # p_K = 1 / (1 + sum(exp(g_k)))
        exp_g = np.exp(distance_features)  # (n_grid, K-1)
        sum_exp_g = exp_g.sum(axis=1, keepdims=True)  # (n_grid, 1)
        p_K = 1.0 / (1.0 + sum_exp_g.ravel())  # (n_grid,)
        p0 = exp_g * p_K[:, np.newaxis]  # (n_grid, K-1)

        # Add baseline category
        p0_full = np.column_stack([p0, p_K])  # (n_grid, K)

        # Assert
        assert p0_full.shape == (n_grid, K_minus_1 + 1)
        # Should sum to 1 for each grid point
        np.testing.assert_allclose(p0_full.sum(axis=1), 1.0, rtol=1e-10)
        # All probabilities should be positive
        assert np.all(p0_full > 0)

    def test_compute_p0_helper_function(self):
        """Test the helper function that computes p0 from dataset."""
        from bayesian_statistics.experiments.run_all import _compute_p0_from_dataset

        # Create mock dataset
        n_grid = 50
        K_minus_1 = 4

        mock_dataset = MagicMock()
        mock_dataset.distance_features_grid = np.random.randn(n_grid, K_minus_1)

        # Act
        p0 = _compute_p0_from_dataset(mock_dataset)

        # Assert: plot_distance_prior expects (n_grid, K) format
        assert p0.shape == (n_grid, K_minus_1 + 1)  # (n_grid, K) format
        np.testing.assert_allclose(p0.sum(axis=1), 1.0, rtol=1e-10)

    def test_compute_p0_handles_none_distance_features(self):
        """Should handle dataset without distance features gracefully."""
        from bayesian_statistics.experiments.run_all import _compute_p0_from_dataset

        mock_dataset = MagicMock()
        mock_dataset.distance_features_grid = None

        # Act
        p0 = _compute_p0_from_dataset(mock_dataset)

        # Assert: Should return None when no distance features
        assert p0 is None


class TestGridSizeValidation:
    """Test grid size validation for model comparison."""

    def test_check_grid_sizes_match_returns_true_when_same(self):
        """Should return True when grid sizes match."""
        from bayesian_statistics.experiments.run_all import _check_grid_sizes_match

        mmcp_probs = np.random.rand(5, 13575)  # (K, n_grid)
        nw_probs = np.random.rand(5, 13575)  # (K, n_grid)

        assert _check_grid_sizes_match(mmcp_probs, nw_probs) is True

    def test_check_grid_sizes_match_returns_false_when_different(self):
        """Should return False when grid sizes don't match."""
        from bayesian_statistics.experiments.run_all import _check_grid_sizes_match

        mmcp_probs = np.random.rand(5, 13575)  # (K, n_grid_subsampled)
        nw_probs = np.random.rand(5, 1357520)  # (K, n_grid_full)

        assert _check_grid_sizes_match(mmcp_probs, nw_probs) is False


class TestMarkedPointProcessDatasetAttributes:
    """Verify MarkedPointProcessDataset has expected attributes."""

    def test_dataset_has_distance_features_grid(self):
        """Dataset should have distance_features_grid (not distance_zscores_grid)."""
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessDataset,
        )

        # Check attribute exists in dataclass definition
        import dataclasses

        field_names = [f.name for f in dataclasses.fields(MarkedPointProcessDataset)]
        assert "distance_features_grid" in field_names
        assert "distance_zscores_grid" not in field_names

    def test_dataset_does_not_have_source_weights_attribute(self):
        """Dataset should NOT have source_weights_full, tau, alpha attributes."""
        from bayesian_statistics.nngp.model.marked_point_process import (
            MarkedPointProcessDataset,
        )

        import dataclasses

        field_names = [f.name for f in dataclasses.fields(MarkedPointProcessDataset)]
        assert "source_weights_full" not in field_names
        assert "tau" not in field_names
        assert "alpha" not in field_names
