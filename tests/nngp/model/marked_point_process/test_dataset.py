"""Tests for MarkedPointProcessDataset."""

from __future__ import annotations

import numpy as np
import pytest


class TestMarkedPointProcessDataset:
    """Tests for the dataset class."""

    def test_dataset_has_site_data(
        self, small_coords: np.ndarray, small_counts: np.ndarray
    ):
        """Dataset should have site coordinates and counts."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=["A", "B", "C"],
        )
        assert dataset.site_coords.shape == (10, 2)
        assert dataset.counts.shape == (10, 3)

    def test_num_sites(self, small_coords: np.ndarray, small_counts: np.ndarray):
        """num_sites should return correct count."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=["A", "B", "C"],
        )
        assert dataset.num_sites() == 10

    def test_num_categories(self, small_coords: np.ndarray, small_counts: np.ndarray):
        """num_categories should return K (including baseline)."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=["A", "B", "C"],
        )
        assert dataset.num_categories() == 3

    def test_total_counts_computed(
        self, small_coords: np.ndarray, small_counts: np.ndarray
    ):
        """total_counts should be computed from counts if not provided."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=["A", "B", "C"],
        )
        # Row sums should match
        expected_totals = small_counts.sum(axis=1)
        np.testing.assert_array_equal(dataset.total_counts, expected_totals)

    def test_design_matrix_defaults_to_intercept(
        self, small_coords: np.ndarray, small_counts: np.ndarray
    ):
        """Default design matrix should be intercept only."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=["A", "B", "C"],
        )
        # Intercept only: (n_sites, 1) matrix of ones
        assert dataset.design_matrix_intensity.shape == (10, 1)
        np.testing.assert_array_equal(dataset.design_matrix_intensity, np.ones((10, 1)))

    def test_design_matrix_marks_defaults_to_intercept(
        self, small_coords: np.ndarray, small_counts: np.ndarray
    ):
        """Default marks design matrix should be intercept only."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=["A", "B", "C"],
        )
        assert dataset.design_matrix_marks.shape == (10, 1)
        np.testing.assert_array_equal(dataset.design_matrix_marks, np.ones((10, 1)))

    def test_region_computed_from_coords(
        self, small_coords: np.ndarray, small_counts: np.ndarray
    ):
        """Region should be inferred from coordinates if not provided."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=["A", "B", "C"],
        )
        # Region should encompass all coordinates
        x_min, x_max = small_coords[:, 0].min(), small_coords[:, 0].max()
        y_min, y_max = small_coords[:, 1].min(), small_coords[:, 1].max()
        assert dataset.region[0][0] <= x_min
        assert dataset.region[0][1] >= x_max
        assert dataset.region[1][0] <= y_min
        assert dataset.region[1][1] >= y_max

    def test_volume_computed_from_region(
        self, small_coords: np.ndarray, small_counts: np.ndarray, region: list
    ):
        """Volume should be area of region."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=["A", "B", "C"],
            region=region,
        )
        x_range = region[0][1] - region[0][0]
        y_range = region[1][1] - region[1][0]
        expected_volume = x_range * y_range
        assert np.isclose(dataset.volume, expected_volume)

    def test_grid_coords_optional(
        self, small_coords: np.ndarray, small_counts: np.ndarray
    ):
        """grid_coords should be optional."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=["A", "B", "C"],
        )
        assert dataset.grid_coords is None

    def test_grid_coords_provided(
        self,
        small_coords: np.ndarray,
        small_counts: np.ndarray,
        small_grid_coords: np.ndarray,
    ):
        """grid_coords can be provided."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=["A", "B", "C"],
            grid_coords=small_grid_coords,
        )
        assert dataset.grid_coords.shape == (25, 2)
        assert dataset.num_grid() == 25

    def test_num_grid_returns_zero_when_no_grid(
        self, small_coords: np.ndarray, small_counts: np.ndarray
    ):
        """num_grid should return 0 when no grid is provided."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        dataset = MarkedPointProcessDataset(
            site_coords=small_coords,
            counts=small_counts,
            origins=["A", "B", "C"],
        )
        assert dataset.num_grid() == 0


class TestDatasetValidation:
    """Tests for input validation."""

    def test_counts_shape_must_match_coords(self, small_coords: np.ndarray):
        """counts rows should match site_coords rows."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        wrong_counts = np.ones((5, 3))  # 5 != 10
        with pytest.raises(ValueError):
            MarkedPointProcessDataset(
                site_coords=small_coords,
                counts=wrong_counts,
                origins=["A", "B", "C"],
            )

    def test_origins_length_must_match_counts_cols(
        self, small_coords: np.ndarray, small_counts: np.ndarray
    ):
        """origins length should match counts columns."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        wrong_origins = ["A", "B"]  # 2 != 3
        with pytest.raises(ValueError):
            MarkedPointProcessDataset(
                site_coords=small_coords,
                counts=small_counts,
                origins=wrong_origins,
            )

    def test_design_matrix_rows_must_match_sites(
        self, small_coords: np.ndarray, small_counts: np.ndarray
    ):
        """design_matrix_intensity rows should match num_sites."""
        from bayesian_statistics.nngp.model.marked_point_process.dataset import (
            MarkedPointProcessDataset,
        )

        wrong_design = np.ones((5, 2))  # 5 != 10
        with pytest.raises(ValueError):
            MarkedPointProcessDataset(
                site_coords=small_coords,
                counts=small_counts,
                origins=["A", "B", "C"],
                design_matrix_intensity=wrong_design,
            )
