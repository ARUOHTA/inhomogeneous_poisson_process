"""Tests for data preprocessor."""

import os
import tempfile
from unittest.mock import MagicMock, patch

import numpy as np
import polars as pl
import pytest


@pytest.fixture
def mock_elevation_df():
    """Create a mock elevation dataframe."""
    return pl.DataFrame(
        {
            "grid_x": [0, 1, 0, 1],
            "grid_y": [0, 0, 1, 1],
            "x": [138.0, 138.1, 138.0, 138.1],
            "y": [35.0, 35.0, 35.1, 35.1],
            "mesh_code_5th": [1, 2, 3, 4],
            "average_elevation": [100.0, 200.0, 150.0, 250.0],
            "is_sea": [False, False, False, False],
        }
    )


@pytest.fixture
def mock_obsidian_df():
    """Create a mock obsidian dataframe."""
    return pl.DataFrame(
        {
            "遺跡ID": [0, 0, 1, 1],
            "緯度": [35.05, 35.05, 35.15, 35.15],
            "経度": [138.05, 138.05, 138.15, 138.15],
            "時期": [0, 0, 1, 1],
            "産地カテゴリ": ["神津島", "信州", "神津島", "信州"],
        }
    )


def test_preprocessor_import():
    """ObsidianDataPreprocessorがインポートできることを確認"""
    from bayesian_statistics.models.preprocessing.data_preprocessor import (
        ObsidianDataPreprocessor,
    )

    assert ObsidianDataPreprocessor is not None


def test_load_tobler_distances_fallback_warning(mock_elevation_df, mock_obsidian_df):
    """Tobler距離ファイルがない場合にwarningが出ることを確認"""
    from bayesian_statistics.models.preprocessing.data_preprocessor import (
        ObsidianDataPreprocessor,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        preprocessor = ObsidianDataPreprocessor(data_dir=tmpdir)
        # Mock internal dataframes
        preprocessor._df_elevation = mock_elevation_df
        preprocessor._df_obsidian = mock_obsidian_df

        with pytest.warns(UserWarning, match="Tobler距離ディレクトリが見つかりません"):
            distances = preprocessor.load_tobler_distances(
                distance_dir="nonexistent_dir"
            )

        # Should return euclidean distances
        assert distances is not None
        assert isinstance(distances, np.ndarray)


def test_load_tobler_distances_fallback_shape(mock_elevation_df, mock_obsidian_df):
    """フォールバック時に正しい形状の距離行列が返されることを確認"""
    from bayesian_statistics.models.preprocessing.data_preprocessor import (
        ObsidianDataPreprocessor,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        preprocessor = ObsidianDataPreprocessor(data_dir=tmpdir)
        preprocessor._df_elevation = mock_elevation_df
        preprocessor._df_obsidian = mock_obsidian_df

        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            distances = preprocessor.load_tobler_distances(
                distance_dir="nonexistent_dir"
            )

        # Shape should be (n_grid, n_sites)
        n_grid = len(mock_elevation_df)
        n_sites = mock_obsidian_df["遺跡ID"].n_unique()
        assert distances.shape == (n_grid, n_sites)


def test_compute_euclidean_distances():
    """ユークリッド距離が正しく計算されることを確認"""
    from bayesian_statistics.models.preprocessing.data_preprocessor import (
        ObsidianDataPreprocessor,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        preprocessor = ObsidianDataPreprocessor(data_dir=tmpdir)

        grid_coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        site_coords = np.array([[0.0, 0.0], [2.0, 0.0]])

        distances = preprocessor._compute_euclidean_distances(grid_coords, site_coords)

        assert distances.shape == (4, 2)
        # Distance from (0,0) to (0,0) should be 0
        np.testing.assert_almost_equal(distances[0, 0], 0.0)
        # Distance from (0,0) to (2,0) should be 2
        np.testing.assert_almost_equal(distances[0, 1], 2.0)
        # Distance from (1,0) to (2,0) should be 1
        np.testing.assert_almost_equal(distances[1, 1], 1.0)


def test_load_data_not_called_error():
    """load_dataを呼ばずにload_tobler_distancesを呼ぶとエラーになることを確認"""
    from bayesian_statistics.models.preprocessing.data_preprocessor import (
        ObsidianDataPreprocessor,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        preprocessor = ObsidianDataPreprocessor(data_dir=tmpdir)

        with pytest.raises(ValueError, match="load_data"):
            preprocessor.load_tobler_distances()
