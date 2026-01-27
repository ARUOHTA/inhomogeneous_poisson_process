"""Tests for NW (Nadaraya-Watson) runner."""

from unittest.mock import MagicMock, patch

import numpy as np
import pytest


@pytest.fixture
def mock_preprocessor():
    """Create a mock preprocessor."""
    preprocessor = MagicMock()
    preprocessor.df_elevation = MagicMock()
    return preprocessor


@pytest.fixture
def mock_config():
    """Create a mock experiment config."""
    from bayesian_statistics.experiments.config import ExperimentConfig

    return ExperimentConfig()


def test_nw_runner_import():
    """NWRunnerがインポートできることを確認"""
    from bayesian_statistics.experiments.models.nw_runner import NWRunner

    assert NWRunner is not None


def test_nw_runner_init(mock_config):
    """NWRunnerが正しく初期化されることを確認"""
    from bayesian_statistics.experiments.models.nw_runner import NWRunner

    runner = NWRunner(config=mock_config)

    assert runner.config == mock_config


def test_nw_runner_has_same_interface_as_mmcp(mock_config):
    """NWRunnerがMMCPRunnerと同じインターフェースを持つことを確認"""
    from bayesian_statistics.experiments.models.nw_runner import NWRunner
    from bayesian_statistics.experiments.models.mmcp_runner import MMCPRunner

    nw_runner = NWRunner(config=mock_config)
    mmcp_runner = MMCPRunner(config=mock_config)

    # Check that both have the same key methods
    assert hasattr(nw_runner, "run_single_period")
    assert hasattr(nw_runner, "run_all_periods")

    assert hasattr(mmcp_runner, "run_single_period")
    assert hasattr(mmcp_runner, "run_all_periods")


def test_run_single_period_returns_dict(mock_config, mock_preprocessor):
    """run_single_periodが辞書を返すことを確認"""
    from bayesian_statistics.experiments.models.nw_runner import NWRunner

    runner = NWRunner(config=mock_config)

    # Mock the internal NW estimator
    with patch.object(runner, "_run_nw") as mock_run:
        mock_result = {
            "site_probs": np.random.rand(5, 10),
            "grid_probs": np.random.rand(5, 100),
            "estimator": MagicMock(),
        }
        mock_run.return_value = mock_result

        result = runner.run_single_period(mock_preprocessor, period=2)

        assert isinstance(result, dict)
        assert "site_probs" in result
        assert "grid_probs" in result
        assert "period" in result


def test_result_probs_shape(mock_config, mock_preprocessor):
    """結果の確率配列が正しい形状を持つことを確認"""
    from bayesian_statistics.experiments.models.nw_runner import NWRunner

    runner = NWRunner(config=mock_config)

    n_origins = 5
    n_sites = 10
    n_grid = 100

    with patch.object(runner, "_run_nw") as mock_run:
        mock_result = {
            "site_probs": np.random.rand(n_origins, n_sites),
            "grid_probs": np.random.rand(n_origins, n_grid),
            "estimator": MagicMock(),
        }
        mock_run.return_value = mock_result

        result = runner.run_single_period(mock_preprocessor, period=2)

        assert result["site_probs"].shape[0] == n_origins
        assert result["grid_probs"].shape[0] == n_origins


def test_probs_are_valid(mock_config, mock_preprocessor):
    """確率が[0, 1]の範囲で、各位置で合計1になることを確認"""
    from bayesian_statistics.experiments.models.nw_runner import NWRunner

    runner = NWRunner(config=mock_config)

    # Create valid probability arrays (sum to 1 along first axis)
    n_origins = 5
    n_sites = 10
    raw_probs = np.random.rand(n_origins, n_sites)
    site_probs = raw_probs / raw_probs.sum(axis=0, keepdims=True)

    raw_grid = np.random.rand(n_origins, 100)
    grid_probs = raw_grid / raw_grid.sum(axis=0, keepdims=True)

    with patch.object(runner, "_run_nw") as mock_run:
        mock_result = {
            "site_probs": site_probs,
            "grid_probs": grid_probs,
            "estimator": MagicMock(),
        }
        mock_run.return_value = mock_result

        result = runner.run_single_period(mock_preprocessor, period=2)

        # Check probability bounds
        assert np.all(result["site_probs"] >= 0)
        assert np.all(result["site_probs"] <= 1)
        assert np.all(result["grid_probs"] >= 0)
        assert np.all(result["grid_probs"] <= 1)

        # Check sum to 1
        np.testing.assert_array_almost_equal(
            result["site_probs"].sum(axis=0), np.ones(n_sites)
        )
