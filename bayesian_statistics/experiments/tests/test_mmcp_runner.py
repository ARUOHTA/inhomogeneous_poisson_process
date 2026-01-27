"""Tests for MMCP runner."""

from pathlib import Path
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

    return ExperimentConfig(
        n_iter=50,
        burn_in=10,
        thinning=2,
    )


def test_mmcp_runner_import():
    """MMCPRunnerがインポートできることを確認"""
    from bayesian_statistics.experiments.models.mmcp_runner import MMCPRunner

    assert MMCPRunner is not None


def test_mmcp_runner_init(mock_config):
    """MMCPRunnerが正しく初期化されることを確認"""
    from bayesian_statistics.experiments.models.mmcp_runner import MMCPRunner

    runner = MMCPRunner(config=mock_config)

    assert runner.config == mock_config


def test_create_mmcp_config(mock_config):
    """MCMCConfigが正しく作成されることを確認"""
    from bayesian_statistics.experiments.models.mmcp_runner import MMCPRunner

    runner = MMCPRunner(config=mock_config)
    mmcp_config = runner._create_mmcp_config()

    assert mmcp_config.n_iter == mock_config.n_iter
    assert mmcp_config.burn_in == mock_config.burn_in
    assert mmcp_config.thinning == mock_config.thinning
    assert mmcp_config.neighbor_count == mock_config.neighbor_count
    assert mmcp_config.tau == mock_config.tau
    assert mmcp_config.alpha == mock_config.alpha


def test_run_single_period_returns_dict(mock_config, mock_preprocessor):
    """run_single_periodが辞書を返すことを確認"""
    from bayesian_statistics.experiments.models.mmcp_runner import MMCPRunner

    runner = MMCPRunner(config=mock_config)

    # Mock the internal methods
    with patch.object(runner, "_run_mcmc") as mock_run:
        mock_results = MagicMock()
        mock_results.predict_probabilities.return_value = np.random.rand(5, 10)
        mock_results.decompose_effects.return_value = {
            "distance": np.random.rand(5, 10),
            "intercept_adjustment": np.random.rand(5, 10),
            "full": np.random.rand(5, 10),
        }
        mock_results.lambda_star_samples = np.random.rand(20)
        mock_dataset = MagicMock()
        mock_dataset.counts = np.random.rand(10, 5)
        mock_run.return_value = (mock_results, mock_dataset)

        result = runner.run_single_period(mock_preprocessor, period=2)

        assert isinstance(result, dict)
        assert "results" in result
        assert "dataset" in result
        assert "site_probs" in result
        assert "grid_probs" in result
        assert "effects" in result


def test_results_have_required_keys(mock_config, mock_preprocessor):
    """結果が必要なキーを持つことを確認"""
    from bayesian_statistics.experiments.models.mmcp_runner import MMCPRunner

    runner = MMCPRunner(config=mock_config)

    with patch.object(runner, "_run_mcmc") as mock_run:
        mock_results = MagicMock()
        mock_results.predict_probabilities.return_value = np.random.rand(5, 10)
        mock_results.decompose_effects.return_value = {
            "distance": np.random.rand(5, 10),
            "intercept_adjustment": np.random.rand(5, 10),
            "full": np.random.rand(5, 10),
        }
        mock_results.lambda_star_samples = np.random.rand(20)
        mock_dataset = MagicMock()
        mock_dataset.counts = np.random.rand(10, 5)
        mock_run.return_value = (mock_results, mock_dataset)

        result = runner.run_single_period(mock_preprocessor, period=2)

        # Check effects keys
        effects = result["effects"]
        assert "distance" in effects
        assert "intercept_adjustment" in effects
        assert "full" in effects


def test_decompose_effects_keys(mock_config, mock_preprocessor):
    """効果分解が正しいキーを持つことを確認"""
    from bayesian_statistics.experiments.models.mmcp_runner import MMCPRunner

    runner = MMCPRunner(config=mock_config)

    with patch.object(runner, "_run_mcmc") as mock_run:
        mock_results = MagicMock()
        mock_results.predict_probabilities.return_value = np.random.rand(5, 10)
        mock_effects = {
            "distance": np.random.rand(5, 10),
            "intercept_adjustment": np.random.rand(5, 10),
            "intercept": np.random.rand(5, 10),
            "full": np.random.rand(5, 10),
        }
        mock_results.decompose_effects.return_value = mock_effects
        mock_results.lambda_star_samples = np.random.rand(20)
        mock_dataset = MagicMock()
        mock_dataset.counts = np.random.rand(10, 5)
        mock_run.return_value = (mock_results, mock_dataset)

        result = runner.run_single_period(mock_preprocessor, period=2)

        effects = result["effects"]
        assert "distance" in effects
        assert "intercept_adjustment" in effects
        assert "full" in effects
