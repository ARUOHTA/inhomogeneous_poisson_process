"""Tests for LOOCV evaluation module."""

from __future__ import annotations

import numpy as np
import pytest

from bayesian_statistics.experiments.config import ExperimentConfig
from bayesian_statistics.experiments.evaluation.loocv_base import (
    LOOCVResult,
    SubsampleLOOCVEvaluator,
)
from bayesian_statistics.models.evaluation.metrics import CompositionalMetrics


class TestCompositionalMetricsIntegration:
    """Test that CompositionalMetrics works correctly for LOOCV."""

    def test_aitchison_distance_positive(self):
        """Aitchison距離が正の値であることを確認."""
        metrics = CompositionalMetrics()
        x = np.array([0.4, 0.3, 0.2, 0.1])
        y = np.array([0.3, 0.3, 0.3, 0.1])

        distance = metrics.aitchison_distance(x, y)

        assert distance > 0
        assert np.isfinite(distance)

    def test_aitchison_distance_zero_for_identical(self):
        """同一ベクトルの距離がゼロであることを確認."""
        metrics = CompositionalMetrics()
        x = np.array([0.4, 0.3, 0.2, 0.1])

        distance = metrics.aitchison_distance(x, x)

        assert np.isclose(distance, 0, atol=1e-10)

    def test_aitchison_distance_symmetric(self):
        """Aitchison距離が対称であることを確認."""
        metrics = CompositionalMetrics()
        x = np.array([0.4, 0.3, 0.2, 0.1])
        y = np.array([0.3, 0.3, 0.3, 0.1])

        d_xy = metrics.aitchison_distance(x, y)
        d_yx = metrics.aitchison_distance(y, x)

        assert np.isclose(d_xy, d_yx)

    def test_aitchison_distance_handles_zeros(self):
        """ゼロ値を含む場合も計算できることを確認."""
        metrics = CompositionalMetrics()
        x = np.array([0.5, 0.5, 0.0, 0.0])
        y = np.array([0.3, 0.3, 0.3, 0.1])

        distance = metrics.aitchison_distance(x, y)

        assert distance > 0
        assert np.isfinite(distance)


class TestLOOCVResult:
    """Test LOOCVResult data structure."""

    def test_loocv_result_creation(self):
        """LOOCVResultが正しく作成されることを確認."""
        result = LOOCVResult(
            period=2,
            site_id=123,
            observed=np.array([0.4, 0.3, 0.2, 0.1]),
            predicted=np.array([0.35, 0.35, 0.2, 0.1]),
            aitchison_distance=0.15,
        )

        assert result.period == 2
        assert result.site_id == 123
        assert result.aitchison_distance == 0.15
        assert len(result.observed) == 4
        assert len(result.predicted) == 4

    def test_loocv_result_fields(self):
        """LOOCVResultが必要なフィールドを持つことを確認."""
        result = LOOCVResult(
            period=0,
            site_id=1,
            observed=np.array([0.5, 0.5]),
            predicted=np.array([0.5, 0.5]),
            aitchison_distance=0.0,
        )

        # Required fields
        assert hasattr(result, "period")
        assert hasattr(result, "site_id")
        assert hasattr(result, "observed")
        assert hasattr(result, "predicted")
        assert hasattr(result, "aitchison_distance")


class TestSubsampleLOOCVEvaluator:
    """Test SubsampleLOOCVEvaluator."""

    def test_evaluator_creation_with_config(self):
        """Configからevaluatorが作成できることを確認."""
        config = ExperimentConfig()
        evaluator = SubsampleLOOCVEvaluator(
            n_samples=config.loocv_n_samples,
            seed=config.loocv_seed,
        )

        assert evaluator.n_samples == 20
        assert evaluator.seed == 42

    def test_evaluator_default_values(self):
        """デフォルト値が正しいことを確認."""
        evaluator = SubsampleLOOCVEvaluator()

        assert evaluator.n_samples == 20
        assert evaluator.seed == 42

    def test_evaluator_custom_values(self):
        """カスタム値が設定できることを確認."""
        evaluator = SubsampleLOOCVEvaluator(n_samples=10, seed=123)

        assert evaluator.n_samples == 10
        assert evaluator.seed == 123

    def test_subsample_reproducibility(self):
        """同じシードで同じサブサンプルが得られることを確認."""
        evaluator1 = SubsampleLOOCVEvaluator(n_samples=5, seed=42)
        evaluator2 = SubsampleLOOCVEvaluator(n_samples=5, seed=42)

        # Mock site_ids list
        site_ids = list(range(100))

        subsample1 = evaluator1._get_subsample(site_ids)
        evaluator2._reset_random_state()  # Reset before second call
        subsample2 = evaluator2._get_subsample(site_ids)

        assert subsample1 == subsample2

    def test_subsample_different_seeds(self):
        """異なるシードで異なるサブサンプルが得られることを確認."""
        evaluator1 = SubsampleLOOCVEvaluator(n_samples=5, seed=42)
        evaluator2 = SubsampleLOOCVEvaluator(n_samples=5, seed=123)

        site_ids = list(range(100))

        subsample1 = evaluator1._get_subsample(site_ids)
        subsample2 = evaluator2._get_subsample(site_ids)

        # Should be different (with high probability)
        assert subsample1 != subsample2

    def test_subsample_count(self):
        """サブサンプル数が設定通りであることを確認."""
        evaluator = SubsampleLOOCVEvaluator(n_samples=15, seed=42)
        site_ids = list(range(100))

        subsample = evaluator._get_subsample(site_ids)

        assert len(subsample) == 15

    def test_subsample_from_smaller_list(self):
        """サイト数がn_samples未満の場合、全サイトが返されることを確認."""
        evaluator = SubsampleLOOCVEvaluator(n_samples=20, seed=42)
        site_ids = list(range(10))  # Only 10 sites

        subsample = evaluator._get_subsample(site_ids)

        assert len(subsample) == 10
        assert set(subsample) == set(site_ids)

    def test_compute_mean_distance(self):
        """平均距離が正しく計算されることを確認."""
        evaluator = SubsampleLOOCVEvaluator()

        results = [
            LOOCVResult(0, 1, np.array([0.5, 0.5]), np.array([0.5, 0.5]), 0.1),
            LOOCVResult(0, 2, np.array([0.5, 0.5]), np.array([0.5, 0.5]), 0.2),
            LOOCVResult(0, 3, np.array([0.5, 0.5]), np.array([0.5, 0.5]), 0.3),
        ]

        mean_distance = evaluator._compute_mean_distance(results)

        assert np.isclose(mean_distance, 0.2)

    def test_compute_std_distance(self):
        """標準偏差が正しく計算されることを確認."""
        evaluator = SubsampleLOOCVEvaluator()

        results = [
            LOOCVResult(0, 1, np.array([0.5, 0.5]), np.array([0.5, 0.5]), 0.1),
            LOOCVResult(0, 2, np.array([0.5, 0.5]), np.array([0.5, 0.5]), 0.2),
            LOOCVResult(0, 3, np.array([0.5, 0.5]), np.array([0.5, 0.5]), 0.3),
        ]

        std_distance = evaluator._compute_std_distance(results)

        expected_std = np.std([0.1, 0.2, 0.3])
        assert np.isclose(std_distance, expected_std)


class TestLOOCVSummary:
    """Test LOOCV summary statistics."""

    def test_summary_statistics_structure(self):
        """サマリー統計が正しい構造を持つことを確認."""
        evaluator = SubsampleLOOCVEvaluator()

        results = [
            LOOCVResult(0, 1, np.array([0.5, 0.5]), np.array([0.5, 0.5]), 0.1),
            LOOCVResult(0, 2, np.array([0.5, 0.5]), np.array([0.5, 0.5]), 0.2),
        ]

        summary = evaluator.compute_summary(results)

        assert "mean_aitchison_distance" in summary
        assert "std_aitchison_distance" in summary
        assert "n_samples" in summary
        assert "period" in summary

    def test_summary_for_empty_results(self):
        """空の結果でもサマリーが作成できることを確認."""
        evaluator = SubsampleLOOCVEvaluator()

        summary = evaluator.compute_summary([])

        assert summary["n_samples"] == 0
        assert np.isnan(summary["mean_aitchison_distance"])
