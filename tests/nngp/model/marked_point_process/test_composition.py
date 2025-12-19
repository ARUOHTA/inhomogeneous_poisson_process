"""Tests for CompositionComponent (mark/composition part)."""

from __future__ import annotations

import numpy as np
import pytest


class TestComputeKappaTildeMark:
    """Tests for mark κ̃ calculation."""

    def test_kappa_tilde_mark_formula(self, small_counts: np.ndarray):
        """κ̃_ik = y_ik - N_i/2 for marks."""
        from bayesian_statistics.nngp.model.marked_point_process.composition import (
            compute_kappa_tilde_mark,
        )

        # For category k
        counts_k = small_counts[:, 0]  # First category counts
        total_counts = small_counts.sum(axis=1)

        kappa_tilde = compute_kappa_tilde_mark(counts_k, total_counts)
        expected = counts_k - total_counts / 2

        np.testing.assert_array_almost_equal(kappa_tilde, expected)

    def test_kappa_tilde_mark_shape(self, small_counts: np.ndarray):
        """κ̃ should have same shape as counts_k."""
        from bayesian_statistics.nngp.model.marked_point_process.composition import (
            compute_kappa_tilde_mark,
        )

        counts_k = small_counts[:, 0]
        total_counts = small_counts.sum(axis=1)

        kappa_tilde = compute_kappa_tilde_mark(counts_k, total_counts)
        assert kappa_tilde.shape == counts_k.shape


class TestSampleXiMark:
    """Tests for PG sampling for marks."""

    def test_sample_xi_mark_uses_b_equals_N(
        self, small_counts: np.ndarray, rng: np.random.Generator
    ):
        """Mark PG uses b = N_i (total count)."""
        from bayesian_statistics.nngp.model.marked_point_process.composition import (
            sample_xi_mark,
        )

        n_sites = small_counts.shape[0]
        total_counts = small_counts.sum(axis=1)
        eta_k = np.zeros(n_sites)  # Linear predictor for one category

        xi_k = sample_xi_mark(eta_k, total_counts, rng)

        assert xi_k.shape == (n_sites,)
        assert np.all(xi_k > 0)  # PG is always positive

    def test_sample_xi_mark_larger_for_larger_N(self):
        """E[PG(N, c)] increases with N."""
        from bayesian_statistics.nngp.model.marked_point_process.composition import (
            sample_xi_mark,
        )

        n = 20
        eta = np.zeros(n)
        N_small = np.ones(n) * 5
        N_large = np.ones(n) * 100

        n_samples = 50
        xi_small_mean = []
        xi_large_mean = []

        for seed in range(n_samples):
            rng = np.random.default_rng(seed)
            xi_small_mean.append(sample_xi_mark(eta, N_small, rng).mean())

            rng = np.random.default_rng(seed + 1000)
            xi_large_mean.append(sample_xi_mark(eta, N_large, rng).mean())

        # Larger N should give larger expected PG value
        assert np.mean(xi_large_mean) > np.mean(xi_small_mean)


class TestSoftmaxWithBaseline:
    """Tests for softmax computation."""

    def test_softmax_sums_to_one(self):
        """Probabilities should sum to 1 at each site."""
        from bayesian_statistics.nngp.model.marked_point_process.composition import (
            softmax_with_baseline,
        )

        K_minus_1 = 2
        n_sites = 10
        eta = np.random.randn(K_minus_1, n_sites)

        probs = softmax_with_baseline(eta)

        assert probs.shape == (K_minus_1 + 1, n_sites)
        np.testing.assert_array_almost_equal(probs.sum(axis=0), np.ones(n_sites))

    def test_softmax_all_positive(self):
        """All probabilities should be positive."""
        from bayesian_statistics.nngp.model.marked_point_process.composition import (
            softmax_with_baseline,
        )

        eta = np.array([[10.0, -10.0, 0.0], [-5.0, 5.0, 0.0]])
        probs = softmax_with_baseline(eta)

        assert np.all(probs > 0)

    def test_softmax_baseline_formula(self):
        """Baseline category should be exp(0) normalized."""
        from bayesian_statistics.nngp.model.marked_point_process.composition import (
            softmax_with_baseline,
        )

        eta = np.array([[0.0], [0.0]])  # Two non-baseline categories
        probs = softmax_with_baseline(eta)

        # With eta = 0 for all, all three categories should be equal
        expected = 1.0 / 3.0
        np.testing.assert_array_almost_equal(probs, np.full((3, 1), expected))


class TestCompositionMatchesExisting:
    """Regression tests against existing multinomial model."""

    def test_softmax_matches_sample_py(self):
        """Our softmax should match sample.py softmax_with_baseline."""
        from bayesian_statistics.nngp.model.marked_point_process.composition import (
            softmax_with_baseline,
        )
        from bayesian_statistics.nngp.model.sample import (
            softmax_with_baseline as original_softmax,
        )

        np.random.seed(42)
        eta = np.random.randn(3, 20)

        our_probs = softmax_with_baseline(eta)
        original_probs = original_softmax(eta)

        np.testing.assert_array_almost_equal(our_probs, original_probs)

    def test_compute_eta_matches_sample_py(self):
        """Our compute_eta should match sample.py compute_eta."""
        from bayesian_statistics.nngp.model.marked_point_process.composition import (
            compute_eta,
        )
        from bayesian_statistics.nngp.model.sample import (
            compute_eta as original_compute_eta,
        )

        K_minus_1 = 3
        p_plus_1 = 2
        n = 10

        np.random.seed(42)
        beta = np.random.randn(K_minus_1, p_plus_1, n)
        W = np.random.randn(n, p_plus_1)
        W[:, 0] = 1.0  # Intercept

        our_eta = compute_eta(beta, W)
        original_eta = original_compute_eta(beta, W)

        np.testing.assert_array_almost_equal(our_eta, original_eta)
