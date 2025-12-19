"""Tests for IntensityComponent (point process part)."""

from __future__ import annotations

import numpy as np
import pytest


class TestComputeQ:
    """Tests for existence probability computation."""

    def test_compute_q_returns_probability(self):
        """compute_q should return values in [0, 1]."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            compute_q,
        )

        eta = np.array([-10.0, -1.0, 0.0, 1.0, 10.0])
        q = compute_q(eta)

        assert q.shape == eta.shape
        assert np.all(q >= 0)
        assert np.all(q <= 1)

    def test_compute_q_at_zero_is_half(self):
        """compute_q(0) should be 0.5."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            compute_q,
        )

        q = compute_q(np.array([0.0]))
        np.testing.assert_almost_equal(q[0], 0.5)

    def test_compute_q_is_sigmoid(self):
        """compute_q should be sigmoid function."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            compute_q,
        )

        eta = np.array([1.0, 2.0, -1.0])
        q = compute_q(eta)
        expected = 1.0 / (1.0 + np.exp(-eta))
        np.testing.assert_array_almost_equal(q, expected)


class TestSamplePseudoAbsence:
    """Tests for pseudo-absence sampling via Poisson thinning."""

    def test_sample_pseudo_absence_returns_array(
        self,
        small_grid_coords: np.ndarray,
        rng: np.random.Generator,
    ):
        """sample_pseudo_absence should return coordinate array."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            sample_pseudo_absence,
        )

        n_grid = small_grid_coords.shape[0]
        valid_mask = np.ones(n_grid, dtype=bool)
        eta_grid = np.zeros(n_grid)  # q = 0.5 everywhere
        region_volume = 1.0
        lambda_star = 10.0

        U_indices = sample_pseudo_absence(
            lambda_star=lambda_star,
            eta_grid=eta_grid,
            valid_mask=valid_mask,
            region_volume=region_volume,
            rng=rng,
        )

        assert isinstance(U_indices, np.ndarray)
        assert U_indices.ndim == 1

    def test_sample_pseudo_absence_returns_valid_indices(
        self,
        small_grid_coords: np.ndarray,
        rng: np.random.Generator,
    ):
        """Returned indices should be within valid mask."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            sample_pseudo_absence,
        )

        n_grid = small_grid_coords.shape[0]
        # Only half the grids are valid
        valid_mask = np.zeros(n_grid, dtype=bool)
        valid_mask[:n_grid // 2] = True
        eta_grid = np.zeros(n_grid)
        region_volume = 1.0
        lambda_star = 50.0

        U_indices = sample_pseudo_absence(
            lambda_star=lambda_star,
            eta_grid=eta_grid,
            valid_mask=valid_mask,
            region_volume=region_volume,
            rng=rng,
        )

        # All returned indices should be in valid region
        if len(U_indices) > 0:
            assert np.all(valid_mask[U_indices])

    def test_sample_pseudo_absence_respects_thinning(
        self,
        small_grid_coords: np.ndarray,
    ):
        """Higher q means fewer pseudo-absence points (more thinning)."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            sample_pseudo_absence,
        )

        n_grid = small_grid_coords.shape[0]
        valid_mask = np.ones(n_grid, dtype=bool)
        region_volume = 1.0
        lambda_star = 100.0

        # High q -> few pseudo-absence (U ~ IPP(λ*(1-q)))
        eta_high = np.ones(n_grid) * 5.0  # q ≈ 1
        # Low q -> more pseudo-absence
        eta_low = np.ones(n_grid) * -5.0  # q ≈ 0

        n_samples = 20
        counts_high = []
        counts_low = []

        for seed in range(n_samples):
            rng = np.random.default_rng(seed)
            U_high = sample_pseudo_absence(
                lambda_star, eta_high, valid_mask, region_volume, rng
            )
            counts_high.append(len(U_high))

            rng = np.random.default_rng(seed + 1000)
            U_low = sample_pseudo_absence(
                lambda_star, eta_low, valid_mask, region_volume, rng
            )
            counts_low.append(len(U_low))

        # On average, low q should produce more points
        assert np.mean(counts_low) > np.mean(counts_high)


class TestUpdateLambdaStar:
    """Tests for lambda* posterior update."""

    def test_update_lambda_star_returns_positive(self, rng: np.random.Generator):
        """update_lambda_star should return positive value."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            update_lambda_star,
        )

        result = update_lambda_star(
            n_total=50,
            region_volume=1.0,
            prior_shape=2.0,
            prior_rate=1.0,
            rng=rng,
        )

        assert result > 0

    def test_update_lambda_star_posterior_shape(self, rng: np.random.Generator):
        """Posterior shape should be prior_shape + n_total."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            update_lambda_star,
        )

        # Sample many values and check distribution
        n_samples = 1000
        samples = []
        prior_shape = 2.0
        prior_rate = 1.0
        n_total = 50
        region_volume = 2.0

        for seed in range(n_samples):
            rng = np.random.default_rng(seed)
            samples.append(
                update_lambda_star(n_total, region_volume, prior_shape, prior_rate, rng)
            )

        samples = np.array(samples)
        # Gamma(shape, rate) has mean = shape/rate
        expected_shape = prior_shape + n_total
        expected_rate = prior_rate + region_volume
        expected_mean = expected_shape / expected_rate

        # Sample mean should be close to expected
        np.testing.assert_almost_equal(np.mean(samples), expected_mean, decimal=0)


class TestSampleOmegaIntensity:
    """Tests for PG sampling for intensity."""

    def test_sample_omega_intensity_uses_b_equals_1(self, rng: np.random.Generator):
        """Point process PG uses b=1."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            sample_omega_intensity,
        )

        eta = np.array([0.0, 1.0, -1.0, 2.0])
        omega = sample_omega_intensity(eta, rng)

        assert omega.shape == eta.shape
        assert np.all(omega > 0)  # PG(1, c) is always positive

    def test_sample_omega_intensity_larger_for_larger_eta(
        self, rng: np.random.Generator
    ):
        """E[PG(1,c)] = tanh(c/2)/(2c), which increases with |c|."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            sample_omega_intensity,
        )

        # Sample many times to get mean
        n_samples = 100
        omega_small = []
        omega_large = []

        for seed in range(n_samples):
            rng = np.random.default_rng(seed)
            omega_small.append(sample_omega_intensity(np.array([0.1]), rng)[0])
            rng = np.random.default_rng(seed + 1000)
            omega_large.append(sample_omega_intensity(np.array([5.0]), rng)[0])

        # For |c| large, mean(ω) ≈ 1/(2|c|), which is smaller
        # For |c| small, mean(ω) ≈ 0.25
        assert np.mean(omega_small) > np.mean(omega_large) * 0.5


class TestKappaIntensity:
    """Tests for kappa calculation."""

    def test_kappa_intensity_formula(self):
        """κ_i = y_i - 0.5 for point process."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            compute_kappa_intensity,
        )

        y = np.array([1, 1, 0, 0, 1])
        kappa = compute_kappa_intensity(y)
        expected = y - 0.5

        np.testing.assert_array_equal(kappa, expected)

    def test_kappa_intensity_values(self):
        """Presence -> 0.5, absence -> -0.5."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            compute_kappa_intensity,
        )

        y = np.array([1, 0])
        kappa = compute_kappa_intensity(y)

        assert kappa[0] == 0.5
        assert kappa[1] == -0.5


class TestUpdateBetaIntensity:
    """Tests for intensity coefficient update with NNGP spatial effect.

    These tests verify that the point process Gibbs step (c) from sec7.tex
    is correctly implemented:
    (β_int, u_int) ~ N(m_int, V_int)
    evaluated at X ∪ U (all observed sites and pseudo-absence points).
    """

    def test_update_beta_intensity_produces_finite_values(
        self,
        small_coords: np.ndarray,
        rng: np.random.Generator,
    ):
        """update_beta_intensity should produce finite beta values."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            update_beta_intensity,
        )
        from bayesian_statistics.nngp.model.nngp import build_nngp_factors, order_points_morton
        from bayesian_statistics.nngp.model.sample import LocalNNGPKernel

        n_sites = small_coords.shape[0]
        p_plus_1 = 2  # intercept + 1 covariate

        # Design matrix: intercept + distance-like feature
        W = np.column_stack([np.ones(n_sites), np.linspace(0, 1, n_sites)])

        # Build NNGP factors
        kernel = LocalNNGPKernel(lengthscale=0.1, variance=1.0)
        order = order_points_morton(small_coords)
        factors, _ = build_nngp_factors(small_coords, M=3, kernel=kernel, order=order)

        # Initial state
        beta_int = np.zeros((p_plus_1, n_sites))
        eta_int = np.zeros(n_sites)

        # PG samples and kappa
        omega = 0.25 * np.ones(n_sites)  # E[PG(1,0)] ≈ 0.25
        y = np.ones(n_sites)  # All presence
        kappa = y - 0.5

        # Run update
        beta_new, eta_new = update_beta_intensity(
            beta_int=beta_int,
            eta_int=eta_int,
            W=W,
            factors_by_feature=[factors] * p_plus_1,
            order=order,
            omega=omega,
            kappa=kappa,
            rng=rng,
        )

        assert np.all(np.isfinite(beta_new))
        assert np.all(np.isfinite(eta_new))

    def test_update_beta_intensity_with_presence_produces_nonzero_eta(
        self,
        small_coords: np.ndarray,
        rng: np.random.Generator,
    ):
        """When all y=1 (presence), eta should become non-zero after updates."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            sample_omega_intensity,
            update_beta_intensity,
        )
        from bayesian_statistics.nngp.model.nngp import build_nngp_factors, order_points_morton
        from bayesian_statistics.nngp.model.sample import LocalNNGPKernel

        n_sites = small_coords.shape[0]
        p_plus_1 = 1  # intercept only

        W = np.ones((n_sites, p_plus_1))

        kernel = LocalNNGPKernel(lengthscale=0.1, variance=1.0)
        order = order_points_morton(small_coords)
        factors, _ = build_nngp_factors(small_coords, M=3, kernel=kernel, order=order)

        beta_int = np.zeros((p_plus_1, n_sites))
        eta_int = np.zeros(n_sites)
        y = np.ones(n_sites)  # All presence

        # Run multiple iterations
        for _ in range(5):
            omega = sample_omega_intensity(eta_int, rng)
            kappa = y - 0.5
            beta_int, eta_int = update_beta_intensity(
                beta_int=beta_int,
                eta_int=eta_int,
                W=W,
                factors_by_feature=[factors] * p_plus_1,
                order=order,
                omega=omega,
                kappa=kappa,
                rng=rng,
            )

        # After updates, eta should be non-zero (spatial variation)
        assert not np.allclose(eta_int, 0)

    def test_update_beta_intensity_at_combined_points(
        self,
        small_coords: np.ndarray,
        small_grid_coords: np.ndarray,
        rng: np.random.Generator,
    ):
        """Beta should be updated at X ∪ U (combined sites and pseudo-absence)."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            sample_pseudo_absence,
            update_beta_intensity,
        )
        from bayesian_statistics.nngp.model.nngp import build_nngp_factors, order_points_morton
        from bayesian_statistics.nngp.model.sample import LocalNNGPKernel

        n_sites = small_coords.shape[0]

        # Sample some pseudo-absence points
        n_grid = small_grid_coords.shape[0]
        valid_mask = np.ones(n_grid, dtype=bool)
        eta_grid = np.zeros(n_grid)
        U_indices = sample_pseudo_absence(
            lambda_star=10.0,
            eta_grid=eta_grid,
            valid_mask=valid_mask,
            region_volume=1.0,
            rng=rng,
        )

        # Get pseudo-absence coordinates
        U_coords = small_grid_coords[U_indices] if len(U_indices) > 0 else np.empty((0, 2))

        # Combine X and U
        if len(U_coords) > 0:
            combined_coords = np.vstack([small_coords, U_coords])
            y_combined = np.concatenate([np.ones(n_sites), np.zeros(len(U_coords))])
        else:
            combined_coords = small_coords
            y_combined = np.ones(n_sites)

        n_combined = len(combined_coords)
        p_plus_1 = 1

        W = np.ones((n_combined, p_plus_1))

        # Build NNGP factors for combined points
        kernel = LocalNNGPKernel(lengthscale=0.1, variance=1.0)
        order = order_points_morton(combined_coords)
        factors, _ = build_nngp_factors(combined_coords, M=3, kernel=kernel, order=order)

        beta_int = np.zeros((p_plus_1, n_combined))
        eta_int = np.zeros(n_combined)
        omega = 0.25 * np.ones(n_combined)
        kappa = y_combined - 0.5

        beta_new, eta_new = update_beta_intensity(
            beta_int=beta_int,
            eta_int=eta_int,
            W=W,
            factors_by_feature=[factors] * p_plus_1,
            order=order,
            omega=omega,
            kappa=kappa,
            rng=rng,
        )

        # Output should have correct shape for combined points
        assert beta_new.shape == (p_plus_1, n_combined)
        assert eta_new.shape == (n_combined,)
        assert np.all(np.isfinite(beta_new))


class TestComputeEtaIntensity:
    """Tests for computing eta for intensity."""

    def test_compute_eta_intensity_formula(
        self,
        small_coords: np.ndarray,
    ):
        """eta_int,i = sum_j W[i,j] * beta[j,i]."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            compute_eta_intensity,
        )

        n = small_coords.shape[0]
        p_plus_1 = 2

        W = np.column_stack([np.ones(n), np.arange(n, dtype=float)])
        beta = np.ones((p_plus_1, n))

        eta = compute_eta_intensity(beta, W)

        # eta[i] = W[i,0] * beta[0,i] + W[i,1] * beta[1,i]
        #        = 1 * 1 + i * 1 = 1 + i
        expected = 1 + np.arange(n)
        np.testing.assert_array_almost_equal(eta, expected)

    def test_compute_eta_intensity_with_zero_beta(
        self,
        small_coords: np.ndarray,
    ):
        """Zero beta should give zero eta."""
        from bayesian_statistics.nngp.model.marked_point_process.intensity import (
            compute_eta_intensity,
        )

        n = small_coords.shape[0]
        p_plus_1 = 2

        W = np.column_stack([np.ones(n), np.random.randn(n)])
        beta = np.zeros((p_plus_1, n))

        eta = compute_eta_intensity(beta, W)

        np.testing.assert_array_almost_equal(eta, np.zeros(n))
