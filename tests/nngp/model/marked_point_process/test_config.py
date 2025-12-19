"""Tests for MarkedPointProcessConfig."""

from __future__ import annotations

import pytest


class TestMarkedPointProcessConfig:
    """Tests for the configuration dataclass."""

    def test_config_has_mcmc_parameters(self):
        """Config should have n_iter, burn_in, thinning, seed."""
        from bayesian_statistics.nngp.model.marked_point_process.config import (
            MarkedPointProcessConfig,
        )

        config = MarkedPointProcessConfig()
        assert hasattr(config, "n_iter")
        assert hasattr(config, "burn_in")
        assert hasattr(config, "thinning")
        assert hasattr(config, "seed")

    def test_config_has_nngp_parameters(self):
        """Config should have neighbor_count and kernel parameters."""
        from bayesian_statistics.nngp.model.marked_point_process.config import (
            MarkedPointProcessConfig,
        )

        config = MarkedPointProcessConfig()
        assert hasattr(config, "neighbor_count")
        assert hasattr(config, "intensity_kernel_lengthscale")
        assert hasattr(config, "intensity_kernel_variance")
        assert hasattr(config, "mark_kernel_lengthscale")
        assert hasattr(config, "mark_kernel_variance")

    def test_config_has_lambda_prior(self):
        """Config should have lambda* prior parameters (Gamma)."""
        from bayesian_statistics.nngp.model.marked_point_process.config import (
            MarkedPointProcessConfig,
        )

        config = MarkedPointProcessConfig()
        assert hasattr(config, "lambda_prior_shape")
        assert hasattr(config, "lambda_prior_rate")

    def test_n_saved_calculation(self):
        """n_saved should correctly compute number of saved samples."""
        from bayesian_statistics.nngp.model.marked_point_process.config import (
            MarkedPointProcessConfig,
        )

        config = MarkedPointProcessConfig(n_iter=1000, burn_in=200, thinning=2)
        # (1000 - 200) / 2 = 400
        assert config.n_saved() == 400

    def test_n_saved_with_zero_thinning_raises(self):
        """n_saved should raise error if thinning <= 0."""
        from bayesian_statistics.nngp.model.marked_point_process.config import (
            MarkedPointProcessConfig,
        )

        config = MarkedPointProcessConfig(thinning=0)
        with pytest.raises(ValueError):
            config.n_saved()

    def test_default_values_are_sensible(self):
        """Default values should be sensible for typical use."""
        from bayesian_statistics.nngp.model.marked_point_process.config import (
            MarkedPointProcessConfig,
        )

        config = MarkedPointProcessConfig()
        assert config.n_iter > 0
        assert config.burn_in >= 0
        assert config.thinning >= 1
        assert config.neighbor_count > 0
        assert config.intensity_kernel_lengthscale > 0
        assert config.intensity_kernel_variance > 0
        assert config.lambda_prior_shape > 0
        assert config.lambda_prior_rate > 0
