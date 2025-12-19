"""Common fixtures for marked point process tests."""

from __future__ import annotations

from typing import List

import numpy as np
import pytest


@pytest.fixture
def rng() -> np.random.Generator:
    """Reproducible random number generator."""
    return np.random.default_rng(42)


@pytest.fixture
def small_coords() -> np.ndarray:
    """Small set of site coordinates for testing (10 sites)."""
    return np.array([
        [139.0, 35.0],
        [139.1, 35.1],
        [139.2, 35.0],
        [139.0, 35.2],
        [139.3, 35.1],
        [138.9, 35.3],
        [139.1, 35.2],
        [139.2, 35.3],
        [138.8, 35.1],
        [139.4, 35.2],
    ])


@pytest.fixture
def small_counts() -> np.ndarray:
    """Small counts matrix for testing (10 sites, 3 categories)."""
    return np.array([
        [5, 3, 2],
        [8, 1, 1],
        [2, 6, 2],
        [4, 4, 2],
        [7, 2, 1],
        [1, 7, 2],
        [3, 3, 4],
        [6, 2, 2],
        [2, 5, 3],
        [4, 3, 3],
    ], dtype=float)


@pytest.fixture
def small_grid_coords() -> np.ndarray:
    """Small grid coordinates for testing (25 points, 5x5 grid)."""
    x = np.linspace(138.7, 139.5, 5)
    y = np.linspace(34.9, 35.4, 5)
    xx, yy = np.meshgrid(x, y)
    return np.column_stack([xx.ravel(), yy.ravel()])


@pytest.fixture
def small_design_matrix(small_coords: np.ndarray) -> np.ndarray:
    """Small design matrix with intercept (10 sites, 2 features)."""
    n = small_coords.shape[0]
    # Intercept + 1 covariate (simulated distance)
    return np.column_stack([
        np.ones(n),
        np.random.default_rng(42).uniform(0, 1, n),
    ])


@pytest.fixture
def origins() -> List[str]:
    """Category names for testing."""
    return ["Origin_A", "Origin_B", "Baseline"]


@pytest.fixture
def region() -> List[List[float]]:
    """Spatial region for testing."""
    return [[138.5, 139.6], [34.8, 35.5]]
