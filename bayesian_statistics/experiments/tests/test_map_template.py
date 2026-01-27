"""Tests for map template (MapPlotter)."""

import tempfile
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def sample_data():
    """Create sample data for testing."""
    np.random.seed(42)
    n_grid = 100
    n_boundary = 20
    n_sites = 10

    # Grid coordinates (2D)
    grid_coords = np.random.rand(n_grid, 2)

    # Boundary points
    boundary = np.random.rand(n_boundary, 2)

    # Land mask
    is_land = np.random.choice([True, False], size=n_grid, p=[0.7, 0.3])

    # Site coordinates
    site_coords = np.random.rand(n_sites, 2)

    # Values for coloring
    grid_values = np.random.rand(n_grid)
    site_values = np.random.rand(n_sites)

    return {
        "grid_coords": grid_coords,
        "boundary": boundary,
        "is_land": is_land,
        "site_coords": site_coords,
        "grid_values": grid_values,
        "site_values": site_values,
    }


def test_map_plotter_init(sample_data):
    """MapPlotterが正しく初期化されることを確認"""
    from bayesian_statistics.experiments.visualization.map_template import MapPlotter

    plotter = MapPlotter(
        boundary=sample_data["boundary"],
        is_land=sample_data["is_land"],
    )

    assert plotter.boundary.shape == sample_data["boundary"].shape
    assert plotter.is_land.shape == sample_data["is_land"].shape


def test_create_figure(sample_data):
    """figure/axesが正しいサイズで作成されることを確認"""
    from bayesian_statistics.experiments.visualization.map_template import MapPlotter

    plotter = MapPlotter(
        boundary=sample_data["boundary"],
        is_land=sample_data["is_land"],
    )

    # Single plot
    fig, ax = plotter.create_figure(1, 1, figsize=(10, 8))
    assert fig is not None
    assert ax is not None
    import matplotlib.pyplot as plt

    plt.close(fig)

    # 2x2 grid
    fig, axes = plotter.create_figure(2, 2, figsize=(12, 10))
    assert fig is not None
    assert axes.shape == (2, 2)
    plt.close(fig)


def test_plot_grid_scatter(sample_data):
    """グリッド散布図が正しく描画されることを確認"""
    from bayesian_statistics.experiments.visualization.map_template import MapPlotter
    import matplotlib.pyplot as plt

    plotter = MapPlotter(
        boundary=sample_data["boundary"],
        is_land=sample_data["is_land"],
    )

    fig, ax = plotter.create_figure(1, 1)

    # Plot grid scatter
    sc = plotter.plot_grid_scatter(
        ax=ax,
        coords=sample_data["grid_coords"],
        values=sample_data["grid_values"],
    )

    assert sc is not None
    plt.close(fig)


def test_plot_site_scatter(sample_data):
    """遺跡点の散布図が正しく描画されることを確認"""
    from bayesian_statistics.experiments.visualization.map_template import MapPlotter
    import matplotlib.pyplot as plt

    plotter = MapPlotter(
        boundary=sample_data["boundary"],
        is_land=sample_data["is_land"],
    )

    fig, ax = plotter.create_figure(1, 1)

    # Plot site scatter
    sc = plotter.plot_site_scatter(
        ax=ax,
        coords=sample_data["site_coords"],
        values=sample_data["site_values"],
    )

    assert sc is not None
    plt.close(fig)


def test_plot_boundary(sample_data):
    """境界線が正しく描画されることを確認"""
    from bayesian_statistics.experiments.visualization.map_template import MapPlotter
    import matplotlib.pyplot as plt

    plotter = MapPlotter(
        boundary=sample_data["boundary"],
        is_land=sample_data["is_land"],
    )

    fig, ax = plotter.create_figure(1, 1)

    # Plot boundary (should not raise)
    plotter.plot_boundary(ax)

    plt.close(fig)


def test_save_figure(sample_data, tmp_path):
    """図がPNG形式で保存されることを確認"""
    from bayesian_statistics.experiments.visualization.map_template import MapPlotter
    import matplotlib.pyplot as plt

    plotter = MapPlotter(
        boundary=sample_data["boundary"],
        is_land=sample_data["is_land"],
    )

    fig, ax = plotter.create_figure(1, 1)
    plotter.plot_grid_scatter(
        ax=ax,
        coords=sample_data["grid_coords"],
        values=sample_data["grid_values"],
    )

    # Save to temp directory
    output_path = tmp_path / "test_figure.png"
    plotter.save_figure(fig, output_path)

    assert output_path.exists()
    assert output_path.stat().st_size > 0

    plt.close(fig)


def test_add_colorbar(sample_data):
    """カラーバーが正しく追加されることを確認"""
    from bayesian_statistics.experiments.visualization.map_template import MapPlotter
    import matplotlib.pyplot as plt

    plotter = MapPlotter(
        boundary=sample_data["boundary"],
        is_land=sample_data["is_land"],
    )

    fig, ax = plotter.create_figure(1, 1)
    sc = plotter.plot_grid_scatter(
        ax=ax,
        coords=sample_data["grid_coords"],
        values=sample_data["grid_values"],
    )

    # Add colorbar (should not raise)
    cbar = plotter.add_colorbar(fig, sc, ax, label="Test Label")
    assert cbar is not None

    plt.close(fig)


def test_set_axis_style(sample_data):
    """軸のスタイル設定が正しく適用されることを確認"""
    from bayesian_statistics.experiments.visualization.map_template import MapPlotter
    import matplotlib.pyplot as plt

    plotter = MapPlotter(
        boundary=sample_data["boundary"],
        is_land=sample_data["is_land"],
    )

    fig, ax = plotter.create_figure(1, 1)

    # Set axis style
    plotter.set_axis_style(ax, title="Test Title", xlabel="X", ylabel="Y")

    assert ax.get_title() == "Test Title"
    assert ax.get_xlabel() == "X"
    assert ax.get_ylabel() == "Y"

    plt.close(fig)


def test_plot_with_mask(sample_data):
    """陸地マスクが正しく適用されることを確認"""
    from bayesian_statistics.experiments.visualization.map_template import MapPlotter
    import matplotlib.pyplot as plt

    plotter = MapPlotter(
        boundary=sample_data["boundary"],
        is_land=sample_data["is_land"],
    )

    fig, ax = plotter.create_figure(1, 1)

    # Plot with mask
    sc = plotter.plot_grid_scatter(
        ax=ax,
        coords=sample_data["grid_coords"],
        values=sample_data["grid_values"],
        mask=sample_data["is_land"],
    )

    assert sc is not None
    plt.close(fig)
