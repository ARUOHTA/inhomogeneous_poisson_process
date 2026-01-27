"""Effect decomposition visualization."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Union

import numpy as np

from .map_template import MapPlotter
from .style import (
    FIGURE_SIZE_GRID_2X2,
    TITLE_FONTSIZE,
)

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


def plot_effect_comparison(
    plotter: MapPlotter,
    grid_coords: np.ndarray,
    effects: Dict[str, np.ndarray],
    origin_index: int,
    origin_name: str,
    effect_names: List[str],
    output_path: Optional[Union[str, Path]] = None,
) -> tuple:
    """Plot comparison of multiple effects for a single origin.

    Parameters
    ----------
    plotter : MapPlotter
        Map plotter instance.
    grid_coords : np.ndarray
        Grid coordinates (n_grid, 2).
    effects : dict
        Dictionary of effect arrays, each (K, n_grid).
    origin_index : int
        Index of the origin to plot.
    origin_name : str
        Name of the origin.
    effect_names : list of str
        Names of effects to plot.
    output_path : str or Path, optional
        Path to save the figure.

    Returns
    -------
    fig : Figure
        Matplotlib figure.
    axes : list
        List of matplotlib axes.
    """
    n_effects = len(effect_names)
    fig, axes = plotter.create_figure(1, n_effects, figsize=(6 * n_effects, 5))
    if n_effects == 1:
        axes = [axes]

    for ax, effect_name in zip(axes, effect_names):
        effect_data = effects[effect_name]

        sc = plotter.plot_grid_scatter(
            ax=ax,
            coords=grid_coords,
            values=effect_data[origin_index],
        )
        plotter.plot_boundary(ax)
        plotter.add_colorbar(fig, sc, ax, label="確率")
        plotter.set_axis_style(
            ax, title=f"{origin_name}: {effect_name}", xlabel="経度", ylabel="緯度"
        )

    if output_path:
        plotter.save_figure(fig, output_path)

    return fig, axes


def plot_all_origins_single_effect(
    plotter: MapPlotter,
    grid_coords: np.ndarray,
    effect_data: np.ndarray,
    effect_name: str,
    origins: List[str],
    output_path: Optional[Union[str, Path]] = None,
) -> tuple:
    """Plot a single effect for all origins in a 2x2 grid.

    Parameters
    ----------
    plotter : MapPlotter
        Map plotter instance.
    grid_coords : np.ndarray
        Grid coordinates (n_grid, 2).
    effect_data : np.ndarray
        Effect data (K, n_grid).
    effect_name : str
        Name of the effect.
    origins : list of str
        Names of origins (first 4 used).
    output_path : str or Path, optional
        Path to save the figure.

    Returns
    -------
    fig : Figure
        Matplotlib figure.
    axes : np.ndarray
        Array of matplotlib axes.
    """
    fig, axes = plotter.create_figure(2, 2, figsize=FIGURE_SIZE_GRID_2X2)
    axes_flat = axes.flatten()

    for idx in range(4):  # 4 main origins
        ax = axes_flat[idx]

        sc = plotter.plot_grid_scatter(
            ax=ax,
            coords=grid_coords,
            values=effect_data[idx],
        )
        plotter.plot_boundary(ax)
        plotter.add_colorbar(fig, sc, ax, label="確率")
        plotter.set_axis_style(ax, title=origins[idx])

    fig.suptitle(f"効果: {effect_name}", fontsize=TITLE_FONTSIZE)

    if output_path:
        plotter.save_figure(fig, output_path)

    return fig, axes


def plot_effect_decomposition_grid(
    plotter: MapPlotter,
    grid_coords: np.ndarray,
    effects: Dict[str, np.ndarray],
    origins: List[str],
    effect_names: List[str] = None,
    output_path: Optional[Union[str, Path]] = None,
) -> tuple:
    """Plot effect decomposition as origins x effects grid.

    Parameters
    ----------
    plotter : MapPlotter
        Map plotter instance.
    grid_coords : np.ndarray
        Grid coordinates (n_grid, 2).
    effects : dict
        Dictionary of effect arrays, each (K, n_grid).
    origins : list of str
        Names of origins.
    effect_names : list of str, optional
        Names of effects to plot. Default: ["distance", "intercept_adjustment", "full"]
    output_path : str or Path, optional
        Path to save the figure.

    Returns
    -------
    fig : Figure
        Matplotlib figure.
    axes : np.ndarray
        Array of matplotlib axes.
    """
    if effect_names is None:
        effect_names = ["distance", "intercept_adjustment", "full"]

    n_origins = min(4, len(origins))  # Use first 4 origins
    n_effects = len(effect_names)

    fig, axes = plotter.create_figure(
        n_origins,
        n_effects,
        figsize=(5 * n_effects, 5 * n_origins),
    )

    for origin_idx in range(n_origins):
        for effect_idx, effect_name in enumerate(effect_names):
            ax = axes[origin_idx, effect_idx]
            effect_data = effects[effect_name]

            sc = plotter.plot_grid_scatter(
                ax=ax,
                coords=grid_coords,
                values=effect_data[origin_idx],
                s=8,
            )
            plotter.plot_boundary(ax)

            # Labels
            if origin_idx == 0:
                ax.set_title(effect_name, fontsize=11)
            if effect_idx == 0:
                ax.set_ylabel(origins[origin_idx], fontsize=11)

            ax.set_xticks([])
            ax.set_yticks([])

    # Common colorbar
    fig.colorbar(sc, ax=axes, label="確率", shrink=0.7, pad=0.02)
    fig.suptitle("効果の分解", fontsize=TITLE_FONTSIZE, y=1.01)

    if output_path:
        plotter.save_figure(fig, output_path)

    return fig, axes
