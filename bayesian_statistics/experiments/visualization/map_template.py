"""Map plot template for consistent visualization."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.figure import Figure

from .style import (
    BOUNDARY_COLOR,
    BOUNDARY_SIZE,
    CMAP_PROBABILITY,
    FIGURE_DPI,
    FIGURE_FORMAT,
    FIGURE_SIZE_SINGLE,
    GRID_SCATTER_ALPHA,
    GRID_SCATTER_SIZE,
    PROB_VMAX,
    PROB_VMIN,
    SITE_EDGE_WIDTH,
    SITE_SCATTER_ALPHA,
    SITE_SCATTER_SIZE,
    SPINE_WIDTH,
    apply_japanese_font,
    set_spine_width,
)


class MapPlotter:
    """Map plotter template for consistent visualization.

    This class provides a common interface for creating map plots with
    consistent styling, supporting grid scatter plots, site markers,
    and boundary visualization.

    Attributes
    ----------
    boundary : np.ndarray
        Boundary points for the study area (n_points, 2).
    is_land : np.ndarray
        Boolean mask indicating land cells (n_grid,).
    """

    def __init__(
        self,
        boundary: np.ndarray,
        is_land: np.ndarray,
    ):
        """Initialize MapPlotter.

        Parameters
        ----------
        boundary : np.ndarray
            Boundary points for the study area (n_points, 2).
        is_land : np.ndarray
            Boolean mask indicating land cells (n_grid,).
        """
        self.boundary = boundary
        self.is_land = is_land

        # Apply Japanese font settings
        try:
            apply_japanese_font()
        except ImportError:
            pass  # japanize_matplotlib not installed

    def create_figure(
        self,
        nrows: int = 1,
        ncols: int = 1,
        figsize: Optional[Tuple[float, float]] = None,
        **kwargs: Any,
    ) -> Tuple[Figure, Union[Axes, np.ndarray]]:
        """Create figure and axes with consistent styling.

        Parameters
        ----------
        nrows : int
            Number of rows in the subplot grid.
        ncols : int
            Number of columns in the subplot grid.
        figsize : tuple, optional
            Figure size (width, height) in inches.
        **kwargs
            Additional keyword arguments passed to plt.subplots.

        Returns
        -------
        fig : Figure
            Matplotlib figure object.
        axes : Axes or np.ndarray
            Matplotlib axes object(s).
        """
        if figsize is None:
            figsize = FIGURE_SIZE_SINGLE

        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=figsize,
            constrained_layout=True,
            squeeze=False if (nrows > 1 or ncols > 1) else True,
            **kwargs,
        )

        # For single plot, return single axes
        if nrows == 1 and ncols == 1:
            return fig, axes

        return fig, axes

    def plot_grid_scatter(
        self,
        ax: Axes,
        coords: np.ndarray,
        values: np.ndarray,
        cmap: str = CMAP_PROBABILITY,
        vmin: float = PROB_VMIN,
        vmax: float = PROB_VMAX,
        norm: Optional[Normalize] = None,
        mask: Optional[np.ndarray] = None,
        s: float = GRID_SCATTER_SIZE,
        alpha: float = GRID_SCATTER_ALPHA,
    ):
        """Plot grid scatter points.

        Parameters
        ----------
        ax : Axes
            Matplotlib axes object.
        coords : np.ndarray
            Coordinates (n_points, 2).
        values : np.ndarray
            Values for coloring (n_points,).
        cmap : str
            Colormap name.
        vmin : float
            Minimum value for colormap.
        vmax : float
            Maximum value for colormap.
        norm : Normalize, optional
            Normalization object (e.g., TwoSlopeNorm).
        mask : np.ndarray, optional
            Boolean mask for selecting points. If None, use is_land.
        s : float
            Marker size.
        alpha : float
            Marker transparency.

        Returns
        -------
        PathCollection
            Scatter plot collection.
        """
        if mask is None:
            mask = self.is_land

        scatter_kwargs = {
            "c": values[mask],
            "cmap": cmap,
            "s": s,
            "alpha": alpha,
        }

        if norm is not None:
            scatter_kwargs["norm"] = norm
        else:
            scatter_kwargs["vmin"] = vmin
            scatter_kwargs["vmax"] = vmax

        sc = ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            **scatter_kwargs,
        )

        return sc

    def plot_site_scatter(
        self,
        ax: Axes,
        coords: np.ndarray,
        values: np.ndarray,
        cmap: str = CMAP_PROBABILITY,
        vmin: float = PROB_VMIN,
        vmax: float = PROB_VMAX,
        norm: Optional[Normalize] = None,
        s: float = SITE_SCATTER_SIZE,
        alpha: float = SITE_SCATTER_ALPHA,
        edgecolors: str = "black",
        linewidths: float = SITE_EDGE_WIDTH,
    ):
        """Plot site scatter points with edge colors.

        Parameters
        ----------
        ax : Axes
            Matplotlib axes object.
        coords : np.ndarray
            Coordinates (n_sites, 2).
        values : np.ndarray
            Values for coloring (n_sites,).
        cmap : str
            Colormap name.
        vmin : float
            Minimum value for colormap.
        vmax : float
            Maximum value for colormap.
        norm : Normalize, optional
            Normalization object.
        s : float
            Marker size.
        alpha : float
            Marker transparency.
        edgecolors : str
            Edge color for markers.
        linewidths : float
            Edge width for markers.

        Returns
        -------
        PathCollection
            Scatter plot collection.
        """
        scatter_kwargs = {
            "c": values,
            "cmap": cmap,
            "s": s,
            "alpha": alpha,
            "edgecolors": edgecolors,
            "linewidths": linewidths,
        }

        if norm is not None:
            scatter_kwargs["norm"] = norm
        else:
            scatter_kwargs["vmin"] = vmin
            scatter_kwargs["vmax"] = vmax

        sc = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            **scatter_kwargs,
        )

        return sc

    def plot_boundary(
        self,
        ax: Axes,
        color: str = BOUNDARY_COLOR,
        s: float = BOUNDARY_SIZE,
    ) -> None:
        """Plot boundary points.

        Parameters
        ----------
        ax : Axes
            Matplotlib axes object.
        color : str
            Color for boundary points.
        s : float
            Marker size.
        """
        ax.scatter(
            self.boundary[:, 0],
            self.boundary[:, 1],
            c=color,
            s=s,
        )

    def add_colorbar(
        self,
        fig: Figure,
        sc,
        ax: Axes,
        label: str = "",
    ) -> Colorbar:
        """Add colorbar to the figure.

        Parameters
        ----------
        fig : Figure
            Matplotlib figure object.
        sc : PathCollection
            Scatter plot collection.
        ax : Axes
            Matplotlib axes object.
        label : str
            Label for the colorbar.

        Returns
        -------
        Colorbar
            Matplotlib colorbar object.
        """
        cbar = plt.colorbar(sc, ax=ax, label=label)
        return cbar

    def set_axis_style(
        self,
        ax: Axes,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        spine_width: float = SPINE_WIDTH,
    ) -> None:
        """Set axis styling.

        Parameters
        ----------
        ax : Axes
            Matplotlib axes object.
        title : str, optional
            Title for the axis.
        xlabel : str, optional
            X-axis label.
        ylabel : str, optional
            Y-axis label.
        spine_width : float
            Width of the axis spines.
        """
        if title is not None:
            ax.set_title(title)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)

        set_spine_width(ax, spine_width)

    def save_figure(
        self,
        fig: Figure,
        filename: Union[str, Path],
        dpi: int = FIGURE_DPI,
        format: str = FIGURE_FORMAT,
    ) -> None:
        """Save figure to file.

        Parameters
        ----------
        fig : Figure
            Matplotlib figure object.
        filename : str or Path
            Output filename.
        dpi : int
            Resolution in dots per inch.
        format : str
            File format (e.g., "png", "pdf").
        """
        fig.savefig(filename, dpi=dpi, format=format, bbox_inches="tight")
