"""Dataset class for the marked point process model."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, List, Optional, Sequence

import numpy as np

if TYPE_CHECKING:
    from bayesian_statistics.models.preprocessing.data_preprocessor import (
        ObsidianDataPreprocessor,
    )


@dataclass
class MarkedPointProcessDataset:
    """Dataset for the marked point process MCMC sampler.

    This class holds all data required for the joint model:
    - Site locations (point process X)
    - Mark data (category counts y)
    - Design matrices for intensity and composition
    - Optional grid data for prediction

    Attributes
    ----------
    site_coords : np.ndarray
        Site coordinates, shape (n_sites, 2).
    counts : np.ndarray
        Category counts per site, shape (n_sites, K).
    origins : List[str]
        Category names (K total, including baseline).
    site_ids : Optional[np.ndarray]
        Site identifiers, shape (n_sites,).
    total_counts : Optional[np.ndarray]
        Total counts per site, shape (n_sites,). Computed if not provided.
    design_matrix_intensity : Optional[np.ndarray]
        Design matrix for intensity, shape (n_sites, p_int).
    design_matrix_marks : Optional[np.ndarray]
        Design matrix for marks, shape (n_sites, p_z).
    grid_coords : Optional[np.ndarray]
        Grid coordinates for prediction, shape (n_grid, 2).
    design_matrix_grid_intensity : Optional[np.ndarray]
        Design matrix for grid (intensity), shape (n_grid, p_int).
    design_matrix_grid_marks : Optional[np.ndarray]
        Design matrix for grid (marks), shape (n_grid, p_z).
    valid_grids : Optional[np.ndarray]
        Boolean mask for valid grid points, shape (n_grid,).
    region : Optional[List[List[float]]]
        Spatial region [[x_min, x_max], [y_min, y_max]].
    period : Optional[int]
        Time period identifier.
    distance_features_sites : Optional[np.ndarray]
        Distance features for sites, shape (n_sites, K-1).
    distance_features_grid : Optional[np.ndarray]
        Distance features for grid, shape (n_grid, K-1).
    prior_mean_intercept_sites : Optional[np.ndarray]
        Prior mean for site intercepts, shape (K-1, n_sites).
    prior_mean_intercept_grid : Optional[np.ndarray]
        Prior mean for grid intercepts, shape (K-1, n_grid).
    """

    site_coords: np.ndarray
    counts: np.ndarray
    origins: List[str]

    site_ids: Optional[np.ndarray] = None
    total_counts: Optional[np.ndarray] = None
    design_matrix_intensity: Optional[np.ndarray] = None
    design_matrix_marks: Optional[np.ndarray] = None
    grid_coords: Optional[np.ndarray] = None
    design_matrix_grid_intensity: Optional[np.ndarray] = None
    design_matrix_grid_marks: Optional[np.ndarray] = None
    valid_grids: Optional[np.ndarray] = None
    region: Optional[List[List[float]]] = None
    period: Optional[int] = None
    distance_features_sites: Optional[np.ndarray] = None
    distance_features_grid: Optional[np.ndarray] = None
    prior_mean_intercept_sites: Optional[np.ndarray] = None
    prior_mean_intercept_grid: Optional[np.ndarray] = None

    volume: float = field(init=False, default=0.0)

    def __post_init__(self):
        """Validate inputs and compute derived quantities."""
        self._validate()
        self._compute_defaults()

    def _validate(self):
        """Validate input data consistency."""
        n_sites = self.site_coords.shape[0]
        n_counts_rows = self.counts.shape[0]
        n_categories = self.counts.shape[1]
        n_origins = len(self.origins)

        if n_counts_rows != n_sites:
            raise ValueError(
                f"counts rows ({n_counts_rows}) must match "
                f"site_coords rows ({n_sites})"
            )

        if n_origins != n_categories:
            raise ValueError(
                f"origins length ({n_origins}) must match "
                f"counts columns ({n_categories})"
            )

        if self.design_matrix_intensity is not None:
            if self.design_matrix_intensity.shape[0] != n_sites:
                raise ValueError(
                    f"design_matrix_intensity rows "
                    f"({self.design_matrix_intensity.shape[0]}) must match "
                    f"num_sites ({n_sites})"
                )

        if self.design_matrix_marks is not None:
            if self.design_matrix_marks.shape[0] != n_sites:
                raise ValueError(
                    f"design_matrix_marks rows "
                    f"({self.design_matrix_marks.shape[0]}) must match "
                    f"num_sites ({n_sites})"
                )

    def _compute_defaults(self):
        """Compute default values for optional fields."""
        n_sites = self.site_coords.shape[0]

        # Compute total_counts if not provided
        if self.total_counts is None:
            self.total_counts = self.counts.sum(axis=1)

        # Default design matrices: intercept only
        if self.design_matrix_intensity is None:
            self.design_matrix_intensity = np.ones((n_sites, 1))

        if self.design_matrix_marks is None:
            self.design_matrix_marks = np.ones((n_sites, 1))

        # Infer region from coordinates if not provided
        if self.region is None:
            x_min, x_max = self.site_coords[:, 0].min(), self.site_coords[:, 0].max()
            y_min, y_max = self.site_coords[:, 1].min(), self.site_coords[:, 1].max()
            # Add small buffer
            buffer_x = 0.05 * (x_max - x_min) if x_max > x_min else 0.01
            buffer_y = 0.05 * (y_max - y_min) if y_max > y_min else 0.01
            self.region = [
                [x_min - buffer_x, x_max + buffer_x],
                [y_min - buffer_y, y_max + buffer_y],
            ]

        # Compute volume (area) from region
        x_range = self.region[0][1] - self.region[0][0]
        y_range = self.region[1][1] - self.region[1][0]
        self.volume = x_range * y_range

    def num_sites(self) -> int:
        """Return the number of sites (n_X)."""
        return self.site_coords.shape[0]

    def num_categories(self) -> int:
        """Return the number of categories (K, including baseline)."""
        return self.counts.shape[1]

    def num_grid(self) -> int:
        """Return the number of grid points."""
        if self.grid_coords is None:
            return 0
        return self.grid_coords.shape[0]


def prepare_marked_point_process_dataset(
    preprocessor: "ObsidianDataPreprocessor",
    period: int,
    origins: Sequence[str],
    grid_subsample_ratio: float = 0.01,
    drop_zero_total_sites: bool = True,
    intensity_variable_names: Optional[Sequence[str]] = None,
) -> MarkedPointProcessDataset:
    """Prepare dataset for the marked point process model.

    This function extracts count data and coordinates from the preprocessor
    and constructs a MarkedPointProcessDataset ready for MCMC sampling.

    Parameters
    ----------
    preprocessor
        ObsidianDataPreprocessor with loaded data.
    period
        Time period to extract (0-4).
    origins
        Category names (K total, including baseline as last element).
    grid_subsample_ratio
        Fraction of grid points to use (0 < ratio <= 1).
    drop_zero_total_sites
        Whether to drop sites with zero total count.
    intensity_variable_names
        Variable names to use as covariates for intensity model.
        If None, only intercept is used.
        Available variables include: "average_elevation", "average_slope_angle",
        "cost_river", etc. (see preprocessor schema).

    Returns
    -------
    MarkedPointProcessDataset
        Dataset ready for MCMC sampling.

    Examples
    --------
    >>> preprocessor = ObsidianDataPreprocessor(data_dir)
    >>> preprocessor.load_data()
    >>> origins = ["神津島", "信州", "箱根", "高原山", "その他"]
    >>> # With covariates
    >>> dataset = prepare_marked_point_process_dataset(
    ...     preprocessor, period=2, origins=origins,
    ...     intensity_variable_names=["average_elevation", "average_slope_angle"],
    ... )
    """
    import polars as pl

    from bayesian_statistics.nngp.model.multinomial_model import _prepare_counts

    # Get counts using the shared utility function
    counts, totals = _prepare_counts(preprocessor, period, origins)

    # Get site coordinates
    site_df = preprocessor.df_sites.sort("遺跡ID")
    coords = site_df.select(["経度", "緯度"]).to_numpy().astype(float)
    site_ids = site_df.select("遺跡ID").to_numpy().astype(int).ravel()

    if counts.shape[0] != site_ids.shape[0]:
        raise ValueError(
            "Mismatch between number of sites in counts and metadata; "
            "ensure site identifiers are contiguous and start at zero."
        )

    # Get intensity covariates from preprocessor
    design_matrix_intensity_sites = None
    design_matrix_intensity_grids_full = None

    if intensity_variable_names is not None and len(intensity_variable_names) > 0:
        W_grids, W_sites = preprocessor.create_explanatory_variables(
            list(intensity_variable_names)
        )
        # Add intercept
        n_sites_full = W_sites.shape[0]
        design_matrix_intensity_sites = np.column_stack([
            np.ones(n_sites_full),
            W_sites,
        ])
        n_grids_full = W_grids.shape[0]
        design_matrix_intensity_grids_full = np.column_stack([
            np.ones(n_grids_full),
            W_grids,
        ])

    # Filter sites with zero total count
    if drop_zero_total_sites:
        mask = totals > 0
        if not np.any(mask):
            raise ValueError("drop_zero_total_sites left no observations")
        coords = coords[mask]
        counts = counts[mask]
        totals = totals[mask]
        site_ids = site_ids[mask]
        if design_matrix_intensity_sites is not None:
            design_matrix_intensity_sites = design_matrix_intensity_sites[mask]

    # Prepare grid data
    elevation_df = preprocessor.df_elevation.sort(["y", "x"])
    n_grid_full = elevation_df.shape[0]

    grid_coords_full = elevation_df.select(["x", "y"]).to_numpy().astype(float)
    is_sea_full = (
        elevation_df.select(pl.col("is_sea").cast(pl.Boolean)).to_numpy().flatten()
    )
    is_valid_full = elevation_df.select("is_valid").to_numpy().flatten().astype(bool)

    # Subsample grid
    if not (0 < grid_subsample_ratio <= 1.0):
        raise ValueError("grid_subsample_ratio must be in (0, 1]")

    if grid_subsample_ratio < 1.0:
        keep = max(int(np.floor(n_grid_full * grid_subsample_ratio)), 1)
        indices = np.linspace(0, n_grid_full - 1, keep, dtype=int)
        grid_coords = grid_coords_full[indices]
        valid_grids = (~is_sea_full[indices]) & is_valid_full[indices]
        if design_matrix_intensity_grids_full is not None:
            design_matrix_grid_intensity = design_matrix_intensity_grids_full[indices]
        else:
            design_matrix_grid_intensity = None
    else:
        grid_coords = grid_coords_full
        valid_grids = (~is_sea_full) & is_valid_full
        design_matrix_grid_intensity = design_matrix_intensity_grids_full

    # Compute region from grid
    region = [
        [float(grid_coords[:, 0].min()), float(grid_coords[:, 0].max())],
        [float(grid_coords[:, 1].min()), float(grid_coords[:, 1].max())],
    ]

    return MarkedPointProcessDataset(
        site_coords=coords,
        counts=counts,
        origins=list(origins),
        site_ids=site_ids,
        total_counts=totals,
        design_matrix_intensity=design_matrix_intensity_sites,
        grid_coords=grid_coords,
        design_matrix_grid_intensity=design_matrix_grid_intensity,
        valid_grids=valid_grids,
        region=region,
        period=period,
    )
