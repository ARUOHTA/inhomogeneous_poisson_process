from typing import List, Optional, Tuple

import numpy as np
import polars as pl
import polars_st as pst
from shapely import contains_xy
from shapely.geometry import Polygon


def prepare_grid_samples(
    polygon: Polygon,
    covariates_df: Optional[pl.DataFrame] = None,
    covariates_names: Optional[List[str]] = None,
    step: float = 0.001,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    ポリゴン内部の一様グリッド座標を生成し、必要なら共変量を空間結合で付与する。

    Returns
    -------
    points : (n_points, 2) ndarray
        [x, y] 座標
    covariates : (n_points, n_covariates) ndarray or None
        指定した共変量列を NumPy 配列で返す（共変量を付けない場合は None）
    """
    minx, miny, maxx, maxy = polygon.bounds

    # グリッド座標を生成
    x_coords = np.arange(minx, maxx, step)
    y_coords = np.arange(miny, maxy, step)

    # メッシュグリッドを作成
    x_grid, y_grid = np.meshgrid(x_coords, y_coords, indexing="xy")

    # ベクトル化 contains_xy でポリゴン内の点を抽出
    mask = contains_xy(polygon, x_grid.ravel(), y_grid.ravel())
    valid_x = x_grid.ravel()[mask]
    valid_y = y_grid.ravel()[mask]
    valid_points = np.column_stack([valid_x, valid_y])

    if covariates_df is None or not covariates_names:
        return valid_points, None

    # 共変量の結合（大量データ時はチャンク処理）
    chunk_size = 100_000
    n_points = len(valid_points)

    def sjoin_cov(points: np.ndarray) -> np.ndarray:
        df = pl.DataFrame({"X": points[:, 0], "Y": points[:, 1]}).with_columns(
            geometry=pst.point(pl.concat_arr(["X", "Y"])).st.set_srid(4326)
        )
        # 左: 点, 右: ポリゴン → 述語は contains（右が左を含む）
        joined = df.st.sjoin(covariates_df, predicate="contains", how="left")
        return joined.select(covariates_names).to_numpy()

    if n_points > chunk_size:
        cov_list = []
        for i in range(0, n_points, chunk_size):
            cov_list.append(sjoin_cov(valid_points[i : i + chunk_size]))
        covariate_values = np.vstack(cov_list) if cov_list else None
    else:
        covariate_values = sjoin_cov(valid_points)

    return valid_points, covariate_values


class PrecomputedGridSampler:
    """事前計算グリッドからの高速サンプリングユーティリティ"""

    def __init__(
        self,
        polygon: Polygon,
        covariates_df: Optional[pl.DataFrame] = None,
        covariates_names: Optional[List[str]] = None,
        step: float = 0.001,
    ) -> None:
        self.points, self.covariates = prepare_grid_samples(
            polygon, covariates_df, covariates_names, step
        )
        self.n_points = len(self.points)
        self.has_covariates = self.covariates is not None
        self.n_covariates = self.covariates.shape[1] if self.has_covariates else 0
        self.indices = np.arange(self.n_points)

    def sample_point(self) -> Tuple[float, float]:
        idx = np.random.randint(0, self.n_points)
        return float(self.points[idx, 0]), float(self.points[idx, 1])

    def sample_point_with_covariates(self) -> Tuple[np.ndarray, np.ndarray]:
        idx = np.random.randint(0, self.n_points)
        point = self.points[idx]
        covariates = self.covariates[idx] if self.has_covariates else np.array([])
        return point, covariates

    def sample_multiple_points(self, n: int, replace: bool = True) -> np.ndarray:
        if replace:
            indices = np.random.randint(0, self.n_points, size=n)
        else:
            indices = np.random.choice(
                self.indices, size=min(n, self.n_points), replace=False
            )
        return self.points[indices]

    def sample_multiple_with_covariates(
        self, n: int, replace: bool = True
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        if replace:
            indices = np.random.randint(0, self.n_points, size=n)
        else:
            indices = np.random.choice(
                self.indices, size=min(n, self.n_points), replace=False
            )
        points = self.points[indices]
        covariates = self.covariates[indices] if self.has_covariates else None
        return points, covariates

    def sample_multiple_with_covariates_idx(
        self, n: int, replace: bool = True
    ) -> Tuple[np.ndarray, Optional[np.ndarray], np.ndarray]:
        """Like sample_multiple_with_covariates but also returns indices."""
        if replace:
            indices = np.random.randint(0, self.n_points, size=n)
        else:
            indices = np.random.choice(
                self.indices, size=min(n, self.n_points), replace=False
            )
        points = self.points[indices]
        covariates = self.covariates[indices] if self.has_covariates else None
        return points, covariates, indices


def create_fast_sampler(
    polygon: Polygon,
    covariates_df: Optional[pl.DataFrame] = None,
    covariates_names: Optional[List[str]] = None,
    step: float = 0.001,
) -> PrecomputedGridSampler:
    """高速サンプラーを作成"""
    return PrecomputedGridSampler(polygon, covariates_df, covariates_names, step)


__all__ = [
    "prepare_grid_samples",
    "PrecomputedGridSampler",
    "create_fast_sampler",
]
