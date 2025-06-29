"""
Nadaraya-Watson推定量クラス
黒曜石の産地構成比を推定する
"""

from typing import Dict, Optional, Tuple

import numpy as np
import polars as pl
from tqdm import tqdm

from .model3_preprocessing import ObsidianDataPreprocessor


class NadarayaWatsonEstimator:
    """Nadaraya-Watson推定量による産地構成比推定"""

    def __init__(self, sigma: float = 500, sigma_for_sites: float = 0.1):
        """
        Parameters
        ----------
        sigma : float
            グリッド間のカーネルバンド幅
        sigma_for_sites : float
            遺跡間のカーネルバンド幅
        """
        self.sigma = sigma
        self.sigma_for_sites = sigma_for_sites
        self._weights: Optional[np.ndarray] = None
        self._weights_sites: Optional[np.ndarray] = None
        self._results: Dict[str, np.ndarray] = {}

    @staticmethod
    def K(x: np.ndarray, sigma: float) -> np.ndarray:
        """
        ガウスカーネル関数

        Parameters
        ----------
        x : np.ndarray
            距離
        sigma : float
            バンド幅

        Returns
        -------
        np.ndarray
            カーネル重み
        """
        return np.exp(-0.5 * (x**2) / (sigma**2)) / (2 * np.pi * sigma**2)

    def calculate_kernel_weights(
        self, distances: np.ndarray, sigma: Optional[float] = None
    ) -> np.ndarray:
        """
        カーネル重みを計算

        Parameters
        ----------
        distances : np.ndarray
            距離行列
        sigma : float, optional
            バンド幅（Noneの場合はself.sigmaを使用）

        Returns
        -------
        np.ndarray
            カーネル重み
        """
        if sigma is None:
            sigma = self.sigma
        return self.K(distances, sigma)

    @staticmethod
    def calculate_distance_W(
        W_grid: np.ndarray,
        W_sites: np.ndarray,
    ) -> np.ndarray:
        """
        説明変数間の距離行列を計算

        Parameters
        ----------
        W_grid : np.ndarray
            グリッドの説明変数 (グリッド数, p)
        W_sites : np.ndarray
            遺跡の説明変数 (遺跡数, p)

        Returns
        -------
        np.ndarray
            説明変数間の距離行列 (グリッド数, 遺跡数, p)
        """
        # W_gridを(グリッド数, 1, p)に、W_sitesを(1, 遺跡数, p)に reshape
        W_grid_expanded = W_grid[:, np.newaxis, :]
        W_sites_expanded = W_sites[np.newaxis, :, :]

        # ユークリッド距離の計算
        distances = np.abs(W_grid_expanded - W_sites_expanded)

        return distances

    def calculate_weighted_ratios(
        self, weights: np.ndarray, counts: np.ndarray, target_counts: np.ndarray
    ) -> np.ndarray:
        """
        重み付き比率を計算

        Parameters
        ----------
        weights : np.ndarray
            重み行列
        counts : np.ndarray
            各遺跡での出土数
        target_counts : np.ndarray
            対象産地の出土数

        Returns
        -------
        ratios : np.ndarray
            各グリッド点での重み付き比率
        """
        # 重み付き合計を計算
        weighted_total = np.sum(weights * counts, axis=1)
        weighted_target = np.sum(weights * target_counts, axis=1)

        # 比率計算（0除算を防ぐ）
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = np.where(weighted_total > 0, weighted_target / weighted_total, 0)

        return ratios

    @staticmethod
    def _create_land_mask_original(grid_coords, df_elevation, lon_mesh, lat_mesh):
        """元のnotebook実装のcreate_land_mask関数"""
        # 地形マスクの作成
        land_points = df_elevation.select(
            ["x", "y", pl.col("is_sea").cast(pl.Boolean)]
        ).to_numpy()

        lons_1d = lon_mesh[0, :]
        lats_1d = lat_mesh[:, 0]
        land_mask = np.full(lon_mesh.shape, False)

        x_indices = np.searchsorted(lons_1d, land_points[:, 0])
        y_indices = np.searchsorted(lats_1d, land_points[:, 1])
        valid_points = (
            (x_indices > 0)
            & (x_indices < len(lons_1d))
            & (y_indices > 0)
            & (y_indices < len(lats_1d))
        )
        is_sea = land_points[valid_points, 2].astype(bool)
        land_mask[y_indices[valid_points], x_indices[valid_points]] = ~is_sea

        # grid_coordsの各点について、対応するland_maskの値を取得
        grid_lons = grid_coords[:, 1] * 180 / np.pi  # ラジアンから度に変換
        grid_lats = grid_coords[:, 0] * 180 / np.pi

        grid_x_indices = np.searchsorted(lons_1d, grid_lons)
        grid_y_indices = np.searchsorted(lats_1d, grid_lats)

        # インデックスが有効範囲内にあることを確認
        valid_grid_points = (
            (grid_x_indices > 0)
            & (grid_x_indices < len(lons_1d))
            & (grid_y_indices > 0)
            & (grid_y_indices < len(lats_1d))
        )

        # 海上の点の重みを0に設定
        grid_is_land = np.zeros(len(grid_coords), dtype=bool)
        grid_is_land[valid_grid_points] = land_mask[
            grid_y_indices[valid_grid_points], grid_x_indices[valid_grid_points]
        ]

        return grid_is_land

    def fit(
        self, preprocessor: ObsidianDataPreprocessor, variable_names: list
    ) -> Dict[str, np.ndarray]:
        """
        モデルを学習

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ
        variable_names : list
            使用する説明変数名

        Returns
        -------
        Dict[str, np.ndarray]
            計算された重み行列
        """
        print("creating weights matrix...")

        # 説明変数の作成
        W_grids, W_sites = preprocessor.create_explanatory_variables(variable_names)

        # 距離行列の読み込み
        distances = preprocessor.load_tobler_distances()

        # 重み行列の計算
        self._weights = self.calculate_kernel_weights(distances, self.sigma)

        # 説明変数間の距離行列の計算
        print("calculating distance_W...")
        distance_W = self.calculate_distance_W(W_grids, W_sites)
        weights_W = self.K(distance_W, self.sigma).prod(axis=2)

        self._weights *= weights_W

        print("updating weights matrix...")

        # 陸地マスクの適用
        grid_coords = preprocessor.create_grid_coords()
        lon_mesh, lat_mesh = preprocessor.create_meshgrid()
        grid_is_land = self._create_land_mask_original(
            grid_coords, preprocessor.df_elevation, lon_mesh, lat_mesh
        )

        # 海上の点からの重みをすべて0に
        self._weights *= grid_is_land[:, np.newaxis]

        # 遺跡間の重み計算
        site_grid_coords = preprocessor.convert_to_grid_coords()
        grid_info = preprocessor.create_grid_info()
        width = grid_info["n_grid_x"]

        def idx(x, y):
            return y * width + x

        site_grid_idx = np.vectorize(idx)(
            site_grid_coords[:, 0], site_grid_coords[:, 1]
        )
        distances_sites = distances[site_grid_idx]
        self._weights_sites = self.calculate_kernel_weights(
            distances_sites, self.sigma_for_sites
        )

        return {"weights": self._weights, "weights_sites": self._weights_sites}

    def predict_single(
        self,
        preprocessor: ObsidianDataPreprocessor,
        target_period: int,
        target_origin: str,
    ) -> Dict[str, np.ndarray]:
        """
        単一の時期・産地について予測

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ
        target_period : int
            対象時期
        target_origin : str
            対象産地

        Returns
        -------
        Dict[str, np.ndarray]
            グリッド上の比率と遺跡での比率
        """
        if self._weights is None or self._weights_sites is None:
            raise ValueError(
                "モデルが学習されていません。fit()を先に実行してください。"
            )

        # データの前処理
        counts, target_counts = preprocessor.preprocess_obsidian_data(
            target_period, target_origin
        )

        # メッシュグリッドの形状を取得
        lon_mesh, lat_mesh = preprocessor.create_meshgrid()

        # グリッド上の比率計算
        ratio_mesh = self.calculate_weighted_ratios(
            self._weights, counts, target_counts
        ).reshape(lon_mesh.shape)

        # 遺跡での比率計算
        ratio_sites = self.calculate_weighted_ratios(
            self._weights_sites, counts, target_counts
        )

        return {
            "ratio_mesh": ratio_mesh,
            "ratio_sites": ratio_sites,
        }

    def predict_all_periods_origins(
        self,
        preprocessor: ObsidianDataPreprocessor,
        time_period_names: Dict[int, str],
        origin_order: list,
    ) -> Tuple[pl.DataFrame, pl.DataFrame]:
        """
        全時期・全産地について予測

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ
        time_period_names : Dict[int, str]
            時期名の辞書
        origin_order : list
            産地のリスト（最後の要素"その他"は除外）

        Returns
        -------
        ratio_df : pl.DataFrame
            グリッド上の比率データフレーム
        ratio_sites_df : pl.DataFrame
            遺跡での比率データフレーム
        """
        if self._weights is None or self._weights_sites is None:
            raise ValueError(
                "モデルが学習されていません。fit()を先に実行してください。"
            )

        # メッシュグリッドの作成
        lon_mesh, lat_mesh = preprocessor.create_meshgrid()

        # 初期データフレーム
        ratio_df = pl.DataFrame({"x": lon_mesh.ravel(), "y": lat_mesh.ravel()})

        # 遺跡の一意なIDを取得
        unique_sites = preprocessor.df_obsidian.unique(subset=["遺跡ID"]).sort("遺跡ID")
        ratio_sites_df = pl.DataFrame({"遺跡ID": unique_sites["遺跡ID"]})

        # 全時期・全産地について計算
        for target_period in tqdm(time_period_names.keys(), desc="時期"):
            for target_origin in origin_order[:-1]:  # "その他"を除外
                print(f"target_period: {target_period}, target_origin: {target_origin}")

                # 予測実行
                results = self.predict_single(
                    preprocessor, target_period, target_origin
                )

                # グリッドデータの追加
                ratio_df = ratio_df.join(
                    pl.DataFrame(
                        {
                            "x": lon_mesh.ravel(),
                            "y": lat_mesh.ravel(),
                            f"ratio_{target_period}_{target_origin}": results[
                                "ratio_mesh"
                            ].ravel(),
                        }
                    ),
                    on=["x", "y"],
                )

                # 遺跡データの追加
                ratio_sites_df = ratio_sites_df.join(
                    pl.DataFrame(
                        {
                            "遺跡ID": unique_sites["遺跡ID"],
                            f"比率_{target_period}_{target_origin}": results[
                                "ratio_sites"
                            ],
                        }
                    ),
                    on="遺跡ID",
                )

        return ratio_df, ratio_sites_df

    @property
    def weights(self) -> np.ndarray:
        """重み行列"""
        if self._weights is None:
            raise ValueError(
                "モデルが学習されていません。fit()を先に実行してください。"
            )
        return self._weights

    @property
    def weights_sites(self) -> np.ndarray:
        """遺跡間の重み行列"""
        if self._weights_sites is None:
            raise ValueError(
                "モデルが学習されていません。fit()を先に実行してください。"
            )
        return self._weights_sites
