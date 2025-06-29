"""
Model 3rd の前処理クラス
データの読み込み、グリッド座標変換、Tobler距離の管理などを行う
"""

import os
import pickle
from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl


class ObsidianDataPreprocessor:
    """黒曜石データの前処理を行うクラス"""

    # データスキーマ定義
    SCHEMA = pl.Schema(
        {
            "grid_x": pl.Int64,
            "x": pl.Float64,
            "grid_y": pl.Int64,
            "y": pl.Float64,
            "mesh_code_5th": pl.Int64,
            "average_elevation": pl.Float64,
            "maximum_elevation": pl.Float64,
            "minimum_elevation": pl.Float64,
            "minimum_elevation_code": pl.Int64,
            "maximum_slope_angle": pl.Float64,
            "maximum_slope_direction": pl.Int64,
            "minimum_slope_angle": pl.Float64,
            "minimum_slope_direction": pl.Int64,
            "average_slope_angle": pl.Float64,
            "geometry": pl.Utf8,
            "is_sea": pl.Boolean,
            "walking_velocity": pl.Float64,
            "travel_time": pl.Float64,
            "elevation_diff_east": pl.Float64,
            "angle_east": pl.Float64,
            "walking_velocity_east": pl.Float64,
            "travel_time_east": pl.Float64,
            "elevation_diff_west": pl.Float64,
            "angle_west": pl.Float64,
            "walking_velocity_west": pl.Float64,
            "travel_time_west": pl.Float64,
            "elevation_diff_north": pl.Float64,
            "angle_north": pl.Float64,
            "walking_velocity_north": pl.Float64,
            "travel_time_north": pl.Float64,
            "elevation_diff_south": pl.Float64,
            "angle_south": pl.Float64,
            "walking_velocity_south": pl.Float64,
            "travel_time_south": pl.Float64,
            "cost_kouzu": pl.Float64,
            "cost_shinshu": pl.Float64,
            "cost_hakone": pl.Float64,
            "cost_takahara": pl.Float64,
            "cost_river": pl.Float64,
            "x_meter": pl.Float64,
            "y_meter": pl.Float64,
        }
    )

    def __init__(self, data_dir: str):
        """
        Parameters
        ----------
        data_dir : str
            データディレクトリのパス
        """
        self.data_dir = data_dir
        self._df_elevation: Optional[pl.DataFrame] = None
        self._df_obsidian: Optional[pl.DataFrame] = None
        self._df_sites: Optional[pl.DataFrame] = None
        self._grid_info: Optional[Dict[str, float]] = None

    def load_data(self) -> Dict[str, pl.DataFrame]:
        """
        データを読み込む

        Returns
        -------
        Dict[str, pl.DataFrame]
            読み込んだデータフレームの辞書
        """
        # 標高データの読み込み
        self._df_elevation = pl.read_csv(
            os.path.join(self.data_dir, "11_gdf_elevation.csv"),
            null_values=["nan"],  # "nan"をnullとして扱う
        )
        self._df_elevation = (
            self._df_elevation.cast(self.SCHEMA)
            .with_columns(
                [
                    pl.all_horizontal(
                        [~pl.col(col).is_null() for col in self._df_elevation.columns]
                    ).alias("is_valid")
                ]
            )
        )

        # 黒曜石データの読み込み
        self._df_obsidian = pl.read_csv(
            os.path.join(self.data_dir, "11_gdf_obsidian.csv")
        )

        # 遺跡データの読み込み
        self._df_sites = pl.read_csv(os.path.join(self.data_dir, "11_gdf_sites.csv"))

        # 遺跡データに標高データをジョイン
        self._df_sites = self._df_sites.join(
            self._df_elevation.drop(["x", "y"]).with_columns(
                [pl.col("mesh_code_5th").cast(pl.Int64).alias("mesh_code_5th")]
            ),
            left_on="メッシュコード",
            right_on="mesh_code_5th",
            how="left",
        )

        return {
            "elevation": self._df_elevation,
            "obsidian": self._df_obsidian,
            "sites": self._df_sites,
        }

    def convert_to_grid_coords(self, df_obsidian: Optional[pl.DataFrame] = None) -> np.ndarray:
        """
        遺跡の座標をグリッド座標に変換

        Parameters
        ----------
        df_obsidian : pl.DataFrame, optional
            黒曜石データフレーム（Noneの場合は内部データを使用）

        Returns
        -------
        np.ndarray
            遺跡のグリッド座標
        """
        if df_obsidian is None:
            df_obsidian = self._df_obsidian
        
        if self._df_elevation is None:
            raise ValueError("データが読み込まれていません。load_data()を先に実行してください。")

        x_spacing = self._df_elevation["x"].unique().sort().diff().mode()[0]
        y_spacing = self._df_elevation["y"].unique().sort().diff().mode()[0]
        x_first_center = self._df_elevation["x"].min()
        y_first_center = self._df_elevation["y"].min()

        coords = (
            df_obsidian.select([(pl.col("遺跡ID")), (pl.col("緯度")), (pl.col("経度"))])
            .unique(subset=["遺跡ID"])
            .sort("遺跡ID")
            .with_columns(
                [
                    ((pl.col("経度") - (x_first_center - (x_spacing / 2))) / x_spacing)
                    .cast(pl.Int64)
                    .alias("grid_x"),
                    ((pl.col("緯度") - (y_first_center - (y_spacing / 2))) / y_spacing)
                    .cast(pl.Int64)
                    .alias("grid_y"),
                ]
            )
        )

        site_grid_coords = np.column_stack(
            [coords["grid_x"].to_numpy(), coords["grid_y"].to_numpy()]
        )

        return site_grid_coords

    def load_tobler_distances(self, distance_dir: str = "16_tobler_distance_with_coast_50_average") -> np.ndarray:
        """
        Tobler距離を読み込む

        Parameters
        ----------
        distance_dir : str
            距離データのディレクトリ名

        Returns
        -------
        np.ndarray
            距離行列 (グリッド数, 遺跡数)
        """
        if self._df_obsidian is None or self._df_elevation is None:
            raise ValueError("データが読み込まれていません。load_data()を先に実行してください。")

        # グリッド座標を取得
        site_coords = self.create_site_coords()
        grid_coords = self.create_grid_coords()

        # 距離行列の読み込み
        distances = np.zeros((len(grid_coords), len(site_coords)))
        for i in range(len(site_coords)):
            with open(
                os.path.join(self.data_dir, distance_dir, f"distance_siteID_{i}"),
                mode="br",
            ) as fi:
                min_costs = pickle.load(fi)
                distances[:, i] = min_costs

        return distances

    def create_explanatory_variables(
        self, variable_names: List[str]
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        説明変数を作成

        Parameters
        ----------
        variable_names : List[str]
            使用する変数名のリスト

        Returns
        -------
        W_grids : np.ndarray
            グリッドの説明変数 (グリッド数, p)
        W_sites : np.ndarray
            遺跡の説明変数 (遺跡数, p)
        """
        if self._df_elevation is None or self._df_sites is None:
            raise ValueError("データが読み込まれていません。load_data()を先に実行してください。")

        # グリッドの説明変数
        W_grids = (
            self._df_elevation.sort(["y", "x"])
            .select(variable_names)
            .to_numpy()
            .astype(np.float64)
        )

        # 遺跡の説明変数
        W_sites = (
            self._df_sites.sort("遺跡ID")
            .select(variable_names)
            .to_numpy()
            .astype(np.float64)
        )

        return W_grids, W_sites

    def create_grid_info(self) -> Dict[str, float]:
        """
        グリッド情報を作成

        Returns
        -------
        Dict[str, float]
            グリッド情報の辞書
        """
        if self._df_elevation is None:
            raise ValueError("データが読み込まれていません。load_data()を先に実行してください。")

        if self._grid_info is not None:
            return self._grid_info

        x_unique = self._df_elevation.select("x").unique().to_series().to_list()
        y_unique = self._df_elevation.select("y").unique().to_series().to_list()
        x_unique.sort()
        y_unique.sort()

        if len(x_unique) < 2 or len(y_unique) < 2:
            raise ValueError("x あるいは y のユニーク値が2つ以上ありません。")

        # 隣り合う中心座標の差からステップ幅を算出
        delta_x = x_unique[1] - x_unique[0]
        delta_y = y_unique[1] - y_unique[0]

        self._grid_info = {
            "x_min_center": x_unique[0],
            "x_max_center": x_unique[-1],
            "y_min_center": y_unique[0],
            "y_max_center": y_unique[-1],
            "delta_x": delta_x,
            "delta_y": delta_y,
            "n_grid_x": len(x_unique),
            "n_grid_y": len(y_unique),
        }

        return self._grid_info

    def create_land_mask(self) -> np.ndarray:
        """
        陸地マスクを作成

        Returns
        -------
        np.ndarray
            陸地マスク（Trueが陸地）
        """
        if self._df_elevation is None:
            raise ValueError("データが読み込まれていません。load_data()を先に実行してください。")

        grid_coords = self.create_grid_coords()
        lon_mesh, lat_mesh = self.create_meshgrid()

        # 地形マスクの作成
        land_points = self._df_elevation.select(
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

    def preprocess_obsidian_data(
        self, target_period: int, target_origin: str
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        黒曜石データを前処理

        Parameters
        ----------
        target_period : int
            対象時期
        target_origin : str
            対象産地カテゴリ

        Returns
        -------
        counts : np.ndarray
            各遺跡での出土数
        target_counts : np.ndarray
            対象産地の出土数
        """
        if self._df_obsidian is None:
            raise ValueError("データが読み込まれていません。load_data()を先に実行してください。")

        # 全遺跡IDのリストを取得
        max_site_id = self._df_obsidian["遺跡ID"].max()

        # 対象時期のデータのみ抽出
        period_df = self._df_obsidian.filter(pl.col("時期") == target_period)

        # 全体のカウント
        counts = (
            period_df.group_by("遺跡ID")
            .agg([pl.len().alias("count")])
            .join(
                pl.DataFrame({"遺跡ID": np.arange(max_site_id + 1)}),
                on="遺跡ID",
                how="right",
            )
            .fill_null(0)
            .sort("遺跡ID")["count"]
            .to_numpy()
        )

        # 対象産地のカウント
        target_counts = (
            period_df.filter(pl.col("産地カテゴリ") == target_origin)
            .group_by("遺跡ID")
            .agg([pl.len().alias("count")])
            .join(
                pl.DataFrame({"遺跡ID": np.arange(max_site_id + 1)}),
                on="遺跡ID",
                how="right",
            )
            .fill_null(0)
            .sort("遺跡ID")["count"]
            .to_numpy()
        )

        return counts, target_counts

    def create_site_coords(self) -> np.ndarray:
        """
        遺跡の座標をラジアンに変換

        Returns
        -------
        np.ndarray
            遺跡の座標（ラジアン）
        """
        if self._df_obsidian is None:
            raise ValueError("データが読み込まれていません。load_data()を先に実行してください。")

        # 座標をラジアンに変換
        coords = (
            self._df_obsidian.select(
                [
                    (pl.col("遺跡ID")),
                    (pl.col("緯度") * np.pi / 180).alias("lat_rad"),
                    (pl.col("経度") * np.pi / 180).alias("lon_rad"),
                ]
            )
            .unique(subset=["遺跡ID"])
            .sort("遺跡ID")
        )

        # 座標と出土数を numpy 配列に変換
        site_coords = np.column_stack(
            [coords["lat_rad"].to_numpy(), coords["lon_rad"].to_numpy()]
        )

        return site_coords

    def create_grid_coords(self) -> np.ndarray:
        """
        グリッドの座標を作成

        Returns
        -------
        np.ndarray
            グリッドの座標（ラジアン）
        """
        lon_mesh, lat_mesh = self.create_meshgrid()
        
        # グリッドの位置座標: (グリッド数, 2)
        grid_coords = np.column_stack(
            [lat_mesh.ravel() * np.pi / 180, lon_mesh.ravel() * np.pi / 180]
        )
        
        return grid_coords

    def create_meshgrid(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        メッシュグリッドを作成

        Returns
        -------
        lon_mesh : np.ndarray
            経度メッシュ
        lat_mesh : np.ndarray
            緯度メッシュ
        """
        if self._df_elevation is None:
            raise ValueError("データが読み込まれていません。load_data()を先に実行してください。")

        # メッシュグリッドを作成
        lon_mesh, lat_mesh = np.meshgrid(
            self._df_elevation["x"].unique().sort(),  # 経度の一意な値
            self._df_elevation["y"].unique().sort(),  # 緯度の一意な値
        )

        return lon_mesh, lat_mesh

    @property
    def df_elevation(self) -> pl.DataFrame:
        """標高データフレーム"""
        if self._df_elevation is None:
            raise ValueError("データが読み込まれていません。load_data()を先に実行してください。")
        return self._df_elevation

    @property
    def df_obsidian(self) -> pl.DataFrame:
        """黒曜石データフレーム"""
        if self._df_obsidian is None:
            raise ValueError("データが読み込まれていません。load_data()を先に実行してください。")
        return self._df_obsidian

    @property
    def df_sites(self) -> pl.DataFrame:
        """遺跡データフレーム"""
        if self._df_sites is None:
            raise ValueError("データが読み込まれていません。load_data()を先に実行してください。")
        return self._df_sites