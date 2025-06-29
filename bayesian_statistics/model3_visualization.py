"""
Model 3rd の可視化クラス
結果の可視化を行う
"""

from typing import Optional, Tuple

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import polars as pl


class ObsidianVisualizer:
    """黒曜石分析結果の可視化"""

    @staticmethod
    def create_grid_visualization(df: pl.DataFrame) -> plt.Figure:
        """
        グリッドの可視化

        Parameters
        ----------
        df : pl.DataFrame
            グリッドデータ

        Returns
        -------
        plt.Figure
            作成した図
        """
        # グリッドの次元を取得
        x_max = df["grid_x"].max() + 1
        y_max = df["grid_y"].max() + 1

        # 2D配列に変換
        grid = np.zeros((y_max, x_max))
        for row in df.iter_rows():
            x, y, cost = row
            grid[y, x] = cost

        # プロットの作成
        plt.figure(figsize=(12, 8))

        # pcolormeshでヒートマップを描画
        im = plt.pcolormesh(grid, cmap="viridis", shading="auto")

        # カラーバーを追加
        plt.colorbar(im, label="Cost(分)")

        # 軸ラベルを設定
        plt.xlabel("Grid X")
        plt.ylabel("Grid Y")
        plt.title("Tobler's hiking function")

        # グリッドを表示
        plt.grid(True)

        return plt.gcf()

    @staticmethod
    def plot_contour(
        df: pl.DataFrame,
        x_col: str = "x",
        y_col: str = "y",
        value_col: str = "cost_kouzu",
        figsize: Tuple[int, int] = (12, 8),
        plot_probability: bool = False,
        n_levels: int = 30,
        cmap: str = "Blues",
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        コンター図を作成

        Parameters
        ----------
        df : pl.DataFrame
            データフレーム
        x_col : str
            x軸のカラム名
        y_col : str
            y軸のカラム名
        value_col : str
            値のカラム名
        figsize : Tuple[int, int]
            図のサイズ
        plot_probability : bool
            確率としてプロット（0-1の範囲）
        n_levels : int
            等高線のレベル数
        cmap : str
            カラーマップ

        Returns
        -------
        fig : plt.Figure
            図オブジェクト
        ax : plt.Axes
            軸オブジェクト
        """
        # まず海陸判定のマスクを作成
        mask = df.with_columns((~pl.col("is_sea")).alias("is_not_sea")).pivot(
            values="is_not_sea", on=x_col, index=y_col
        )

        # 値のピボットテーブルを作成
        grid_data = df.pivot(values=value_col, on=x_col, index=y_col)

        # マスクを適用（y_col列は保持）
        grid_data = grid_data.with_columns(
            [
                pl.col(col) * mask.get_column(col)
                for col in grid_data.columns
                if col != y_col
            ]
        )

        # メッシュグリッドの作成
        x_mesh = np.array(grid_data.columns[1:], dtype=float)
        y_mesh = np.array(grid_data.to_numpy()[:, 0], dtype=float)
        values_mesh = grid_data.to_numpy()[:, 1:]

        # プロット作成
        fig, ax = plt.subplots(figsize=figsize)

        if plot_probability:
            # 確率表示モード（0-1の範囲）
            contour = ax.contourf(
                x_mesh,
                y_mesh,
                values_mesh,
                levels=np.linspace(0, 1, n_levels + 1),
                cmap=cmap,
                alpha=0.7,
                vmin=0,
                vmax=1,
            )
        else:
            # 通常モード（データをそのまま使用）
            contour = ax.contourf(
                x_mesh, y_mesh, values_mesh, levels=n_levels, cmap=cmap, alpha=0.7
            )

        # カラーバーの追加
        if plot_probability:
            # カラーバー（0-1の範囲に固定）
            plt.colorbar(contour, ax=ax, label="Ratio", ticks=np.linspace(0, 1, 6))
        else:
            plt.colorbar(contour, ax=ax)

        # ラベルの設定
        ax.set_xlabel("経度")
        ax.set_ylabel("緯度")

        return fig, ax

    @staticmethod
    def plot_ratio_map(
        df_elevation: pl.DataFrame,
        df_sites: pl.DataFrame,
        period: int,
        origin: str,
        time_period_names: dict,
        figsize: Tuple[int, int] = (12, 8),
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        産地構成比の地図をプロット

        Parameters
        ----------
        df_elevation : pl.DataFrame
            標高データフレーム
        df_sites : pl.DataFrame
            遺跡データフレーム
        period : int
            時期
        origin : str
            産地
        time_period_names : dict
            時期名の辞書
        figsize : Tuple[int, int]
            図のサイズ

        Returns
        -------
        fig : plt.Figure
            図オブジェクト
        ax : plt.Axes
            軸オブジェクト
        """
        # コンター図の作成
        fig, ax = ObsidianVisualizer.plot_contour(
            df_elevation,
            value_col=f"ratio_{period}_{origin}",
            plot_probability=True,
            figsize=figsize,
        )

        # 境界の描画
        boundary_df = df_elevation.filter(
            pl.col("is_sea") == 0, pl.col("average_elevation").is_null()
        )
        ax.scatter(boundary_df["x"], boundary_df["y"], c="black", s=0.001)

        # 遺跡のプロット
        ax.scatter(
            df_sites["経度"],
            df_sites["緯度"],
            c=df_sites[f"比率_{period}_{origin}"],
            cmap="Blues",
            edgecolors="black",
            linewidths=0.5,
            vmin=0,
            vmax=1,
        )

        ax.set_title(f"黒曜石の産地構成比 ({time_period_names[period]}, {origin})")

        return fig, ax

    @staticmethod
    def plot_ratio_map_from_result(
        df_elevation: pl.DataFrame,
        df_sites: pl.DataFrame,
        ratio_mesh: np.ndarray,
        ratio_sites: np.ndarray,
        period: int,
        origin: str,
        time_period_names: dict,
        figsize: Tuple[int, int] = (12, 8),
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        予測結果から直接産地構成比の地図をプロット

        Parameters
        ----------
        df_elevation : pl.DataFrame
            標高データフレーム
        df_sites : pl.DataFrame
            遺跡データフレーム
        ratio_mesh : np.ndarray
            グリッド上の比率（1次元配列）
        ratio_sites : np.ndarray
            遺跡での比率（1次元配列）
        period : int
            時期
        origin : str
            産地
        time_period_names : dict
            時期名の辞書
        figsize : Tuple[int, int]
            図のサイズ

        Returns
        -------
        fig : plt.Figure
            図オブジェクト
        ax : plt.Axes
            軸オブジェクト
        """
        # ratio_meshを2次元に戻す
        lon_mesh, lat_mesh = np.meshgrid(
            df_elevation["x"].unique().sort(),  # 経度の一意な値
            df_elevation["y"].unique().sort(),  # 緯度の一意な値
        )
        ratio_2d = ratio_mesh.reshape(lon_mesh.shape)

        # プロット作成
        fig, ax = plt.subplots(figsize=figsize)

        # コンター図の作成
        contour = ax.contourf(
            lon_mesh,
            lat_mesh,
            ratio_2d,
            levels=np.linspace(0, 1, 31),
            cmap="Blues",
            alpha=0.7,
            vmin=0,
            vmax=1,
        )

        # カラーバーの追加
        plt.colorbar(contour, ax=ax, label="Ratio", ticks=np.linspace(0, 1, 6))

        # 境界の描画
        boundary_df = df_elevation.filter(
            pl.col("is_sea") == 0, pl.col("average_elevation").is_null()
        )
        ax.scatter(boundary_df["x"], boundary_df["y"], c="black", s=0.001)

        # 遺跡のプロット
        ax.scatter(
            df_sites["経度"],
            df_sites["緯度"],
            c=ratio_sites,
            cmap="Blues",
            edgecolors="black",
            linewidths=0.5,
            vmin=0,
            vmax=1,
        )

        # ラベルの設定
        ax.set_xlabel("経度")
        ax.set_ylabel("緯度")
        ax.set_title(f"黒曜石の産地構成比 ({time_period_names[period]}, {origin})")

        return fig, ax

    @staticmethod
    def plot_site_probability(
        df_elevation: pl.DataFrame,
        df_sites: pl.DataFrame,
        figsize: Tuple[int, int] = (12, 8),
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        遺跡存在確率の地図をプロット

        Parameters
        ----------
        df_elevation : pl.DataFrame
            標高データフレーム
        df_sites : pl.DataFrame
            遺跡データフレーム
        figsize : Tuple[int, int]
            図のサイズ

        Returns
        -------
        fig : plt.Figure
            図オブジェクト
        ax : plt.Axes
            軸オブジェクト
        """
        # コンター図の作成
        fig, ax = ObsidianVisualizer.plot_contour(
            df_elevation,
            value_col="site_probability",
            plot_probability=True,
            cmap="magma",
            figsize=figsize,
        )

        # 境界の描画
        boundary_df = df_elevation.filter(
            pl.col("is_sea") == 0, pl.col("average_elevation").is_null()
        )
        ax.scatter(boundary_df["x"], boundary_df["y"], c="black", s=0.001)

        # 遺跡のプロット
        ax.scatter(
            df_sites["経度"],
            df_sites["緯度"],
            c="white",
            edgecolors="black",
            linewidths=0.5,
            vmin=0,
            vmax=1,
        )

        ax.set_title("遺跡の存在確率")

        return fig, ax

    @staticmethod
    def plot_mcmc_diagnostics(
        idata: az.InferenceData,
        figsize: Optional[Tuple[int, int]] = None,
    ) -> None:
        """
        MCMC診断プロットを作成

        Parameters
        ----------
        idata : az.InferenceData
            InferenceDataオブジェクト
        figsize : Optional[Tuple[int, int]]
            図のサイズ
        """
        # ArviZのスタイルを設定
        az.style.use("arviz-doc")

        # トレースプロット
        az.plot_trace(idata, figsize=figsize)
        plt.show()

        # 事後分布のプロット
        az.plot_posterior(idata, figsize=figsize)
        plt.show()

        # 自己相関プロット
        az.plot_autocorr(idata, figsize=figsize)
        plt.show()

        # 統計的な概要情報の表示
        summary_df = az.summary(idata)
        print(summary_df)

        # スタイルをリセット
        az.style.use("default")

    @staticmethod
    def plot_all_periods_origins(
        df_elevation: pl.DataFrame,
        df_sites: pl.DataFrame,
        time_period_names: dict,
        origin_order: list,
        save_dir: Optional[str] = None,
    ) -> None:
        """
        全時期・全産地の図を作成

        Parameters
        ----------
        df_elevation : pl.DataFrame
            標高データフレーム
        df_sites : pl.DataFrame
            遺跡データフレーム
        time_period_names : dict
            時期名の辞書
        origin_order : list
            産地のリスト
        save_dir : Optional[str]
            保存ディレクトリ（Noneの場合は保存しない）
        """
        for period in time_period_names.keys():
            for origin in origin_order[:-1]:  # "その他"を除外
                fig, ax = ObsidianVisualizer.plot_ratio_map(
                    df_elevation,
                    df_sites,
                    period,
                    origin,
                    time_period_names,
                )
                
                if save_dir is not None:
                    filename = f"ratio_{period}_{origin}.png"
                    fig.savefig(f"{save_dir}/{filename}", dpi=300, bbox_inches="tight")
                    plt.close(fig)
                else:
                    plt.show()