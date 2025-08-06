"""
黒曜石分析結果の可視化クラス
既存のmodel3_visualization.pyから移行
"""

from typing import Any, Dict, Optional, Tuple

import arviz as az
import japanize_matplotlib
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from .base_visualizer import BaseVisualizer

japanize_matplotlib.japanize()


class ObsidianVisualizer(BaseVisualizer):
    """黒曜石分析結果の可視化"""

    def plot_prediction_map(
        self,
        data: pl.DataFrame,
        df_elevation: pl.DataFrame,
        df_sites: pl.DataFrame,
        period: int,
        origin: str,
        **kwargs,
    ) -> plt.Figure:
        """
        予測結果の地図プロット（BaseVisualizer実装）

        Parameters
        ----------
        data : pl.DataFrame
            予測結果データ
        df_elevation : pl.DataFrame
            標高データ
        df_sites : pl.DataFrame
            遺跡データ
        period : int
            時期
        origin : str
            産地名
        **kwargs
            追加のプロットパラメータ

        Returns
        -------
        plt.Figure
            作成された図
        """
        return self.plot_ratio_map(
            df_elevation=df_elevation,
            df_sites=df_sites,
            period=period,
            origin=origin,
            **kwargs,
        )

    def plot_model_diagnostics(self, model_results: Dict[str, Any]) -> plt.Figure:
        """
        モデル診断プロット（BaseVisualizer実装）

        Parameters
        ----------
        model_results : Dict[str, Any]
            モデルの学習結果（'trace'キーを含む）

        Returns
        -------
        plt.Figure
            診断プロット
        """
        if "trace" in model_results:
            return self.plot_mcmc_diagnostics(model_results["trace"])
        else:
            # 簡単な診断情報を表示
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(
                0.5,
                0.5,
                "No MCMC trace available for diagnostics",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title("Model Diagnostics")
            return fig

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
        # データを取得
        x = df[x_col].to_numpy()
        y = df[y_col].to_numpy()
        z = df[value_col].to_numpy()

        # 海のマスクを作成（zが0またはnanの場合は海とみなす）
        mask = (z == 0) | np.isnan(z)

        # グリッドのサイズを取得
        x_unique = np.unique(x)
        y_unique = np.unique(y)
        X, Y = np.meshgrid(x_unique, y_unique)

        # Zを2D配列に変換
        Z = np.full(X.shape, np.nan)
        for i, row in enumerate(df.iter_rows()):
            x_val, y_val, z_val = (
                row[df.get_column_index(x_col)],
                row[df.get_column_index(y_col)],
                row[df.get_column_index(value_col)],
            )
            xi = np.where(x_unique == x_val)[0][0]
            yi = np.where(y_unique == y_val)[0][0]
            Z[yi, xi] = z_val

        # プロット作成
        fig, ax = plt.subplots(figsize=figsize)

        # 確率プロットの場合は0-1の範囲を使用
        if plot_probability:
            levels = np.linspace(0, 1, n_levels)
        else:
            valid_z = z[~mask]
            if len(valid_z) > 0:
                levels = np.linspace(valid_z.min(), valid_z.max(), n_levels)
            else:
                levels = np.linspace(0, 1, n_levels)

        # コンター図を描画
        contour = ax.contourf(X, Y, Z, levels=levels, cmap=cmap, extend="both")

        # 海の部分を灰色で塗りつぶし
        sea_mask = np.isnan(Z) | (Z == 0)
        ax.contourf(
            X,
            Y,
            sea_mask.astype(float),
            levels=[0.5, 1.5],
            colors=["lightgray"],
            alpha=0.7,
        )

        # カラーバーを追加
        cbar = plt.colorbar(contour, ax=ax)
        if plot_probability:
            cbar.set_label("Probability")
        else:
            cbar.set_label("Value")

        # 軸ラベル
        ax.set_xlabel("X coordinate")
        ax.set_ylabel("Y coordinate")
        ax.set_aspect("equal")

        return fig, ax

    @staticmethod
    def plot_ratio_map(
        df_elevation: pl.DataFrame,
        df_sites: pl.DataFrame,
        period: int,
        origin: str,
        figsize: Tuple[int, int] = (12, 10),
        n_levels: int = 20,
        cmap: str = "Reds",
        alpha: float = 0.7,
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        産地構成比のマップをプロット

        Parameters
        ----------
        df_elevation : pl.DataFrame
            標高データ
        df_sites : pl.DataFrame
            遺跡データ
        period : int
            時期のインデックス
        origin : str
            産地名
        figsize : Tuple[int, int]
            図のサイズ
        n_levels : int
            等高線のレベル数
        cmap : str
            カラーマップ
        alpha : float
            透明度

        Returns
        -------
        fig : plt.Figure
            図オブジェクト
        ax : plt.Axes
            軸オブジェクト
        """
        # 該当する時期・産地の比率カラムを特定
        ratio_col = f"ratio_{period}_{origin}"

        if ratio_col not in df_elevation.columns:
            raise ValueError(f"Column {ratio_col} not found in df_elevation")

        # コンター図を作成
        fig, ax = ObsidianVisualizer.plot_contour(
            df_elevation,
            x_col="x",
            y_col="y",
            value_col=ratio_col,
            figsize=figsize,
            plot_probability=True,
            n_levels=n_levels,
            cmap=cmap,
        )

        # 遺跡データをプロット
        # 該当時期の遺跡のみフィルタ
        period_sites = df_sites.filter(pl.col(f"時期_{period}") == 1)

        if len(period_sites) > 0:
            # 遺跡の位置をプロット
            ax.scatter(
                period_sites["経度"].to_numpy(),
                period_sites["緯度"].to_numpy(),
                c="black",
                s=50,
                marker="o",
                alpha=alpha,
                edgecolors="white",
                linewidth=1,
            )

        # タイトルを設定
        time_periods = {0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"}
        period_name = time_periods.get(period, f"期間{period}")
        ax.set_title(f"{period_name}・{origin}の構成比")

        return fig, ax

    @staticmethod
    def plot_ratio_map_from_result(
        df_elevation: pl.DataFrame,
        df_sites: pl.DataFrame,
        ratio_mesh: np.ndarray,
        ratio_sites: np.ndarray,
        period: int,
        origin: str,
        grid_x: np.ndarray,
        grid_y: np.ndarray,
        figsize: Tuple[int, int] = (12, 10),
        n_levels: int = 20,
        cmap: str = "Reds",
        alpha: float = 0.7,
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        予測結果から産地構成比のマップをプロット

        Parameters
        ----------
        df_elevation : pl.DataFrame
            標高データ
        df_sites : pl.DataFrame
            遺跡データ
        ratio_mesh : np.ndarray
            メッシュ上の比率
        ratio_sites : np.ndarray
            遺跡での比率
        period : int
            時期のインデックス
        origin : str
            産地名
        grid_x : np.ndarray
            グリッドX座標
        grid_y : np.ndarray
            グリッドY座標
        figsize : Tuple[int, int]
            図のサイズ
        n_levels : int
            等高線のレベル数
        cmap : str
            カラーマップ
        alpha : float
            透明度

        Returns
        -------
        fig : plt.Figure
            図オブジェクト
        ax : plt.Axes
            軸オブジェクト
        """
        # プロット作成
        fig, ax = plt.subplots(figsize=figsize)

        # メッシュ上の予測結果をプロット
        levels = np.linspace(0, 1, n_levels)
        contour = ax.contourf(
            grid_x, grid_y, ratio_mesh, levels=levels, cmap=cmap, extend="both"
        )

        # 海の部分を灰色で塗りつぶし（比率が0の部分）
        sea_mask = ratio_mesh == 0
        ax.contourf(
            grid_x,
            grid_y,
            sea_mask.astype(float),
            levels=[0.5, 1.5],
            colors=["lightgray"],
            alpha=0.7,
        )

        # カラーバーを追加
        cbar = plt.colorbar(contour, ax=ax)
        cbar.set_label("Composition Ratio")

        # 遺跡データをプロット
        period_sites = df_sites.filter(pl.col(f"時期_{period}") == 1)

        if len(period_sites) > 0 and len(ratio_sites) > 0:
            # 予測された比率で色付け
            scatter = ax.scatter(
                period_sites["経度"].to_numpy(),
                period_sites["緯度"].to_numpy(),
                c=ratio_sites,
                s=80,
                marker="o",
                alpha=alpha,
                edgecolors="white",
                linewidth=1,
                cmap=cmap,
                vmin=0,
                vmax=1,
            )

        # タイトルを設定
        time_periods = {0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"}
        period_name = time_periods.get(period, f"期間{period}")
        ax.set_title(f"{period_name}・{origin}の構成比（予測結果）")

        # 軸ラベル
        ax.set_xlabel("X coordinate")
        ax.set_ylabel("Y coordinate")
        ax.set_aspect("equal")

        return fig, ax

    @staticmethod
    def plot_site_probability(
        df_elevation: pl.DataFrame,
        df_sites: pl.DataFrame,
        probability_mesh: np.ndarray,
        grid_x: np.ndarray,
        grid_y: np.ndarray,
        period: int,
        figsize: Tuple[int, int] = (12, 10),
        n_levels: int = 20,
        cmap: str = "Blues",
        alpha: float = 0.7,
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        遺跡存在確率のマップをプロット

        Parameters
        ----------
        df_elevation : pl.DataFrame
            標高データ
        df_sites : pl.DataFrame
            遺跡データ
        probability_mesh : np.ndarray
            メッシュ上の確率
        grid_x : np.ndarray
            グリッドX座標
        grid_y : np.ndarray
            グリッドY座標
        period : int
            時期のインデックス
        figsize : Tuple[int, int]
            図のサイズ
        n_levels : int
            等高線のレベル数
        cmap : str
            カラーマップ
        alpha : float
            透明度

        Returns
        -------
        fig : plt.Figure
            図オブジェクト
        ax : plt.Axes
            軸オブジェクト
        """
        # プロット作成
        fig, ax = plt.subplots(figsize=figsize)

        # 確率の等高線プロット
        levels = np.linspace(0, 1, n_levels)
        contour = ax.contourf(
            grid_x, grid_y, probability_mesh, levels=levels, cmap=cmap, extend="both"
        )

        # 海の部分を灰色で塗りつぶし
        sea_mask = np.isnan(probability_mesh) | (probability_mesh == 0)
        ax.contourf(
            grid_x,
            grid_y,
            sea_mask.astype(float),
            levels=[0.5, 1.5],
            colors=["lightgray"],
            alpha=0.7,
        )

        # カラーバーを追加
        cbar = plt.colorbar(contour, ax=ax)
        cbar.set_label("Site Existence Probability")

        # 実際の遺跡をプロット
        period_sites = df_sites.filter(pl.col(f"時期_{period}") == 1)

        if len(period_sites) > 0:
            ax.scatter(
                period_sites["経度"].to_numpy(),
                period_sites["緯度"].to_numpy(),
                c="red",
                s=60,
                marker="o",
                alpha=alpha,
                edgecolors="white",
                linewidth=1,
                label="Observed Sites",
            )
            ax.legend()

        # タイトルを設定
        time_periods = {0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"}
        period_name = time_periods.get(period, f"期間{period}")
        ax.set_title(f"{period_name}の遺跡存在確率")

        # 軸ラベル
        ax.set_xlabel("X coordinate")
        ax.set_ylabel("Y coordinate")
        ax.set_aspect("equal")

        return fig, ax

    @staticmethod
    def plot_mcmc_diagnostics(
        idata,
        var_names: Optional[list] = None,
        figsize: Tuple[int, int] = (15, 10),
        max_vars: int = 6,
    ) -> plt.Figure:
        """
        MCMC診断プロットを作成

        Parameters
        ----------
        idata : arviz.InferenceData
            MCMCサンプルデータ
        var_names : list, optional
            プロットする変数名のリスト
        figsize : Tuple[int, int]
            図のサイズ
        max_vars : int
            最大表示変数数

        Returns
        -------
        plt.Figure
            診断プロット
        """
        # 変数名が指定されていない場合は自動選択
        if var_names is None:
            all_vars = list(idata.posterior.data_vars)
            var_names = all_vars[:max_vars]  # 最初の数個を選択

        # 診断プロットを作成
        fig, axes = plt.subplots(len(var_names), 2, figsize=figsize)

        if len(var_names) == 1:
            axes = axes.reshape(1, -1)

        for i, var_name in enumerate(var_names):
            # トレースプロット
            az.plot_trace(idata, var_names=[var_name], axes=axes[i])

        plt.tight_layout()
        return fig

    @staticmethod
    def plot_all_periods_origins(
        df_elevation: pl.DataFrame,
        df_sites: pl.DataFrame,
        time_periods: Dict[int, str],
        origins: list,
        output_dir: str = "output",
        figsize: Tuple[int, int] = (12, 10),
        n_levels: int = 20,
        alpha: float = 0.7,
    ) -> Dict[str, plt.Figure]:
        """
        全時期・全産地の比率マップを一括作成

        Parameters
        ----------
        df_elevation : pl.DataFrame
            標高データ
        df_sites : pl.DataFrame
            遺跡データ
        time_periods : Dict[int, str]
            時期の辞書
        origins : list
            産地のリスト
        output_dir : str
            出力ディレクトリ
        figsize : Tuple[int, int]
            図のサイズ
        n_levels : int
            等高線のレベル数
        alpha : float
            透明度

        Returns
        -------
        Dict[str, plt.Figure]
            作成した図の辞書（キー: period_origin）
        """
        figures = {}

        for period, period_name in time_periods.items():
            for origin in origins:
                if origin == "その他":
                    continue

                try:
                    fig, ax = ObsidianVisualizer.plot_ratio_map(
                        df_elevation=df_elevation,
                        df_sites=df_sites,
                        period=period,
                        origin=origin,
                        figsize=figsize,
                        n_levels=n_levels,
                        alpha=alpha,
                    )

                    key = f"{period}_{origin}"
                    figures[key] = fig

                    # ファイル保存
                    plt.savefig(
                        f"{output_dir}/ratio_map_{key}.png",
                        dpi=300,
                        bbox_inches="tight",
                    )
                    plt.close(fig)

                except Exception as e:
                    print(f"Error plotting {period_name} - {origin}: {e}")
                    continue

        return figures

    @staticmethod
    def plot_ratio_map_from_grid_data(
        df_elevation: pl.DataFrame,
        df_sites: pl.DataFrame,
        grid_ratios: pl.DataFrame,
        period: int,
        origin: str,
        figsize: Tuple[int, int] = (12, 10),
        n_levels: int = 20,
        cmap: str = "Reds",
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        grid_ratiosデータから産地構成比のマップをプロット

        Parameters
        ----------
        df_elevation : pl.DataFrame
            標高データ
        df_sites : pl.DataFrame
            遺跡データ
        grid_ratios : pl.DataFrame
            predict_grid_ratios()の結果
        period : int
            時期のインデックス
        origin : str
            産地名
        figsize : Tuple[int, int]
            図のサイズ
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
        # 該当する比率カラムを特定
        ratio_col = f"ratio_{period}_{origin}"

        if ratio_col not in grid_ratios.columns:
            raise ValueError(f"Column {ratio_col} not found in grid_ratios")

        # df_elevationとgrid_ratiosを結合
        df_combined = df_elevation.join(
            grid_ratios.select(["x", "y", ratio_col]), on=["x", "y"], how="left"
        )

        # コンター図を作成
        fig, ax = ObsidianVisualizer.plot_contour(
            df_combined,
            x_col="x",
            y_col="y",
            value_col=ratio_col,
            figsize=figsize,
            plot_probability=True,
            n_levels=n_levels,
            cmap=cmap,
        )

        # 遺跡データをプロット
        # 該当時期の遺跡のみフィルタ
        time_col = f"時期_{period}"
        if time_col in df_sites.columns:
            period_sites = df_sites.filter(pl.col(time_col) == 1)
        else:
            # 時期カラムがない場合は全遺跡を表示
            period_sites = df_sites

        if len(period_sites) > 0:
            # df_sitesは経度・緯度カラムを使用
            ax.scatter(
                period_sites["経度"].to_numpy(),
                period_sites["緯度"].to_numpy(),
                c="black",
                s=20,
                alpha=0.8,
                label="Sites",
            )
            ax.legend()

        ax.set_title(f"Ratio Map: {origin} (Period {period})")

        return fig, ax
