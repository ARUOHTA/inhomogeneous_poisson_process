#!/usr/bin/env python3
"""
修正版ベイズ空間多項分布回帰モデル
可視化の問題（縦線と均一予測）を修正
"""

from typing import List, Tuple

import arviz as az
import japanize_matplotlib
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pymc as pm

from bayesian_statistics.model3_config import Model3Config
from bayesian_statistics.model3_preprocessing import ObsidianDataPreprocessor
from bayesian_statistics.model3_visualization import ObsidianVisualizer

japanize_matplotlib.japanize()


class FixedBayesianSpatialMultinomial:
    """修正版ベイズ空間多項ロジット回帰モデル"""

    def __init__(self, n_origins: int, prior_sigma: float = 0.5):
        self.n_origins = n_origins
        self.prior_sigma = prior_sigma
        self.model = None
        self.trace = None
        self.origins = None
        self.X_mean = None
        self.X_std = None

    def prepare_data(
        self,
        preprocessor: ObsidianDataPreprocessor,
        target_period: int,
        time_periods: dict,
        origins_filtered: List[str],
        config,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, pl.DataFrame]:
        """データ準備（修正版）"""
        print("修正版ベイズモデル用データ準備")

        # まず単一変数でテスト（標高のみ）
        spatial_variables = ["average_elevation"]

        # 重要："その他"を含む完全なoriginsを使用
        data = preprocessor.create_X_Y(spatial_variables, time_periods, config.origins)
        X_all = data["X"]
        Y_all = data["Y"]

        print(f"全データ形状: X={X_all.shape}, Y={Y_all.shape}")

        # 対象時期のデータを正しく抽出
        period_start_idx = target_period * len(origins_filtered)
        period_end_idx = period_start_idx + len(origins_filtered)
        Y_period = Y_all[:, period_start_idx:period_end_idx]

        print(f"対象時期のデータ形状: Y_period={Y_period.shape}")
        print(
            f"期間2の産地別合計: {[Y_period[:, i].sum() for i in range(len(origins_filtered))]}"
        )

        # 空間変数のみ使用（遺跡IDと緯度経度を除く）
        # X_all構造: [site_id, 緯度, 経度, 変数1, 変数2, ...]
        print(f"X_all形状: {X_all.shape}")
        print(f"spatial_variables: {spatial_variables}")

        # 最初の2列は経度(lon), 緯度(lat)なので、空間変数のみ抽出
        X_spatial = X_all[:, 2 : 2 + len(spatial_variables)]  # 指定した変数数だけ取得
        # データがある遺跡のみ
        has_data_mask = Y_period.sum(axis=1) > 0
        X_sites = X_spatial[has_data_mask]
        Y_sites = Y_period[has_data_mask]

        # インターセプト追加
        X_sites_with_intercept = np.column_stack([np.ones(len(X_sites)), X_sites])

        print(
            f"学習用遺跡データ形状: X={X_sites_with_intercept.shape}, Y={Y_sites.shape}"
        )
        print(
            f"学習用データの産地別合計: {[Y_sites[:, i].sum() for i in range(len(origins_filtered))]}"
        )

        # グリッドデータ（陸地のみ、サンプリングして高速化）
        df_viz = preprocessor.df_elevation.clone()
        land_mask = (~df_viz["is_sea"]) & (~df_viz["average_elevation"].is_null())
        df_land = df_viz.filter(land_mask)

        # 大きすぎる場合はサンプリング
        if len(df_land) > 50000:
            df_land = df_land.sample(n=50000, seed=42)
            print(f"グリッドデータをサンプリング: {len(df_land)}点")

        # グリッド用特徴量（遺跡データと同じ構造にする）
        grid_features = []

        # インターセプト追加
        grid_features.append(np.ones(len(df_land)))

        # 指定した空間変数を順次追加
        for var_name in spatial_variables:
            values = df_land[var_name].to_numpy()
            values = np.where(np.isnan(values), np.nanmean(values), values)
            grid_features.append(values)

        X_grid = np.column_stack(grid_features)

        print(f"グリッドデータ形状: {X_grid.shape}")

        return X_sites_with_intercept, Y_sites, X_grid, df_land

    def build_model(self, X: np.ndarray, Y: np.ndarray) -> None:
        """ベイズモデル構築（修正版）"""
        n_sites, n_vars = X.shape
        n_origins = Y.shape[1]

        # データ標準化
        self.X_mean = X.mean(axis=0)
        self.X_std = X.std(axis=0)
        self.X_std[self.X_std == 0] = 1
        X_scaled = (X - self.X_mean) / self.X_std

        print(
            f"ベイズモデル構築: 観測数={n_sites}, 変数数={n_vars}, 産地数={n_origins}"
        )

        with pm.Model() as model:
            # より弱い事前分布で空間変動を許可
            beta = pm.Normal(
                "beta", mu=0, sigma=self.prior_sigma, shape=(n_vars, n_origins - 1)
            )

            # 線形予測子
            eta = pm.math.dot(X_scaled, beta)
            eta_padded = pm.math.concatenate([eta, pm.math.zeros((n_sites, 1))], axis=1)

            # 多項分布
            p = pm.math.softmax(eta_padded, axis=1)
            n_total = pm.math.sum(Y, axis=1)
            likelihood = pm.Multinomial("likelihood", n=n_total, p=p, observed=Y)

        self.model = model

    def fit(self, X: np.ndarray, Y: np.ndarray) -> None:
        """モデル学習"""
        print("ベイズ推定実行中...")
        self.build_model(X, Y)

        with self.model:
            self.trace = pm.sample(
                draws=1000,
                tune=500,
                return_inferencedata=True,
                random_seed=42,
                progressbar=True,
            )

        rhat = az.rhat(self.trace).to_array()
        print(f"収束診断: 最大R-hat = {float(rhat.max()):.4f}")

        # 係数の要約
        summary = az.summary(self.trace)
        print("係数の事後分布:")
        print(summary)

    def predict_probabilities(self, X_new: np.ndarray) -> np.ndarray:
        """事後予測（修正版）"""
        print(f"事後予測実行: {len(X_new)}地点")

        X_scaled = (X_new - self.X_mean) / self.X_std

        beta_samples = self.trace.posterior["beta"].values
        beta_flat = beta_samples.reshape(
            -1, beta_samples.shape[-2], beta_samples.shape[-1]
        )

        # より多くのサンプルで予測精度向上
        n_samples = min(500, len(beta_flat))
        idx = np.random.choice(len(beta_flat), n_samples, replace=False)
        beta_selected = beta_flat[idx]

        predictions = []
        for beta_sample in beta_selected:
            eta = X_scaled @ beta_sample
            eta_padded = np.column_stack([eta, np.zeros(len(X_new))])
            p = self._softmax(eta_padded)
            predictions.append(p)

        mean_predictions = np.array(predictions).mean(axis=0)

        # 予測の統計情報
        print("予測結果の統計:")
        for i, origin in enumerate(self.origins):
            values = mean_predictions[:, i]
            print(
                f"  {origin}: 範囲=[{values.min():.4f}, {values.max():.4f}], 平均={values.mean():.4f}"
            )

        return mean_predictions

    def _softmax(self, x: np.ndarray) -> np.ndarray:
        """数値安定なソフトマックス"""
        x_clipped = np.clip(x, -50, 50)
        x_max = np.max(x_clipped, axis=1, keepdims=True)
        exp_x = np.exp(x_clipped - x_max)
        return exp_x / np.sum(exp_x, axis=1, keepdims=True)


def create_fixed_visualization_data(
    df_land: pl.DataFrame,
    probabilities_grid: np.ndarray,
    target_period: int,
    origins_filtered: List[str],
) -> pl.DataFrame:
    """修正版可視化データ作成"""
    print("修正版可視化データ作成中...")

    df_viz = df_land.clone()

    for i, origin in enumerate(origins_filtered):
        col_name = f"ratio_{target_period}_{origin}"
        df_viz = df_viz.with_columns(
            pl.Series(col_name, probabilities_grid[:, i]).alias(col_name)
        )

        values = probabilities_grid[:, i]
        print(f"  {origin}: 範囲=[{values.min():.4f}, {values.max():.4f}]")

        # 空間変動があるかチェック
        spatial_var = np.var(values)
        print(f"    空間変動（分散）: {spatial_var:.6f}")

    return df_viz


def visualize_fixed_origin(
    df_elevation_viz: pl.DataFrame,
    preprocessor: ObsidianDataPreprocessor,
    target_period: int,
    origin: str,
    config: Model3Config,
) -> None:
    """修正版可視化"""
    print(f"  {origin}の修正版可視化中...")

    # 遺跡データ準備
    spatial_variables = ["average_elevation"]
    origins_filtered = [org for org in config.origins if org != "その他"]
    data = preprocessor.create_X_Y(
        spatial_variables, config.time_periods, config.origins
    )
    Y_all = data["Y"]

    period_start_idx = target_period * len(origins_filtered)
    period_end_idx = period_start_idx + len(origins_filtered)
    Y_period = Y_all[:, period_start_idx:period_end_idx]

    has_data_mask = Y_period.sum(axis=1) > 0
    import numpy as np

    sites_indices = np.where(has_data_mask)[0]

    sites_with_data = preprocessor.df_sites.filter(
        pl.int_range(pl.len()).is_in(sites_indices.tolist())
    )

    # 実際の構成比
    Y_sites_subset = Y_period[has_data_mask]
    total_counts = Y_sites_subset.sum(axis=1)
    origin_idx = origins_filtered.index(origin)
    actual_ratios = Y_sites_subset[:, origin_idx] / (total_counts + 1e-12)

    sites_with_data = sites_with_data.with_columns(
        pl.Series(f"比率_{target_period}_{origin}", actual_ratios).alias(
            f"比率_{target_period}_{origin}"
        )
    )

    try:
        # contour plot使用、ただし縦線問題を回避するためデータを事前整理

        # グリッドデータを規則的な格子に変換
        x_unique = sorted(df_elevation_viz["x"].unique())
        y_unique = sorted(df_elevation_viz["y"].unique())

        print(f"    グリッド範囲: x={len(x_unique)}, y={len(y_unique)}")

        # 規則的な格子データに変換（欠損値は補間）
        import numpy as np
        from scipy.interpolate import griddata

        # 既存データポイント
        points = np.column_stack(
            [df_elevation_viz["x"].to_numpy(), df_elevation_viz["y"].to_numpy()]
        )
        values = df_elevation_viz[f"ratio_{target_period}_{origin}"].to_numpy()

        # 規則的なグリッド作成
        x_grid, y_grid = np.meshgrid(x_unique, y_unique)

        # 補間でデータ補完（縦線の原因となる欠損を除去）
        values_interpolated = griddata(
            points, values, (x_grid, y_grid), method="linear", fill_value=0
        )

        # 海のマスク処理を追加 - 元の標高データから海域情報を取得
        df_full_elevation = preprocessor.df_elevation

        # 全標高データから海域を特定
        all_points_full = np.column_stack(
            [df_full_elevation["x"].to_numpy(), df_full_elevation["y"].to_numpy()]
        )
        is_sea_values_full = df_full_elevation["is_sea"].to_numpy()

        # 海域の情報を同じグリッドに補間
        sea_interpolated = griddata(
            all_points_full,
            is_sea_values_full,
            (x_grid, y_grid),
            method="nearest",
            fill_value=1,  # 欠損部分は海とする
        )

        # 海域部分（is_sea=1）は値を0にして白く表示
        values_interpolated = np.where(sea_interpolated >= 0.5, 0, values_interpolated)

        # contour plot作成
        fig, ax = plt.subplots(figsize=(12, 8))

        contour = ax.contourf(
            x_grid,
            y_grid,
            values_interpolated,
            levels=np.linspace(0, 1, 31),
            cmap="Blues",
            alpha=0.7,
            vmin=0,
            vmax=1,
        )

        # カラーバー
        plt.colorbar(contour, ax=ax, label="Ratio", ticks=np.linspace(0, 1, 6))

        # 海岸線の境界描画（海と陸の境界）
        boundary_df = df_full_elevation.filter(
            pl.col("is_sea") == 0, pl.col("average_elevation").is_null()
        )
        if len(boundary_df) > 0:
            ax.scatter(boundary_df["x"], boundary_df["y"], c="black", s=0.001)

        # 遺跡プロット
        ax.scatter(
            sites_with_data["経度"],
            sites_with_data["緯度"],
            c=sites_with_data[f"比率_{target_period}_{origin}"],
            cmap="Blues",
            edgecolors="black",
            linewidths=0.5,
            vmin=0,
            vmax=1,
            s=30,
        )

        ax.set_xlabel("経度")
        ax.set_ylabel("緯度")
        ax.set_title(f"線形回帰: {config.time_periods[target_period]} - {origin}")

        output_path = f"output/fixed_bayesian_map_{target_period}_{origin}.png"
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"    保存完了: {output_path}")

        plt.show()

    except Exception as e:
        print(f"    エラー: {e}")
        import traceback

        traceback.print_exc()


def main():
    """修正版ベイズモデル可視化のメイン実行"""
    print("=== 修正版ベイズ空間多項分布回帰モデルの可視化 ===")

    config = Model3Config(data_dir="/home/ohta/dev/bayesian_statistics/data/")
    preprocessor = ObsidianDataPreprocessor(config.data_dir)
    preprocessor.load_data()

    target_period = 2
    origins_filtered = [origin for origin in config.origins if origin != "その他"]

    print("分析設定:")
    print(f"  対象時期: {target_period} ({config.time_periods[target_period]})")
    print(f"  対象産地: {origins_filtered}")

    # 修正版ベイズモデル初期化
    model = FixedBayesianSpatialMultinomial(n_origins=len(origins_filtered))
    model.origins = origins_filtered

    # データ準備
    X_sites, Y_sites, X_grid, df_land = model.prepare_data(
        preprocessor, target_period, config.time_periods, origins_filtered, config
    )

    # モデル学習（インターセプト付きデータを使用）
    model.fit(X_sites, Y_sites)

    # 予測実行
    probabilities_grid = model.predict_probabilities(X_grid)

    # 可視化データ準備
    df_elevation_viz = create_fixed_visualization_data(
        df_land, probabilities_grid, target_period, origins_filtered
    )

    # 全産地を可視化
    print("\n全産地の修正版可視化:")
    for origin in origins_filtered:
        visualize_fixed_origin(
            df_elevation_viz, preprocessor, target_period, origin, config
        )

    print("\n修正版ベイズモデル可視化完了!")


if __name__ == "__main__":
    main()
