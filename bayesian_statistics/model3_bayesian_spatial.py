"""
Fixed Bayesian Spatial Multinomial モデルのmodel3パイプライン統合版
"""

from typing import Any, Dict, List, Optional, Tuple

import arviz as az
import numpy as np
import polars as pl
import pymc as pm

from .model3_preprocessing import ObsidianDataPreprocessor


class BayesianSpatialMultinomialModel:
    """ベイズ空間多項ロジット回帰モデル"""

    def __init__(
        self,
        variable_names: List[str],
        prior_sigma: float = 0.5,
        n_draws: int = 1000,
        n_tune: int = 500,
        random_seed: int = 42,
    ):
        """
        Parameters
        ----------
        variable_names : List[str]
            使用する説明変数名のリスト
        prior_sigma : float
            事前分布の標準偏差
        n_draws : int
            MCMC描画数
        n_tune : int
            チューニング回数
        random_seed : int
            乱数シード
        """
        self.variable_names = variable_names
        self.prior_sigma = prior_sigma
        self.n_draws = n_draws
        self.n_tune = n_tune
        self.random_seed = random_seed

        self.model = None
        self.trace = None
        self.X_mean = None
        self.X_std = None
        self._fit_results: Optional[Dict[str, Any]] = None

    def _prepare_data_for_period(
        self, preprocessor: ObsidianDataPreprocessor, target_period: int = 0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """特定時期のデータを準備"""

        # 時期設定
        time_periods = {0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"}
        origins = ["神津島", "信州", "箱根", "高原山", "その他"]
        origins_filtered = [org for org in origins if org != "その他"]  # "その他"を除外

        # データ作成
        data = preprocessor.create_X_Y(self.variable_names, time_periods, origins)
        X_all = data["X"]
        Y_all = data["Y"]

        # 対象時期のデータを抽出
        period_start_idx = target_period * len(origins_filtered)
        period_end_idx = period_start_idx + len(origins_filtered)
        Y_period = Y_all[:, period_start_idx:period_end_idx]

        # 空間変数のみ使用（最初の2列は経度・緯度なので除外）
        X_spatial = X_all[:, 2 : 2 + len(self.variable_names)]

        # データがある遺跡のみ
        has_data_mask = Y_period.sum(axis=1) > 0
        X_sites = X_spatial[has_data_mask]
        Y_sites = Y_period[has_data_mask]

        # インターセプト追加
        X_sites_with_intercept = np.column_stack([np.ones(len(X_sites)), X_sites])

        print(f"学習用データ形状: X={X_sites_with_intercept.shape}, Y={Y_sites.shape}")
        print(f"産地別合計: {[Y_sites[:, i].sum() for i in range(Y_sites.shape[1])]}")

        return X_sites_with_intercept, Y_sites

    def _build_model(self, X: np.ndarray, Y: np.ndarray) -> None:
        """ベイズモデル構築"""
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
            # 弱い事前分布で空間変動を許可
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

    def fit(self, preprocessor: ObsidianDataPreprocessor) -> Dict[str, Any]:
        """
        モデルを学習

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ

        Returns
        -------
        Dict[str, Any]
            学習結果
        """
        print("=== ベイズ空間多項回帰学習を開始 ===")

        # データ準備（period=0のみ）
        X_sites, Y_sites = self._prepare_data_for_period(preprocessor, target_period=0)

        # モデル構築
        self._build_model(X_sites, Y_sites)

        # MCMC実行
        print("MCMC推定実行中...")
        with self.model:
            self.trace = pm.sample(
                draws=self.n_draws,
                tune=self.n_tune,
                return_inferencedata=True,
                random_seed=self.random_seed,
                progressbar=True,
            )

        # 収束診断
        rhat = az.rhat(self.trace).to_array()
        print(f"収束診断: 最大R-hat = {float(rhat.max()):.4f}")

        # 係数の要約
        summary = az.summary(self.trace)
        print("係数の事後分布:")
        print(summary)

        # 結果を保存
        self._fit_results = {
            "X_train": X_sites,
            "Y_train": Y_sites,
            "summary": summary,
            "rhat_max": float(rhat.max()),
        }

        print("=== ベイズ空間多項回帰学習完了 ===")

        return self._fit_results

    def predict_probabilities(self, X_new: np.ndarray) -> np.ndarray:
        """事後予測"""
        if self.trace is None:
            raise ValueError(
                "モデルが学習されていません。fit()を先に実行してください。"
            )

        print(f"事後予測実行: {len(X_new)}地点")

        X_scaled = (X_new - self.X_mean) / self.X_std

        beta_samples = self.trace.posterior["beta"].values
        beta_flat = beta_samples.reshape(
            -1, beta_samples.shape[-2], beta_samples.shape[-1]
        )

        # サンプル数を制限
        n_samples = min(200, len(beta_flat))
        idx = np.random.choice(len(beta_flat), n_samples, replace=False)
        beta_selected = beta_flat[idx]

        predictions = []
        for beta_sample in beta_selected:
            eta = X_scaled @ beta_sample
            eta_padded = np.column_stack([eta, np.zeros(len(X_new))])
            p = self._softmax(eta_padded)
            predictions.append(p)

        mean_predictions = np.array(predictions).mean(axis=0)

        return mean_predictions

    def _softmax(self, x: np.ndarray) -> np.ndarray:
        """数値安定なソフトマックス"""
        x_clipped = np.clip(x, -50, 50)
        x_max = np.max(x_clipped, axis=1, keepdims=True)
        exp_x = np.exp(x_clipped - x_max)
        return exp_x / np.sum(exp_x, axis=1, keepdims=True)

    def predict_site_ratios(
        self, preprocessor: ObsidianDataPreprocessor
    ) -> Dict[str, pl.DataFrame]:
        """遺跡の産地構成比を予測"""
        if self._fit_results is None:
            raise ValueError("モデルが学習されていません。")

        # 学習データと同じ形式で予測用データを作成
        X_train = self._fit_results["X_train"]
        predictions = self.predict_probabilities(X_train)

        # 学習に使用した遺跡IDを取得
        # データがある遺跡のみを対象としている
        time_periods = {0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"}
        origins = ["神津島", "信州", "箱根", "高原山", "その他"]
        origins_filtered = [org for org in origins if org != "その他"]

        data = preprocessor.create_X_Y(self.variable_names, time_periods, origins)
        Y_all = data["Y"]
        Y_period = Y_all[:, 0 : len(origins_filtered)]  # period=0

        has_data_mask = Y_period.sum(axis=1) > 0
        site_indices = np.where(has_data_mask)[0]

        # 遺跡IDを取得
        all_site_ids = preprocessor.df_sites["遺跡ID"].sort().to_list()
        training_site_ids = [all_site_ids[i] for i in site_indices]

        # 結果をDataFrameに変換
        ratio_sites_df = pl.DataFrame({"遺跡ID": training_site_ids})

        for i, origin in enumerate(origins_filtered):
            ratio_sites_df = ratio_sites_df.with_columns(
                pl.Series(f"比率_0_{origin}", predictions[:, i])
            )

        return {"0": ratio_sites_df}

    def predict_grid_ratios(
        self, preprocessor: ObsidianDataPreprocessor
    ) -> pl.DataFrame:
        """グリッドの産地構成比を予測"""
        if self._fit_results is None:
            raise ValueError("モデルが学習されていません。")

        # グリッドデータの準備
        df_viz = preprocessor.df_elevation.clone()
        land_mask = (~df_viz["is_sea"]) & (~df_viz["average_elevation"].is_null())
        df_land = df_viz.filter(land_mask)

        # 大きすぎる場合はサンプリング
        if len(df_land) > 30000:
            df_land = df_land.sample(n=30000, seed=42)
            print(f"グリッドデータをサンプリング: {len(df_land)}点")

        # グリッド用特徴量
        grid_features = [np.ones(len(df_land))]  # インターセプト

        for var_name in self.variable_names:
            values = df_land[var_name].to_numpy()
            values = np.where(np.isnan(values), np.nanmean(values), values)
            grid_features.append(values)

        X_grid = np.column_stack(grid_features)

        # 予測実行
        predictions = self.predict_probabilities(X_grid)

        # 結果をDataFrameに変換
        ratio_df = pl.DataFrame({"x": df_land["x"], "y": df_land["y"]})

        origins_filtered = ["神津島", "信州", "箱根", "高原山"]
        for i, origin in enumerate(origins_filtered):
            ratio_df = ratio_df.with_columns(
                pl.Series(f"ratio_0_{origin}", predictions[:, i])
            )

        return ratio_df

    @property
    def results(self) -> Dict[str, Any]:
        """学習結果"""
        if self._fit_results is None:
            raise ValueError("モデルが学習されていません。")
        return self._fit_results
