"""
Fixed Bayesian Spatial Multinomial モデルのmodel3パイプライン統合版
"""

from typing import Any, Dict, List, Optional, Tuple

import arviz as az
import numpy as np
import polars as pl
import pymc as pm
from tqdm import tqdm

from ..base.base_model import BaseCompositionModel
from ..preprocessing.data_preprocessor import ObsidianDataPreprocessor


class BayesianSpatialMultinomialModel(BaseCompositionModel):
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
        self._fit_results: Optional[Dict[int, Dict[str, Any]]] = None

    def _prepare_data_for_period(
        self, preprocessor: ObsidianDataPreprocessor, target_period: int = 0
    ) -> Tuple[np.ndarray, np.ndarray, list]:
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
        site_indices = np.where(has_data_mask)[0]
        X_sites = X_spatial[has_data_mask]
        Y_sites = Y_period[has_data_mask]

        # 遺跡IDを取得
        all_site_ids = preprocessor.df_sites["遺跡ID"].sort().to_list()
        training_site_ids = [all_site_ids[i] for i in site_indices]

        # インターセプト追加
        X_sites_with_intercept = np.column_stack([np.ones(len(X_sites)), X_sites])

        print(f"学習用データ形状: X={X_sites_with_intercept.shape}, Y={Y_sites.shape}")
        print(f"産地別合計: {[Y_sites[:, i].sum() for i in range(Y_sites.shape[1])]}")

        return X_sites_with_intercept, Y_sites, training_site_ids

    def _fit_single_period(
        self, preprocessor: ObsidianDataPreprocessor, target_period: int
    ) -> Dict[str, Any]:
        """単一時期のモデル学習"""
        print(f"=== 時期{target_period}の学習開始 ===")

        # データ準備
        x_sites, y_sites, site_ids = self._prepare_data_for_period(
            preprocessor, target_period
        )

        # データがない場合はスキップ
        if len(x_sites) == 0:
            print(f"時期{target_period}のデータがありません")
            return {}

        # モデル構築（各時期で独立）
        n_sites, n_vars = x_sites.shape
        n_origins = y_sites.shape[1]

        # データ標準化（時期別）
        x_mean = x_sites.mean(axis=0)
        x_std = x_sites.std(axis=0)
        x_std[x_std == 0] = 1
        x_scaled = (x_sites - x_mean) / x_std

        with pm.Model() as model:
            # 係数の事前分布
            beta = pm.Normal(
                "beta", mu=0, sigma=self.prior_sigma, shape=(n_vars, n_origins - 1)
            )

            # 線形予測子
            eta = pm.math.dot(x_scaled, beta)
            eta_padded = pm.math.concatenate([eta, pm.math.zeros((n_sites, 1))], axis=1)

            # 多項分布
            p = pm.math.softmax(eta_padded, axis=1)
            n_total = pm.math.sum(y_sites, axis=1)
            _ = pm.Multinomial("likelihood", n=n_total, p=p, observed=y_sites)

        # MCMC実行
        with model:
            trace = pm.sample(
                draws=self.n_draws,
                tune=self.n_tune,
                return_inferencedata=True,
                random_seed=self.random_seed,
                progressbar=True,
            )

        # 結果保存
        summary = az.summary(trace, var_names=["beta"])
        result = {
            "X_train": x_sites,
            "Y_train": y_sites,
            "site_ids": site_ids,
            "X_mean": x_mean,
            "X_std": x_std,
            "trace": trace,
            "summary": summary,
            "rhat_max": float(summary["r_hat"].max()),
        }

        print(f"時期{target_period}: 学習完了 (Rhat最大={result['rhat_max']:.3f})")
        return result

    def fit(self, preprocessor: ObsidianDataPreprocessor) -> Dict[str, Any]:
        """全時期のモデル学習"""
        print("=== ベイズ空間多項回帰学習開始（全時期） ===")

        time_periods = {0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"}
        self._fit_results = {}

        for period in time_periods.keys():
            result = self._fit_single_period(preprocessor, period)
            if result:  # データがある時期のみ保存
                self._fit_results[period] = result

        print(f"=== 学習完了: {len(self._fit_results)}時期のモデル ===")
        return self._fit_results

    def _predict_probabilities_single_period(
        self, x_new: np.ndarray, period_result: Dict[str, Any]
    ) -> np.ndarray:
        """単一時期の事後予測"""
        x_scaled = (x_new - period_result["X_mean"]) / period_result["X_std"]

        beta_samples = period_result["trace"].posterior["beta"].values
        beta_flat = beta_samples.reshape(
            -1, beta_samples.shape[-2], beta_samples.shape[-1]
        )

        # サンプル数を制限
        n_samples = min(200, len(beta_flat))
        idx = np.random.choice(len(beta_flat), n_samples, replace=False)
        beta_selected = beta_flat[idx]

        predictions = []
        for beta_sample in beta_selected:
            eta = x_scaled @ beta_sample
            eta_padded = np.column_stack([eta, np.zeros(len(x_new))])
            p = self._softmax(eta_padded)
            predictions.append(p)

        return np.array(predictions).mean(axis=0)

    def _softmax(self, x: np.ndarray) -> np.ndarray:
        """数値安定なソフトマックス"""
        x_clipped = np.clip(x, -50, 50)
        x_max = np.max(x_clipped, axis=1, keepdims=True)
        exp_x = np.exp(x_clipped - x_max)
        return exp_x / np.sum(exp_x, axis=1, keepdims=True)

    def predict_site_ratios(
        self, preprocessor: ObsidianDataPreprocessor
    ) -> Dict[str, pl.DataFrame]:
        """遺跡の産地構成比を予測（BaseCompositionModel準拠・標準化済み）"""
        if self._fit_results is None:
            raise ValueError("モデルが学習されていません。")

        # 新しいpredict_all_periods_originsを使用
        time_periods = {0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"}
        origins = ["神津島", "信州", "箱根", "高原山", "その他"]
        _, ratio_sites_dict = self.predict_all_periods_origins(
            preprocessor, time_periods, origins
        )

        # BayesianNadarayaWatsonと同じ形式に変換
        result_dict = {}
        if isinstance(ratio_sites_dict, dict):
            for period_idx, period_df in tqdm(ratio_sites_dict.items()):
                # 主要4産地のカラムのみ選択
                select_columns = ["遺跡ID"]
                for origin in ["神津島", "信州", "箱根", "高原山"]:
                    col_name = f"ratio_{period_idx}_{origin}"
                    if col_name in period_df.columns:
                        select_columns.append(col_name)

                if len(select_columns) > 1:  # 遺跡ID以外のカラムがある場合
                    filtered_df = period_df.select(select_columns)
                    result_dict[str(period_idx)] = filtered_df

        return result_dict

    def predict_grid_ratios(
        self, preprocessor: ObsidianDataPreprocessor
    ) -> pl.DataFrame:
        """グリッドの産地構成比を予測"""
        # 新しいpredict_all_periods_originsを使用
        time_periods = {0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"}
        origins = ["神津島", "信州", "箱根", "高原山", "その他"]
        ratio_df, _ = self.predict_all_periods_origins(
            preprocessor, time_periods, origins
        )
        return ratio_df

    @property
    def model_name(self) -> str:
        """モデル名"""
        return "BayesianSpatialRegression"

    def _predict_single_period(
        self, preprocessor: ObsidianDataPreprocessor, target_period: int
    ) -> Tuple[pl.DataFrame, pl.DataFrame]:
        """単一時期の予測（グリッド・遺跡両方）"""
        if target_period not in self._fit_results:
            return pl.DataFrame(), pl.DataFrame()

        period_result = self._fit_results[target_period]

        # グリッド予測の準備
        df_viz = preprocessor.df_elevation.clone()
        land_mask = (~df_viz["is_sea"]) & (~df_viz["average_elevation"].is_null())
        df_land = df_viz.filter(land_mask)

        # グリッド用特徴量
        grid_features = [np.ones(len(df_land))]  # インターセプト
        for var_name in self.variable_names:
            values = df_land[var_name].to_numpy()
            values = np.where(np.isnan(values), np.nanmean(values), values)
            grid_features.append(values)

        x_grid = np.column_stack(grid_features)
        grid_predictions = self._predict_probabilities_single_period(
            x_grid, period_result
        )

        # グリッド結果
        origins_filtered = ["神津島", "信州", "箱根", "高原山"]
        grid_df = pl.DataFrame({"x": df_land["x"], "y": df_land["y"]})
        for i, origin in enumerate(origins_filtered):
            grid_df = grid_df.with_columns(
                pl.Series(f"ratio_{target_period}_{origin}", grid_predictions[:, i])
            )

        # 遺跡予測
        x_sites = period_result["X_train"]
        site_predictions = self._predict_probabilities_single_period(
            x_sites, period_result
        )

        # 遺跡結果
        sites_df = pl.DataFrame({"遺跡ID": period_result["site_ids"]})
        for i, origin in enumerate(origins_filtered):
            sites_df = sites_df.with_columns(
                pl.Series(f"ratio_{target_period}_{origin}", site_predictions[:, i])
            )

        return grid_df, sites_df

    def predict_all_periods_origins(
        self,
        preprocessor: ObsidianDataPreprocessor,
        time_period_names: Dict[int, str],
        origin_order: list,
    ) -> Tuple[pl.DataFrame, Dict[int, pl.DataFrame]]:
        """全時期・全産地について予測"""
        if self._fit_results is None:
            raise ValueError("モデルが学習されていません。")

        # origin_orderは互換性のために受け取るが、固定の産地順序を使用
        _ = origin_order  # 明示的に使用しないことを示す

        # 基本グリッドDataFrame
        df_viz = preprocessor.df_elevation.clone()
        land_mask = (~df_viz["is_sea"]) & (~df_viz["average_elevation"].is_null())
        df_land = df_viz.filter(land_mask)
        ratio_df = pl.DataFrame({"x": df_land["x"], "y": df_land["y"]})

        ratio_sites_dict = {}

        for period in time_period_names.keys():
            if period in self._fit_results:
                grid_result, sites_result = self._predict_single_period(
                    preprocessor, period
                )

                # グリッド結果をマージ
                for col in grid_result.columns:
                    if col.startswith("ratio_"):
                        ratio_df = ratio_df.with_columns(grid_result[col])

                # 遺跡結果を保存
                ratio_sites_dict[period] = sites_result

        return ratio_df, ratio_sites_dict

    @property
    def results(self) -> Dict[str, Any]:
        """学習結果"""
        if self._fit_results is None:
            raise ValueError("モデルが学習されていません。")
        return self._fit_results
