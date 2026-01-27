"""
Leave-One-Out Cross Validation (LOOCV) システム
Nadaraya-Watson推定量の評価用
"""

import random
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl
from tqdm import tqdm

from ..base.base_model import BaseCompositionModel
from ..composition.nadaraya_watson import NadarayaWatsonEstimator
from ..config.model_config import Model3Config
from ..preprocessing.data_preprocessor import ObsidianDataPreprocessor
from .metrics import CompositionalMetrics


@dataclass
class LOOCVConfig:
    """LOOCV実行のための設定"""

    n_trials: int = 100
    zero_replacement: float = 1e-6
    random_seed: int = 42
    save_intermediate: bool = True
    verbose: bool = True


class LOOCVEvaluator:
    """Leave-One-Out Cross Validation評価器"""

    def __init__(
        self,
        preprocessor: ObsidianDataPreprocessor,
        model_config: Model3Config,
        loocv_config: LOOCVConfig,
    ):
        """
        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みデータ
        model_config : Model3Config
            モデル設定
        loocv_config : LOOCVConfig
            LOOCV設定
        """
        self.preprocessor = preprocessor
        self.model_config = model_config
        self.loocv_config = loocv_config
        self.metrics = CompositionalMetrics(
            zero_replacement=loocv_config.zero_replacement
        )

        # 乱数シードを設定
        if loocv_config.random_seed is not None:
            random.seed(loocv_config.random_seed)
            np.random.seed(loocv_config.random_seed)

        # 遺跡IDのリストを取得
        self.site_ids = self.preprocessor.df_sites["遺跡ID"].unique().sort().to_list()

        # 各時期・産地での観測構成比を事前に計算
        self._precompute_observed_compositions()

    def _precompute_observed_compositions(self):
        """各遺跡・時期での観測構成比を事前計算"""
        self.observed_compositions = {}

        for period in self.model_config.time_periods.keys():
            self.observed_compositions[period] = {}

            # 対象時期のデータを抽出
            period_data = self.preprocessor.df_obsidian.filter(pl.col("時期") == period)

            for site_id in self.site_ids:
                # 各遺跡での産地別出土数を集計
                site_data = period_data.filter(pl.col("遺跡ID") == site_id)

                if site_data.height == 0:
                    # データがない場合は均等分布で初期化
                    n_origins = len(self.model_config.origins)
                    composition = np.ones(n_origins) / n_origins
                else:
                    composition = np.zeros(len(self.model_config.origins))
                    total_count = site_data.height

                    for i, origin in enumerate(self.model_config.origins):
                        if origin == "その他":
                            # "その他"は残りの全て
                            other_count = site_data.filter(
                                ~pl.col("産地カテゴリ").is_in(
                                    self.model_config.origins[:-1]
                                )
                            ).height
                            composition[i] = (
                                other_count / total_count if total_count > 0 else 0
                            )
                        else:
                            count = site_data.filter(
                                pl.col("産地カテゴリ") == origin
                            ).height
                            composition[i] = (
                                count / total_count if total_count > 0 else 0
                            )

                    # 構成比の合計が1になるように正規化
                    if np.sum(composition) > 0:
                        composition = composition / np.sum(composition)
                    else:
                        # 全て0の場合は均等分布
                        composition = np.ones(len(self.model_config.origins)) / len(
                            self.model_config.origins
                        )

                self.observed_compositions[period][site_id] = composition

    def _create_excluded_dataset(
        self, excluded_site_ids: List[int]
    ) -> ObsidianDataPreprocessor:
        """
        指定された遺跡を除外したデータセットを作成

        Parameters
        ----------
        excluded_site_ids : List[int]
            除外する遺跡IDのリスト

        Returns
        -------
        ObsidianDataPreprocessor
            除外後のデータセット
        """
        # 元のデータをコピー
        excluded_preprocessor = ObsidianDataPreprocessor(self.preprocessor.data_dir)

        # データを読み込み
        data = excluded_preprocessor.load_data()

        # 遺跡データから指定IDを除外
        excluded_obsidian = excluded_preprocessor.df_obsidian.filter(
            ~pl.col("遺跡ID").is_in(excluded_site_ids)
        )
        excluded_sites = excluded_preprocessor.df_sites.filter(
            ~pl.col("遺跡ID").is_in(excluded_site_ids)
        )

        # データを更新
        excluded_preprocessor._df_obsidian = excluded_obsidian
        excluded_preprocessor._df_sites = excluded_sites

        return excluded_preprocessor

    def _predict_composition_at_coordinates(
        self,
        nw_estimator: NadarayaWatsonEstimator,
        preprocessor: ObsidianDataPreprocessor,
        target_lon: float,
        target_lat: float,
        period: int,
    ) -> np.ndarray:
        """
        指定された座標での産地構成比を予測

        座標に基づく直接的な予測を行い、配列サイズの問題を回避する

        Parameters
        ----------
        nw_estimator : NadarayaWatsonEstimator
            学習済みのNW推定量
        preprocessor : ObsidianDataPreprocessor
            前処理済みデータ
        target_lon : float
            対象経度
        target_lat : float
            対象緯度
        period : int
            対象時期

        Returns
        -------
        np.ndarray
            予測された産地構成比ベクトル
        """
        predicted_composition = np.zeros(len(self.model_config.origins))

        # メッシュグリッドの取得
        lon_mesh, lat_mesh = preprocessor.create_meshgrid()

        # 対象座標に最も近いグリッド点を探す
        distances = np.sqrt((lon_mesh - target_lon) ** 2 + (lat_mesh - target_lat) ** 2)
        min_idx = np.unravel_index(np.argmin(distances), distances.shape)
        grid_flat_idx = min_idx[0] * lon_mesh.shape[1] + min_idx[1]

        # 各産地について個別に予測
        for i, origin in enumerate(self.model_config.origins):
            if origin == "その他":
                # "その他"は最後に計算
                continue

            try:
                # 直接的にNW推定を行う
                predicted_ratio = self._predict_single_origin_at_grid(
                    nw_estimator, preprocessor, period, origin, grid_flat_idx
                )
                predicted_composition[i] = predicted_ratio

            except Exception as e:
                if self.loocv_config.verbose:
                    print(
                        f"警告: 座標({target_lon:.3f}, {target_lat:.3f}), 産地{origin}の予測でエラー: {e}"
                    )
                predicted_composition[i] = 0.0

        # "その他"の比率を計算
        other_ratio = 1.0 - np.sum(predicted_composition[:-1])
        predicted_composition[-1] = max(0.0, other_ratio)

        # 構成比の正規化
        total = np.sum(predicted_composition)
        if total > 0:
            predicted_composition = predicted_composition / total
        else:
            # 全て0の場合は均等分布
            predicted_composition = np.ones(len(self.model_config.origins)) / len(
                self.model_config.origins
            )

        return predicted_composition

    def _predict_single_origin_at_grid(
        self,
        nw_estimator: NadarayaWatsonEstimator,
        preprocessor: ObsidianDataPreprocessor,
        period: int,
        origin: str,
        grid_idx: int,
    ) -> float:
        """
        特定のグリッド点での単一産地の比率を予測

        Parameters
        ----------
        nw_estimator : NadarayaWatsonEstimator
            学習済みのNW推定量
        preprocessor : ObsidianDataPreprocessor
            前処理済みデータ
        period : int
            対象時期
        origin : str
            対象産地
        grid_idx : int
            グリッドのフラットインデックス

        Returns
        -------
        float
            予測された比率
        """
        # 除外後のデータセットでのpreprocess_obsidian_dataを実行
        counts, target_counts = preprocessor.preprocess_obsidian_data(period, origin)

        # NW推定量の重み行列を取得
        weights = nw_estimator.weights  # (グリッド数, 残存遺跡数)

        # 指定グリッド点での重み付き比率を計算
        grid_weights = weights[grid_idx, :]  # (残存遺跡数,)

        # 残存遺跡のIDを取得
        remaining_site_ids = preprocessor.df_sites["遺跡ID"].sort().to_list()

        # countsとtarget_countsから必要な部分だけを抽出
        valid_counts = np.array(
            [
                counts[site_id] if site_id < len(counts) else 0
                for site_id in remaining_site_ids
            ]
        )
        valid_target_counts = np.array(
            [
                target_counts[site_id] if site_id < len(target_counts) else 0
                for site_id in remaining_site_ids
            ]
        )

        # 重み付き比率の計算
        weighted_total = np.sum(grid_weights * valid_counts)
        weighted_target = np.sum(grid_weights * valid_target_counts)

        if weighted_total > 0:
            return weighted_target / weighted_total
        else:
            return 0.0

    def run_single_trial(self, period: int, test_site_id: int) -> dict:
        """
        単一のLOOCV試行を実行

        Parameters
        ----------
        period : int
            対象時期
        test_site_id : int
            テスト対象の遺跡ID

        Returns
        -------
        dict
            試行結果
        """
        try:
            # 1. テスト遺跡を除外したデータセットを作成
            excluded_preprocessor = self._create_excluded_dataset([test_site_id])

            # 2. NW推定量を再学習
            nw_estimator = NadarayaWatsonEstimator(
                sigma=self.model_config.nw_sigma,
                sigma_for_sites=self.model_config.nw_sigma_for_sites,
                variable_names=self.model_config.nw_variable_names,
            )
            nw_estimator.fit(excluded_preprocessor)

            # 3. テスト遺跡での構成比を予測
            # まず対象遺跡の座標を取得
            target_site_info = self.preprocessor.df_sites.filter(
                pl.col("遺跡ID") == test_site_id
            )

            if target_site_info.height == 0:
                raise ValueError(f"遺跡ID {test_site_id} が見つかりません")

            target_lon = target_site_info["経度"].item()
            target_lat = target_site_info["緯度"].item()

            predicted_composition = self._predict_composition_at_coordinates(
                nw_estimator, excluded_preprocessor, target_lon, target_lat, period
            )

            # 4. 観測値を取得
            observed_composition = self.observed_compositions[period][test_site_id]

            # 5. 評価指標を計算
            metrics = self.metrics.compute_all_metrics(
                observed_composition, predicted_composition
            )

            return {
                "success": True,
                "period": period,
                "test_site_id": test_site_id,
                "observed_composition": observed_composition,
                "predicted_composition": predicted_composition,
                "aitchison_distance": metrics["aitchison_distance"],
                "total_variation": metrics["total_variation"],
                "observed_sum": metrics["observed_sum"],
                "predicted_sum": metrics["predicted_sum"],
                "error": None,
            }

        except Exception as e:
            if self.loocv_config.verbose:
                print(f"エラー: 遺跡{test_site_id}, 時期{period}の試行でエラー: {e}")

            return {
                "success": False,
                "period": period,
                "test_site_id": test_site_id,
                "observed_composition": None,
                "predicted_composition": None,
                "aitchison_distance": np.nan,
                "total_variation": np.nan,
                "observed_sum": np.nan,
                "predicted_sum": np.nan,
                "error": str(e),
            }

    def run_period_evaluation(
        self, period: int, n_trials: Optional[int] = None
    ) -> dict:
        """
        特定時期での全体評価を実行

        Parameters
        ----------
        period : int
            対象時期
        n_trials : int, optional
            試行回数（Noneの場合は設定値を使用）

        Returns
        -------
        dict
            評価結果
        """
        if n_trials is None:
            n_trials = self.loocv_config.n_trials

        # 該当時期にデータがある遺跡のみを対象
        available_sites = []
        for site_id in self.site_ids:
            if site_id in self.observed_compositions[period]:
                # 全て均等分布でない場合のみ対象とする
                composition = self.observed_compositions[period][site_id]
                if not np.allclose(composition, 1.0 / len(composition)):
                    available_sites.append(site_id)

        if len(available_sites) == 0:
            return {
                "period": period,
                "n_trials": 0,
                "trial_results": [],
                "summary_statistics": {},
                "error": f"時期{period}にデータのある遺跡がありません",
            }

        # 試行回数を調整（利用可能な遺跡数が少ない場合）
        actual_trials = min(n_trials, len(available_sites) * 3)  # 最大でも各遺跡3回まで

        trial_results = []

        desc = f"Period {period} LOOCV"
        for trial in tqdm(
            range(actual_trials), desc=desc, disable=not self.loocv_config.verbose
        ):
            # ランダムに遺跡を選択
            test_site_id = random.choice(available_sites)

            # 試行を実行
            result = self.run_single_trial(period, test_site_id)
            result["trial_id"] = trial
            trial_results.append(result)

        # 統計量を計算
        summary_stats = self._compute_summary_statistics(trial_results)

        return {
            "period": period,
            "n_trials": actual_trials,
            "available_sites": available_sites,
            "trial_results": trial_results,
            "summary_statistics": summary_stats,
        }

    def run_all_periods_evaluation(self, n_trials: Optional[int] = None) -> dict:
        """
        全時期での評価を実行

        Parameters
        ----------
        n_trials : int, optional
            各時期での試行回数

        Returns
        -------
        dict
            全時期の評価結果
        """
        all_results = {}

        for period in self.model_config.time_periods.keys():
            if self.loocv_config.verbose:
                period_name = self.model_config.time_periods[period]
                print(f"\n=== 時期 {period} ({period_name}) の評価を開始 ===")

            period_result = self.run_period_evaluation(period, n_trials)
            all_results[period] = period_result

        return all_results

    def _compute_summary_statistics(self, trial_results: List[dict]) -> dict:
        """
        試行結果から統計量を計算

        Parameters
        ----------
        trial_results : List[dict]
            試行結果のリスト

        Returns
        -------
        dict
            統計量
        """
        successful_trials = [r for r in trial_results if r["success"]]

        if len(successful_trials) == 0:
            return {
                "n_successful_trials": 0,
                "success_rate": 0.0,
                "mean_aitchison_distance": np.nan,
                "std_aitchison_distance": np.nan,
                "mean_total_variation": np.nan,
                "std_total_variation": np.nan,
            }

        aitchison_distances = [
            r["aitchison_distance"]
            for r in successful_trials
            if not np.isnan(r["aitchison_distance"])
        ]
        total_variations = [
            r["total_variation"]
            for r in successful_trials
            if not np.isnan(r["total_variation"])
        ]

        return {
            "n_successful_trials": len(successful_trials),
            "success_rate": len(successful_trials) / len(trial_results)
            if len(trial_results) > 0
            else 0.0,
            "mean_aitchison_distance": np.mean(aitchison_distances)
            if len(aitchison_distances) > 0
            else np.nan,
            "std_aitchison_distance": np.std(aitchison_distances)
            if len(aitchison_distances) > 0
            else np.nan,
            "median_aitchison_distance": np.median(aitchison_distances)
            if len(aitchison_distances) > 0
            else np.nan,
            "mean_total_variation": np.mean(total_variations)
            if len(total_variations) > 0
            else np.nan,
            "std_total_variation": np.std(total_variations)
            if len(total_variations) > 0
            else np.nan,
            "median_total_variation": np.median(total_variations)
            if len(total_variations) > 0
            else np.nan,
            "aitchison_distances": aitchison_distances,
            "total_variations": total_variations,
        }
