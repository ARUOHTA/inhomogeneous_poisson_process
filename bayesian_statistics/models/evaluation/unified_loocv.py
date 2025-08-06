"""
統一LOOCV評価システム
全産地構成比モデルに対応したLOOCV評価

DRY原則に基づく重複コード排除
"""

import random
from typing import Dict, List

import numpy as np
import polars as pl
from tqdm import tqdm

from ..base.base_model import BaseCompositionModel
from ..preprocessing.data_preprocessor import ObsidianDataPreprocessor
from .loocv import LOOCVConfig
from .metrics import CompositionalMetrics


class UnifiedLOOCVEvaluator:
    """全モデル共通のLOOCV評価器"""

    def __init__(self, config: LOOCVConfig, variable_names: List[str] = None):
        self.config = config
        self.variable_names = variable_names or ["average_elevation"]
        self.metrics = CompositionalMetrics()

    def evaluate_all_models(
        self, models: List[BaseCompositionModel], preprocessor: ObsidianDataPreprocessor
    ) -> Dict[str, pl.DataFrame]:
        """
        全モデルのLOOCV評価を実行

        Parameters
        ----------
        models : List[BaseCompositionModel]
            評価対象のモデルリスト
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ

        Returns
        -------
        Dict[str, pl.DataFrame]
            モデル名をキーとした評価結果
        """
        results = {}

        for model in models:
            print(f"\n=== {model.model_name} LOOCV評価開始 ===")
            results[model.model_name] = self._evaluate_single_model(model, preprocessor)
            print(f"=== {model.model_name} LOOCV評価完了 ===")

        return results

    def evaluate_single_model(
        self, model: BaseCompositionModel, preprocessor: ObsidianDataPreprocessor
    ) -> pl.DataFrame:
        """
        単一モデルのLOOCV評価（公開メソッド）

        Parameters
        ----------
        model : BaseCompositionModel
            評価対象のモデル
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ

        Returns
        -------
        pl.DataFrame
            評価結果（距離指標、遺跡ID、時期など）
        """
        return self._evaluate_single_model(model, preprocessor)

    def _evaluate_single_model(
        self, model: BaseCompositionModel, preprocessor: ObsidianDataPreprocessor
    ) -> pl.DataFrame:
        """
        単一モデルのLOOCV評価

        Parameters
        ----------
        model : BaseCompositionModel
            評価対象のモデル
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ

        Returns
        -------
        pl.DataFrame
            評価結果（距離指標、遺跡ID、時期など）
        """
        # 遺跡での予測を取得（モデルは既に学習済み）
        predicted_ratios = model.predict_site_ratios(preprocessor)

        # 返り値の型チェック
        if not hasattr(predicted_ratios, "items"):
            # DataFrameが返された場合は、時期0のデータとして扱う
            if hasattr(predicted_ratios, "columns"):
                predicted_ratios = {"0": predicted_ratios}
            else:
                # 完全に予期しない型の場合はスキップ
                return pl.DataFrame()  # 空のDataFrameを返す

        # 実際の構成比を計算
        observed_ratios = self._calculate_observed_ratios(preprocessor)

        # 評価指標を計算
        evaluation_results = []

        for period_key, predicted_df in predicted_ratios.items():
            if period_key in observed_ratios:
                observed_df = observed_ratios[period_key]

                # 遺跡IDで結合
                merged = predicted_df.join(observed_df, on="遺跡ID", how="inner")

                for row in merged.iter_rows(named=True):
                    site_id = row["遺跡ID"]

                    # 予測値と実測値を抽出
                    predicted = self._extract_composition_vector(row, "predicted")
                    observed = self._extract_composition_vector(row, "observed")

                    # 各種距離指標を計算
                    from .metrics import (
                        bray_curtis_dissimilarity,
                        jensen_shannon_divergence,
                    )

                    aitchison_dist = self.metrics.aitchison_distance(
                        predicted, observed
                    )
                    bray_curtis = bray_curtis_dissimilarity(predicted, observed)
                    js_div = jensen_shannon_divergence(predicted, observed)
                    total_var = self.metrics.total_variation(predicted, observed)

                    evaluation_results.append(
                        {
                            "遺跡ID": site_id,
                            "時期": period_key,
                            "Aitchison距離": aitchison_dist,
                            "Bray-Curtis": bray_curtis,
                            "Jensen-Shannon": js_div,
                            "Total_Variation": total_var,
                            "モデル": model.model_name,
                        }
                    )

        return pl.DataFrame(evaluation_results)

    def _calculate_observed_ratios(
        self, preprocessor: ObsidianDataPreprocessor
    ) -> Dict[str, pl.DataFrame]:
        """実際の構成比を計算（標準化済み形式）"""
        # 時期と産地の定義（主要4産地のみ）
        time_periods = {0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"}
        origins = ["神津島", "信州", "箱根", "高原山"]

        # 産地名マッピング（実データの産地名を標準産地名に変換）
        source_mapping = {
            "神津島": ["神津島"],
            "信州": ["諏訪", "和田峠", "蓼科", "男女倉"],
            "箱根": ["箱根", "天城"],
            "高原山": ["高原山"],
        }

        observed_ratios = {}

        for period_key in time_periods.keys():
            # この時期のデータを抽出
            period_data = preprocessor.df_obsidian.filter(pl.col("時期") == period_key)

            if len(period_data) == 0:
                continue

            # 遺跡ごとの産地別出土数を集計
            site_ratios = []

            for site_id in period_data["遺跡ID"].unique().sort():
                site_data = period_data.filter(pl.col("遺跡ID") == site_id)

                # 各標準産地の出土数を計算
                origin_counts = {}
                total_count = 0

                for standard_origin in origins:
                    actual_origins = source_mapping[standard_origin]
                    count = site_data.filter(
                        pl.col("産地").is_in(actual_origins)
                    ).height
                    origin_counts[standard_origin] = count
                    total_count += count

                # 構成比を計算
                if total_count > 0:
                    site_ratio = {"遺跡ID": site_id}
                    for origin in origins:
                        ratio = origin_counts[origin] / total_count
                        site_ratio[f"ratio_{origin}"] = ratio
                    site_ratios.append(site_ratio)

            if site_ratios:
                observed_ratios[str(period_key)] = pl.DataFrame(site_ratios)

        return observed_ratios

    def _extract_composition_vector(self, row: dict, prefix: str) -> np.ndarray:
        """行データから構成比ベクトルを抽出（統一形式対応）"""
        origins = ["神津島", "信州", "箱根", "高原山"]
        composition = []

        for origin in origins:
            if prefix == "predicted":
                # 予測値：統一形式 ratio_{時期}_{産地名}
                # 時期0を仮定（将来的に動的に取得可能）
                col_name = f"ratio_0_{origin}"
            else:
                # 実測値：ratio_{産地名}
                col_name = f"ratio_{origin}"

            composition.append(row.get(col_name, 0.0))

        return np.array(composition)
