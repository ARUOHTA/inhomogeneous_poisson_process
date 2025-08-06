"""
モデル比較フレームワーク

複数モデルの学習・評価・比較を統一的に実行
DRY原則に基づく重複コード排除
"""

from typing import Dict, List

import matplotlib.pyplot as plt
import polars as pl

from ..base.base_model import BaseCompositionModel
from ..preprocessing.data_preprocessor import ObsidianDataPreprocessor
from .unified_loocv import LOOCVConfig, UnifiedLOOCVEvaluator


class ModelComparison:
    """複数モデルの比較実験フレームワーク"""

    def __init__(
        self,
        models: List[BaseCompositionModel],
        loocv_config: LOOCVConfig = None,
        variable_names: List[str] = None,
    ):
        """
        Parameters
        ----------
        models : List[BaseCompositionModel]
            比較対象のモデルリスト
        loocv_config : LOOCVConfig, optional
            LOOCV設定
        variable_names : List[str], optional
            使用する説明変数名
        """
        self.models = models
        self.loocv_config = loocv_config or LOOCVConfig()
        self.variable_names = variable_names or ["average_elevation"]
        self.evaluator = UnifiedLOOCVEvaluator(self.loocv_config, self.variable_names)
        self.results = {}

    def run_comparison(
        self, preprocessor: ObsidianDataPreprocessor
    ) -> "ComparisonResults":
        """
        全モデルの学習・評価・比較を実行

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ

        Returns
        -------
        ComparisonResults
            比較結果
        """
        print("=== モデル比較実験開始 ===")

        # 1. 各モデルの学習
        print("\n1. モデル学習フェーズ")
        trained_models = []
        for model in self.models:
            print(f"  学習中: {model.model_name}")
            model.fit(preprocessor)
            trained_models.append(model)

        # 2. LOOCV評価
        print("\n2. LOOCV評価フェーズ")
        loocv_results = self.evaluator.evaluate_all_models(trained_models, preprocessor)

        # 3. 予測結果の比較
        print("\n3. 予測結果比較フェーズ")
        prediction_results = self._compare_predictions(trained_models, preprocessor)

        # 4. 結果の統合
        comparison_results = ComparisonResults(
            loocv_results=loocv_results,
            prediction_results=prediction_results,
            model_names=[model.model_name for model in trained_models],
        )

        print("=== モデル比較実験完了 ===")
        return comparison_results

    def _compare_predictions(
        self, models: List[BaseCompositionModel], preprocessor: ObsidianDataPreprocessor
    ) -> Dict[str, Dict[str, pl.DataFrame]]:
        """各モデルの予測結果を比較"""
        prediction_results = {}

        for model in models:
            prediction_results[model.model_name] = {
                "site_ratios": model.predict_site_ratios(preprocessor),
                "grid_ratios": model.predict_grid_ratios(preprocessor),
            }

        return prediction_results


class ComparisonResults:
    """モデル比較結果を格納するクラス"""

    def __init__(
        self,
        loocv_results: Dict[str, pl.DataFrame],
        prediction_results: Dict[str, Dict[str, pl.DataFrame]],
        model_names: List[str],
    ):
        self.loocv_results = loocv_results
        self.prediction_results = prediction_results
        self.model_names = model_names

    def get_loocv_summary(self) -> pl.DataFrame:
        """LOOCV結果の要約統計を取得"""
        summary_data = []

        for model_name, results_df in self.loocv_results.items():
            summary = {
                "モデル": model_name,
                "Aitchison距離_平均": results_df["Aitchison距離"].mean(),
                "Aitchison距離_標準偏差": results_df["Aitchison距離"].std(),
                "Bray-Curtis_平均": results_df["Bray-Curtis"].mean(),
                "Bray-Curtis_標準偏差": results_df["Bray-Curtis"].std(),
                "Jensen-Shannon_平均": results_df["Jensen-Shannon"].mean(),
                "Jensen-Shannon_標準偏差": results_df["Jensen-Shannon"].std(),
                "評価件数": len(results_df),
            }
            summary_data.append(summary)

        return pl.DataFrame(summary_data)

    def plot_loocv_comparison(self) -> plt.Figure:
        """LOOCV結果の比較プロット"""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        metrics = ["Aitchison距離", "Bray-Curtis", "Jensen-Shannon", "Total_Variation"]

        for i, metric in enumerate(metrics):
            ax = axes[i // 2, i % 2]

            for model_name, results_df in self.loocv_results.items():
                ax.hist(results_df[metric], alpha=0.6, label=model_name, bins=20)

            ax.set_xlabel(metric)
            ax.set_ylabel("頻度")
            ax.set_title(f"{metric}の分布比較")
            ax.legend()

        plt.tight_layout()
        return fig

    def get_best_model(self, metric: str = "Aitchison距離") -> str:
        """指定された評価指標で最良のモデルを取得"""
        summary = self.get_loocv_summary()
        metric_col = f"{metric}_平均"

        if metric_col in summary.columns:
            best_idx = summary[metric_col].arg_min()
            return summary[best_idx, "モデル"]
        else:
            raise ValueError(f"指定された評価指標 '{metric}' が見つかりません")
