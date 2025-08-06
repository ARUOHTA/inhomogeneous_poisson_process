"""
Model3リファクタリング - モデル統合パッケージ

DRY原則に基づく重複コード排除と共通インターフェース化
"""

# 設定とパイプライン
# 基底クラス
from .base import BaseCompositionModel, BaseIntensityModel

# 構成比モデル
from .composition import (
    BayesianNadarayaWatson,
    BayesianSpatialMultinomialModel,
    KSBPModel,
    NadarayaWatsonEstimator,
)
from .config import Model3Config, Model3Pipeline

# 評価モジュール
from .evaluation import (
    ComparisonResults,
    CompositionalMetrics,
    LOOCVConfig,
    LOOCVEvaluator,
    ModelComparison,
    UnifiedLOOCVEvaluator,
    bray_curtis_dissimilarity,
    jensen_shannon_divergence,
)

# 強度モデル
from .intensity import InhomogeneousPoissonProcess, IntensityFunction

# データ前処理
from .preprocessing import ObsidianDataPreprocessor

# 可視化
from .visualization import BaseVisualizer, ObsidianVisualizer

__all__ = [
    # 設定とパイプライン
    "Model3Config",
    "Model3Pipeline",
    # 基底クラス
    "BaseCompositionModel",
    "BaseIntensityModel",
    # 構成比モデル
    "NadarayaWatsonEstimator",
    "BayesianNadarayaWatson",
    "KSBPModel",
    "BayesianSpatialMultinomialModel",
    # 強度モデル
    "IntensityFunction",
    "InhomogeneousPoissonProcess",
    # データ前処理
    "ObsidianDataPreprocessor",
    # 評価モジュール
    "UnifiedLOOCVEvaluator",
    "ModelComparison",
    "ComparisonResults",
    "CompositionalMetrics",
    "LOOCVEvaluator",
    "LOOCVConfig",
    "bray_curtis_dissimilarity",
    "jensen_shannon_divergence",
    # 可視化
    "BaseVisualizer",
    "ObsidianVisualizer",
]
