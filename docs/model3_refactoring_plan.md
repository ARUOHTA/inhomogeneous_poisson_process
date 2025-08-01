# Model3 リファクタリング計画書

## 1. 現状分析

### 1.1 プロジェクトの目的
黒曜石の産地構成比を予測する空間統計モデルの実装と比較

### 1.2 現在の問題点
- 複数のモデル間で共通インターフェースが統一されていない
- LOOCV評価がNadaraya-Watsonモデル専用になっている
- 各モデルの比較実験が困難
- 重複コードが存在する

### 1.3 既存モデルの種類

#### 産地構成比モデル（今回のリファクタリング対象）
1. **Nadaraya-Watson (NW)**: 古典的ノンパラメトリックカーネル回帰
2. **Bayesian NW (bayes_NW)**: ベイズ拡張版NW
3. **KSBP**: Kernel Stick-Breaking Process（ベイズノンパラメトリック）
4. **Spatial Regression**: ベイズ空間多項ロジット回帰

#### 遺跡存在確率モデル
- **IPP**: Inhomogeneous Poisson Process（今回は対象外）

## 2. 依存関係の整理

### 2.1 基底モジュール（他に依存しない）
- `model3_preprocessing.py`: データ前処理
- `model3_compositional_metrics.py`: 評価指標

### 2.2 モデル実装モジュール
```
model3_nadaraya_watson.py
    └── model3_bayes_nw.py (継承関係)

model3_ksbp.py (独立)

model3_bayesian_spatial.py (独立)
```

### 2.3 統合・評価モジュール
- `model3_config.py`: 全モデルを統合するパイプライン
- `model3_loocv.py`: LOOCV評価（現在NW専用）

## 3. リファクタリング方針

### 3.1 共通インターフェースの設計

#### BaseCompositionModel（抽象基底クラス）
```python
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional
import numpy as np

class BaseCompositionModel(ABC):
    """産地構成比モデルの共通インターフェース"""

    @abstractmethod
    def fit(self, preprocessor: ObsidianDataPreprocessor, **kwargs) -> Dict[str, Any]:
        """モデルの学習"""
        pass

    @abstractmethod
    def predict_site_ratios(self, preprocessor: ObsidianDataPreprocessor) -> np.ndarray:
        """遺跡での産地構成比予測"""
        pass

    @abstractmethod
    def predict_grid_ratios(self, preprocessor: ObsidianDataPreprocessor,
                          grid_points: Optional[np.ndarray] = None) -> np.ndarray:
        """グリッド点での産地構成比予測"""
        pass

    @property
    @abstractmethod
    def model_name(self) -> str:
        """モデル名（評価結果の識別用）"""
        pass
```

### 3.2 各モデルの修正内容

1. **NadarayaWatsonEstimator**
   - BaseCompositionModelを継承
   - 既存メソッド名の統一（predict_all_periods_origins → predict_site_ratios）

2. **BayesianNadarayaWatson**
   - 親クラスの変更対応

3. **KSBPModel**
   - BaseCompositionModelを継承
   - インターフェースに合わせたメソッド名の調整

4. **BayesianSpatialMultinomialModel**
   - BaseCompositionModelを継承
   - predict_probabilitiesをpredict_site_ratiosにリネーム

### 3.3 LOOCV評価の共通化

#### UnifiedLOOCVEvaluator（新クラス）
```python
class UnifiedLOOCVEvaluator:
    """全モデル共通のLOOCV評価器"""

    def __init__(self, models: List[BaseCompositionModel],
                 config: LOOCVConfig):
        self.models = models
        self.config = config

    def evaluate_all_models(self, preprocessor: ObsidianDataPreprocessor) -> Dict[str, pd.DataFrame]:
        """全モデルのLOOCV評価を実行"""
        results = {}
        for model in self.models:
            results[model.model_name] = self._evaluate_single_model(model, preprocessor)
        return results

    def _evaluate_single_model(self, model: BaseCompositionModel,
                             preprocessor: ObsidianDataPreprocessor) -> pd.DataFrame:
        """単一モデルのLOOCV評価"""
        # 既存のLOOCV実装を汎用化
        pass
```

### 3.4 モデル比較フレームワーク

#### ModelComparison（新クラス）
```python
class ModelComparison:
    """複数モデルの比較実験フレームワーク"""

    def __init__(self, models: List[BaseCompositionModel]):
        self.models = models
        self.evaluator = UnifiedLOOCVEvaluator(models)

    def run_comparison(self, preprocessor: ObsidianDataPreprocessor) -> ComparisonResults:
        """全モデルの学習・評価・比較を実行"""
        # 1. 各モデルの学習
        # 2. LOOCV評価
        # 3. 予測結果の比較
        # 4. 可視化
        pass
```

## 4. 実装手順

### Phase 1: 共通インターフェースの作成（優先度: 高）
1. BaseCompositionModelクラスの作成
2. 各モデルクラスへの継承実装

### Phase 2: LOOCV評価の共通化（優先度: 高）
1. UnifiedLOOCVEvaluatorの実装
2. 既存のLOOCVコードのリファクタリング

### Phase 3: モデル比較フレームワーク（優先度: 中）
1. ModelComparisonクラスの実装
2. 比較結果の可視化機能

### Phase 4: ドキュメント整備（優先度: 低）
1. 新しいAPIドキュメント
2. 使用例のノートブック作成

## 5. 期待される効果

1. **開発効率の向上**
   - 新しいモデルの追加が容易に
   - 共通評価の自動化

2. **実験の再現性**
   - 統一的な評価フレームワーク
   - 公平な比較が可能

3. **コードの保守性**
   - 重複コードの削減
   - 明確な責任分離

## 6. 注意事項

- 既存のノートブックとの互換性を考慮
- 段階的な移行を実施（既存コードと新コードの共存期間を設ける）
- 各段階でテストを実施

## 7. 新しいディレクトリ構造

### 7.1 提案するディレクトリ構造

```
bayesian_statistics/
├── __init__.py
├── models/                      # モデル関連の全コード
│   ├── __init__.py
│   ├── base/                   # 基底クラスとインターフェース
│   │   ├── __init__.py
│   │   ├── base_model.py       # BaseCompositionModel, BaseIntensityModel
│   │   └── interfaces.py       # その他のインターフェース定義
│   │
│   ├── preprocessing/          # 前処理関連
│   │   ├── __init__.py
│   │   ├── data_preprocessor.py  # 現model3_preprocessing.py
│   │   └── distance_calculator.py # 現10_2_calculate_distance.py
│   │
│   ├── composition/            # 産地構成比モデル
│   │   ├── __init__.py
│   │   ├── nadaraya_watson.py    # 現model3_nadaraya_watson.py
│   │   ├── bayesian_nw.py        # 現model3_bayes_nw.py
│   │   ├── ksbp.py               # 現model3_ksbp.py
│   │   └── spatial_regression.py # 現model3_bayesian_spatial.py
│   │
│   ├── intensity/              # 遺跡存在確率モデル
│   │   ├── __init__.py
│   │   ├── ipp.py                # 現model3_ipp.py
│   │   └── bayesian_ipp.py       # 将来の拡張用
│   │
│   ├── evaluation/             # 評価関連
│   │   ├── __init__.py
│   │   ├── loocv.py              # 汎用化されたLOOCV
│   │   ├── metrics.py            # 現model3_compositional_metrics.py
│   │   └── comparison.py         # モデル比較フレームワーク
│   │
│   ├── visualization/          # 可視化関連
│   │   ├── __init__.py
│   │   ├── base_visualizer.py    # 共通可視化インターフェース
│   │   ├── composition_plots.py  # 産地構成比モデル用プロット
│   │   ├── intensity_plots.py    # 遺跡存在確率モデル用プロット
│   │   ├── comparison_plots.py   # モデル比較用プロット
│   │   └── utils.py              # 現model3_visualization.pyの一部
│   │
│   ├── config/                 # 設定関連
│   │   ├── __init__.py
│   │   ├── model_config.py       # 現model3_config.py（設定部分）
│   │   └── pipeline.py           # 現model3_config.py（パイプライン部分）
│   │
│   └── utils/                  # 共通ユーティリティ
│       ├── __init__.py
│       ├── data_utils.py         # データ処理ユーティリティ
│       └── spatial_utils.py      # 空間データ処理
│
├── legacy/                     # 移行期間中の旧コード
│   ├── __init__.py
│   ├── utils_2.py               # IPPの旧実装
│   ├── utils_3.py               # その他の旧実装
│   └── model3_fixed_bayesian_spatial_multinomial.py
│
└── notebooks/                  # ノートブックは現状維持
```

### 7.2 各ファイルの移動マッピング

| 現在のファイル | 新しい場所 |
|--------------|-----------|
| model3_preprocessing.py | models/preprocessing/data_preprocessor.py |
| model3_nadaraya_watson.py | models/composition/nadaraya_watson.py |
| model3_bayes_nw.py | models/composition/bayesian_nw.py |
| model3_ksbp.py | models/composition/ksbp.py |
| model3_bayesian_spatial.py | models/composition/spatial_regression.py |
| model3_ipp.py | models/intensity/ipp.py |
| model3_loocv.py | models/evaluation/loocv.py |
| model3_compositional_metrics.py | models/evaluation/metrics.py |
| model3_visualization.py | models/visualization/utils.py (分割) |
| model3_config.py | models/config/model_config.py + pipeline.py |
| model3_fixed_bayesian_spatial_multinomial.py | legacy/ |
| utils_2.py | legacy/ |
| utils_3.py | legacy/ |

### 7.3 可視化機能の統合設計

#### BaseVisualizer（抽象基底クラス）
```python
class BaseVisualizer(ABC):
    """可視化の共通インターフェース"""
    
    @abstractmethod
    def plot_prediction_map(self, data: pl.DataFrame, **kwargs) -> plt.Figure:
        """予測結果の地図プロット"""
        pass
    
    @abstractmethod
    def plot_model_diagnostics(self, model_results: Dict[str, Any]) -> plt.Figure:
        """モデル診断プロット"""
        pass
```

#### CompositionVisualizer
```python
class CompositionVisualizer(BaseVisualizer):
    """産地構成比モデル用可視化"""
    
    def plot_origin_distribution(self, df: pl.DataFrame, origin: str, period: int) -> plt.Figure:
        """特定産地・時期の分布図"""
        pass
    
    def plot_composition_comparison(self, models: List[BaseCompositionModel]) -> plt.Figure:
        """複数モデルの予測結果比較"""
        pass
```

### 7.4 パイプラインの再設計

```python
# models/config/pipeline.py
class UnifiedPipeline:
    """統一的なモデル実行パイプライン"""
    
    def __init__(self, config: ModelConfig):
        self.config = config
        self.preprocessor = DataPreprocessor(config.data_dir)
        self.models = {}
        self.visualizers = {}
    
    def register_model(self, name: str, model: BaseCompositionModel):
        """モデルの登録"""
        self.models[name] = model
        
    def run_experiments(self) -> ExperimentResults:
        """全モデルの実験実行"""
        # 1. データ前処理
        # 2. 各モデルの学習
        # 3. LOOCV評価
        # 4. 結果の可視化
        # 5. 比較レポート生成
        pass
```

## 8. 実装の優先順位（更新）

### Phase 0: ディレクトリ構造の作成（優先度: 最高）
1. 新しいディレクトリ構造の作成
2. __init__.pyファイルの配置
3. インポートパスの設計

### Phase 1: 基底クラスとインターフェース（優先度: 高）
1. base/base_model.py の作成
2. visualization/base_visualizer.py の作成
3. 各インターフェースの定義

### Phase 2: 既存コードの移動とリファクタリング（優先度: 高）
1. 前処理モジュールの移動
2. 各モデルの移動と共通インターフェース実装
3. 可視化機能の分割と整理

### Phase 3: 統合機能の実装（優先度: 中）
1. 汎用LOOCV評価器
2. モデル比較フレームワーク
3. 統一パイプライン

### Phase 4: ドキュメントとテスト（優先度: 低）
1. 新しいAPIドキュメント
2. 移行ガイド
3. サンプルノートブック

## 9. 今後の拡張可能性

- 新しい評価指標の追加
- ハイパーパラメータ最適化の共通化
- 並列処理による高速化
- Webアプリケーション対応
- より高度な可視化（インタラクティブ地図など）
