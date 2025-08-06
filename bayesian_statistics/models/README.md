# Bayesian Statistics Models API Reference

## 目次

1. [概要](#概要)
2. [基底クラス](#基底クラス)
3. [産地構成比モデル](#産地構成比モデル)
4. [遺跡存在確率モデル](#遺跡存在確率モデル)
5. [前処理](#前処理)
6. [評価・比較](#評価比較)
7. [設定・パイプライン](#設定パイプライン)
8. [可視化](#可視化)
9. [完全な使用例](#完全な使用例)

---

## 概要

このディレクトリには、考古学的な黒曜石分布データを分析するための統一されたベイズ統計モデル群が含まれています。
---

## 基底クラス

### BaseCompositionModel

産地構成比を予測する全モデルの抽象基底クラス

```python
from bayesian_statistics.models.base import BaseCompositionModel

class BaseCompositionModel(ABC):
    @abstractmethod
    def fit(self, preprocessor, **kwargs) -> Dict[str, Any]:
        """モデル学習"""

    @abstractmethod
    def predict_site_ratios(self, preprocessor) -> Dict[str, pl.DataFrame]:
        """遺跡での産地構成比予測"""

    @abstractmethod
    def predict_grid_ratios(self, preprocessor) -> pl.DataFrame:
        """グリッド点での産地構成比予測"""

    @property
    @abstractmethod
    def model_name(self) -> str:
        """モデル名"""

    @property
    @abstractmethod
    def results(self) -> Dict[str, Any]:
        """学習結果"""
```

### BaseIntensityModel

遺跡存在確率を予測する全モデルの抽象基底クラス

```python
from bayesian_statistics.models.base import BaseIntensityModel

class BaseIntensityModel(ABC):
    @abstractmethod
    def fit(self, preprocessor, **kwargs) -> Dict[str, Any]:
        """モデル学習"""

    @abstractmethod
    def predict_probability(self, preprocessor) -> pl.DataFrame:
        """グリッド上での存在確率予測"""

    @abstractmethod
    def predict_intensity(self, coordinates: np.ndarray) -> np.ndarray:
        """指定座標での強度予測"""

    @property
    @abstractmethod
    def model_name(self) -> str:
        """モデル名"""

    @property
    @abstractmethod
    def results(self) -> Dict[str, Any]:
        """学習結果"""
```

---

## 産地構成比モデル

### NadarayaWatsonEstimator

クラシカルなNadaraya-Watsonカーネル回帰による産地構成比推定

```python
from bayesian_statistics.models.composition import NadarayaWatsonEstimator

# 初期化
nw_model = NadarayaWatsonEstimator(
    sigma=1000,              # グリッド用バンド幅
    sigma_for_sites=0.1,     # 遺跡用バンド幅
    zero_replacement=1e-6    # ゼロ値置換
)

# 学習
results = nw_model.fit(preprocessor)

# 予測
site_ratios = nw_model.predict_site_ratios(preprocessor)
# 返り値: Dict[str, pl.DataFrame] = {"0": DataFrame, "1": DataFrame, ...}

grid_ratios = nw_model.predict_grid_ratios(preprocessor)
# 返り値: pl.DataFrame with columns ['x', 'y', '比率_0_神津島', ...]

print(nw_model.model_name)  # "NadarayaWatson"
```

### BayesianNadarayaWatson

Nadaraya-Watsonのベイズ拡張版

```python
from bayesian_statistics.models.composition import BayesianNadarayaWatson

# NadarayaWatsonEstimatorを継承、同じAPIで不確実性定量化を追加
bayesian_nw = BayesianNadarayaWatson(
    sigma=1000,
    sigma_for_sites=0.1
)

results = bayesian_nw.fit(preprocessor)
site_ratios = bayesian_nw.predict_site_ratios(preprocessor)

print(bayesian_nw.model_name)  # "BayesianNadarayaWatson"
```

### KSBPModel

Kernel Stick-Breaking Processによるノンパラメトリックベイズモデル

```python
from bayesian_statistics.models.composition import KSBPModel

ksbp_model = KSBPModel(
    kappa_coords=0.002,      # 座標カーネルパラメータ
    kappa_costs=2.0,         # コストカーネルパラメータ
    kappa_elevation=10000.0, # 標高カーネルパラメータ
    kappa_angle=3.0,         # 角度カーネルパラメータ
    kappa_river=2.0,         # 河川カーネルパラメータ
    gamma=0.1,               # 集中パラメータ
    alpha0=0.1,              # ディリクレパラメータ
    j_max=80,                # 最大クラスター数
    n_iter=2000,             # イテレーション数
    burnin=100               # バーンイン
)

results = ksbp_model.fit(preprocessor)
site_ratios = ksbp_model.predict_site_ratios(preprocessor)

print(ksbp_model.model_name)  # "KSBP"
```

### BayesianSpatialMultinomialModel

ベイズ空間多項ロジット回帰モデル（全時期対応）

```python
from bayesian_statistics.models.composition import BayesianSpatialMultinomialModel

spatial_model = BayesianSpatialMultinomialModel(
    variable_names=['average_elevation'],  # 使用する説明変数
    prior_sigma=0.5,                       # 事前分布の標準偏差
    n_draws=1000,                          # MCMC描画数
    n_tune=500,                            # チューニング回数
    random_seed=42                         # 乱数シード
)

# 全時期（0-4）のモデル学習
results = spatial_model.fit(preprocessor)
print(f"学習完了: {len(results)}時期のモデル")

# 遺跡での予測（時期別に分かれたDict形式）
site_ratios = spatial_model.predict_site_ratios(preprocessor)
# 返り値: {"0": DataFrame, "1": DataFrame, ...}

# グリッドでの予測（全時期・全産地のカラムを含む）
grid_ratios = spatial_model.predict_grid_ratios(preprocessor)
# 返り値: pl.DataFrame with columns ['x', 'y', 'ratio_0_神津島', 'ratio_1_神津島', ...]

print(spatial_model.model_name)  # "BayesianSpatialRegression"

# 学習結果の確認
for period, period_result in spatial_model.results.items():
    print(f"時期{period}: Rhat最大値 = {period_result['rhat_max']:.3f}")
```

**重要な変更点（v2.0）**：

- **全時期対応**: 従来の時期0のみから、全時期（0-4: 早期・早々期〜晩期）に対応
- **時期別学習**: 各時期で独立したMCMCサンプリングを実行
- **統一インターフェース**: 他の構成比モデル（`BayesianNadarayaWatson`等）と同じAPI

---

## 遺跡存在確率モデル

### InhomogeneousPoissonProcess

非斉次ポアソン過程による遺跡存在確率の推定

```python
from bayesian_statistics.models.intensity import InhomogeneousPoissonProcess

# 初期化
ipp_model = InhomogeneousPoissonProcess(
    variable_names=['average_elevation'],
    region=[[138.0, 141.0], [34.0, 37.0]]  # [x_min, x_max], [y_min, y_max]
)

# 学習（MCMC）
results = ipp_model.fit(preprocessor)

# 予測
probability_df = ipp_model.predict_probability(preprocessor)
# 返り値: pl.DataFrame with columns ['x', 'y', 'probability']

# 任意座標での強度予測
coordinates = np.array([[139.0, 35.0], [140.0, 36.0]])  # [x, y] in degrees
intensities = ipp_model.predict_intensity(coordinates)

# 診断
idata = ipp_model.create_inference_data()

print(ipp_model.model_name)  # "InhomogeneousPoissonProcess"
```

---

## 前処理

### ObsidianDataPreprocessor

全データの中央集権的前処理

```python
from bayesian_statistics.models.preprocessing import ObsidianDataPreprocessor

# 初期化
preprocessor = ObsidianDataPreprocessor(
    data_dir="/path/to/data",
    x_min=138.0, x_max=141.0,
    y_min=34.0, y_max=37.0
)

# データ読み込み
preprocessor.load_data()

# 学習用データ作成
X, Y = preprocessor.create_X_Y(
    variable_names=['average_elevation'],
    time_periods={0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"},
    origins=["神津島", "信州", "箱根", "高原山"]
)

# グリッド座標取得
grid_coords = preprocessor.create_grid_coords()

# 遺跡座標取得
site_coords = preprocessor.convert_to_grid_coords()

# グリッド情報
grid_info = preprocessor.create_grid_info()
# 返り値: {"n_grid_x": int, "n_grid_y": int, "delta_x": float, ...}

# データアクセス
elevation_df = preprocessor.df_elevation  # 標高データ
sites_df = preprocessor.df_sites          # 遺跡データ
obsidian_df = preprocessor.df_obsidian    # 黒曜石データ
```

---

## 評価・比較

### UnifiedLOOCVEvaluator

全産地構成比モデル対応の統一LOOCV評価

```python
from bayesian_statistics.models.evaluation import UnifiedLOOCVEvaluator, LOOCVConfig
from bayesian_statistics.models.composition import NadarayaWatsonEstimator

# モデル準備
nw_model = NadarayaWatsonEstimator(sigma=1000)

# 評価器初期化
loocv_config = LOOCVConfig(n_trials=100, random_seed=42)
evaluator = UnifiedLOOCVEvaluator(
    config=loocv_config,
    variable_names=['average_elevation']
)

# 評価実行
evaluation_results = evaluator.evaluate_single_model(nw_model, preprocessor)
# 返り値: pl.DataFrame with columns ['遺跡ID', '時期', '産地', 'observed_ratio',
#                                   'predicted_ratio', 'aitchison_distance',
#                                   'bray_curtis', 'jensen_shannon', 'total_variation']

# 要約統計
summary_stats = evaluation_results.group_by(['時期', '産地']).agg([
    pl.col('aitchison_distance').mean().alias('mean_aitchison'),
    pl.col('bray_curtis').mean().alias('mean_bray_curtis')
])
```

### ModelComparison

複数モデルの自動比較フレームワーク

```python
from bayesian_statistics.models.evaluation import ModelComparison, LOOCVConfig
from bayesian_statistics.models.composition import (
    NadarayaWatsonEstimator, BayesianNadarayaWatson, KSBPModel
)

# モデルリスト準備
models = [
    NadarayaWatsonEstimator(sigma=1000),
    BayesianNadarayaWatson(sigma=1000),
    KSBPModel(n_iter=1000)
]

# 比較実行
comparison = ModelComparison(
    models=models,
    variable_names=['average_elevation']
)

results = comparison.run_comparison(preprocessor)

# 結果可視化
results.plot_loocv_comparison()

# 要約統計取得
summary = results.get_summary_stats()
print(summary)
```

---

## 設定・パイプライン

### Model3Config

統一設定管理

```python
from bayesian_statistics.models.config import Model3Config

config = Model3Config(
    # データパス
    data_dir="/path/to/data",

    # グリッド設定
    x_min=138.0, x_max=141.0,
    y_min=34.0, y_max=37.0,

    # NW推定設定
    nw_sigma=1000,
    nw_sigma_for_sites=0.1,

    # IPP設定
    mcmc_iterations=30000,
    burn_in=5000,
    ipp_variable_names=['average_elevation'],

    # KSBP設定
    ksbp_n_iter=2000,
    ksbp_j_max=80,

    # その他
    time_periods={0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"},
    origins=["神津島", "信州", "箱根", "高原山"]
)
```

### Model3Pipeline

エンドツーエンド実行パイプライン

```python
from bayesian_statistics.models.config import Model3Pipeline, Model3Config

# 設定とパイプライン初期化
config = Model3Config(data_dir="/path/to/data")
pipeline = Model3Pipeline(config)

# 段階的実行
preprocessor = pipeline.run_preprocessing()
nw_model = pipeline.run_nadaraya_watson(preprocessor)
ipp_model = pipeline.run_ipp(preprocessor)

# 一括実行
results = pipeline.run_full_pipeline()
# 返り値: {
#   "preprocessor": ObsidianDataPreprocessor,
#   "nw_estimator": NadarayaWatsonEstimator,
#   "ipp_model": InhomogeneousPoissonProcess,
#   "results": Dict[str, Any]
# }

# 結果保存
pipeline.save_results("/path/to/output")
```

---

## 可視化

### ObsidianVisualizer

包括的可視化システム

```python
from bayesian_statistics.models.visualization import ObsidianVisualizer
import matplotlib.pyplot as plt

# 産地構成比マップ
# 予測結果を取得
grid_ratios = model.predict_grid_ratios(preprocessor)

fig, ax = ObsidianVisualizer.plot_ratio_map_from_grid_data(
    preprocessor.df_elevation,
    preprocessor.df_sites,
    grid_ratios,
    period=0,           # 早期・早々期
    origin="神津島",      # 神津島産地
    figsize=(12, 10),
    n_levels=20,
    cmap="Reds"
)
plt.show()

# 予測結果から直接可視化
nw_results = nw_model.predict_all_periods_origins(
    preprocessor, time_periods, origins
)
fig, ax = ObsidianVisualizer.plot_ratio_map_from_result(
    preprocessor.df_elevation,
    preprocessor.df_sites,
    ratio_mesh=nw_results["ratio_mesh"],
    ratio_sites=nw_results["ratio_sites"],
    period=0,
    origin="神津島",
    grid_x=nw_results["grid_x"],
    grid_y=nw_results["grid_y"]
)

# 遺跡存在確率マップ
probability_df = ipp_model.predict_probability(preprocessor)
grid_info = preprocessor.create_grid_info()
probability_mesh = probability_df["probability"].to_numpy().reshape(
    grid_info["n_grid_y"], grid_info["n_grid_x"]
)

x_coords = np.linspace(config.x_min, config.x_max, grid_info["n_grid_x"])
y_coords = np.linspace(config.y_min, config.y_max, grid_info["n_grid_y"])
grid_x, grid_y = np.meshgrid(x_coords, y_coords)

fig, ax = ObsidianVisualizer.plot_site_probability(
    preprocessor.df_elevation,
    preprocessor.df_sites,
    probability_mesh,
    grid_x, grid_y,
    period=0
)

# MCMC診断
idata = ipp_model.create_inference_data()
fig = ObsidianVisualizer.plot_mcmc_diagnostics(idata)

# 一括可視化
figures = ObsidianVisualizer.plot_all_periods_origins(
    preprocessor.df_elevation,
    preprocessor.df_sites,
    time_periods={0: "早期・早々期", 1: "前期"},
    origins=["神津島", "信州"],
    output_dir="output"
)
```

---

## 完全な使用例

### シナリオ1: 単一モデルでの分析

```python
import numpy as np
import matplotlib.pyplot as plt
from bayesian_statistics.models.config import Model3Config, Model3Pipeline
from bayesian_statistics.models.composition import NadarayaWatsonEstimator
from bayesian_statistics.models.evaluation import UnifiedLOOCVEvaluator, LOOCVConfig
from bayesian_statistics.models.visualization import ObsidianVisualizer

# 1. 設定
config = Model3Config(
    data_dir="/path/to/data",
    nw_sigma=1000,
    x_min=138, x_max=141, y_min=34, y_max=37
)

# 2. 前処理
pipeline = Model3Pipeline(config)
preprocessor = pipeline.run_preprocessing()

# 3. モデル学習
nw_model = NadarayaWatsonEstimator(sigma=config.nw_sigma)
results = nw_model.fit(preprocessor)

# 4. 予測
site_ratios = nw_model.predict_site_ratios(preprocessor)
grid_ratios = nw_model.predict_grid_ratios(preprocessor)

# 5. 評価
loocv_config = LOOCVConfig(n_trials=100, random_seed=42)
evaluator = UnifiedLOOCVEvaluator(
    config=loocv_config,
    variable_names=['average_elevation']
)
evaluation_results = evaluator.evaluate_single_model(nw_model, preprocessor)

# 6. 可視化
# 予測結果を取得
grid_ratios = nw_model.predict_grid_ratios(preprocessor)

# 新しいメソッドで可視化
fig, ax = ObsidianVisualizer.plot_ratio_map_from_grid_data(
    preprocessor.df_elevation,
    preprocessor.df_sites,
    grid_ratios,
    period=0,
    origin="神津島"
)
plt.show()

# 7. 結果分析
print(f"Model: {nw_model.model_name}")
print(f"Mean Aitchison Distance: {evaluation_results['aitchison_distance'].mean():.3f}")
```

### シナリオ2: 複数モデル比較

```python
from bayesian_statistics.models.composition import (
    NadarayaWatsonEstimator, BayesianNadarayaWatson,
    KSBPModel, BayesianSpatialMultinomialModel
)
from bayesian_statistics.models.evaluation import ModelComparison, LOOCVConfig

# 1. 設定・前処理（同上）
config = Model3Config(data_dir="/path/to/data")
pipeline = Model3Pipeline(config)
preprocessor = pipeline.run_preprocessing()

# 2. モデル群準備
models = [
    NadarayaWatsonEstimator(sigma=1000),
    BayesianNadarayaWatson(sigma=1000),
    KSBPModel(n_iter=1000, j_max=50),
    BayesianSpatialMultinomialModel(
        variable_names=['average_elevation'],
        n_draws=1000,
        n_tune=500
    )
]

# 3. 一括比較
comparison = ModelComparison(models, variable_names=['average_elevation'])
results = comparison.run_comparison(preprocessor)

# 4. 結果可視化・分析
results.plot_loocv_comparison()
plt.show()

summary = results.get_summary_stats()
print("Model Performance Summary:")
print(summary)

# 5. ベストモデル選択
best_model_name = summary.filter(
    pl.col('metric') == 'aitchison_distance'
).sort('mean_value').row(0)[0]
print(f"Best model: {best_model_name}")
```

### シナリオ3: IPPモデルとの統合分析

```python
from bayesian_statistics.models.intensity import InhomogeneousPoissonProcess

# 1. 産地構成比 + 存在確率の統合分析
config = Model3Config(data_dir="/path/to/data")
pipeline = Model3Pipeline(config)

# 2. 一括実行
results = pipeline.run_full_pipeline()
preprocessor = results["preprocessor"]
nw_model = results["nw_estimator"]
ipp_model = results["ipp_model"]

# 3. 統合可視化
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# 産地構成比（神津島・早期）
ObsidianVisualizer.plot_ratio_map(
    preprocessor.df_elevation, preprocessor.df_sites,
    period=0, origin="神津島"
)

# 遺跡存在確率
probability_df = ipp_model.predict_probability(preprocessor)
grid_info = preprocessor.create_grid_info()
probability_mesh = probability_df["probability"].to_numpy().reshape(
    grid_info["n_grid_y"], grid_info["n_grid_x"]
)
x_coords = np.linspace(config.x_min, config.x_max, grid_info["n_grid_x"])
y_coords = np.linspace(config.y_min, config.y_max, grid_info["n_grid_y"])
grid_x, grid_y = np.meshgrid(x_coords, y_coords)

ObsidianVisualizer.plot_site_probability(
    preprocessor.df_elevation, preprocessor.df_sites,
    probability_mesh, grid_x, grid_y, period=0
)

plt.tight_layout()
plt.show()

# 4. 結果保存
pipeline.save_results("/path/to/output")
```

---

## 拡張ガイド

### 新しい産地構成比モデルの追加

1. `BaseCompositionModel`を継承
2. 必須メソッドの実装: `fit()`, `predict_site_ratios()`, `predict_grid_ratios()`, `model_name`, `results`
3. `composition/`ディレクトリにファイル配置
4. 自動的に`UnifiedLOOCVEvaluator`と`ModelComparison`で使用可能

### 新しい遺跡存在確率モデルの追加

1. `BaseIntensityModel`を継承
2. 必須メソッドの実装: `fit()`, `predict_probability()`, `predict_intensity()`, `model_name`, `results`
3. `intensity/`ディレクトリにファイル配置

### カスタム評価指標の追加

1. `evaluation/metrics.py`に新しい関数を追加
2. `UnifiedLOOCVEvaluator`の評価ループに組み込み

---

## 注意事項

- 全てのDataFrameはPolars形式を使用
- 座標系: 度単位（経度、緯度）
- 時期インデックス: 0-4（早期・早々期〜晩期）
- 産地: ["神津島", "信州", "箱根", "高原山"]が標準
