# Model 3rd リファクタリング設計

## 概要
`notebooks/16_model_3rd.ipynb`のコードを整理し、前処理クラスとモデルクラスに分離して保守性を向上させる。

## 全体的なアーキテクチャ

### 1. モジュール構成
```
bayesian_statistics/
├── model3_preprocessing.py  # 前処理クラス
├── model3_nadaraya_watson.py  # Nadaraya-Watson推定量クラス
├── model3_ipp.py  # IPP (Inhomogeneous Poisson Process) クラス
└── model3_visualization.py  # 可視化関数
```

### 2. クラス設計

#### 2.1 前処理クラス (`ObsidianDataPreprocessor`)

**責務:**
- 生データの読み込みと検証
- グリッド座標系への変換
- Tobler距離の計算と読み込み
- 説明変数の作成と標準化

**主要メソッド:**
```python
class ObsidianDataPreprocessor:
    def __init__(self, data_dir: str)
    def load_data(self) -> Dict[str, pl.DataFrame]
    def convert_to_grid_coords(self, df_obsidian: pl.DataFrame) -> np.ndarray
    def load_tobler_distances(self) -> np.ndarray
    def create_explanatory_variables(self, variable_names: List[str]) -> Tuple[np.ndarray, np.ndarray]
    def create_grid_info(self) -> Dict[str, float]
    def create_land_mask(self) -> np.ndarray
    def preprocess_obsidian_data(self, target_period: int, target_origin: str) -> Tuple[np.ndarray, np.ndarray]
```

#### 2.2 Nadaraya-Watson推定量クラス (`NadarayaWatsonEstimator`)

**責務:**
- カーネル重みの計算
- 重み付き比率の計算
- 全時期・全産地の推定実行

**主要メソッド:**
```python
class NadarayaWatsonEstimator:
    def __init__(self, sigma: float = 500, sigma_for_sites: float = 0.1)
    def calculate_kernel_weights(self, distances: np.ndarray) -> np.ndarray
    def calculate_weighted_ratios(self, weights: np.ndarray, counts: np.ndarray, target_counts: np.ndarray) -> np.ndarray
    def fit(self, preprocessor: ObsidianDataPreprocessor) -> Dict[str, np.ndarray]
    def predict_single(self, target_period: int, target_origin: str) -> np.ndarray
    def predict_all_periods_origins(self) -> pl.DataFrame
```

#### 2.3 IPPクラス (`InhomogeneousPoissonProcess`)

**責務:**
- 強度関数の定義と更新
- MCMC推論の実行
- 事後分布の計算

**主要メソッド:**
```python
class InhomogeneousPoissonProcess:
    def __init__(self, variable_names: List[str], region: List[List[float]])
    def create_design_matrix(self, xy: np.ndarray) -> np.ndarray
    def poisson_thinning_U(self, valid_grids: np.ndarray, volume: float) -> Tuple[np.ndarray, np.ndarray]
    def mcmc_sampler(self, X_obs: np.ndarray, y_obs: np.ndarray, num_iterations: int) -> Tuple[np.ndarray, np.ndarray]
    def fit(self, preprocessor: ObsidianDataPreprocessor, num_iterations: int = 30000) -> Dict[str, np.ndarray]
    def predict_site_probability(self, grid_coords: np.ndarray) -> np.ndarray
```

#### 2.4 可視化クラス (`ObsidianVisualizer`)

**責務:**
- コンター図の作成
- 結果の可視化
- MCMC診断プロット

**主要メソッド:**
```python
class ObsidianVisualizer:
    @staticmethod
    def plot_contour(df: pl.DataFrame, value_col: str, **kwargs) -> Tuple[plt.Figure, plt.Axes]
    @staticmethod
    def plot_ratio_map(df_elevation: pl.DataFrame, df_sites: pl.DataFrame, period: int, origin: str)
    @staticmethod
    def plot_mcmc_diagnostics(idata: az.InferenceData)
    @staticmethod
    def create_grid_visualization(df: pl.DataFrame) -> plt.Figure
```

### 3. データフロー

```
1. データ読み込み (ObsidianDataPreprocessor)
   ↓
2. 前処理実行 (ObsidianDataPreprocessor)
   ↓
3a. NW推定 (NadarayaWatsonEstimator) → 産地構成比
   ↓
3b. IPP推定 (InhomogeneousPoissonProcess) → 遺跡存在確率
   ↓
4. 結果可視化 (ObsidianVisualizer)
```

### 4. 設定管理

#### 4.1 設定クラス (`Model3Config`)
```python
@dataclass
class Model3Config:
    # データパス
    data_dir: str

    # グリッド設定
    x_min: float = 138
    x_max: float = 141
    y_min: float = 34
    y_max: float = 37

    # NW推定設定
    nw_sigma: float = 500
    nw_sigma_for_sites: float = 0.1

    # IPP設定
    mcmc_iterations: int = 30000
    burn_in: int = 5000

    # 説明変数
    variable_names: List[str] = field(default_factory=lambda: [
        "average_elevation",
        "average_slope_angle",
        "cost_shinshu",
        "cost_river"
    ])

    # 時期・産地設定
    time_periods: Dict[int, str] = field(default_factory=lambda: {
        0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"
    })
    origins: List[str] = field(default_factory=lambda: [
        "神津島", "信州", "箱根", "高原山"
    ])
```

### 5. メインの実行クラス (`Model3Pipeline`)

```python
class Model3Pipeline:
    def __init__(self, config: Model3Config)
    def run_preprocessing(self) -> ObsidianDataPreprocessor
    def run_nadaraya_watson(self, preprocessor: ObsidianDataPreprocessor) -> NadarayaWatsonEstimator
    def run_ipp(self, preprocessor: ObsidianDataPreprocessor) -> InhomogeneousPoissonProcess
    def run_full_pipeline(self) -> Dict[str, Any]
    def save_results(self, output_dir: str)
```

### 6. 使用例

```python
# 設定
config = Model3Config(data_dir="/path/to/data")

# パイプライン実行
pipeline = Model3Pipeline(config)
results = pipeline.run_full_pipeline()

# 個別実行の例
preprocessor = pipeline.run_preprocessing()
nw_estimator = pipeline.run_nadaraya_watson(preprocessor)
ipp_model = pipeline.run_ipp(preprocessor)

# 結果保存
pipeline.save_results("output/")
```

### 7. 実装上の考慮点

#### 7.1 パフォーマンス
- 重い計算（距離行列、カーネル重み）のキャッシュ機能
- 並列処理の活用（可能な箇所）
- メモリ効率的な実装

#### 7.2 エラーハンドリング
- データ検証機能の強化
- 適切な例外処理とログ出力
- 中間結果の保存機能

#### 7.3 テスタビリティ
- 各クラスの単体テスト
- モックデータでの統合テスト
- 回帰テスト用のベンチマーク

#### 7.4 拡張性
- 新しい説明変数の追加が容易
- 異なるカーネル関数の選択肢
- 他の推定手法への拡張

### 8. 既存コードからの移行

#### 8.1 段階的移行
1. 前処理部分のクラス化
2. NW推定部分のクラス化
3. IPP推定部分のクラス化
4. 可視化部分の整理
5. 統合とテスト

#### 8.2 後方互換性
- 既存の関数をラッパーとして残す
- ノートブックからの段階的移行
- 結果の一致確認

この設計により、コードの保守性、再利用性、テスタビリティが大幅に向上し、今後の拡張も容易になります。
