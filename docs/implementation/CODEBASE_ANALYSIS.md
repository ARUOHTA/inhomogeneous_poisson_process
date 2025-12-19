# コードベース分析ドキュメント

このドキュメントは、マーク付き点過程モデル実装に向けた現行コードの詳細分析です。

## 1. 全体構造

```
bayesian_statistics/
├── bayesian_statistics/           # メインPythonパッケージ
│   ├── models/                    # 旧世代モデル群（統一API）
│   │   ├── base/                  # 抽象基底クラス
│   │   ├── composition/           # 産地構成比モデル
│   │   ├── intensity/             # 遺跡存在確率モデル（IPP）
│   │   ├── preprocessing/         # データ前処理
│   │   ├── evaluation/            # モデル評価
│   │   └── visualization/         # 可視化
│   │
│   └── nngp/                      # 新世代NNGPモジュール
│       ├── model/                 # NNGPモデル実装
│       │   ├── nngp.py            # Vecchia因子計算の核
│       │   ├── sample.py          # サンプリングユーティリティ
│       │   ├── multinomial_model.py  # 多項NNGP
│       │   └── distance_prior_model.py  # 距離事前付きNNGP（最新）
│       ├── mesh/                  # メッシュ処理
│       ├── poi/                   # POI処理
│       └── utils/                 # ユーティリティ
│
├── data/                          # データファイル
├── notebooks/                     # 実験ノートブック
└── docs/                          # ドキュメント
```

## 2. 二つのモデルエコシステム

現在、**2つの異なるモデル設計パラダイム**が共存しています。

### 2.1 旧世代: `models/` ディレクトリ

**設計思想**: 統一インターフェースによるモデル比較・評価

**基底クラス**:
- `BaseCompositionModel` - 産地構成比モデル用
- `BaseIntensityModel` - 遺跡存在確率モデル用

**特徴**:
- `fit()`, `predict_*()`, `model_name`, `results` の統一API
- `UnifiedLOOCVEvaluator` による自動評価
- Polarsベースのデータ入出力

**実装モデル**:
| モデル | ファイル | 説明 |
|--------|----------|------|
| NadarayaWatsonEstimator | composition/nadaraya_watson.py | カーネル回帰 |
| BayesianNadarayaWatson | composition/bayesian_nw.py | ベイズNW |
| KSBPModel | composition/ksbp.py | Kernel Stick-Breaking |
| BayesianSpatialMultinomialModel | composition/spatial_regression.py | 空間多項ロジット |
| InhomogeneousPoissonProcess | intensity/ipp.py | 非斉次ポアソン過程 |

### 2.2 新世代: `nngp/` ディレクトリ

**設計思想**: 高性能なNNGP実装に特化

**核となるコンポーネント**:
- `nngp.py` - Vecchia因子の構築（`NNGPFactors`, `FactorCache`）
- `sample.py` - Polya-Gamma サンプリング、ギブス更新
- `multinomial_model.py` - 多項NNGPモデル
- `distance_prior_model.py` - 距離依存事前付きモデル（**最新・推奨**）

**特徴**:
- 疎精度行列によるスケーラブルな推論
- Morton順序によるキャッシュ効率
- per-feature カーネルパラメータ

## 3. 核心コンポーネントの詳細

### 3.1 データ前処理: `ObsidianDataPreprocessor`

**場所**: `models/preprocessing/data_preprocessor.py`

**役割**: 全データの中央集権的管理

**主要メソッド**:
```python
class ObsidianDataPreprocessor:
    def load_data(self) -> Dict[str, pl.DataFrame]
        # df_elevation, df_obsidian, df_sites を読み込み

    def create_explanatory_variables(self, variable_names) -> Tuple[W_grids, W_sites]
        # グリッド・遺跡の説明変数行列

    def preprocess_obsidian_data(self, period, origin) -> Tuple[counts, target_counts]
        # 時期・産地別カウントデータ

    def create_grid_coords(self) -> np.ndarray  # (n_grid, 2) ラジアン
    def create_site_coords(self) -> np.ndarray  # (n_sites, 2) ラジアン
```

**依存データファイル**:
- `11_gdf_elevation.csv` - 標高・地形データ
- `11_gdf_obsidian.csv` - 黒曜石出土データ
- `11_gdf_sites.csv` - 遺跡メタデータ

### 3.2 NNGP核: `nngp.py`

**場所**: `nngp/model/nngp.py`

**主要データ構造**:
```python
@dataclass
class NNGPFactors:
    neighbor_idx: List[np.ndarray]  # 各点の近傍インデックス
    a_rows: List[np.ndarray]        # 条件付き平均係数 a_i
    d: np.ndarray                   # 条件付き分散 d_i
```

**主要関数**:
```python
def build_nngp_factors(coords, M, kernel, order=None) -> (NNGPFactors, order)
    # 観測点に対するVecchia因子構築

def build_cross_factors(prior_coords, target_coords, M, kernel) -> NNGPFactors
    # グリッド点への予測用因子

class FactorCache:
    # per-feature でNNGP因子をキャッシュ
    # Morton順序で近傍探索を効率化
    def factors_S_all(self) -> List[NNGPFactors]
    def A_grid_all(self) -> List[sp.csr_matrix]
```

### 3.3 サンプリングユーティリティ: `sample.py`

**場所**: `nngp/model/sample.py`

**Polya-Gammaサンプラー**:
```python
class PolyaGammaSampler:
    def sample(self, b, c) -> np.ndarray
        # PG(b, c) からサンプリング
```

**カーネル**:
```python
@dataclass
class LocalNNGPKernel:
    lengthscale: float
    variance: float
    def K(self, X1, X2) -> np.ndarray  # RBFカーネル
```

**ギブス更新**:
```python
def update_beta_category(
    beta_k, eta_k, W, factors_by_feature, order,
    omega, kappa_tilde, rng, prior_mean_by_feature=None
) -> (beta_new, eta_new)
    # カテゴリkの係数場を1サイト・1特徴ずつ更新
```

**多項ロジット関連**:
```python
def compute_eta(beta, W) -> np.ndarray          # 線形予測子
def softmax_with_baseline(eta) -> np.ndarray    # 確率変換
def compute_log_rest_terms(eta) -> np.ndarray   # log(1 + Σexp(η_ℓ))
```

### 3.4 多項NNGPモデル: `multinomial_model.py`

**場所**: `nngp/model/multinomial_model.py`

**データ構造**:
```python
@dataclass
class MultinomialNNGPConfig:
    n_iter: int = 2000
    burn_in: int = 500
    thinning: int = 1
    neighbor_count: int = 25
    kernel_lengthscale: float = 0.05
    kernel_variance: float = 1.0
    seed: Optional[int] = None

@dataclass
class MultinomialDataset:
    period: int
    origins: List[str]
    coords: np.ndarray        # (n_sites, 2)
    counts: np.ndarray        # (n_sites, K)
    total_counts: np.ndarray  # (n_sites,)
    design_matrix_sites: np.ndarray  # (n_sites, p+1)
    grid_points: np.ndarray   # (n_grid, 2)
    design_matrix_grid: np.ndarray   # (n_grid, p+1)

@dataclass
class MultinomialNNGPResults:
    dataset: MultinomialDataset
    config: MultinomialNNGPConfig
    factor_cache: FactorCache
    kernels: List[LocalNNGPKernel]
    beta_samples: np.ndarray  # (n_save, K-1, p+1, n_sites)
```

**MCMCループ構造**:
```python
for iteration in range(n_iter):
    for k in range(K-1):  # 各カテゴリ
        # 1. log_rest = log(1 + Σ_{ℓ≠k} exp(η_ℓ))
        # 2. psi = η_k - log_rest
        # 3. ω ~ PG(N_i, psi)
        # 4. κ̃ = (y_k - N/2) + ω * log_rest
        # 5. β_k 更新 (update_beta_category)
```

### 3.5 距離事前付きモデル: `distance_prior_model.py`

**場所**: `nngp/model/distance_prior_model.py`（**最新・推奨**）

**設計思想**: 距離情報を加法項ではなく**事前平均**として組み込み

```
β_0k(s) = λ_k * g_k(s) + δ_k(s)
          ↑距離事前      ↑データ駆動調整
```

**追加データ構造**:
```python
@dataclass
class DistancePriorDataset(MultinomialDataset):
    # 距離Z-score
    distance_zscores_sites: np.ndarray  # (n_sites, K)
    distance_zscores_grid: np.ndarray   # (n_grid, K)
    # 対数比特徴量
    distance_features_sites: np.ndarray # (n_sites, K-1)
    distance_features_grid: np.ndarray  # (n_grid, K-1)
    # 事前平均
    prior_mean_intercept_sites: np.ndarray  # (K-1, n_sites)
    prior_mean_intercept_grid: np.ndarray   # (K-1, n_grid)
```

**距離特徴量計算**:
```python
def compute_distance_features(Z, w, tau, alpha) -> np.ndarray
    # g_k = log(p_0k) - log(p_0K)
    # p_0 = weighted_inverse_softmax(Z, w, tau, alpha)
```

### 3.6 IPPモデル: `ipp.py`

**場所**: `models/intensity/ipp.py`

**クラス構造**:
```python
class InhomogeneousPoissonProcess(BaseIntensityModel):
    def __init__(self, variable_names, region)

    def poisson_thinning_U(self, intensity_func, valid_grids, volume)
        # U ~ IPP(λ*(1-q)) をthinningでサンプル

    def mcmc_sampler(self, ...)
        # 1. U のサンプリング (Poisson thinning)
        # 2. データ結合 (X ∪ U)
        # 3. ω ~ PG(1, W@β)
        # 4. β のサンプリング（正規完全条件付き）
        # 5. λ* ~ Gamma(m0+n, r0+|D|)
```

**重要**:
- 現在のIPPはNNGPを使っていない（通常のベイズロジスティック回帰）
- 空間ランダム効果なし
- マーク情報（産地構成比）なし

## 4. モジュール依存グラフ

```
ObsidianDataPreprocessor
         │
         ▼
┌────────────────────────────────────────────────────┐
│                                                    │
│  ┌──────────────────┐    ┌──────────────────────┐ │
│  │   models/        │    │   nngp/model/        │ │
│  │                  │    │                      │ │
│  │  IPP ─────────┐  │    │  nngp.py             │ │
│  │  (intensity/) │  │    │    │ NNGPFactors     │ │
│  │               │  │    │    │ FactorCache     │ │
│  │               │  │    │    ▼                 │ │
│  │               │  │    │  sample.py           │ │
│  │               │  │    │    │ PolyaGammaSampler│
│  │               │  │    │    │ update_beta_*   │ │
│  │               │  │    │    ▼                 │ │
│  │               │  │    │  multinomial_model.py│ │
│  │               │  │    │    ▼                 │ │
│  │               │  │    │  distance_prior_model│ │
│  │               │  │    │  (最新)              │ │
│  └───────────────┘  │    └──────────────────────┘ │
│                                                    │
│         ★ 現在は分離している                        │
│         → マーク付き点過程で統合が必要              │
└────────────────────────────────────────────────────┘
```

## 5. 現行コードとマーク付き点過程理論の対応

| 理論セクション | 現行実装 | 状況 |
|---------------|---------|------|
| Sec2: モデル定式化 | | |
| - IPP強度 λ*q(s,t) | ipp.py | ✓ 実装済み（NNGPなし） |
| - 多項ロジット構成比 | multinomial_model.py | ✓ 実装済み |
| - 距離事前NNGP | distance_prior_model.py | ✓ 実装済み |
| Sec3: 階層ベイズ構造 | | |
| - β事前分布 | sample.py | ✓ |
| - u事前分布（NNGP） | nngp.py | ✓ |
| Sec4: 擬似尤度・PG拡張 | | |
| - 偽不在点過程 U | ipp.py (poisson_thinning_U) | ✓ |
| - PG変換（二項） | ipp.py | ✓ |
| - PG変換（多項） | sample.py | ✓ |
| Sec5: 点過程完全条件付き | | |
| - (β_int, u_int) 更新 | ipp.py（NNGPなし） | △ 要拡張 |
| - λ* 更新 | ipp.py | ✓ |
| - U 更新 | ipp.py | ✓ |
| - ω 更新 | ipp.py | ✓ |
| Sec6: マーク完全条件付き | | |
| - (β_k, u_k) 更新 | multinomial_model.py | ✓ |
| - ξ 更新 | sample.py | ✓ |
| Sec7: ギブスアルゴリズム | | |
| - 点過程のみ | ipp.py | ✓ |
| - マークのみ | multinomial_model.py | ✓ |
| - **統合版** | **未実装** | ✗ 新規必要 |
| Sec8-9: DGP時間構造 | **未実装** | ✗ 将来拡張 |

## 6. 実装上の課題

### 6.1 IPPとNNGPの分離

現在のIPP (`ipp.py`) は:
- 通常のベイズロジスティック回帰
- 空間ランダム効果 u_int なし
- Vecchia/NNGP近似なし

→ **u_int(s,t) のNNGP実装が必要**

### 6.2 二つのPG拡張の統一

- 点過程側: `ipp.py` で `PG(1, ψ)`
- マーク側: `sample.py` で `PG(N_i, η_k)`

→ **同一フレームワーク内で両方を扱う設計が必要**

### 6.3 偽不在点とマーク評価

現在の偽不在 U は:
- `ipp.py` で IPP(λ*(1-q)) からサンプル
- マークモデルでは使用されていない

→ **U の座標でマークモデルは評価しない（U には y がない）**
→ **点過程部分は X ∪ U、マーク部分は X のみ**

### 6.4 座標系の不整合

- `ipp.py`: 度単位の座標
- `nngp/`: 度単位の座標（preprocessor経由）
- 内部計算: 一部でラジアン変換

→ **座標系の統一が必要**

## 7. 推奨モジュール構造（新規実装用）

```
bayesian_statistics/nngp/model/
├── nngp.py                    # 既存：Vecchia因子
├── sample.py                  # 既存：PGサンプラー、ギブス更新
├── multinomial_model.py       # 既存：マーク部分のみ
├── distance_prior_model.py    # 既存：距離事前付きマーク
│
└── marked_point_process/      # 新規：マーク付き点過程
    ├── __init__.py
    ├── config.py              # 統合設定
    ├── dataset.py             # 統合データセット
    ├── intensity.py           # 点過程部分（NNGP付きIPP）
    ├── composition.py         # マーク部分（既存を再利用）
    ├── sampler.py             # 統合ギブスサンプラー
    └── results.py             # 結果・予測
```

## 8. 次のステップ

1. **intensity.py の設計**
   - `ipp.py` を参考に NNGP 付き IPP を実装
   - `u_int(s,t)` の空間ランダム効果を追加
   - `FactorCache` を再利用

2. **sampler.py の設計**
   - 点過程更新とマーク更新を交互に実行
   - 偽不在 U のサンプリングを統合
   - 共通の iteration ループ

3. **dataset.py の設計**
   - `DistancePriorDataset` を拡張
   - 点過程用の追加フィールド（valid_grids, region など）

4. **既存コードの再利用方針**
   - `nngp.py`, `sample.py` はそのまま利用
   - `distance_prior_model.py` のデータ準備ロジックを再利用
   - `ipp.py` の `poisson_thinning_U` を移植

---

*作成日: 2024-12-19*
*目的: マーク付き点過程モデル実装の準備*
