# マーク付き点過程モデル実装計画

## 概要

このドキュメントは、`docs/sections/sec2.tex` - `sec9.tex` で定式化したマーク付き点過程モデルを、既存コードベースに統合実装するための計画です。

---

## 1. 理論と実装の対応表

### 1.1 モデル構造

```
【理論】                          【実装】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
遺跡位置 X ~ IPP(λ*q)            marked_point_process/intensity.py
  ├─ 強度 λ* ~ Gamma             - lambda_star パラメータ
  ├─ 存在確率 q(s,t)             - logistic(η_int)
  │    └─ η_int = W_int'β_int + u_int
  └─ 空間効果 u_int ~ NNGP       - FactorCache 利用

マーク y_i ~ Multinomial(N_i, π)  marked_point_process/composition.py
  ├─ 構成比 π_k(s,t)             - softmax(η_k)
  │    └─ η_k = W_z'β_k + u_k
  └─ 空間効果 u_k ~ NNGP         - FactorCache 利用（共有可）

偽不在 U ~ IPP(λ*(1-q))          marked_point_process/intensity.py
  └─ Poisson thinning            - 既存 ipp.py から移植
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### 1.2 ギブスサンプラーのステップ（sec7.tex対応）

```
【Sec7: ギブスアルゴリズム】      【実装場所】                     【評価点】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
(a) U ~ IPP(λ*(1-q))             intensity.py: sample_pseudo_absence()   グリッド全体
(b) ω_i ~ PG(1, η_int,i)         sampler.py: sample_omega_intensity()    X ∪ U（n点）
(c) (β_int, u_int) ~ Normal      sampler.py: update_intensity_params()   X ∪ U（n点）
(d) λ* ~ Gamma(m0+n, r0+|D|)     sampler.py: update_lambda_star()        -
(e) ξ_ik ~ PG(N_i, η_ik)         sampler.py: sample_xi_marks()           X のみ（n_X点）
(f) (β_k, u_k) ~ Normal          sampler.py: update_mark_params()        X のみ（n_X点）
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### 1.3 重要：点過程とマークの評価点の違い

```
┌─────────────────────────────────────────────────────────────────┐
│ 【点過程部分】                                                    │
│   評価点: X ∪ U = {(s_i, t_i) : i=1,...,n}                       │
│   - n = n_X + n_U （反復ごとにn_Uは変化）                         │
│   - 観測遺跡 X: y_i = 1                                          │
│   - 偽不在 U: y_i = 0                                            │
│   - NNGP因子は n 点で構築                                         │
├─────────────────────────────────────────────────────────────────┤
│ 【マーク部分】                                                    │
│   評価点: X のみ = {(s_i, t_i) : i=1,...,n_X}                     │
│   - 偽不在 U にはマークデータ（産地出土数）がない                  │
│   - NNGP因子は n_X 点で構築                                       │
│   - 既存の multinomial_model.py の構造をそのまま利用可能           │
└─────────────────────────────────────────────────────────────────┘
```

---

## 2. 各ギブスステップの数式と実装詳細

### 2.1 ステップ(a): 偽不在のサンプリング（sec5.tex §5.3）

**数式:**
```
U | λ*, β_int, u_int ~ IPP(λ*(1-q))

where q(s,t) = exp(η_int) / (1 + exp(η_int))
      η_int(s,t) = w_int(s,t)' β_int + u_int(s,t)
```

**Poisson Thinning アルゴリズム:**
1. N ~ Poisson(λ* |D|) をサンプル
2. j=1,...,N について (s_j, t_j) ~ Uniform(D) をサンプル
3. u_j ~ Uniform(0,1) をサンプルし、u_j < 1-q(s_j,t_j) なら採用

**実装:**
```python
def sample_pseudo_absence(
    lambda_star: float,
    eta_int_grid: np.ndarray,  # グリッド上の η_int
    grid_coords: np.ndarray,
    valid_mask: np.ndarray,
    region_volume: float,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns:
        U_coords: (n_U, 2) 偽不在点の座標
        U_design: (n_U, p_int+1) 偽不在点の設計行列
    """
    q_grid = 1.0 / (1.0 + np.exp(-eta_int_grid))

    # 1. 候補点数をサンプル
    n_candidates = rng.poisson(lambda_star * region_volume)

    # 2. 有効グリッド上でランダムに候補点を選択
    valid_indices = np.where(valid_mask)[0]
    candidate_indices = rng.choice(valid_indices, size=n_candidates)

    # 3. Thinning: 1-q(s,t) の確率で採用
    u = rng.uniform(size=n_candidates)
    accept_mask = u < (1.0 - q_grid[candidate_indices])

    U_coords = grid_coords[candidate_indices[accept_mask]]
    # ... 設計行列の構築
    return U_coords, U_design
```

### 2.2 ステップ(b): 点過程側PGサンプリング（sec5.tex §5.4）

**数式:**
```
ω_i | η_int,i ~ PG(1, η_int,i)    for i = 1,...,n (X ∪ U)

where κ_i = y_i - 1/2
      y_i = 1 (i ∈ X), 0 (i ∈ U)
```

**実装:**
```python
def sample_omega_intensity(
    eta_int: np.ndarray,  # (n,) X ∪ U 全点
    pg_sampler: PolyaGammaSampler,
) -> np.ndarray:
    """
    Returns:
        omega: (n,) PG(1, η_int,i) からのサンプル
    """
    n = len(eta_int)
    b = np.ones(n)  # 点過程は b=1
    return pg_sampler.sample(b, eta_int)
```

### 2.3 ステップ(c): 点過程パラメータの更新（sec5.tex §5.1）

**数式:**
```
θ_int = [β_int; u_int]' | ω, X, U ~ N(m_int, V_int)

where:
  V_int^{-1} = H' Ω H + R_0
  m_int = V_int (H' Ω z + R_0 μ_0)

  H = [W_int  I_n]
  Ω = diag(ω_1, ..., ω_n)
  z = (κ_1/ω_1, ..., κ_n/ω_n)'
  κ_i = y_i - 1/2

  R_0 = [B_0^{-1}  0; 0  Q_int]  (事前精度)
  μ_0 = [b_0; 0]  (事前平均)
```

**注意:** Q_int はNNGP精度行列。n点で構築するため、反復ごとに再構築が必要。

**実装:**
```python
def update_intensity_params(
    omega: np.ndarray,         # (n,) PG変数
    y: np.ndarray,             # (n,) 0/1
    W_int: np.ndarray,         # (n, p_int+1) 設計行列
    factors_int: NNGPFactors,  # NNGP因子
    beta_prior_mean: np.ndarray,
    beta_prior_prec: np.ndarray,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns:
        beta_int: (p_int+1,) 固定効果
        u_int: (n,) 空間効果
    """
    n = len(omega)
    kappa = y - 0.5
    z = kappa / omega
    Omega = np.diag(omega)

    # NNGP精度行列
    Q_int = factors_int.precision_matrix()

    # ブロック行列の構築と解法
    # ... 効率的な実装が必要（スパース行列利用）
```

### 2.4 ステップ(d): λ* の更新（sec5.tex §5.2）

**数式:**
```
λ* | X, U ~ Gamma(m_0 + n, r_0 + |D|)

where n = n_X + n_U
```

**実装:**
```python
def update_lambda_star(
    n_total: int,          # n_X + n_U
    region_volume: float,  # |D|
    prior_shape: float,    # m_0
    prior_rate: float,     # r_0
    rng: np.random.Generator,
) -> float:
    shape = prior_shape + n_total
    rate = prior_rate + region_volume
    return rng.gamma(shape, 1.0 / rate)
```

### 2.5 ステップ(e): マーク側PGサンプリング（sec6.tex §6.3）

**数式:**
```
ξ_ik | η_ik, N_i ~ PG(N_i, η_ik)    for i = 1,...,n_X, k = 1,...,K-1

where κ̃_ik = y_ik - N_i/2
      N_i = Σ_k y_ik（総出土数）
```

**注意:** 点過程とは異なり、b = N_i（総出土数）

**実装:**
```python
def sample_xi_marks(
    eta_marks: np.ndarray,     # (K-1, n_X) 各カテゴリの線形予測子
    total_counts: np.ndarray,  # (n_X,) 各遺跡の総出土数 N_i
    pg_sampler: PolyaGammaSampler,
) -> np.ndarray:
    """
    Returns:
        xi: (K-1, n_X) PG(N_i, η_ik) からのサンプル
    """
    K_minus_1, n_X = eta_marks.shape
    xi = np.zeros_like(eta_marks)
    for k in range(K_minus_1):
        xi[k] = pg_sampler.sample(total_counts, eta_marks[k])
    return xi
```

### 2.6 ステップ(f): マークパラメータの更新（sec6.tex §6.2）

**数式:**
```
θ_k = [β_k; u_k]' | Ξ, Y, X ~ N(m_k, V_k)    for k = 1,...,K-1

where:
  V_k^{-1} = H_z' Ω_k H_z + R_{0,k}
  m_k = V_k (H_z' Ω_k z_k + R_{0,k} μ_{0,k})

  H_z = [W_z  I_{n_X}]
  Ω_k = diag(ξ_1k, ..., ξ_{n_X,k})
  z_k = (κ̃_1k/ξ_1k, ..., κ̃_{n_X,k}/ξ_{n_X,k})'
  κ̃_ik = y_ik - N_i/2
```

**注意:**
- n_X 点のみで評価（Uは含まない）
- 既存の `update_beta_category` をほぼそのまま利用可能

**実装:**
```python
def update_mark_params(
    xi_k: np.ndarray,          # (n_X,) カテゴリkのPG変数
    counts_k: np.ndarray,      # (n_X,) カテゴリkの出土数
    total_counts: np.ndarray,  # (n_X,) 総出土数
    W_z: np.ndarray,           # (n_X, p_z+1) 設計行列
    factors_k: NNGPFactors,    # カテゴリkのNNGP因子
    beta_prior_mean: np.ndarray,
    beta_prior_prec: np.ndarray,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns:
        beta_k: (p_z+1,) 固定効果
        u_k: (n_X,) 空間効果
    """
    kappa_tilde = counts_k - 0.5 * total_counts
    z_k = kappa_tilde / xi_k
    # ... 既存の update_beta_category と同様の処理
```

---

## 3. モジュール設計

### 3.1 ディレクトリ構造

```
bayesian_statistics/nngp/model/marked_point_process/
├── __init__.py
├── config.py           # 設定クラス
├── dataset.py          # データセットクラス
├── intensity.py        # 点過程コンポーネント
├── composition.py      # マークコンポーネント（既存ラッパー）
├── sampler.py          # 統合ギブスサンプラー
└── results.py          # 結果クラス・予測
```

### 3.2 各モジュールの責務

#### `config.py`

```python
@dataclass
class MarkedPointProcessConfig:
    """マーク付き点過程モデルの統合設定"""

    # MCMC設定
    n_iter: int = 5000
    burn_in: int = 1000
    thinning: int = 1
    seed: Optional[int] = None

    # NNGP設定（点過程・マーク共通）
    neighbor_count: int = 25

    # 点過程側カーネル
    intensity_kernel_lengthscale: float = 0.05
    intensity_kernel_variance: float = 1.0

    # マーク側カーネル
    mark_kernel_lengthscale: float = 0.05
    mark_kernel_variance: float = 1.0

    # λ* 事前分布
    lambda_prior_shape: float = 2.0
    lambda_prior_rate: float = 1.0

    # 距離事前（マーク側）
    use_distance_prior: bool = True
    distance_tau: float = 1.0
    distance_alpha: float = 1.0

    # 領域設定
    region: List[List[float]] = None  # [[x_min, x_max], [y_min, y_max]]
```

#### `dataset.py`

```python
@dataclass
class MarkedPointProcessDataset:
    """マーク付き点過程の統合データセット"""

    # 時期・産地情報
    period: int
    origins: List[str]  # K 個のカテゴリ名（基準含む）

    # 遺跡データ（点過程 X）
    site_coords: np.ndarray          # (n_sites, 2) 座標
    site_ids: np.ndarray             # (n_sites,) ID

    # マークデータ（構成比 y）
    counts: np.ndarray               # (n_sites, K) カウント
    total_counts: np.ndarray         # (n_sites,) 総数 N_i

    # 設計行列
    design_matrix_intensity: np.ndarray  # (n_sites, p_int+1) 点過程用
    design_matrix_marks: np.ndarray      # (n_sites, p_z+1) マーク用

    # グリッドデータ（予測用）
    grid_coords: np.ndarray          # (n_grid, 2)
    design_matrix_grid_intensity: np.ndarray  # (n_grid, p_int+1)
    design_matrix_grid_marks: np.ndarray      # (n_grid, p_z+1)
    valid_grids: np.ndarray          # (n_grid,) bool

    # 距離事前（オプション）
    distance_features_sites: Optional[np.ndarray] = None  # (n_sites, K-1)
    distance_features_grid: Optional[np.ndarray] = None   # (n_grid, K-1)
    prior_mean_intercept_sites: Optional[np.ndarray] = None  # (K-1, n_sites)
    prior_mean_intercept_grid: Optional[np.ndarray] = None   # (K-1, n_grid)

    # 領域情報
    region: List[List[float]] = None
    volume: float = 0.0  # 有効領域の面積

    def num_sites(self) -> int: ...
    def num_categories(self) -> int: ...
    def num_grid(self) -> int: ...
```

#### `intensity.py`

```python
class IntensityComponent:
    """点過程の強度部分"""

    def __init__(self, config: MarkedPointProcessConfig):
        self.config = config
        self.kernel = LocalNNGPKernel(
            lengthscale=config.intensity_kernel_lengthscale,
            variance=config.intensity_kernel_variance
        )

    def build_factors(self, coords: np.ndarray) -> NNGPFactors:
        """遺跡座標に対するNNGP因子構築"""
        # nngp.py の build_nngp_factors を利用
        pass

    def sample_pseudo_absence(
        self, lambda_star: float, q: np.ndarray,
        valid_grids: np.ndarray, region: List[List[float]]
    ) -> Tuple[np.ndarray, np.ndarray]:
        """偽不在点 U のサンプリング"""
        # ipp.py の poisson_thinning_U を参考に実装
        # U ~ IPP(λ*(1-q))
        pass

    def compute_q(self, eta_int: np.ndarray) -> np.ndarray:
        """存在確率 q(s,t) = sigmoid(η_int)"""
        return 1 / (1 + np.exp(-eta_int))
```

#### `composition.py`

```python
class CompositionComponent:
    """マーク（構成比）部分"""

    def __init__(self, config: MarkedPointProcessConfig):
        self.config = config
        self.kernel = LocalNNGPKernel(
            lengthscale=config.mark_kernel_lengthscale,
            variance=config.mark_kernel_variance
        )

    def build_factors(self, coords: np.ndarray) -> List[NNGPFactors]:
        """マーク用NNGP因子構築"""
        # 既存の FactorCache ロジックを利用
        pass

    def compute_probs(self, eta: np.ndarray) -> np.ndarray:
        """構成確率 π = softmax(η)"""
        return softmax_with_baseline(eta)
```

#### `sampler.py`

```python
class MarkedPointProcessSampler:
    """統合ギブスサンプラー"""

    def __init__(
        self,
        dataset: MarkedPointProcessDataset,
        config: MarkedPointProcessConfig,
    ):
        self.dataset = dataset
        self.config = config

        # コンポーネント初期化
        self.intensity = IntensityComponent(config)
        self.composition = CompositionComponent(config)

        # NNGP因子構築
        self._build_all_factors()

        # Polya-Gamma サンプラー
        self.pg = PolyaGammaSampler()

    def run(self) -> 'MarkedPointProcessResults':
        """MCMCを実行"""
        # 初期化
        lambda_star = self._init_lambda_star()
        beta_int, u_int = self._init_intensity_params()
        beta_marks, u_marks = self._init_mark_params()

        # サンプル保存用
        samples = self._init_storage()

        for iteration in tqdm(range(self.config.n_iter)):
            # === 点過程部分 ===
            # (a) 偽不在のサンプリング
            eta_int = self._compute_eta_intensity(beta_int, u_int)
            q = self.intensity.compute_q(eta_int)
            U, U_design = self.intensity.sample_pseudo_absence(
                lambda_star, q, self.dataset.valid_grids, self.dataset.region
            )

            # (b) ω_i ~ PG(1, η_int,i) for i in X ∪ U
            omega_int = self._sample_omega_intensity(eta_int, U)

            # (c) (β_int, u_int) の更新
            beta_int, u_int = self._update_intensity_params(omega_int, U, U_design)

            # (d) λ* の更新
            n_total = self.dataset.num_sites() + len(U)
            lambda_star = self._update_lambda_star(n_total)

            # === マーク部分 ===
            # (e) ξ_ik ~ PG(N_i, η_ik)
            eta_marks = self._compute_eta_marks(beta_marks, u_marks)
            xi = self._sample_xi_marks(eta_marks)

            # (f) (β_k, u_k) の更新
            beta_marks, u_marks = self._update_mark_params(xi)

            # サンプル保存
            self._store_samples(iteration, samples, ...)

        return MarkedPointProcessResults(...)

    # === 点過程更新メソッド ===
    def _sample_omega_intensity(self, eta_int, U):
        """PG(1, η_int) からサンプル"""
        # X ∪ U の全点について
        pass

    def _update_intensity_params(self, omega, U, U_design):
        """(β_int, u_int) の正規完全条件付きからサンプル"""
        # Sec5.1 の導出に従う
        # 既存の update_beta_category を参考に
        pass

    def _update_lambda_star(self, n_total):
        """λ* ~ Gamma(m0 + n, r0 + |D|)"""
        shape = self.config.lambda_prior_shape + n_total
        rate = self.config.lambda_prior_rate + self.dataset.volume
        return np.random.gamma(shape, 1/rate)

    # === マーク更新メソッド ===
    def _sample_xi_marks(self, eta_marks):
        """PG(N_i, η_ik) からサンプル"""
        # Sec6 の導出に従う
        pass

    def _update_mark_params(self, xi):
        """(β_k, u_k) の更新"""
        # 既存の update_beta_category をそのまま利用可能
        pass
```

#### `results.py`

```python
@dataclass
class MarkedPointProcessResults:
    """MCMC結果"""

    dataset: MarkedPointProcessDataset
    config: MarkedPointProcessConfig

    # 事後サンプル
    lambda_star_samples: np.ndarray  # (n_save,)
    beta_int_samples: np.ndarray     # (n_save, p_int+1, n_sites)
    u_int_samples: np.ndarray        # (n_save, n_sites)
    beta_mark_samples: np.ndarray    # (n_save, K-1, p_z+1, n_sites)
    u_mark_samples: np.ndarray       # (n_save, K-1, n_sites)

    # NNGP因子（予測用）
    intensity_factors: NNGPFactors
    mark_factor_cache: FactorCache

    def predict_site_probability(self) -> np.ndarray:
        """遺跡での存在確率 q(s_i, t)"""
        pass

    def predict_grid_probability(self) -> np.ndarray:
        """グリッドでの存在確率"""
        # NNGP条件付き分布による外挿
        pass

    def predict_site_composition(self) -> np.ndarray:
        """遺跡での産地構成比 π(s_i, t)"""
        pass

    def predict_grid_composition(self) -> np.ndarray:
        """グリッドでの産地構成比"""
        pass

    def decompose_effects(self):
        """効果の分解（距離事前 vs データ駆動）"""
        # distance_prior_model.py の decompose_effects を参考に
        pass
```

---

## 4. 理論文書（sec*.tex）との整合性

### 4.1 修正された数式の反映

以下の修正が理論文書で行われており、実装時はこれに従う：

| セクション | 修正内容 | 実装への影響 |
|-----------|---------|-------------|
| sec4.tex | IPP尤度に階乗 $n_X! n_U!$ **なし** | 尤度計算で階乗を含めない |
| sec6.tex | マークPG: $\xi_{ik} \sim \mathrm{PG}(N_i, \eta_{ik})$, $\tilde{\kappa}_{ik} = y_{ik} - N_i/2$ | 点過程(b=1)とマーク(b=N_i)で異なる |
| sec8.tex | AR(1)精度行列: $Q = Q_t \otimes Q_s$（積の形） | 時間構造導入時はKronecker積 |

### 4.2 PG変数のパラメータ比較

```
┌─────────────────────────────────────────────────────────────────┐
│ 【点過程】sec5.tex                                               │
│   ω_i ~ PG(1, η_int,i)                                          │
│   κ_i = y_i - 1/2  where y_i ∈ {0, 1}                           │
│   z_i = κ_i / ω_i                                                │
├─────────────────────────────────────────────────────────────────┤
│ 【マーク】sec6.tex                                               │
│   ξ_ik ~ PG(N_i, η_ik)                                          │
│   κ̃_ik = y_ik - N_i/2  where y_ik ∈ {0,1,...,N_i}              │
│   z_k,i = κ̃_ik / ξ_ik                                           │
└─────────────────────────────────────────────────────────────────┘
```

### 4.3 完全条件付き分布の共通構造

点過程（sec5.tex）とマーク（sec6.tex）の更新は同じ構造：

```
θ | ω, データ ~ N(m, V)

V^{-1} = H' Ω H + R_0
m = V (H' Ω z + R_0 μ_0)

where:
  H = [W  I_n]      # 設計行列と単位行列の結合
  Ω = diag(ω)       # PG変数の対角行列
  R_0 = [B_0^{-1}  0; 0  Q]  # 事前精度（NNGPスパース）
```

この共通構造により、`update_beta_category`関数を点過程・マーク両方で再利用可能。

---

## 5. 既存コードの再利用

### 5.1 そのまま利用

| モジュール | 利用部分 |
|-----------|---------|
| `nngp.py` | `NNGPFactors`, `build_nngp_factors`, `build_cross_factors`, `FactorCache` |
| `sample.py` | `PolyaGammaSampler`, `LocalNNGPKernel`, `compute_eta`, `softmax_with_baseline`, `update_beta_category` |

### 5.2 参考にして新規実装

| 既存 | 新規 | 内容 |
|-----|-----|-----|
| `ipp.py: poisson_thinning_U` | `intensity.py: sample_pseudo_absence` | 偽不在サンプリング |
| `ipp.py: mcmc_sampler` | `sampler.py: run` | 点過程部分のループ構造 |
| `multinomial_model.py: run_mcmc` | `sampler.py: run` | マーク部分のループ構造 |
| `distance_prior_model.py: prepare_*` | `dataset.py` | 距離事前の準備 |

### 5.3 拡張が必要

| モジュール | 拡張内容 |
|-----------|---------|
| `sample.py: update_beta_category` | `u_int` 単独の更新にも対応（現在は β と u を一緒に更新） |

---

## 6. 実装の順序

### Phase 1: 基盤（1-2日）
1. `config.py` - 設定クラス
2. `dataset.py` - データセットクラス
3. 既存 `ObsidianDataPreprocessor` との連携関数

### Phase 2: コンポーネント（2-3日）
4. `intensity.py` - 点過程部分
   - NNGP因子構築
   - 偽不在サンプリング
   - q(s,t) 計算
5. `composition.py` - マーク部分
   - 既存ロジックのラッパー

### Phase 3: サンプラー（3-4日）
6. `sampler.py` - 統合ギブスサンプラー
   - 点過程更新ステップ
   - マーク更新ステップ
   - 交互実行ループ

### Phase 4: 結果・検証（2-3日）
7. `results.py` - 結果クラス
8. 予測メソッド
9. 単体テスト
10. ノートブックでの検証

## 7. 注意点

### 7.1 座標系
- 全て度単位で統一
- preprocessor から取得する座標は度単位であることを確認

### 7.2 点過程とマークの評価点（重要）
- **点過程**: X ∪ U の全点で η_int を評価（n = n_X + n_U 点）
- **マーク**: X のみで η_k を評価（n_X 点、U にはマーク y がない）
- NNGP因子は点過程用（n点）とマーク用（n_X点）で**別々に構築**

### 7.3 NNGP因子の構築タイミング
- **マーク側**: 遺跡座標 X は固定なので、MCMCループ開始前に1回だけ構築
- **点過程側**: X ∪ U の座標は反復ごとに変化するため、毎反復で再構築が必要
  - これが計算のボトルネックになる可能性あり

### 7.4 計算量の考慮
- 偽不在 U のサイズは反復ごとに変化
- U が大きくなると点過程部分の計算が重くなる
- 対策案:
  - U のサイズに上限を設ける
  - グリッド解像度を下げる
  - NNGP近傍数を減らす

---

## 8. テスト駆動開発（TDD）計画

### 8.1 テストファイル構成

```
tests/nngp/model/marked_point_process/
├── __init__.py
├── conftest.py              # 共通フィクスチャ
├── test_config.py           # 設定クラスのテスト
├── test_dataset.py          # データセットクラスのテスト
├── test_intensity.py        # 点過程コンポーネントのテスト
├── test_sampler.py          # サンプラーのテスト
├── test_consistency.py      # 既存モデルとの一致テスト【最重要】
└── test_integration.py      # 統合テスト
```

### 8.2 単体テスト詳細

#### `test_intensity.py`
```python
def test_sample_pseudo_absence_returns_valid_coords():
    """偽不在点が有効なグリッド内に収まることを確認"""

def test_sample_pseudo_absence_respects_thinning():
    """1-q(s,t)の確率でthinningされることを確認"""

def test_compute_q_is_sigmoid():
    """q(s,t) = sigmoid(η)が正しく計算されることを確認"""
```

#### `test_sampler.py`
```python
def test_sample_omega_intensity_uses_b_equals_1():
    """点過程用PG: b=1でサンプルされることを確認"""

def test_sample_xi_marks_uses_b_equals_N_i():
    """マーク用PG: b=N_iでサンプルされることを確認"""

def test_update_lambda_star_posterior_params():
    """λ* ~ Gamma(m0+n, r0+|D|)のパラメータが正しいことを確認"""

def test_kappa_intensity_formula():
    """κ_i = y_i - 1/2 が正しく計算されることを確認"""

def test_kappa_tilde_marks_formula():
    """κ̃_ik = y_ik - N_i/2 が正しく計算されることを確認"""
```

### 8.3 回帰テスト（既存モデルとの一致）【最重要】

#### `test_consistency.py`
```python
def test_mark_update_matches_multinomial_model():
    """
    【最重要テスト】
    点過程部分を無効化（U=空, λ*=固定）して、
    マーク部分だけ更新した場合、
    既存のmultinomial_model.run_mcmc()と同じ結果になることを確認。

    手順:
    1. 同じデータセットを準備
    2. 同じシード、同じ設定で両モデルを実行
    3. β_samples が数値誤差内で一致することを確認
    """

def test_point_process_nngp_factors_correct():
    """
    X∪U座標でNNGP因子を構築した場合、
    既存のbuild_nngp_factorsと同じ結果になることを確認
    """
```

### 8.4 統合テスト

#### `test_integration.py`
```python
def test_full_mcmc_runs_without_error():
    """小さいデータ（10遺跡、3カテゴリ）でMCMCが完走"""

def test_posterior_samples_are_finite():
    """全事後サンプルにNaN/Infがないことを確認"""

def test_posterior_probabilities_sum_to_one():
    """各地点で Σπ_k = 1 を確認"""

def test_site_probability_in_valid_range():
    """存在確率 q ∈ [0, 1] を確認"""
```

### 8.5 共通フィクスチャ（conftest.py）

```python
@pytest.fixture
def small_dataset():
    """テスト用の小規模データセット"""
    # 10遺跡、3カテゴリ、100グリッド点

@pytest.fixture
def preprocessor_fixture():
    """実データのpreprocessor（統合テスト用）"""

@pytest.fixture
def rng():
    """再現可能な乱数生成器"""
    return np.random.default_rng(42)
```

---

## 9. TDD実装フェーズ

### Phase 1: データ構造（テスト先行）
1. `test_config.py` を書く
2. `config.py` を実装してテストを通す
3. `test_dataset.py` を書く
4. `dataset.py` を実装してテストを通す

### Phase 2: 点過程コンポーネント（テスト先行）
5. `test_intensity.py` を書く
6. `intensity.py` を実装してテストを通す

### Phase 3: マークコンポーネント（回帰テスト重視）
7. `test_consistency.py` の `test_mark_update_matches_multinomial_model` を書く
8. `composition.py` を実装（既存コードのラッパー）
9. 回帰テストが通ることを確認

### Phase 4: 統合サンプラー
10. `test_sampler.py` を書く
11. `sampler.py` を実装
12. `test_integration.py` で全体動作確認

### Phase 5: 結果クラスと予測
13. `results.py` を実装
14. ノートブックでの動作確認

---

*作成日: 2024-12-19*
*更新日: 2024-12-19（テスト戦略追加、理論文書修正の反映）*
*参照: CODEBASE_ANALYSIS.md, docs/sections/sec2-9.tex*
