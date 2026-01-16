# Møller et al. (1998) - Log Gaussian Cox過程

## 基本情報
- **タイトル**: Log Gaussian Cox Processes
- **著者**: Jesper Møller, Anne Randi Syversveen, Rasmus Plenge Waagepetersen
- **所属**: Aalborg University (Denmark), Norwegian University of Science and Technology, University of Aarhus (Denmark)
- **掲載**: Scandinavian Journal of Statistics
- **年**: 1998
- **DOI**: 記載なし
- **関連論点**: 論点1（空間点過程）、論点3（クラスタリング）

## 主な貢献
- **Log Gaussian Cox Process (LGCP)の体系的な理論**: 1次・2次・3次モーメント特性の完全な導出
- **Pair correlation functionとの直接的関係**: $g(s) = \exp\{\sigma^2 r(s)\}$ （定理1）
- **Bayesian推論の枠組み**: Empirical Bayesian methodsによる強度曲面の予測（Section 8）
- **シミュレーション手法**: 辺縁効果なしの正確なシミュレーション（Section 6）
- **多変量LGCPへの拡張**: Multivariate log Gaussian Cox processesの定義（Section 5）

## 手法の概要

### Cox過程の定義

**Cox過程** (doubly stochastic Poisson process):
- ランダムな強度測度$\Lambda$で決定される非斉次Poisson過程
- 条件付き分布: $X | \Lambda \sim \text{Poisson}(\Lambda)$

$$
\text{card}(X \cap B) | \Lambda \sim \text{Poisson}\left(\int_B \Lambda(s) ds\right)
$$

for bounded Borel sets $B \subset \mathbb{R}^2$

**特徴**:
- Doubly stochastic（二重確率性）
- 集積（aggregation）を確率的環境の不均質性でモデル化
- クラスタリングの自然なモデル

### Log Gaussian Cox Process (LGCP)

**定義**: 強度過程を対数ガウス過程でモデル化：

$$
\Lambda(s) = \exp\{Y(s)\}
$$

ここで$Y = \{Y(s): s \in \mathbb{R}^2\}$はガウス過程

**定常LGCPのパラメータ化**:
- 平均: $\mu = \mathrm{E}[Y(s)]$
- 分散: $\sigma^2 = \text{Var}[Y(s)]$
- 相関関数: $r(s_1 - s_2) = \text{Cov}(Y(s_1), Y(s_2)) / \sigma^2$

**存在条件**:
1. **Positive semi-definiteness**: $\sum_{i,j} a_i a_j r(s_i - s_j) \geq 0$ for all $a_i \in \mathbb{R}$, $s_i \in \mathbb{R}^2$
2. **Continuity**: ほぼ確実に連続な修正の存在

$$
1 - r(s) < \frac{K}{(-\log(\|s\|))^{1+\varepsilon}} \quad \text{for } \|s\| < 1
$$

または（より強い条件）$1 - r(s) < \tilde{K} |s|^\alpha$ for some $\tilde{K} > 0$, $\alpha > 0$

### パラメータの解釈

**スケールとシェイプパラメータ**:

$$
\Lambda_{\mu, \sigma} \overset{\mathscr{D}}{=} e^\mu \Lambda_{0,1}^\sigma
$$

ここで$\ln \Lambda_{\mu, \sigma}$は平均$\mu$、分散$\sigma^2$、相関関数$r(\cdot)$の定常ガウス過程

- $e^\mu > 0$: スケールパラメータ
- $\sigma > 0$: シェイプパラメータ（クラスタリングの程度）

**極端ケース**:
- $\sigma \to 0$: 斉次Poisson過程（クラスタリングなし）
- $r(\cdot) = 1$: Mixed Poisson過程（$\Lambda(\cdot) = \lambda$が定数でlog-Gaussian分布）

### 等方的相関モデル（Table 1）

| モデル | 相関関数$r(a)$ | パラメータ |
|--------|---------------|------------|
| 1. Gaussian | $\exp(-(a/\beta)^2)$ | スケール$\beta$ |
| 2. Exponential | $\exp(-a/\beta)$ | スケール$\beta$ |
| 3. Cardinal sine | $\sin(a/\beta) / (a/\beta)$ | スケール$\beta$ |
| 4. Stable | $\exp(-\sqrt{a/\beta})$ | スケール$\beta$ |

ここで$a = \|s\|$ (distance)

**性質**:
- 全て$a \to \infty$で$r(a) \to 0$
- 全て連続修正の存在条件を満たす
- パラメータ$\beta$: 空間的スケールを制御

## 理論的結果

### 定理1: 積密度と対相関関数

**$n$次積密度**:

$$
\rho^{(n)}(s_1, \ldots, s_n) = \mathrm{E}\left[\prod_{i=1}^n \Lambda(s_i)\right]
$$

**定常LGCPの場合**:

$$
\rho^{(n)}(s_1, \ldots, s_n) = \exp\left\{n\mu + \sigma^2\left[\frac{n}{2} + \sum_{1 \leq i < j \leq n} r(s_i - s_j)\right]\right\}
$$

$$
= \rho^n \prod_{1 \leq i < j \leq n} g(s_i - s_j)
$$

ここで：

**強度** (intensity):

$$
\rho = \rho^{(1)}(s) = \exp\{\mu + \sigma^2/2\}
$$

**対相関関数** (pair correlation function):

$$
g(s_1 - s_2) = \frac{\rho^{(2)}(s_1, s_2)}{\rho^2} = \exp\{\sigma^2 r(s_1 - s_2)\}
$$

**重要な関係**:
- $g(0) = \exp\{\sigma^2\}$ （完全な重なり）
- $r(s) \to 0$ as $s \to \infty$ ⇒ $g(s) \to 1$（独立性）
- $g(s) > 1$ ⇔ クラスタリング（正の相関）

### 1次・2次統計量

**期待点数** (単位正方形$[0,1]^2$内):

$$
\mathrm{E}[\text{card}(X \cap [0,1]^2)] = \rho
$$

**分散**:

$$
\text{Var}[\text{card}(X \cap [0,1]^2)] = \rho + \rho^2 \left(\int_{[0,1]^2} \int_{[0,1]^2} \exp\{\sigma^2 r(s-t)\} ds dt - 1\right)
$$

**性質**:
- 分散 ≥ 平均（over-dispersion）
- $\beta$（相関範囲）↑ ⇒ 分散↑

### 3次統計量

**3次密度**:

$$
\rho^{(3)}(x_1, x_2, x_3) = \rho^3 g(x_1 - x_2) g(x_2 - x_3) g(x_1 - x_3)
$$

**3次関数$z$の定義**:

$$
z(u, v) = \frac{\rho^{(3)}(0, u, v)}{\rho^{(3)}(0, u, 0) \cdot \rho^{(3)}(0, 0, v)} = g(u-v)
$$

**性質**: $z(u, v)$は$g(\cdot)$のみで決定される（定理2）

**応用**: モデル診断（Section 7）
- 推定量$\hat{z}(u, v)$と理論値$g(u-v)$の比較
- 3次モーメントに基づくモデル検証

### 定理3: エルゴード性

**定理**: 等方的な定常LGCPは$r(s) \to 0$ as $\|s\| \to \infty$ならば**ergodic**

**意味**: サンプル平均 → 期待値（大数の法則）

$$
\frac{1}{|W|} \int_W h(X \cap (W + s)) ds \to \mathrm{E}[h(X)] \quad \text{as } W \nearrow \mathbb{R}^2
$$

**応用**: パラメータの一致推定量の構成

## Neyman-Scott過程との比較（Section 4）

**Neyman-Scott過程**: 親点過程 + 子点の散布

**Thomas過程**（特殊ケース）:
- 親点: Poisson過程、強度$\kappa$
- 子点: 各親の周りにGaussian散布、標準偏差$\omega$

**対相関関数** (Thomas):

$$
g_{\text{Thomas}}(s) = 1 + \frac{1}{4\pi\kappa\omega^2} \exp\left\{-\frac{\|s\|^2}{4\omega^2}\right\}
$$

**LGCPとの関係** (Figure 1(e), (f)):
- Thomas過程 ≈ LGCP（Gaussian相関関数）
- パラメータの対応関係で近似的に一致
- **違い**: LGCPは理論的に扱いやすい（モーメント、推論）

## 多変量LGCP（Section 5）

**定義**: $m$種類の点過程$X_1, \ldots, X_m$

$$
\Lambda_k(s) = \exp\{Y_k(s)\}, \quad k=1, \ldots, m
$$

ここで$(Y_1(s), \ldots, Y_m(s))$は多変量ガウス過程

**Cross pair correlation function**:

$$
g_{jk}(s) = \exp\{\text{Cov}(Y_j(0), Y_k(s))\}
$$

**応用**: 複数種の生物の空間分布、異なる産地の遺跡分布など

## シミュレーション手法（Section 6）

### 方法1: Circulant embedding

**ステップ**:
1. ガウス場$Y$を周期的拡張してFFTで高速シミュレーション
2. $\Lambda(s) = \exp\{Y(s)\}$を計算
3. 条件付きPoisson過程をシミュレーション

**利点**: 高速、正確、辺縁効果なし

### 方法2: 条件付きシミュレーション

**条件**: 点数$n = \text{card}(X \cap W)$を固定

**手順**: Accept-reject法により$n$を条件とするLGCPをサンプル

**応用**: Figure 3の全シミュレーション（$n=148$固定）

## パラメータ推定とモデル診断（Section 7）

### 推定方法

**1. Minimum contrast estimation**:

$$
(\hat{\mu}, \hat{\sigma}, \hat{\beta}) = \arg\min_{\mu, \sigma, \beta} \int_0^{t^*} (\hat{K}(t) - K(t; \mu, \sigma, \beta))^2 dt
$$

ここで$K(t)$はRipley's K関数、$\hat{K}(t)$はその推定量

**2. Maximum pseudo-likelihood**:

ガウス近似を用いた疑似尤度最大化

### モデル診断

**1. 対相関関数の比較**:
- 推定$\hat{g}(s)$ vs 理論$g(s; \hat{\theta})$

**2. 3次統計量の比較**:
- 推定$\hat{z}(u, v)$ vs 理論$g(u-v; \hat{\theta})$

**3. Residual analysis**:
- Transformed residuals のPoisson性の検定

### 実データ例

**Example 1: Japanese black pines**（79点）:
- Exponential相関モデルが適合
- $\hat{\mu} = 2.70$, $\hat{\sigma} = 1.36$, $\hat{\beta} = 0.027$

**Example 2: Tropical rain forest** (bivariate, 種A: 65点, 種B: 73点):
- 2種間の正の相関を検出
- Cross pair correlation function によるモデル化

## Empirical Bayesian推論（Section 8）

### 目的

観測されたLGCP $X$から未観察のガウス場$Y$と強度曲面$\Lambda$を予測

### ベイズ的枠組み

**事後分布**:

$$
p(Y | X) \propto p(X | Y) p(Y)
$$

ここで：
- $p(X | Y)$: Poisson尤度
- $p(Y)$: ガウス過程事前分布

**Empirical Bayesアプローチ**:
1. データから$(\mu, \sigma, \beta)$を推定
2. 推定値を事前分布に代入
3. MCMCで事後分布をサンプル

### Metropolis-adjusted Langevin algorithm (MALA)

**提案分布**:

$$
Y^* = Y + \frac{\epsilon^2}{2} \nabla \log p(Y | X) + \epsilon Z
$$

ここで$Z \sim N(0, I)$、$\nabla \log p(Y | X)$は対数事後密度の勾配

**利点**:
- Gradient情報を利用 → 効率的
- Gibbs samplingより実装が容易
- 辺縁効果の問題なし

## 実証結果

### シミュレーション研究（Figure 2, 3）

**Figure 2**: ガウス場$Y$のシミュレーション（4種類の相関関数）
- Gaussian: 滑らかな変動
- Exponential: 中程度の滑らかさ
- Cardinal sine: 周期的構造
- Stable: 小スケールで変動

**Figure 3**: 対応するLGCP（条件付き$n=148$）
- **Row 1**: $\sigma=1$, 大$\beta$ → 大きく疎なクラスタ
- **Row 2**: $\sigma=2.4$, 中$\beta$ → 中程度の大きさのクラスタ
- **Row 3**: $\sigma=2.4$, 小$\beta$ → 多数の小クラスタ

**視覚的比較**:
- Gaussian と Cardinal sine: 類似パターン
- Stable: 他よりクラスタリングが弱い

### 実データ解析結果

**Example 1 (Japanese black pines)**:
- $\hat{K}(t)$とExponential LGCP理論値が良好に一致
- $\hat{z}(u, v)$もモデルを支持

**Example 2 (Tropical rain forest)**:
- 2種間の空間的相互作用を定量化
- Bivariate LGCPが適合

## 限界・残された課題

### 理論的限界

1. **相関関数の選択**:
   - Table 1の7モデルのどれを選ぶか?
   - **課題**: モデル選択基準（AIC, BIC）の導出

2. **非定常性**:
   - 本論文は定常LGCPに焦点
   - **課題**: 非定常ガウス場への拡張（トレンド + 定常成分）

3. **Inhibition（反発）の扱い**:
   - LGCPは$g(s) \geq 1$ → クラスタリングのみ
   - **課題**: $g(s) < 1$となる反発過程への拡張

4. **高次モーメント**:
   - 3次まで導出済み
   - **課題**: 4次以上のモーメント特性、Cumulantsの利用

### 実装上の課題

5. **計算コスト**:
   - Circulant embedding: FFT $O(n \log n)$だが、大規模データで重い
   - MALA: 各イテレーションで勾配計算が必要
   - **課題**: 数千〜数万点のデータへのスケーラビリティ

6. **パラメータ推定の不安定性**:
   - Minimum contrast推定は初期値に敏感
   - **課題**: Robust推定法、ベイズ推定への完全移行

7. **辺縁効果**:
   - 論文は「辺縁効果なし」と主張
   - **現実**: 有界窓$W$での観測 → edge correction必要
   - **課題**: Ripley's edge correctionとの整合性

### 応用上の課題

8. **ゼロ過剰データ**:
   - 実データでポイントが全くない領域が多い場合
   - **課題**: Zero-inflated LGCP

9. **Marked point processへの拡張**:
   - 本論文は位置のみ
   - **課題**: 連続マークを持つLGCP

10. **時空間LGCP**:
    - 本論文は空間のみ
    - **課題**: $Y(s, t)$の時空間ガウス過程

11. **外れ値**:
    - 極端なクラスタ（例: 数百点の密集）
    - **課題**: Robust LGCP

## 本研究（MMCP）との関係

### 借用すべき手法

1. **LGCPの基本的枠組み**:
   - **直接的な適用**: MMCPの点過程部分
   - $\Lambda(s) = \exp\{Y(s)\}$
   - 黒曜石遺跡の分布をLGCPでモデル化

2. **Pair correlation functionの関係**:
   - $g(s) = \exp\{\sigma^2 r(s)\}$
   - 相関関数$r(\cdot)$の選択 → クラスタリングパターンの制御
   - **応用**: トブラー距離に基づく相関関数

3. **理論的性質の活用**:
   - 1次・2次・3次モーメントの明示的公式
   - モデル診断統計量（$\hat{K}(t)$, $\hat{z}(u,v)$）
   - エルゴード性による一致推定

4. **シミュレーション手法**:
   - Circulant embeddingによる高速シミュレーション
   - 条件付きシミュレーション（点数固定）
   - **応用**: 予測分布のサンプリング

5. **Bayesian推論**:
   - Empirical Bayes + MALA
   - 強度曲面$\Lambda(s)$の事後分布
   - **応用**: 遺跡密度の空間的予測

6. **多変量LGCP**:
   - 複数の産地別遺跡分布を同時モデル化
   - Cross pair correlation function $g_{jk}(s)$
   - **応用**: 産地間の空間的相互作用

### 拡張すべき点

1. **マーク付き点過程への統合**:
   - Møller et al. (1998): 位置のみ
   - MMCP: 位置 + 組成マーク
   - **課題**: $\Lambda(s, \mathbf{y})$の共同モデリング

2. **組成データのマーク分布**:
   - 各遺跡$s_i$に組成マーク$\mathbf{y}_i \in \mathbb{S}^{d-1}$
   - 強度関数が組成に依存: $\Lambda(s, \mathbf{y}) = \Lambda_1(s) f(\mathbf{y} | s)$
   - **課題**: $f(\mathbf{y} | s)$のLogistic-normalモデル

3. **空間的に変化するマーク分布**:
   - 産地からの距離で組成分布が変化
   - $f(\mathbf{y} | s; \boldsymbol{\mu}(s), \boldsymbol{\Sigma}(s))$
   - **課題**: Spatially varying coefficients (SVC)との統合

4. **トブラー距離の組み込み**:
   - ユークリッド距離ではなく地形考慮距離
   - 相関関数$r(d_{\text{Tobler}}(s_1, s_2))$
   - **課題**: 非ユークリッド距離に対するガウス過程

5. **Preferential samplingへの拡張**:
   - 遺跡位置と組成が依存（産地近くに高濃度）
   - Shared latent process $Y(s)$
   - **課題**: Gelfand & Shirota (2019)との統合

6. **階層モデルへの拡張**:
   - 時代レベル（縄文早期〜後期）
   - 産地レベル
   - **課題**: Hierarchical LGCP with compositional marks

7. **ベイズ推定の完全実装**:
   - Empirical Bayesではなくフルベイズ
   - $(\mu, \sigma, \beta)$の事前分布
   - **課題**: MALAの収束診断、事後予測チェック

### MMCPにおける位置づけ

**第3章での役割**:
- Cox過程の基礎理論としてLGCPを紹介
- Clusteringの確率的モデリングの主要アプローチ
- 対相関関数の理論的導出（定理1）の引用

**実装の基礎**:
- `bayesian_statistics/`モジュールでのLGCPシミュレーション
- Pair correlation functionの推定と可視化
- MALA for posterior samplingの実装

**方法論的貢献**:
- Doubly stochasticの概念 → 2段階モデリング
- ガウス場の対数変換 → 非負性の自然な保証
- Empirical Bayes → 計算効率とフレキシビリティのバランス

**理論的保証**:
- エルゴード性 → 大標本での推定の一致性
- 積密度の明示的公式 → モデル診断の理論的基盤
- 3次統計量 → 高次の依存構造の検証

## 引用すべき箇所

### Cox過程の定義と動機

> "Cox processes provide useful and frequently applied models for aggregated spatial point patterns where the aggregation is due to a stochastic environmental heterogeneity." (Section 1)

→ Cox過程の応用上の重要性

> "A Cox process is 'doubly stochastic' as it arises as an inhomogeneous Poisson process with a random intensity measure." (Section 1)

→ Doubly stochasticの定義

### LGCPの優位性

> "We show that the class of stationary log Gaussian Cox processes possesses various appealing properties. (i) The distribution is completely characterized by the intensity and the pair correlation function of the Cox process. This makes parametric models easy to interpret and simple methods are available for parameter estimation and model checking." (Section 1)

→ LGCPの解釈の容易さ

> "(ii) Theoretical properties are easily derived. Higher-order properties are for instance simply expressed by the intensity and the pair correlation function of the log Gaussian Cox process." (Section 1)

→ 理論的扱いやすさ

### 定理1の中核

> "$\rho^{(n)}(s_1, \ldots, s_n) = \rho^n \prod_{1 \leq i < j \leq n} g(s_i - s_j)$ where [...] $g(s_1 - s_2) = \rho^{(2)}(s_1, s_2) / \rho^2 = \exp\{\sigma^2 r(s_1 - s_2)\}$ are the intensity and the pair correlation function of the process, respectively." (Theorem 1)

→ 積密度と対相関関数の明示的公式

### パラメータの解釈

> "The parameters $e^\mu > 0$ and $\sigma > 0$ have a clear interpretation as a scale and shape parameter, respectively." (Section 2)

→ $\mu$と$\sigma$の解釈

### Rathbun & Cressie (1994)との違い

> "The log Gaussian Cox processes studied in the present paper are in contrast to those [Rathbun & Cressie 1994] specified by such characteristics, and discretized versions of our log Gaussian Cox processes can be simulated exactly without any problem with edge effects." (Section 1)

→ 本論文の優位性（辺縁効果なし）

### エルゴード性

> "[Theorem 3] An isotropic stationary log Gaussian Cox process is ergodic if $r(s) \to 0$ as $\|s\| \to \infty$." (Theorem 3)

→ 大数の法則の成立

### 3次統計量の応用

> "This [third-order property] can be used to check our model assumptions as demonstrated in section 7." (After Theorem 2)

→ モデル診断への利用

### MALAの利点

> "Also the Metropolis-adjusted Langevin algorithm (Besag, 1994; Roberts & Tweedie, 1997) for simulating from the posterior of the intensity surface as studied in section 8 is both easy to specify and implement." (Section 1)

→ Bayesian推論の実装の容易さ
