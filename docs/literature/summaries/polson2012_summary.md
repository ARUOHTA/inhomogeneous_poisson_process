# Polson et al. (2013) - Pólya-Gamma Augmentation

## 基本情報
- **タイトル**: Bayesian inference for logistic models using Pólya-Gamma latent variables
- **著者**: Nicholas G. Polson, James G. Scott, Jesse Windle
- **所属**: University of Chicago (Polson), University of Texas at Austin (Scott, Windle)
- **年**: 2013（First Draft: August 2011, This Draft: July 2013）
- **ジャーナル**: Journal of the American Statistical Association
- **関連論点**: 論点2（組成データ - 技術的基盤）

## 主な貢献

1. **Pólya-Gamma分布の導入**
   - Gamma分布の無限畳み込みとして定義される新しい分布族
   - Laplace変換: $\mathbb{E}\{\exp(-\omega t)\} = \cosh^{-b}(\sqrt{t/2})$
   - 指数tilting: $PG(b, c)$クラスの構成

2. **ロジスティック回帰の厳密なデータ拡張**
   - Albert & Chib (1993)のprobitモデルと同様の単純さ
   - 1層の潜在変数のみ（従来手法は2層以上）
   - 厳密推論（近似なし）

3. **Theorem 1: 基本的な積分恒等式**
   - 二項尤度とPólya-Gamma分布の関係を確立
   - Gaussianの尺度混合（scale mixture）表現
   - Gibbsサンプラーの理論的基盤

4. **効率的なサンプリング手法**
   - Alternating-series method（Devroye 1986）の拡張
   - Accept-reject確率: 0.99919以上（$PG(1, c)$の場合）
   - 自動調整不要（チューニングフリー）

5. **Rパッケージ BayesLogit の開発**
   - 高速・正確なPólya-Gammaサンプラー
   - 様々なモデルへの拡張（混合モデル、空間モデル等）

## 手法の概要

### Pólya-Gamma分布の定義

**Definition 1**: 確率変数$X$がパラメータ$b > 0$, $c \in \mathbb{R}$のPólya-Gamma分布に従う、$X \sim PG(b, c)$、とは：
$$
X \stackrel{D}{=} \frac{1}{2\pi^2} \sum_{k=1}^\infty \frac{g_k}{(k-1/2)^2 + c^2/(4\pi^2)}
$$
ここで$g_k \sim \text{Ga}(b, 1)$は独立なGamma分布。

**Laplace変換**（$PG(b, 0)$の場合）:
$$
\mathbb{E}\{\exp(-\omega t)\} = \frac{1}{\cosh^b(\sqrt{t/2})}
$$

**一般的な$PG(b, c)$**（指数tilting）:
$$
p(\omega | b, c) = \frac{\exp(-\frac{c^2}{2}\omega) p(\omega | b, 0)}{\mathbb{E}_\omega\{\exp(-\frac{c^2}{2}\omega)\}} = \frac{\exp(-\frac{c^2}{2}\omega) p(\omega | b, 0)}{\cosh^{-b}(c/2)}
$$

**Laplace変換**（$PG(b, c)$の場合）:
$$
\mathbb{E}_\omega\{\exp(-\omega t)\} = \frac{\cosh^b(c/2)}{\cosh^b(\sqrt{(c^2/2 + t)/2})}
$$

### Theorem 1: 基本的な積分恒等式

**核心的な結果**: $p(\omega)$を$\omega \sim PG(b, 0)$の密度とするとき、全ての$a \in \mathbb{R}$に対して：
$$
\frac{(e^\psi)^a}{(1+e^\psi)^b} = 2^{-b} e^{\kappa\psi} \int_0^\infty e^{-\omega\psi^2/2} p(\omega) d\omega
$$
ここで$\kappa = a - b/2$。

**重要な帰結**: 条件付き分布
$$
p(\omega | \psi) = \frac{e^{-\omega\psi^2/2} p(\omega)}{\int_0^\infty e^{-\omega\psi^2/2} p(\omega) d\omega}
$$
もPólya-Gamma分布: $(\omega | \psi) \sim PG(b, \psi)$。

**証明の要点**:
$$
\frac{(e^\psi)^a}{(1+e^\psi)^b} = \frac{2^{-b} \exp\{\kappa\psi\}}{\cosh^b(\psi/2)} = 2^{-b} e^{\kappa\psi} \mathbb{E}_\omega\{\exp(-\omega\psi^2/2)\}
$$

### ロジスティック回帰へのデータ拡張

**モデル**:
- $y_i \sim \text{Binom}(n_i, 1/(1+e^{-\psi_i}))$
- $\psi_i = x_i^T\beta$（対数オッズ）
- $\beta \sim N(b, B)$（Gaussian事前分布）

**Gibbsサンプラー**（2ステップ）:

1. **潜在変数のサンプリング**:
$$
(\omega_i | \beta) \sim PG(n_i, x_i^T\beta)
$$

2. **回帰係数のサンプリング**:
$$
(\beta | y, \omega) \sim N(m_\omega, V_\omega)
$$
ここで：
$$
\begin{aligned}
V_\omega &= (X^T\Omega X + B^{-1})^{-1} \\
m_\omega &= V_\omega(X^T\kappa + B^{-1}b)
\end{aligned}
$$
- $\kappa = (y_1 - n_1/2, \ldots, y_N - n_N/2)$
- $\Omega = \text{diag}(\omega_1, \ldots, \omega_N)$

**尤度の表現**（観測$i$の寄与）:
$$
L_i(\beta) = \frac{\{\exp(x_i^T\beta)\}^{y_i}}{1+\exp(x_i^T\beta)} \propto \exp(\kappa_i x_i^T\beta) \int_0^\infty \exp\{-\omega_i(x_i^T\beta)^2/2\} p(\omega_i | n_i, 0) d\omega_i
$$

**条件付き事後分布**（$\beta$の完全条件付き分布）:
$$
\begin{aligned}
p(\beta | \omega, y) &\propto p(\beta) \prod_{i=1}^N \exp\{\kappa_i x_i^T\beta - \omega_i(x_i^T\beta)^2/2\} \\
&\propto p(\beta) \exp\left\{-\frac{1}{2}(z - X\beta)^T\Omega(z - X\beta)\right\}
\end{aligned}
$$
ここで$z = (\kappa_1/\omega_1, \ldots, \kappa_N/\omega_N)$（working responses）。

### Pólya-Gammaのサンプリング手法

**Alternating-series method**（Devroye 1986の拡張）:

**基本アイデア**: 密度を交互符号級数で表現：
$$
f(x) = \sum_{n=0}^\infty (-1)^n a_n(x)
$$
ここで$a_n(x)$は減少数列（固定$x$に対して）。

**部分和**:
$$
S_n(x) = \sum_{i=0}^n (-1)^i a_i(x)
$$
は次を満たす：
$$
S_0(x) > S_2(x) > \cdots > f(x) > \cdots > S_3(x) > S_1(x)
$$

**Accept-reject手順**:
1. 提案分布$g(x)$から$X$を生成
2. $U \sim \mathcal{U}(0, cg(X))$を生成（$\|f/g\|_\infty \leq c$）
3. $U \leq S_i(X)$が奇数$i$で成立すればaccept
4. $U > S_i(X)$が偶数$i$で成立すればreject

**PG(1, z)の場合の性質**:
- Accept確率: 0.99919以上（非常に効率的）
- 必要な分布: Exponential, Inverse-Gaussian（標準的）
- チューニング不要（自動的に最適性能）

**PG(n, z)の場合**:
- $n$個の独立な$PG(1, z)$の和として表現
- $PG(b_1, z) + PG(b_2, z) = PG(b_1 + b_2, z)$（畳み込み性）

### Pólya-Gammaの性質

**期待値**（Laplace変換の微分から）:
$$
\mathbb{E}(\omega) = \frac{b}{2c}\tanh(c/2) = \frac{b}{2c}\left(\frac{e^c - 1}{1 + e^c}\right)
$$

**密度の表現**（Inverse-Gaussianの交互和）:
$$
f(x | b, c) = \{\cosh^b(c/2)\} \frac{2^{b-1}}{\Gamma(b)} \sum_{n=0}^\infty (-1)^n \frac{\Gamma(n+b)}{\Gamma(n+1)} \frac{(2n+b)}{\sqrt{2\pi x^3}} e^{-\frac{(2n+b)^2}{8x} - \frac{c^2}{2}x}
$$

**全ての有限モーメントが閉形式**:
- EM algorithmでの応用可能
- 完全データ十分統計量として利用

### 従来手法との比較

**Albert & Chib (1993) - Probit**:
- 潜在効用$z_i = x_i^T\beta + \epsilon_i$、$\epsilon_i \sim N(0, 1)$
- 1層の潜在変数（truncated normal）
- 単純で効率的

**Holmes & Held (2006) - Logit**:
- 潜在効用$z_i = x_i^T\beta + \epsilon_i$、$\epsilon_i \sim \text{Logistic}(1)$
- 2層の潜在変数: Kolmogorov-Smirnov分布を追加
- 複雑で効率低い

**Frühwirth-Schnatter & Frühwirth (2010) - Logit**:
- 離散混合正規分布による近似（$K = 10$）
- 2層の潜在変数
- 近似だが実用的

**Pólya-Gamma（本手法）**:
- 1層の潜在変数のみ
- 厳密推論（近似なし）
- 単純かつ効率的

**比較図**（Figure 1）:
- 左: Holmes & Held / Frühwirth-Schnattor（2層構造）
- 右: Pólya-Gamma（1層構造、直接的なデータ拡張）

### 混合モデルへの拡張

**例: バングラデシュ避妊使用データ**（mlmRevパッケージ）

**モデル**:
$$
\begin{aligned}
y_{ij} &\sim \text{Binom}(1, p_{ij}), \quad p_{ij} = \frac{e^{\psi_{ij}}}{1 + e^{\psi_{ij}}} \\
\psi_{ij} &= m + \delta_j + x_{ij}'\beta \\
\delta_j &\sim N(0, 1/\phi) \\
m &\sim N(0, \kappa^2/\phi), \quad \kappa \to \infty
\end{aligned}
$$

**結果**:
- 地区ごとのランダム切片を適切に推定
- 中央ESS = 8,168、中央ESR = 59.88
- 非階層モデル（26.1秒）vs 混合モデル（27.3秒）- ほぼ同じ計算時間

## 実証結果

### 効率性の評価（Section 5）

**主な主張**:
1. **単純なロジットモデル**（データ豊富、階層なし）: Independence Metropolis-Hastingsとほぼ同等（提案分布が適切な場合）
2. **他の全てのケース**: Pólya-Gamma法が最も効率的

**例外**: 負の二項回帰で観測値あたりのカウント数が大きい場合
- ESS（有効サンプルサイズ）は最良
- ESR（有効サンプリングレート）が低下
- 理由: $PG(n, c)$を$n$個の$PG(1, c)$の和で表現するため

### 一様エルゴード性

**Choi & Hobert (2013)の理論結果**:
- Pólya-Gamma Gibbsサンプラーは一様エルゴード
- Monte Carlo平均の中心極限定理が保証される
- 他のMCMC手法には同様の結果なし

## 限界・残された課題

### 著者が述べる限界

1. **大きな形状パラメータの計算コスト**
   - $PG(n, c)$を$n$個の$PG(1, c)$の和で計算
   - 負の二項回帰で$n$が大きい場合にボトルネック
   - Windle et al. (2013b)で改良版を提案

2. **多項分布への拡張**
   - Technical supplementで議論
   - Multinomial logitへの適用は可能だが詳細は別文献

3. **EM algorithmへの適用**
   - 潜在変数$\omega$の期待値が閉形式
   - 完全データ十分統計量として利用可能（未実装）

### 本研究の観点からの限界

1. **組成データ特有の制約**
   - 本論文はbinomial/negative-binomialに焦点
   - 単体制約（sum-to-one）の明示的扱いなし

2. **点過程への適用**
   - ロジスティック回帰やGLMが主な対象
   - 点過程の強度関数への適用は別研究

3. **空間過程との統合**
   - 空間モデルへの言及はあるが詳細は限定的
   - NNGPとの統合は扱っていない

## 本研究（MMCP）との関係

### 借用する要素

1. **Pólya-Gamma augmentationの基本手法**
   - Theorem 1の積分恒等式
   - Gibbsサンプリングの枠組み
   - 尺度混合としてのGaussian表現

2. **効率的なサンプリング手法**
   - Alternating-series method
   - BayesLogit Rパッケージの実装
   - Accept-reject確率の高さ（0.99919以上）

3. **階層モデルへの自然な拡張**
   - 混合モデル、空間モデルへの応用
   - 1層の潜在変数による単純性

### MMCPによる拡張

1. **Multinomial logitへの適用**
   - Polsonはbinomial（2値）が主
   - MMCPは多項選択（組成）への拡張
   - 各組成成分に対するPólya-Gamma拡張

2. **点過程の強度関数への応用**
   - Polsonはロジスティック回帰
   - MMCPは点過程の強度関数モデリング
   - Moreira & Gamerman (2022)がpresence-only点過程にPólya-Gammaを適用

3. **空間過程（NNGP）との統合**
   - PolsonではGPは言及程度
   - MMCPはPólya-Gamma + NNGP + 点過程の統合
   - Moreira et al. (2024)が既に実装

4. **組成データへの適用**
   - 単体制約の下でのPólya-Gamma拡張
   - $D$-部分組成の各成分へのmultinomial logit

### 位置づけ

Polson et al. (2013)は：
- ロジスティック回帰のベイズ推論における画期的な手法
- Pólya-Gamma分布という新しい分布族の導入
- Albert & Chib (1993)のprobitモデルと同等の単純さを実現
- 厳密推論・高効率・拡張性の三拍子

MMCPは：
- この技術的基盤を継承
- Multinomial logit（組成データ）への拡張
- 点過程の強度関数への適用
- NNGPとの統合（計算効率化）
- Presence-onlyデータへの対応

特にMoreira & Gamerman (2022)がPólya-Gammaをpresence-only点過程に初めて適用し、Moreira et al. (2024)がNNGPとの統合を実現。MMCPはこれらに組成マークを追加する形。

## 引用すべき箇所

### データ拡張の新規性（Abstract）
> "We propose a new data-augmentation strategy for fully Bayesian inference in models with binomial likelihoods. The approach appeals to a new class of Pólya-Gamma distributions, which are constructed in detail."

> "In each case, our data-augmentation strategy leads to simple, effective methods for posterior inference that: (1) circumvent the need for analytic approximations, numerical integration, or Metropolis-Hastings; and (2) outperform other known data-augmentation strategies, both in ease of use and in computational efficiency."

### ロジスティック回帰の困難さ（Section 1）
> "Bayesian inference for the logistic regression model has long been recognized as a hard problem, due to the analytically inconvenient form of the model's likelihood function. By comparison, Bayesian inference for the probit model is much easier, owing to the simple latent-variable method of Albert and Chib (1993) for posterior sampling."

### 本手法の単純性（Section 1）
> "In this paper, we present a new data-augmentation algorithm for Bayesian logistic regression. Although our method involves a different missing-data mechanism from that of Albert and Chib (1993), it is nonetheless a direct analogue of their construction, in that it is both exact and simple."

### Theorem 1（核心的結果）（Section 3.1）
> "Let $p(\omega)$ denote the density of the random variable $\omega \sim PG(b, 0), b > 0$. Then the following integral identity holds for all $a \in \mathbb{R}$:
> $$\frac{(e^\psi)^a}{(1+e^\psi)^b} = 2^{-b} e^{\kappa\psi} \int_0^\infty e^{-\omega\psi^2/2} p(\omega) d\omega$$
> where $\kappa = a - b/2$."

### 従来手法との違い（Section 3.2）
> "Our data-augmentation scheme differs from each of these approaches in several ways. First, it does not appeal directly to the random-utility interpretation of the logit model. Instead, it represents the logistic CDF as a mixture with respect to an infinite convolution of gammas. Second, the method is exact, in the sense of making draws from the correct joint posterior distribution, rather than an approximation to the posterior that arises out of an approximation to the link function. Third, like the Albert and Chib (1993) method, it requires only a single layer of latent variables."

### サンプリング手法の効率性（Section 4.1）
> "For the basic $\mathrm{PG}(1, c)$ case, the sampler is very efficient: it requires only exponential and inverse-Gaussian draws, and the probability of accepting a proposed draw is uniformly bounded below at 0.99919. The method is also fully automatic, with no tuning needed to get optimal performance."

### 一様エルゴード性（Section 1）
> "Furthermore, in a recent paper based on an early technical report of our method, Choi and Hobert (2013) have proven that the Pólya-Gamma Gibbs sampler for Bayesian logistic regression is uniformly ergodic. This result has important practical consequences; most notably, it guarantees the existence of a central limit theorem for Monte Carlo averages of posterior draws. We are aware of no similar result for any other MCMC-based approach to the Bayesian logit model."

### 混合モデルへの拡張の容易さ（Section 3.3）
> "But the real advantage of data augmentation, and the Pólya-Gamma technique in particular, is that it becomes easy to construct and fit more complicated models. For instance, the Pólya-Gamma method trivially accommodates mixed models, factor models, and models with a spatial or dynamic structure."

### 計算効率の実証（Section 3.3）
> "Together these changes require just a few lines of code and a few extra seconds of runtime compared to the non-hierarchical logit model. A posterior draw of 2,000 samples for this data set takes 26.1 seconds for a binomial logistic regression, versus 27.3 seconds for a binomial logistic mixed model."
