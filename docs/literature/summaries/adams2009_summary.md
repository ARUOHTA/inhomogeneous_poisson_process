# Adams et al. (2009) - Tractable Inference in Gaussian Cox Processes

## 基本情報
- **タイトル**: Tractable Nonparametric Bayesian Inference in Poisson Processes with Gaussian Process Intensities
- **著者**: Ryan Prescott Adams, Iain Murray, David J.C. MacKay
- **会議**: International Conference on Machine Learning (ICML)
- **年**: 2009
- **関連論点**: 論点1（Gaussian Cox processes, Computational methods）、論点3（Point process models）

## 主な貢献
1. **Sigmoidal Gaussian Cox Process (SGCP)**: Gaussian Cox processにおいて、approximationやfinite-dimensional proxy distributionを導入せずに、完全ノンパラメトリックベイズ推論を可能にする初めての手法
2. **Tractable data generation via thinning**: GP-based intensity functionから exact Poisson dataを生成する方法（Algorithm 1）
3. **Data augmentation for tractable inference**: Thinned events（拒否されたイベント）をlatent variablesとして導入することで、intractable integralを回避
4. **Hamiltonian Monte Carlo integration**: Infinite-dimensional function $g(s)$の効率的な事後サンプリング

## 手法の概要

### Inhomogeneous Poisson Process (IPP)

Domain $\mathcal{S} = \mathbb{R}^D$上のPoisson processで、intensity function $\lambda(s): \mathcal{S} \to \mathbb{R}^+$でパラメータ化：
- Subregion $\mathcal{T} \subset \mathcal{S}$内のイベント数：$N(\mathcal{T}) \sim \text{Poisson}(\lambda_{\mathcal{T}})$
- $\lambda_{\mathcal{T}} = \int_{\mathcal{T}} \lambda(s) ds$

### Gaussian Cox Process (GCP)

**Log Gaussian Cox Process（従来）**:
$$\log \lambda(s) \sim \mathcal{GP}(m(\cdot), C(\cdot, \cdot))$$

**Likelihood（intractable）**:
$$p(\{s_k\}_{k=1}^K | \lambda(s)) = \exp\left\{-\int_{\mathcal{T}} \lambda(s) ds\right\} \prod_{k=1}^K \lambda(s_k)$$

積分$\int_{\mathcal{T}} \lambda(s) ds$が計算不可能（infinite-dimensional random function）。

### Sigmoidal Gaussian Cox Process (SGCP)

**Key transformation**:
$$\lambda(s) = \lambda^* \sigma(g(s))$$

ここで：
- $g(s) \sim \mathcal{GP}(m(\cdot), C(\cdot, \cdot))$: Gaussian process
- $\sigma(z) = (1 + e^{-z})^{-1}$: Logistic function（sigmoid）
- $\lambda^*$: Upper bound on $\lambda(s)$

**利点**: $0 \leq \lambda(s) \leq \lambda^*$が保証され、thinning algorithmが適用可能。

### Exact Data Generation via Thinning (Algorithm 1)

**手順**:

1. **Homogeneous Poisson process sampling**:
   - $J \sim \text{Poisson}(\lambda^* \mu(\mathcal{T}))$（イベント数）
   - $\{\hat{s}_j\}_{j=1}^J \sim \text{Uniform}(\mathcal{T})$（イベント位置）

2. **GP sampling at proposed locations**:
   - $\{g(\hat{s}_j)\}_{j=1}^J \sim \mathcal{GP}(m(\cdot), C(\cdot, \cdot), \theta, \{\hat{s}_j\}_{j=1}^J)$

3. **Thinning（rejection sampling）**:
   - 各$j$について：$r_j \sim \text{Uniform}(0,1)$
   - Accept $\hat{s}_j$ if $r_j < \sigma(g(\hat{s}_j))$
   - Accepted events: $\mathcal{E} = \{s_k\}_{k=1}^K$

**重要な性質**: $\{s_k\}_{k=1}^K$はintensity $\lambda(s) = \lambda^* \sigma(g(s))$を持つinhomogeneous Poisson processからのexact sample。

### Inference via Data Augmentation

**Latent variables**:
1. $M$: Number of thinned（rejected）events
2. $\{\tilde{s}_m\}_{m=1}^M$: Locations of thinned events
3. $\boldsymbol{g}_M = \{g(\tilde{s}_m)\}_{m=1}^M$: Function values at thinned events
4. $\boldsymbol{g}_K = \{g(s_k)\}_{k=1}^K$: Function values at observed events

**Joint distribution（tractable）**:
$$p(\{s_k\}_{k=1}^K, M, \{\tilde{s}_m\}_{m=1}^M, \boldsymbol{g}_{M+K} | \lambda^*, \mathcal{T}, \theta) = (\lambda^*)^{K+M} \exp\{-\lambda^* \mu(\mathcal{T})\} \prod_{k=1}^K \sigma(g(s_k)) \prod_{m=1}^M \sigma(-g(\tilde{s}_m)) \times \mathcal{GP}(\boldsymbol{g}_{M+K} | \{s_k\}_{k=1}^K, \{\tilde{s}_m\}_{m=1}^M, \theta)$$

**重要な点**: $\sigma(-z) = 1 - \sigma(z)$なので、thinned eventsはaccept probabilityの補集合で重み付け。

### MCMC Sampling Scheme

#### 1. Sampling the number of thinned events $M$

**Insertion move**: $M \to M+1$
- Propose new location: $\tilde{s}' \sim \text{Uniform}(\mathcal{T})$
- Sample GP value: $g(\tilde{s}') | \boldsymbol{g}_{M+K} \sim \mathcal{GP}$
- Acceptance ratio:
$$a_{\text{ins}} = \frac{(1-b(K, M+1)) \mu(\mathcal{T}) \lambda^*}{(M+1) b(K, M) (1 + \exp\{g(\tilde{s}')\})}$$

**Deletion move**: $M \to M-1$
- Randomly select thinned event $m$ to remove
- Acceptance ratio:
$$a_{\text{del}} = \frac{M b(K, M-1) (1 + \exp\{g(\tilde{s}_m)\})}{(1-b(K, M)) \mu(\mathcal{T}) \lambda^*}$$

ここで$b(K, M)$はinsertion probabilityで、typically $b(K, M) = 0.5$。

#### 2. Sampling the locations of thinned events $\{\tilde{s}_m\}_{m=1}^M$

各thinned event $m$について：
- Propose new location: $\tilde{s}_m' \sim q(\cdot \leftarrow \tilde{s}_m)$（perturbative proposal）
- Sample GP value: $g(\tilde{s}_m') | \boldsymbol{g}_{M+K} \sim \mathcal{GP}$
- Acceptance ratio:
$$a_{\text{loc}} = \frac{q(\tilde{s}_m \leftarrow \tilde{s}_m') (1 + \exp\{g(\tilde{s}_m)\})}{q(\tilde{s}_m' \leftarrow \tilde{s}_m) (1 + \exp\{g(\tilde{s}_m')\})}$$

#### 3. Sampling the function values $\boldsymbol{g}_{M+K}$

**Hamiltonian Monte Carlo（HMC）**を使用：

Log conditional posterior:
$$\ln p(\boldsymbol{g}_{M+K} | M, \{s_k\}_{k=1}^K, \{s_m\}_{m=1}^M, \theta) = -\frac{1}{2} \boldsymbol{g}_{M+K}^T \boldsymbol{\Sigma}^{-1} \boldsymbol{g}_{M+K} - \sum_{k=1}^K \ln(1 + \exp\{-g(s_k)\}) - \sum_{m=1}^M \ln(1 + \exp\{g(\tilde{s}_m)\}) + \text{const.}$$

Gradientが解析的に計算可能なため、HMCで効率的にサンプリング。

**Whitening transformation**: Cholesky decomposition $\boldsymbol{\Sigma} = \boldsymbol{L}\boldsymbol{L}^T$を使い、$\boldsymbol{L}^{-1}\boldsymbol{g}_{M+K}$の空間でHMCを実行（better-conditioned）。

#### 4. Hyperparameter inference

**GP hyperparameters $\theta$**:
- HMC for continuous hyperparameters（lengthscale、variance等）
- Conditional on data、thinned events、function values

**Upper bound $\lambda^*$**:
- Conjugate Gamma prior: $\lambda^* \sim \text{Gamma}(\alpha, \beta)$
- Conditional posterior（given $M, K$）:
$$\lambda^* | M, K \sim \text{Gamma}(\alpha + K + M, \beta + \mu(\mathcal{T}))$$

$\{\hat{s}_j\}_{j=1}^{K+M} = \{s_k\}_{k=1}^K \cup \{\tilde{s}_m\}_{m=1}^M$はhomogeneous Poisson process（rate $\lambda^*$）からのsampleなので、conjugacy成立。

### Predictive Distribution

**新しいイベント列の生成**:
1. Current MCMC state: $\boldsymbol{g}_{M+K}, \{\tilde{s}_m\}_{m=1}^M, \{s_k\}_{k=1}^K, \theta, \lambda^*$
2. Algorithm 1を実行するが、Line 4でcurrent functionにcondition:
   - $\{g(\hat{s}_j)\}_{j=1}^J \sim \mathcal{GP}(\cdot | \{s_k, g(s_k)\}_{k=1}^K, \{\tilde{s}_m, g(\tilde{s}_m)\}_{m=1}^M, \theta)$

これにより、posterior over intensity functionsを積分した predictive distributionからのsampleが得られる。

## 実証結果

### Synthetic Data（3つのシナリオ）

**設定**:
1. $\lambda_1(s) = 2\exp\{-s/15\} + \exp\{-((s-25)/10)^2\}$ on $[0,50]$（53 events）
2. $\lambda_2(s) = 5\sin(s^2) + 6$ on $[0,5]$（29 events）
3. $\lambda_3(s)$: Piecewise linear on $[0,100]$（235 events）

**比較手法**:
- **SGCP**: 本論文の手法
- **KS**: Kernel smoothing（Diggle 1985）
- **LGCP10/25/100**: Log Gaussian Cox Process with discretization（Møller et al. 1998）、10/25/100 bins

**評価指標**:
- $\ell_2$ distance: $\|\hat{\lambda}(s) - \lambda_{\text{true}}(s)\|_2$
- Mean log predictive probability（lp）: 10個のhold-out time seriesでの平均対数尤度

**結果（Table 1）**:

| Method | $\lambda_1$ $\ell_2$ | $\lambda_1$ lp | $\lambda_2$ $\ell_2$ | $\lambda_2$ lp | $\lambda_3$ $\ell_2$ | $\lambda_3$ lp |
|--------|---------------------|----------------|---------------------|----------------|---------------------|----------------|
| **SGCP** | **4.20** | **-45.11** | **38.38** | 24.45 | 11.41 | **-43.39** |
| KS | 6.65 | -46.41 | 73.71 | **28.19** | 30.56 | -46.47 |
| LGCP100 | 5.44 | -45.24 | 43.51 | 25.29 | **10.79** | -47.16 |

- **SGCP**: $\lambda_1$と$\lambda_3$で最良、$\lambda_2$でも competitive
- **Kernel smoothing**: $\lambda_2$でのみ良好（sinusoidは滑らかで KSに有利）
- **LGCP**: Discretizationに依存し、bins数の選択が critical

### Coal Mining Disaster Data（Real data）

**データ**:
- 1875年3月15日～1962年3月22日のイギリスの炭鉱爆発事故（死者10人以上）
- 191 events

**結果（Figure 5）**:
- Inferred intensity $\lambda(s)$: 19世紀末に高く、20世紀半ばに減少
- Upper bound $\lambda^*$: Posterior mean $\approx 3.5$ events/year
- Number of thinned events $M$: Posterior mean $\approx 100$

### Redwoods Data（Spatial data）

**データ**:
- Ripley (1977)のRedwood treesの位置データ
- 195 points on unit square

**結果（Figure 4）**:
- Mean intensity estimate: Spatial clusteringを捉える
- Thinned events histogram: Intensity が低い領域に集中（"peg down"効果）

**Thinned eventsの役割**:
- Low-intensity regionsでGPをanchor
- Intensity functionが0に近い場所を"pegging down"

## 限界・残された課題

### 著者が指摘する限界
1. **Computational complexity**: $O((K+M)^3)$ per MCMC step（GaussianProcessの逆行列計算）
   - Several thousand eventsを超えると計算不可能
   - 大規模データへの拡張が必要（sparse GP、NNGP等）

2. **Choice of $\lambda^*$**: Upper bound $\lambda^*$の選択がpriorとして重要
   - Too large: 多数のthinned events（$M$大）、計算コスト増
   - Too small: モデルが真のintensityを表現できない
   - Hierarchical prior（Gamma）でinferenceは可能だが、initial choiceは依然として重要

3. **Choice of sigmoid $\sigma(\cdot)$**: Logistic function以外の選択肢
   - Normal CDF: Marginal uniformity of $\lambda(s)$（if zero-mean stationary GP）
   - 他のsigmoids: 異なるprior beliefsを表現

4. **Mixing of Markov chain**: Thinned eventsの数$M$と位置$\{\tilde{s}_m\}$のmixing
   - HMCで function values $\boldsymbol{g}_{M+K}$のmixingは改善
   - しかし$M$と$\{\tilde{s}_m\}$のMetropolis-Hastingsは依然として遅い可能性

### 本研究の視点からの限界
1. **NNGPとの統合なし**: 大規模データへのscalability issues（$O(n^3)$）を解決する近似手法（NNGP）の議論なし
2. **Spatial covariatesの欠如**: Intensity $\lambda(s)$がcovariatesに依存するモデルへの拡張（例：$\log \lambda(s) = \beta^T x(s) + g(s)$）が未検討
3. **Marksの考慮なし**: Marked point processへの拡張なし
4. **Temporal dynamicsのみ**: Spatial-temporal modelsへの拡張は示唆されているが未実装
5. **Preferential samplingの無視**: Sampling locationとintensityの相関（Gelfand & Shirota 2019）は考慮されていない

## 本研究(MMCP)との関係

### 直接的な関連
- **Gaussian Cox processの実装**: 本研究でLog GCPを使用する場合、Adams et al. (2009)のthinning-based data augmentationが直接適用可能
- **Tractable inference**: Intractable integralを回避する data augmentation strategy
- **Latent variablesの役割**: Thinned eventsの考え方は、他のlatent variable models（例：augmented likelihood）にも示唆的

### 本研究での具体的適用シナリオ

**Scenario 1: Unmarked point process（遺跡分布）**
- $\{s_i\}_{i=1}^n$: 遺跡の空間位置
- Intensity: $\lambda(s) = \lambda^* \sigma(g(s))$
- GP prior: $g(s) \sim \mathcal{GP}(m(s), C(s, s'))$
- Covariate-dependent mean: $m(s) = \beta^T x(s)$（産地からの距離等）

SGCPを直接適用し、遺跡密度の空間変化を推定。

**Scenario 2: Marked point process with compositional marks**
- Point process: SGCPで遺跡の空間分布
- Mark distribution: $Y_i | s_i \sim f(y | s_i, \theta)$（compositional、Wehrhahn et al. 2022）

Hierarchical model:
1. Spatial process: SGCP with thinning（Adams et al. 2009）
2. Marks given location: DMBPP（Wehrhahn et al. 2022）

**Computational challenge**: $O((K+M)^3)$がボトルネック。本研究で数千遺跡を扱う場合、NNGPへの拡張が必須。

### Adams et al. (2009) vs Møller et al. (1998)

| 側面 | Adams et al. (2009) | Møller et al. (1998) |
|------|---------------------|----------------------|
| Approach | Thinning + Data augmentation | Discretization |
| Dimensionality | Infinite-dimensional（true GP） | Finite-dimensional proxy |
| Approximation | None（exact inference） | Bin-based approximation |
| Computational cost | $O((K+M)^3)$ | $O(B^3)$（$B$ = bins） |
| Flexibility | GP prior fully respected | GP prior approximated |
| Scalability | Limited（several thousand events） | Better for large data |

**本研究での選択**:
- 小規模データ（< 1000 events）: Adams et al. (2009)が理論的に優位
- 大規模データ（> 10,000 events）: Discretization（Møller et al. 1998）またはNNGP-based approximation

### 本研究での拡張の方向性

**Extension 1: Covariate-dependent intensity**
$$\lambda(s) = \lambda^* \sigma(\beta^T x(s) + g(s))$$

ここで$\beta$はcoefficients、$g(s)$はGP residual。

**Extension 2: NNGP for scalability**
- Replace full GP $\mathcal{GP}(\boldsymbol{g}_{M+K})$ with NNGP approximation
- $O((K+M)m^3)$ where $m \ll K+M$（nearest neighbors）
- Thinning algorithmはそのまま適用可能

**Extension 3: Marked SGCP**
- Intensity: $\lambda_k(s) = \lambda^* \sigma(g_k(s))$ for mark $k$
- Joint GP: $(g_1(s), \ldots, g_K(s))^T \sim \mathcal{GP}(\boldsymbol{m}, \boldsymbol{C})$
- Cross-covariance $\boldsymbol{C}$でmarks間の依存関係をモデル化

## 引用すべき箇所

### Gaussian Cox processの定義と問題
> "The inhomogeneous Poisson process is a point process that has varying intensity across its domain (usually time or space). For nonparametric Bayesian modeling, the Gaussian process is a useful way to place a prior distribution on this intensity. The combination of a Poisson process and GP is known as a Gaussian Cox process, or doubly-stochastic Poisson process. Likelihood-based inference in these models requires an intractable integral over an infinite-dimensional random function." (Abstract)

Gaussian Cox processの定義とintractabilityの問題を簡潔に述べた引用。

### 本論文の貢献
> "In this paper we present the first approach to Gaussian Cox processes in which it is possible to perform inference without introducing approximations or finite-dimensional proxy distributions." (Abstract)

本論文の主要貢献（approximation-free inference）を明確に述べた引用。

### Thinning algorithmの原理
> "We use the transformation of Equation 1 because it allows us to simulate exact Poisson data from a random intensity function drawn from the prior provided by the Gaussian process. By exact we mean that the data are not biased by, for example, the starting state of a finite Markov chain. We generate these exact data via thinning, which is a point-process variant of rejection sampling..." (Section 2.3)

Thinning algorithmの動機と"exact"の意味を説明する引用。

### Data augmentationの利点
> "By considering the procedure as a latent variable model, we inherit these convenient properties for inference. The generative procedure did not require integrating an infinite-dimensional random function, nor did it require knowledge of $g(s)$ or $\lambda(s)$ at more than a finite number of locations." (Section 3.1)

Data augmentationによる tractabilityの獲得を説明する引用。

### Computational limitations
> "Gaussian processes have significant computational demands: they have $O(n^3)$ time complexity for $n$ input points, and $O(n^2)$ space complexity. When performing inference in the SGCP model, this means that each MCMC step costs $O((K+M)^3)$ as the thinned events must be included in the GP. Thus the approach we present is infeasible for data sets that have more than several thousand events." (Section 5.1)

Computational complexityの限界を明確に述べた引用（本研究での適用可能性を議論する際に重要）。

### Thinned eventsの役割（spatial context）
> "As expected, it is approximately a 'negative' of the mean intensity; the thinned events are moving to places where it is necessary to 'peg down' the intensity function." (Section 4.3、Redwoods data)

Thinned eventsがlow-intensity regionsでGPを"peg down"する役割を説明する引用（intuitionとして有用）。

### 従来手法との比較
> "The Sigmoidal Gaussian Cox Process is superior to the frequentist kernel density approach (Diggle, 1985) in several ways. First, we obtain samples from the posterior rather than a point estimate of the unknown intensity function. Second, we are able to perform bandwidth selection in a principled way by sampling from the hyperparameters of the Gaussian process." (Section 5.2)

Bayesian approachの利点を説明する引用。

### Model variationsの可能性
> "There are several ways in which the Sigmoidal Gaussian Cox Process we have presented could be modified for different modeling situations. For example, to arrive at bounded random intensities we used a constant dominating function $\lambda^*$, but other tractable parametric forms would be suitable and potentially more efficient." (Section 5.3)

モデルの拡張可能性を述べた引用（本研究での応用を検討する際の示唆）。
