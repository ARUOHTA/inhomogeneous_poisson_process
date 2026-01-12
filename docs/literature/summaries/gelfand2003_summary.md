# Gelfand et al. (2003) - Spatially Varying Coefficient Processesによる空間モデリング

## 基本情報
- **タイトル**: Spatial Modeling With Spatially Varying Coefficient Processes
- **著者**: Alan E. Gelfand, Hyon-Jung Kim, C. F. Sirmans, Sudipto Banerjee
- **ジャーナル**: Journal of the American Statistical Association
- **年**: 2003
- **DOI**: 10.2307/30045248
- **関連論点**: 論点1（Spatial modeling, Gaussian processes）

## 主な貢献
1. **Spatially Varying Coefficients (SVC) model**: 回帰係数が空間的に変化することを明示的にモデル化し、Gaussian processを用いて係数の空間サーフェスを表現
2. **Bayesian hierarchical framework**: 完全なBayesian推論により、observedおよびunobserved locationsでの係数プロセスの事後分布を取得
3. **Multivariate spatial processes**: 複数の係数プロセス間の依存性をseparable covariance構造でモデル化
4. **Spatio-temporal extensions**: 静的モデルから時空間動的モデルへの4つの拡張を提案
5. **Slice Gibbs sampler**: チューニング不要で効率的なMCMCアルゴリズム

## 手法の概要

### 従来の空間モデル（背景）

**Standard spatial process model**（Cressie 1993）:
$$Y(s) = \mu(s) + W(s) + \epsilon(s)$$

ここで：
- $\mu(s) = x(s)^T \beta$: Mean structure（fixed coefficients）
- $W(s) \sim \mathcal{GP}(0, \sigma^2, \rho(\cdot; \phi))$: Spatial random effects
- $\epsilon(s) \sim N(0, \tau^2)$: White noise（i.i.d.）

**問題点**: 係数$\beta$が region全体で constant（場所によって関係が変わる場合に不適）。

### Single Covariate SVC Model

**Model 1（Intercept process）**:

$W(s)$を random intercept process $\beta_0(s)$と解釈：
$$Y(s) = \beta_0 + \beta_1 x(s) + \beta_0(s) + \epsilon(s)$$

ここで：
- $\beta_0$: Overall intercept
- $\beta_0(s) \sim \mathcal{GP}(0, \sigma_0^2, \rho_0(\cdot; \phi_0))$: Random spatial adjustment to intercept
- $\tilde{\beta}_0(s) = \beta_0 + \beta_0(s)$: Random intercept process

**Marginalized likelihood**（$\beta_0$を積分消去）:
$$L(\beta_0, \beta_1, \tau^2, \sigma_0^2, \phi_0; y) = |\sigma_0^2 H_0(\phi_0) + \tau^2 I|^{-1/2} \exp\left\{-\frac{1}{2}(y - \beta_0 \mathbf{1} - \beta_1 x)^T (\sigma_0^2 H_0(\phi_0) + \tau^2 I)^{-1} (y - \beta_0 \mathbf{1} - \beta_1 x)\right\}$$

ここで$H_0(\phi_0)$は correlation matrix: $(H_0(\phi_0))_{ij} = \rho_0(s_i - s_j; \phi_0)$。

**Model 2（Slope process）**:
$$Y(s) = \beta_0 + \beta_1 x(s) + \beta_1(s) x(s) + \epsilon(s)$$

- $\beta_1(s) \sim \mathcal{GP}(0, \sigma_1^2, \rho_1(\cdot; \phi_1))$: Random spatial adjustment to slope
- $\tilde{\beta}_1(s) = \beta_1 + \beta_1(s)$: Random slope process

**Nonstationary property**:
$$\text{var}(Y(s) | \beta_0, \beta_1, \tau^2, \sigma_1^2, \phi_1) = x^2(s) \sigma_1^2 + \tau^2$$
$$\text{cov}(Y(s), Y(s') | \cdots) = \sigma_1^2 x(s) x(s') \rho_1(s - s'; \phi_1)$$

**重要な注意**: $x(s) > 0$が必要（centeringは不適切）。$x(s) \approx 0$では$Y(s)$と$Y(s')$の相関が消失。

**Model 3（Both intercept and slope processes）**:
$$Y(s) = \beta_0 + \beta_1 x(s) + \beta_0(s) + \beta_1(s) x(s) + \epsilon(s)$$

Bivariate process specification for $(\beta_0(s), \beta_1(s))$が必要（Section 3で詳述）。

Independence assumptionの下での marginal likelihood:
$$L(\beta_0, \beta_1, \tau^2, \sigma_0^2, \sigma_1^2, \phi_0, \phi_1; y) = |\sigma_0^2 H_0(\phi_0) + \sigma_1^2 D_x H_1(\phi_1) D_x + \tau^2 I|^{-1/2} \times \exp\{\cdots\}$$

ここで$D_x$は diagonal with $(D_x)_{ii} = x(s_i)$。

### Multivariate SVC Model（Multiple Covariates）

**General model**（$p$個のcovariates）:
$$Y(s) = X^T(s) \tilde{\beta}(s) + \epsilon(s)$$

ここで：
- $X(s) = (1, x_1(s), \ldots, x_{p-1}(s))^T$: $p \times 1$ covariate vector（intercept含む）
- $\tilde{\beta}(s) = \mu_\beta + \beta(s)$: Spatially varying coefficient vector
- $\beta(s) \sim \mathcal{GP}(\mathbf{0}, C(\cdot, \cdot))$: $p$-variate Gaussian process

**Separable covariance structure**（Mardia & Goodall 1993）:
$$C(s, s')_{lm} = \rho(s - s'; \phi) \tau_{lm}$$

ここで：
- $\rho(\cdot; \phi)$: Scalar correlation function（spatial dependence）
- $T = [\tau_{lm}]_{p \times p}$: Positive-definite symmetric matrix（cross-covariance at same location）

**Distribution of $\tilde{\beta}$**（$n$個の観測地点）:
$$\tilde{\beta} \sim N(\mathbf{1}_{n \times 1} \otimes \mu_\beta, H(\phi) \otimes T)$$

ここで：
- $\tilde{\beta}$: $np \times 1$ vector（concatenation of $\tilde{\beta}(s_1), \ldots, \tilde{\beta}(s_n)$）
- $H(\phi)$: $n \times n$ spatial correlation matrix
- $\otimes$: Kronecker product

**Marginalized likelihood**:
$$L(\mu_\beta, \tau^2, T, \phi; y) = |X(H(\phi) \otimes T)X^T + \tau^2 I|^{-1/2} \times \exp\left\{-\frac{1}{2}(y - X(\mathbf{1} \otimes \mu_\beta))^T (X(H(\phi) \otimes T)X^T + \tau^2 I)^{-1} (y - X(\mathbf{1} \otimes \mu_\beta))\right\}$$

ここで$X^T$は$n \times np$ block diagonal with blocks $X^T(s_i)$。

### Matérn Correlation Function

全てのモデルでMatérn classを使用：
$$\rho(h; \phi) \propto (\gamma \|h\|)^\nu K_\nu(\gamma \|h\|)$$

ここで：
- $K_\nu$: Modified Bessel function
- $\phi = (\gamma, \nu)$
- $\gamma$: Decay parameter（大きいほど相関が速く減衰）
- $\nu$: Smoothness parameter（$\nu > 1$でmean-squared differentiable）

**Range**: $\rho(r; \phi) = 0.05$となる距離$r$（practical range）。

### Inference

**Prior specification**（multivariate model）:
- $\mu_\beta \sim N(\mathbf{0}, 10^5 I)$（vague）
- $T \sim \text{IW}(p, \text{diag}(0.001))$（Inverted Wishart）
- $\tau^2 \sim \text{IG}(2, 1)$（Inverse Gamma）
- $\gamma, \nu \sim \text{Gamma}(2, 0.1)$

**Slice Gibbs sampler**（Agarwal & Gelfand 2001, Neal 2002）:
1. Auxiliary variable $U \sim \text{Uniform}(0, L(\mu_\beta, \tau^2, T, \phi; y))$を導入
2. Joint posterior: $f(\mu_\beta, \tau^2, T, \phi, U | y) \propto I(U < L) \pi(\mu_\beta, \tau^2, T, \phi)$
3. Gibbs sampler:
   - Update $U$: Uniform draw
   - Update parameters: Prior draw subject to $I(U < L)$

**利点**: Tuning不要、Metropolis alternativesより速い収束、autocorrelation問題の回避。

**Prediction**:

1. **Coefficient process at observed locations**:
$$f(\tilde{\beta} | \mu_\beta, \tau^2, T, \phi, y) = N(B b, B)$$
   - $B = (X^T X / \tau^2 + H^{-1}(\phi) \otimes T^{-1})^{-1}$
   - $b = X^T y / \tau^2 + (H^{-1}(\phi) \otimes T^{-1})(\mathbf{1} \otimes \mu_\beta)$

2. **Coefficient at new location $s_{\text{new}}$**:
$$f(\tilde{\beta}(s_{\text{new}}) | \tilde{\beta}, \mu_\beta, T, \phi) = N(\mu_\beta + (h_{\text{new}}^T(\phi) H^{-1}(\phi) \otimes I)(\tilde{\beta} - \mathbf{1} \otimes \mu_\beta), (1 - h_{\text{new}}^T(\phi) H^{-1}(\phi) h_{\text{new}}(\phi))T)$$
   - $h_{\text{new}}(\phi)$: $n \times 1$ vector with entries $\rho(s_i - s_{\text{new}}; \phi)$

3. **Response at new location**:
$$f(Y(s_{\text{new}}) | y) = \int f(Y(s_{\text{new}}) | \beta_0, \beta_1, \beta_0(s_{\text{new}}), \tau^2) f(\beta_0(s_{\text{new}}) | \beta_0, \sigma_0^2, \phi_0) f(\beta_0, \beta_1, \tau^2, \sigma_0^2, \phi_0 | y)$$

### Spatio-Temporal Extensions（Section 4）

**General framework**:
$$Y(s, t) = X^T(s, t) \tilde{\beta}(s, t) + \epsilon(s, t), \quad t = 1, 2, \ldots, M$$

**4つのモデル**:

1. **Model 1**: $\beta(s, t) = \beta(s)$（Constant over time、longitudinal growth curve analog）
2. **Model 2**: $\beta(s, t) = \beta(s) + \alpha(t)$（Additive space-time）
   - 2a: $\alpha_k(t)$ i.i.d. time dummies
   - 2b: $\alpha(t)$ follows random walk/AR process
3. **Model 3**: $\beta(s, t) = \beta^{(t)}(s)$（Nested within time）
   - Independent spatial processes at each $t$
   - 3a: Common $\mu_\beta$ across time
   - 3b: Time-varying $\mu_\beta^{(t)}$
4. **Model 4**: Separable covariance in space and time
   - $\Sigma_{[\beta(s, t), \beta(s', t')]} = \rho^{(1)}(s - s'; \phi) \times \rho^{(2)}(t - t'; \gamma) T$

**Computational note**: $n$ sites × $T$ time points → $nT \times nT$ matrices。

### Model Comparison（Posterior Predictive Loss）

**Gelfand & Ghosh (1998) criterion**:
$$D_k = \sum_{(s,t)} (Y_{\text{obs}}(s, t) - \mu(s, t))^2 + k \sum_{(s,t)} \sigma^2(s, t)$$

- First term（$G$）: Goodness-of-fit
- Second term（$P$）: Penalty for model complexity
- Typically $k = 1$
- Smallest $D$ → best model

## 実証結果

### Real Estate Data（Baton Rouge, LA, 1985-1992）

**データ**: Log selling price of single-family homes

**Covariates**:
1. Age of house
2. Square feet of living area
3. Square feet of other area（garage、carport等）
4. Number of bathrooms

**Static spatial models（1992年、$n=237$）**:

| Model | $G$ | $P$ | $D$ |
|-------|-----|-----|-----|
| **2D models**（intercept + 1 SVC）: ||||
| Living area | 69.87 | 46.24 | **116.11** |
| Age | 74.52 | 44.58 | 119.10 |
| Other area | 70.24 | 49.87 | 120.11 |
| Bathrooms | 78.02 | 52.93 | 130.95 |
| **3D models**（intercept + 2 SVCs）: ||||
| Age, living area | 61.38 | 47.83 | **109.21** |
| Age, other area | 63.80 | 48.45 | 112.25 |
| Living area, other area | 72.35 | 50.75 | 123.10 |
| **5D models**（intercept + 4 SVCs）: ||||
| Dependent process | 42.21 | 36.01 | **78.22** |
| Independent process | 94.36 | 59.34 | 153.70 |

**重要な発見**:
- Dependent process model（separable covariance）が圧倒的に優れる
- Independence assumptionは inadequate
- 全ての係数を空間的に変化させる5D modelが最良

**5D Dependent Model posterior summary**（Table 4）:

- $\beta_0$（intercept）: 9.917（median）
- $\beta_1$（age）: -0.005（負の効果、期待通り）
- $\beta_2$（living area）: 0.341（正の効果）
- $\beta_3$（other area）: 0.313
- $\beta_4$（bathrooms）: 0.292
- Range: 4.17 km（median、parish全体は約22 km × 33 km）
- Smoothness $\nu$: 1.47（mean-squared differentiable realizations）

**Spatial variance contributions**:
- $T_{11}$（intercept process）: 0.322（largest contribution）
- $\bar{x}_4^2 T_{55}$（bathrooms process）: 0.183（second）
- $\tau^2$（pure error）: 0.049（smallest）

**Correlations between processes**（$T_{lm} / \sqrt{T_{ll} T_{mm}}$）:
- Intercept-Age: -0.203（negative、expected）
- Intercept-Living area: -0.186
- Intercept-Bathrooms: -0.583（強い負の相関）
- Other area-Bathrooms: -0.839（非常に強い負の相関）

**Cross-validation（20 holdout sites、Table 5）**:
- 5D model: 95% predictive intervals（19/20 coverage）、最も狭い区間
- Independence model: 18/20 coverage、非常に広い区間

### Spatio-Temporal Models（1989-1992、$n=120$ distinct locations）

**Model comparison**（Table 6）:

| Model | $G$ | $P$ | $D$ |
|-------|-----|-----|-----|
| **Dependent process**: ||||
| Model 1（constant $\beta(s)$） | 54.54 | 29.11 | 83.65 |
| Model 2a（time dummies） | 47.92 | 26.95 | 74.87 |
| Model 2b（temporal process） | 43.38 | 29.10 | 72.48 |
| Model 3a（nested、common $\mu_\beta$） | 43.74 | 20.63 | 64.37 |
| Model 3b（nested、$\mu_\beta^{(t)}$） | 42.35 | 21.04 | **63.39** |
| Model 4（separable） | 37.84 | 26.47 | 64.31 |

**結果**: Model 3b（space nested within time）が最良、Model 4も competitive。

**Temporal evolution（Model 3b、Table 7）**:
- $\mu_\beta^{(t)}$: 時間的に安定（わずかな変動）
- Spatial range: 年によって変化（1989: 2.95 km、1991: 4.55 km）

### Generalized Linear Models（Section 7）

**GLM extension**:
$$f(y(s_i) | \theta(s_i)) = h(y(s_i)) \exp(\theta(s_i) y(s_i) - b(\theta(s_i)))$$

Canonical link: $\theta(s_i) = X^T(s_i) \tilde{\beta}(s_i)$

**Challenges**:
- Componentwise updating（Gibbs sampler with adaptive rejection sampling）
- Hierarchically centered parameterization（still slow mixing）
- Langevin diffusion for blocked Metropolis（computationally intensive）

**結論**: GLMケースでの効率的なfitting strategyは今後の課題。

## 限界・残された課題

### 著者が指摘する限界
1. **GLM settingでの computational efficiency**: Section 7で指摘されているように、GLMへの拡張では componentwise updating が遅く、blocked Metropolis も tractable でない
2. **Parametric correlation function choice**: Matérn classを使用しているが、他の選択肢（exponential、Gaussian等）との比較は未検討
3. **Nonstationarity in spatial correlation**: $\rho(\cdot; \phi)$はstationaryで、spatial rangeが場所によって変わることは許容されない
4. **Measurement error**: Covariates $X(s)$の測定誤差は考慮されていない

### 本研究の視点からの限界
1. **Large-scale computational challenges**: $n$が数千を超えるとGP計算（$O(n^3)$）が prohibitive。NNGPのような近似手法の議論なし
2. **Cross-covariance structure choice**: Separable covariance（$\rho(s - s'; \phi) \tau_{lm}$）は computationally convenientだが制約的。Non-separable structureへの拡張は未検討
3. **Spatial misalignment**: Covariates $X(s)$が異なる空間解像度で観測される場合の change-of-support問題は未対応
4. **Point pattern dataへの拡張なし**: 本論文はresponse $Y(s)$が point-referenced dataを想定。Point process intensityへの適用は未検討

## 本研究(MMCP)との関係

### 直接的な関連
- **Spatially varying coefficients**: 本研究で距離から遺跡密度への影響（coefficient）が場所によって異なる可能性を考慮
- **Gaussian process modeling**: Coefficient processをGPでモデル化する考え方は、本研究のintensity functionのGP modelingと共通
- **Bayesian hierarchical framework**: 完全なBayesian推論により不確実性を定量化

### 本研究での具体的適用シナリオ

**Scenario 1: Distance effectの空間変化**

本研究でobsidian産地からの距離$d_k(s)$が遺跡密度$\lambda(s)$に与える影響が場所によって異なる場合：

**Standard IPP model**（constant coefficients）:
$$\log \lambda(s) = \beta_0 + \sum_{k=1}^K \beta_k d_k(s)$$

**SVC-IPP model**（spatially varying coefficients）:
$$\log \lambda(s) = \beta_0 + \sum_{k=1}^K (\beta_k + \beta_k(s)) d_k(s)$$

ここで$(\beta_1(s), \ldots, \beta_K(s))^T \sim \mathcal{GP}(\mathbf{0}, H(\phi) \otimes T)$。

**期待される効果**:
- 産地近傍: $\beta_k(s)$が大きく負（強い距離減衰効果）
- 遠方: $\beta_k(s)$が0に近い（距離効果が弱い）
- Separable covariance $T$でdistance coefficients間の依存性をモデル化

**Scenario 2: Marked point process with SVC**

Point pattern: 遺跡の空間分布
Mark: Obsidian source（discrete）+ Composition（continuous）

**Intensity for source $k$**:
$$\lambda_k(s) = \lambda^* \exp\{X^T(s) \tilde{\beta}_k(s)\}$$

ここで：
- $X(s)$: Spatial covariates（elevation、distance to water等）
- $\tilde{\beta}_k(s) \sim \mathcal{GP}(\mu_{\beta_k}, C_k(\cdot, \cdot))$: Source-specific SVC process

Multivariate GP（$K$個のsources）で$\tilde{\beta}_1(s), \ldots, \tilde{\beta}_K(s)$間の依存性をモデル化。

### Gelfand et al. (2003) vs NNGPの統合

**Challenge**: SVC modelは$np \times np$行列（$n$ locations × $p$ coefficients）の逆行列計算が必要→ large $n$でprohibitive。

**Solution direction**（本研究での適用）:
- Replace full GP with NNGP（Datta et al. 2016）
- Sparse precision matrix: $O(nm^3)$ where $m \ll n$
- Separable covariance構造を維持したままNNGPを適用

**Implementation**:
$$\tilde{\beta}(s) | \tilde{\beta}(N(s)) \sim N(\mu_\beta + C(s, N(s)) C(N(s), N(s))^{-1} (\tilde{\beta}(N(s)) - \mu_\beta), C(s, s) - C(s, N(s)) C(N(s), N(s))^{-1} C(N(s), s))$$

ここで$N(s)$は$s$の$m$個のnearest neighbors。

### 本研究での限界として明記すべき点
第3章で引用する際、以下を明確に述べるべき：
- SVC modelは computationally intensive（本研究で数千遺跡を扱う場合、NNGPへの拡張が必須）
- Separable covarianceは tractableだが、異なる係数プロセス間の空間相関構造が場所によって変わることは許容されない
- Coefficients $\beta_k(s)$の解釈には注意が必要：$x(s) > 0$の条件、centering不適切等
- GLM case（例：Poisson regression for counts）での実装は computationally challenging

## 引用すべき箇所

### SVC modelの動機
> "In nearly all of this work, the regression coefficients are assumed to be constant across the region. In certain applications, this would not be appropriate. The coefficients may be expected to vary at the local or subregion level." (Introduction)

Constant coefficientsの限界とSVCの必要性を説明する引用。

### GP-based approach vs parametric approach
> "One possible approach would be to model the spatial surface for the coefficient parametrically...However, this requires selection of a spline function and determination of the number of and locations of the knots in the space...The approach that we adopt here is arguably more natural and at least as flexible. We model the spatially varying coefficient surface as a realization from a spatial process." (Introduction)

GP-basedアプローチの優位性を説明する引用。

### Bayesian frameworkの利点
> "We adopt a Bayesian approach for our modeling framework. This is attractive in the proposed setting, because we are specifically interested in inference for the random spatial effects. In particular, we obtain an entire posterior for the spatial coefficient process at both observed and unobserved locations...Interpolation for a process that is neither observed nor arising as a residual seems inaccessible in any other framework." (Introduction)

Bayesian推論の必要性（係数プロセスの予測）を説明する引用。

### Nonstationarity property
> "Note that (8) provides a heterogeneous, nonstationary process for the data regardless of the choice of covariance function for the $\beta_1(s)$ process. Here, $\text{var}(Y(s) | \cdots) = x^2(s) \sigma_1^2 + \tau^2$ and $\text{cov}(Y(s), Y(s') | \cdots) = \sigma_1^2 x(s) x(s') \rho_1(s - s'; \phi_1)$." (Section 2)

SVC modelが自動的にnonstationaryになる性質を説明する引用。

### Centering issueの警告
> "As a result, we observe that in practice, (8) is sensible only if we have $x(s) > 0$. In fact, centering and scaling, usually advocated for better-behaved model fitting, is inappropriate here." (Section 2)

Covariatesのcentering/scalingに関する重要な注意（実装時に critical）。

### Dependence between processes
> "In practice, to assume that the component processes of $\tilde{\beta}(s)$ are independent is likely inappropriate. That is, in the simpler case of simple linear regression, negative association between slope and intercept is usually seen." (Section 3)

Multivariate GP（依存性モデル化）の必要性を説明する引用。

### Separable covarianceの定義
> "We require a valid $p$-variate choice. In the sequel we work with a computationally convenient separable choice following Mardia and Goodall (1993)...Specifically, let...$$C(s, s')_{lm} = \rho(s - s'; \phi) \tau_{lm}$$" (Section 3)

Separable covariance structureの定義（computational tractabilityとのトレードオフ）。

### Model comparison result
> "From Table 3 the five-dimensional model using (13) is far superior, and the independence model is clearly worst, supporting our earlier intuition." (Section 6)

Dependent multivariate processの優位性を示す実証結果。

### GLM settingの課題
> "Clearly, more work is needed to provide efficient model-fitting strategies for these models." (Section 7)

GLM extensionにおける computational challengesを認める引用（将来の研究課題）。
