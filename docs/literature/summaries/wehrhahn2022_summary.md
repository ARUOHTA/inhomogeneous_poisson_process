# Wehrhahn et al. (2022) - Compositional DataのためのDependent Bayesian Nonparametric Modeling

## 基本情報
- **タイトル**: Dependent Bayesian nonparametric modeling of compositional data using random Bernstein polynomials
- **著者**: Claudia Wehrhahn, Andrés F. Barrientos, Alejandro Jara
- **ジャーナル**: Annals of Applied Statistics（査読中）
- **年**: 2022
- **DOI**: 未公開（受理済み：Received August 2021）
- **関連論点**: 論点2（Compositional data, Stick-breaking, Bayesian nonparametrics）

## 主な貢献
1. **Dependent Multivariate Bernstein Polynomial Process (DMBPP)**: Compositional responses（$\Delta_m$ simplexで支持されるデータ）のための完全ノンパラメトリック回帰モデルを提案
2. **Bernstein polynomialsとdependent stick-breaking processの統合**: Modified Bernstein polynomials（Barrientos et al. 2015）をdependent stick-breaking processに拡張
3. **Zero valuesへの対応**: Dirichlet distributionベースのモデルで、compositional dataがzero valuesを含む場合でも適切に定義される（log-ratio transformationの問題を回避）
4. **Spike-and-slab priorsによるmodel selection**: Predictor dependencyをweights、atoms、または両方のどこに導入すべきかを自動選択
5. **理論的性質の証明**: Continuity、association structure、support、posterior consistencyを確立

## 手法の概要

### Compositional dataの定義

$m$-dimensional simplex:
$$\Delta_m = \left\{ (y_1, \ldots, y_m) \in [0,1]^m : \sum_{i=1}^m y_i \leq 1 \right\}$$

観測データ: $\{(y_i, x_i): i=1,\ldots,n\}$
- $y_i \in \Delta_m$: Compositional response（例：固形廃棄物の組成比率）
- $x_i \in \mathscr{X} \subseteq \mathbb{R}^p$: Predictors（共変量）

### Modified Multivariate Bernstein Polynomials (MBP)

**Barrientos et al. (2015)のMBP**（本論文の基礎）:

Degree $k \in \mathbb{N}$のMBPは、CDF $G$に対して：
$$B(y | k, G) = \sum_{j \in \mathscr{H}_{k,m}} G\left(\frac{j_1}{k}, \ldots, \frac{j_m}{k}\right) \text{Mult}(j | k+m-1, y)$$

ここで：
- $\mathscr{H}_{k,m} = \{(j_1, \ldots, j_m) \in \{0,\ldots,k\}^m : \sum_{l=1}^m j_l \leq k+m-1\}$
- $\text{Mult}(\cdot | k+m-1, y)$: Multinomial pmf

**重要な性質**:
- $G$がCDFなら、$B(\cdot | k, G)$も有効なCDF
- Pointwise convergence: $B(\cdot | k, G) \to G|_{\Delta_m}$ as $k \to \infty$
- Uniform convergence: $G|_{\Delta_m}$が連続なら一様収束

**Density representation**（Mixture of Dirichlet）:

$G$が$\Delta_m^0 = \{y \in \Delta_m : y_j > 0, j=1,\ldots,m\}$上の確率測度のCDFなら：
$$b(y | k, G) = \sum_{j \in \mathscr{H}_{k,m}^0} G\left(\left(\frac{j_1-1}{k}, \frac{j_1}{k}\right] \times \cdots \times \left(\frac{j_m-1}{k}, \frac{j_m}{k}\right]\right) \text{dir}(y | \alpha(k, j))$$

ここで：
- $\mathscr{H}_{k,m}^0 = \{(j_1, \ldots, j_m) \in \{1,\ldots,k\}^m : \sum_{l=1}^m j_l \leq k+m-1\}$
- $\alpha(k, j) = (j, k+m-\|j\|_1)$（Dirichlet parameters）

### Dependent MBP Process (DMBPP)

**基本的なアイデア**: Mixing measure $G$をpredictor-dependent mixing measure $G_x$に置き換え。

**Random conditional densities**:
$$f_x(y | k, G_x) = \int_{\Delta_m} \text{dir}(y | \alpha(k, \lceil k\theta \rceil)) G_x(d\theta)$$

**Dependent stick-breaking representation**:
$$G_x(\cdot) = \sum_{j=1}^\infty w_j(x) \delta_{\theta_j(x)}(\cdot)$$

ここで：
$$w_j(x) = V_j(x) \prod_{l<j} [1 - V_l(x)]$$

### 3つのモデルバージョン

#### 1. General DMBPP（最も一般的）

**Definition 1**の完全版:
$$f_x(\cdot) = \sum_{j=1}^\infty w_j(x) \text{dir}(\cdot | \alpha(k, \lceil k\theta_j(x) \rceil))$$

ここで：
- $V_j(x) = v_x(\eta_j(x))$: Predictor-dependent weights
- $\theta_j(x) = h_x(z_j(x))$: Predictor-dependent atoms
- $\eta_j = \{\eta_j(x): x \in \mathscr{X}\}$: i.i.d. real-valued stochastic processes（law indexed by $\Psi_\eta$）
- $z_j = \{z_j(x): x \in \mathscr{X}\}$: i.i.d. real-valued stochastic processes（law indexed by $\Psi_z$）
- $v_x: \mathbb{R} \to [0,1]$, $h_x: \mathbb{R}^m \to \Delta_m^0$: Bijective continuous functions

**記号**: $\text{DMBPP}(\lambda, \Psi_\eta, \Psi_z, \mathscr{V}, \mathscr{H})$

#### 2. Single-weights DMBPP (wDMBPP)

Weightsが共通、atomsがpredictor-dependent:
$$f_x(\cdot) = \sum_{j=1}^\infty w_j \text{dir}(\cdot | \alpha(k, \lceil k\theta_j(x) \rceil))$$

- $v_j \sim$ i.i.d. $[0,1]$-valued random variables（indexed by $\Psi_v$）
- $w_j = v_j \prod_{l<j} (1 - v_l)$（共通のweights）

**記号**: $w\text{DMBPP}(\lambda, \Psi_v, \Psi_z, \mathscr{H})$

#### 3. Single-atoms DMBPP ($\theta$DMBPP)

Atomsが共通、weightsがpredictor-dependent:
$$f_x(\cdot) = \sum_{j=1}^\infty w_j(x) \text{dir}(\cdot | \alpha(k, \lceil k\theta_j \rceil))$$

- $\theta_j \sim$ i.i.d. $\Delta_m^0$-valued random vectors（indexed by $\Psi_\theta$）

**記号**: $\theta\text{DMBPP}(\lambda, \Psi_\eta, \mathscr{V}, \Psi_\theta)$

### Spike-and-slab Priorsによる Model Selection

**問題**: 3つのモデルバージョンのどれが最適か？

**解決策**: Binary indicators $(\gamma^\eta, \gamma^z)$を導入：
$$(\gamma^\eta, \gamma^z) \sim \pi_1 \delta_{(1,1)} + \pi_2 \delta_{(0,1)} + \pi_3 \delta_{(1,0)} + \pi_4 \delta_{(0,0)}$$

- $(1,1)$: General DMBPP（weights + atoms dependent）
- $(0,1)$: wDMBPP（atomsのみdependent）
- $(1,0)$: $\theta$DMBPP（weightsのみdependent）
- $(0,0)$: i.i.d. model（predictor independenceを表す退化ケース）

**Prior specifications**:
- **Prior I**: $\pi_1 = \pi_2 = \pi_3 = \pi_4 = 0.25$（等重み）
- **Prior II**: $\pi_2 = \pi_3 = 0.45$, $\pi_1 = \pi_4 = 0.05$（wDMBPPと$\theta$DMBPPを優遇）

### 実装

**Prior specifications**:
- $k | \lambda \sim \text{Poisson}(\lambda) \mathbb{I}_{\{k \geq 1\}}$（degree）
- $\lambda \sim \text{Gamma}(a_\lambda, b_\lambda)$
- $\eta_j(x) \sim \mathcal{N}(\mu_\eta(x), \sigma_\eta^2)$（weights process）
- $z_j(x) \sim \mathcal{N}_m(\mu_z(x), \Sigma_z)$（atoms process）
- Link functions: $v_x(a) = \Phi(a)$（probit link）, $h_x(b) = \text{softmax}(b)$

**MCMC**:
- Slice sampler for truncated stick-breaking（Kalli et al. 2011）
- Gibbs sampling for conditional posteriors
- Reversible jump MCMC for spike-and-slab indicators

## 実証結果

### Simulation Study（Section 6.1）

**4つのScenarios**:
- **Scenario I**: wDMBPP（atoms dependent、weights constant）
- **Scenario II**: wDMBPP（3 components）
- **Scenario III**: $\theta$DMBPP（weights dependent、atoms constant）
- **Scenario IV**: General DMBPP（both dependent）

**評価指標**:
- Integrated $L_1$ distance: $\hat{L}_1 = \int \int |\hat{f}(y|x) - f_0(y|x)| dy dx$
- Supremum $L_\infty$ distance: $\hat{L}_\infty = \max_i \max_j |\hat{f}(y_i|x_j) - f_0(y_i|x_j)|$

**結果（Table 2）**:

| Scenario | $n=250$ | $n=500$ | $n=1000$ |
|----------|---------|---------|----------|
| I        | 0.413   | 0.345   | 0.326    |
| II       | 0.479   | 0.426   | 0.411    |
| III      | 0.368   | 0.278   | 0.228    |
| IV       | 0.417   | 0.350   | 0.311    |

（Prior I使用時のintegrated $L_1$ distance）

- サンプルサイズ増加とともに精度向上
- Scenario IIIが最も良い性能（weightsのみdependentが最もシンプル）

**Model selection performance（Table 3）**:

| Scenario | $n=250$ | $n=500$ | $n=1000$ |
|----------|---------|---------|----------|
| I        | 1.000   | 1.000   | 1.000    |
| II       | 0.440   | 0.670   | 0.870    |
| III      | 0.500   | 0.960   | 0.990    |
| IV       | 0.460   | 0.530   | 0.530    |

（Prior Iでの真のモデル構造選択の正答率）

- Scenario I: 常に正しく選択（最もシンプル）
- Scenario II-III: サンプルサイズとともに改善
- Scenario IV: 複雑なモデルの識別は困難

### Application: Solid Waste Data from Colombia（Section 6.2）

**データ**:
- コロンビアの都市における住宅廃棄物の組成データ
- $n = 303$ households
- 5 components: Organic、Paper、Plastic、Others、Recyclables
- Predictor: Socioeconomic stratum（$x \in \{1,2,3\}$）

**主な発見**:
- Spike-and-slab priorはwDMBPP（atoms dependent）を選択
- Stratum 1（低所得）: Organic廃棄物の割合が高い
- Stratum 3（高所得）: Plasticの割合が高い
- モデルは条件付き密度の複雑な形状（multimodality）を捉えられる

## 限界・残された課題

### 著者が指摘する限界
1. **Computational cost**: Slice samplerと truncated stick-breakingの組み合わせは計算コストが高い（特に$m$が大きい場合）
2. **Link functions $v_x, h_x$の選択**: Probit linkとsoftmax transformationをdefaultとして使用しているが、他の選択肢（例：logit link）の影響は未検討
3. **Prior sensitivity**: $\lambda$（Poissonのparameter for $k$）やspike-and-slab prior $\{\pi_i\}$の選択にsensitiveの可能性
4. **Large $m$への拡張**: $m$が大きい（例：$m > 10$）場合のscalabilityが不明
5. **Zero-inflated data**: Zero valuesが大量に存在する場合（例：microbiome data）への対応は議論されているが、実証は限定的

### 本研究の視点からの限界
1. **Point pattern dataへの適用なし**: 本論文はregressionに焦点を置いており、marked point processへの直接適用は議論されていない
2. **Spatial dependenceの欠如**: Predictors $x$に空間位置が含まれる場合、空間相関の明示的モデル化がない
3. **Discrete vs continuous marks**: 本論文は continuous compositional marksを扱うが、discrete marksとの併用（例：obsidian source + composition）は未検討
4. **Computational efficiency for large $n$**: MCMC-basedで、大規模データ（$n > 10,000$）へのscalabilityは明示されていない

## 本研究(MMCP)との関係

### 直接的な関連
- **Compositional marksのモデリング**: 本研究の黒曜石組成データ（各産地からの黒曜石の化学組成）はcompositional dataであり、DMBPPの直接的な適用対象
- **Stick-breaking representation**: Linderman et al. (2015)のstick-breaking multinomialとの理論的関連性
- **Bayesian nonparametrics**: ノンパラメトリックベイズの枠組みでmarksの分布を柔軟にモデル化

### 本研究での具体的適用シナリオ

**Scenario 1: Obsidian composition as marks**
- Point pattern: 遺跡の空間分布（IPPまたはCox process）
- Mark: 各遺跡で発見された黒曜石の化学組成（$y \in \Delta_m$）
- Predictor: 各産地からの距離 $x = (d_1, d_2, \ldots, d_K)$

Conditional distribution:
$$f(y | s, d(s)) \sim \text{DMBPP}(\lambda, \Psi_\eta, \Psi_z, \mathscr{V}, \mathscr{H})$$

ここで$d(s) = (d_1(s), \ldots, d_K(s))$は遺跡位置$s$から各産地への距離。

**期待される効果**:
- 産地に近い遺跡：その産地由来の黒曜石組成が支配的
- 中間地点：複数産地の混合（mixture of Dirichlets）
- wDMBPPで組成の空間変化（atoms dependent）をモデル化

**Scenario 2: Joint modeling of discrete source + continuous composition**
- Discrete mark: Obsidian source $C \in \{$信州, 神津島, 箱根, 高原山$\}$（multinomial）
- Continuous mark: Composition $Y | C \in \Delta_m$（conditionally compositional）

Hierarchical model:
1. Source selection: $C | s \sim \text{Multinomial}(\pi(s))$（stick-breaking multinomial, Linderman et al. 2015）
2. Composition given source: $Y | C=c, s \sim f_c(y | s)$

ここで各$f_c(\cdot | s)$をDMBPPでモデル化。

### Wehrhahn et al. (2022) vs Eckardt et al. (2025)

本研究で両論文を統合する可能性:

| 側面 | Wehrhahn et al. (2022) | Eckardt et al. (2025) | 統合の可能性 |
|------|------------------------|----------------------|-------------|
| Focus | Regression（$y|x$） | Marked point process | 結合モデル |
| Compositional marks | ✓（Dirichlet mixtures） | ✓（Aitchison transformations） | 相補的 |
| Spatial dependence | × | ✓ | Eckardtから借用 |
| Stick-breaking | ✓（dependent） | × | Wehrhahnから借用 |
| Zero values | ✓（自然に対応） | △（alr/clr/ilrに問題） | Wehrhahnが優位 |

**統合モデルの構想**:
1. Spatial point process（Eckardt et al. 2025）: $N(A) \sim \text{Poisson}(\Lambda(A))$
2. Marks given location（Wehrhahn et al. 2022）: $Y | s \sim \text{DMBPP}(\lambda, \Psi_\eta, \Psi_z, \mathscr{V}, \mathscr{H})$
3. Intensity $\Lambda$とmark distribution $f_s(y)$の spatial covariance（Eckardtの枠組み）

### 本研究での限界として明記すべき点
第3章で引用する際、以下を明確に述べるべき：
- 本研究はpoint pattern dataであり、Wehrhahn et al. (2022)のregressionフレームワークへの直接適用は困難
- DMBPPの計算コストが高く、大規模な考古学データ（数千遺跡）へのscalabilityが懸念
- Spatial dependenceを明示的にモデル化するには、Eckardt et al. (2025)との統合が必要
- Zero-inflated compositional data（例：特定の産地由来が完全に欠如）への対応は理論的には可能だが、実装は未検証

## 引用すべき箇所

### Compositional dataの定義
> "From a mathematical point of view, compositional data can be defined as multivariate data supported on the $m$-dimensional simplex, $\Delta_m$..." (Introduction)

Compositional dataの数学的定義を明確にする引用。

### Log-ratio transformationの問題点
> "Although those approaches can be applied to compositional responses, by transforming the responses from $\Delta_m$ to $\mathbb{R}^m$, i) the resulting density in the simplex could not be well defined at the edges or ii) the resulting density in the simplex could be equal to zero at the edges. This can cause important problems if zeros are observed in the data..." (Introduction)

Log-ratio transformation（Aitchison変換）の問題点を指摘する引用（本研究での限界を議論する際に有用）。

### Bernstein polynomialsの近似性質
> "This class of MBP retains the appealing approximation properties of univariate BP and the standard class of MBP given in Equation (1). Specifically, if $G$ is a real-valued function defined on $\mathbb{R}^m$ and $G|_{\Delta_m}$ is its restriction on $\Delta_m$, then $B(\cdot | k, G)$ converges pointwise to $G|_{\Delta_m}$, as $k$ goes to infinity, and the relation holds uniformly on $\Delta_m$ if $G|_{\Delta_m}$ is a continuous function." (Section 2)

Bernstein polynomialsの理論的性質（一様収束）を説明する引用。

### Mixture of Dirichlet representation
> "Furthermore, if $G$ is the CDF of a probability measure defined on $\Delta_m^0$, then $B(\cdot | k, G)$ is the restriction of the CDF of a probability measure with density function given by the following mixture of Dirichlet distributions..." (Section 2, Equation 2)

Mixture of Dirichletとしての表現（本研究での実装の基礎）を説明する引用。

### Dependent stick-breaking process
> "To this end, we replace the mixing measure $G$ by a predictor-dependent mixing measure $G_x$. Under this approach, the random conditional densities are given by...where the set of mixing distributions $\{G_x : x \in \mathscr{X}\}$ follows a dependent stick-breaking process..." (Section 3)

Dependent stick-breaking processの導入を説明する引用（Linderman et al. 2015との関連）。

### Spike-and-slab priorsの動機
> "The use of the dependent stick-breaking process raises the question of where to introduce the predictor dependency: on weights, atoms, or both. Each selection leading to a different version of the model. Rather than fitting all versions of the model...we use spike-and-slab mixtures (George & McCulloch, 1993) to define a prior that automatically chooses the version of the model that best accommodates to the complexity of the underlying data-generating mechanism." (Introduction)

Model selectionの自動化（spike-and-slab）の動機を説明する引用。

### Simulation resultのまとめ
> "The model is able to choose the version of the model that is in agreement with the predictor dependency structure of the true model..." (Section 6.1)

Spike-and-slab priorsの有効性を示す実証結果の要約。

### Solid waste applicationの発見
> "The proposed approach is also applied for the analysis of solid waste data from Colombia." (Introduction)

Real data applicationの紹介（compositional dataの実例として有用）。

### Zero valuesへの対応
> "An important property of the considered model class is that the densities are well-defined in scenarios with compositional data containing zero values. When the response vector contains zero values, either models based on the Dirichlet distribution without restrictions on the parameter space or approaches based on the log-ratio transformation are not properly defined and cannot be employed unless a zero value imputation is applied first..." (Introduction)

Zero valuesの扱い（本論文の重要な利点）を説明する引用（本研究での適用可能性を論じる際に重要）。
