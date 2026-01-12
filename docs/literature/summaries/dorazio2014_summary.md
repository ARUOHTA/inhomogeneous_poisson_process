# Dorazio (2014) - Presence-Only Dataにおける不完全検出と調査バイアスの統計的対処

## 基本情報
- **タイトル**: Accounting for imperfect detection and survey bias in statistical analysis of presence-only data
- **著者**: Robert M. Dorazio (US Geological Survey)
- **ジャーナル**: Global Ecology and Biogeography
- **年**: 2014
- **DOI**: 未記載
- **関連論点**: 論点3（Presence-only data modeling）

## 主な貢献
1. **不完全検出と調査バイアスの明示的モデル化**: 従来のpresence-only data分析が無視していた、不完全な検出（imperfect detection）と調査バイアス（survey bias：調査地点の偏った選択）を考慮した階層モデルを提案
2. **データ統合フレームワーク**: Presence-onlyデータ（opportunistic surveys）とplanned surveysのデータ（例：replicated point counts）を統一的にモデル化し、少数の高品質データで大規模presence-onlyデータを「leverage」
3. **パラメータ識別可能性の理論的検討**: Fisher information matrixを用いて、presence-onlyデータ単独ではSDMパラメータが識別不可能になる条件を明示
4. **シミュレーション実証**: わずか50 quadrats（0.5%のサンプリング）のplanned surveysを追加するだけで、推定バイアスと不確実性が劇的に削減されることを実証

## 手法の概要

### 階層モデル構造

**Component 1: Spatial point process model**
個体の活動中心（activity centers）がPoisson point processに従うと仮定：
$$\lambda(s) = \exp(\beta^T x(s))$$

領域$B$内の個体数：
$$N(B) \sim \text{Poisson}(\mu(B)), \quad \mu(B) = \int_B \lambda(s) ds$$

**Component 2a: Detections in opportunistic surveys**
検出確率$p(s)$による位置依存的thinning：
$$\text{logit}(p(s)) = \alpha^T w(s)$$

検出された個体のpoint processは：
$$M(B) \sim \text{Poisson}(\nu(B)), \quad \nu(B) = \int_B \lambda(s) p(s) ds$$

Likelihood（presence-only data）：
$$L(\beta, \alpha) = \frac{\exp(-\nu(B))}{m!} \prod_{i=1}^m \lambda(s_i) p(s_i)$$

**Component 2b: Detections in planned surveys**
Sample unit $C_k$内の個体数：
$$N(C_k) \sim \text{Poisson}(\mu(C_k))$$

$J_k$回の独立したpoint count surveys（product-binomial model）：
$$\Pr(Y_k = y_k | N(C_k) = n_k) = \prod_{j=1}^{J_k} \binom{n_k}{y_{kj}} p_{kj}^{y_{kj}} (1 - p_{kj})^{n_k - y_{kj}}$$

ここで$p_{kj}$は検出確率：
$$\text{logit}(p_{kj}) = \gamma^T v(C_k)$$

Unconditional likelihood（point counts）：
$$\Pr(Y_k = y_k) = \sum_{n_k = \max(y_k)}^{\infty} \frac{\exp(-\mu(C_k)) (\mu(C_k))^{n_k}}{n_k!} \prod_{j=1}^{J_k} \binom{n_k}{y_{kj}} p_{kj}^{y_{kj}} (1 - p_{kj})^{n_k - y_{kj}}$$

### 統合likelihood

Opportunistic surveysとplanned surveysは独立に取得されるため：
$$L(\beta, \alpha, \gamma) = L(\beta, \alpha) \times L(\beta, \gamma)$$

### パラメータ識別可能性

**Presence-only dataのみの場合の問題**:

1. **検出確率が一定の場合**: $p(s) = p$（定数）とすると、$\beta_0$と$\log(p)$は和としてのみ識別可能：
   $$\log L(\beta, p) = -\exp(\beta_0 + \log(p)) \int_B \exp(\tilde{\beta}^T \tilde{x}(s)) ds + m(\beta_0 + \log(p)) + \sum_{i=1}^m \tilde{\beta}^T \tilde{x}(s_i)$$

2. **検出確率が低い場合**: $p(s) < 0.2$のとき、
   $$\lambda(s) p(s) \approx \exp(\beta^T x(s) + \alpha^T w(s))$$
   となり、$x(s)$と$w(s)$が線形従属だとパラメータ識別不可能

**識別可能性の診断**:
Fisher information matrix $I(\theta)$（$\theta = (\beta^T, \alpha^T)^T$）の条件数（最大固有値/最小固有値）を計算。
- $I(\theta)$が正定値（full rank）ならパラメータ識別可能
- 条件数の逆数（最小/最大固有値比）が0に近いと識別困難

**必要条件**:
- $\lambda(s)$の共変量と$p(s)$の共変量が線形独立であること
- 検出確率が全ての場所でextreme（ほぼ0またはほぼ1）でないこと

### データ統合による改善

Planned surveysを追加すると：
- $J_k \geq 2$の反復調査があれば、abundanceとdetectabilityが分離可能（Royle 2004）
- $J_k = 1$でも、共変量が線形独立なら識別可能（Sólymos et al. 2012）

## 実証結果

### シミュレーション設定
- 領域$B$: 正方形領域、10,000 quadratsに分割
- 真のパラメータ: $\log(\lambda(s)) = \log(8000) + 0.5 x(s)$
- 検出モデル: $\text{logit}(p(s)) = \alpha_0 - 1.0 w(s)$（または$x(s)$）
- Planned surveys: $K = 50, 100, 200, 400, 800$ quadrats, $J = 4$ replicates
- 期待値: $E(N(B)) = 35,857$, $E(M(B)) \approx 7,922$～$11,251$

### 結果1: 独立共変量の場合（$w(s) \perp x(s)$）

Presence-onlyのみの推定:
- $\hat{\beta}$のバイアスは無視できるほど小さい
- ただし不確実性は大きい

Presence-only + planned surveys:
- $\hat{\beta_0}$の不確実性が大幅に削減
- $\hat{\beta_1}$の不確実性は変わらず（Figure 5）

### 結果2: 同一共変量の場合（$w(s) = x(s)$）

Presence-onlyのみの推定:
- $\hat{\beta}$に**強いバイアス**と**高い分散**
- パラメータ識別不可能（Fisher information matrixの条件数が極端）

Presence-only + planned surveys:
- わずか50 quadrats（0.5%）を追加するだけで**バイアスと不確実性が劇的に削減**
- 200 quadrats以上でほぼバイアスなし（Figure 6）

### 重要な発見
- $J = 4$ replicatesが必要（$J = 1$では改善なし）
- Planned surveysのサンプルサイズは小さくて良い（全体の0.5～8%）
- 高品質データ（abundanceとdetectabilityの両方に情報的）が少数あればpresence-onlyデータを有効活用できる

## 限界・残された課題

### 著者が指摘する限界
1. **Spatial autocorrelationの無視**: Poisson processは空間独立性を仮定しており、実際のデータに存在する空間相関を考慮していない
2. **個体の移動**: 活動中心での検出を仮定しており、高度に移動性の種や大きなテリトリーを持つ種には不適
3. **Planned surveysの必要性**: パラメータ識別のためにはplanned surveysが実質的に必須
4. **Observer heterogeneityの簡略化**: 観測者間の検出能力の違いを単一のinterceptで処理（多数の観測者がいる場合の限界）

### 本研究の視点からの限界
1. **頻度論的アプローチ**: Bayesian frameworkではないため、事前情報の活用や不確実性の完全な定量化が困難
2. **Spatial random effectsなし**: Log Gaussian Cox Processのような空間的潜在構造を持たない
3. **Marksの考慮なし**: Marked point processへの拡張がない
4. **Computational efficiency**: 大規模データへの適用可能性が不明（NNGPのような近似手法の議論なし）
5. **実データ適用の限界**: 実際にplanned surveysを実施することは費用と労力がかかり、多くのpresence-only dataには適用困難

## 本研究(MMCP)との関係

### 直接的な関連
- **Presence-only dataの扱い**: 本研究も黒曜石遺跡データという形でpresence-onlyデータを扱う
- **階層モデル構造**: 点過程モデル（真の分布）+ 観測プロセス（検出メカニズム）という二段階構造
- **パラメータ識別の課題**: 観測プロセスと真の分布の分離が困難という問題意識を共有

### 本研究での課題
- **Planned surveysの欠如**: 本研究では考古学的遺跡データのみで、planned surveysに相当するデータがない
- **検出バイアスの扱い**: Dorazioのモデルは検出確率$p(s)$を明示的にモデル化するが、本研究では完全な確率的観測（$p(s) = 1$または調査範囲内での完全検出）を暗黙的に仮定している可能性
- **時間的側面**: 考古学データは時代を跨ぐが、検出プロセスは現代の調査努力に依存する（時代による検出バイアス）

### 本研究で借用する要素
- **Thinned point processの考え方**: 観測された遺跡 = 真の分布 × 検出プロセス、という概念的枠組み
- **空間的検出バイアス**: 調査努力が高い地域（例：都市近郊、アクセスの良い地域）でより多くの遺跡が発見される可能性
- **Fisher information matrixによる診断**: パラメータ識別可能性の事前チェック（特にBayesian推定でも有用）

### 本研究で拡張する点
- **Bayesian framework**: Dorazioの頻度論的モデルをBayesian IPP/LGCPに拡張
- **Spatial random effects**: NNGPやVecchia approximationによる空間相関の導入
- **Marked point process**: Obsidian sourceというdiscrete marksの追加
- **Compositional marks**: さらにEckardt et al. (2025)のように組成値marksへの拡張
- **Planned surveysなしでの推定**: Dorazioが「必須」とするplanned surveysなしで、どこまで推定可能かを検討（限界を明示的に議論）

### 本研究の限界として明記すべき点
第3章で引用する際、以下を明確に述べるべき：
- 本研究は考古学的調査努力の空間的偏りを明示的にモデル化していない
- Dorazio (2014)が示すように、検出バイアスを無視すると推定にバイアスが生じる可能性
- Planned surveys相当のデータがないため、パラメータ識別可能性に課題が残る
- 将来の拡張として、調査努力データ（例：各地域の発掘調査回数）を共変量として組み込む必要性

## 引用すべき箇所

### 不完全検出の重要性
> "Nearly all surveys of natural populations, including opportunistic surveys that produce presence-only observations, are prone to detection errors." (Introduction)

Presence-only dataにおける不完全検出の普遍性を述べる際の引用。

### 検出バイアスがSDM推定に与える影響
> "Failure to account for imperfect detectability in models of presence-only data induces bias in estimates of SDMs when the covariates of abundance are not distinct and stochastically independent of the covariates of detectability." (Introduction)

検出バイアスを無視した場合の問題を説明する引用。

### パラメータ識別不可能性の条件
> "If the probability of detecting an individual is identical at all locations, the parameters of a SDM are not identifiable." (Discussion)

検出確率が一定の場合の識別問題。

> "If the detection probability differs among locations, the parameters of a SDM can be estimated if some of the regressors used to specify spatial differences in individual density are linearly independent of the regressors of detection probability." (Discussion)

識別可能性のための必要条件。

### データ統合の利点
> "A relatively small number of high-quality data (from planned surveys) can be used to leverage the information in presence-only observations, which usually have broad spatial coverage but may not be informative of both occurrence and detectability of individuals." (Abstract)

Planned surveysとpresence-onlyデータの統合の効果を述べる引用（本研究の限界を議論する際、このアプローチが取れていないことを明記する文脈で使用）。

### Point processモデルの利点
> "The parameters of a point-process model are invariant to spatial scale, so the model can be used to predict the abundance or occurrence of individuals for any subregion located within the region of interest." (Introduction)

Point processモデルの空間スケール不変性の利点。

### Occurrence probabilityの定義
> "$$\Pr(N(C) > 0) = 1 - \exp(-\mu(C))$$
> where $\mu(C) = \int_C \lambda(s) ds$ is a function of the model's parameters $\beta$ and the values of the spatially varying regressors $x$ in $C$." (Discussion)

Occurrence probabilityと intensityの関係（本研究でも同様の定義を使用）。

### Site occupancy modelの限界
> "Conventional site-occupancy models (MacKenzie et al., 2002; Tyre et al., 2003) are limited for use as SDMs, i.e. the interpretation of occurrence probability in site-occupancy models depends on the spatial resolution used in the analysis." (Discussion)

Site occupancy modelの空間解像度依存性の問題（point processモデルの優位性を論じる際の引用）。

### 本研究の限界として引用可能な箇所
> "One limitation of the point-process model of Warton & Shepherd (2010) is that it fails to account for the effects of errors in the detection of individuals." (Introduction)

本研究もWarton & Shepherd (2010)と同様に検出エラーを考慮していないことを認める際の引用。
