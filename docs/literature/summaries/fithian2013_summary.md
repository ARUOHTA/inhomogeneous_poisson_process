# Fithian & Hastie (2013) - Presence-Only Data統計モデルの有限サンプル等価性

## 基本情報
- **タイトル**: Statistical Models for Presence-Only Data: Finite-Sample Equivalence and Addressing Observer Bias
- **著者**: William Fithian, Trevor Hastie (Stanford University)
- **ジャーナル**: Annals of Applied Statistics
- **年**: 2013
- **DOI**: 10.1214/13-AOAS667
- **関連論点**: 論点3（Presence-only data modeling）

## 主な貢献
1. **有限サンプルでの等価性証明**: Inhomogeneous Poisson Process (IPP)、Maxent、infinitely weighted logistic regressionが有限サンプルで完全に等価であることを示した
2. **Infinitely weighted logistic regression**: IPP推定値を復元する新しい手法を提案し、標準的な統計ソフトウェアで実装可能にした
3. **Observer bias対処**: Thinned IPP modelを用いて、観測バイアス（調査努力の空間的偏り）を明示的にモデル化
4. **データ統合**: Presence-onlyデータとpresence-absenceデータを統一的に扱う枠組みを提示

## 手法の概要

### IPP、Maxent、logistic regressionの等価性

**Occurrence rate vs occurrence probability**:
- Occurrence rate $\lambda(z)$: 単位面積あたりの期待出現数（無制限）
- Occurrence probability $\theta(z)$: グリッドセル内での出現確率（0-1に制約）
- 関係式: $\theta(z) \approx 1 - e^{-a\lambda(z)} \approx a\lambda(z)$（セルサイズ$a$が小さい時）

**IPPモデル**:
$$\lambda(z) = \lambda_0 e^{\beta^T z}$$

対数尤度の数値近似（quadrature points $z_1, \ldots, z_m$を使用）:
$$\ell(\beta) \approx \sum_{i=1}^n \beta^T z_i - \sum_{j=1}^m a_j \lambda_0 e^{\beta^T z_j}$$

**Maxentとの関係**:
- Maxentは制約付き最大エントロピー分布として定式化される
- 双対問題を解くと、IPPの対数尤度推定と同一の解が得られる（Proposition 1）

**Infinitely weighted logistic regression**:
- Background points（quadrature points）に weight $W \to \infty$ を割り当てる
- Proposition 2: $\lim_{W \to \infty} \hat{\beta}_W = \hat{\beta}_{\text{IPP}}$
- 実装上は $W = 10^4$ 程度で十分な近似が得られる

### Observer biasのモデル化

**Thinned IPP**:
観測プロセスを2段階で考える：
1. 真の出現プロセス: $\tilde{\lambda}(z)$ （潜在的な分布）
2. 観測プロセス: $\lambda(z) = \tilde{\lambda}(z) \cdot s(z)$

ここで $s(z)$ は観測確率（sampling effort）を表す。

**推定手順**:
- $s(z)$ を別データ（調査努力データ）から推定、またはcovariateとしてモデル化
- $\tilde{\lambda}(z) = \lambda(z) / s(z)$ を復元
- これにより調査努力が低い地域での真の分布を推定可能

### Presence-onlyとpresence-absenceの統合

Case-control sampling interpretation:
- Presence-only: cases from $P(z | \text{presence})$
- Presence-absence: cases and controls from $P(z | \text{presence})$ and $P(z | \text{absence})$

統合推定量:
$$\ell_{\text{pooled}}(\beta) = \ell_{\text{PO}}(\beta) + \ell_{\text{PA}}(\beta)$$

複数種のデータプーリング:
$$\ell_{\text{all}}(\beta) = \sum_{k=1}^K \left[ \sum_{i=1}^{n_k} \beta^T z_{ki} - \lambda_{0k} \sum_{j=1}^m a_j e^{\beta^T z_j} \right]$$

## 実証結果

### Bradypus variegatus（ミユビナマケモノ）データ分析
- 116 presence-only records（GBIF）
- 6 environmental covariates
- Maxent vs infinitely weighted logistic regression: 推定値がほぼ一致（Table 2）
- Cross-validationによるモデル選択でも同一のモデルを選択

### シミュレーション研究
- 異なる重み $W$ での収束を確認
- $W = 10^4$ で IPP推定値とほぼ同一の結果
- Background pointsの数 $m$ を増やすと近似精度向上

### 複数種のデータ統合
- Willow Warblerのpresence-onlyデータ（約20,000点）
- 他の鳥種のpresence-absenceデータと統合
- 共通のhabitatパラメータ $\beta$ を推定し、種間での環境選好の違いを定量化

## 限界・残された課題

### 著者が指摘する限界
1. **Spatial autocorrelationの無視**: IPPモデルは空間独立性を仮定しており、実際のデータにある空間相関を考慮していない
2. **Observer biasの推定困難**: $s(z)$ を直接観測することは困難で、別途調査努力データが必要
3. **Background pointsの選択**: Quadrature pointsの数と配置が推定精度に影響するが、最適選択の理論は未発達
4. **計算コスト**: Background pointsを多数使用すると計算量が増大

### 本研究の視点からの限界
1. **Bayesian frameworkの欠如**: 頻度論的アプローチのみで、事後分布や不確実性の完全な定量化が困難
2. **Spatial random effectsなし**: Log Gaussian Cox Processのような空間的潜在構造を持たない
3. **Marksの考慮なし**: 点パターンのみでmarked point processへの拡張がない
4. **Static model**: 時空間動態や共変量の時間変化を扱えない

## 本研究(MMCP)との関係

### 直接的な関連
- **Presence-only data framework**: 本研究もpresence-only dataを扱うため、IPP-logistic regression等価性は基礎理論として重要
- **Observer bias認識**: Thinned IPP modelの考え方は、本研究の観測プロセスモデリングに示唆を与える
- **Quadrature scheme**: Background pointsの使用は、本研究の数値計算でも採用している

### 本研究での拡張
- **Bayesian IPP**: Fithianらの頻度論的IPPをBayesian frameworkに拡張（Moreira et al. 2022の流れ）
- **Spatial random effects**: Log Gaussian Cox Processによる空間相関の導入
- **Marked process**: Obsidian sourceという離散marksの追加
- **Compositional marks**: さらにEckardt et al. (2025)のように組成値marksへ拡張

### 本研究で借用する要素
- **Infinitely weighted logistic regression**: 計算上の利便性から、同様の近似手法を検討可能
- **Quadrature points**: 積分近似のためのbackground pointsの配置戦略
- **Case-control interpretation**: データ拡張におけるpresence/absenceの解釈

### 本研究で拡張する点
- Spatial GPによる $\lambda(s)$ の柔軟なモデル化
- NNGPによる大規模データへの対応
- Pólya-Gamma augmentationによる効率的な事後サンプリング
- Marked point process（特に組成値marks）への拡張

## 引用すべき箇所

### 理論的背景（IPP-logistic regression等価性）
> "We show that under certain conditions, presence-only data analysis using inhomogeneous Poisson process (IPP) models, maximum entropy (Maxent) models, and logistic regression models all lead to the same parameter estimates." (Abstract)

IPP、Maxent、logistic regressionの等価性を説明する際の基礎引用。

### Occurrence rate vs probability の区別
> "It is important to distinguish between the occurrence rate λ(z) (also called intensity), which can be arbitrarily large, and the occurrence probability θ(z), which is constrained to lie in [0,1]." (Section 2.1)

Presence-only dataのモデル化において、occurrence rateとoccurrence probabilityの違いを明確にする際の引用。

### Infinitely weighted logistic regression の提案
> "We propose a simple modification of the weighted logistic regression approach: assign infinitely large weights to the background points. This yields parameter estimates identical to those from IPP likelihood maximization." (Section 2.3)

実装上の工夫を説明する際の引用。

### Observer bias の定式化
> "We model the observed point pattern as a thinned version of the true occurrence process: λ(z) = λ̃(z)s(z), where λ̃(z) is the true intensity and s(z) is the probability of observing a specimen at location z." (Section 4)

Observer biasの明示的なモデル化の必要性を論じる際の引用（本研究では完全な確率的観測を仮定しているが、限界として言及可能）。

### データ統合の枠組み
> "Presence-only and presence-absence data can be combined in a unified likelihood framework, yielding more efficient parameter estimates when both data types are available." (Section 5)

複数のデータソース統合の重要性を論じる際の引用。

### 関連研究レビュー
Section 1.1 "Related Work" に、species distribution modeling、Maxent、point process modelsの文献レビューが簡潔にまとめられており、関連研究の位置づけに有用。
