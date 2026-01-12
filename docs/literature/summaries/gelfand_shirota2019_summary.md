# Gelfand & Shirota (2019) - Presence/absenceとPresence-onlyデータの融合

## 基本情報
- **タイトル**: Preferential sampling for presence/absence data and for fusion of presence/absence data with presence-only data
- **著者**: Alan E. Gelfand, Shinichiro Shirota
- **ジャーナル**: Ecological Monographs
- **年**: 2019
- **DOI**: 10.1002/ecm.1372
- **関連論点**: 論点3（Presence-onlyデータ）、論点4（統合の困難さ）

## 主な貢献

1. **Presence/absenceとPresence-onlyの根本的な非互換性の指摘**
   - 従来研究が無視してきた2つのデータタイプの確率論的非互換性を明確化
   - "probability of presence"の定義が両者で異なることを詳細に論証

2. **点レベルPresence/absenceモデリングの提唱**
   - 領域単位（areal unit）ではなく点レベルでのモデリングの必要性
   - 真のpresence/absence表面は局所的に定数（locally constant）であるべき

3. **Preferential Samplingによる改善**
   - Shared process（$\eta(\mathbf{s})$または$\omega(\mathbf{s})$）を用いたpreferential samplingモデル
   - サンプリングバイアスを考慮した推論の改善

4. **確率論的に一貫したデータ融合の提案**
   - 従来の融合手法の問題点を指摘
   - Preferential samplingフレームワークによる一貫した融合

## 手法の概要

### 根本的な問題（Section 1.1, 4）

**Presence/absence（点レベル）**:
- $Y(\mathbf{s}) \sim \text{Bernoulli}(p(\mathbf{s}))$
- $p(\mathbf{s})$: 点$\mathbf{s}$でのpresence確率
- 領域$A$への拡張: $E(Y(A)) = \int_A p(\mathbf{s})d\mathbf{s}/|A|$（ランダム選択された点での確率）

**Presence-only（点過程）**:
- $N(A) \sim \text{Poisson}(\lambda(A))$、$\lambda(A) = \int_A \lambda(\mathbf{s})d\mathbf{s}$
- $P(\text{presence in } A) = P(N(A) \geq 1) = 1 - e^{-\lambda(A)}$
- 領域$A$のサイズに依存（点レベルに縮退しない）

> **非互換性**: 2つの定義は確率論的に異なる。$\lambda(A) \to 0$（$A \to 0$）で意味を失う。

### Presence/absenceの局所定常性（Section 4.1）

**パッチ概念**:
- 種の存在は有限個の「パッチ」（密集した個体群の領域）
- パッチ内の任意の点でpresence = 1、パッチ外でabsence = 0
- **結果**: 実現された表面は局所的に定数（locally constant）

**モデリングへの含意**:
- 条件付き独立Bernoulli試行は不適切（everywhere discontinuous）
- 潜在GP $Z(\mathbf{s})$を用いた第1段階モデル: $Y(\mathbf{s}) = 1(Z(\mathbf{s}) > 0)$
- $Z(\mathbf{s})$が平滑なら、$Y(\mathbf{s})$も局所定数

### Preferential Samplingモデル（Section 5）

**基本アイデア**: サンプリング位置$\mathcal{S}$とデータ$\mathcal{Y}$が独立か？

#### サンプリング位置$\mathcal{S}$のモデル

(i) **NHPP**: $\lambda(\mathbf{s}) = \mathbf{w}^T(\mathbf{s})\boldsymbol{\beta}$

(ii) **LGCP**: $\lambda(\mathbf{s}) = \mathbf{w}^T(\mathbf{s})\boldsymbol{\beta} + \eta(\mathbf{s})$

#### Presence/absenceデータ$\mathcal{Y}$のモデル

$Y(\mathbf{s}) = 1(Z(\mathbf{s}) > 0)$として：

(a) $Z(\mathbf{s}) = \mathbf{x}^T(\mathbf{s})\boldsymbol{\alpha} + \epsilon(\mathbf{s})$（空間回帰）

(b) $Z(\mathbf{s}) = \mathbf{x}^T(\mathbf{s})\boldsymbol{\alpha} + \omega(\mathbf{s}) + \epsilon(\mathbf{s})$（地球統計モデル）

(c) $Z(\mathbf{s}) = \mathbf{x}^T(\mathbf{s})\boldsymbol{\alpha} + \delta\eta(\mathbf{s}) + \omega(\mathbf{s}) + \epsilon(\mathbf{s})$（**Shared process**）

(d) $Z(\mathbf{s}) = \mathbf{x}^T(\mathbf{s})\boldsymbol{\alpha} + \delta\eta(\mathbf{s}) + \epsilon(\mathbf{s})$（Shared process、$\omega$なし）

**Preferentiality parameter $\delta$の解釈**:
- $\delta > 0$: presenceがoversampled（高いpresence確率領域でサンプリング）
- $\delta < 0$: presenceがundersampled
- $\delta = 0$: non-preferential sampling

**Shared processの役割**:
- $\eta(\mathbf{s})$が$\mathcal{S}$と$\mathcal{Y}$の両方に影響
- $[\mathcal{S} | \boldsymbol{\eta}_D] \neq [\mathcal{S}]$となり、確率的依存が生じる

### データ融合（Section 6）

**問題**: Presence-only ($\mathcal{S}_{PO}$)とPresence/absence ($\mathcal{S}_{PA}, \mathcal{Y}_{PA}$)をどう統合？

**従来手法の問題**:
- Presence-onlyの点過程フレームワークにPresence/absenceを無理やり押し込む
- 確率論的に非一貫

**提案**: 2つのshared processを用いた融合

$$
\begin{aligned}
Z(\mathbf{s}) &= \mathbf{x}^T(\mathbf{s})\boldsymbol{\alpha} + \delta_{PA}\eta_{PA}(\mathbf{s}) + \delta_{PO}\eta_{PO}(\mathbf{s}) + \omega(\mathbf{s}) + \epsilon(\mathbf{s}) \\
\mathcal{S}_{PA} &\sim \text{LGCP}(\lambda_{PA}(\mathbf{s}) = \mathbf{w}^T(\mathbf{s})\boldsymbol{\beta}_{PA} + \eta_{PA}(\mathbf{s})) \\
\mathcal{S}_{PO} &\sim \text{LGCP}(\lambda_{PO}(\mathbf{s}) = \mathbf{w}^T(\mathbf{s})\boldsymbol{\beta}_{PO} + \eta_{PO}(\mathbf{s}))
\end{aligned}
$$

**完全モデル**:
$$
[\mathcal{Y} | \mathcal{S}_{PA}, \boldsymbol{\alpha}, \boldsymbol{\eta}_{PA,Y}, \delta_{PA}, \boldsymbol{\eta}_{PO,Y}, \delta_{PO}][\mathcal{S}_{PA} | \boldsymbol{\beta}_{PA}, \boldsymbol{\eta}_{PA,D}][\mathcal{S}_{PO} | \boldsymbol{\beta}_{PO}, \boldsymbol{\eta}_{PO,D}]
$$

**期待**: $\delta_{PO} > 0$（presence-only位置は強くpresenceを示唆）

### Presence-onlyの劣化（Degradation）（Section 6.2）

**3つの表面**:
- **Potential intensity** $\lambda(\mathbf{s})$: 劣化がない場合
- **Availability** $U(\mathbf{s})$: 土地利用変化で利用可能か
- **Sampling effort** $T(\mathbf{s})$: 実際にサンプリングされたか

**Realized intensity**: $\lambda(\mathbf{s})U(\mathbf{s})T(\mathbf{s})$

グリッドセル$A_i$レベルで：
- $u_i = \int_{A_i} U(\mathbf{s})d\mathbf{s}/|A_i|$: 利用可能割合
- $q_i = \int_{A_i} T(\mathbf{s})U(\mathbf{s})d\mathbf{s}/|A_i|$: 利用可能かつサンプリングされた割合
- $p_i = q_i/u_i$: 利用可能な場所がサンプリングされた条件付き確率

## 実証結果

### データ
- **Presence/absence**: IPANE（4,314地点）、ニューイングランドの侵入植物7種
- **Presence-only**: GBIF、同じ7種
- **共変量**: WorldClimから7変数（気温・降水量関連）

### Preferential Samplingの効果（Section 5.3）

**モデル比較**（Tjur $R^2$、値が大きいほど良い）:
- モデル(a)（空間回帰）< モデル(d)（$\delta\eta_{PA}$追加）
- モデル(b)（$\omega$追加）が最良
- モデル(c)（$\omega + \delta\eta_{PA}$）はモデル(b)と同等以上

**$\delta$の推定**:
- 4種（MR, JB, GB, AO）で$\delta$が0から有意に離れる
- $\delta > 0$: presenceがoversampledされている証拠

### データ融合の効果（Section 6.3）

**モデル(e), (f)の結果**:
- $\delta_{PO}$は全種で0から大きく離れる（強いpreferentiality）
- $\delta_{PA}$は有意でなくなる（$\delta_{PO}$が主要）
- Presence-only位置情報がPresence/absence推論を改善

## 限界・残された課題

### 著者が述べる限界

1. **計算コスト**
   - NNGP（$k=15$近傍）を使用したが、依然として計算負荷が高い
   - MHステップ（$\eta$のサンプリング）が時間を要する

2. **識別可能性**
   - 3つのGP（$\omega, \eta_{PA}, \eta_{PO}$）を同時に扱うと識別問題
   - 強い事前情報（range制約等）が必要

3. **将来の拡張**
   - 複数種の同時分布モデル（Joint species distribution models）
   - 個体数データ（abundance）への拡張

### 本研究の観点からの限界

1. **マークの未対応**
   - 点の位置のみでマーク（組成等）は扱っていない

2. **Pólya-Gamma拡張の未使用**
   - ロジスティック回帰部分はMHサンプリング
   - Pólya-Gammaによる効率化は適用されていない

3. **空間変化係数の未導入**
   - 共変量効果は空間的に一定

## 本研究（MMCP）との関係

### 借用する要素

1. **点レベルモデリングの哲学**
   - Presence/absenceを点レベルで扱う
   - 局所定常性の考え方

2. **Shared processによるpreferential sampling**
   - $\delta\eta(\mathbf{s})$による観測バイアスの考慮
   - データタイプ間の依存構造のモデリング

3. **LGCP + NNGPの計算効率化**
   - 点過程のベイズ推論での空間依存構造

4. **データ融合の概念**
   - 異なるデータタイプの確率論的に一貫した統合

### MMCPによる拡張

1. **組成マークの追加**
   - Gelfand & Shirotaは点の位置のみ
   - MMCPは各点に組成（黒曜石産地比率）を付与

2. **Pólya-Gamma拡張**
   - multinomial logitリンクをGibbsサンプリング
   - 計算効率の向上

3. **空間変化係数**
   - 共変量効果が空間的に変化するモデル

### 位置づけ

Gelfand & Shirota (2019)は：
- Presence/absenceとPresence-onlyの確率論的非互換性を明確化
- Preferential samplingによる改善と融合の枠組みを提供

MMCPは：
- この枠組みを継承しつつ、組成マークを追加
- Presence-only点過程における組成分布のモデリング
- ただしMMCPはpresence/absenceとの融合は範囲外（presence-onlyのみ）

## 引用すべき箇所

### 非互換性の主張（Abstract）
> "We illuminate the fundamental modeling differences between the two types of data. Most simply, locations are considered as fixed under presence/absence data; locations are random under presence-only data. The definition of 'probability of presence' is incompatible between the two."

### 点レベルの必要性（Section 1.1）
> "This would seem to be at odds with how presence/absence arises in nature. In fact, it argues that, if presence/absence is to be viewed coherently, it must be viewed at point level."

### 局所定常性（Section 4.1）
> "Second, presence/absence is a neighborhood phenomenon. If there is a presence at $\mathbf{s}$ then there is presence everywhere in a sufficiently small neighborhood, $\partial\mathbf{s}$, of $\mathbf{s}$. [...] As a result, the realized presence absence surface is locally constant."

### Preferential samplingの定義（Section 5.2）
> "Working with (b) and (iii), if $\psi = 0$, then, following Diggle et al. (2010), we have non-preferential sampling while if $\psi \neq 0$, we have strong preferential sampling."

### データ融合の問題（Section 6）
> "Since we argue that a point pattern specification is inappropriate for presence/absence data, a different type of fusion is required. We have a point pattern model for the presence-only data and a binary map model for the presence/absence data."

### 将来の方向性（Section 7）
> "Future work offers much opportunity. More experience is needed with regard to the rich set of modeling specifications that we have presented in Sections 5 and 6."
