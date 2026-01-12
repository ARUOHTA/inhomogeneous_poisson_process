# Renner et al. (2015) - Point Process Models for Presence-only Analysis

## 基本情報
- **タイトル**: Point process models for presence-only analysis
- **著者**: Ian W. Renner, Jane Elith, Adrian Baddeley, William Fithian, Trevor Hastie, Steven J. Phillips, Gordana Popovic, David I. Warton
- **所属**: 複数機関（Newcastle, Melbourne, Curtin, Stanford, NSW等）
- **ジャーナル**: Methods in Ecology and Evolution
- **年**: 2015
- **DOI**: 確認中
- **関連論点**: 論点3（Presence-onlyデータ）

## 主な貢献

1. **Presence-only dataに対する点過程モデル（PPM）のレビュー**
   - 種分布モデリング（Species Distribution Modeling: SDM）の文脈で点過程モデルの理論と実践を体系化
   - 複数のアプローチ（Maxent, GLM, spatstat等）の統一的理解

2. **Quadrature points（背景点）選択の客観的枠組み**
   - 従来の「pseudo-absences」「background points」を数値積分の問題として定式化
   - 選択基準とモデル安定性のチェック手法

3. **モデル仮定の明確化とチェック手法**
   - 応答変数が何か（強度関数）の明確化
   - 空間独立性の仮定チェック
   - 空間依存性への対処（Gibbs, Cox過程）

4. **サンプリングバイアスへの対処**
   - Observer biasを共変量として組み込む方法
   - 観測努力（sampling effort）のモデリング

5. **複数ソフトウェアツールの実用的比較**
   - spatstat, ppmlasso, Maxent, R-INLA, lgCP等の比較
   - 各ツールの長所・短所の整理

## 手法の概要

### Presence-only dataと点過程

**Presence-only data**: 種が観測された地点の集合$\{s_1, \ldots, s_n\}$（absenceデータなし）

**点過程モデル（PPM）**: 点の位置$s \in \mathcal{D}$（$\mathcal{D}$: 研究領域）をランダム事象として扱う

**強度関数** $\lambda(s)$: 微小領域$ds$あたりの期待点数
$$
\mathbb{E}[\text{# points in } A] = \int_A \lambda(s) ds
$$

**重要な注意**: $\lambda(s)$は通常、種の**真の豊度**（abundance）を反映せず、**報告の期待数**（expected abundance of species reportings）を反映。相対パターンの推論のみ可能（Fithian & Hastie 2013）。

### Poisson点過程（空間独立性を仮定）

**仮定**:
(a) 点の位置は互いに独立
(b) 強度関数$\lambda(s)$が環境共変量$\mathbf{x}(s)$の関数として表現可能

**対数線形モデル**:
$$
\lambda(s) = \exp\{\beta_0 + \mathbf{x}(s)^T\boldsymbol{\beta}\}
$$

### 対数尤度関数

**Poisson点過程の対数尤度**（Cressie 1993）:
$$
\log L(\boldsymbol{\beta}) = \sum_{i=1}^n \log\lambda(s_i) - \int_\mathcal{D} \lambda(s) ds
$$

**積分項の解釈**: 研究領域$\mathcal{D}$全体での期待点数

**数値積分による近似**（Quadrature points $\{q_1, \ldots, q_m\}$と重み$\{w_1, \ldots, w_m\}$）:
$$
\int_\mathcal{D} \lambda(s) ds \approx \sum_{j=1}^m w_j \lambda(q_j)
$$

### Quadrature pointsの選択

**従来の問題設定**: 「Pseudo-absences」や「background points」をどう選ぶか？

**点過程の視点**: 数値積分の問題として定式化
- 目的: 尤度（積分）の正確な推定
- 基準: モデルが安定し、Quadrature pointsの再サンプリングに対して不変

**推奨手順**:
1. 十分な数のQuadrature pointsを選択
2. Quadrature pointsを変えてモデルを再フィッティング
3. パラメータ推定値が安定していることを確認
4. 不十分な場合は数を増やす

**典型的な選択**:
- ランダムサンプリング（一様分布）
- 格子点
- 共変量空間での層別サンプリング

### Presence-background (PB) 回帰モデルとの関係

**PB回帰**: Presence点を$y_i = 1$、Background点を$y_i = 0$としてロジスティック回帰

**点過程尤度との関係**: Quadrature近似を用いると、点過程の対数尤度はPB回帰の対数尤度に近似的に等しい（Warton & Shepherd 2010）

**結論**: GLMソフトウェアで点過程モデルをフィッティング可能

### Maxentとの関係

**Maxent**: 最大エントロピー原理に基づく種分布モデリング手法（Phillips et al. 2006）

**等価性**: Maxentは特定の正則化を持つPoisson点過程モデルと等価（Fithian & Hastie 2013）

**結論**: Maxentソフトウェアも点過程モデルをフィッティングしている

### 空間依存性への対処

**問題**: Presence-onlyデータはしばしば空間的にクラスター化（独立性の仮定違反）

**2つの主なアプローチ**:

#### 1. Gibbs (Interaction) Processes

**基本アイデア**: 点間の相互作用を明示的にモデリング

**例: Area-interaction process**（Widom & Rowlinson 1970; Baddeley & van Lieshout 1995）:
- 距離$2r$以内の全ての点間に相互作用
- 強度関数に追加項を導入（観測点を条件とする）

**解釈**: 点間の誘引（clustering）または反発（regularity）をモデリング

#### 2. Cox Processes

**Log-Gaussian Cox Process (LGCP)**（Møller et al. 1998）:

**定義**: 強度関数自体がランダム
$$
\lambda(s) = \exp\{\beta_0 + \mathbf{x}(s)^T\boldsymbol{\beta} + \xi(s)\}
$$
ここで$\xi(s)$は平均0のGaussian process。

**解釈**: GLMMのランダム切片の点過程版

**利点**:
- 測定されていない共変量の効果を捉える
- 空間的クラスタリングを自然にモデリング
- ベイズ推論フレームワーク（INLA, MCMC）

### フィッティング手法とソフトウェア

**複数のソフトウェアツール**（Table 1で比較）:

1. **spatstat** (Rパッケージ):
   - Poisson PPMとGibbs processesに対応
   - モデル診断ツールが充実
   - Quadrature pointsの選択が柔軟

2. **ppmlasso** (Rパッケージ):
   - LASSO正則化を用いたPPM
   - 変数選択が可能
   - Area-interaction modelsに対応

3. **iwlr / dwpr**:
   - 重み付けロジスティック回帰
   - GLMソフトウェアで実装可能
   - 計算が高速

4. **Maxent**:
   - 広く使われている
   - 正則化が自動
   - モデル仮定のチェック機能が限定的

5. **R-INLA** (Rパッケージ):
   - LGCPのフィッティング
   - 高速な近似ベイズ推論
   - 空間依存性を考慮

6. **lgCP** (Rパッケージ):
   - LGCPのMCMC推論
   - 完全なベイズ推論
   - 計算コストが高い

### サンプリングバイアスへの対処

**問題**: 観測努力が空間的に不均一（道路沿い、都市近郊に偏る等）

**解決策**: バイアス変数を共変量として組み込む

**例**:
- 道路からの距離
- 都市からの距離
- 人口密度
- 観測努力の代理変数

**モデル**:
$$
\lambda(s) = \exp\{\underbrace{\mathbf{x}_{\text{env}}(s)^T\boldsymbol{\beta}_{\text{env}}}_{\text{環境効果}} + \underbrace{\mathbf{x}_{\text{bias}}(s)^T\boldsymbol{\beta}_{\text{bias}}}_{\text{バイアス効果}}\}
$$

## 実証例: Eucalyptus sparsifolia

### データ
- 230 presence点（オーストラリア・Blue Mountains地域）
- 環境共変量: 降水量、気温、標高等

### モデル比較（Figure 5）

1. **Poisson PPM** (spatstat / Maxent):
   - 類似した予測強度
   - 環境変数のみで説明

2. **Area-interaction model** (ppmlasso):
   - 空間相互作用を考慮
   - 西方の追加的な高強度領域を検出

3. **Log-Gaussian Cox process** (R-INLA):
   - 空間依存性を考慮
   - ppmlassoと類似した追加領域を検出

### 主要な発見

**最低年間気温の効果**（Figure 6a）:
- 二次項が有意に負（ppmlasso, R-INLA）
- Maxentの応答曲線とも整合
- 気候変動シナリオ下での分布縮小を示唆

## 限界・残された課題

### 著者が述べる限界

1. **計算コスト**
   - LGCPのベイズ推論（MCMC）は大規模データで困難
   - INLAは高速だが近似

2. **モデル選択**
   - 複数のモデル（Poisson, Gibbs, Cox）から選ぶ基準が不明確
   - データ駆動的な選択手法が必要

3. **Quadrature pointsの最適化**
   - 数値積分の精度を定量的に評価する手法が限定的
   - 適応的Quadrature手法の開発が必要

4. **サンプリングバイアスの推定**
   - バイアス変数の選択が主観的
   - バイアス関数の不確実性の定量化が困難

### 本研究の観点からの限界

1. **厳密なベイズ推論の欠如**
   - Poisson PPMは最尤推定が主
   - データ拡張によるGibbsサンプリングは扱っていない

2. **マーク付き点過程の未対応**
   - 点の位置のみ（マークなし）
   - 組成データ等の複雑なマークは範囲外

3. **計算効率化手法の限定的議論**
   - NNGPへの言及なし
   - 大規模データへのスケーラビリティは課題

## 本研究（MMCP）との関係

### 借用する要素

1. **Presence-onlyデータの点過程モデリング**
   - 強度関数$\lambda(s)$による定式化
   - Quadrature pointsの概念
   - サンプリングバイアスへの対処方法

2. **Poisson点過程の理論的基盤**
   - 対数尤度関数の構造
   - 数値積分による近似
   - GLMとの関係

3. **空間依存性の扱い**
   - Cox過程（特にLGCP）の枠組み
   - Gaussian processによる空間相関のモデリング

4. **モデル仮定のチェック**
   - 空間独立性の診断
   - 残差分析

### MMCPによる拡張

1. **厳密なベイズ推論**
   - Rennerは最尤推定・INLAが主
   - MMCPはデータ拡張+Gibbsサンプリング（Moreira & Gamerman 2022）
   - Pólya-Gamma augmentationの活用

2. **組成マークの追加**
   - Rennerは点の位置のみ
   - MMCPは各点に組成値マーク（黒曜石産地比率）
   - Multinomial logitとの統合

3. **計算効率化**
   - RennerはフルGP（LGCPでの計算困難を指摘）
   - MMCPはNNGP（Datta et al. 2016）による効率化
   - $O(nm^3)$の計算複雑度

4. **Thinningモデルの統合**
   - Rennerは観測過程を暗黙的に扱う
   - MMCPは観測確率を明示的にモデリング（Moreira & Gamerman 2022）

### 位置づけ

Renner et al. (2015)は：
- Presence-onlyデータに対する点過程モデルの包括的レビュー
- 種分布モデリングの文脈での理論と実践の橋渡し
- 複数のソフトウェアツールの統一的理解
- Quadrature pointsとサンプリングバイアスへの体系的アプローチ

MMCPは：
- この理論的基盤を継承
- ベイズ推論（データ拡張+Pólya-Gamma）による厳密推論
- 組成マークの追加
- NNGP による計算効率化
- Thinningモデルによる観測過程の明示的モデリング

Rennerは**点過程モデリングの基本概念**を提供し、後続研究（Moreira & Gamerman 2022, Moreira et al. 2024）が**厳密なベイズ推論手法**を開発。MMCPはそこに**組成データ**を統合する。

## 引用すべき箇所

### Presence-onlyデータの課題（Abstract）
> "Presence-only data are widely used for species distribution modelling, and point process regression models are a flexible tool that has considerable potential for this problem, when data arise as point events."

### 点過程モデルの利点（Abstract）
> "Advantages include (and are not limited to) clarification of what the response variable is that is modelled; a framework for choosing the number and location of quadrature points (commonly referred to as pseudoabsences or 'background points') objectively; clarity of model assumptions and tools for checking them; models to handle spatial dependence between points when it is present; and ways forward regarding difficult issues such as accounting for sampling bias."

### 強度関数の解釈（Point process models）
> "It should be emphasised that in most instances, the intensity $\lambda(s)$ does not reflect the expected abundance per unit area of a species; rather, it reflects the expected abundance of species reportings. It can typically only be used to make inferences about relative patterns in species abundance (Fithian & Hastie 2013)."

### 空間独立性の仮定（The Poisson case）
> "Assumption (a), in which point locations are independent, is a restrictive assumption which often is not satisfied by presence-only data. Methods to check the independence assumption are described in Section 'Software for fitting point process models' and methods to fit models that account for dependence are described in Section 'Spatial dependence in point processes'."

### Log-Gaussian Cox process（Spatial dependence）
> "An alternative way to deal with clustering and the effects of unmeasured covariates is by fitting a Cox process, the most common example of which is the spatial log-Gaussian Cox process (LGCP) (Møller, Syversveen & Waagepetersen 1998). This can be understood as a point process analogue of a generalised linear mixed model with a random intercept that is normally distributed."

### Quadrature pointsの選択（How to choose quadrature points）
> "A first step in fitting a PPM is selection of quadrature points, to allow estimation of the PPM likelihood. This is an equivalent issue to that of background selection (Section 'Regression models of presence-background data'). However, our choice of the term 'quadrature points' in this Section reflects a desire to pose the question of their choice as a quadrature problem, which clarifies their role in analysis and provides a framework for their selection."

### GLMとMaxentの問題点（Software comparison）
> "A key issue, however, with implementations of Poisson PPM using GLM and maxent software is the lack of assumption checking tools. One should not take 'on faith' the assumption that there is no spatial dependence in the data beyond that explained by environmental variables included in the model."

### Quadrature問題としての定式化（Discussion）
> "Our hope is that approaching the 'pseudo-absence problem' via numerical quadrature will shift attention of analysis away from quadrature point choice and towards where it belongs - developing and interpreting a plausible model for intensity as a function of environment and possibly observer bias variables."
