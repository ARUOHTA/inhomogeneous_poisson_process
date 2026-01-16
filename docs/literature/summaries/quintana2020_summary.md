# Quintana et al. (2020) - 依存ディリクレ過程とその関連モデル

## 基本情報
- **タイトル**: The Dependent Dirichlet Process and Related Models
- **著者**: Fernand A. Quintana, Peter Müller, Alejandro Jara, Steven N. MacEachern
- **所属**: Pontificia Universidad Católica de Chile, University of Texas at Austin, Ohio State University
- **年**: 2020
- **ジャーナル/掲載**: レビュー論文（掲載誌の詳細は本文に記載なし）
- **DOI**: 記載なし
- **関連論点**: 論点2（ベイズ推論手法）、論点3（共変量依存の分布モデリング）

## 主な貢献
- **依存ディリクレ過程（DDP）の包括的レビュー**: 過去25年の主要な構成法を体系的に整理
- **完全ノンパラメトリック回帰の定式化**: 応答変数の分布全体が共変量$\boldsymbol{x}$で柔軟に変化するモデル
- **多様な構成法の比較**: MacEachern DDP, KSBP, PSB, HDP, NDP等の特性と利点・欠点を議論
- **理論的性質の解説**: 完全サポート、軌跡の滑らかさ、パーティション構造の共変量依存性
- **応用分野の提示**: 自己回帰モデル、金融時系列、メタ分析等への適用例

## 手法の概要

### 完全ノンパラメトリック回帰問題
応答変数$\boldsymbol{y}_i | F_{\boldsymbol{x}_i} \stackrel{ind}{\sim} F_{\boldsymbol{x}_i}$, $i=1, \ldots, n$を仮定し、共変量依存の確率測度の族を推定：

$$
\mathscr{F} = \{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\}
$$

ここで$F_{\boldsymbol{x}}$は応答空間$\mathscr{Y}$上の確率測度であり、分布の形状全体が$\boldsymbol{x}$で変化可能。

**従来の回帰との違い**:
- 一般化線形モデル (GLM): 平均のみが共変量で変化
- 分位点回帰: 特定の分位点のみが変化
- **DDP**: 分布全体（平均、分散、形状、多峰性など）が変化

### ディリクレ過程混合モデル（DPM）の拡張
標準的なDPM:

$$
f(\boldsymbol{y}|G) = \int_{\Theta} \psi(\boldsymbol{y}, \boldsymbol{\theta}) G(d\boldsymbol{\theta}), \quad G \sim \text{DP}(M, G_0)
$$

スティックブレーキング表現（Sethuraman 1994）:

$$
G(B) = \sum_{h=1}^{\infty} w_h \delta_{\boldsymbol{\theta}_h}(B)
$$

ここで$w_h = V_h \prod_{\ell < h}(1-V_\ell)$, $V_h \stackrel{iid}{\sim} \text{Be}(1, M)$, $\boldsymbol{\theta}_h \stackrel{iid}{\sim} G_0$

DDPは$G$を$G_{\boldsymbol{x}}$に置き換え、共変量依存性を導入。

## 主要な構成法

### 1. MacEachern DDP（基本定義）

**一般形**（MacEachern 1999, 2000）:

$$
G_{\boldsymbol{x}}(\bullet) = \sum_{h=1}^{\infty} w_h(\boldsymbol{x}) \delta_{\boldsymbol{\theta}_h(\boldsymbol{x})}(\bullet)
$$

ここで：
- $V_h(\boldsymbol{x})$: $[0,1]$値の独立確率過程、周辺分布$\text{Be}(1, M_{\boldsymbol{x}})$
- $\boldsymbol{\theta}_h(\boldsymbol{x})$: 独立確率過程、周辺分布$G_{\boldsymbol{x}}^0$
- $w_h(\boldsymbol{x}) = V_h(\boldsymbol{x}) \prod_{\ell < h}[1-V_\ell(\boldsymbol{x})]$

**正準構成**: 2つの独立な確率過程$Z_h^V(\boldsymbol{x})$, $Z_h^{\boldsymbol{\theta}}(\boldsymbol{x})$を変換：

$$
V_h(\boldsymbol{x}) = B^{-1}(\Phi(\sigma^{-1} Z_h^V(\boldsymbol{x})))
$$

ここで$\Phi$: 標準正規CDF、$B$: Beta$(1, M)$ CDF

**代替定義**（Barrientos, Jara & Quintana 2012）: コピュラ$\mathcal{C}_{\mathscr{X}}^V$, $\mathcal{C}_{\mathscr{X}}^{\theta}$を用いて有限次元分布を定義し、Kolmogorovの整合性条件を満たす。

#### (a) Single-weights DDP
重みが共通、原子が共変量依存：

$$
G_{\boldsymbol{x}}(\bullet) = \sum_{h=1}^{\infty} w_h \delta_{\boldsymbol{\theta}_h(\boldsymbol{x})}(\bullet)
$$

**特徴**:
- 標準的なDP推論アルゴリズムが小修正で適用可能
- パーティション$\rho_{\boldsymbol{x}}$の事前分布が$\boldsymbol{x}$で不変
- 計算が比較的容易

#### (b) Single-atoms DDP
原子が共通、重みが共変量依存：

$$
G_{\boldsymbol{x}}(\bullet) = \sum_{h=1}^{\infty} w_h(\boldsymbol{x}) \delta_{\boldsymbol{\theta}_h}(\bullet)
$$

**特徴**:
- パーティション$\rho_{\boldsymbol{x}}$の事前分布が$\boldsymbol{x}$で変化（柔軟性）
- 重み過程$w_h(\boldsymbol{x})$の計算が複雑
- 実装が困難

#### (c) Linear DDP (LDDP)（De Iorio et al. 2004）
原子が線形関数：

$$
\boldsymbol{\theta}_h(\boldsymbol{x}) = \boldsymbol{\beta}_{0h} + \mathbf{B}_h \boldsymbol{x}
$$

ここで$(\boldsymbol{\beta}_{0h}, \mathbf{B}_h) \stackrel{iid}{\sim} G_0$

**応用**: 条件付き密度推定、生存時間解析、生物統計

#### (d) ANCOVA-DDP（De Iorio et al. 2004）
カテゴリカル変数$c \in \{1, \ldots, C\}$と連続変数$\boldsymbol{z}$の混合：

$$
\boldsymbol{\theta}_{ch}(\boldsymbol{z}) = \boldsymbol{\beta}_{0ch} + \mathbf{B}_{ch}\boldsymbol{z}
$$

**特徴**: カテゴリー間で異なる線形関数、群間の情報借用

### 2. Kernel Stick-Breaking Process (KSBP)（Dunson & Park 2008）

$$
G_{\boldsymbol{x}}(\bullet) = \sum_{h=1}^{\infty} \left\{ W(\boldsymbol{x}; \boldsymbol{\Gamma}_h, V_h) \prod_{\ell < h} [1 - W(\boldsymbol{x}; \boldsymbol{\Gamma}_\ell, V_\ell)] \right\} G_h(\bullet)
$$

ここで：
- $W(\boldsymbol{x}; \boldsymbol{\Gamma}_h, V_h) = V_h K(\boldsymbol{x}, \boldsymbol{\Gamma}_h)$
- $K: \mathscr{X} \times \mathscr{X} \to [0,1]$: カーネル関数（例: ガウシアンカーネル）
- $V_h \stackrel{ind}{\sim} \text{Be}(a_h, b_h)$
- $\boldsymbol{\Gamma}_h \stackrel{iid}{\sim} H$: ランダムカーネル位置
- $G_h \stackrel{iid}{\sim} \mathcal{G}$: ランダム確率測度

**簡略版**（単一原子）:

$$
G_{\boldsymbol{x}}(\bullet) = \sum_{h=1}^{\infty} \left\{ W(\boldsymbol{x}; \boldsymbol{\Gamma}_h, V_h) \prod_{\ell < h} [1 - W(\boldsymbol{x}; \boldsymbol{\Gamma}_\ell, V_\ell)] \right\} \delta_{\boldsymbol{\theta}_h}(\bullet)
$$

**特徴**:
- ランダム位置$\boldsymbol{\Gamma}_h$からの距離に基づく重み
- 柔軟な共変量依存構造
- モデル複雑性の低減（簡略版）

**カーネル例**（ガウシアンカーネル）:

$$
K(\boldsymbol{x}, \boldsymbol{\Gamma}_h) = \exp\left\{-\frac{\|\boldsymbol{x} - \boldsymbol{\Gamma}_h\|^2}{2\tau^2}\right\}
$$

### 3. Probit Stick-Breaking Process (PSB)（Chung & Dunson 2009）

$$
G(\bullet) = \sum_{h=1}^{\infty} \left\{ \Phi(\eta_h) \prod_{\ell < h} [1 - \Phi(\eta_\ell)] \right\} \delta_{\boldsymbol{\theta}_h}(\bullet)
$$

ここで$\eta_h \stackrel{iid}{\sim} N(\mu, 1)$, $\Phi$: 標準正規CDF

**共変量依存版**:

$$
G_{\boldsymbol{x}}(\bullet) = \sum_{h=1}^{\infty} \left\{ \Phi(\eta_h(\boldsymbol{x})) \prod_{\ell < h} [1 - \Phi(\eta_\ell(\boldsymbol{x}))] \right\} \delta_{\boldsymbol{\theta}_h}(\bullet)
$$

ここで$\eta_h(\boldsymbol{x}) = \alpha_h + f_h(\boldsymbol{x})$, $\alpha_h \sim N(\mu, 1)$

**特徴**:
- Beta変数をnormalで置き換え
- 変数選択への応用（Chung & Dunson 2009）
- **Logit Stick-Breaking**（Ren et al. 2011）: probitをlogitで置き換えた類似構成

### 4. Hierarchical Mixture of DP（Müller, Quintana & Rosner 2004）

離散的共変量$\mathscr{X} = \{1, \ldots, J\}$の場合：

$$
G_j(\bullet) = \epsilon H_0(\bullet) + (1-\epsilon) H_j(\bullet)
$$

ここで：
- $\epsilon \in [0,1]$: 依存度のパラメータ
- $H_0, H_1, \ldots, H_J$: 独立なDP

**極端ケース**:
- $\epsilon = 1$: 完全プーリング（全群で共通）
- $\epsilon = 0$: 独立（情報借用なし）

**連続共変量への拡張**:

$$
G_{j,\boldsymbol{z}}(\bullet) = \epsilon H_{0,\boldsymbol{z}}(\bullet) + (1-\epsilon) H_{j,\boldsymbol{z}}(\bullet)
$$

ここで$H_{0,\boldsymbol{z}}, H_{j,\boldsymbol{z}}$は独立なMacEachern DDP

**応用**: メタ分析、複数研究間の情報借用

### 5. Hierarchical DP (HDP)（Teh et al. 2006）

$$
G_j | M_j, G \stackrel{ind}{\sim} \text{DP}(M_j, G), \quad j=1, \ldots, J, \quad G | M, G_0 \sim \text{DP}(M, G_0)
$$

**特徴**:
- 階層的構造（条件付きDP）
- 群間で原子（クラスタ）を共有、重みは異なる
- **Chinese Restaurant Franchise**: 多階層クラスタリング
- **応用**: テキスト解析（複数文書・コーパス間でクラスタ共有）

### 6. Nested DP (NDP)（Rodríguez, Dunson & Gelfand 2008）

$$
G_j \stackrel{ind}{\sim} \sum_{h=1}^{\infty} \pi_h \delta_{G_h^*}(\bullet), \quad j=1, \ldots, J, \quad G_h^* | M_2, H \stackrel{iid}{\sim} \text{DP}(M_2, H)
$$

ここで$\pi_h = V_h \prod_{\ell < h}(1-V_\ell)$, $V_h \stackrel{iid}{\sim} \text{Be}(1, M_1)$

**特徴**:
- DPの軌跡の無限混合
- 多層クラスタリング（個体内 + 分布間）
- HDPとの違い: NDPは原子と重みの両方を共有または完全に独立

**制約**: $G_{j_1} = G_{j_2}$または完全独立の2択のみ（中間構成なし）

**代替モデル**:
- **Latent nested process**（Camerlenghi et al. 2019）
- **Semi-hierarchical DP**（Beraha, Guglielmi & Quintana 2020）

### 7. Product of Independent DPs（Gelfand & Kottas 2001）

$$
G_j(\bullet) \equiv H_j(\bullet) \prod_{\ell < j} H_\ell(\bullet), \quad j=1, \ldots, J
$$

ここで$H_j | M_j, H_{0j} \stackrel{ind}{\sim} \text{DP}(M_j, H_{0j})$

**特徴**: 独立DPの積による依存構造の導入

### 8. Normalized Priors（Barrientos et al. 2017）
一般的な完全ランダム測度（CRM）$\tilde{\mu}$の正規化：

$$
P_{\boldsymbol{x}}(A) = \frac{\tilde{\mu}_{\boldsymbol{x}}(A)}{\tilde{\mu}_{\boldsymbol{x}}(\mathscr{Y})}, \quad A \subset \mathscr{Y}
$$

**例**:
- **Normalized generalized gamma process**
- **Normalized inverse Gaussian process**
- **Dependent normalized inverse Gaussian process** (DNIGP)

## 応用例

### 1. 自己回帰モデル (Section 6)
**Old Faithful Geyser data**: 噴火時間の時系列

$$
y_t | y_{t-1} \sim G_{y_{t-1}}
$$

**モデル**: Conditional DDP（条件付き密度を直接モデル化）

**結果**: Figure 4に事後推定密度$G_x$の共変量$x=y_{t-1}$依存性を図示

### 2. 金融時系列
**S&P500 data**: 株価リターンの時系列

**モデル**: LDDP vs 条件付き密度アプローチの比較

**結果**: Figure 5にLDDPと条件付き密度アプローチの違いを図示

## 理論的性質

### Full Support（Barrientos, Jara & Quintana 2012）
**定理**: 確率過程$V_h(\boldsymbol{x})$, $\boldsymbol{\theta}_h(\boldsymbol{x})$が完全サポートを持つ条件下で、DDPは$\mathscr{F} = \{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\}$上で完全サポートを持つ。

**意味**: DDPは十分に柔軟で、任意の共変量依存分布族を近似可能

### 滑らかさ (Smoothness)
**性質**: $\boldsymbol{x} \to \boldsymbol{x}_0$のとき、$F_{\boldsymbol{x}} \to F_{\boldsymbol{x}_0}$（ある位相で）

**実装**: $\text{Corr}\{F_{\boldsymbol{x}}(A), F_{\boldsymbol{x}_0}(A)\} \to 1$ as $\boldsymbol{x} \to \boldsymbol{x}_0$

### パーティション構造
$\rho_{\boldsymbol{x}}$: $G_{\boldsymbol{x}}$からのサンプルサイズ$n$が誘導するパーティション

- **Single-weights DDP**: $p(\rho_{\boldsymbol{x}}|\mathcal{G})$が$\boldsymbol{x}$で不変
- **Single-atoms DDP**: $p(\rho_{\boldsymbol{x}}|\mathcal{G})$が$\boldsymbol{x}$で変化（より柔軟）

## 限界・残された課題

### 実装上の課題

1. **Single-atoms vs Single-weights のトレードオフ**:
   - Single-weights: 計算容易だが、パーティション構造が$\boldsymbol{x}$で不変
   - Single-atoms: パーティション構造が柔軟だが、重み過程の計算が複雑
   - 両方に依存: 最も柔軟だが、計算的に最も困難

2. **計算複雑性**: 多くの構成法で標準的なMCMCアルゴリズムの適用が非自明
   - LDDP, single-weights: 比較的容易
   - Single-atoms, KSBP: 複雑

### 理論的限界

3. **事前情報の導入困難**:
   - 任意の関数（平均、分位点など）の事前分布を直接指定できない
   - 関数の誘導分布の導出が困難
   - **対策**: Kessler, Hoff & Dunson (2015)のMarginal specificationアプローチ（単一分布）
   - **課題**: 共変量依存分布族への拡張が未解決

4. **局所交換可能性** (Local Exchangeability)（Campbell et al. 2019）:
   - 通常の交換可能性の緩和版を導入
   - DDPがde Finettiの定理の対応する測度となる条件
   - **課題**: 拡張と応用が今後の課題

5. **縦断的・空間的・関数的データへの拡張**:
   - 現状のDDPは周辺分布$G_{\boldsymbol{x}}$に焦点
   - 観測間の依存構造（時系列相関、空間相関）を直接モデル化できない
   - **例**: Xu, MacEachern & Xu (2015)の金融データ解析
   - **課題**: 周辺分布と依存構造の分離

6. **正規化による依存構造の導入**（Section 3.8）:
   - より一般的なCRMへの拡張
   - DDP以外の枠組みへの適用

### モデル選択の課題

7. **どのDDPが最良か？**: 一般的な答えはない
   - 応用分野、データの性質、計算資源による
   - 両方に依存するモデルは少ない（計算複雑性のため）
   - 条件付きアプローチ（Section 4）は例外

8. **パーティション構造の評価**:
   - Single-weightsとsingle-atomsで異なる事前分布
   - 実データでどちらが適切かの判断基準が不明確

## 本研究（MMCP）との関係

### 借用すべき手法

1. **Kernel Stick-Breaking Process (KSBP)**:
   - **直接的な関連**: プロジェクトには`20_model_KSBP.ipynb`が存在
   - MMCPでの適用: 黒曜石組成分布が空間的に変化する構造をモデル化
   - カーネル$K(\boldsymbol{x}, \boldsymbol{\Gamma}_h)$で空間的な滑らかさを導入
   - ランダム位置$\boldsymbol{\Gamma}_h$が黒曜石産地に対応する可能性

2. **Linear DDP (LDDP)**:
   - 組成係数$\boldsymbol{\theta}_h(\boldsymbol{x}) = \boldsymbol{\beta}_{0h} + \mathbf{B}_h \boldsymbol{x}$
   - $\boldsymbol{x}$: 距離、標高などの空間共変量
   - 黒曜石組成（多変量連続）の空間的変化をモデル化

3. **完全ノンパラメトリック回帰の枠組み**:
   - 組成分布全体（平均だけでなく）が空間的に変化
   - 黒曜石産地ごとに異なる組成分布の形状を許容
   - 多峰性、非対称性などの複雑な分布形状への対応

4. **情報借用 (Borrowing Strength)**:
   - 近隣遺跡間での組成分布の類似性を活用
   - Hierarchical mixture of DPのアイデア: 共通部分と個別部分の混合
   - 産地ごとの共通分布と遺跡固有の分布

5. **理論的基盤**:
   - Full supportの性質: DDPが十分に柔軟
   - 滑らかさの保証: 空間的に近い位置で分布が類似
   - パーティション構造: 遺跡のクラスタリングが空間的に変化

### 拡張すべき点

1. **点過程への統合**:
   - 本論文: 応答変数$\boldsymbol{y}$の分布が共変量依存
   - MMCP: 点過程（遺跡の位置）+ マーク（組成）の同時モデル
   - **課題**: DDPをCox過程の強度関数に組み込む方法

2. **マーク付き点過程への適用**:
   - DDPは単変量/多変量応答に焦点
   - MMCPはマーク付き点過程 → 位置とマークの同時分布をモデル化
   - **課題**: 点過程の空間構造とDDPの統合

3. **組成データ制約の導入**:
   - DDPは一般的な$\mathbb{R}^d$上の分布
   - 黒曜石組成: シンプレックス$\{\boldsymbol{y}: y_i \geq 0, \sum_i y_i = 1\}$上
   - **課題**: ilr変換後の$\mathbb{R}^{d-1}$上でDDPを適用、または直接シンプレックス上で定義

4. **空間的な距離の導入**:
   - KSBPのカーネル$K(\boldsymbol{x}, \boldsymbol{\Gamma}_h)$にトブラー距離を使用
   - 黒曜石プロジェクトの地形考慮距離行列との統合
   - **課題**: 地形ベースの距離をカーネル関数に組み込む方法

5. **計算効率の向上**:
   - DDPは一般に計算コストが高い
   - 大規模遺跡データ（数千点）への適用
   - **課題**: KSBPの打ち切り、近似アルゴリズムの開発

6. **Preferential Samplingとの統合**:
   - DDPは応答の分布のみをモデル化
   - MMCPは遺跡位置自体が選好的（産地に近い場所に集中）
   - **課題**: DDPとCox過程の強度関数を共有成分で結合

7. **事前分布の導出**:
   - 黒曜石組成の平均・分散の考古学的先行知識
   - DDPの関数の誘導分布を導出
   - **課題**: Kessler et al. (2015)の marginal specification をDDP族に拡張

### MMCPにおける位置づけ

**第3章での役割**:
- 共変量依存のノンパラメトリック分布モデリングの主要アプローチとしてDDPを紹介
- KSBPを中心に解説（プロジェクトで実装済み）
- 他のアプローチ（空間Cox過程、Gelfand-Schliep SVC）との比較材料

**実装の選択肢**:
- Stan/NUTSでのKSBP実装（`20_model_KSBP.ipynb`）
- PyroやPyMC4でのDDP実装の可能性
- LDDPのシンプルな実装と計算効率

**トレードオフの議論**:
- パラメトリック（SVC） vs ノンパラメトリック（DDP/KSBP）
- 計算コスト vs モデル柔軟性
- 解釈可能性 vs 予測性能

**理論的貢献**:
- 完全サポートの性質により、MMCPの柔軟性を理論的に保証
- 滑らかさの概念が空間的な連続性の仮定と整合

## 引用すべき箇所

### 完全ノンパラメトリック回帰の動機

> "Standard regression approaches assume that some finite number of the response distribution characteristics, such as location and scale, change as a (parametric or nonparametric) function of predictors. However, it is not always appropriate to assume a location/scale representation, where the error distribution has unchanging shape over the predictor space." (p.1, Abstract)

→ 従来の回帰の限界とDDPの必要性

> "In fact, it often happens in applied research that the distribution of responses under study changes with predictors in ways that cannot be reasonably represented by a finite dimensional functional form." (p.1, Abstract)

→ 分布の形状全体が共変量で変化するケースの実在性

### DDPの定義

> "The key idea behind the DDP construction is to define a set of random measures that are marginally (i.e. for every possible predictor value $\boldsymbol{x} \in \mathscr{X}$) DP-distributed random measures." (p.4, Section 2.1)

→ DDPの基本アイデア（周辺的にDPを保つ）

> "Intuitively, the constructed DDP can be thought of as taking an ordinary DP and modifying some of its components (i.e. weights and atoms) according to the type of desired indexing or functional dependence of predictors $\boldsymbol{x} \in \mathscr{X}$." (p.4, Section 2.1)

→ DDPの直感的理解

### KSBPの定義

> "The KSBP thus begins with an infinite sequence of basis random distributions $\{G_h\}$ and then constructs covariate-dependent random measures by mixing according to distance from the random locations $\Gamma_h$, with stick-breaking probabilities that are defined as a kernel multiplied by Beta-distributed weights." (p.11, Section 3.2)

→ KSBPの構成原理

### 理論的性質

> "The results in Barrientos, Jara and Quintana (2012) show that under full support of the stochastic processes that are used to convey covariate dependence, the resulting DDP has full support in the space $\mathscr{F} = \{F_{\boldsymbol{x}}: \boldsymbol{x} \in \mathscr{X}\}$." (p.24, Section 7)

→ 完全サポートの性質（柔軟性の理論的保証）

### Single-atoms vs Single-weights

> "The single-weights models are typically easier to fit, as the standard algorithms designed to implement posterior simulation in the context of DPs can be applied with minor adjustments. [...] The single-atoms models are typically less attractive from a computational viewpoint, mainly due to how covariate dependence is encoded in the definition of the weight processes $\{w_h(\boldsymbol{x}): \boldsymbol{x} \in \mathscr{X}\}$." (p.24, Section 7)

→ 計算複雑性のトレードオフ

> "However, the single-atoms DDP allows for the prior probability distribution on the partitions to change with $\boldsymbol{x}$, a feature that is not supported by the single-weights DDP." (p.24, Section 7)

→ パーティション構造の柔軟性の違い

### 事前情報の導入困難性

> "Models for dependent probability distributions do not easily allow for the incorporation of existing prior information about arbitrary functionals. A modeler is unlikely to have prior knowledge about all aspects of a collection of probability measures, but could have real historical prior information about specific functionals (such as the mean or quantile functions)." (p.24-25, Section 7)

→ 実用上の重要な限界

### 情報借用の重要性

> "A main aim here is that subjects under study $j_1$ should inform inference about subjects enrolled in a different but related study $j_2 \neq j_1$." (p.2, Section 1)

→ 関連研究/遺跡間での情報借用の動機

> "Two extreme modeling choices would be (i) to pool all patients and assume one common effects distribution, or (ii) to assume $J$ distinct distributions with independent priors. [...] In most applications, the desired level of borrowing strength is somewhere in-between these two extremes." (p.2, Section 1)

→ 完全プーリングと独立の中間を取るモデルの必要性

### HDPとNDPの違い

> "In the HDP of Teh et al. (2006), the random measures in $\mathcal{G} = \{G_1, \ldots, G_J\}$ share the same atoms but assign them different weights, while in the NDP two distributions $G_{j_1}$ and $G_{j_2}$ either share both atoms and weights (i.e. they are identical) or share nothing at all." (p.15, Section 3.6)

→ 階層的DDPモデル間の重要な違い

### 滑らかさの概念

> "For modeling, one important property is the notion of distributions changing smoothly with respect to $\boldsymbol{x} \in \mathscr{X}$, just as is the case of generalized linear models in the scale of the transformed mean. The smoothness could be expressed as continuity of $F_{\boldsymbol{x}}$ [...] or as the notion that $F_{\boldsymbol{x}}$ 'approaches' $F_{\boldsymbol{x}_0}$ as $\boldsymbol{x} \to \boldsymbol{x}_0$, for instance, $\text{Corr}\{F_{\boldsymbol{x}}(A), F_{\boldsymbol{x}_0}(A)\} \to 1$ as $\boldsymbol{x} \to \boldsymbol{x}_0$ for any event $A$." (p.3, Section 1)

→ 空間的な連続性の数学的定式化

### 今後の課題

> "When data are longitudinal, spatial or functional, the observations may be considered to have dependence that cannot be captured by the marginal distributions $G_{\boldsymbol{x}}$. [...] Many open questions remain in this direction." (p.26, Section 7)

→ 空間/時系列データへの拡張の必要性（MMCPに直接関連）

> "The idea of introducing dependence through normalization [...] can be further exploited and extended to more general cases, including going beyond the context of DDPs." (p.26, Section 7)

→ 正規化によるアプローチの今後の可能性
