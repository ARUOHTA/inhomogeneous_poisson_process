# Aitchison (1983) - 組成データの主成分分析

## 基本情報
- **タイトル**: Principal component analysis of compositional data
- **著者**: J. Aitchison
- **所属**: University of Hong Kong, Department of Statistics
- **掲載**: Biometrika
- **年**: 1983
- **DOI**: 記載なし
- **関連論点**: 論点2（組成データの統計手法、次元削減）

## 主な貢献
- **Log linear contrast PCAの提案**: 組成データの曲率（curvature）と定和制約の両方に対処する主成分分析
- **Centered log-ratio (clr)変換の導入**: 幾何平均を基準とするlog ratio変換 $\log\{x_i / g(x)\}$
- **等方的共分散構造の特定**: シンプレックス上の球対称性に対応する共分散行列$G_{d+1}$
- **不変性の確立**: 任意の成分を基準に選んでも同じ主成分が得られることを証明
- **従来手法の批判的検討**: Le Maitre (1968), Webb & Briggs (1966)の限界を指摘

## 手法の概要

### 組成データPCAの2つの難しさ

**1. 曲率の問題** (Curvature difficulty):
- 組成データはしばしば顕著な曲線パターンを示す
- 標準的な線形PCAでは不十分（Figure 1(b)のAphyric Skye lavasの例）
- **解決策**: 非線形変換（対数変換）により曲率を捉える

**2. 定和制約の問題** (Constraint difficulty):
- シンプレックス $\mathbb{S}^d = \{(x_1, \ldots, x_d): x_i > 0, x_1 + \ldots + x_d < 1\}$
- 生の比率の共分散・相関は解釈困難（Pearson 1897以来の問題）
- 負の相関へのバイアス、「null correlation」の概念（Chayes & Kruskal 1966）
- **解決策**: Log ratio変換により$\mathbb{R}^d$の標準的多変量手法を適用

### 従来手法の限界

**方法(i)**: $\text{cov}(x_{-j})$ （1成分を除外）
- **問題**: 残りの成分も制約$\sum_{i \neq j} x_i < 1$を受ける
- どの成分$x_j$を除外するかで結果が変わる（不変性の欠如）

**方法(ii)**: $\text{cov}(x^{(d+1)})$ （Le Maitre 1968）
- **問題**: 特異行列（ゼロ固有値1つ）、生の比率の線形結合
- 解釈困難な相関構造
- 曲率を捉えられない（線形主成分のみ）

**方法(iii)**: $\text{cov}(x_{-j} / x_j)$ （Webb & Briggs 1966）
- **問題**: 依然として線形、$x_j$の選択に依存
- 正しい方向だが不十分

### Log linear contrast PCA

#### 基本アイデア

**Aitchison (1981, 1982)の変換**を出発点：

$$
y^{(d)} = \log(x_{-j} / x_j)
$$

共分散行列$\Omega_j = \text{cov}\{\log(x_{-j} / x_j)\}$を考えるが、$x_j$の選択依存性が問題。

**対称化の鍵**: 主成分は**log linear contrast**の形：

$$
\sum_{i \neq j} a_i \log(x_i / x_j) = \sum_{i=1}^{d+1} a_i \log x_i, \quad \text{where } a_1 + \ldots + a_{d+1} = 0
$$

**幾何平均への正規化**:

$$
g(x) = (x_1 \cdots x_{d+1})^{1/(d+1)}
$$

を用いて、

$$
\sum_{i=1}^{d+1} a_i \log x_i = \sum_{i=1}^{d+1} a_i \log\{x_i / g(x)\}
$$

#### 対称的定式化

**Centered log-ratio (clr)変換**:

$$
y_i = \log\{x_i / g(x)\}, \quad i=1, \ldots, d+1
$$

**共分散行列**:

$$
\Omega = \text{cov}[\log\{x^{(d+1)} / g(x)\}]
$$

**性質**:
- $\Omega$は半正定値（positive-semidefinite）
- ゼロ固有値1つ：対応する固有ベクトル$u_{d+1}$ = $(1, \ldots, 1)^T$
- 正の固有値$\lambda_1 \geq \ldots \geq \lambda_d > 0$（降順）

**主成分**:

固有値問題:

$$
(\Omega - \lambda I_{d+1}) a = 0
$$

の$d$個の非ゼロ固有値$\lambda_1, \ldots, \lambda_d$と対応する固有ベクトル$a_1, \ldots, a_d$を求める。

固有ベクトルは$u_{d+1}$と直交 → $a_1 + \ldots + a_{d+1} = 0$ (log linear contrast)

**第$k$主成分**:

$$
PC_k = a_{k1} \log(x_1 / g(x)) + \ldots + a_{k,d+1} \log(x_{d+1} / g(x)) = \sum_{i=1}^{d+1} a_{ki} \log x_i
$$

#### 非対称的定式化（計算上の等価性）

**$\Omega$と$\Omega_j$の関係**:

$$
\Omega_j = B_j \Omega B_j^T
$$

ここで$B_j$は$d \times (d+1)$行列で、$j$列目に$-1$のベクトルを挿入した単位行列：

$$
B_j = \begin{bmatrix} I_{j-1} & -u_{j-1} & 0 \\ 0 & -u_{d-j+1} & I_{d-j+1} \end{bmatrix}
$$

**重要な性質**:

$$
B_j B_j^T = H_d
$$

ここで$H_d$は$d \times d$行列で、対角要素2、非対角要素1の定数行列。

**非対称固有値問題**:

$$
(\Omega_j - \mu H_d) b = 0
$$

**対称版との関係**:
- 固有値は同一: $\mu = \lambda$
- 固有ベクトル: $a = B_j^T b$ または $b = a_{-j}$ （$a$から$a_j$を除いたベクトル）
- **log linear principal componentは同一**（不変性の証明）

固有ベクトル$b_1, \ldots, b_d$は正規直交ではなく、$H_d$-正規直交:

$$
b_i^T H_d b_k = \begin{cases} 1 & (i=k) \\ 0 & (i \neq k) \end{cases}
$$

### 等方的共分散構造

**$\mathbb{R}^d$の場合**: 球対称分布の共分散行列 $\sigma^2 I_d$ は直交変換で不変

**シンプレックスの対応物**: $\Omega$の形で全成分に対称、かつ$u_{d+1}$が固有ベクトルとなるもの：

$$
\sigma^2 G_{d+1}, \quad G_{d+1} = I_{d+1} - \frac{1}{d+1} U_{d+1}
$$

ここで$U_{d+1}$は全要素1の$(d+1) \times (d+1)$行列。

**非対称版**:

$$
\sigma^2 H_d = \sigma^2 B_j B_j^T
$$

**性質**:
- $d$個の固有値がすべて$\sigma^2$ （等しい）
- 任意の正規直交系（$u_{d+1}$に直交）を固有ベクトルとして選べる
- **解釈**: 等方的構造では主成分分析で次元削減の利得なし

### 総変動量

**測度**:

$$
\text{Total variability} = \sum_{i=1}^d \lambda_i = \sum_{i=1}^{d+1} \text{var}[\log\{x_i / g(x)\}]
$$

**第1〜第$c$主成分で説明される割合**:

$$
\frac{\lambda_1 + \ldots + \lambda_c}{\lambda_1 + \ldots + \lambda_d}
$$

## 実装

### 計算手順

1. データ: $x_{ri}$ ($r=1, \ldots, n$; $i=1, \ldots, d+1$)

2. 対数変換:
   $$
   z_{ri} = \log x_{ri}
   $$

3. Centering（幾何平均からの偏差）:
   $$
   z_{r\cdot} = \frac{1}{d+1} \sum_{i=1}^{d+1} z_{ri}, \quad y_{ri} = z_{ri} - z_{r\cdot}
   $$

4. 共分散行列の推定:
   $$
   S_y = \text{sample covariance matrix of } \{(y_{r1}, \ldots, y_{r,d+1}): r=1, \ldots, n\}
   $$

5. 固有値・固有ベクトル計算: $S_y$の固有値問題を解く

6. **注意**: ゼロ固有値1つを持つため、$d$個の非ゼロ固有値のみを使用

## 応用例

### Example 1: Steroid metabolites（ステロイド代謝物）
- **データ**: 37健常成人の尿中ステロイド代謝物、3成分
- **パターン**: 楕円形、ほぼ線形（Figure 1(a)）
- **結果**:
  - 第1主成分で変動の大部分を説明
  - clr変換PCAが適切にフィット
  - 線形パターンへの柔軟な対応を実証

### Example 2: Aphyric Skye lavas（スカイ島玄武岩AFM組成）
- **データ**: 23標本、A=Al₂O₃, F=Fe₂O₃, M=MgO（Figure 1(b)）
- **パターン**: 顕著な曲率、ほぼ1次元の曲線変動
- **従来手法（Le Maitre 1968）の失敗**:
  - 線形主成分軸がデータパターンを捉えられない（Figure 1(b)の直線）
- **clr変換PCAの成功**:
  - 曲線パターンを効果的に捉える
  - 第1主成分で変動のほとんどを説明
  - **解釈**: 地質学的トレンド vs 統計的変動パターンの再考（Butler 1979の議論）

## 実証結果

### 固有値と変動説明率

| 例 | 固有値$\lambda_1$ | 固有値$\lambda_2$ | 第1PC説明率 | 第1-2PC説明率 |
|----|------------------|------------------|------------|--------------|
| Steroid metabolites | （大）| （小） | 高い | ほぼ100% |
| Aphyric Skye lavas | （非常に大）| （非常に小） | 極めて高い | ほぼ100% |

### 視覚的検証

**Figure 1**の比較:
- Le Maitre (1968)の線形軸: 曲率を無視、データの散らばりを誤表現
- clr変換PCA: 曲線変動を適切に捉え、楕円・曲線パターンの両方に対応

## 限界・残された課題

### 理論的限界

1. **高次元での解釈困難**:
   - 主成分がlog linear contrastの形 $\sum a_i \log x_i$ with $\sum a_i = 0$
   - 係数$a_i$の地質学的・生物学的解釈が不明確
   - 正負の係数の意味（どの成分が増加/減少するか）

2. **等方的構造の判定**:
   - $H_d$構造の検定手法が未発達
   - 「主成分分析が無益」と「有益」の境界が不明確

3. **変換の選択**:
   - clr変換の唯一性・最適性の理論的根拠が弱い
   - 他の変換（alr, ilr）との比較が不足

4. **外れ値への感度**:
   - 対数変換は外れ値に敏感
   - 極端に小さい$x_i$（例: 検出限界以下）の扱い

### 実装上の課題

5. **ゼロ成分の問題**:
   - $\log(0)$は未定義 → データの前処理が必要
   - 小さな正の値で置換（例: $0.05$）の恣意性
   - 欠測データとの区別

6. **次元削減の判断基準**:
   - スクリープロット（scree plot）の解釈
   - 何個の主成分を保持すべきか?
   - 固有値の有意性検定（Bartlett検定等）の適用可能性

7. **計算の安定性**:
   - $g(x)$の計算時のアンダーフロー/オーバーフロー
   - 共分散行列の数値的安定性

### 応用上の課題

8. **地質学的解釈**:
   - AFM図の曲率 = 地質学的トレンド? 統計的変動?（Butler 1979の疑問）
   - 主成分の物理的意味（化学プロセス、結晶分化など）

9. **生物統計への適用**:
   - 食事組成、遺伝子発現などへの拡張
   - 時系列組成データ（longitudinal compositional data）への拡張

10. **予測への利用**:
    - 主成分スコアによる分類・判別
    - 新標本の組成予測

## 本研究（MMCP）との関係

### 借用すべき手法

1. **Centered log-ratio (clr)変換**:
   - **直接的な関連**: MMCPで黒曜石組成データの前処理
   - $y_i = \log(x_i / g(x))$ where $g(x) = (x_1 \cdots x_d)^{1/d}$
   - プロジェクトで使用しているilr変換の基礎

2. **次元削減**:
   - 高次元組成データ（多数の元素）の主成分分析
   - 第1-2主成分での可視化（2次元散布図）
   - 遺跡クラスタリングの前処理

3. **等方的構造の概念**:
   - 黒曜石組成が産地間で等方的か?（すべての元素が同程度に変動）
   - 特定の元素比が産地識別に重要か?（非等方的）
   - **応用**: 産地特有の元素の選択

4. **曲率への対応**:
   - 組成空間での非線形パターン
   - 混合産地の影響（曲線的な変動）
   - トブラー距離と組成の非線形関係

5. **不変性の理論**:
   - どの元素を基準（divisor）に選んでも同じ結果
   - alr変換での基準元素選択の恣意性を解消

### 拡張すべき点

1. **空間構造との統合**:
   - 本論文: 独立サンプルの組成PCA
   - MMCP: 空間的に相関する遺跡の組成 → 空間PCAの必要性
   - **課題**: clr変換 + 空間相関構造のモデリング

2. **点過程への組み込み**:
   - Aitchison (1983): 組成の変動構造のみ
   - MMCP: 点過程（遺跡位置）+ 組成の同時モデル
   - **課題**: 主成分スコアを点過程の強度関数に組み込む

3. **ベイズ推論への拡張**:
   - 本論文: 標本主成分分析（頻度論的）
   - MMCP: ベイズ推論による不確実性の定量化
   - **課題**: 主成分の事後分布、予測的主成分

4. **多群比較**:
   - 産地ごとの主成分の違い
   - 共通主成分 vs 群特有の主成分（Flury 1988の共通主成分分析）
   - **課題**: 階層的主成分モデル

5. **欠測値・検出限界以下の値**:
   - 微量元素の検出限界
   - ゼロ値の適切な処理（imputation）
   - **課題**: EM algorithm for compositional data

6. **Dynamic PCA**:
   - 時代（縄文早期〜後期）での組成変化
   - 主成分の時間的変動
   - **課題**: Time-varying covariance structure

7. **Sparse PCA**:
   - 多数の元素から少数の重要な元素を選択
   - 産地識別に寄与する元素の特定
   - **課題**: $\ell_1$正則化 in compositional space

### MMCPにおける位置づけ

**第3章での役割**:
- 組成データの次元削減手法として紹介
- clr変換の理論的基盤（Aitchison 1982とセット）
- 可視化・探索的データ解析の手段

**実装の基礎**:
- `bayesian_statistics/`モジュールでの主成分分析実装
- 組成データの散布図行列（pairs plot）の前処理
- 産地別の主成分スコアの比較

**方法論的貢献**:
- 曲率問題への認識 → 非線形空間過程モデルの動機づけ
- 不変性の重要性 → 座標系に依存しないモデリング
- 等方的構造の概念 → 産地特異性の定量化

**可視化への応用**:
- 黒曜石組成の主成分biplot
- 遺跡の主成分スコアの地図上プロット
- 産地クラスタの判別

## 引用すべき箇所

### 組成データPCAの難しさ

> "Compositional data, consisting of vectors of proportions, have proved difficult to handle statistically because of the awkward constraint that the components of each vector must sum to unity. Moreover such data sets frequently display marked curvature so that linear techniques such as standard principal component analysis are likely to prove inadequate." (Abstract)

→ 組成データPCAの2つの本質的困難

### 曲率の問題

> "The second [data set] has a decidedly nonlinear pattern but with little variation about a conceptual curved line. In purely geometrical terms standard principal component analysis, being a linear reduction technique, may fail to capture successfully the essentially one-dimensional curved variability of the second set." (Section 1.2)

→ 線形PCAの失敗例（Skye lavas）

### 生の比率の相関の問題

> "It is now well known, perhaps more so in geological than in either biometrical or statistical circles, that covariances of raw proportions do not have the simple interpretations which can be placed on their counterparts used in describing variability in $\mathbb{R}^d$." (Section 1.3)

→ 制約による相関構造の歪み（地質学者には周知）

### 等方的共分散構造

> "In principal component analysis in $\mathbb{R}^d$ scalar multiples of the identity matrix $I_d$ play a central role as the only covariance structures which remain invariant under the group of orthogonal transformations. [...] The only form of $\Omega$ which is symmetric in all the components [...] is a scalar multiple of $G_{d+1} = I_{d+1} - \{1/(d+1)\} U_{d+1}$." (Section 2.2)

→ シンプレックスにおける球対称性の対応物

### 不変性の証明

> "Thus the asymmetric eigenvalue problem [...] and the symmetric problem (2.3) produce identical eigenvalues and identical log linear principal components, through the relations $a = B_j^T b$ and $b = a_{-j}$." (Section 2.2)

→ どの成分を基準に選んでも同じ主成分

### 総変動量の測度

> "As in standard principal component analysis the sum of the eigenvalues can be used as a measure of the total variability [...] We note that the measure of total variability for compositional data can be expressed in the form $\sum_{i=1}^{d+1} \text{var}[\log\{x_i / g(x)\}]$." (Section 2.2)

→ 組成データの総変動量の定義

### Le Maitre手法の批判

> "Le Maitre (1968) advocates the use of principal components based on the $d$ nonzero eigenvalues and their corresponding eigenvectors. The principal components are then linear contrasts in the raw proportions and thus also subject to the difficulty of any correlation interpretation caused by the constant-sum constraint." (Section 2.1)

→ 既存手法の根本的問題

### 地質学的解釈の再考

> "Instead of investigating the possibility of regarding nonlinearity in the simplex as a form of natural variability he [Le Maitre] gives the departure from the linear an interpretation of a geological trend, in contrast to more recent views that visual trends in compositional data may be a matter of fantasy rather than fact (Butler, 1979)." (Section 2.1)

→ 曲率を「トレンド」とみなすか「自然な変動」とみなすかの論争

### 変換アプローチの正当化

> "The success of this device [log ratio transformation] in a number of applications, and in particular the successful capture of a curved data set by a predictive region (Aitchison, 1982, Fig. 1), encourages the view that it may also prove effective in principal component analysis of data sets with and without curvature." (Section 1.3)

→ Aitchison (1982)からの自然な拡張
