# Pawlowsky-Glahn & Egozcue (2001) - シンプレックス上の統計分析への幾何学的アプローチ

## 基本情報
- **タイトル**: Geometric approach to statistical analysis on the simplex
- **著者**: V. Pawlowsky-Glahn, J. J. Egozcue
- **所属**: Universitat de Girona (Spain), Universitat Politècnica de Catalunya (Spain)
- **掲載**: Stochastic Environmental Research and Risk Assessment (SERRA)
- **年**: 2001
- **DOI**: 記載なし
- **関連論点**: 論点2（組成データの幾何学、統計的推論）

## 主な貢献
- **Aitchison geometryの数学的基盤**: シンプレックスを$(d-1)$次元Hilbert空間として定式化
- **Metric centerとmetric varianceの一般化**: 任意のHilbert空間における中心と分散の定義
- **Closed geometric meanの正当化**: 組成データの自然な中心測度としての理論的根拠
- **不偏推定量の定義**: Metric unbiasedness（$M_\theta$-不偏性）の概念
- **最良線形不偏推定量（BLUE）の拡張**: シンプレックス上のBLU推定量（$M_a$-best estimator）

## 手法の概要

### 動機と背景

**Aitchison (1982)のlogratio approachへの抵抗**:
- 実務家は新理論を受け入れにくい
- **問題**: 逆変換後の推定量が古典的性質（不偏性、最小分散）を持たない

**本論文の答え**:
- Metric center と metric variance（Aitchison 2001, Pawlowsky-Glahn & Egozcue 2001）
- Hilbert空間の幾何学的構造に基づく定義
- 古典的性質をシンプレックス上で再構築

### 一般的なHilbert空間における定義

#### $m$次元Hilbert空間 $\mathscr{E}$

**演算**:
- 内部演算: $\oplus$ (加法)
- 外部演算: $\otimes$ (スカラー倍)
- 逆元: $\ominus$ (減法), $\oslash$ (除法)
- 中性元: $\mathbf{e}$

**内積と距離**:
- 内積: $\langle \mathbf{x}, \mathbf{y} \rangle$
- ノルム: $\|\mathbf{x}\| = \sqrt{\langle \mathbf{x}, \mathbf{x} \rangle}$
- 距離: $d(\mathbf{x}, \mathbf{y})$

**距離の不変性**:

$$
d(\mathbf{x} \oplus \mathbf{z}, \mathbf{y} \oplus \mathbf{z}) = d(\mathbf{x}, \mathbf{y})
$$

$$
d(\alpha \otimes \mathbf{x}, \alpha \otimes \mathbf{y}) = |\alpha| \cdot d(\mathbf{x}, \mathbf{y})
$$

**重要な性質**: 全ての$m$次元Hilbert空間は$\mathbb{R}^m$と**等長同型** (isometric)

$$
h: \mathscr{E} \to \mathbb{R}^m, \quad d(\mathbf{x}, \mathbf{y}) = d_e(h(\mathbf{x}), h(\mathbf{y}))
$$

### Metric centerとmetric variance

#### 定義1: Metric variance around $\xi$

**Dispersion（分散）**:

$$
\text{Mvar}[\mathbf{X}, \xi] = \mathrm{E}[d^2(\mathbf{X}, \xi)]
$$

幾何学的解釈: ランダムベクトル$\mathbf{X}$から固定点$\xi$までの期待2乗距離

#### 定義2: Metric center

**中心** $\text{Mcen}[\mathbf{X}]$:

$$
\text{Mcen}[\mathbf{X}] = \arg\min_{\xi \in \mathscr{E}} \text{Mvar}[\mathbf{X}, \xi]
$$

#### 定義3: Metric variance (around the center)

**Metric variance**:

$$
\text{Mvar}[\mathbf{X}] = \text{Mvar}[\mathbf{X}, \text{Mcen}[\mathbf{X}]] = \mathrm{E}[d^2(\mathbf{X}, \text{Mcen}[\mathbf{X}])]
$$

**Metric standard deviation**:

$$
\text{Mstd}[\mathbf{X}] = \sqrt{\text{Mvar}[\mathbf{X}]}
$$

### 基本的な命題

**命題1**: 等長同型$h: \mathscr{E} \to \mathbb{R}^m$に対して：

$$
\text{Mcen}[\mathbf{X}] = h^{-1}(\mathrm{E}[h(\mathbf{X})])
$$

**命題2**: $h(\mathbf{X}) = \mathbf{Y} \in \mathbb{R}^m$のとき：

$$
\text{Mvar}[\mathbf{X}] = \sum_{i=1}^m \text{Var}[Y_i]
$$

**命題3**: 線形性（平行移動と拡大縮小）:

$$
\text{Mcen}[(\alpha \otimes \mathbf{X}) \oplus \mathbf{b}] = (\alpha \otimes \text{Mcen}[\mathbf{X}]) \oplus \mathbf{b}
$$

**命題4**: Centering:

$$
\text{Mcen}[\mathbf{X} \ominus \text{Mcen}[\mathbf{X}]] = \mathbf{e}
$$

**命題5**: 加法的線形性:

$$
\text{Mcen}[\mathbf{X} \oplus \mathbf{Y}] = \text{Mcen}[\mathbf{X}] \oplus \text{Mcen}[\mathbf{Y}]
$$

**命題6**: Varianceのスケール不変性:

$$
\text{Mvar}[(\alpha \otimes \mathbf{X}) \oplus \mathbf{b}] = \alpha^2 \text{Mvar}[\mathbf{X}]
$$

**命題7**: 分散の分解:

$$
\text{Mvar}[\mathbf{X}] = \text{Mvar}[\mathbf{X}, \mathbf{b}] - d^2(\text{Mcen}[\mathbf{X}], \mathbf{b})
$$

**命題8**: 独立ランダムベクトルの分散:

$$
\text{Mvar}[\mathbf{X} \oplus \mathbf{Y}] = \text{Mvar}[\mathbf{X}] + \text{Mvar}[\mathbf{Y}] \quad (\mathbf{X} \perp \mathbf{Y})
$$

**命題9**: Chebyshevの不等式:

$$
P[d(\mathbf{X}, \text{Mcen}[\mathbf{X}]) \geq k \text{Mstd}[\mathbf{X}]] \leq \frac{1}{k^2}
$$

## 推定理論

### パラメータ空間$\Theta$における推定

**定義4**: $M_\theta$-不偏推定量（$M_\theta$-centered/unbiased estimator）:

$$
\text{Mcen}_\theta[\hat{\theta}] = \theta
$$

**定義5**: $M_\theta$-bias:

$$
\text{Mcen}_\theta[\hat{\theta} \oplus_\theta \theta^{-1}] = \text{Mcen}_\theta[\hat{\theta}] \ominus_\theta \theta
$$

**定義6**: $M_\theta$-quadratic error:

$$
\text{Mvar}_\theta[\hat{\theta}, \theta] = \mathrm{E}[d_\theta^2(\hat{\theta}, \theta)]
$$

**命題12**: Bias-variance分解:

$$
\text{Mvar}_\theta[\hat{\theta}, \theta] = \text{Mvar}_\theta[\hat{\theta}] + d_\theta^2(\text{Mcen}_\theta[\hat{\theta}], \theta)
$$

**定義7**: $M_\theta$-efficiency:

$$
\hat{\theta}_1 \text{ is more } M_\theta\text{-efficient than } \hat{\theta}_2 \Leftrightarrow \text{Mvar}_\theta[\hat{\theta}_1, \theta] < \text{Mvar}_\theta[\hat{\theta}_2, \theta]
$$

**定義8**: $M_\theta$-best estimator (within class $\hat{\Theta}$):

$M_\theta$-不偏かつ最小$M_\theta$-varianceを持つ推定量

**定義9**: Linear estimator:

$$
\hat{\theta} = \bigoplus_{n=1}^N (\alpha_n \otimes_\theta g(\mathbf{X}_n))
$$

## シンプレックスへの応用

### シンプレックスの定義

**$d$-part composition** $\mathbf{x} = (x_1, \ldots, x_d)'$:

$$
\mathscr{S}_c^d = \{\mathbf{x} = (x_1, \ldots, x_d)' : x_i > 0, \sum_{i=1}^d x_i = c\}
$$

ここで$c$は定数（$c=1$ for parts per unit, $c=100$ for percent）

### Aitchison演算

**Perturbation（摂動）**: $\mathbf{x}, \mathbf{y} \in \mathscr{S}_c^d$に対して

$$
\mathbf{x} \circ \mathbf{y} = \mathscr{C}(x_1 y_1, \ldots, x_d y_d)'
$$

**Power transformation（べき変換）**: $\mathbf{x} \in \mathscr{S}_c^d$, $\alpha \in \mathbb{R}$に対して

$$
\alpha \diamond \mathbf{x} = \mathscr{C}(x_1^\alpha, \ldots, x_d^\alpha)'
$$

**Closure operation（閉包操作）**: $\mathbf{z} = (z_1, \ldots, z_d)'$に対して

$$
\mathscr{C}(\mathbf{z}) = \left(\frac{c \cdot z_1}{\sum_{i=1}^d z_i}, \ldots, \frac{c \cdot z_d}{\sum_{i=1}^d z_i}\right)'
$$

**Neutral element（中性元）**:

$$
\mathbf{e} = \left(\frac{c}{d}, \ldots, \frac{c}{d}\right)' \quad \text{(baricenter)}
$$

### Aitchison geometry

**内積**（Aitchison inner product）:

$$
\langle \mathbf{x}, \mathbf{y} \rangle_a = \frac{1}{d} \sum_{i<j} \ln\frac{x_i}{x_j} \ln\frac{y_i}{y_j}
$$

**ノルム**（Aitchison norm）:

$$
\|\mathbf{x}\|_a = \sqrt{\frac{1}{d} \sum_{i<j} \left(\ln\frac{x_i}{x_j}\right)^2}
$$

**距離**（Aitchison distance）:

$$
d_a(\mathbf{x}, \mathbf{y}) = \sqrt{\frac{1}{d} \sum_{i<j} \left(\ln\frac{x_i}{x_j} - \ln\frac{y_i}{y_j}\right)^2}
$$

### CLR変換（等長同型）

**Centered log-ratio transformation**:

$$
\text{clr}(\mathbf{x}) = \left(\ln\frac{x_1}{g(\mathbf{x})}, \ldots, \ln\frac{x_d}{g(\mathbf{x})}\right)
$$

ここで$g(\mathbf{x}) = (\prod_{i=1}^d x_i)^{1/d}$（幾何平均）

**性質**:
- $\text{clr}: \mathscr{S}_c^d \to \mathbb{R}^d$の$(d-1)$次元超平面（原点を通り、シンプレックスに平行）
- 等長同型: $d_a(\mathbf{x}, \mathbf{y}) = d_e(\text{clr}(\mathbf{x}), \text{clr}(\mathbf{y}))$
- 逆変換: 指数関数 + closure operation

### シンプレックス上の統計量

**命題16**: Metric center of random composition:

$$
\text{Mcen}_a[\mathbf{X}] = \mathscr{C}(\exp\{\mathrm{E}[\ln(X_1)]\}, \ldots, \exp\{\mathrm{E}[\ln(X_d)]\})'
$$

= **Closed geometric mean**（Aitchison 1997の定義）

**命題17**: Metric variance:

$$
\text{Mvar}_a[\mathbf{X}] = \frac{1}{d} \sum_{i<j} \text{Var}\left[\ln\frac{X_i}{X_j}\right] = \sum_{i=1}^d \text{Var}\left[\ln\frac{X_i}{g(\mathbf{X})}\right]
$$

= **Total variance**（Aitchison 1997）

**系**:

$$
\mathrm{E}\left[\ln\frac{X_i}{X_j}\right] = \ln\frac{\gamma_i}{\gamma_j}, \quad \gamma = \text{Mcen}_a[\mathbf{X}]
$$

### シンプレックス上の推定

**$M_a$-best linear unbiased estimator**:

$$
\overline{\mathbf{X}}_a = \bigcirc_{n=1}^N \left(\frac{1}{N} \diamond \mathbf{X}_n\right)
$$

**性質**:

$$
\text{Mvar}_a[\overline{\mathbf{X}}_a] = \frac{\text{Mvar}_a[\mathbf{X}]}{N}
$$

### 標準化組成

**Standardized random composition**:

$$
\mathbf{U} = \frac{1}{\text{Mstd}_a[\mathbf{X}]} \diamond \left(\mathbf{X} \circ (\text{Mcen}_a[\mathbf{X}])^{-1}\right)
$$

**性質**:
- $\text{Mcen}_a[\mathbf{U}] = \mathbf{e}$ (baricenter)
- $\text{Mvar}_a[\mathbf{U}] = 1$ (unit variance)

### Perturbationの影響

**命題3の適用**:
- Perturbation: $\mathbf{X} \circ \mathbf{b}$ → Metric centerも摂動される
- Centering: $\mathbf{X} \circ (\text{Mcen}_a[\mathbf{X}])^{-1}$ （Buccianti et al. 1999等の理論的裏付け）

**命題6の適用**:
- Perturbationは metric varianceに影響しない
- Power transformation $\alpha \diamond \mathbf{X}$ → $\text{Mvar}_a[\alpha \diamond \mathbf{X}] = \alpha^2 \text{Mvar}_a[\mathbf{X}]$

## 実証結果

本論文は理論的枠組みを提供する純粋数学論文のため、実データの実証結果はなし。

ただし、以下の既存研究への理論的基盤を提供：
- **Buccianti et al. (1999)**: Centering による三角図の可視化
- **Martín-Fernández et al. (1999)**: 組成データの発散測度
- **Eynatten et al. (2001)**: 三角図での摂動の理解

## 限界・残された課題

### 理論的限界

1. **Hilbert空間構造の一意性**:
   - シンプレックスに複数のHilbert空間構造を定義可能
   - Aitchison geometryの選択の一意性・最適性は未証明
   - **課題**: 他の距離（例: Fisher-Rao metric）との比較

2. **非連続データへの拡張**:
   - 離散分布、混合分布への適用が未発達
   - ゼロ成分を持つ組成データの扱いが不明確

3. **高次モーメントの定義**:
   - Metric centerとmetric varianceのみ定義
   - **課題**: Skewness, kurtosisの対応物

4. **複素Hilbert空間への拡張**:
   - 論文で言及されているが詳細な展開なし
   - **応用**: 周期的組成データ（時系列）

### 実装上の課題

5. **CLR変換の特異性**:
   - CLRは$(d-1)$次元超平面へ写像 → 特異共分散行列
   - ilr変換（Egozcue et al. 2003）がより実用的
   - **課題**: CLRとilrの使い分け基準

6. **数値的安定性**:
   - 幾何平均$g(\mathbf{x})$の計算時のアンダーフロー
   - 極端に小さい成分の対数変換
   - **対策**: log-sum-exp trick

7. **大規模データへのスケーラビリティ**:
   - $d$が大きい場合（例: 遺伝子発現データ、$d > 10000$）
   - Aitchison距離の計算コスト$O(d^2)$
   - **課題**: 疎な組成データへの適用

### 応用上の課題

8. **統計検定への拡張**:
   - Metric varianceの分布理論が未発達
   - **課題**: 仮説検定、信頼区間の構成

9. **多群比較**:
   - 複数の組成分布の比較（ANOVA的手法）
   - **課題**: $M_a$-distance based MANOVA

10. **回帰モデルへの組み込み**:
    - 応答変数・説明変数が組成の場合
    - **課題**: Compositional regression with metric unbiasedness

11. **ベイズ推論への拡張**:
    - 事前分布のmetric centerとmetric variance
    - **課題**: Dirichlet分布のAitchison geometry上の解釈

## 本研究（MMCP）との関係

### 借用すべき手法

1. **Aitchison geometryの理論的基盤**:
   - **直接的な応用**: MMCPで黒曜石組成データの距離測度
   - Aitchison distance $d_a(\mathbf{x}, \mathbf{y})$を用いた産地間距離
   - 空間的な組成変動のモデリング

2. **Closed geometric mean**:
   - 黒曜石組成の平均値推定
   - 産地ごとの代表的組成の計算
   - **性質**: Metric unbiasedness, minimum metric variance

3. **Total variance（Metric variance）**:
   - 組成データの変動性の測度
   - 産地内変動 vs 産地間変動の定量化
   - **応用**: 産地識別の判別力評価

4. **Perturbation演算**:
   - Centering: $\mathbf{X} \circ (\text{Mcen}_a[\mathbf{X}])^{-1}$
   - 産地ごとの組成を標準化
   - 三角図での可視化（Buccianti et al. 1999）

5. **Power transformation**:
   - 組成データのスケーリング
   - 標準化組成$\mathbf{U}$の構成
   - **応用**: 異なる産地のデータの統合

6. **CLR変換による等長同型**:
   - シンプレックス → $\mathbb{R}^{d-1}$超平面
   - ユークリッド幾何学の手法が適用可能
   - **基盤**: Aitchison (1983)のPCAの理論的根拠

### 拡張すべき点

1. **空間過程への統合**:
   - 本論文: 独立同分布のランダムベクトル
   - MMCP: 空間的に相関する組成データ → 空間相関 + Aitchison geometry
   - **課題**: Metric varianceの空間的構造

2. **点過程への組み込み**:
   - Pawlowsky-Glahn & Egozcue (2001): 組成データのみ
   - MMCP: 点過程（遺跡位置）+ 組成マーク
   - **課題**: Aitchison距離を用いたCox過程の強度関数

3. **Marked point processの統計理論**:
   - Metric centerとmetric varianceをマーク付き点過程に拡張
   - **課題**: 点過程の強度関数 + マークの同時分布

4. **空間的に変化する組成分布**:
   - 産地からの距離で組成分布のmetric centerが変化
   - **課題**: Spatially varying compositional parameters

5. **ベイズ推論への適用**:
   - Metric unbiasedness in Bayesian framework
   - 事前分布・事後分布のmetric center
   - **課題**: Aitchison geometry上のDirichlet分布

6. **階層モデルへの拡張**:
   - 遺跡レベル + 産地レベルの階層
   - Hierarchical metric center
   - **課題**: Random effects in Aitchison geometry

7. **時空間モデルへの適用**:
   - 時代（縄文早期〜後期）での組成変化
   - 時空間的なmetric variance
   - **課題**: Spatio-temporal Aitchison geometry

### MMCPにおける位置づけ

**第3章での役割**:
- 組成データの幾何学的構造の理論的基盤
- Aitchison (1982, 1983)の数学的厳密化
- Closed geometric meanとtotal varianceの正当化

**実装の基礎**:
- `bayesian_statistics/`モジュールでのCLR変換実装の理論
- Aitchison距離の計算アルゴリズム
- Perturbation演算の実装

**方法論的貢献**:
- 不偏推定量の新定義 → 逆変換後の不偏性
- シンプレックスを「制約付き空間」ではなく「固有の幾何学を持つHilbert空間」として扱う
- 古典的統計手法の自然な拡張

**理論的保証**:
- CLRによる等長同型 → ユークリッド空間の理論が適用可能
- Metric unbiasedness → 変換後の推定量の最適性
- Chebyshevの不等式 → 確率的な挙動の保証

## 引用すべき箇所

### Aitchison approachへの抵抗

> "This approach makes it possible to perform classical statistical analysis on transformed data and to back transform the results, which is a clear advantage due to the large amount of methods available for multivariate normally distributed phenomena and the robustness of those. But there has been a certain reluctance in using the new approach by practitioners, which, besides the usual resistance to new theories, is due to the lack of classical properties of backtransformed estimators and models, like unbiasedness and minimum variance." (Introduction)

→ Log ratio approachの利点と実務家の抵抗の理由

### 幾何学的アプローチの動機

> "Nevertheless, the center can be defined as that value $\mu$ which minimizes the expected squared Euclidean distance $\mathrm{E}[d_e(X, \xi)^2]$, and the variance $\sigma^2$ can be defined as the expected value of the squared Euclidean distance around $\mu$, $\sigma^2 = \mathrm{E}[d_e(X, \mu)^2]$. [...] To our understanding, this geometric approach gives its real meaning to the center as a measure of central tendency and to the variance as a measure of dispersion." (Section 1)

→ 中心と分散の幾何学的定義

### Hilbert空間の等長同型

> "Note that an $m$-Hilbert space is always isometric to a real Euclidean space, both having the same dimension $m$. This property is essential for subsequent developments." (Section 2)

→ Hilbert空間 = ユークリッド空間（等長同型）の重要性

### Metric centerの定義

> "The metric center of the distribution of $\mathbf{X}$ is that element $\xi \in \mathscr{E}$ which minimizes $\text{Mvar}[\mathbf{X}, \xi]$." (Definition 2)

→ Metric centerの最適化による定義

### 基本命題

> "If $h: \mathscr{E} \to \mathbb{R}^m$ is an isometry, then $\text{Mcen}[\mathbf{X}] = h^{-1}(\mathrm{E}[h(\mathbf{X})])$." (Proposition 1)

→ 等長同型による期待値の対応

### Aitchison geometry

> "To refer to the properties of $(\mathscr{S}_c^d, \circ, \diamond)$ as a $(d-1)$-Hilbert space, we shall talk globally about the Aitchison geometry on the simplex, and in particular about the Aitchison distance, norm and inner product." (Section "Estimation on the simplex")

→ Aitchison geometryの命名と定義

### Closed geometric meanの正当化

> "The metric center $\text{Mcen}_a[\mathbf{X}]$ of $\mathbf{X}$ is the center or closed geometric mean of $\mathbf{X}$, [...] This result is actually the original definition of center of a random composition given by Aitchison (1997)." (Proposition 16)

→ Closed geometric mean = Metric center（理論的正当化）

### Total varianceとの等価性

> "The first equality states that the metric variance with respect to the Aitchison distance is identical to the total variance defined by Aitchison (1997)." (Proposition 17)

→ Aitchison (1997)の定義との整合性

### Metric varianceの最小化

> "As a result of these statements, we can say that the closed geometric mean of random compositions minimizes the metric variance on the simplex with respect to the Aitchison distance. We can also say that the total variance defined by Aitchison (1997) is an appropriate measure of compositional variability within the simplex, as it coincides with the expected value of the squared Aitchison distance to the metric center of the distribution." (After Proposition 17)

→ Closed geometric meanとtotal varianceの最適性

### Perturbationの影響

> "Perturbation of a random composition affects the metric center in that it leads to a perturbed metric center (Proposition 3), whereas it has no effect on the metric variance (Proposition 6)." (After Proposition 17)

→ Perturbation演算の性質

### Centeringの理論的支持

> "As a consequence, we can center the random composition by perturbing it with the inverse of the metric center (proposition 4), thus giving theoretical support to the approach presented in (Buccianti et al., 1999; Martín-Fernández et al., 1999; Eynatten et al., 2001)." (After Proposition 17)

→ 既存手法への理論的根拠の提供

### $M_a$-best estimator

> "Taking [...] the identity function $g(\mathbf{X}_n) = \text{id}(\mathbf{X}_n) = \mathbf{X}_n$ we obtain that $\overline{\mathbf{X}}_a = \bigcirc_{n=1}^N (\frac{1}{N} \diamond \mathbf{X}_n)$ is the $M_a$-best id-estimator of $\text{Mcen}_a[\mathbf{X}]$ within the class of $M_a$-linear $M_a$-unbiased id-estimators of $\text{Mcen}_a[\mathbf{X}]$." (After Proposition 17)

→ Sample closed geometric mean = BLU推定量

### 結論

> "The existence of an appropriate $m$-Hilbert space structure in the simplex suggests a different approach to the statistical analysis of compositional data based on geometric reasoning. Based on this approach, which is completely parallel to the usual one in Euclidean space, it is straightforward to define reasonable properties for estimators of compositional parameters." (Conclusions)

→ 幾何学的アプローチの全体像

> "In particular, the closed geometric mean is a linear, unbiased estimator that minimizes the metric variance with respect to the Aitchison geometry on the simplex." (Conclusions)

→ Closed geometric meanの最適性のまとめ
