# Gneiting & Raftery (2007) - Strictly Proper Scoring Rules（厳密適正スコアリング規則）

## 基本情報
- **タイトル**: Strictly Proper Scoring Rules, Prediction, and Estimation
- **著者**: Tilmann Gneiting, Adrian E. Raftery
- **ジャーナル**: Journal of the American Statistical Association, 102(477), 359-378
- **年**: 2007
- **DOI**: 10.1198/016214506000001437
- **関連論点**: 論点4（モデル診断・評価）

## 主な貢献

確率的予測の評価と推定のための**厳密適正スコアリング規則（strictly proper scoring rules）**の包括的理論を展開。

1. **一般理論**: 一般確率空間上でのproper scoring rulesの特徴付け定理
2. **凸関数との関係**: Bregman divergence、情報測度、エントロピー関数との関連
3. **連続変数のスコア**: Continuous Ranked Probability Score (CRPS)、energy score等の新規提案
4. **カーネル表現**: 負定値関数に基づくkernel scoreの一般的構成
5. **最適スコア推定**: 厳密適正スコアを損失関数とする点推定・区間推定手法
6. **実証評価**: 天気予報データでの不適正スコアの問題点を実証

## 手法の概要

### Proper Scoring Ruleの定義

**確率的予測**: $P \in \mathcal{P}$（確率測度のクラス）

**スコアリング規則**: $S: \mathcal{P} \times \Omega \rightarrow \bar{\mathbb{R}}$
- 予測 $P$ を発表、事象 $\omega$ が実現 → 報酬 $S(P, \omega)$

**期待スコア**:
$$
S(P, Q) = \int S(P, \omega) dQ(\omega)
$$

**Proper scoring rule**: 全ての $P, Q \in \mathcal{P}$ に対し
$$
S(Q, Q) \geq S(P, Q)
$$

**Strictly proper**: 等号成立 ⇔ $P = Q$

**意味**: 予測者の真の信念が $Q$ なら、$Q$ を正直に報告することが期待スコアを最大化

### 特徴付け定理

**定理1**: 正則スコアリング規則 $S$ がproperである ⇔ ある凸関数 $G: \mathcal{P} \to \mathbb{R}$ が存在して
$$
S(P, \omega) = G(P) - \int G^*(P, \omega) dP(\omega) + G^*(P, \omega)
$$

ここで $G^*(P, \cdot)$ は $G$ の点 $P$ における劣接線（subtangent）

**厳密適正**: $G$ が厳密凸 ⇔ $S$ が厳密適正

**情報測度（generalized entropy）**:
$$
G(P) = \sup_{Q \in \mathcal{P}} S(Q, P)
$$

最大達成可能効用

**Divergence function**:
$$
d(P, Q) = S(Q, Q) - S(P, Q) \geq 0
$$

（$P = Q$ のとき等号、厳密適正なら $P \neq Q$ で厳密正）

### 代表的なScoring Rules

#### 1. Logarithmic Score（対数スコア）

$$
S_{\log}(P, \omega) = \log p(\omega)
$$

- **適用**: 予測密度 $p$
- **情報測度**: Shannon entropy $H(P) = -\int p(\omega) \log p(\omega) d\omega$
- **Divergence**: Kullback-Leibler divergence $\text{KL}(Q \| P)$

#### 2. Brier Score（Brierスコア）

二値イベント $A$ に対し、$P(A) = q$ の予測:
$$
S_{\text{Brier}}(q, \omega) = -(\mathbb{1}_A(\omega) - q)^2
$$

（ $-$ 符号は損失 → 報酬への変換）

#### 3. Spherical Score（球面スコア）

$$
S_{\text{sph}}(P, \omega) = \frac{p(\omega)}{\|p\|_{L^2(Q)}} = \frac{p(\omega)}{\sqrt{\int p^2 d\mu}}
$$

- 予測密度の $L^2$ ノルムで正規化
- 計算が容易、有界

#### 4. Quadratic Score（二次スコア）

$$
S_{\text{quad}}(P, \omega) = 2p(\omega) - \int p^2 d\mu
$$

#### 5. Continuous Ranked Probability Score (CRPS)

実数値変数 $X \in \mathbb{R}$ の予測CDF $F$ に対し:
$$
\text{CRPS}(F, x) = -\int_{-\infty}^\infty (F(y) - \mathbb{1}\{x \leq y\})^2 dy
$$

**別表現**:
$$
\text{CRPS}(F, x) = -E_F|X - x| + \frac{1}{2} E_F|X - X'|
$$

ここで $X, X'$ は $F$ からのiid標本

**特性**:
- 絶対誤差（MAE）の一般化
- 点予測 $\delta_x$ なら $\text{CRPS}(\delta_x, y) = -|x - y|$
- 予測の sharpness（分散）と calibration（バイアス）を同時評価

#### 6. Energy Score

多変量 $\mathbf{X} \in \mathbb{R}^d$ の予測分布 $F$ に対し:
$$
\text{ES}(F, \mathbf{x}) = -E_F\|\mathbf{X} - \mathbf{x}\| + \frac{1}{2} E_F\|\mathbf{X} - \mathbf{X}'\|
$$

ここで $\|\cdot\|$ はユークリッドノルム

**特性**:
- CRPSの多変量版
- 任意の次元 $d$ で厳密適正
- 分布の全構造（共分散行列等）を評価

### Kernel Score（カーネルスコア）

**一般構成**: 負定値関数 $k: \Omega \times \Omega \to \mathbb{R}$ を用いて
$$
S_k(P, \omega) = -\int k(\omega, y) dP(y) + \frac{1}{2} \iint k(x, y) dP(x) dP(y)
$$

**負定値関数の定義**: 全ての $n$, $\{\omega_1, \ldots, \omega_n\} \subset \Omega$, $\{a_1, \ldots, a_n\} \subset \mathbb{R}$ with $\sum a_i = 0$ に対し
$$
\sum_{i,j} a_i a_j k(\omega_i, \omega_j) \leq 0
$$

**例**:
- $k(\omega, y) = |\omega - y|$ → CRPS
- $k(\mathbf{x}, \mathbf{y}) = \|\mathbf{x} - \mathbf{y}\|$ → Energy score
- 一般のMinkowski距離 $\|・\|_\beta$ ($0 < \beta \leq 2$) → kernel score

**Hoeffding型不等式との関係**:
$$
E|X - Y| \leq \sqrt{2} \sqrt{E(X - Y)^2}
$$

### Quantile and Interval Scores

**分位点予測**: $\alpha$-quantile $\xi_\alpha$ の予測スコア
$$
S_\alpha(\xi, x) = -(\mathbb{1}\{x < \xi\} - \alpha)(\xi - x)
$$

**Interval Score**: $\alpha$-prediction interval $[\ell, u]$ に対し
$$
\text{IS}_\alpha([\ell, u], x) = -(u - \ell) - \frac{2}{\alpha}(\ell - x)\mathbb{1}\{x < \ell\} - \frac{2}{\alpha}(x - u)\mathbb{1}\{x > u\}
$$

**解釈**:
- 区間幅 $(u - \ell)$ にペナルティ（narrow is better）
- カバレッジ違反に追加ペナルティ
- $\alpha$-信頼区間の幅とカバレッジをバランス

### Bayes Factorとの関係

**ベイズ因子**（モデル $M_1$ vs $M_2$）:
$$
\text{BF}_{12} = \frac{p(x | M_1)}{p(x | M_2)} = \frac{\int p(x | \theta_1, M_1) \pi(\theta_1 | M_1) d\theta_1}{\int p(x | \theta_2, M_2) \pi(\theta_2 | M_2) d\theta_2}
$$

**対数スコアとの関係**:
$$
\log \text{BF}_{12} = S_{\log}(P_1, x) - S_{\log}(P_2, x)
$$

ここで $P_i$ は予測密度

### Cross-Validation

**Leave-one-out CV**:
$$
\text{CV}_{\text{LOO}} = \frac{1}{n} \sum_{i=1}^n S(\hat{P}_{-i}, X_i)
$$

$\hat{P}_{-i}$: $i$ 番目を除いたデータで推定した予測分布

**Random-fold CV**（本論文で提案）:
- データをランダムに $K$ 個に分割
- 各fold $k$ で $\text{CV}_k = \frac{1}{n_k} \sum_{i \in \text{fold } k} S(\hat{P}_{-k}, X_i)$
- 分割を $M$ 回繰り返し平均

**利点**: LOOより計算効率的、$K$-fold CVより分散小

### Optimum Score Estimation（最適スコア推定）

**点推定**:
$$
\hat{\theta}_n = \arg\max_\theta \frac{1}{n} \sum_{i=1}^n S(P_\theta, X_i)
$$

- 厳密適正スコアを損失関数として使用
- 最尤推定（対数スコア）の一般化
- $M$-estimation の特殊ケース

**分位点推定**:
$$
\hat{\xi}_\alpha = \arg\max_\xi \frac{1}{n} \sum_{i=1}^n S_\alpha(\xi, X_i)
$$

古典的分位点推定（絶対損失最小化）と一致

**区間推定**:
$$
\hat{I}_\alpha = [\hat{\ell}, \hat{u}] = \arg\max_{[\ell, u]} \frac{1}{n} \sum_{i=1}^n \text{IS}_\alpha([\ell, u], X_i)
$$

幅とカバレッジを同時最適化

## 実証結果

### 気象予測評価（Pacific Northwest）

**データ**: 米国太平洋北西部、2002-2003年の12時間先気温予測
- $n = 365$ 日
- 予測: アンサンブル予測（50メンバー）から予測CDF作成

**評価スコア**:
1. **Mean Absolute Error (MAE)**: 予測中央値の絶対誤差
2. **CRPS**: 予測CDF全体を評価
3. **Log Score**: 予測密度（カーネル密度推定）

**結果**（Table 3の概要）:

| モデル | MAE | CRPS | Log Score |
|--------|-----|------|-----------|
| Naive | 4.52 | -2.45 | -3.78 |
| Model A | 3.87 | -2.12 | -3.21 |
| Model B | 3.95 | -2.08 | -3.15 |

Model B: CRPSで最良、sharpnessとcalibrationのバランス良

**不適正スコアの問題**: Figure 4
- 不適正スコア（例: RMSEの平方根）は誤ったモデル選択を導く
- 特に過度に"sharp"（過信）な予測を過大評価

### Skill Score（相対評価）

$$
\text{SS} = 1 - \frac{\bar{S}_{\text{forecast}}}{\bar{S}_{\text{reference}}}
$$

- 参照予測（climatology等）に対する改善度
- $\text{SS} > 0$: 参照より良い、$\text{SS} < 0$: 悪い
- 気象予報で標準的

## 限界・残された課題

### 著者が述べる限界

1. **計算コスト**: CRPSの数値積分、energy scoreの多重積分
   - 高次元データで計算困難
   - モンテカルロ近似が必要

2. **カーネル選択**: Kernel scoreの負定値関数 $k$ の選択基準が不明確
   - $k(\mathbf{x}, \mathbf{y}) = \|\mathbf{x} - \mathbf{y}\|_\beta$ で $\beta$ をどう選ぶか

3. **時系列データ**: 独立性の仮定
   - 自己相関を持つデータへの拡張は今後の課題

4. **離散-連続混合**: カテゴリカル変数と連続変数が混在する場合のスコア未開発

5. **空間データ**: 空間的に配置された予測点への拡張
   - 各点のスコアをどう集約するか

### 私たちが認識する限界

6. **多変量のscalability**: Energy scoreは $d$ 次元で $O(d)$ 計算
   - 組成データ（$d = 100$等）では実用的か不明

7. **ベイズ推論の不確実性**: パラメータ $\theta$ の事後分布は考慮していない
   - 予測分布 $P_\theta$ の $\theta$ 不確実性の伝播

8. **スパース予測**: 高次元で多くの変数がゼロの場合（例: 産地比率）
   - $L^1$ノルムベースのスコアが適切か

9. **階層モデル**: 階層ベイズモデルの各層での評価
   - hyperparameterのスコアは未議論

10. **時空間予測**: 時間・空間の両方で変動する予測
    - 時空間CRPSの定義

### 実装上の問題

11. **CRPSの数値計算**:
$$
\text{CRPS}(F, x) = \int_{-\infty}^\infty (F(y) - \mathbb{1}\{x \leq y\})^2 dy
$$
- CDFが閉形式でない場合（MCMC標本等）の近似
- 論文にも言及あり（経験CDF使用）

12. **Energy scoreのモンテカルロ近似**:
$$
\widehat{\text{ES}}(F, \mathbf{x}) \approx -\frac{1}{M} \sum_{j=1}^M \|\mathbf{X}_j - \mathbf{x}\| + \frac{1}{2M^2} \sum_{j,j'} \|\mathbf{X}_j - \mathbf{X}_{j'}\|
$$
- $M$ の選択（収束速度 $O(M^{-1/2})$）

13. **Interval scoreの最適化**:
$$
\arg\max_{[\ell, u]} \text{IS}_\alpha([\ell, u], X_i)
$$
- 区間幅 $(u - \ell)$ と外れ値ペナルティのバランス
- 数値最適化の安定性

## 本研究（MMCP）との関係

### 直接的な関連

1. **モデル診断・評価**: MMCPの事後予測チェック
   - 論点4（モデル診断）の中核
   - 予測分布 vs 観測データのスコアリング

2. **組成データの評価**: 産地マーク $\mathbf{w}_i \in \mathbb{S}^{d-1}$
   - Energy scoreの適用: $\mathbf{w}$ は simplex上の点
   - Aitchison距離に基づくkernel score

3. **点過程強度の評価**: LGCP $\lambda(\mathbf{s})$ の事後予測
   - CRPS for点過程: イベント数の予測CDF評価
   - 空間強度場全体のenergy score

### MMCPでの応用可能性

4. **産地比率の予測分布評価**:
   - 予測: $\hat{F}(\mathbf{w} | \mathbf{s})$（ある位置の産地比率）
   - 観測: 実際の黒曜石マーク $\mathbf{w}_{\text{obs}}$
   - スコア: $\text{ES}(F, \mathbf{w}_{\text{obs}})$ on simplex

5. **距離場の予測評価**:
   - トブラー距離 $d_T(\mathbf{s}_i, \mathbf{s}_j)$ の予測
   - 実測距離 vs 予測距離のCRPS

6. **イベント数の分布予測**:
   - 領域 $\mathcal{A}$ 内の遺跡数 $N(\mathcal{A})$
   - 予測分布 $\hat{F}_N$ vs 観測 $n_{\text{obs}}$
   - $\text{CRPS}(\hat{F}_N, n_{\text{obs}})$

7. **クロスバリデーションによるモデル選択**:
   - Random-fold CVでMMCPのハイパーパラメータ選択
   - 産地数 $K$、相関関数パラメータ $\theta$ 等

### MMCPで借用すべき要素

8. **Interval scoreの区間推定**:
   - トブラー距離の95%信頼区間評価
   - 幅 $(u - \ell)$ とカバレッジを同時最適化
   - Bayesian credible intervalの評価基準

9. **Skill score**:
   - MMCPのベースラインモデル（等方性ポアソン過程）に対する改善度
$$
\text{SS}_{\text{MMCP}} = 1 - \frac{\text{CRPS}_{\text{MMCP}}}{\text{CRPS}_{\text{baseline}}}
$$

10. **Optimum score estimation**:
    - 事後モードの代替: 厳密適正スコア最大化推定量
    - 計算効率的、解釈容易

### MMCPで拡張すべき点

11. **Simplex上のEnergy score**:
    - Gneiting論文はユークリッド空間 $\mathbb{R}^d$
    - MMCP: Aitchison距離 $d_a(\mathbf{w}, \mathbf{w}')$ on $\mathbb{S}^{d-1}$
$$
\text{ES}_{\text{simplex}}(F, \mathbf{w}) = -E_F d_a(\mathbf{W}, \mathbf{w}) + \frac{1}{2} E_F d_a(\mathbf{W}, \mathbf{W}')
$$

12. **Marked point processのスコア**:
    - イベント位置 $\mathbf{s}_i$ + マーク $\mathbf{w}_i$ の同時評価
    - 既存研究なし → 新規開発必要

13. **時空間スコア**:
    - 考古学データは時代情報あり
    - 時空間CRPS: $\text{CRPS}(F_{s,t}, (x, y, t))$

14. **階層的スコア**:
    - MMCPの多層構造（強度場 → マーク分布）
    - 各層のスコアをどう統合するか

### 技術的な注意点

15. **MCMCサンプルからのCRPS計算**:
    - 事後サンプル $\{\mathbf{w}^{(m)}\}_{m=1}^M$ → 経験CDF $\hat{F}$
    - 経験CRPS: $\widehat{\text{CRPS}}(\hat{F}, \mathbf{w}_{\text{obs}}) = \frac{1}{M} \sum_{m=1}^M d(\mathbf{w}^{(m)}, \mathbf{w}_{\text{obs}}) - \frac{1}{2M^2} \sum_{m,m'} d(\mathbf{w}^{(m)}, \mathbf{w}^{(m')})$

16. **Log scoreの発散**:
    - 予測密度 $p(\mathbf{w})$ が観測点で0 → $\log p(\mathbf{w}_{\text{obs}}) = -\infty$
    - simplex境界（成分=0）で問題
    - CRPSは境界で有界（ロバスト）

17. **計算量**:
    - CRPS: $O(M)$ (Mはサンプル数)
    - Energy score: $O(M^2)$ (全ペア距離)
    - 大規模MCMCでは間引き必要

## 引用すべき箇所

### Proper scoring ruleの定義

> "A scoring rule is proper if the forecaster maximizes the expected score for an observation drawn from the distribution $F$ if he or she issues the probabilistic forecast $F$, rather than $G \neq F$. It is strictly proper if the maximum is unique. In prediction problems, proper scoring rules encourage the forecaster to make careful assessments and to be honest." (p. 359, Abstract)

### 特徴付け定理

> "Theorem 1. A regular scoring rule $S: \mathcal{P} \times \Omega \rightarrow \bar{\mathbb{R}}$ is proper relative to the class $\mathcal{P}$ if and only if there exists a convex, real-valued function $G$ on $\mathcal{P}$ such that $S(P, \omega) = G(P) - \int G^*(P, \omega) dP(\omega) + G^*(P, \omega)$ for $P \in \mathcal{P}$ and $\omega \in \Omega$, where $G^*(P, \cdot)$ is a subtangent of $G$ at the point $P \in \mathcal{P}$." (p. 362)

### CRPS の定義と解釈

> "The continuous ranked probability score (CRPS) is defined by $\text{CRPS}(F, x) = -\int_{-\infty}^\infty (F(y) - \mathbb{1}\{x \leq y\})^2 dy$. ... The CRPS generalizes the absolute error and can be expressed as $\text{CRPS}(F, x) = E_F|X - x| - \frac{1}{2} E_F|X - X'|$ where $X$ and $X'$ are independent random variables with distribution $F$." (p. 368-369)

### Energy scoreの提案

> "We propose a novel score, the energy score, which applies to probabilistic forecasts in the form of multivariate predictive distributions. ... The energy score can be written as $\text{ES}(F, \mathbf{x}) = E_F\|\mathbf{X} - \mathbf{x}\| - \frac{1}{2} E_F\|\mathbf{X} - \mathbf{X}'\|$ where $\mathbf{X}$ and $\mathbf{X}'$ are independent random vectors with distribution $F$, and $\|\cdot\|$ denotes the Euclidean norm." (p. 371)

### Kernel scoreの一般構成

> "A general construction gives rise to kernel scores based on negative definite functions. ... For a negative definite function $k: \Omega \times \Omega \rightarrow \mathbb{R}$, the kernel score is $S_k(P, \omega) = -\int k(\omega, y) dP(y) + \frac{1}{2} \iint k(x, y) dP(x) dP(y)$." (p. 372)

### Interval scoreの提案

> "We propose the interval score as a utility function that addresses width as well as coverage. For a $\alpha$-prediction interval $[\ell, u]$, the interval score is $\text{IS}_\alpha([\ell, u], x) = (u - \ell) + \frac{2}{\alpha}(\ell - x)\mathbb{1}\{x < \ell\} + \frac{2}{\alpha}(x - u)\mathbb{1}\{x > u\}$." (p. 374)

### Propriety（適正性）の重要性

> "Propriety is essential in scientific and operational forecast evaluation; and we present a case study that provides a striking example of the potential issues that result from the use of intuitively appealing but improper scoring rules." (p. 360)

### Optimum score estimation

> "To fix the idea, suppose that we wish to fit a parametric model $P_\theta$ based on a sample $X_1, \ldots, X_n$. To estimate $\theta$, we might measure the goodness-of-fit by the mean score $\mathcal{S}_n(\theta) = \frac{1}{n} \sum_{i=1}^n S(P_\theta, X_i)$, where $S$ is a strictly proper scoring rule. ... This suggests a general approach to estimation: Choose a strictly proper scoring rule that is tailored to the problem at hand and use $\hat{\theta}_n = \arg\max_\theta \mathcal{S}_n(\theta)$ as the optimum score estimator." (p. 360)

### SharpnessとCalibration

> "The goal of probabilistic forecasting is to maximize the sharpness of the predictive distributions subject to calibration. Calibration refers to the statistical consistency between the distributional forecasts and the observations, and is a joint property of the forecasts and the events or values that materialize. Sharpness refers to the concentration of the predictive distributions and is a property of the forecasts only." (p. 359)

### Bregman divergence との関係

> "If the sample space is finite and the entropy function is sufficiently smooth, then the divergence function becomes the Bregman divergence associated with the convex function $G$. Bregman divergences play major roles in optimization and have recently attracted the attention of the machine learning community." (p. 363)

### 実証評価での教訓

> "A case study on probabilistic weather forecasts in the North American Pacific Northwest illustrates the importance of propriety. ... Improper scoring rules can lead to incorrect ranking of competing forecast procedures and encourage forecasters to issue forecasts that are overly sharp (overconfident)." (Case study section, p. 375-376)
