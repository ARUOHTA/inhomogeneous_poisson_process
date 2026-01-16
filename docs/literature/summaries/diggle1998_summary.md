# Diggle et al. (1998) - Model-based geostatistics

## 基本情報
- **タイトル**: Model-based geostatistics (with discussion)
- **著者**: Peter J. Diggle, Jonathan A. Tawn, Rosemary A. Moyeed
- **ジャーナル**: Journal of the Royal Statistical Society: Series C (Applied Statistics), 47(3), 299-350
- **年**: 1998
- **DOI**: 10.1111/1467-9876.00113
- **関連論点**: 論点1（空間統計モデル）、論点3（presence-onlyデータとIPP）

## 主な貢献

クリギング（kriging）の枠組みを非ガウシアンデータに拡張し、一般化線形モデル（GLM）と空間統計を統合した**モデルベース地球統計学**を提案。

1. **一般化線形予測（Generalized Linear Prediction）**: 標準的なガウシアンクリギングを、ポアソン、二項、その他の指数型分布族に拡張
2. **ベイズMCMC推論**: パラメータ不確実性を完全に組み込んだ予測・推論手法
3. **バリオグラムの拡張**: 非ガウシアンデータに対するバリオグラムの理論的導出
4. **実データへの応用**: 放射能汚染（Rongelap島）、感染症疫学（Campylobacter）への適用

## 手法の概要

### モデル構造

**(a) 潜在ガウシアン過程**:
$$
S(\mathbf{x}) \sim \text{GP}(0, \sigma^2 \rho(\mathbf{x} - \mathbf{x}'))
$$

**(b) 条件付き独立性**:
$$
Y_i | S(\mathbf{x}_i) \sim f_i\{y \mid S(\mathbf{x}_i)\} \text{ (mutually independent)}
$$

**(c) リンク関数**:
$$
h\{M_i\} = S(\mathbf{x}_i) + \mathbf{d}_i^T \boldsymbol{\beta}
$$
ここで $M_i = E[Y_i | S(\mathbf{x}_i)]$

### 一般化線形予測

予測量は以下の多重積分で表される：
$$
\hat{S}(\mathbf{x}_j^*) = \frac{E_{1,j}}{f(\mathbf{y})}
$$
$$
E_{1,j} = \int s_{n+j} \left\{ \prod_{i=1}^n f_i(y_i | s_i) \right\} g_{n+m}(\mathbf{s}) \, d\mathbf{s} \, ds_{n+j}
$$

### MCMC推論アルゴリズム

**Step 0**: 初期値設定（GLM推定値を使用）

**Step 1**: パラメータ $\boldsymbol{\theta}$ の更新
- $\pi(\boldsymbol{\theta} | \mathbf{S}) \propto p(\mathbf{S} | \boldsymbol{\theta}) p(\boldsymbol{\theta})$

**Step 2**: 潜在過程 $\mathbf{S}$ の更新
- $\pi(S_i | \mathbf{S}_{-i}, \mathbf{Y}, \boldsymbol{\theta}, \boldsymbol{\beta}) \propto f(y_i | s_i, \boldsymbol{\beta}) \cdot p(S_i | \mathbf{S}_{-i}, \boldsymbol{\theta})$
- Metropolis-Hastings受理確率: $\Delta(S_i, S_i') = \min\{f(y_i|s_i', \boldsymbol{\beta}) / f(y_i|s_i, \boldsymbol{\beta}), 1\}$

**Step 3**: 回帰パラメータ $\boldsymbol{\beta}$ の更新
- $\pi(\boldsymbol{\beta} | \mathbf{Y}, \mathbf{S}) \propto \prod_{j=1}^n f(y_j | s_j, \boldsymbol{\beta}) \cdot p(\boldsymbol{\beta})$

**Step 4**: 予測点 $\mathbf{S}^*$ のサンプリング
- $\mathbf{S}^* | (\mathbf{S}, \boldsymbol{\theta}) \sim \text{MVN}(\Sigma_{12}^T \Sigma_{11}^{-1} \mathbf{S}, \Sigma_{22} - \Sigma_{12}^T \Sigma_{11}^{-1} \Sigma_{12})$

### バリオグラム

非ガウシアンデータに対するバリオグラムの一般形：
$$
\begin{aligned}
C(\mathbf{u}) = & \frac{1}{2} [E_S\{\tau^2(\mathbf{x}) + \tau^2(\mathbf{x}+\mathbf{u})\} + \text{var}_S\{M(\mathbf{x})\} + \text{var}_S\{M(\mathbf{x}+\mathbf{u})\}] \\
& - \text{cov}_S\{M(\mathbf{x}), M(\mathbf{x}+\mathbf{u})\}
\end{aligned}
$$

**ポアソンの場合**: $M(\mathbf{x}) = \exp\{\beta + S(\mathbf{x})\}$, $\tau^2(\mathbf{x}) = M(\mathbf{x})$
$$
C(\mathbf{u}) = \exp(\beta + \sigma^2/2) + \exp(2\beta + \sigma^2) [\exp(\sigma^2) - \exp\{\sigma^2 \rho(\mathbf{u})\}]
$$

小さい $\sigma$ に対する近似:
$$
C(\mathbf{u}) \approx \tau_*^2 + \sigma_*^2 \{1 - \rho(\mathbf{u})\}
$$
ここで $\tau_*^2 = \exp(\beta + \sigma^2/2)$, $\sigma_*^2 = \sigma^2 \tau_*^4$

### 相関関数の仮定

すべての応用例で以下を使用：
$$
\rho(u) = \exp\{-(\alpha u)^\delta\}, \quad \alpha > 0, \, 0 < \delta < 2
$$

## 実証結果

### 6.1 シミュレーション研究（1次元）

- **データ**: 150地点、ポアソン分布 $\log\{M(t)\} = \beta_0 + \beta_1 d(t) + S(t)$
- **真値**: $\boldsymbol{\theta} = (1, 20, 1.8)$, $\boldsymbol{\beta} = (-1, 1)$
- **結果**:
  - 事後平均が真のシグナルを良好に推定
  - GLMのみ（$S=0$を仮定）では信頼区間が狭すぎ、推定値がバイアス
  - $\beta_1$ の事後平均 ≈ 1.0（真値と一致）
  - 収束のために再パラメータ化が必要: $\beta_0^* = \beta_0 + \bar{s}$

### 6.2 Rongelap島の放射能汚染

- **データ**: 157地点の $^{137}$Cs 計数（ポアソン）
- **モデル**: $M_i = t_i \exp\{\beta + S(\mathbf{x}_i)\}$
- **事後モード**: $\boldsymbol{\theta} = (1.7, 0.65, 4.7, 0.7)$ for $(\beta, \sigma, \alpha, \delta)$
- **結果**:
  - $\delta < 1$ により、対数ガウシアンクリギング（$\delta=2$）より空間変動が大きい
  - パラメータ不確実性を無視すると予測標準偏差が約80%過小評価
  - 最大強度の事後分布: 平均15-20 counts/time（点推定11では過小）
  - 汚染域の95%信頼領域を地図化（Figure 10）

### 6.3 Campylobacter感染（イングランド北部）

- **データ**: 248地点の二項データ（234 Campylobacter / 399全腸内感染）
- **モデル**: $\text{logit}\{P(\mathbf{x})\} = \beta + S(\mathbf{x})$
- **デクラスタリング**: 同一郵便番号・5日以内のケースを1ケースとして扱う
- **事後モード**: $\boldsymbol{\theta} = (0.62, 1.0, 6.5)$ for $(\beta, \sigma, \alpha)$（$\delta=1$に固定）
- **結果**:
  - 残差逸脱度 361.4 (df=247, $p \approx 3 \times 10^{-6}$) → 二項過分散有意
  - 推定対数オッズ範囲: (-0.46, 1.42) → Campylobacter比率 0.39-0.81
  - 農村部で高リスク（農場が環境貯蔵庫の可能性）

## 限界・残された課題

### 著者が述べる限界

1. **計算コスト**: $n$ 個の多重積分の評価が必要、MCMCによる近似が不可欠
2. **収縮問題**: パラメータ空間で0が吸収状態となる可能性（実際には問題なし）
3. **収束判定**: 視覚的・時系列診断に依存、より形式的な方法の必要性
4. **再パラメータ化の必要性**: $\beta_0$ と $\bar{s}$ の直交化が収束に不可欠
5. **非定常性の問題**: 短距離空間依存 $S_1$ と長距離 $S_2$ の混合モデルでバイアスの可能性

### 私たちが認識する限界

6. **単変量データのみ**: 多変量（マルチスピーシーズ）への拡張は議論されていない
7. **定常性・等方性の仮定**: 強い仮定、非定常性への拡張が必要
8. **離散化の問題**: 連続空間過程を有限点で推論、空間分解能の変換が非自明
9. **共変量なしの応用**: Rongelap、Campylobacterとも環境共変量を含まない
10. **時空間への拡張**: Campylobacterは本来時空間感染過程だが空間マージナルのみ分析

## 本研究（MMCP）との関係

### 直接的な関連

1. **潜在ガウシアン過程の枠組み**: MMCPも潜在強度過程 $\Lambda(\mathbf{x})$ を想定
   - Diggle: $\Lambda = \exp\{Y\}$ で $Y \sim$ GP
   - MMCP: 対数ガウシアンコックス過程（LGCP）の枠組みと一致

2. **ベイズMCMC推論**: 本研究のMetropolis-Hastingsアルゴリズムは、MMCPのPolyá-Gamma augmentation以前の標準的手法
   - Step 2のシグナル更新: $S_i | \mathbf{S}_{-i}, \mathbf{Y}, \boldsymbol{\theta}$ のサンプリング
   - 条件付きガウシアン $p(S_i | \mathbf{S}_{-i}, \boldsymbol{\theta})$ の利用

3. **パラメータ不確実性**: $\boldsymbol{\theta}$ と $\boldsymbol{\beta}$ の同時推論
   - Rongelap例で示された「固定vs変動パラメータ」の比較はMMCPでも重要

4. **バリオグラム診断**: 非ガウシアンデータのバリオグラム理論
   - MMCPのマーク付き点過程でも類似の二次モーメント診断が必要

### 拡張すべき点

5. **点過程データへの適用**: Diggleは格子データ（固定点での観測）
   - MMCPは点パターン自体が確率的（イベント位置もランダム）
   - 強度関数 $\lambda(\mathbf{x})$ の推定が主目的

6. **マルチスピーシーズへの拡張**: Diggleは単変量
   - MMCPは黒曜石産地マーク $\mathbf{w} \in \mathbb{S}^{d-1}$（組成データ）
   - 依存Dirichlet過程による産地混合の扱いが必要

7. **presence-onlyデータへの応用**: Diggleは完全なcount/prevalenceデータ
   - MMCPは考古学遺跡（発見バイアスあり）
   - weighted likelihoodやcase-control approachの統合が課題

### 本研究が借用できる要素

8. **相関関数パラメータ化**: $\rho(u) = \exp\{-(\alpha u)^\delta\}$
   - $\delta < 2$ により滑らかさ制御（mean square continuity）
   - $\delta=2$ では数値不安定（condition number $\approx 10^{60}$）

9. **診断プロット**:
   - Figure 5のバリオグラム許容区間（シミュレーション包絡線）
   - パラメータ時系列プロット（Figure 2）による収束診断

10. **予測量の柔軟性**: 最大値、超過確率、excursion setなどの関数形
    - MMCPでも「ある産地が主要供給源である確率」などに応用可能

## 引用すべき箇所

### GLMとkrigingの統合

> "The distinctive feature of kriging methodology derived from the Gaussian assumptions (a) and (b) in Section 2 is that the predictor $\hat{S}(\mathbf{x})$ is linear in the data Y. We now extend this methodology in exactly the same way that generalized linear models (Nelder and Wedderburn, 1972; McCullagh and Nelder, 1989) extend the classical Gaussian linear model." (p. 307)

### ベイズ枠組みの利点

> "The Bayesian framework also provides a convenient way of incorporating parameter uncertainty into predictive inferences for the process $S$. We adopt independent proper uniform priors throughout, a consequence being that the resulting joint posterior distribution is proportional to the likelihood surface over the specified region." (p. 308)

### MCMCの必然性

> "Since $n$ will be very large, it is clear that special methods will be needed for the approximate evaluation of $f(\mathbf{y})$, $E_{1,j}$ and $E_{2,j}$. ... In Section 4, we attack this problem by using MCMC methods (Smith and Roberts, 1993; Besag and Green, 1993; Gilks et al., 1993), and adopting a Bayesian approach to inference." (pp. 307-308)

### 条件付き vs 周辺パラメータ

> "The explanation for this is that, within a single realization, there is partial confounding between the deterministic trend, $\beta_0 + \beta_1 t$, and the smooth stochastic variation about the trend, $S(t)$. This emphasizes that the $\boldsymbol{\beta}$-parameters must be interpreted conditionally on $S(t)$, rather than marginally." (p. 315)

### パラメータ不確実性の重要性

> "For the Rongelap data it is difficult to distinguish visually the maps of the posterior mean of $\lambda(\mathbf{x})$ produced with fixed and with varying parameters. However, potentially important differences do exist. ... However, for the prediction sites, at which data are not observed, estimates obtained with varying parameters are on average approximately 5% larger, rarely being smaller, and have prediction standard deviations which are 80% larger on average." (p. 320)

### モデル仮定の両刃性

> "As with all model-based approaches to inference, generalized linear prediction is a two-edged sword, requiring its users to address the assumptions made more critically than in the case of nonparametric smoothing methods." (p. 327)

### 空間vs離散モデル

> "This follows the standard practice of classical geostatistics and allows us to make predictions at arbitrary locations within the study region without redefining our underlying model. When the spatial index set can be represented by a finite number of locations, Markov random field models would be computationally more convenient..." (p. 326)

### バリオグラムと非ガウシアン性

> "Despite the similarities between our model for generalized linear prediction and the generalized linear model, the diagnostic checks for the latter (McCullagh and Nelder, 1989) are inappropriate here owing to the spatial nature of the problem (Handcock and Wallis, 1994). However, a diagnostic for the second-moment structure of the fitted model is valuable and this is given by a comparison of the fitted and empirical variograms." (p. 312)

### $\delta$ パラメータの解釈

> "The change from $\delta > 1$ to $\delta < 1$ is most important as the estimates of the underlying process $S(\mathbf{x})$ no longer smooth the observed spatial variations in the data so heavily." (p. 319)

### Discussionより（John Kent）

> "The parameter $\delta$, $0 < \delta \leq 2$, determines the fractal dimension of the path of $X(t)$ by $D = 2 - \frac{1}{2}\delta$, through the behaviour of the covariance function near the origin, $\rho(u) = 1 - O(|u|^\delta)$ as $u \to 0$." (p. 335, Discussion)
