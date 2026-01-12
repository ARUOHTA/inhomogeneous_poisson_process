# Gonçalves & Gamerman (2018) - 厳密ベイズ推論によるCox過程

## 基本情報
- **タイトル**: Exact Bayesian inference in spatiotemporal Cox processes driven by multivariate Gaussian processes
- **著者**: Flávio B. Gonçalves, Dani Gamerman
- **ジャーナル**: Journal of the Royal Statistical Society: Series B (Statistical Methodology)
- **年**: 2018（Received March 2015. Revised April 2017）
- **DOI**: 記載なし
- **関連論点**: 論点1（空間的依存性）、論点3（点過程）

## 主な貢献

1. **厳密なベイズ推論の実現**
   - 離散化誤差なし（discretization-free）のベイズ推論
   - Adams et al. (2009)の単変量空間モデルを多変量時空間へ一般化
   - MCMC誤差のみが唯一の不確実性源（離散化誤差を完全に回避）

2. **多変量動的GPによる強度関数のモデリング**
   - 複数のGP成分による共変量効果・時間依存性の表現
   - 空間変化係数（spatially varying covariate effects）への対応
   - 動的GP（difference equation）による時間進化

3. **Poisson thinning + augmentationの理論的基盤**
   - データ拡張により尤度の非可解性を回避
   - Gibbs samplerによる直接サンプリング
   - Skew normal分布からのサンプリング手法の開発

4. **計算アルゴリズムの洗練**
   - 完全条件付き分布からの直接サンプリング
   - Metropolis-Hastingsの調整不要（高速収束・良好なmixing）
   - 強度関数の任意位置での推定手法

## 手法の概要

### モデル構造

**基本モデル**（空間・時空間Cox過程）:
$$
\begin{aligned}
Y &\sim \text{Cox process}(\lambda_t(s)) \\
\lambda_t(s) &= \lambda_t^* \Phi\{W_t(s)\beta_t(s)\}
\end{aligned}
$$

ここで：
- $\Phi$: 標準正規分布の累積分布関数
- $\beta_t(s) = (\beta_{0,t}(s), \ldots, \beta_{q,t}(s))$: $q+1$個の独立GP（多変量GP）
- $W_t(s)$: 共変量行列（$f\{\beta_t(s), W_t(s)\} = W_t(s)\beta_t(s)$と線形化）
- $\lambda_t^*$: 強度の上限（supremum）

**多変量GPの例**（共変量モデル）:
$$
f\{\beta_t(s), W_t(s)\} = \beta_{0,t}(s) + \beta_{1,t}(s)W_{1,t}(s) + \ldots + \beta_{q,t}(s)W_{q,t}(s)
$$

**動的GP**（離散時間）:
$$
\beta_{t'}(\cdot) = G_{t',t}\beta_t(\cdot) + w_{t',t}(\cdot), \quad w_{t',t} \sim \text{GP}(\text{zero mean, time independent})
$$
- $G$: 時間遷移行列
- $w_{t',t}$: 動的誤差GP（恒等行列の場合、独立な時間進化）

### データ拡張モデル

Adams et al. (2009)のアイデアを時空間・多変量へ拡張：

1. **潜在過程の導入**:
   - $X_t$: 均質ポアソン過程（強度$\lambda_t^*$、$S$上）
   - $Y_t$: 観測された点過程（thinning後）
   - 各点$s_k \in X_t$に対して、$Z_{t,k} \sim \text{Ber}(\Phi\{W_t(s_k)\beta_t(s_k)\})$でthinning

2. **Poisson thinningアルゴリズム**（Algorithm 1）:
   - Step 1: $K_t \sim \text{Poisson}\{\lambda_t^*\mu(S)\}$を生成
   - Step 2: $K_t$個の点を$S$上に一様分布
   - Step 3: $\beta_t$をサンプリング
   - Step 4: 各点を確率$\lambda_t(s_k)/\lambda_t^*$で保持
   - Step 5: 保持された点が$Y_t$

3. **拡張尤度**:
$$
\begin{aligned}
\pi(\psi | W, S) &= \Phi_N(W_N\beta_N; I_N) \Phi_{K-N}(-W_M\beta_M; I_{K-N}) \pi_{\text{GP}}(\beta_K | \theta) \\
&\quad \times \exp\{-\lambda^*\mu(S)\}(\lambda^*)^K \frac{1}{K!} \pi(\lambda^*)\pi(\theta)
\end{aligned}
$$

ここで：
- $N$: 観測点数（$Y$）
- $K$: 潜在過程の点数（$X$）
- $M = K - N$: thinningされた点数
- $\psi := (\{s_n\}_{n=1}^N, \{s_m\}_{m=1}^{K-N}, \beta_K, \{s_k\}_{k=1}^K, K, \lambda^*, \theta)$

### Gibbsサンプラー

**ブロック構造**（空間モデル）:

1. **$(\{s_m\}, \beta_M, K)$のサンプリング**（Algorithm 3）:
   - $K \sim \text{Poisson}(\lambda^*\mu(S))$ truncated to $\{N, N+1, \ldots\}$
   - 各$m = 1, \ldots, K-N$について：
     - $s_m \sim \mathcal{U}(S)$を提案
     - $\beta(s_m) \sim \pi_{\text{GP}}(\beta(s_m) | \beta_N, \beta_{1:m-1}, \theta)$を生成
     - 確率$\Phi\{-W(s_m)\beta(s_m)\}$でaccept（pointwise rejection sampling）

2. **$\beta_K$のサンプリング**（Algorithm 4）:
   - Skew normal分布$\text{SN}(\mu_K, \Sigma_K, W_K)$からサンプリング
   - $\mu_K, \Sigma_K$: GP事前分布の平均・共分散
   - $W_K$: 観測・未観測点の共変量を統合した行列

3. **$\lambda^*$のサンプリング**:
   - 共役Gamma事前分布$\mathcal{G}(\alpha_\lambda, \beta_\lambda)$の場合
   - 完全条件付き分布: $\mathcal{G}\{\alpha_\lambda + K, \beta_\lambda + \mu(S)\}$

4. **$\theta$のサンプリング**:
   - GPのハイパーパラメータ（共分散関数のパラメータ）
   - Metropolis-Hastingsまたは適応的Gaussian random walk

**Skew Normal分布**:

Arellano-Valle & Azzalini (2006)の一般クラス：
$$
f(z) = \frac{1}{\Phi_m(\gamma; \Gamma)} \phi_d(z-\xi; \Sigma) \Phi_m(Wz; I_m)
$$

Algorithm 2によりサンプリング：
1. Truncated multivariate normalからサンプリング
2. Cholesky分解を用いた変換
3. 条件付き正規分布からサンプリング

**時空間モデルへの拡張**（Section 5）:

- ブロック: $\{(K_t, \{s_{t,m}\}, \beta_{M_t})\}_{t=0}^T$、$\{\beta_{K_t}\}_{t=0}^T$
- 時間依存性の因数分解:
$$
\pi(\{\beta_{K_t}\}_{t=0}^T | \cdot) \propto \prod_{t=0}^T \left\{\Phi_{N_t}(W_{N_t}\beta_{N_t}; I_{N_t}) \Phi_{K_t-N_t}(-W_{M_t}\beta_{M_t}; I_{K_t-N_t}) \pi_{\text{GP}}(\beta_{K_t} | \beta_{K_{(t-1)}}, \theta)\right\}
$$

### 強度関数の推定

**Lemma 3**: 任意の有限集合$S_0 = \{\tilde{s}_1, \ldots, \tilde{s}_G\}$での事後分布：
$$
\pi(\beta_{S_0} | \{s_n\}_{n=1}^N, W, S) = \int \pi(\beta_{S_0} | \beta_K, \theta) \pi(\beta_K, \theta | \{s_n\}_{n=1}^N, W, S) d\beta_K d\theta
$$

**積分強度の推定**（離散化誤差なし）:
$$
\hat{\Lambda}(R) = \mu(R) \frac{1}{J}\sum_{j=1}^J \lambda^{(j)}(U^{(j)}), \quad U^{(j)} \sim \mathcal{U}(R)
$$

これは強法則により強一致推定量。

## 実証結果

### 1. Lansing Woods data
- **データ**: 448本のwhite oak trees（924 ft × 924 ft）
- **設定**:
  - GP: $\gamma = 3/2$指数共分散、$\sigma^2 = 2$、$\tau^2 = 2$
  - 事前分布: $\lambda^* \sim \text{Gamma}(76, 6)$ truncated below 15
- **計算**:
  - 500 MCMC iterations（burn-in 100）
  - 各イテレーション約30秒
  - Monte Carlo誤差: 0.47%（$\Lambda([0,4]^2)$の推定）
- **結果**:
  - 滑らかな強度関数の推定
  - 高速収束・低自己相関（Fig. 3 of supplementary）

**比較**: 離散化手法（lgcpパッケージ）との比較
- 質的に類似だが、離散化レベル（1600, 10000, 40000セル）に依存
- 本手法は真の連続的IF（離散化誤差なし）

### 2. New Brunswick fires
- **データ**: 1987-2003年の森林火災（年ごとの点過程、1988年を除く）
- **モデル**: 動的GP
  - $f\{\beta_t(s), W_t(s)\} = \beta_{0,t}(s)$
  - $\beta_{0,t}(s) = \beta_{0,t-1}(s) + w_t(s)$
  - $\beta_{0,0} \sim \text{GP}(0, 1.75^2, 10)$
  - $w_t \sim \text{GP}(0, 0.5^2, 15)$
- **結果**:
  - 空間的・時間的平滑性
  - 連続する年での類似した空間パターン
  - U字型の高強度領域（時間を通じて持続）

## 限界・残された課題

### 著者が述べる限界

1. **計算コスト**
   - GP共分散行列の逆行列・Cholesky分解が$O(K^3)$
   - $K$（潜在点数）が大きいと計算困難
   - 近似手法（lower dimension, NNGP）が必要

2. **識別可能性**
   - $\lambda^*$と$\beta$の切片の識別問題
   - $\lambda^*$は強度の上限（supremum）として識別すべき
   - 事前分布の慎重な設定が必要

3. **将来の拡張**
   - マーク付き点過程への拡張（言及のみ）
   - 非時空間共変量の導入（Pinto et al. 2015に委ねる）
   - 計算効率化（GP近似、HMC、MALA等）

### 本研究の観点からの限界

1. **Presence-only未対応**
   - 完全観測データを前提
   - Presence-onlyデータのthinningは扱っていない

2. **組成マークの未対応**
   - マークへの拡張は将来の方向性として言及
   - 組成データ特有の扱いは範囲外

3. **Pólya-Gamma拡張の未使用**
   - ロジスティック部分（$\Phi$リンク）はskew normalサンプリング
   - Pólya-Gammaによる効率化は適用されていない

4. **計算スケーラビリティ**
   - フルGPでは大規模データに対応困難
   - NNGP等の近似が必須

## 本研究（MMCP）との関係

### 借用する要素

1. **データ拡張による厳密ベイズ推論の枠組み**
   - Poisson thinningとaugmented modelの構造
   - 潜在過程による尤度の可解化
   - Gibbs samplerによる直接サンプリング

2. **多変量GPによる柔軟なモデリング**
   - 複数成分（共変量・時間効果等）の同時モデリング
   - 動的GPによる時間依存構造
   - 空間変化係数の概念

3. **階層ベイズフレームワーク**
   - 第1段：データ | 潜在過程、パラメータ
   - 第2段：潜在過程（GP）| パラメータ
   - 第3段：パラメータの事前分布

4. **強度関数の推定手法**
   - 任意位置でのGPサンプリング（Lemma 3の応用）
   - 積分強度のMonte Carlo推定（離散化誤差なし）

### MMCPによる拡張

1. **Presence-only + thinningモデルの統合**
   - Gonçalvesはthinningモデル（完全観測）
   - MMCPはpresence-onlyデータへの対応（Moreira & Gamermanと統合）
   - 観測バイアスの明示的モデリング

2. **組成マークの追加**
   - Gonçalvesは点の位置のみ
   - MMCPは組成値マーク（黒曜石産地比率）を追加
   - Multinomial logitリンク

3. **Pólya-Gamma拡張**
   - Gonçalvesはskew normal sampling（embedded Gibbs）
   - MMCPはmultinomial logit + Pólya-Gamma
   - 計算効率の向上

4. **計算効率化**
   - GonçalvesはフルGP（$O(K^3)$の計算コスト）
   - MMCPはNNGPによる効率化（$O(nm^3)$、$m$=近傍数）

### 位置づけ

Gonçalves & Gamerman (2018)は：
- Cox過程の厳密ベイズ推論の理論的基盤を確立
- データ拡張とPoisson thinningによる計算可能な推論
- 多変量GP・動的GP・空間変化係数の概念を導入
- Adams et al. (2009)の単変量空間モデルを一般化

MMCPは：
- この厳密推論の枠組みを継承
- Presence-onlyデータへの対応（Moreira & Gamerman 2022と統合）
- 組成値マークの追加（Eckardt et al. 2025と統合）
- Pólya-GammaとNNGPによる計算効率化
- Moreira et al. (2024)のマーク付きpresence-onlyと概念的に近い

## 引用すべき箇所

### 厳密推論の動機（Section 1）
> "Solutions for the inference problem have required, until recently, the use of discrete approximations [...]. These represent a considerable source of error and, therefore, ought to be used with care [...]. This motivates the development of exact methodologies, i.e. free from discretization errors, and helps to understand its advantages."

> "In this paper, the term exact refers to the fact that no discretization-based approximation is used. In particular, the methodology proposed here has Markov chain Monte Carlo (MCMC) error (Markov chain convergence plus Monte Carlo error) as the only source of inaccuracy, which is generally well understood and controlled."

### データ拡張の目的（Section 1）
> "The aim of this work is to propose an exact inference methodology for spatiotemporal Cox processes in which the IF dynamics are driven by a GP. The exactness feature stems from an augmented model approach as in Adams et al. (2009). However, we generalize their point pattern models by firstly considering spatiotemporal models and, secondly, by using multivariate (possibly dynamic) GPs to allow the inclusion of different model components (regression and temporal effects) in a flexible manner."

### Poisson thinningの役割（Section 2.1）
> "The first advantage of the formulation in expressions (1)-(4) is that it allows exact simulation of data from the model, which is the key to develop exact inference methods. Exact simulation of the model is based on a key result from Poisson processes called Poisson thinning."

### augmented modelの構造（Section 2.2）
> "The crucial step to develop exact methods is to avoid dealing with likelihood (6). One possible solution is to define an augmented model for $Y$ and some additional variable $X$, such that the joint (pseudo)likelihood based on $(X, Y)$ is tractable. [...] Most importantly, this approach leads to a tractable likelihood when the joint distribution of $X$ and $Y$ is considered."

### 動的GPの定義（Section 2.3）
> "A process $\beta$ follows a dynamic GP in discrete time if it can be described by a difference equation $\beta_{t'}(\cdot) = G_{t',t}\beta_t(\cdot) + w_{t',t}(\cdot), \quad w_{t',t} \sim \text{GP}$, where the multivariate GP disturbances $w_{t',t}(\cdot)$ are zero mean and time independent."

### 識別可能性の議論（Section 3.3）
> "The model proposed may suffer from identifiability problems concerning parameter $\lambda^*$. The natural way to identify it is to have this parameter as the supremum of the IF which, under the Bayesian approach, should be achieved by an appropriate specification of the prior distribution."

### 計算コストの問題と対策（Section 3.3）
> "Despite their great flexibility in a variety of statistical modelling problems, GPs have a considerable practical limitation when it comes to computational cost. More specifically, simulating an $n$-dimensional GP has a cost which is typically of the order of $n^3$. [...] Nevertheless, this issue is mitigated by several reasons."

> "Furthermore, for the cases where the cost is still too high, some (approximating) strategies may be employed to reduce it. Most importantly, none of these defy the exactness of our methodology. Lower dimension approximations [...] may be used to deal with the GP prior and to speed the computation of covariance matrices."

### Gibbsサンプラーの性質（Section 4.2）
> "The blocking and sampling schemes of our Gibbs sampling result in good convergence properties. In fact, the examples that are presented here suggest that convergence is attained after a few iterations."

### 将来の方向性（Section 7）
> "An immediate extension of our models involves consideration of marks to the Poisson events. These marks may be described with a variety of components, whose effects are allowed to vary smoothly, in line with the models that are used for the IF."
