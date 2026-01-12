# Rue et al. (2009) - INLA: 高速ベイズ近似推定

## 基本情報
- **タイトル**: Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations
- **著者**: Håvard Rue, Sara Martino, Nicolas Chopin
- **ジャーナル**: Journal of the Royal Statistical Society: Series B (Statistical Methodology)
- **年**: 2009
- **DOI**: 10.1111/j.1467-9868.2008.00700.x
- **関連論点**: 論点1（空間統計の計算効率化）

## 主な貢献
1. **INLA (Integrated Nested Laplace Approximations) の提案**: MCMCに代わる決定論的近似手法で、潜在ガウシアンモデルの事後周辺分布を高速・高精度に計算
2. **圧倒的な計算効率**: MCMCが数時間〜数日要する推定を、数秒〜数分で完了（精度はMCMCより高い）
3. **広範な適用範囲**: GLM、時系列、空間モデル、Log-Gaussian Cox過程など、構造化加法回帰モデル全般に適用可能
4. **自動化とブラックボックス化**: チューニング不要、並列計算対応、汎用プログラム `inla` として実装

## 手法の概要

### 潜在ガウシアンモデル (Latent Gaussian Models)
構造化加法回帰モデルのサブクラス：

$$
\eta_i = \alpha + \sum_{j=1}^{n_f} f^{(j)}(u_{ji}) + \sum_{k=1}^{n_\beta} \beta_k z_{ki} + \varepsilon_i
$$

- 潜在場 $\mathbf{x}$: $n$ 次元ガウシアン変数（$n = 10^2$〜$10^5$）
- ハイパーパラメータ $\boldsymbol{\theta}$: 少数（$m \leq 6$）
- 応答変数: 非ガウシアン（Poisson、Binomial、Student-$t$ など）

### GMRF (Gaussian Markov Random Fields)
- **条件付き独立性**: 精度行列 $\mathbf{Q}(\boldsymbol{\theta})$ がスパース
- **Cholesky分解**: $\mathbf{Q} = \mathbf{L}\mathbf{L}^T$ の計算量は次元に依存
  - 1次元: $\mathcal{O}(n)$
  - 2次元: $\mathcal{O}(n^{3/2})$
  - 3次元: $\mathcal{O}(n^2)$
- 周辺分散の計算も効率的: 空間データで $\mathcal{O}\{n \log(n)^2\}$

### INLAの3段階アプローチ

#### Step 1: $\tilde{\pi}(\boldsymbol{\theta} \mid \mathbf{y})$ の計算
Tierney & Kadane (1986) のLaplace近似：

$$
\tilde{\pi}(\boldsymbol{\theta} \mid \mathbf{y}) \propto \left. \frac{\pi(\mathbf{x}, \boldsymbol{\theta}, \mathbf{y})}{\tilde{\pi}_G(\mathbf{x} \mid \boldsymbol{\theta}, \mathbf{y})} \right|_{\mathbf{x} = \mathbf{x}^*(\boldsymbol{\theta})}
$$

- $\tilde{\pi}_G$: ガウス近似
- $\mathbf{x}^*(\boldsymbol{\theta})$: 条件付き分布のモード
- 誤差率: $\mathcal{O}(n_d^{-3/2})$（規準化後）

**探索戦略**:
- モード $\boldsymbol{\theta}^*$ を準ニュートン法で特定
- ヘッセ行列から標準化変数 $\mathbf{z}$ へ変換（スケール・回転補正）
- 対数密度が閾値以下に落ちるまで各座標方向を探索

#### Step 2: $\tilde{\pi}(x_i \mid \boldsymbol{\theta}, \mathbf{y})$ の近似（3つの手法）

**A. ガウス近似 (Gaussian)**:
$$
\tilde{\pi}_G(x_i \mid \boldsymbol{\theta}, \mathbf{y}) = \mathcal{N}\{x_i; \mu_i(\boldsymbol{\theta}), \sigma_i^2(\boldsymbol{\theta})\}
$$

- 最も高速だが、位置誤差や歪度の欠如

**B. Laplace近似 (Laplace)**:
$$
\tilde{\pi}_{LA}(x_i \mid \boldsymbol{\theta}, \mathbf{y}) \propto \left. \frac{\pi(\mathbf{x}, \boldsymbol{\theta}, \mathbf{y})}{\tilde{\pi}_{GG}(\mathbf{x}_{-i} \mid x_i, \boldsymbol{\theta}, \mathbf{y})} \right|_{\mathbf{x}_{-i} = \mathbf{x}_{-i}^*(x_i, \boldsymbol{\theta})}
$$

- モード最適化の代わりに条件付き期待値で近似: $\mathbf{x}_{-i}^* \approx E_{\tilde{\pi}_G}(\mathbf{x}_{-i} \mid x_i)$
- "Region of interest" $R_i(\boldsymbol{\theta})$: 相関が強いノードのみ考慮

$$
R_i(\boldsymbol{\theta}) = \{j: |a_{ij}(\boldsymbol{\theta})| > 0.001\}
$$

**C. 簡略化Laplace近似 (Simplified Laplace)**:
Laplace近似を $x_i = \mu_i(\boldsymbol{\theta})$ 周りでテイラー展開：

$$
\log\{\tilde{\pi}_{SLA}(x_i^{(s)} \mid \boldsymbol{\theta}, \mathbf{y})\} = c - \frac{1}{2}(x_i^{(s)})^2 + \gamma_i^{(1)}(\boldsymbol{\theta}) x_i^{(s)} + \frac{1}{6}(x_i^{(s)})^3 \gamma_i^{(3)}(\boldsymbol{\theta}) + \ldots
$$

ここで、
- $\gamma_i^{(1)}$: 位置補正項（1次）
- $\gamma_i^{(3)}$: 歪度補正項（3次）
- Skew normal分布にフィット（対称分布には別処理）

**計算量**: 全 $n$ 個の周辺分布を計算する総コストは $\boldsymbol{\theta}$ の次元に指数的 × $\mathcal{O}\{n^2 \log(n)\}$（空間データの場合）

#### Step 3: 数値積分
$$
\tilde{\pi}(x_i \mid \mathbf{y}) = \sum_k \tilde{\pi}(x_i \mid \boldsymbol{\theta}_k, \mathbf{y}) \tilde{\pi}(\boldsymbol{\theta}_k \mid \mathbf{y}) \Delta_k
$$

選択した $\{\boldsymbol{\theta}_k\}$ 点での値を重み付き和で積分

### 誤差評価

**有効パラメータ数**:
$$
p_D(\boldsymbol{\theta}) \approx n - \mathrm{tr}\{\mathbf{Q}(\boldsymbol{\theta}) \mathbf{Q}^*(\boldsymbol{\theta})^{-1}\}
$$

- $p_D$ が小さい → 近似誤差が小さい
- 残差 $r(\mathbf{x}; \boldsymbol{\theta}, \mathbf{y})$ を1000サンプルで評価

## 実証結果

### 1. シミュレーション例（Section 5.1）
**データ**: AR(1) 潜在場（$n=50$）+ Student-$t_3$ / Bernoulli観測
- Simplified Laplace vs Gaussian の SKLD: 0.20 (Student-$t$), 0.05 (Bernoulli)
- Simplified Laplace vs (Full) Laplace の SKLD: 0.001, 0.0004（ほぼ一致）
- **計算時間**: 全近似計算 < 0.08秒 vs MCMC 25秒

### 2. 縦断データGLMM（Section 5.2, Epilてんかんデータ）
**データ**: 59患者 × 4訪問 = 236観測、Poisson応答、$n=301$
- 過分散モデル: 個人ランダム効果 + 訪問ランダム効果
- 最大SKLD: 0.23（$\beta_0$ でGaussian vs Simplified Laplace）
- 有効パラメータ数: $p_D = 121.1$（約 2 samples/parameter）
- **計算時間**: 約1.5秒 vs OpenBUGS 数時間

### 3. 確率的ボラティリティモデル（Section 5.3, ポンド・ドル為替レート）
**データ**: $n_d = 945$ 日次対数差分
- AR(1) + log-variance、Student-$t_\nu$ 拡張あり
- 部分データ（$n=50$）での比較: INLA 0.3秒 vs EP法（MATLAB）40分（かつ精度低い）
- 全データ（$n=945$）: $p_D \approx 63$、**計算時間**: 約11秒

### 4. 疾病マッピング（Section 5.4, 旧東ドイツ子宮頸がん）
**データ**: 6990症例、216地区 × 15年齢層、$n=3241$
- Logistic回帰: 年齢効果（2次RW）+ 空間効果（intrinsic CAR）+ 地区ランダム効果
- 最大SKLD: 0.058
- $p_D \approx 101$（$\ll n_d$）
- **計算時間**: 約34秒

### 5. Log-Gaussian Cox過程（Section 5.5, 熱帯雨林樹木データ）
**データ**: 3605本の木、50ヘクタール、200×100格子（$n_d = 20000$）
- 共変量: 標高、勾配 + 空間効果（2次多項式GMRF）+ 非構造効果
- 大規模離散化でも高速推定
- **スケーラビリティ**: 数万観測点でも実用的

### 比較結果まとめ
| 例 | $n$ | 計算時間（INLA） | 計算時間（MCMC） | 精度 |
|----|-----|------------------|------------------|------|
| シミュレーション | 50 | < 0.08秒 | 25秒 | ほぼ完全一致 |
| Epilデータ | 301 | 1.5秒 | 数時間 | 完全一致 |
| 為替レート | 945 | 11秒 | - | 完全一致 |
| 疾病マッピング | 3241 | 34秒 | 数時間 | 完全一致 |
| Cox過程 | 20000 | 数分 | 数日（推定） | 高精度 |

## 限界・残された課題

### 著者が述べる限界
1. **ハイパーパラメータ数 $m$ への依存**: 計算コストが $m$ に対して指数的
   - $m \leq 6$ では実用的だが、$m = 10$ では困難
   - 対策: CCDアプローチ（composite conditional design）、経験ベイズ近似

2. **非ガウス観測への依存**: 対数尤度が極端に非ガウシアンな場合、近似精度が低下
   - 厚い裾を持つ分布（Student-$t$ with low d.f.）では追加工夫が必要

3. **事前分布の影響**: Gaussianity の保持度合いが $p_D$ に反映される
   - 非情報事前分布では $p_D \to 0$ で誤差なし
   - 情報的データでは $p_D \ll n_d$ が理想

### 追加の課題
1. **モデル選択の自動化**: どの効果を空間変動させるか、どの共変量を含めるかの選択支援
2. **非定常共分散の扱い**: GMRFの枠組みを超える一般的な共分散構造への拡張
3. **動的モデルへの拡張**: 時空間モデルでの逐次推定（filtering/smoothing）
4. **ベイズ因子の計算コスト**: 周辺尤度の計算は可能だが、モデル比較の自動化が必要

## 本研究(MMCP)との関係

### 直接的な関係
1. **Log-Gaussian Cox過程への応用**: 黒曜石遺跡の点過程モデルはLGCPの特殊ケース
   - INLA論文のSection 5.5が直接的な先行研究
   - 空間強度 $\lambda(s) = \exp\{\eta(s)\}$ のモデル化手法

2. **計算効率化の必然性**: MMCPは多変量応答 + マーク付き点過程で $n$ が大きい
   - 黒曜石データ: 数百〜数千遺跡 × 複数産地 → 数万パラメータ
   - MCMCでは実行不可能な規模でも、INLA的アプローチなら実用的

3. **空間構造のモデル化**: GMRF（intrinsic CAR、RW priors）が直接利用可能
   - 距離産地効果の空間変動をGMRFでモデル化
   - トブラー距離などの地形情報を共変量として組み込み

4. **ハイパーパラメータの推定**: 空間range、smoothness、ナゲット効果の同時推定
   - MMCPでも $\boldsymbol{\theta}$ = (空間decay, 分散成分, Dirichlet濃度) が典型的に $m \leq 6$

### 本研究への示唆
- **MCMCからの脱却可能性**: マーク付き点過程でもINLA的アプローチで高速化
  - 点位置のLGCP + 組成のDirichlet回帰を階層的に近似
  - Simplified Laplace approximationで位置・歪度補正

- **モデル診断の実用性**: $p_D$、残差評価で近似誤差を定量化
  - MMCPの複雑なモデルでも、事後分布の妥当性を数値的に検証可能

- **並列計算の活用**: OpenMPによるマルチコア並列化
  - 各遺跡 $i$ の周辺分布 $\tilde{\pi}(x_i \mid \mathbf{y})$ は独立に計算可能

### 借用すべき手法
1. **GMRF理論とスパース行列計算**: Rue & Held (2005) の実装
2. **Laplace近似の入れ子構造**: $\boldsymbol{\theta}$ の周辺化 → $x_i$ の周辺化
3. **Simplified Laplace の skew normal fitting**: 計算量 $\mathcal{O}(n^2 \log n)$ の効率
4. **R-INLAパッケージ**: 2013年以降、`INLA` パッケージとして公開

### 拡張すべき点
1. **マーク（組成）の扱い**: 論文はスカラー応答だが、MMCPは多変量
   - 組成ベクトル $\mathbf{c}_i$ への拡張（Aitchison幾何学との統合）
   - Stick-breaking + Pólya-Gamma augmentationとの併用

2. **点位置とマークの同時モデリング**: 2段階階層構造
   - Stage 1: 点位置 $\mathbf{s}_i \sim \mathrm{LGCP}(\lambda)$
   - Stage 2: マーク $\mathbf{c}_i \mid \mathbf{s}_i \sim \mathrm{Dirichlet}(\boldsymbol{\alpha}(\mathbf{s}_i))$

3. **距離依存効果**: トブラー距離 $d_{ij}$ を共変量として組み込み
   - $\eta_i = \beta_0 + \sum_j \beta_j f_j(d_{ij}) + w(s_i)$

## 引用すべき箇所

### INLA手法の優位性（Abstract, p.1）
> "We show that, by using an integrated nested Laplace approximation and its simplified version, we can directly compute very accurate approximations to the posterior marginals. The main benefit of these approximations is computational: where Markov chain Monte Carlo algorithms need hours or days to run, our approximations provide more precise estimates in seconds or minutes."

**用途**: MCMC代替手法としてのINLAの圧倒的計算効率を示す際に引用

### 潜在ガウシアンモデルの定義（Section 1.3, p.4）
> "Latent Gaussian models are a subset of all Bayesian additive models with a structured additive predictor, namely those which assign a Gaussian prior to α, {f^(j)(·)}, {β_k} and {ε_i}."

**用途**: MMCPが潜在ガウシアンモデルの枠組みに属することを明示する際に引用

### GMRFの計算量（Section 2.1, p.6）
> "The typical cost of factorizing Q into LL^T depends on the dimension of the GMRF, e.g. O(n) for one dimension, O(n^{3/2}) for two dimensions and O(n^2) for three dimensions."

**用途**: 空間モデルの計算量を議論する際に引用

### Simplified Laplace approximationの効率（Section 3.2.3, p.11）
> "The simplified Laplace approximation appears to be highly accurate for many observational models. The computational cost is dominated by the calculation of vector a_i.(θ), for each i; thus the 'region of interest' strategy is unhelpful here. Most of the other terms do not depend on i and thus are computed only once. The cost for computing equation (22), for a given i, is of the same order as the number of non-zero elements of the Cholesky triangle, e.g. O{n log(n)} in the spatial case."

**用途**: Simplified Laplace近似の計算効率を具体的に説明する際に引用

### 誤差評価の重要性（Section 4.2, p.13）
> "Obviously, there is only one way to assess with certainty the approximation error of our approach, which is to run an MCMC sampler for an infinite time. However, we propose to use the following two strategies to assess the approximation error, which should be reasonable in most situations."

**用途**: INLAの近似誤差をどう評価すべきかを議論する際に引用

### Log-Gaussian Cox過程への応用（Section 5.5, p.22）
> "A log-Gaussian Cox process is a hierarchical Poisson process: Y in W ⊂ R^d is a Poisson point process with a random-intensity function λ(ξ) = exp{Z(ξ)}, where Z(ξ) is a Gaussian field at ξ ∈ R^d. In this way, the dependence in the point pattern is modelled through a common latent Gaussian variable Z(·)."

**用途**: MMCPの点過程部分がLGCPであることを説明する際に引用

### 実用性と展望（Discussion, p.27）
> "Near instant inference will make latent Gaussian models more applicable, useful and appealing for the end user, who has no time or patience to wait for the results of an MCMC algorithm, notably if he or she must analyse many different data sets with the same model."

**用途**: 高速推定の実践的意義を強調する際に引用

### ソフトウェア実装（Discussion, p.27）
> "A prototype of such a program, inla, is already available (Martino and Rue, 2008) and all the latent Gaussian models in Section 5 were specified and analysed by using this program. inla is built on the GMRFLib library, which is open source and available from the first author's Web page. (An interface to the inla program from R is soon to come.)"

**用途**: R-INLAパッケージの起源と実装を説明する際に引用

### 主要な欠点（Discussion, p.28）
> "The main disadvantage of the INLA approach is that its computational cost is exponential with respect to the number of hyperparameters m. In most applications m is small, but applications where m goes up to 10 do exist."

**用途**: INLA手法の限界を正直に認める際に引用、MMCPでの $m$ の設定を正当化
