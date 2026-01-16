# Blei & Jordan (2006) - ディリクレ過程混合モデルに対する変分推論

## 基本情報
- **タイトル**: Variational Inference for Dirichlet Process Mixtures
- **著者**: David M. Blei, Michael I. Jordan
- **ジャーナル**: Bayesian Analysis
- **年**: 2006
- **DOI**: 記載なし
- **関連論点**: 論点2（ベイズ推論手法、計算効率）

## 主な貢献
- **ディリクレ過程（DP）混合モデルに対する変分推論アルゴリズムの開発**: MCMCサンプリングの決定論的代替手法を提供
- **スティックブレーキング表現の活用**: 無限次元のDPを有限次元の変分分布で近似
- **モデルと変分分布の分離**: モデル自体は打ち切らず、変分分布のみを打ち切ることで真のDPを保持
- **計算効率の実証**: 高次元データにおいてGibbsサンプリングより高速で収束時間の分散が小さい

## 手法の概要

### ディリクレ過程混合モデル
ディリクレ過程（Ferguson 1973）は測度上の測度であり、以下のように定義される：

$$
\begin{aligned}
G | \{G_0, \alpha\} &\sim \text{DP}(G_0, \alpha) \\
\eta_n &\sim G, \quad n \in \{1, \ldots, N\}
\end{aligned}
$$

スティックブレーキング表現（Sethuraman 1994）により、$G$は離散測度として明示的に表現される：

$$
\begin{aligned}
\pi_i(\mathbf{v}) &= v_i \prod_{j=1}^{i-1}(1-v_j), \quad v_i \sim \text{Beta}(1, \alpha) \\
G &= \sum_{i=1}^{\infty} \pi_i(\mathbf{v}) \delta_{\eta_i^*}, \quad \eta_i^* \sim G_0
\end{aligned}
$$

### 変分推論アプローチ
真の事後分布$p(\mathbf{w}|\mathbf{x}, \theta)$を変分分布$q_\nu(\mathbf{w})$で近似し、KLダイバージェンスを最小化：

$$
\text{D}(q_\nu(\mathbf{w}) \| p(\mathbf{w}|\mathbf{x}, \theta)) = \text{E}_q[\log q_\nu(\mathbf{W})] - \text{E}_q[\log p(\mathbf{W}, \mathbf{x}|\theta)] + \log p(\mathbf{x}|\theta)
$$

これは対数周辺尤度の下界の最大化と等価：

$$
\log p(\mathbf{x}|\theta) \geq \text{E}_q[\log p(\mathbf{W}, \mathbf{x}|\theta)] - \text{E}_q[\log q_\nu(\mathbf{W})]
$$

### 打ち切りスティックブレーキング表現
変分分布を打ち切りレベル$T$で有限次元化（$q(v_T=1)=1$）：

$$
q(\mathbf{v}, \boldsymbol{\eta}^*, \mathbf{z}) = \prod_{t=1}^{T-1} q_{\gamma_t}(v_t) \prod_{t=1}^{T} q_{\tau_t}(\eta_t^*) \prod_{n=1}^{N} q_{\phi_n}(z_n)
$$

ここで：
- $q_{\gamma_t}(v_t)$: ベータ分布
- $q_{\tau_t}(\eta_t^*)$: 指数型分布族（自然パラメータ$\tau_t$）
- $q_{\phi_n}(z_n)$: 多項分布

**重要**: 打ち切りは変分分布にのみ適用され、モデル自体は完全なディリクレ過程のまま。

### 座標降下アルゴリズム
指数型分布族の条件付き分布$p(w_i|\mathbf{w}_{-i}, \mathbf{x}, \theta)$に対して、変分パラメータの更新式は：

$$
\nu_i = \text{E}_q[g_i(\mathbf{W}_{-i}, \mathbf{x}, \theta)]
$$

DP混合モデルの具体的な更新式：

$$
\begin{aligned}
\gamma_{t,1} &= 1 + \sum_n \phi_{n,t} \\
\gamma_{t,2} &= \alpha + \sum_n \sum_{j=t+1}^{T} \phi_{n,j} \\
\tau_{t,1} &= \lambda_1 + \sum_n \phi_{n,t} x_n \\
\tau_{t,2} &= \lambda_2 + \sum_n \phi_{n,t} \\
\phi_{n,t} &\propto \exp(S_t)
\end{aligned}
$$

ここで、

$$
S_t = \text{E}_q[\log V_t] + \sum_{i=1}^{t-1} \text{E}_q[\log(1-V_i)] + \text{E}_q[\eta_t^*]^T X_n - \text{E}_q[a(\eta_t^*)]
$$

### 予測分布
変分事後分布を用いた予測分布の近似：

$$
p(x_{N+1}|\mathbf{x}, \alpha, \lambda) \approx \sum_{t=1}^{T} \text{E}_q[\pi_t(\mathbf{V})] \text{E}_q[p(x_{N+1}|\eta_t^*)]
$$

## 実証結果

### シミュレーション実験（Gaussian DP混合モデル）
- **データ**: 5次元から50次元、各次元100データポイント + 100 held-out points
- **共分散構造**: AR(1)型（$\rho=0.9$）、強い依存構造
- **打ち切りレベル**: $T=K=20$

**収束時間の比較**:
- 変分推論が最も高速で、分散が小さい
- collapsed Gibbs sampler > TDP Gibbs sampler の順
- 変分推論の収束時間は次元に対してほぼ一定（5次元〜50次元）

**Held-out対数確率** (平均と標準誤差):

| 次元 | 変分推論 | Collapsed Gibbs | TDP Gibbs |
|------|----------|-----------------|-----------|
| 5    | -147.96 (4.12) | -148.08 (3.93) | -147.93 (3.88) |
| 10   | -266.59 (7.69) | -266.29 (7.64) | -265.89 (7.66) |
| 20   | -494.12 (7.31) | -492.32 (7.54) | -491.96 (7.59) |
| 30   | -721.55 (8.18) | -720.05 (7.92) | -720.02 (7.96) |
| 40   | -943.39 (10.65) | -941.04 (10.15) | -940.71 (10.23) |
| 50   | -1151.01 (15.23) | -1148.51 (14.78) | -1147.48 (14.55) |

→ 変分推論の予測性能はMCMCとほぼ同等

### 大規模画像データ解析
- **データ**: Associated Press画像5000枚、192次元（8×8グリッドのRGB平均値）
- **モデル**: DP mixture of Gaussians（平均$\mu$、共分散$\sigma^2 I$）
- **打ち切りレベル**: $T=150$
- **事前分布**: $\alpha \sim \text{Gamma}(1,1)$、$1/\sigma^2 \sim \text{Gamma}(4,2)$、$\mu \sim \mathcal{N}(0, 5\sigma^2)$

**結果**:
- 収束時間: 約4時間
- 使用コンポーネント数: 79個（150個中）
- collapsed Gibbs: 1反復15分 → 4時間で16反復のみ（収束不十分）

発見されたクラスタ例：
- バスケットボールのシュート
- 曇天の屋外シーン
- 顔写真
- 青背景の画像

## 限界・残された課題

### 理論的限界
1. **局所最適解の問題**: 変分パラメータ空間に複数の局所最適解が存在する可能性
   - 対策: 複数回の初期化とリスタート
   - トレードオフ: 計算コストの増加

2. **近似の限界**: 変分分布は真の事後分布の近似にとどまる
   - 完全因子化仮定により事後依存構造を無視
   - より複雑な変分分布（階層的変分法など）でギャップを縮小可能だが計算コスト増

3. **理論的保証の弱さ**: MCMCの収束保証に比べて近似精度の保証が限定的

### 実装上の課題
4. **打ち切りレベル$T$の選択**:
   - $T$が小さすぎる: 近似精度の低下
   - $T$が大きすぎる: 計算コストの増加
   - 論文では$T$を固定したが、KLダイバージェンスに基づく最適化も可能

5. **非共役基底分布への拡張**: $G_0$が共役でない場合、座標降下アップデートが複雑化
   - 共役分布の混合の場合は単純なアップデートが可能

6. **初期化の重要性**: sequential importance samplingによる初期化を推奨
   - データのランダム順列による逐次アップデート
   - 最良の下界を与える初期化を選択

### 収束診断の課題
7. **MCMCとの比較**:
   - MCMC: 収束診断が困難（定常分布への到達判定）
   - 変分推論: 最適化基準（Equation 16）で収束評価可能だが、局所最適解の可能性

## 本研究（MMCP）との関係

### 借用すべき手法
1. **変分推論フレームワーク**: MMCP推論における計算効率向上の選択肢
   - MCMCが収束しない/遅い場合の代替手段
   - 決定論的アルゴリズムによる再現性の向上

2. **スティックブレーキング表現**: 無限混合モデルの扱い
   - コンポーネント数を事前に固定しないノンパラメトリックアプローチ
   - データから自動的にクラスタ数を推定

3. **座標降下最適化**: 指数型分布族に対する効率的な更新式
   - Equation 15の一般形を他のモデルに適用可能

### 拡張すべき点
1. **空間過程への適用**:
   - 本論文はi.i.d.混合モデル
   - MMCPは空間相関を持つCox過程 → 空間構造の導入が必要

2. **マーク過程への拡張**:
   - 本論文は単変量観測（画像ピクセル、Gaussianデータ）
   - MMCPは多変量マーク（黒曜石組成） → マルチタスク変分推論の検討

3. **大規模データへの適用**:
   - 遺跡データが数千規模になった場合、変分推論が有効
   - 特に高次元共変量（標高、距離行列など）を扱う際に優位性

4. **打ち切りレベルの自動選択**:
   - MMCPにおける最適な打ち切りレベル$T$の決定法
   - クロスバリデーション、情報量規準（WAIC、WBIC）の活用

### MMCPにおける位置づけ
- **第3章での役割**: ノンパラメトリックベイズ推論の計算手法の1つとして紹介
- **実装の選択肢**: Stan/NUTSの代替として変分推論を検討する際の理論的基盤
- **トレードオフの理解**: MCMC vs 変分推論の速度・精度トレードオフを議論する材料

## 引用すべき箇所

### MCMCの限界と代替手法の必要性
> "However, while an unquestioned success story, MCMC is not an unqualified one—MCMC methods can be slow to converge and their convergence can be difficult to diagnose." (p.1, Introduction)

→ MMCP推論でMCMCの問題点を述べる際の引用

> "MCMC sampling can be prohibitively slow, and it is important to explore alternatives." (p.1, Abstract)

→ 大規模データに対する変分推論の動機づけ

### 変分推論の基本原理
> "The basic idea of variational inference is to formulate the computation of a marginal or conditional probability in terms of an optimization problem." (p.1, Introduction)

→ 変分推論の概念説明

### 打ち切りの扱い
> "Note that our use of truncation is rather different [from Ishwaran & James 2001]. In our case, the model is a full Dirichlet process and is not truncated; only the variational distribution is truncated." (p.7, Section 3.2)

→ モデルと変分分布の打ち切りの違いを強調（重要な設計思想）

### 実証結果
> "The variational algorithm was faster and exhibited significantly less variance in its convergence time." (p.17, Section 5)

→ 変分推論の計算効率の実証的優位性

> "Though it is based on an approximation to the posterior, the resulting predictive distributions are very accurate for this class of DP mixtures." (p.16, Section 5)

→ 近似にもかかわらず予測性能が高いことの強調

### 大規模データへの適用
> "The variational algorithm required approximately four hours to converge. [...] For a rough comparison to Gibbs sampling, an iteration of collapsed Gibbs takes 15 minutes with this data set. In the same four hours, one could perform only 16 iterations. This is not enough for a chain to converge to its stationary distribution." (p.18, Section 6)

→ 大規模データ（5000画像、192次元）における実用的優位性

### 座標降下アルゴリズムの導出
> "The optimization of KL divergence with respect to a single variational parameter $\nu_i$ is achieved by computing the following expectation: $\nu_i = \text{E}_q[g_i(\mathbf{W}_{-i}, \mathbf{x}, \theta)]$." (p.6, Section 3.1, Equation 15)

→ 平均場変分推論の一般的更新式（Appendix A で詳細導出）

### 方法論の比較
> "Both variational and MCMC methods have strengths and weaknesses, and it is unlikely that one methodology will dominate the other in general." (p.19, Conclusions)

→ MMCP推論手法の選択を議論する際のバランスの取れた視点
