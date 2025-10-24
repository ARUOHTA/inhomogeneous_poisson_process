<script>
MathJax = {
  tex: {
    inlineMath: [['$', '$']],
    displayMath: [['$$', '$$']]
  }
};
</script>
<script type="text/javascript" id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
</script>

# 距離先行情報付き NNGP 多項比率モデル：階層的GP事前アプローチ

本ドキュメントでは、`docs/nngp_multinomial_model.md` で導入した多項ロジット + NNGP モデルに「産地までの距離」を明示的な先行情報として組み込む新しいアプローチを定式化する。距離情報を**GPの事前平均として組み込む階層的定式化**により、以下を実現する：

1. **データがある場所**: 距離事前から出発し、観測データに基づいてGPが局所的に調整
2. **データがない場所**: GPの不確実性が大きく、距離事前分布が支配的
3. **解釈性**: 「距離効果をregress outした残差をGPで推定」という直感的な構造

---

# 1. データと記法

- 観測領域: $\mathcal{D} \subset \mathbb{R}^2$。
- 遺跡（観測地点）: $s_i \in \mathcal{D}$（$i=1,\ldots,n$）。
- 産地（距離の基準点）: $c_k \in \mathcal{D}$（$k=1,\ldots,K$）。
- カウントデータ: $\mathbf{y}_i = (y_{ik})_{k=1}^K$、$y_{ik} \ge 0$、$N_i = \sum_{k=1}^K y_{ik}$。
- 共変量ベクトル: $\mathbf{W}(s) = (1, W_1(s),\ldots,W_p(s))^{\top}$（先頭は切片1）。
- 距離Zスコア行列: $\mathbf{Z} \in \mathbb{R}^{n \times K}$、$Z_{ik}$ は地点 $s_i$ から産地 $k$ への距離のZスコア（小さいほど近い）。

---

# 2. 重み付き逆ソフトマックスによる距離ベース事前確率

## 2.1 距離ベース事前確率の定義

距離のZスコア行列 $\mathbf{Z}$ と産地重要度 $\mathbf{w} = (w_1,\ldots,w_K)^{\top}$（各成分 $> 0$）から、各地点 $s_i$ での**距離ベース事前確率**を次のように定義する：

$$
p_{0k}(s_i) = \frac{\exp(-Z_{ik}/\tau + \alpha \log w_k)}{\sum_{\ell=1}^K \exp(-Z_{i\ell}/\tau + \alpha \log w_\ell)},
\qquad k=1,\ldots,K.
$$

ここで：

- $\tau > 0$: **温度パラメータ**。小さいほど最近接産地への集中度が高まる。
- $\alpha \ge 0$: **重要度の効き具合**。$\alpha=0$ なら重要度を無視し、純粋に距離のみ。
- $w_k$: 産地 $k$ の**重要度**。考古学的に重要な産地に高い重みを付与できる。

## 2.2 対数比（log-ratio）特徴量への変換

多項ロジットの識別性を保つため、基準確率 $p_0(s)$ を**基準カテゴリ $K$ との対数比**に変換する：

$$
g_k(s) := \log p_{0k}(s) - \log p_{0K}(s), \qquad k=1,\ldots,K-1.
$$

この $g_k(s)$ は「距離ベース確率で見たカテゴリ $k$ 対 $K$ の**相対優位度**」を表す。負の値は「$K$ の方が近い」、正の値は「$k$ の方が近い」ことを意味する。

## 2.3 事前確率の性質

1. **正規化**: $\sum_{k=1}^K p_{0k}(s) = 1$ かつ $p_{0k}(s) \in (0,1)$。
2. **相対性・スケール不変性**: 対数比 $g_k$ を通じて使われるため、距離の**全体スケールや一様オフセット**にほぼ不変。
3. **最近接産地優位**: $Z_{ik}$ が最小（最も近い）の産地 $k$ で $g_k(s) > 0$ となる傾向。
4. **固定ハイパーパラメータ**: 本モデルでは $\tau, \alpha, \mathbf{w}$ を**固定**とする（必要に応じて階層化も可能）。

---

# 3. 階層的モデルの定式化

## 3.1 線形予測子の構成

各カテゴリ $k=1,\ldots,K-1$ について、線形予測子を次のように定義する：

$$
\eta_k(s) = \mathbf{W}(s)^{\top} \boldsymbol{\beta}_k(s),
$$

ここで $\boldsymbol{\beta}_k(s) = (\beta_{0k}(s), \beta_{1k}(s), \ldots, \beta_{pk}(s))^{\top}$ は空間的に変化する係数ベクトル。

基準カテゴリ $K$ のロジットはゼロ。全カテゴリの確率は softmax により：

$$
\begin{aligned}
\pi_k(s)
  &= \frac{\exp(\eta_k(s))}
           {1 + \sum_{\ell=1}^{K-1} \exp(\eta_\ell(s))},
  &&k = 1,\ldots,K-1, \\
\pi_K(s)
  &= \frac{1}
           {1 + \sum_{\ell=1}^{K-1} \exp(\eta_\ell(s))}.
\end{aligned}
$$

## 3.2 距離事前を組み込んだGP事前分布

**鍵となるアイデア**：距離情報を切片 $\beta_{0k}(s)$ の**事前平均**として組み込む。

### 切片項の事前分布（距離情報を含む）

$$
\beta_{0k}(\cdot) \sim \operatorname{GP}\bigl(\mu_k(\cdot), C_{\theta_k}(\cdot, \cdot)\bigr),
$$

ここで事前平均関数を：

$$
\mu_k(s) = \lambda_k \cdot g_k(s)
$$

と定義する。$\lambda_k$ は**距離事前の強さ**を制御するスケーリングパラメータ（固定またはサンプリング）。

### 共変量係数の事前分布（標準的なGP）

共変量係数 $\beta_{jk}(s)$（$j=1,\ldots,p$）は従来通りゼロ平均のGP：

$$
\beta_{jk}(\cdot) \sim \operatorname{GP}\bigl(0, C_{\theta_k}(\cdot, \cdot)\bigr), \quad j=1,\ldots,p.
$$

### 共分散関数

$C_{\theta}$ は典型的には RBF カーネル：

$$
C_{\theta}(s, s') = \sigma^2 \exp\left(-\frac{\|s - s'\|^2}{2\ell^2}\right).
$$

## 3.3 モデルの解釈

この定式化により：

1. **切片 $\beta_{0k}(s)$ の事前平均** = $\lambda_k \cdot g_k(s)$（距離ベースの期待値）
2. **データがある場所**: 事後分布がデータに引っ張られ、$\beta_{0k}(s)$ が距離事前から局所的に調整される
3. **データがない場所**: 事後分布が事前平均 $\lambda_k \cdot g_k(s)$ に留まり、距離事前が支配的
4. **直感的解釈**: 「距離効果を baseline とし、その周りをGPが局所的に変動」

従来の加算モデル $\eta_k(s) = \lambda_k g_k(s) + \mathbf{W}(s)^{\top} \boldsymbol{\beta}_k(s)$ では、データのない場所でもGP項 $\boldsymbol{\beta}_k$ の事前分布（平均0）が影響し、距離事前が希釈されてしまう問題があった。本アプローチではこの問題が解決される。

---

# 4. 観測尤度

遺跡 $i$ のカウントは多項分布に従う：

$$
\mathbf{y}_i \mid N_i, \boldsymbol{\pi}(s_i)
  \sim \operatorname{Multinomial}\bigl(N_i, \boldsymbol{\pi}(s_i)\bigr).
$$

対数尤度は（定数項を除き）：

$$
\log p(\mathbf{y} \mid \boldsymbol{\eta})
= \sum_{i=1}^n \sum_{k=1}^{K-1} y_{ik} \eta_k(s_i)
   - \sum_{i=1}^n N_i \log\left(1 + \sum_{\ell=1}^{K-1} \exp(\eta_\ell(s_i))\right),
$$

ここで $\eta_k(s_i) = \mathbf{W}(s_i)^{\top} \boldsymbol{\beta}_k(s_i)$ である。

---

# 5. 事前分布の詳細

## 5.1 距離事前強度パラメータ $\lambda_k$

距離事前をどれだけ信頼するかを制御する $\lambda_k$ について、以下の選択肢がある：

### オプション A: 固定値（ハイパーパラメータ）

$\lambda_k$ を固定値として扱う（例: $\lambda_k = 1.0$）。これは最もシンプルで、事前知識を直接反映できる。

### オプション B: サンプリング（弱情報事前）

全カテゴリで共通の係数 $\lambda$ を仮定し、弱情報事前を設定：

$$
\lambda_k \equiv \lambda \sim \mathcal{N}(1, \sigma_\lambda^2),
$$

ここで $\sigma_\lambda^2$ は比較的小さく設定し、$\lambda \approx 1$（距離事前を全面的に採用）に緩やかに誘導する。

### オプション C: カテゴリ別階層事前

各カテゴリごとに $\lambda_k$ を設定し、階層事前で安定化：

$$
\begin{aligned}
\lambda_k \mid \mu_\lambda, \tau_\lambda &\sim \mathcal{N}(\mu_\lambda, \tau_\lambda^2), \quad k=1,\ldots,K-1, \\
\mu_\lambda &\sim \mathcal{N}(1, s^2), \\
\tau_\lambda &\sim \operatorname{Half-Normal}(\tilde{\sigma}).
\end{aligned}
$$

**本実装では**: まずオプションAで固定値として実装し、後にサンプリングに拡張可能とする。

## 5.2 空間変化係数の事前分布（NNGP近似）

観測点集合 $S = \{s_1,\ldots,s_n\}$ 上での切片ベクトル：

$$
\boldsymbol{\beta}_{0k,S} = (\beta_{0k}(s_1),\ldots,\beta_{0k}(s_n))^{\top}
  \sim \mathcal{N}(\boldsymbol{\mu}_k, \Sigma_{\theta_k}),
$$

ここで：

- $\boldsymbol{\mu}_k = (\mu_k(s_1),\ldots,\mu_k(s_n))^{\top} = \lambda_k (g_k(s_1),\ldots,g_k(s_n))^{\top}$（距離ベース事前平均）
- $\Sigma_{\theta_k}$ は共分散行列（RBFカーネルから構成）

共変量係数については従来通り：

$$
\boldsymbol{\beta}_{jk,S} = (\beta_{jk}(s_1),\ldots,\beta_{jk}(s_n))^{\top}
  \sim \mathcal{N}(\mathbf{0}, \Sigma_{\theta_k}), \quad j=1,\ldots,p.
$$

大規模データでは Vecchia/NNGP 近似を用いる：

$$
p(\boldsymbol{\beta}_{jk,S})
  \approx \prod_{i=1}^n p\bigl(\beta_{jk}(s_i) \mid \boldsymbol{\beta}_{jk}(S_{N_i})\bigr),
$$

ここで $S_{N_i}$ は $s_i$ の近傍点集合（Morton順序で前方の点のみ）。

---

# 6. Pólya–Gamma 補助変数

## 6.1 恒等式の適用

多項ロジット尤度の非線形項を線形化するため、Polson et al. (2013) の Pólya–Gamma 恒等式を用いる：

$$
\frac{(\exp(\psi))^a}{(1+\exp(\psi))^b}
  = 2^{-b} \exp(\kappa \psi)
    \int_0^\infty \exp\left(-\frac{\omega \psi^2}{2}\right)
                 p_{\mathrm{PG}}(\omega \mid b,0)\,d\omega,
$$

ここで $\kappa = a - b/2$、$p_{\mathrm{PG}}$ は Pólya–Gamma 分布密度。

## 6.2 補助変数の導入

各遺跡 $i$、カテゴリ $k$ について補助変数を導入：

$$
\omega_{ik} \sim \operatorname{PG}(N_i, \eta_k(s_i)), \quad k=1,\ldots,K-1.
$$

すると条件付き対数尤度は：

$$
\log p(\mathbf{y}_i \mid \boldsymbol{\eta}(s_i), \boldsymbol{\omega}_i)
= \sum_{k=1}^{K-1}
   \left(
     \kappa_{ik} \eta_k(s_i)
     - \tfrac{1}{2}\,\omega_{ik}\,\eta_k(s_i)^2
   \right)
   + \text{const},
$$

ここで：

$$
\kappa_{ik} = y_{ik} - \tfrac{N_i}{2}.
$$

この形は $\eta_k(s_i)$ について二次形式であり、ガウス過程事前と組み合わせることでガウス型の完全条件付き分布が得られる。

---

# 7. ギブスサンプラーと完全条件付き分布

## 7.1 事後分布の分解

補助変数を含む完全データの事後分布は：

$$
p(\boldsymbol{\beta}, \boldsymbol{\omega} \mid \mathbf{y}, \mathbf{W}, \boldsymbol{\mu})
\propto
  p(\mathbf{y} \mid \boldsymbol{\beta}, \boldsymbol{\omega}, \mathbf{W})
  p(\boldsymbol{\omega} \mid \boldsymbol{\beta}, \mathbf{W})
  p(\boldsymbol{\beta} \mid \boldsymbol{\mu}),
$$

ここで $\boldsymbol{\mu}$ は距離ベース事前平均（$\lambda_k g_k$ から構成、固定）。

ギブスサンプラーは以下の2ステップを反復する（$\lambda_k$ 固定の場合）：

1. $\boldsymbol{\omega}$ のサンプリング（独立）
2. $\boldsymbol{\beta}$ のサンプリング（Vecchia順序付き）

## 7.2 ステップ 1: $\omega_{ik}$ の完全条件付き分布

Pólya–Gamma 分布の性質より：

$$
\omega_{ik} \mid \mathbf{y}, \boldsymbol{\beta}
\sim \operatorname{PG}\bigl(N_i, \eta_k(s_i)\bigr),
$$

ここで：

$$
\eta_k(s_i) = \mathbf{W}(s_i)^{\top} \boldsymbol{\beta}_k(s_i).
$$

各 $\omega_{ik}$ は他の補助変数と独立にサンプル可能。

## 7.3 ステップ 2: $\beta_{jk}(s_i)$ の完全条件付き分布

カテゴリ $k$ の共変量 $j$ の係数 $\beta_{jk}(s_i)$ について、他を固定したとき：

$$
\eta_k(s_i) = w_{ij} \beta_{jk}(s_i) + \eta_k^{(-j)}(s_i),
$$

ここで：
- $w_{ij} = W_j(s_i)$ は共変量 $j$ の値（$j=0$ なら $w_{i0}=1$）
- $\eta_k^{(-j)}(s_i) = \eta_k(s_i) - w_{ij} \beta_{jk}(s_i)$ は「成分 $j$ を除いた線形予測子」

### 7.3.1 尤度項の寄与

Pólya–Gamma 補助後の条件付き対数尤度は：

$$
\log p(\mathbf{y}_i \mid \beta_{jk}(s_i), \text{rest})
\propto \kappa_{ik} \cdot w_{ij} \beta_{jk}(s_i)
        - \tfrac{1}{2} \omega_{ik} \left(w_{ij} \beta_{jk}(s_i) + \eta_k^{(-j)}(s_i)\right)^2.
$$

他のカテゴリのロジットとの関係から生じる補正項：

$$
C_{ik} = \log\left(1 + \sum_{\ell \neq k} \exp(\eta_\ell(s_i))\right),
$$

を導入すると、$\tilde{\kappa}_{ik} = \kappa_{ik} + \omega_{ik} C_{ik}$ として尤度を整理できる。

### 7.3.2 事前の寄与（Vecchia 近似）

Vecchia 近似により、$\beta_{jk}(s_i)$ の事前は近傍点に条件付けられる。

**切片の場合**（$j=0$）:

$$
\beta_{0k}(s_i) \mid \boldsymbol{\beta}_{0k}(S_{N_i})
  \sim \mathcal{N}\Bigl(
    \mu_k(s_i) + a_i^{\top} \bigl[\boldsymbol{\beta}_{0k}(S_{N_i}) - \boldsymbol{\mu}_k(S_{N_i})\bigr],
    d_i
  \Bigr),
$$

ここで：
- $\mu_k(s_i) = \lambda_k g_k(s_i)$ は事前平均
- $\boldsymbol{\mu}_k(S_{N_i}) = \lambda_k (g_k(s_{i_1}),\ldots,g_k(s_{i_m}))^{\top}$ は近傍点での事前平均
- $a_i$ は NNGP 係数ベクトル
- $d_i$ は条件付き分散

**共変量係数の場合**（$j \ge 1$）:

$$
\beta_{jk}(s_i) \mid \boldsymbol{\beta}_{jk}(S_{N_i})
  \sim \mathcal{N}\bigl(a_i^{\top} \boldsymbol{\beta}_{jk}(S_{N_i}),\, d_i\bigr).
$$

### 7.3.3 完全条件付き分布（平方完成）

尤度 + 事前を $\beta_{jk}(s_i)$ について平方完成すると：

$$
\beta_{jk}(s_i) \mid \text{rest}
  \sim \mathcal{N}(\mu_{\text{post}}, \sigma_{\text{post}}^2),
$$

ここで：

$$
\sigma_{\text{post}}^2 = \frac{1}{\omega_{ik} w_{ij}^2 + d_i^{-1}}.
$$

**切片の場合**（$j=0$）:

$$
\mu_{\text{post}} = \sigma_{\text{post}}^2 \left[
  \tilde{\kappa}_{ik} w_{i0} - \omega_{ik} w_{i0} \eta_k^{(0)}(s_i)
  + \frac{\mu_k(s_i) + a_i^{\top} (\boldsymbol{\beta}_{0k}(S_{N_i}) - \boldsymbol{\mu}_k(S_{N_i}))}{d_i}
\right],
$$

ここで $w_{i0} = 1$、$\eta_k^{(0)}(s_i)$ は切片を除いた線形予測子。

**共変量係数の場合**（$j \ge 1$）:

$$
\mu_{\text{post}} = \sigma_{\text{post}}^2 \left[
  \tilde{\kappa}_{ik} w_{ij} - \omega_{ik} w_{ij} \eta_k^{(-j)}(s_i)
  + \frac{a_i^{\top} \boldsymbol{\beta}_{jk}(S_{N_i})}{d_i}
\right].
$$

---

# 8. アルゴリズムのまとめ

## 8.1 前処理

1. **距離Zスコア計算**: 各地点 $s$ から各産地 $k$ への距離を計算し、Zスコア化して $\mathbf{Z}$ を得る。
2. **距離ベース事前平均の計算**:
   - 重み付き逆ソフトマックスにより基準確率 $p_0$ を計算
   - 対数比 $g_k(s) = \log p_{0k}(s) - \log p_{0K}(s)$ を算出
   - 事前平均 $\mu_k(s) = \lambda_k g_k(s)$ を計算（全地点で）
3. **NNGP ファクター構築**: Morton 順序付けにより Vecchia 近似のための $(a_i, d_i)$ を全共変量・全地点で前計算（`FactorCache` に保存）。

## 8.2 MCMC 反復

### 初期化

- $\boldsymbol{\beta}_k$ を距離ベース事前平均 $\boldsymbol{\mu}_k$ で初期化（または全カテゴリの総計に基づくロジット）。
- 線形予測子を計算：$\eta_k(s_i) = \mathbf{W}(s_i)^{\top} \boldsymbol{\beta}_k(s_i)$。

### 各反復で以下を実行

**ステップ 1: Pólya–Gamma 補助変数のサンプリング**

$$
\omega_{ik} \sim \operatorname{PG}\bigl(N_i, \eta_k(s_i)\bigr),
\quad \forall i=1,\ldots,n,\; k=1,\ldots,K-1.
$$

**ステップ 2: 空間変化係数のサンプリング**

各カテゴリ $k=1,\ldots,K-1$ について：

- 各共変量成分 $j=0,\ldots,p$ を Morton 順序で更新。
- 各地点 $s_i$ で：
  1. $C_{ik} = \log\bigl(1 + \sum_{\ell \neq k} \exp(\eta_\ell(s_i))\bigr)$ を計算。
  2. $\tilde{\kappa}_{ik} = \kappa_{ik} + \omega_{ik} C_{ik}$ を計算。
  3. 完全条件付き分布 $\mathcal{N}(\mu_{\text{post}}, \sigma_{\text{post}}^2)$ からサンプル。
     - $j=0$ なら距離事前平均 $\mu_k(s_i)$ を含む形で更新
     - $j \ge 1$ なら標準的な更新
  4. $\eta_k(s_i)$ を逐次更新。

### サンプル保存

バーンイン期間後、間引き (thinning) に従って $\boldsymbol{\beta}$ を保存。

## 8.3 予測

### 遺跡位置での構成比

保存したサンプルごとに $\eta_k(s_i) = \mathbf{W}(s_i)^{\top} \boldsymbol{\beta}_k(s_i)$ を計算し、softmax を適用。事後平均を取る。

### グリッド点での予測

グリッド点 $u$ について：

1. 距離ベース事前平均 $\mu_k(u) = \lambda_k g_k(u)$ を計算。
2. Vecchia クロスファクター $(a_*(u), d_*(u))$ を用いて $\boldsymbol{\beta}_k(u)$ を条件付きサンプル：
   - 切片: $\beta_{0k}(u) \sim \mathcal{N}(\mu_k(u) + a_*^{\top}(\boldsymbol{\beta}_{0k}(S_*) - \boldsymbol{\mu}_k(S_*)), d_*)$
   - 共変量: $\beta_{jk}(u) \sim \mathcal{N}(a_*^{\top}\boldsymbol{\beta}_{jk}(S_*), d_*)$
3. $\eta_k(u) = \mathbf{W}(u)^{\top} \boldsymbol{\beta}_k(u)$ を計算し、softmax で確率に変換。

**重要**: データから遠いグリッド点では、$\boldsymbol{\beta}_{0k}(u)$ が事前平均 $\mu_k(u) = \lambda_k g_k(u)$ に近づき、距離事前が支配的になる。

---

# 9. 効果の分解

事後推定結果から、各効果を以下のように分解できる：

## 9.1 基本的な分解

事後平均 $\bar{\boldsymbol{\beta}}_k(s)$ を用いて：

1. **距離事前のみ**: $\eta_k^{\text{dist}}(s) = \lambda_k g_k(s)$
   - これは切片の事前平均に相当

2. **切片のみ**: $\eta_k^{\text{int}}(s) = \bar{\beta}_{0k}(s)$
   - 距離事前 + データによる調整の結果

3. **個別共変量効果**: $\eta_k^{(j)}(s) = W_j(s) \bar{\beta}_{jk}(s)$
   - 各共変量の寄与

4. **全共変量効果**: $\eta_k^{\text{cov}}(s) = \sum_{j=1}^p W_j(s) \bar{\beta}_{jk}(s)$
   - 切片を除く全ての共変量

5. **完全モデル**: $\eta_k^{\text{full}}(s) = \mathbf{W}(s)^{\top} \bar{\boldsymbol{\beta}}_k(s)$
   - 切片 + 全共変量

## 9.2 距離効果の抽出

「純粋な距離効果」を見たい場合：

$$
\eta_k^{\text{dist-only}}(s) = \bar{\beta}_{0k}(s) - \lambda_k g_k(s)
$$

これは「事後切片から事前平均を引いた残差」であり、データによって距離事前からどれだけずれたかを示す。

## 9.3 解釈上の注意

- これらの分解は**非加算的**（softmaxの非線形性のため）
- 各効果の確率への寄与は独立ではなく、相互作用がある
- あくまで「各項を単独で見た場合の傾向」を示す探索的分析

---

# 10. 既存実装との統合

## 10.1 再利用可能なコンポーネント

本モデルは既存の NNGP 多項モデル実装（`bayesian_statistics/nngp/model/`）と高い互換性を持つ：

1. **FactorCache**: Vecchia 近似の前計算機構をそのまま利用。
2. **PolyaGammaSampler**: $\omega_{ik}$ のサンプリングに既存実装を流用。
3. **compute_eta, softmax_with_baseline**: 線形予測子計算と確率変換をそのまま利用。
4. **weighted_inverse_softmax, compute_distance_features**: 距離事前確率計算（既に `sample.py` に実装済み）。

## 10.2 必要な新規実装

### update_beta_category の拡張

既存の `update_beta_category` 関数を以下のように拡張：

1. **引数に事前平均を追加**: `prior_mean_by_feature: Optional[Sequence[np.ndarray]]`
   - 切片（j=0）の場合のみ非ゼロの事前平均 $\mu_k(s)$ を渡す
   - 共変量（j≥1）の場合は None または zeros

2. **NNGP事前の修正**:
   - 切片の場合: $\mu_{\text{prior},i} = \mu_k(s_i) + a_i^{\top}(\beta_{0k}(S_{N_i}) - \mu_k(S_{N_i}))$
   - 共変量の場合: $\mu_{\text{prior},i} = a_i^{\top}\beta_{jk}(S_{N_i})$（従来通り）

### データセット準備の拡張

`prepare_distance_prior_dataset` 関数を新規作成：

1. 既存の `prepare_multinomial_dataset` を呼び出し
2. 距離Zスコアから $g_k(s)$ を計算
3. 事前平均 $\mu_k(s) = \lambda_k g_k(s)$ を計算
4. 拡張された `DistancePriorDataset` に格納

## 10.3 コード量の最小化

- `update_beta_category` の変更は最小限（事前平均オフセットの追加のみ）
- MCMC メインループは `run_mcmc` とほぼ同一（$\lambda$ サンプリングなし）
- 既存の `run_mcmc` を内部で呼び出し、前処理のみ `prepare_distance_prior_dataset` で実施

---

# 11. まとめと比較

## 11.1 従来モデルとの違い

### 旧モデル（加算型）

$$
\eta_k(s) = \lambda_k g_k(s) + \mathbf{W}(s)^{\top} \boldsymbol{\beta}_k(s), \quad \beta_{jk} \sim \operatorname{GP}(0, C)
$$

**問題点**:
- データのない場所で $\beta_{jk}(s) \to 0$（事前平均）となり、距離事前が希釈される
- 産地に極めて近い未観測地点でも $\pi_k \not\to 1$

### 新モデル（階層型）

$$
\eta_k(s) = \mathbf{W}(s)^{\top} \boldsymbol{\beta}_k(s), \quad \beta_{0k} \sim \operatorname{GP}(\lambda_k g_k, C),\; \beta_{jk} \sim \operatorname{GP}(0, C)
$$

**利点**:
- データのない場所で $\beta_{0k}(s) \to \lambda_k g_k(s)$（距離事前）となる
- 産地に近い未観測地点で自然に $\pi_k \to 1$
- 「距離効果をregress out」という直感的解釈

## 11.2 モデルの特徴

1. **理論的クリーンさ**: GP事前の標準的な階層構造
2. **自動的なバランス調整**: データがある場所ではGPが学習、ない場所では距離事前が支配
3. **計算効率**: 加算モデルと同等（事前平均の計算は前処理で完了）
4. **拡張性**: $\lambda_k$ をサンプリング対象にすることも容易

## 11.3 実装上の利点

- **最小限の変更**: `update_beta_category` への小さな拡張のみ
- **コード再利用**: 既存の NNGP インフラを最大限活用（DRY 原則）
- **保守性**: 既存モデルとコードパスを共有するため、バグ修正・機能追加が容易

---

# 参考文献

## 主要文献

- **Polson, N., Scott, J., & Windle, J. (2013)**. *Bayesian inference for logistic models using Pólya–Gamma latent variables*. Journal of the American Statistical Association, 108(504), 1339-1349.
  - Pólya–Gamma 補助変数法の原論文。

- **Finley, A. O., Datta, A., Banerjee, S., & Gelfand, A. E. (2019)**. *Nearest Neighbor Gaussian Processes for massive spatial data sets*. Journal of the American Statistical Association, 114(525), 15-30.
  - Vecchia 近似と NNGP の詳細。

- **Banerjee, S., Carlin, B. P., & Gelfand, A. E. (2014)**. *Hierarchical modeling and analysis for spatial data*. CRC press.
  - 空間統計の階層モデリングの基礎。

## 実装関連

- **既存実装**: `bayesian_statistics/nngp/model/`
  - `multinomial_model.py`: ベースとなる多項 NNGP モデル。
  - `nngp.py`: FactorCache と Vecchia 近似の実装。
  - `sample.py`: Pólya–Gamma サンプラーと係数更新ルーチン。

---

# 12. 拡張: 距離依存分散カーネルによる空間構造の制御

本節では、産地からの距離に応じてGPの分散を変化させる**距離依存分散カーネル**を導入し、黒曜石の運搬経路の連続性をより適切にモデル化する手法を説明する。

## 12.1 動機：運搬経路の連続性と空間的影響範囲

### 問題設定

現在のモデルでは、データがない場所で距離事前が支配的になることを保証しているが、以下の状況が生じる可能性がある：

1. **観測点から離れた場所での孤立した高確率領域**: データがない場所では、距離事前のみが支配的になるため、観測点と空間的に離れた場所でも高い確率を持つ「島」が生じる可能性がある。
2. **運搬経路の連続性が表現されない**: 考古学的には、産地から遠い場所で黒曜石が出土した場合、その運搬経路が必ず存在するが、現在のモデルでは経路の連続性が空間構造に反映されない。

### 考古学的動機

黒曜石の分布には以下の性質がある：

- **産地から遠い場所で黒曜石が出土したら、その運搬経路が必ず存在する**: 例えば、産地から100km離れた場所で黒曜石が出土した場合、その間のどこかで人々が運搬したはずである。
- **遠方の観測点はより広い空間領域に影響を与えるべき**: 産地から遠い場所での観測は、その周辺のより広い範囲で黒曜石が存在する可能性を示唆する。
- **すべてのデータは等しく尊重されるべき**: 距離による重み付けではなく、観測点の空間的影響範囲を変化させることでこの性質を実現する。

この性質をモデルに組み込むため、**距離依存分散カーネル**を導入する。

## 12.2 距離依存分散カーネルの定義

### 基本アイデア

標準的な等方性 RBF カーネルは、一定の分散と長さスケールを持つ：

$$
k(\mathbf{s}_i, \mathbf{s}_j) = \sigma^2 \exp\left(-\frac{\|\mathbf{s}_i - \mathbf{s}_j\|^2}{2\ell^2}\right)
$$

これを拡張し、**産地 $\mathbf{o}$ からの距離に応じて分散が変化する**カーネルを定義する：

$$
k_{\text{dist}}(\mathbf{s}_i, \mathbf{s}_j; \mathbf{o}) = \sigma(\mathbf{s}_i, \mathbf{o}) \cdot \sigma(\mathbf{s}_j, \mathbf{o}) \cdot \exp\left(-\frac{\|\mathbf{s}_i - \mathbf{s}_j\|^2}{2\ell^2}\right)
$$

ここで、$\sigma(\mathbf{s}, \mathbf{o})$ は地点 $\mathbf{s}$ における**局所標準偏差**であり、産地 $\mathbf{o}$ からの距離 $d(\mathbf{s}, \mathbf{o})$ に依存する。

### 局所標準偏差関数の設計

産地からの距離に応じて標準偏差を増加させる関数を定義する：

#### オプション A: 線形スケーリング

$$
\sigma(\mathbf{s}, \mathbf{o}) = \sigma_0 \cdot \bigl(1 + \gamma \cdot d(\mathbf{s}, \mathbf{o})\bigr)
$$

ここで：
- $\sigma_0 > 0$: 基準標準偏差（産地近傍での標準偏差）
- $\gamma \geq 0$: 距離スケーリング係数
- $d(\mathbf{s}, \mathbf{o})$: 産地 $\mathbf{o}$ から地点 $\mathbf{s}$ への距離

#### オプション B: 指数スケーリング

$$
\sigma(\mathbf{s}, \mathbf{o}) = \sigma_0 \cdot \exp\bigl(\gamma \cdot d(\mathbf{s}, \mathbf{o})\bigr)
$$

より急激な距離依存性を表現したい場合に適している。

### カーネルの性質

このカーネルは以下の性質を持つ：

1. **対称性**: $k_{\text{dist}}(\mathbf{s}_i, \mathbf{s}_j; \mathbf{o}) = k_{\text{dist}}(\mathbf{s}_j, \mathbf{s}_i; \mathbf{o})$（GP のカーネルとして有効）
2. **正定値性**: 適切なパラメータ選択により正定値性が保証される
3. **距離依存分散**: 産地から遠い地点ほど、GP の事前分散が大きくなる

## 12.3 距離の標準化

距離 $d(\mathbf{s}, \mathbf{o})$ は、カテゴリ間での比較可能性を保つため、標準化することが望ましい。本モデルでは、既に計算済みの**距離Zスコア** $Z_{ik}$ を用いる：

$$
\sigma(\mathbf{s}_i, \mathbf{o}_k) = \sigma_0 \cdot \bigl(1 + \gamma \cdot Z_{ik}\bigr)
$$

あるいは、元の距離を観測領域のスケールで正規化する：

$$
d_{\text{norm}}(\mathbf{s}, \mathbf{o}) = \frac{d(\mathbf{s}, \mathbf{o})}{D_{\max}}
$$

ここで、$D_{\max}$ は観測領域の直径（例：すべての産地間距離の最大値）。

## 12.4 解釈：大きな分散は重み付けではない

### よくある誤解

「産地から遠い地点で分散を大きくする = そのデータを down-weight する」という解釈は**誤り**である。

### 正しい解釈

GP における分散の意味：

1. **事前不確実性の大きさ**: 分散が大きい地点では、事前分布の不確実性が高い。
2. **空間的影響範囲**: 分散が大きい地点は、より広い範囲の近傍点と強く相関する。
3. **データによる調整の柔軟性**: 分散が大きい地点では、データによって事前分布から大きく調整される可能性がある。

### 考古学的意味

産地から遠い地点で黒曜石が観測された場合：

- **運搬経路の存在を示唆**: その地点への運搬経路が必ず存在する。
- **広い領域への影響**: その観測は、周辺のより広い範囲で黒曜石が存在する可能性を高める。
- **柔軟な空間構造**: 遠方のデータは、距離事前から大きくずれる可能性があるため、より柔軟なフィッティングが必要。

距離依存分散カーネルは、**データの重みを変えるのではなく、データの空間的影響範囲を調整する**。

## 12.5 実装方針

### カーネルクラスの拡張

現在の `LocalNNGPKernel` クラスを拡張し、距離依存分散をサポートする：

```python
@dataclass
class DistanceDependentNNGPKernel:
    """距離依存分散を持つ RBF カーネル。

    産地からの距離に応じて局所標準偏差が変化するカーネル。
    産地から遠い地点ほど分散が大きくなり、より広い空間的影響範囲を持つ。

    Parameters
    ----------
    lengthscale : float
        RBF カーネルの長さスケール（空間相関の範囲）
    base_variance : float
        基準分散（産地近傍での分散）
    distance_scaling : float
        距離スケーリング係数 gamma。
        gamma=0 なら標準的な等方性カーネルに等しい。
    scaling_type : str
        スケーリング関数の種類。"linear" または "exponential"。
    """
    lengthscale: float = 0.15
    base_variance: float = 1.0
    distance_scaling: float = 0.0
    scaling_type: str = "linear"

    def local_std(self, distances: np.ndarray) -> np.ndarray:
        """産地からの距離に基づく局所標準偏差を計算。

        Parameters
        ----------
        distances : np.ndarray
            産地からの距離（Zスコアまたは正規化距離）、shape (n,)

        Returns
        -------
        np.ndarray
            局所標準偏差、shape (n,)
        """
        sigma_0 = np.sqrt(self.base_variance)
        if self.distance_scaling == 0.0:
            return np.full_like(distances, sigma_0)

        if self.scaling_type == "linear":
            return sigma_0 * (1.0 + self.distance_scaling * distances)
        elif self.scaling_type == "exponential":
            return sigma_0 * np.exp(self.distance_scaling * distances)
        else:
            raise ValueError(f"Unknown scaling_type: {self.scaling_type}")

    def K(
        self,
        X1: np.ndarray,
        X2: np.ndarray,
        distances1: np.ndarray,
        distances2: np.ndarray,
    ) -> np.ndarray:
        """距離依存分散カーネル行列を計算。

        Parameters
        ----------
        X1 : np.ndarray
            第1の地点座標、shape (n1, 2)
        X2 : np.ndarray
            第2の地点座標、shape (n2, 2)
        distances1 : np.ndarray
            X1 の各点から産地への距離、shape (n1,)
        distances2 : np.ndarray
            X2 の各点から産地への距離、shape (n2,)

        Returns
        -------
        np.ndarray
            カーネル行列、shape (n1, n2)
        """
        # 基本 RBF カーネル（分散=1で正規化）
        X1 = np.atleast_2d(X1)
        X2 = np.atleast_2d(X2)
        a2 = np.sum(X1**2, axis=1, keepdims=True)
        b2 = np.sum(X2**2, axis=1, keepdims=True).T
        sqd = a2 + b2 - 2.0 * (X1 @ X2.T)
        K_base = np.exp(-0.5 * sqd / (self.lengthscale**2))

        # 局所標準偏差
        sigma1 = self.local_std(distances1)  # (n1,)
        sigma2 = self.local_std(distances2)  # (n2,)

        # 距離依存分散の適用
        sigma_outer = sigma1[:, np.newaxis] * sigma2[np.newaxis, :]  # (n1, n2)
        return K_base * sigma_outer
```

### FactorCache への統合

`build_nngp_factors` 関数を拡張し、カーネルが `DistanceDependentNNGPKernel` の場合は距離情報を渡す：

```python
def build_nngp_factors(
    coords: np.ndarray,
    neighbor_idx: List[np.ndarray],
    kernel: Union[LocalNNGPKernel, DistanceDependentNNGPKernel],
    origin_distances: Optional[np.ndarray] = None,
) -> NNGPFactors:
    """
    NNGP ファクターを構築（距離依存カーネル対応）。

    Parameters
    ----------
    ...
    origin_distances : np.ndarray, optional
        各地点から産地への距離、shape (n,)。
        DistanceDependentNNGPKernel の場合は必須。
    """
    if isinstance(kernel, DistanceDependentNNGPKernel) and origin_distances is None:
        raise ValueError("origin_distances required for DistanceDependentNNGPKernel")

    # ... カーネル行列計算時に distances を渡す
```

## 12.6 パラメータ選択と解釈

### 距離スケーリング係数 $\gamma$ の選択

- **$\gamma = 0$**: 距離依存なし（標準等方性カーネル）
- **$\gamma \in [0.1, 0.5]$**: 穏やかな距離依存（産地から遠いほど若干広い影響範囲）
- **$\gamma > 1.0$**: 強い距離依存（産地から遠い地点は非常に広い空間的影響を持つ）

実際の値は、事前予測チェック (prior predictive check) やクロスバリデーションで選択する。

### 考古学的解釈

- **$\gamma$ が大きい場合**: 産地から遠い場所での観測は、その周辺のより広い範囲に影響を与える。運搬経路の連続性を強く仮定する場合に適切。
- **$\gamma$ が小さい場合**: 距離によらず、観測点の影響範囲はほぼ一定。局所的な交易パターンを仮定する場合に適切。

## 12.7 期待される効果

距離依存分散カーネルの導入により、以下の効果が期待される：

1. **運搬経路の連続性の表現**: 産地から遠い場所で黒曜石が観測された場合、その周辺のより広い範囲で確率が上昇し、経路の連続性が表現される。
2. **孤立した島の抑制**: 遠方の観測点が広い空間的影響を持つため、データと空間的に連続しない高確率領域が生じにくくなる。
3. **すべてのデータの等しい尊重**: データの重みを変えるのではなく、影響範囲を調整するため、すべての観測が適切に考慮される。
4. **柔軟な空間構造**: 産地から遠い地点では、距離事前から大きくずれることが許容されるため、複雑な運搬パターンに適応できる。

## 12.8 実装の注意点

### 計算コスト

距離依存分散カーネルは、標準 RBF カーネルに局所標準偏差の外積を乗じるだけなので、計算コストの増加は最小限である。NNGP 近似を用いる場合、影響はほぼ無視できる。

### 数値安定性

距離スケーリング係数 $\gamma$ が大きすぎる場合、分散が過度に大きくなり数値的不安定性を招く可能性がある。実装では、局所標準偏差の上限を設定することを推奨する：

```python
sigma = np.clip(sigma_0 * (1.0 + gamma * distances), sigma_0, sigma_max)
```

### 複数産地への対応

多項ロジスティックモデルでは、各カテゴリ $k$ が異なる産地 $\mathbf{o}_k$ に対応する。カテゴリごとに異なる距離依存分散カーネルを用いることで、各産地からの距離に応じた空間構造を個別にモデル化できる。

---

**最終更新**: 2025年10月24日
**関連ドキュメント**: `nngp_multinomial_model.md`（ベースモデルの詳細導出）
