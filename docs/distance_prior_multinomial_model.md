# 距離先行情報付き NNGP 多項比率モデル

本ドキュメントでは、`docs/nngp_multinomial_model.md` で導入した多項ロジット + NNGP モデルに「産地までの距離」を明示的な先行情報として組み込み、距離に裏付けられた決定論的成分と NNGP 残差の役割分担を数式レベルで整理する。観測密度が低い地域では距離による信念を優先し、境界付近や混在領域では NNGP 残差により柔軟に補正する構成を目指す。

---

# 1. データと記法
\[
\begin{aligned}
&\mathcal{D} \subset \mathbb{R}^2:&& \text{観測領域}, \\
&s_i \in \mathcal{D},\ i=1,\ldots,n:&& \text{遺跡（観測地点）}, \\
&c_k \in \mathcal{D},\ k=1,\ldots,K:&& \text{産地（距離の基準点）}, \\
&\mathbf{y}_i = (y_{ik})_{k=1}^K,\quad y_{ik}\ge 0,\quad N_i = \sum_{k=1}^K y_{ik}:&& \text{カウントデータ}, \\
&\mathbf{W}(s) = \bigl(1, W_1(s),\ldots,W_p(s)\bigr)^\top:&& \text{共変量ベクトル}.
\end{aligned}
\]

距離関数はカテゴリごとに
\[
  d_k(s) = \mathrm{dist}(s, c_k)
  \quad
  (\text{ユークリッド距離、歩行コストなど})
\]
で定義する。必要に応じて $d_k(s)$ を $[0,1]$ に射影する（例：$d_k(s) \leftarrow d_k(s)/\max_i d_k(s)$）。

---

# 2. 距離ベースの決定論的基底項
距離から誘導される基底項を
\[
\begin{aligned}
\mu_k(s)
  &= \alpha_k - \gamma_k\, d_k(s),
  &&k = 1,\ldots,K-1, \\
\mu_K(s)
  &= 0
  &&\text{(基準カテゴリ)}.
\end{aligned}
\]
ここで
\[
\begin{aligned}
\alpha_k &= \mu_k(c_k) &&\text{(距離ゼロでのロジット)},\\
\gamma_k &> 0 &&\text{(減衰係数: 距離 1 の増分でロジットを $\gamma_k$ だけ減少)}.
\end{aligned}
\]
距離を正規化する場合は $d_k(s) \in [0,1]$ を前提にすれば、$\alpha_k$ と $\gamma_k$ の大きさの意味付けが容易になる。

---

# 3. 距離基底 + 残差による線形予測子
距離基底項に加え、共変量係数と NNGP 残差を加える。
\[
\begin{aligned}
\mathbf{b}_k(s) &= (\beta_{0k}(s),\ldots,\beta_{pk}(s))^\top, \\
\delta_k(s) &\sim \text{Residual GP}, \\
\eta_k(s)
  &= \mu_k(s) + \mathbf{W}(s)^\top \mathbf{b}_k(s) + \delta_k(s),
     &&k = 1,\ldots,K-1.
\end{aligned}
\]
基準カテゴリ $K$ のロジットはゼロ（softmax での基準）とし、全体の softmax：
\[
\begin{aligned}
\pi_k(s)
  &= \frac{\exp\{\eta_k(s)\}}
           {1 + \sum_{\ell=1}^{K-1} \exp\{\eta_\ell(s)\}},
  &&k = 1,\ldots,K-1, \\
\pi_K(s)
  &= \frac{1}
           {1 + \sum_{\ell=1}^{K-1} \exp\{\eta_\ell(s)\}}.
\end{aligned}
\]
$\eta_k(s)$ の中心を距離基底 $\mu_k(s)$ が決め、$\mathbf{W}(s)^\top \mathbf{b}_k(s)$ と $\delta_k(s)$ が局所的な偏りを補正する。

---

# 4. 観測尤度
遺跡 $i$ のカウントは
\[
\mathbf{y}_i \mid N_i, \boldsymbol{\pi}(s_i)
  \sim \operatorname{Multinomial}\bigl(N_i, \boldsymbol{\pi}(s_i)\bigr)
\]
で、尤度は
\[
\begin{aligned}
p(\mathbf{y}_i \mid N_i, \{\boldsymbol{\beta}_k\})
  &= \frac{N_i!}{\prod_{k=1}^K y_{ik}!}
     \prod_{k=1}^{K-1}
       \left(
         \frac{\exp\{\eta_k(s_i)\}}
              {1 + \sum_{\ell=1}^{K-1} \exp\{\eta_\ell(s_i)\}}
       \right)^{\!y_{ik}}
     \left(
       \frac{1}{1 + \sum_{\ell=1}^{K-1} \exp\{\eta_\ell(s_i)\}}
     \right)^{\!y_{iK}} \\
  &= c_i
     \prod_{k=1}^{K-1}
       \exp\{y_{ik}\,\eta_k(s_i)\}
       \left(1 + \sum_{\ell=1}^{K-1} \exp\{\eta_\ell(s_i)\}\right)^{-N_i},
\end{aligned}
\]
ここで $c_i = \dfrac{N_i!}{\prod_{k=1}^K y_{ik}!}$ はデータだけに依存する多項係数。

---

# 5. 距離基底に対する事前
距離基底のパラメータに弱情報事前を置く例：
\[
\begin{aligned}
\alpha_k &\sim \mathcal{N}(m_\alpha, s_\alpha^2), \\
\gamma_k &\sim \operatorname{LogNormal}(m_\gamma, s_\gamma^2)
          \quad (\text{あるいは } \gamma_k \sim \mathcal{HN}).
\end{aligned}
\]
この事前によって「距離ゼロでほぼ 1」「距離が増えると急激に減衰」等の信念を柔らかく調整できる。
距離を正規化した場合は $\alpha_k$ をロジットでの基準値（例：$\alpha_k = \log 9$ なら 90 %）、$\gamma_k$ を距離 1（=最大距離）でどれだけ減衰させたいかで設定する。

---

# 6. 残差項への NNGP 事前
距離基底で大域的構造を抑えた後、残差 $\delta_k(s)$ にのみ Vecchia 近似を適用する。
\[
\delta_k(s_i) \mid \delta_k(S_{N_i})
  \sim \mathcal{N}\bigl(a_i^\top \delta_k(S_{N_i}),\, d_i\bigr),
\quad
\delta_k(s_*) \mid \delta_k(S)
  \sim \mathcal{N}\bigl(a_*(s_*)^\top \delta_k(S),\, d_*(s_*)\bigr).
\]
これにより距離基底が説明し切れない局所的な偏りだけを NNGP で表現できる。

---

# 7. Pólya–Gamma 補助変数の導入
\[
\begin{aligned}
\log p(\mathbf{y}_i \mid \boldsymbol{\eta}(s_i))
  &= \sum_{k=1}^{K-1} y_{ik} \eta_k(s_i)
     - N_i \log\!\Bigl(1 + \sum_{\ell=1}^{K-1} \exp\{\eta_\ell(s_i)\}\Bigr).
\end{aligned}
\]
Polson et al. (2013) の恒等式
\[
\frac{(e^\psi)^a}{(1+e^\psi)^b}
  = 2^{-b} e^{\kappa \psi}
    \int_0^\infty e^{-\frac{\omega \psi^2}{2}}
                 p_{\mathrm{PG}}(\omega \mid b,0)\,d\omega,
    \quad \kappa = a - \tfrac{b}{2},
\]
を各カテゴリに適用すると、条件付き尤度は
\[
\begin{aligned}
p(\mathbf{y}_i \mid \boldsymbol{\eta}(s_i), \boldsymbol{\omega}_i)
  &\propto
     \prod_{k=1}^{K-1}
       \exp\left\{
         \left(y_{ik}-\tfrac{N_i}{2}\right)\eta_k(s_i)
         - \tfrac{1}{2}\omega_{ik}\,\eta_k(s_i)^2
       \right\}
\end{aligned}
\]
となり、$\omega_{ik} \sim \mathrm{PG}(N_i,\eta_k(s_i))$。

---

# 8. オフセットつき残差更新
距離基底をオフセットとして扱うため
\[
\psi_{ik} = \eta_k(s_i) - \mu_k(s_i)
          = \mathbf{W}(s_i)^\top \mathbf{b}_k(s_i) + \delta_k(s_i)
\]
と書き換える。平方完成から
\[
\beta_{jk}(s_i) \mid \text{rest}
  \sim \mathcal{N}\!\left(
      \frac{
        \omega_{ik} w_{ij} (\psi_{ik}^{(-j)} - C_{ik})
        + \kappa_{ik} w_{ij}
        + d_i^{-1}\mu_{ik}^{\text{prior}}
      }{
        \omega_{ik} w_{ij}^2 + d_i^{-1}
      },
      \;
      \frac{1}{\omega_{ik} w_{ij}^2 + d_i^{-1}}
    \right),
\]
\[
\begin{aligned}
\psi_{ik}^{(-j)} &= \psi_{ik} - w_{ij} \beta_{jk}(s_i), \\
C_{ik}
  &= \log\left(1 + \sum_{\ell \neq k} \exp\{\mu_\ell(s_i) + \psi_{i\ell}\}\right), \\
\mu_{ik}^{\text{prior}}
  &= a_i^\top \beta_{jk}(S_{N_i}).
\end{aligned}
\]
距離基底 $\mu_k(s_i)$ は常に既知のオフセットとして加算され、残差更新では $\psi_{ik}$ にのみランダム性が残る。

---

# 9. ギブスサンプラー
1. 距離基底 $\mu_k(s)$ を全点で前計算。
2. 残差 $\delta_k(s)$ と共変量係数をゼロで初期化（初期ロジットは $\mu_k(s)$）。
3. 各反復で
   \[
     \omega_{ik} \sim \mathrm{PG}\!\left(N_i,\mu_k(s_i) + \psi_{ik}\right).
   \]
4. 距離基底をオフセットとして残したまま $\psi_{ik}$ の平方完成 → $\beta_{jk}(s_i)$ 更新。
5. 必要なら距離パラメータ $\alpha_k,\gamma_k$ をメトロポリス等で更新。
6. MCMC サンプルから構成比 $\boldsymbol{\pi}(s)$ を取得。

距離基底を固定したまま更新できるため、従来モデルと同程度の計算コストで学習可能。

---

# 10. 予測と境界の推定
グリッド点 $u$ のロジット：
\[
\eta_k(u) = \mu_k(u) + \mathbf{W}(u)^\top \mathbf{b}_k(u) + \delta_k(u).
\]
$\mu_k(u)$ は距離に基づく滑らかなグラデーション、$\delta_k(u)$ は NNGP 条件付き分布で予測。
\[
\delta_k(u) \mid \delta_k(S) \sim \mathcal{N}\bigl(a_*(u)^\top \delta_k(S),\, d_*(u)\bigr)
\]
を用いてサンプルし、softmax を通して構成比マップと境界推定を行う。

---

# 11. まとめ
- 距離基底 $\mu_k(s)$ は「産地に近い → ロジットが大きい」を明示的に実装。
- NNGP 残差 $\delta_k(s)$ は境界付近の局所構造を補正。
- Pólya–Gamma 補助により、距離基底をオフセットとしたままガウス型のギブス更新が可能。
- 未観測域でも距離 prior に従う構成比が得られ、境界線は NNGP により自動調整される。

---

# 参考文献
- Polson, N., Scott, J., & Windle, J. (2013). *Bayesian inference for logistic models using Pólya–Gamma latent variables*. JASA.
- Linderman, S., Johnson, M., & Adams, R. (2015). *Dependent multinomial models made easy: Stick-breaking with the Pólya-Gamma augmentation*. NIPS.
- Finley, A. O., Datta, A., Banerjee, S., & Gelfand, A. E. (2019). *Nearest Neighbor Gaussian Processes for massive spatial data sets*. JASA.
