# 多項カウントデータに対する NNGP モデル

本ドキュメントでは、黒曜石産地構成比のような多カテゴリカウントデータを扱うために構築した **Nearest-Neighbor Gaussian Process (NNGP) 多項ロジットモデル** を、モデル生成過程のレベルから詳細に記述する。`docs/model_description.md` の点過程モデルと同じ筆致で、生成モデル → 尤度 → 事前分布 → ギブスサンプラーの完全条件付き分布を順を追って導出する。

---

# 1. 生成モデル

## 1.1 観測データと共変量

- 観測領域：$\mathcal{D} \subset \mathbb{R}^2$。
- 遺跡（観測地点）：$s_i \in \mathcal{D}$, $i=1,\ldots,n$。
- 産地カテゴリ：$K$ 種（例：神津島・信州・箱根・高原山・その他）。
- 観測カウント：
  \[
  \mathbf{y}_i = (y_{i1},\ldots,y_{iK})^\top, \qquad N_i = \sum_{k=1}^K y_{ik}.
  \]
- 共変量ベクトル（切片を含む）：
  \[
  \mathbf{W}(s) = \bigl(1, W_1(s),\ldots,W_p(s)\bigr)^\top.
  \]

基準カテゴリを $K$ 番目とし、$K-1$ 個のロジットを明示的にモデル化する。実装では `prepare_multinomial_dataset` がこれらの配列を生成し、`MultinomialDataset` に格納する。

## 1.2 ロジットとソフトマックス

カテゴリ $k=1,\ldots,K-1$ について、空間的に変化する係数ベクトル
\[
\boldsymbol{\beta}_k(s) = \bigl(\beta_{0k}(s),\ldots,\beta_{pk}(s)\bigr)^\top
\]
を導入し、線形予測子を
\[
\eta_k(s) = \mathbf{W}(s)^\top \boldsymbol{\beta}_k(s)
\]
と定義する。softmax 写像により、確率ベクトル
\[
\begin{aligned}
\pi_k(s) &= \frac{\exp\{\eta_k(s)\}}{1 + \sum_{\ell=1}^{K-1} \exp\{\eta_\ell(s)\}}, && k=1,\ldots,K-1,\\
\pi_K(s) &= \frac{1}{1 + \sum_{\ell=1}^{K-1} \exp\{\eta_\ell(s)\}}
\end{aligned}
\]
を得る。$\boldsymbol{\pi}(s) = (\pi_1(s),\ldots,\pi_K(s))$ は単体上の確率ベクトルである。

## 1.3 生成手順

基準カテゴリを含む多項分布に従う生成過程は以下のとおり。

1. 各遺跡 $i$ について、総出土数 $N_i$ を観測する（既知の外生データ）。
2. 共変量 $\mathbf{W}(s_i)$ を構築する。
3. 係数場 $\{\boldsymbol{\beta}_k(s)\}$ をガウス過程事前から生成する（後述）。
4. softmax より構成比 $\boldsymbol{\pi}(s_i)$ を計算。
5. 観測カウントを
   \[
   \mathbf{y}_i \mid N_i, \boldsymbol{\pi}(s_i) \sim \operatorname{Multinomial}\!\bigl(N_i, \boldsymbol{\pi}(s_i)\bigr)
   \]
   から得る。

---

# 2. 尤度の導出

## 2.1 多項ロジットの尤度

多項分布の確率質量関数は
\[
p(\mathbf{y}_i \mid N_i, \boldsymbol{\pi}(s_i))
= \frac{N_i!}{\prod_{k=1}^K y_{ik}!} \prod_{k=1}^K \pi_k(s_i)^{y_{ik}}
\]
である。softmax を代入し、ログを取ると
\[
\begin{aligned}
\log p(\mathbf{y}_i \mid \boldsymbol{\eta}(s_i))
&= \log \frac{N_i!}{\prod_{k=1}^K y_{ik}!}
  + \sum_{k=1}^{K-1} y_{ik} \eta_k(s_i)
  - N_i \log\!\left(1 + \sum_{\ell=1}^{K-1} \exp\{\eta_\ell(s_i)\}\right).
\end{aligned}
\]
第1項はデータのみの関数であるため定数として扱う。

## 2.2 全データの対数尤度

全遺跡分を総和すると
\[
\begin{aligned}
\log p(\mathbf{y} \mid \boldsymbol{\eta})
&= \sum_{i=1}^n \sum_{k=1}^{K-1} y_{ik} \eta_k(s_i)
   - \sum_{i=1}^n N_i \log\!\left(1 + \sum_{\ell=1}^{K-1} \exp\{\eta_\ell(s_i)\}\right)
   + \text{const}.
\end{aligned}
\]
softmax の分母が非線形のため、事前と組み合わせても閉じた形にはならない。次節で Pólya–Gamma 補助変数を用いてガウス化する。

---

# 3. Pólya–Gamma 補助変数

## 3.1 恒等式

Polson, Scott & Windle (2013) の恒等式を用いると、
\[
\frac{\exp\{a \psi\}}{(1 + \exp\{\psi\})^b}
= 2^{-b} \exp\!\bigl\{\kappa \psi\bigr\}
  \int_0^\infty \exp\!\left(-\frac{\omega \psi^2}{2}\right)
  p_{\text{PG}}(\omega \mid b, 0)\, d\omega,
\]
ただし $\kappa = a - b/2$、$p_{\text{PG}}$ は Pólya–Gamma 分布密度を表す。この恒等式をカテゴリごとに適用する。

## 3.2 補助変数の導入

各 $i,k$ について
\[
\omega_{ik} \mid - \sim \operatorname{PG}\bigl(N_i, \eta_k(s_i)\bigr)
\]
を導入すると、条件付き対数尤度は
\[
\begin{aligned}
\log p(\mathbf{y}_i \mid \boldsymbol{\eta}(s_i), \boldsymbol{\omega}_i)
&= \sum_{k=1}^{K-1}
   \left(
     \kappa_{ik} \eta_k(s_i)
     - \tfrac{1}{2}\,\omega_{ik}\,\eta_k(s_i)^2
   \right)
   + \text{const},
\end{aligned}
\]
となる。ここで
\[
\kappa_{ik} = y_{ik} - \tfrac{N_i}{2}.
\]
この形は $\eta_k(s_i)$ について二次形式であり、ガウス過程事前と組み合わせるとガウス分布になる。

---

# 4. ガウス過程事前

## 4.1 ベクトル表現

係数をベクトル化する。観測点集合 $S = \{s_1,\ldots,s_n\}$ に対し、
\[
\boldsymbol{\beta}_{jk,S}
  = \bigl(\beta_{jk}(s_1),\ldots,\beta_{jk}(s_n)\bigr)^\top.
\]
ゼロ平均ガウス過程を仮定すると
\[
\boldsymbol{\beta}_{jk,S} \sim \mathcal{N}(\mathbf{0}, K_\theta),
\]
ただし $(K_\theta)_{ab} = C_\theta(s_a, s_b)$ である。典型的には RBF カーネル
\[
C_\theta(s,s') = \sigma^2 \exp\!\left(-\frac{\|s - s'\|^2}{2\ell^2}\right)
\]
を用いるが、Matérn 等でもよい。

## 4.2 Vecchia / NNGP 近似

大規模データでは $K_\theta$ の逆行列を扱うのが困難なため、Vecchia 近似を適用する。Morton 順 $o$ に基づき、各点 $s_i$ の先行点から $M$ 個の近傍集合 $N_i$ を選ぶ。条件付き分布を
\[
\beta_{jk}(s_i) \mid \boldsymbol{\beta}_{jk}(S_{N_i})
  \sim \mathcal{N}\bigl(a_i^\top \boldsymbol{\beta}_{jk}(S_{N_i}),\, d_i\bigr)
\]
と近似し、$a_i$・$d_i$ をカーネル行列の部分行列から計算する。`build_nngp_factors` がこの分解を実装する。

全点をまとめると、
\[
\log p(\boldsymbol{\beta}_{jk,S})
= -\tfrac{1}{2} \sum_{i=1}^n
   \left[
     \frac{\bigl(\beta_{jk}(s_i) - a_i^\top \boldsymbol{\beta}_{jk}(S_{N_i})\bigr)^2}{d_i}
     + \log d_i
   \right]
   + \text{const}.
\]
行列表現では疎行列 $A$（行 $i$ に $a_i$ を配置）と対角行列 $D=\operatorname{diag}(d_1,\ldots,d_n)$ を用いる。

## 4.3 グリッド点への外挿

グリッド点 $u$ に対しても同様にクロス因子
\[
\beta_{jk}(u) \mid \boldsymbol{\beta}_{jk,S}
  \sim \mathcal{N}\bigl(a_*(u)^\top \boldsymbol{\beta}_{jk,S},\, d_*(u)\bigr)
\]
を前計算し、予測時に利用する。`build_cross_factors` がこれを担う。

---

# 5. 事後分布とギブスサンプラー

## 5.1 事後分布の分解

補助変数を含む完全データの事後分布は
\[
p(\boldsymbol{\beta}, \boldsymbol{\omega}
 \mid \mathbf{y}, \mathbf{W})
\propto
  p(\mathbf{y} \mid \boldsymbol{\beta}, \boldsymbol{\omega}, \mathbf{W})
  p(\boldsymbol{\omega} \mid \boldsymbol{\beta}, \mathbf{W})
  p(\boldsymbol{\beta})
\]
に比例する。PG 変数は $\boldsymbol{\beta}$ に条件付きで独立にサンプルできるため、ギブスサンプラーを構成できる。

\begin{frame}{$\omega_{ik}$ Conditional}
\[
\begin{aligned}
p(\omega_{ik} \mid \mathbf{y}, \boldsymbol{\beta})
&\propto p(\mathbf{y}_i \mid \boldsymbol{\beta}, \omega_{ik})\,
p(\omega_{ik} \mid N_i,0) \\
&\propto
\exp\!\left\{
  \left(y_{ik}-\tfrac{N_i}{2}\right)\mathbf{W}_i^\top\boldsymbol{\beta}_k
\right\}
\exp\!\left\{
  -\tfrac{\omega_{ik}}{2}(\mathbf{W}_i^\top\boldsymbol{\beta}_k)^2
\right\}
\cdot p(\omega_{ik} \mid N_i,0) \\
&\propto
\exp\!\left\{
  -\tfrac{\omega_{ik}}{2}(\mathbf{W}_i^\top\boldsymbol{\beta}_k)^2
\right\}
\cdot p(\omega_{ik} \mid N_i,0).
\end{aligned}
\]
Polson et al. (2013), Theorem 1:
\[
p(\omega \mid b, \psi)
\propto
\exp\!\left(-\tfrac{\omega}{2}\psi^2\right)p(\omega \mid b,0)
\Rightarrow
\omega \mid \psi \sim \mathrm{PG}(b, \psi).
\]
したがって
\[
\omega_{ik} \mid \mathbf{y}, \boldsymbol{\beta}
\sim \mathrm{PG}\!\left(N_i,\,\mathbf{W}_i^\top\boldsymbol{\beta}_k\right).
\]
\end{frame}

## 5.3 $\boldsymbol{\beta}$ の完全条件付き分布

カテゴリ $k$ を固定する。PG 補助を導入した条件付き対数尤度は
\[
\begin{aligned}
\log p(\mathbf{y}_{\cdot k} \mid \boldsymbol{\beta}_k, \boldsymbol{\omega}_{\cdot k})
= \sum_{i=1}^n
  \left(
    \kappa_{ik} \eta_k(s_i)
    - \tfrac{1}{2}\,\omega_{ik}\,\eta_k(s_i)^2
  \right)
\end{aligned}
\]
であった。設計行列を $W = [w_{ij}]$、係数場を行列 $B_k = [\beta_{jk}(s_i)]_{j=0,\ldots,p}^{i=1,\ldots,n}$ とすると
\[
\eta_k = W \odot B_k^\top \mathbf{1}_{p+1},
\]
ここで $\odot$ は要素積、$\mathbf{1}_{p+1}$ は単位ベクトルである。計算を簡潔にするため、共変量ごとに更新する。成分 $j$ のみを取り出すと、
\[
\eta_k(s_i) = w_{ij} \beta_{jk}(s_i) + \eta_k^{(-j)}(s_i),
\]
を得る。これを条件付き対数尤度に代入し、$\beta_{jk}(s_i)$ について平方完成を行うと
\[
\beta_{jk}(s_i) \mid \text{他} \sim
\mathcal{N}\!\left(
  \frac{ \omega_{ik} w_{ij} \left(\eta_k^{(-j)}(s_i) - C_{ik}\right)
         + \kappa_{ik} w_{ij}
         + \mu_{ik}^{\text{prior}} / d_i }{\omega_{ik} w_{ij}^2 + 1/d_i},
  \;
  \frac{1}{\omega_{ik} w_{ij}^2 + 1/d_i}
\right),
\]
となる。ここで
\[
\begin{aligned}
\mu_{ik}^{\text{prior}} &= a_i^\top \boldsymbol{\beta}_{jk}(S_{N_i}),\\
C_{ik} &= \log\!\left(1 + \sum_{\ell \ne k} \exp\{\eta_\ell(s_i)\}\right),
\end{aligned}
\]
である。`update_beta_category` はこの一変量正規を Morton 順に沿ってサンプルし、$\eta_k$ を逐次更新する。

## 5.4 アルゴリズムのまとめ

1. 初期化：$\boldsymbol{\beta}_k$ をゼロ（または総計に基づくロジット）で初期化し、$\eta_k = W \beta_k$ を計算。
2. 反復：
   - $\omega_{ik} \sim \operatorname{PG}(N_i, \eta_k(s_i))$ を全 $i,k$ でサンプル。
   - 各カテゴリ $k$ について、成分 $j=0,\ldots,p$ を順に更新。
   - 必要に応じてハイパーパラメータ（長さ尺度・分散）をメトロポリス–ヘイスティング等で更新（現状は固定）。
3. バーンイン後、所定の間隔で $\boldsymbol{\beta}$ サンプルを保存。

`run_mcmc` はこの流れを実装し、`MultinomialNNGPResults` が結果を保持する。

---

# 6. 予測と後モデルチェック

## 6.1 遺跡位置での構成比

保存した $\boldsymbol{\beta}$ サンプルごとに
\[
\eta_k(s_i) = \mathbf{W}(s_i)^\top \boldsymbol{\beta}_k(s_i)
\]
を計算し、softmax を通すことで事後分布を得る。`posterior_site_mean` はサンプル平均を返す。

## 6.2 グリッド点での予測

グリッド点 $u$ について、`FactorCache` に事前計算されたクロス因子 $(a_*(u), d_*(u))$ を用い、
\[
\boldsymbol{\beta}_k(u) \mid \boldsymbol{\beta}_k(S)
  \sim \mathcal{N}\bigl(a_*(u)^\top \boldsymbol{\beta}_k(S),\, d_*(u)\bigr)
\]
を計算する。`posterior_grid_mean` はこれを softmax に通した平均を返し、`sample_conditional=True` とすると Vecchia 分散を用いたノイズを加えて不確実性を反映する。

---

# 7. 実装上の注意

- **共変量の標準化**：`ObsidianDataPreprocessor(scale_variables=True)` を用いると、遺跡側の平均・標準偏差で共変量を中心化・スケーリングできる。
- **ゼロカウント遺跡の除外**：`prepare_multinomial_dataset(..., drop_zero_total_sites=True)` により、$N_i = 0$ の遺跡を学習から除外可能。
- **グリッド間引き**：`grid_subsample_ratio` を設定すると、グリッド点を等間隔に間引き NNGP 前処理時間を削減できる。
- **カーネルの列別設定**：`kernel_lengthscale_by_feature` や `kernel_variance_by_feature` で共変量ごとにハイパーパラメータを設定できるが、FactorCache 初期化時間は列数に比例して増える。
- **FactorCache の再利用**：同一データセットで複数回サンプラーを走らせる場合、FactorCache を外部で保持すると初期化コストを節約できる。

---

# 8. 点過程モデルとの比較

`docs/model_description.md` の非斉次ポアソン過程モデルでは
\[
\lambda_t(s) = \lambda_t^* \,\Phi\bigl(\mathbf{W}_t(s)^\top \boldsymbol{\beta}_t(s)\bigr)
\]
を推定し、点の強度を扱った。多項モデルは
\[
\pi_k(s) = \frac{\exp\{\mathbf{W}(s)^\top \boldsymbol{\beta}_k(s)\}}
                 {1 + \sum_{\ell=1}^{K-1} \exp\{\mathbf{W}(s)^\top \boldsymbol{\beta}_\ell(s)\}}
\]
を推定し、構成比に焦点を当てる。両モデルとも Pólya–Gamma 補助変数と Vecchia/NNGP 近似を共有しており、実装上の要素（FactorCache、KDTree 近傍探索など）も共通である。

---

# 9. 距離先行情報との接続

観測が乏しい地域で地理的常識を反映させるため、距離に基づく基底項を導入する拡張を `docs/distance_prior_multinomial_model.md` に整理した。線形予測子を
\[
\eta_k(s) = \underbrace{\alpha_k - \gamma_k d_k(s)}_{\text{距離先行}}
            + \mathbf{W}(s)^\top \tilde{\boldsymbol{\beta}}_k(s)
            + \delta_k(s)
\]
と分解し、距離先行を既知オフセットとして扱えば、PG 補助および Vecchia 近似の実装をそのまま流用できる。

---

# 参考文献

- Finley, A. O., Datta, A., Banerjee, S., & Gelfand, A. E. (2019). *Nearest Neighbor Gaussian Processes for massive spatial data sets*. Journal of the American Statistical Association.
- Linderman, S., Johnson, M., & Adams, R. (2015). *Dependent multinomial models made easy: Stick-breaking with the Pólya-Gamma augmentation*. Advances in Neural Information Processing Systems.
- Polson, N., Scott, J., & Windle, J. (2013). *Bayesian inference for logistic models using Pólya–Gamma latent variables*. Journal of the American Statistical Association.
