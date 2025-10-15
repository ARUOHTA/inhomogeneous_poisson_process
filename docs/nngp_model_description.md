https://onlinelibrary.wiley.com/doi/10.1111/jtsa.12457

# GCPの定式化

### 設定

**観測データ**：$\mathcal{S} = \set{\mathcal{S}_1, ..., \mathcal{S}_T}$*、*ただし $\mathcal{S}_t = \set{s_{t,1}, ..., s_{t,n_t}}$

**パラメータ**：

- $\boldsymbol{\beta} = {\boldsymbol{\beta}_1(s), ..., \boldsymbol{\beta}_T}(s)$：回帰係数　$*\boldsymbol{\beta}_t(s) = (\beta_{t,0}(s), \beta_{t,1}(s), ..., \beta_{t,p}(s))^\top*$
- $\boldsymbol{\lambda}^{*} = {\lambda_1^*, ..., \lambda_T^*}$：強度の上限

**共変量**：

- $W_t(s) = (1, W_{t,1}(s), ..., W_{t,p}(s))^\top$：位置s、時刻tでのp個の共変量 + 切片

### 強度関数の定義

時刻t、位置sでの強度関数を

$$
\lambda_t(s) = \lambda_t^* \cdot \Phi[f(W_t(s), \boldsymbol{\beta}_t(s))]
$$

とする。ここで、$\Phi$ は標準正規分布の累積分布関数で線型予測子 $f$ の値を $[0, 1]$ に制限する。また、

$$
\begin{aligned}
f(W_t(s), \boldsymbol{\beta}_t(s)) &= W_t(s)^\top\boldsymbol{\beta}_t(s) \\
&= \beta_{t,0}(s) + \sum_{j=1}^p \beta_{t,j}(s)W_{t,j}(s) \\

\end{aligned}
$$

である。

### 条件付き尤度

$\boldsymbol{\beta}, \boldsymbol{\lambda}^*$ が与えられた下での $\mathcal{S}$ の尤度：

$$
\begin{aligned}
\mathcal{L}(\mathcal{S}|\boldsymbol{\beta}, \boldsymbol{\lambda}^*) &= \prod_{t=1}^T \mathcal{L}_t(\mathcal{S}_t|\boldsymbol{\beta}_t, \lambda_t^*)
\end{aligned}
$$

各時刻tについて考える。

任意の部分集合$A \in \mathcal{D}$について、強度関数の積分を

$$
\begin{aligned}

\Lambda_t(A)&=\int_A \lambda_t(s) ds \\

&=\int_A \lambda_t^* \cdot \Phi[f(W_t(s), \boldsymbol{\beta}_t(s))]ds \\
\end{aligned}
$$

とすると、観測データ $\mathcal{S}_t, n_t$ が強度 $\lambda_t$ の非斉次ポアソン過程に従うことは

$$
\mathcal{S_t}, n_t \sim \text{IPP}(\lambda_t(\cdot)) \quad \Leftrightarrow \quad n_t \sim \operatorname{Poisson}\left(\Lambda_t(\mathcal{D})\right), \quad s_{t,i} \stackrel{\mathrm{iid}}{\sim} \frac{\lambda_t\left(s\right)}{\Lambda_t(\mathcal{D})}, \quad i=1, \ldots, n_t
$$

と定義される。つまり、個数$n_t$はポアソン分布から生成され、位置$s_{t, i}$は$\frac{\lambda_t\left(s\right)}{\Lambda_t(\mathcal{D})}$を確率密度関数とする分布から独立同分布として生成される。この定義からモデルの尤度を導くことができ、

$$
\begin{aligned}
\mathcal{L}_t(\mathcal{S}_t, n_t|\boldsymbol{\beta}_t, \lambda_t^*)

&= P\left(n_t \mid \boldsymbol{\beta}_t, \lambda_t^*\right) \cdot \prod_{i=1}^{n_t} P\left(s_{t,i}  \mid \boldsymbol{\beta}_t, \lambda_t^*, n_t\right) \\

&= \frac{\left[\Lambda_t(\mathcal{D})\right]^{n_t} e^{-\Lambda_t(\mathcal{D})}}{n_t!} \cdot \prod_{i=1}^{n_t} \frac{\lambda_t(s_{t,i} )}{\Lambda_t(\mathcal{D})}\\

&= \frac{e^{-\Lambda_t(\mathcal{D})}}{n_t!} \cdot \prod_{i=1}^{n_t} \lambda_t(s_{t,i} )
\\
&= \frac{1}{n_t !} \exp\left(-\int_D \lambda_t(s)ds\right) \prod_{i=1}^{n_t} \lambda_t(s_{t,i}) \\
&= \frac{(\lambda_t^*)^{n_t}}{n_t !} \exp\left(-\int_D \lambda_t^* \Phi[f(W_t(s), \boldsymbol{\beta}_t(s))]ds\right) \prod_{i=1}^{n_t} \Phi[f(W_t(s_{t,i}), \boldsymbol{\beta}_t(s_{t,i}))] \\
&= \frac{(\lambda_t^*)^{n_t}}{n_t !} \exp\left(-\lambda_t^* \int_D \Phi[W_t(s)^\top\boldsymbol{\beta}_t(s) ]ds\right) \prod_{i=1}^{n_t} \Phi[W_t(s_{t,i})^\top\boldsymbol{\beta}_t(s_{t, i})]
\end{aligned}
$$

となる。

## 潜在変数 $\mathcal{U}$ によるData Augmentation

### 1. 拡張データの同時分布

潜在データ $\mathcal{U}_t, m_t$ は、端的に

$$
\mathcal{U}_t, m_t \sim \text{IPP}(\lambda_t^*-\lambda_t(\cdot))
$$

として定義される。よって尤度は

$$
\begin{aligned}
\mathcal{L}_t(\mathcal{U}_t|\boldsymbol{\beta}_t, \lambda_t^*) &= \frac{1}{m_t !} \exp\left(-\int_D \lambda_t^* - \lambda_t(u)du\right) \prod_{i=1}^{m_t} (\lambda_t^* - \lambda_t(u_{t,i})) \\
&= \frac{1}{m_t !} \exp\left(-\int_D \lambda_t^*-\lambda_t^* \Phi[f(W_t(u), \boldsymbol{\beta}_t(u))]du\right) \prod_{i=1}^{m_t} \lambda_t^* -\lambda_t^* \Phi[f(W_t(u_{t,i}), \boldsymbol{\beta}_t(u_{t,i}))] \\
&= \frac{(\lambda_t^*)^{m_t}}{m_t !}\exp\left(-\lambda_t^* \int_D 1-\Phi[W_t(u)^\top\boldsymbol{\beta}_t(u)]du\right) \prod_{i=1}^{m_t}  \left(1- \Phi[W_t(u_{t,i})^\top\boldsymbol{\beta}_t(u_{t,i})]\right)
\end{aligned}
$$

よって、潜在変数 $\mathcal{U}$ も含めた場合の尤度は、元々の強度$\lambda_t$からサンプルされたSと、その反対の強度$\lambda_t^* - \lambda_t$からサンプルされたUが打ち消し合うことで、積分が消去される：

$$

\begin{aligned}
&\mathcal{L}(\mathcal{S}_{aug}| \boldsymbol{\beta}, \boldsymbol{\lambda}^*) \\
&=\mathcal{L}(\mathcal{S}, \mathcal{U} | \boldsymbol{\beta}, \boldsymbol{\lambda}^*) \\
&= \prod_{t=1}^T \mathcal{L}_t(\mathcal{S}_t|\boldsymbol{\beta}_t, \lambda_t^*) \cdot \mathcal{L}_t(\mathcal{U}_t|\boldsymbol{\beta}_t, \lambda_t^*)
\\
&= \prod_{t=1}^T  \left[\frac{(\lambda_t^*)^{n_t}}{n_t !} \exp\left(-\lambda_t^* \int_D \Phi[W_t(s)^\top\boldsymbol{\beta}_t(s)]ds\right) \prod_{i=1}^{n_t} \Phi[W_t(s_{t,i})^\top\boldsymbol{\beta}_t(s_{t,i})]
\times \frac{(\lambda_t^*)^{m_t}}{m_t !}\exp\left(-\lambda_t^* \int_D 1-\Phi[W_t(u)^\top\boldsymbol{\beta}_t(u)]du\right) \prod_{i=1}^{m_t}  \left(1- \Phi[W_t(u_{t,i})^\top\boldsymbol{\beta}_t(u_{t,i})]\right)\right]
\\
&= \prod_{t=1}^T  \left[\frac{(\lambda_t^*)^{n_t + m_t}}{n_t!m_t!} \exp\left(-\lambda_t^* \int_D \Phi[W_t(s)^\top\boldsymbol{\beta}_t(s)] + (1 - \Phi[W_t(s)^\top\boldsymbol{\beta}_t(s)) ds\right) \prod_{i=1}^{n_t} \Phi[W_t(s_{t,i})^\top\boldsymbol{\beta}_t(s_{t,i})]
\prod_{i=1}^{m_t}  \left(1- \Phi[W_t(u_{t,i})^\top\boldsymbol{\beta}_t(u_{t,i})]\right)\right]
\\
&= \prod_{t=1}^T  \left[\frac{(\lambda_t^*)^{n_t + m_t}}{n_t!m_t!} \exp\left(-\lambda_t^* \int_D  du\right) \prod_{i=1}^{n_t} \Phi[W_t(s_{t,i})^\top\boldsymbol{\beta}_t(s_{t,i})]
\prod_{i=1}^{m_t}  \left(1- \Phi[W_t(u_{t,i})^\top\boldsymbol{\beta}_t(u_{t,i})]\right)\right]
\\
&= \prod_{t=1}^T  \left[\frac{(\lambda_t^*)^{n_t + m_t}}{n_t!m_t!} \exp\left(-\lambda_t^* |\mathcal{D}|\right) \prod_{i=1}^{n_t} \Phi[W_t(s_{t,i})^\top\boldsymbol{\beta}_t(s_{t,i})]
\prod_{i=1}^{m_t}\Phi[-W_t(u_{t,i})^\top\boldsymbol{\beta}_t(u_{t,i})]\right]
\end{aligned}
$$

最後の行は、累積分布関数の性質：$1 - \Phi[x] = \Phi[-x]$ を使った。

# 事前分布

推定が必要な変数は

- $\boldsymbol{\beta} = ({\boldsymbol{\beta}_1(s), ..., \boldsymbol{\beta}_T}(s))$
- $\boldsymbol{\lambda}* = (\lambda_1^*, \ldots, \lambda_T^*)$
- 潜在変数 $\mathcal{U}$

の3つになる。

## $\boldsymbol{\beta}$ の事前分布：Gaussian Process

隣り合う時刻同士のガウス過程は $AR(1)$ 過程として記述される。以下では、共変量については独立に考える。

$$

\boldsymbol{\beta}_{t}(s) = G\boldsymbol{\beta}_{t-1}(s) + \boldsymbol{\eta}_t(s), \quad \boldsymbol{\eta}_t(s) \sim \mathcal{GP}(0, C_{\boldsymbol{\theta}_t}), \quad \text{for } t = 2, ..., T

$$

$$

\boldsymbol{\beta}_1(s) = \boldsymbol{\eta}_1(s), \quad \boldsymbol{\eta}_1(s) \sim \mathcal{GP}(0, C_{\boldsymbol{\theta}_1})

$$

ここで：

- $\boldsymbol{\eta}_t$*：共分散 $C_{\boldsymbol{\theta}_t}$* を持つ独立なガウス過程
    - $\boldsymbol{\theta}_2 = \cdots = \boldsymbol{\theta}_T = \boldsymbol{\theta}$ と設定
    - $\boldsymbol{\theta}_1$ *は $C_{\boldsymbol{\theta}_1}$* がより大きな分散と強い空間依存性を持つように設定
- $G$：時間遷移行列（例：自己回帰係数、単位行列）
    - $G = I$（単位行列）とすることが多い

関数 $\boldsymbol{\beta}_t(s)$ の非可算無限個の値のうち、尤度に関連するのは $\mathcal{S}$ と $\mathcal{U}$ に含まれる点における有限個の値のみである。それらをそれぞれ：

$$

\boldsymbol{\beta}_{\mathcal{S},t} = (\boldsymbol{\beta}_t(s_{t, 1}), \ldots, \boldsymbol{\beta}_t(s_{t, n_t}))^\top, \quad t = 1, \ldots, T

$$

$$

\boldsymbol{\beta}_{\mathcal{U},t} = (\boldsymbol{\beta}_t(u_{t, 1}), \ldots, \boldsymbol{\beta}_t(u_{t, m_t}))^\top, \quad t = 1, \ldots, T

$$

とする。以下では簡単のため $G = I$ とする。

### 時系列全体の同時分布

時刻 t における全ての点での値を：

$$

\boldsymbol{\beta}_{t} = \begin{pmatrix} \boldsymbol{\beta}_{\mathcal{S},t} \ \boldsymbol{\beta}_{\mathcal{U},t}\end{pmatrix}

$$

とする。$AR(1)$ 過程の再帰的表現から：

$$

\begin{align}
\boldsymbol{\beta}_1 &= \boldsymbol{\eta}_1 \\
\boldsymbol{\beta}_2 &= \boldsymbol{\beta}_1 + \boldsymbol{\eta}_2 = \boldsymbol{\eta}_1 + \boldsymbol{\eta}_2 \\
\boldsymbol{\beta}_3 &= \boldsymbol{\beta}_2 + \boldsymbol{\eta}_3 = \boldsymbol{\eta}_1 + \boldsymbol{\eta}_2 + \boldsymbol{\eta}_3 \\
&\vdots \\
\boldsymbol{\beta}_t &= \sum_{i=1}^t \boldsymbol{\eta}_i
\end{align}

$$

となる。同時分布の導出のためには、任意の2つの時刻 $t, t' （t \leq t'）$と位置 $s, s'$ について、共分散を導出する。

$$

\begin{aligned}
\text{Cov}[\boldsymbol{\beta}_t(s), \boldsymbol{\beta}_{t'}(s')] &= \text{Cov}\left[\sum_{i=1}^t \boldsymbol{\eta}_i(s), \sum_{j=1}^{t'} \boldsymbol{\eta}_j(s')\right]
\end{aligned}

$$

ここで、$\boldsymbol{\eta}_i$ は異なる時刻間で独立なガウス過程である。したがって、$i \neq j$ のとき：

$$
\text{Cov}[\boldsymbol{\eta}_i(s), \boldsymbol{\eta}_j(s')] = 0

$$

この独立性により、共分散の和は $i = j$ となる項のみが残る：

$$

\begin{aligned}
\text{Cov}[\boldsymbol{\beta}_t(s), \boldsymbol{\beta}_{t'}(s')]
&= \sum_{i=1}^t \sum_{j=1}^{t'} \text{Cov}[\boldsymbol{\eta}_i(s), \boldsymbol{\eta}_j(s')]
\\
&= \sum_{i=1}^{\min(t,t')} \text{Cov}[\boldsymbol{\eta}_i(s), \boldsymbol{\eta}_i(s')]
\\
&= \sum_{i=1}^{t} \text{Cov}[\boldsymbol{\eta}_i(s), \boldsymbol{\eta}_i(s')]
\end{aligned}

$$

次に、各 $\boldsymbol{\eta}_i$ の共分散関数を考慮する：

- $\boldsymbol{\eta}_1 \sim \mathcal{GP}(0, C_{\boldsymbol{\theta}_1})$ より、$\text{Cov}[\boldsymbol{\eta}_1(s), \boldsymbol{\eta}_1(s')] = C_{\boldsymbol{\theta}_1}(s, s')$
- $\boldsymbol{\eta}_i \sim \mathcal{GP}(0, C_{\boldsymbol{\theta}}) \ \forall i \geq 2$ より、$\text{Cov}[\boldsymbol{\eta}_i(s), \boldsymbol{\eta}_i(s')] = C_{\boldsymbol{\theta}}(s, s')$

したがって：

$$

\begin{aligned}
\text{Cov}[\boldsymbol{\beta}_t(s), \boldsymbol{\beta}_{t'}(s')] &= \sum_{i=1}^{t} \text{Cov}[\boldsymbol{\eta}_i(s), \boldsymbol{\eta}_i(s')] \\
&= C_{\boldsymbol{\theta}_1}(s, s') + \sum_{i=2}^{t} C_{\boldsymbol{\theta}}(s, s') \\
&= C_{\boldsymbol{\theta}_1}(s, s') + (t-1) C_{\boldsymbol{\theta}}(s, s')
\end{aligned}

$$

と書くことができる。

### 条件付き分布の導出

MCMCサンプリングにおいて、新しい潜在点 $u_{t,j}$ における $\boldsymbol{\beta}_t(u_{t,j})$ の値を既知の値を条件としてサンプルする必要がある。

既知の値を：

- $\boldsymbol{\beta}_{\mathcal{S}} = {\boldsymbol{\beta}_{\mathcal{S},1}, \ldots, \boldsymbol{\beta}_{\mathcal{S},T}}$：全時刻の観測点での値
- $\boldsymbol{\beta}_{\mathcal{U} \setminus u{t,j}}$：$u_{t,j}$ を除く全ての潜在点での値

とすると、多変量正規分布の条件付き分布の性質から：

$$

\boldsymbol{\beta}_t(u_{t,j}) \mid \boldsymbol{\beta}_{\mathcal{S}}, \boldsymbol{\beta}_{\mathcal{U} \setminus u_{t,j}} \sim \mathcal{N}(\mu_{u_{t,j}|rest}, \Sigma_{u_{t,j}|rest})

$$

となる。ここで、条件付き平均と分散は：

$$

\begin{aligned}
\mu_{u_{t,j}|rest} &= K_{u_{t,j},rest} K_{rest,rest}^{-1} \boldsymbol{\beta}_{rest} \\
\Sigma_{u_{t,j}|rest} &= K_{u_{t,j},u_{t,j}} - K_{u_{t,j},rest} K_{rest,rest}^{-1} K_{rest,u_{t,j}}
\end{aligned}

$$

となる（参考：https://zenn.dev/m4thphobia/articles/9ec50e8f0d513d）。ただし、$\boldsymbol{\beta}_{rest}$ は既知の全ての値を並べたベクトル、$K$ は対応する共分散行列のブロックを表す。

導出した共分散構造$\text{Cov}[\boldsymbol{\beta}_t(s), \boldsymbol{\beta}_{t'}(s')] = C_{\boldsymbol{\theta}_1}(s, s') + (t-1)C_{\boldsymbol{\theta}}(s, s')$ を用いて、$K_{u_{t,j},rest}$ の各要素を計算すれば良い。

また、$u_{t,j}$ 自身の分散は：

$$
K_{u_{t,j},u_{t,j}} = \text{Var}[\boldsymbol{\beta}_t(u_{t,j})] = C_{\boldsymbol{\theta}_1}(u_{t,j}, u_{t,j}) + (t-1)C_{\boldsymbol{\theta}}(u_{t,j}, u_{t,j})

$$

となる。

これらの計算によって、理論上は任意の観測点が得られたもとでの新たな点での値をサンプルすることが可能になった。

ただ、計算量が点数の３乗に比例するため、大規模なデータになると、この共分散行列の計算の効率化は必須になる。

以下で、時間方向の近似と、空間方向の近似の2つの手法を導入する。

### 時間方向の近似

逐次フィルタリングによって、完全条件付き分布の時間方向での近似を行う

### 近似の内容

完全条件付き分布では、時刻 t の点 u_{t,j} をサンプリングする際に：

$$

\boldsymbol{\beta}t(u{t,j}) \mid {\text{全時刻の既知の値}} \sim \mathcal{N}(\mu_{full}, \sigma^2_{full})

$$

という分布を用いる。これに対して、逐次フィルタリングでは：

$$
\boldsymbol{\beta}_t(u_{t,j}) = \boldsymbol{\beta}_{t-1}(u_{t,j}) + \boldsymbol{\eta}_t(u_{t,j})
$$

$$
\boldsymbol{\eta}_t(u_{t,j}) \mid \text{同時刻の観測値のみ} \sim \mathcal{N}(\mu_{\text{filter}}, \sigma^2_{\text{filter}})
$$

という分解を行い、**innovation $\boldsymbol{\eta}_t$ の条件付き分布のみを考慮**する。

### 何を無視しているか

逐次フィルタリングが無視している情報は主に2つある：

1. **未来の情報の無視**

    完全条件付き分布では、未来の時刻 $\tau > t$ の観測も条件に含まれる：

    $\text{Cov}[\boldsymbol{\beta}t(u{t,j}), \boldsymbol{\beta}\tau(s)] = C{\boldsymbol{\theta}1}(u{t,j}, s) + (t-1)C_{\boldsymbol{\theta}}(u_{t,j}, s)$

    この共分散により、未来の観測が現在の推定に影響を与える（スムージング効果）。逐次フィルタリングではこれを完全に無視する。

2. **過去の情報の部分的利用**

    逐次フィルタリングでは、過去の情報は \boldsymbol{\beta}_{t-1} を通じてのみ伝播する。つまり：

    - 時刻 t-1 での推定値は固定値として扱われる
    - 時刻 t-2 以前の観測との直接的な共分散関係は考慮されない

### 近似の正当化

この近似が実用的に有効な理由：

1. **マルコフ性の仮定**：$\boldsymbol{\beta}t$ が$*\boldsymbol{\beta}{t-1}*$のみに依存するAR(1)構造を持つため、過去の情報の多くは$\boldsymbol{\beta}_{t-1}$に集約されている
2. **計算効率**：
    - 完全条件付き：O(N^3) の計算量（N は全時空間点数）
    - 逐次フィルタリング：O(n^3) の計算量（n は各時刻の点数）
3. **オンライン推定**：データが逐次的に得られる場合、未来の情報を待たずに推定できる

### 近似の影響

逐次フィルタリングと完全条件付き分布の主な違い：

- **不確実性の過大評価**：未来の情報を使わないため、推定の不確実性が大きくなる
- **時間遅れ**：観測の影響が過去に遡及しない（バックワードパスがない）
- **局所的な推定**：各時刻で独立に innovation を推定するため、時系列全体の一貫性が弱い

### 空間方向の近似：Nearest Neighbor Gaussian Process

### Nearest Neighbor Gaussian Processes（最近傍ガウス過程）

https://andrewcharlesjones.github.io/journal/nngp.html

### 基本概念

Datta et al. (2016a)によるNNGPは、空間的近傍の観測データの影響のみを考えることで計算量の大幅な削減を考える近似手法。

**同時密度の分解**：

$$
\pi(z(\mathcal{S})) = \pi(z(s_1))\pi(z(s_2)|z(s_1)) \cdots \pi(z(s_i)|z_{<i}) \cdots \pi(z(s_n)|z_{<n})
$$

ここで $z_{<i} = {z(s_1), ..., z(s_{i-1})}$

**NNGP近似**：

$$
\pi(z(\mathcal{S})) \approx \pi(z(s_1))\pi(z(s_2)|z(s_1)) \cdots \pi(z(s_i)|z_{N_i}) \cdots \pi(z(s_n)|z_{N_n}) = \tilde{\pi}(z(\mathcal{S}))
$$

ここで：

- $N_i$：$s_i$の近傍のインデックス集合
- $z_{N_i} \subseteq z_{<i}$（$z_{<i}$の部分集合）

**条件付き分布**：

$z(s_i)|z_{N_i} \sim \mathcal{N}(\mu_i, \sigma_i^2)$

ここで：

$$
\mu_i = C_{\boldsymbol{\theta}}(s_i, N_i)C_{\boldsymbol{\theta}}(N_i, N_i')^{-1}z(N_i)
$$

$$
\sigma_i^2 = C_{\boldsymbol{\theta}}(s_i, s_i) - C_{\boldsymbol{\theta}}(s_i, N_i)C_{\boldsymbol{\theta}}(N_i, N_i')^{-1}C_{\boldsymbol{\theta}}(N_i, s_i)
$$

### NNGPの行列表現

以下のように書ける：

$$
z(s_1) = \eta_1
$$

$$
z(s_i) = a_{i,1}z(s_1) + a_{i,2}z(s_2) + \cdots + a_{i,i-1}z(s_{i-1}) + \eta_i, \quad i = 2, ..., n
$$

行列形式で：

$$
z(\mathcal{S}) = Az(\mathcal{S}) + \boldsymbol{\eta}
$$

ここで：

- $A$：$n \times n$の厳密下三角行列（$a_{i,j} = 0$  when  $j \geq i$）
- $\boldsymbol{\eta} \sim \mathcal{N}(0, D)$：$D$は対角行列で $d_{i,i} = \text{Var}[z(s_i)|z_{<i}]$

**共分散行列の近似**：

$$
C = (I - A)^{-1}D(I - A)^{-T}
$$

近似共分散行列：

$$
\tilde{C} = (I - \tilde{A})^{-1}\tilde{D}(I - \tilde{A})^{-T}
$$

ここで：

- $\tilde{A}$は疎行列（$a_{i,j} \neq 0$ only when $j \in N_i$）
- $\tilde{d}_{i,i} = \text{Var}[z_i|z_{N_i}]$

### 計算複雑性

- **NNGP**：$O(nM^3)$ flops（Mは近傍の数、通常小さい値に固定）
- **通常のGP**：$O(n^3)$ flops
- **INLA**：$O(n^{3/2})$ flops、$O(n\log(n))$のメモリ

M = 10または20程度で、数百万の位置に対するGP実現を十分に近似できることが示されている。

**具体例**

Dynamic GPにおいて、t=2のみで観測値（赤丸）を得た際の各時刻での事後分布からの1つのサンプル例。（過去➡️未来方向の、時間窓1の隣しか見ていない時間近似をした。）

t-0, 1においては情報がないので事前分布からサンプルしているだけだが、t=2では観測データによって更新されており、t=3, 4では、AR構造によってt=2の観測データの情報が伝播している。

![image.png](attachment:69b2a812-a98d-437a-9c9a-d268927cd18c:image.png)

2次元での同様の例

![image.png](attachment:16e25b49-8f0e-4ca0-8320-5da5f0cfb7b1:image.png)

**近似手法の適用**

以下のように、時間方向、あるいは空間方向に近似を行うことで、共分散行列のうち計算が必要な要素を大幅に削減できる。

![image.png](attachment:95456d3c-d0d3-4621-92b4-43b9b262dae6:image.png)

これら4つの近似手法ごとの推定結果の比較。

```python
Grid: 40 locations × 12 times = 480 total points
Observations: 12 points
```

```python
Gibbs sampling: 100%|██████████| 30/30 [02:12<00:00,  4.42s/it]
Full Gibbs          : 132.56s, mean=  0.056, std= 0.832
Gibbs sampling: 100%|██████████| 30/30 [00:08<00:00,  3.60it/s]
Time window (w=1)   :   8.34s, mean=  0.009, std= 0.859
Gibbs sampling: 100%|██████████| 30/30 [00:02<00:00, 14.38it/s]
NNGP (k=5)          :   2.09s, mean=  0.156, std= 0.785
Gibbs sampling: 100%|██████████| 30/30 [00:05<00:00,  5.42it/s]
Combined (w=1, k=3) :   5.54s, mean=  0.072, std= 0.889
```

近似、特にNNGPを組み込むと相当な高速化が見込まれる。

また、推定結果も以下のようになった↓

![Cursor 2025-08-07 19.25.17.png](attachment:f010cc2a-1835-4bf3-86de-a114d0424187:Cursor_2025-08-07_19.25.17.png)

今後の実装では、NNGP (近傍数は調整する)＋時間窓w=1で計算を行う。

## $\boldsymbol{\lambda}^*$の事前分布：Gamma分布

時系列的に独立にガンマ分布に従うとしても良いし、依存させても良い。ここでは簡単のため独立とする：

$$
p(\boldsymbol{\lambda^*}) = \prod_{t=0}^T(\lambda_t^*)^{m_0 - 1} e^{- r_0 \lambda_t^*}
$$

（あるいは以下のようにしても良い）

時間依存モデル：

$$
p(\boldsymbol{\lambda}^*) = p(\lambda_1^*) \prod_{t=2}^T p(\lambda_t^*|\lambda_{t-1}^*)
$$

ここで：

- $p(\lambda_1^*) = \text{Gamma}(\lambda_1^*|a_1, b_1)$
- $p(\lambda_t^*|\lambda_{t-1}^*) = w^{-1}\lambda_{t-1}^{*-1} \text{Beta}(\lambda_t^*/\lambda_{t-1}^*|wa_t, (1-w)a_t)$

- memo

    事後推論の後に、空間内の任意の点の有限集合から強度関数の値をサンプルする必要がある。

    4.3

    IFの上限値の識別可能性について。IFの事前分布の選択がとても重要。

    オンライン補足資料のSection１。

    MCMCの高速化について

    - *Lower dimension approximations Banerjee et al. (2013)*
    - GPからのサンプルをHMCなどでやる
    - GPからのサンプルをHMCなどでやる

    Skew Nからのサンプリング：4.1

    潜在変数について

    - Homogeneous PPからサンプリングした際のremainingをMとして定義することもできる（そちらの方がいい？）→結局ほぼ等価らしい（2023のcorrigendum）

    thetaからのサンプルは、できるなら直接、そうでないならMH。thetaをブロック化してサンプルできるならそうする。adaptive MH（Roberts and Rosenthal (2009)）を使って適応的にやるのもあり。


# 完全条件付き分布の導出

ベイズの定理より、パラメータと潜在変数の事後分布は

$$
\begin{aligned}
&p(\mathcal{U}, \boldsymbol{\lambda}^*, \boldsymbol{\beta}, \boldsymbol{\theta} \mid \mathcal{S}) \\
&\propto \mathcal{L}(\mathcal{S}, \mathcal{U}| \boldsymbol{\lambda}^*, \boldsymbol{\beta}) \times p(\boldsymbol{\beta} | \boldsymbol{\theta}) \times p(\boldsymbol{\lambda}^*) \times p(\boldsymbol{\theta})
\\
&= \prod_{t=1}^T  \left[\frac{(\lambda_t^*)^{n_t + m_t}}{n_t!m_t!} \exp\left(-\lambda_t^* |\mathcal{D}|\right) \prod_{i=1}^{n_t} \Phi[W_t(s_{t,i})^\top\boldsymbol{\beta}_t(s_{t,i})]
\prod_{i=1}^{m_t}\Phi[-W_t(u_{t,i})^\top\boldsymbol{\beta}_t(u_{t,i})]\right] \times \pi_{GP}(\boldsymbol{\beta}_{\mathcal{S}}, \boldsymbol{\beta}_{\mathcal{U}} \mid \boldsymbol{\theta}) \times \prod_{t=0}^T(\lambda_t^*)^{m_0 - 1} e^{- r_0 \lambda_t^*}\times\pi(\boldsymbol{\theta})
\end{aligned}
$$

となる。ここから各パラメータの完全条件付き分布を導出する。

### $\mathcal{U}$の完全条件付き分布

$$
\begin{aligned}
&p(\mathcal{U}, \boldsymbol{\beta}_{\mathcal{U}} | \cdot) \\&\propto \prod_{t=0}^T \left[\frac{(\lambda_t^*)^{m_t}}{m_t !}\exp\left(-\lambda_t^* \int_D 1-\Phi[W_t(u)^\top\boldsymbol{\beta}_t(u)]du\right) \prod_{i=1}^{m_t}  \left(1- \Phi[W_t(u_{t,i})^\top\boldsymbol{\beta}_t(u_{t,i})]\right)\right]
\times \pi_{GP}(\boldsymbol{\beta}_{\mathcal{U}} \mid \boldsymbol{\beta}_{\mathcal{S}}, \boldsymbol{\theta})
\\
&
\end{aligned}
$$

これは、IPPの尤度そのものである。つまり、現在の$\boldsymbol{\beta}_{\mathcal{S}}$の値に基づいてGPの条件付き分布から$\boldsymbol{\beta_{\mathcal{U}}}$の値をサンプルし、それを元にIPPから$\mathcal{U}$をサンプリングすれば良い。

IPPのサンプリングの仕方については、以下のPoisson thinningを用いる：

### IPPからのサンプリング：Poisson thinning

（Thinningによるポアソン過程の実現：[https://www.neuralengine.org/res/book/node8.html#:~:text=特徴である．-,希薄化による数値計算法,-最後に簡潔](https://www.neuralengine.org/res/book/node8.html#:~:text=%E7%89%B9%E5%BE%B4%E3%81%A7%E3%81%82%E3%82%8B%EF%BC%8E-,%E5%B8%8C%E8%96%84%E5%8C%96%E3%81%AB%E3%82%88%E3%82%8B%E6%95%B0%E5%80%A4%E8%A8%88%E7%AE%97%E6%B3%95,-%E6%9C%80%E5%BE%8C%E3%81%AB%E7%B0%A1%E6%BD%94)）

すべての $t = 1, \ldots, T$ について、以下を繰り返す：

1. $K \sim \text{Poisson}(\lambda_t^*)$ をサンプリング
2. 以下を $K$ 回繰り返す：
    1. 候補点
    $u \sim \text{Uniform}(\mathcal{D})$ をサンプリング
    2. 候補点 $u$ における共変量 $W_t(u)$ を用意し、回帰係数 $\beta_t(u) \sim \pi_{GP}(\cdot \mid \boldsymbol{\beta}_{\mathcal{S}, 1:t}, \boldsymbol{\beta}_{\mathcal{U}, 1:t-1}, \boldsymbol{\theta})$ をサンプリング
    3. 候補点 $u$ を、確率 $\Phi[-W_t(u)\beta_t(u)]$ で保持する
3. 最終的に保持された点の数を $m_t$ とし、その集合を $\mathcal{U}_t = \set{u_1, \ldots, u_{m_t}}$ とする。

こうして得られた $\mathcal{U} = \set{\mathcal{U}_1, \ldots, \mathcal{U_T}}$ は、[$\mathcal{U}$の完全条件付き分布](https://www.notion.so/mathcal-U-2412929953a980e0b5eed22af1fd80b6?pvs=21) からのモンテカルロサンプルになっている。

### $\boldsymbol{\beta}$の完全条件付き分布

$$
\begin{aligned}
p(\boldsymbol{\beta} \mid \cdot) &\propto \prod_{t=0}^T \left[\prod_{i=1}^{n_t} \Phi[W_t(s_{t,i})^\top\boldsymbol{\beta}_t(s_{t,i})]
\prod_{i=1}^{m_t}\Phi[-W_t(u_{t,i})^\top\boldsymbol{\beta}_t(u_{t,i})]\right] \times \pi_{\text{GP}}(\boldsymbol{\beta}_{\mathcal{S}}, \boldsymbol{\beta}_{\mathcal{U}}\mid \boldsymbol{\theta})
\\
&= \prod_{t=0}^T \left[\Phi_{n_t}[\boldsymbol{W}_{\mathcal{S}, t}^\top \boldsymbol{\beta}_{\mathcal{S},t}]\Phi_{m_t}[-\boldsymbol{W}_{\mathcal{U}, t}^\top \boldsymbol{\beta}_{\mathcal{U},t}]\pi_{\text{GP}}(\boldsymbol{\beta}_{\mathcal{S}, t}, \boldsymbol{\beta}_{\mathcal{U},t} \mid \boldsymbol{\beta}_{\mathcal{S}, 1:t-1}, \boldsymbol{\beta}_{\mathcal{U},1:t-1}, \boldsymbol{\theta})\right]
\\
&= \prod_{t=0}^T \left[\Phi_{n_t+m_t}[\boldsymbol{W}_{t}^\top \boldsymbol{\beta}_{t}]\pi_{\text{GP}}(\boldsymbol{\beta}_t \mid \boldsymbol{\beta}_{\mathcal{S}, 1:t-1}, \boldsymbol{\beta}_{\mathcal{U},1:t-1}, \boldsymbol{\theta})\right]
\\
&= \prod_{t=0}^T \left[\Phi_{n_t+m_t}[\boldsymbol{W}_{t}^\top \boldsymbol{\beta}_{t}]
\phi_{n_t+m_t}\left(\boldsymbol{\beta}_t-\boldsymbol{\mu}_t ; \boldsymbol{\Sigma}_t\right)
\right]
\end{aligned}
$$

1. 1行目から2行目について、多変量正規分布にまとめた。 $\Phi_{n_t}[\cdot]$ と $\Phi_{m_t}[\cdot]$ はそれぞれ、次元が $n_t$ と $m_t$ の多変量正規分布の累積分布関数である。ガウス過程の事前分布については、時系列で分解をした。
2. 2行目から3行目について、以下の変換をした：

    $$
    \begin{aligned}
    \boldsymbol{W}_{t} &= \left(\begin{array}{cc}
    \boldsymbol{W}_{\mathcal{S, t}} & \boldsymbol{O}_{n_t \times m_t} \\
    \boldsymbol{O}_{m_t \times n_t} & -\boldsymbol{W}_{\mathcal{U, t}}
    \end{array}\right)
    \\
    \boldsymbol{\beta}_t &= [\boldsymbol{\beta}_{\mathcal{S}, t}, \boldsymbol{\beta}_{\mathcal{U},t}]^\top
    \end{aligned}
    $$

    この時、多変量標準正規分布の累積分布関数の性質により、以下が成り立つ：

    $$
    \Phi_{n_t+m_t}[\boldsymbol{W}_{t}^\top \boldsymbol{\beta}_{t}] = \Phi_{n_t}[\boldsymbol{W}_{\mathcal{S}, t}^\top \boldsymbol{\beta}_{\mathcal{S},t}]\Phi_{m_t}[-\boldsymbol{W}_{\mathcal{U}, t}^\top \boldsymbol{\beta}_{\mathcal{U},t}]
    $$

3. 3行目から4行目まで、ガウス過程の条件付き分布で書き換えた。

以上より、$\boldsymbol{\beta}$ の完全条件付き分布は、Skewed Normal Distribution $S N(\boldsymbol{\mu}_t, \boldsymbol{\Sigma}_t, \boldsymbol{W}_t)$ の形である。このサンプリングについては以下で述べる：

### Skewed Normal Distributionからのサンプリング

[Scandinavian J Statistics - 2006 - ARELLANO‐VALLE - On the Unification of Families of Skew‐normal Distributions.pdf](attachment:9039248e-0c1e-40b5-a16e-2cdcbef76ac5:Scandinavian_J_Statistics_-_2006_-_ARELLANOVALLE_-_On_the_Unification_of_Families_of_Skewnormal_Distributions.pdf)

m次元の多変量歪正規分布 $S N(\xi, \Sigma, W)$ の確率密度関数は以下のようになる。

$$
f(z)=\frac{1}{\Phi_m(\gamma ; \Gamma)} \phi_d(z-\xi ; \Sigma) \Phi_m\left(W z ; I_m\right)
$$

ここで、

$$
\begin{aligned}& \gamma=W \xi \\& \Gamma=I_m+W \Sigma W^T\end{aligned}
$$

である。

- Skew-Normalからのサンプリングコード

    本論文のSection 4.1、または以下の論文をもとに作成

    [Scandinavian J Statistics - 2006 - ARELLANO‐VALLE - On the Unification of Families of Skew‐normal Distributions.pdf](attachment:722cdcc1-c752-4445-9bff-eb743656cb56:Scandinavian_J_Statistics_-_2006_-_ARELLANOVALLE_-_On_the_Unification_of_Families_of_Skewnormal_Distributions.pdf)

    ```python
    import numpy as np
    from scipy.stats import truncnorm
    from scipy.linalg import cholesky, solve_triangular

    def sample_skew_normal(xi, Sigma, W, n_gibbs_iter=5):
       """
       Sample from Skew Normal distribution SN(xi, Sigma, W)

       Parameters:
       -----------
       xi : array_like, shape (d,)
           Location parameter
       Sigma : array_like, shape (d, d)
           Scale matrix (positive definite)
       W : array_like, shape (m, d)
           Skewness matrix
       n_gibbs_iter : int
           Number of Gibbs sampling iterations for truncated normal (default: 5)

       Returns:
       --------
       z : array_like, shape (d,)
           Sample from SN(xi, Sigma, W)
       """
       xi = np.asarray(xi)
       Sigma = np.asarray(Sigma)
       W = np.asarray(W)

       m, d = W.shape

       # Step 0: Compute necessary matrices
       # Gamma = I_m + W @ Sigma @ W.T
       Gamma = np.eye(m) + W @ Sigma @ W.T

       # Delta = W @ Sigma
       Delta = W @ Sigma

       # gamma = W @ xi
       gamma = W @ xi

       # Cholesky decomposition of Gamma
       A = cholesky(Gamma, lower=True)

       # Step 1: Sample u* from truncated multivariate normal
       u_star = sample_truncated_mvn(A, gamma, n_gibbs_iter)

       # Step 2: Transform to u = A @ u*
       u = A @ u_star

       # Step 3: Sample z* from conditional normal distribution
       # Mean: Delta.T @ inv(Gamma) @ u
       # Covariance: Sigma - Delta.T @ inv(Gamma) @ Delta

       # Efficient computation using Cholesky decomposition
       # Solve A @ y = u for y
       y = solve_triangular(A, u, lower=True)
       # Solve A.T @ x = y for x (so x = inv(Gamma) @ u)
       x = solve_triangular(A.T, y, lower=False)

       # Conditional mean
       mu_cond = Delta.T @ x

       # Conditional covariance
       # First solve A @ Y = Delta for Y
       Y = solve_triangular(A, Delta, lower=True)
       # Then Sigma_cond = Sigma - Y.T @ Y
       Sigma_cond = Sigma - Y.T @ Y

       # Sample from conditional normal
       z_star = np.random.multivariate_normal(mu_cond, Sigma_cond)

       # Step 4: Shift by xi
       z = z_star + xi

       return z

    def sample_truncated_mvn(A, gamma, n_iter=5):
       """
       Sample from truncated multivariate standard normal with constraint A @ u* > -gamma
       using Gibbs sampling

       Parameters:
       -----------
       A : array_like, shape (m, m)
           Lower triangular matrix from Cholesky decomposition
       gamma : array_like, shape (m,)
           Constraint vector
       n_iter : int
           Number of Gibbs iterations

       Returns:
       --------
       u_star : array_like, shape (m,)
           Sample from truncated distribution
       """
       m = len(gamma)
       u_star = np.zeros(m)

       # Initialize with a valid point
       for i in range(m):
           # For lower triangular A, A[i,j] = 0 for j > i
           # So constraint i is: sum(A[i,j] * u_star[j] for j <= i) > -gamma[i]
           lower_bound = (-gamma[i] - np.sum(A[i, :i] * u_star[:i])) / A[i, i]
           # Sample from truncated normal(0, 1) with lower bound
           u_star[i] = truncnorm.rvs(a=lower_bound, b=np.inf, loc=0, scale=1)

       # Gibbs sampling iterations
       for _ in range(n_iter):
           for i in range(m):
               # Calculate bounds for u_star[i] from all constraints
               bounds = []

               # Constraint j involves u_star[i] only if A[j,i] != 0
               # Since A is lower triangular, this means j >= i
               for j in range(i, m):
                   # Constraint: sum(A[j,k] * u_star[k] for k <= j) > -gamma[j]
                   # Isolate u_star[i]: A[j,i] * u_star[i] > -gamma[j] - sum(...)

                   other_sum = np.sum(A[j, :i] * u_star[:i])
                   if i < j:
                       other_sum += np.sum(A[j, i+1:j+1] * u_star[i+1:j+1])

                   if A[j, i] > 0:
                       lower = (-gamma[j] - other_sum) / A[j, i]
                       bounds.append((lower, np.inf))
                   elif A[j, i] < 0:
                       upper = (-gamma[j] - other_sum) / A[j, i]
                       bounds.append((-np.inf, upper))

               # Find the tightest bounds
               lower_bound = max([b[0] for b in bounds])
               upper_bound = min([b[1] for b in bounds])

               # Sample from truncated normal
               u_star[i] = truncnorm.rvs(a=lower_bound, b=upper_bound, loc=0, scale=1)

       return u_star
    ```

- Skew-Normalの例

    正規分布を、ベクトル$W$ によって歪ませている

    ![image.png](attachment:86112b56-6e9b-4bf6-bf0d-5387700f7c7e:image.png)


### $\boldsymbol{\lambda^*}$の完全条件付き分布

$$
\begin{aligned}
p(\boldsymbol{\lambda_t^*} \mid \cdot)
&\propto \prod_{t=0}^{T}\frac{(\lambda_t^*)^{n_t + m_t}}{n_t!m_t!} \exp\left(-\lambda_t^* |\mathcal{D}|\right) \times \prod_{t=0}^T(\lambda_t^*)^{m_0 - 1} e^{- r_0 \lambda_t^*} \\
&\propto \prod_{t=0}^{T}\left[ (\lambda^*_t)^{m_0 - 1} e^{- r_0 \lambda^*_t} \cdot e^{- \lambda^*_t |\mathcal{D}|} \cdot \frac{ (\lambda^*_t)^{n_t + m_t} }{ n_t! m_t! } \right] \\
&\propto \prod_{t=0}^{T}\left[ (\lambda^*_t)^{m_0 - 1} e^{- r_0 \lambda^*_t} \cdot e^{- \lambda^*_t |\mathcal{D}|} \cdot (\lambda^*_t)^{n_t+m_t} \right]\\
&\propto \prod_{t=0}^{T}\left[ (\lambda^*_t)^{m_0 - 1 + n_t+m_t} e^{- (r_0 + |\mathcal{D}|) \lambda^*_t} \right]\\

\end{aligned}
$$

より、再びガンマ分布 $Ga(m_0 + n_t + m_t, r_0 + |\mathcal{D}|)$ となる。

### ハイパーパラメータ $\boldsymbol{\theta}$ の完全条件付き分布

ハイパーパラメータも同時にサンプルする場合は、以下のようにして行う。

$$
p(\boldsymbol{\theta} \mid \cdot) \propto \pi_{GP}(\boldsymbol{\beta}_{\mathcal{S}}, \boldsymbol{\beta}_{\mathcal{U}} \mid \boldsymbol{\theta}) \cdot \pi(\boldsymbol{\theta})
$$

これは一般的に解析的な形にならないため、MH法、HMCなどでサンプルする。

分割される領域ごとに推定する

河川で断絶された部分を障害物で入れる

カテゴリ間の相関を考える

共変量は最初は考えない
