# Walker (2007) - スライスサンプリングによるディリクレ混合モデル

## 基本情報
- **タイトル**: Sampling the Dirichlet Mixture Model with Slices
- **著者**: Stephen G. Walker
- **ジャーナル**: Working Paper Series, University of Kent (Working Paper No. 16/2006)
- **年**: 2006 (Working Paper)、2007年出版
- **DOI**: 記載なし
- **関連論点**: 論点2（ベイズ推論手法、サンプリングアルゴリズム）

## 主な貢献
- **潜在変数（スライス変数）による無限和の有限化**: 補助変数$u$を導入し、各イテレーションで必要なコンポーネント数$k^*$を正確に決定
- **シンプルなGibbsサンプラー**: 受理・棄却ステップや複雑なretrospective sampling不要
- **既存手法の改良**: Ishwaran & Zarepour (2000)の近似打ち切りやPapaspiliopoulos & Roberts (2005)のdetailed balance基準より実装が容易
- **一般的なスティックブレーキング事前分布への拡張可能性**: $v_j \sim \text{Beta}(a_j, b_j)$の一般形にも適用可能

## 手法の概要

### ディリクレ過程混合モデル（MDP）
ディリクレ過程$D(c, P_0)$からの確率測度$P$はスティックブレーキング表現（Sethuraman & Tiwari 1982, Sethuraman 1994）により：

$$
P = \sum_{j=1}^{\infty} w_j \delta_{\theta_j}
$$

ここで、
- $v_j \sim \text{Beta}(1, c)$ （独立）
- $\theta_j \sim P_0$ （独立、密度$g_0$）
- $w_1 = v_1$, $w_j = v_j \prod_{l<j}(1-v_l)$ （$j > 1$）

MDP密度関数：

$$
f_P(y) = \int \text{N}(y|\theta) \, dP(\theta) = \sum_{j=1}^{\infty} w_j \text{N}(y|\theta_j)
$$

### 潜在変数の導入（スライスサンプリングのアイデア）
補助変数$u$を導入し、$(y, u)$の同時密度を定義：

$$
f_{w,\theta}(y, u) = \sum_{j=1}^{\infty} \mathbf{1}(u < w_j) \text{N}(y|\theta_j)
$$

これは以下のように書き換えられる：

$$
f_{w,\theta}(y, u) = \sum_{j=1}^{\infty} w_j \text{U}(u|0, w_j) \text{N}(y|\theta_j)
$$

**重要な性質**: 集合$A_w(u) = \{j : w_j > u\}$は全ての$u > 0$に対して**有限**

$u$の周辺密度：

$$
f_w(u) = \sum_{j=1}^{\infty} w_j \text{U}(u|0, w_j) = \sum_{j=1}^{\infty} \mathbf{1}(u < w_j) = |A_w(u)|
$$

$u$を条件とする$y$の条件付き密度：

$$
f_{w,\theta}(y|u) = \frac{1}{f_w(u)} \sum_{j \in A_w(u)} \text{N}(y|\theta_j)
$$

→ **有限混合モデル**（等しい重み$1/f_w(u)$）

### 指示変数の導入
さらに混合コンポーネントを識別する指示変数$\delta$を導入：

$$
f_{w,\theta}(y, \delta=k, u) = \text{N}(y|\theta_k) \mathbf{1}(k \in A_w(u))
$$

完全データ尤度（サンプルサイズ$n$）：

$$
l_{w,\theta}(\{y_i, u_i, \delta_i=k_i\}_{i=1}^n) = \prod_{i=1}^n \text{N}(y_i|\theta_{k_i}) \mathbf{1}(u_i < w_{k_i})
$$

### Gibbsサンプラーの条件付き分布

**A. 潜在変数$u_i$のサンプリング**:

$$
u_i | \cdots \sim \text{Uniform}(0, w_{k_i})
$$

**B. パラメータ$\theta_j$のサンプリング**:

$$
f(\theta_j | \cdots) \propto g_0(\theta_j) \prod_{k_i=j} \text{N}(y_i|\theta_j)
$$

$k_i = j$となる$i$が存在しない場合、$f(\theta_j|\cdots) = g_0(\theta_j)$

**C. スティック長$v_j$のサンプリング**:

$$
f(v_j | v_{-j}, \cdots) \propto \text{Beta}(v_j|1, c) \mathbf{1}(\alpha_j < v_j < \beta_j)
$$

ここで、

$$
\alpha_j = \max_{k_i=j} \left\{ \frac{u_i}{\prod_{l<j}(1-v_l)} \right\}
$$

$$
\beta_j = 1 - \max_{k_i>j} \left\{ \frac{u_i}{v_{k_i} \prod_{l<k_i, l \neq j}(1-v_l)} \right\}
$$

$j \leq k^*$（$k^* = \max\{k_1, \ldots, k_n\}$）に対してのみサンプリングが必要。$j > k^*$に対しては$v_j \sim \text{Beta}(1, c)$

分布関数（逆CDF法によるサンプリング）：

$$
F(v_j) = \frac{(1-\alpha_j)^c - (1-v_j)^c}{(1-\alpha_j)^c - (1-\beta_j)^c}
$$

**D. 指示変数$\delta_i$のサンプリング**:

$$
\Pr(\delta_i = k | \cdots) \propto \mathbf{1}(k \in A_w(u_i)) \text{N}(y_i|\theta_k)
$$

**重要**: $A_w(u_i)$は有限集合であり、必要なコンポーネント数$k^*$は以下の条件で決定：

$$
\sum_{j=1}^{k^*} w_j > 1 - u^*
$$

ここで$u^* = \min\{u_1, \ldots, u_n\}$

事前分布の下での$k^*$の分布：

$$
k^* \sim 1 + \text{Poisson}(-c \log u^*)
$$

（Muliere & Tardella 1998参照）

**E. スケールパラメータ$c$のサンプリング**（事前分布$\pi(c)$を設定した場合）:

$$
f(c | d, n) \propto c^d \Gamma(c) \pi(c) / \Gamma(c+n)
$$

ここで$d$はクラスタ数（異なる$k_i$の個数）、$n$はサンプルサイズ

### 予測分布のサンプリング
予測分布$f(y_{n+1}|y_1, \ldots, y_n)$からのサンプリング：

1. 各イテレーションで$(w_j, \theta_j)$を得る
2. 一様乱数$r \sim \text{Uniform}(0, 1)$をサンプル
3. $w_{j-1} < r < w_j$となる$\theta_j$を選択（$w_0 = 0$）
4. 必要なら追加の$v_j \sim \text{Beta}(1, c)$と$\theta_j \sim g_0$をサンプル
5. $y_{n+1} \sim \text{N}(\cdot|\theta_j)$

## 実証結果

### シミュレーション例（Gaussian MDP）
**データ生成モデル**:

$$
f(y) = \frac{1}{3} \text{N}(y|-4, 1) + \frac{1}{3} \text{N}(y|0, 1) + \frac{1}{3} \text{N}(y|8, 1)
$$

**設定**:
- サンプルサイズ: $n = 50$
- パラメータ: $\epsilon = 0.5$, $s = 0.1$（非情報的事前分布）
- $c$の事前分布: $\text{Ga}(0.1, 0.1)$
- イテレーション: 20,000回
- バーンイン: 10,000回
- 予測サンプル: 10,000個（イテレーション10,000以降）

**結果**:
- Figure 1: データのヒストグラムと予測密度推定が良好に一致
- Figure 2: クラスタ数の移動平均が10,000イテレーション以内に定常状態に到達
- 密度推定: R density routineでバンド幅0.3を使用

**実装**: Scilab（フリーソフト）でコーディング、シンプルで高速

### 正規分布の場合の条件付き分布
$\theta = (\mu, \sigma^2)$、$\lambda = \sigma^{-2}$として：

**$\mu_j$の条件付き分布**:

$$
f(\mu_j | \cdots) = \text{N}\left(\frac{\xi_j \lambda_j}{m_j \lambda_j + s}, \frac{1}{m_j \lambda_j + s}\right)
$$

ここで$\xi_j = \sum_{k_i=j} y_i$、$m_j = \sum_{k_i=j} 1$

**$\lambda_j$の条件付き分布**:

$$
f(\lambda_j | \cdots) = \text{Ga}\left(\epsilon + m_j/2, \epsilon + d_j/2\right)
$$

ここで$d_j = \sum_{k_i=j}(y_i - \theta_j)^2$

## 限界・残された課題

### 理論的・実装上の限界

1. **非共役モデルへの拡張**: $\text{N}(y|\theta)$と$g_0(\theta)$が非共役対の場合
   - 条件付きサンプリングが困難
   - **解決策**: Damien et al. (1999, Sections 4 & 5)の補助変数法を適用

2. **スティックブレーキング事前分布の一般化**: 本論文は$v_j \sim \text{Beta}(1, c)$に焦点
   - $v_j \sim \text{Beta}(a_j, b_j)$への拡張は可能だが詳細は省略
   - Ishwaran & James (2001)の条件: $\sum_{j=1}^{\infty} \log(1 + a_j/b_j) = +\infty$

3. **$k^*$の決定の計算コスト**: 各イテレーションで$k^*$を見つける必要
   - $\sum_{j=1}^{k^*} w_j > 1 - u^*$の判定
   - ただし、$k^* \sim 1 + \text{Poisson}(-c \log u^*)$は比較的小さい

4. **既存手法との比較**: 定量的な比較実験が限定的
   - collapsed Gibbs samplerやretrospective samplingとの速度・精度比較が不十分
   - 単一のシミュレーション例のみ提示

5. **高次元データへのスケーラビリティ**: 実証例は1次元データのみ
   - 多変量データや高次元での性能が不明

### Pólya-urnモデルとの関係

6. **依存構造の扱い**:
   - Pólya-urn scheme（Blackwell & MacQueen 1973）: $\theta_{k_i}$間に依存
   - 本手法: ランダム分布関数$P$を保持することで依存を除去
   - トレードオフ: 無限次元表現への対処が必要

7. **周辺化との比較**: Escobar (1988, 1994)のcollapsed Gibbsは$P$を周辺化
   - 本手法はより多くの変数（$v_j$, $\theta_j$, $u_i$, $\delta_i$）をサンプリング
   - 各条件付き分布が単純化されるメリット

## 本研究（MMCP）との関係

### 借用すべき手法

1. **潜在変数による有限化のアイデア**:
   - 無限混合モデルを各イテレーションで有限化
   - MMCPでスティックブレーキング表現を使う場合に適用可能
   - 特に空間Cox過程でDirichlet processを使う場合に有用

2. **シンプルなGibbs sampler設計**:
   - 受理・棄却不要
   - 全ての条件付き分布が標準的な分布または簡単にサンプリング可能
   - 実装の簡便性とデバッグの容易さ

3. **$k^*$の正確な決定法**:
   - 近似打ち切りではなく、必要なコンポーネント数を正確に決定
   - $k^* \sim 1 + \text{Poisson}(-c \log u^*)$の分布を利用した事前予測

4. **スライスサンプリングの一般原理**:
   - Damien et al. (1999)の補助変数法の枠組み
   - 他の無限次元モデルへの適用可能性

### 拡張すべき点

1. **空間構造への適応**:
   - 本論文はi.i.d. MDP
   - MMCPは空間的に構造化されたCox過程 → 空間相関の導入が課題

2. **マーク付き点過程への拡張**:
   - 本論文は単変量観測（1次元正規分布）
   - MMCPは多変量マーク（黒曜石組成）→ マルチタスク拡張が必要

3. **空間的に変化する係数（SVC）との統合**:
   - MMCPは組成係数が空間的に変化
   - スライスサンプリングをSVCモデルに組み込む方法の検討

4. **大規模データへの適用性の評価**:
   - 遺跡データが数千規模になった場合の計算効率
   - $k^*$の増加による計算コスト増への対策

5. **非共役モデルへの拡張**:
   - MMCPで非共役事前分布を使う場合（例：ロジスティック回帰）
   - Damien et al. (1999)の補助変数法の適用

### MMCPにおける位置づけ

**第3章での役割**:
- ノンパラメトリックベイズのサンプリング手法の1つとして紹介
- Blei & Jordan (2006)の変分推論との対比材料（MCMC vs 変分）
- スティックブレーキング表現の実装例

**実装の選択肢**:
- Stan/NUTSやcollapsed Gibbsの代替
- 特に無限混合コンポーネントを明示的に扱いたい場合に有用
- Retrospective samplingより実装が容易

**トレードオフの理解**:
- シンプルさ（本手法）vs 周辺化による効率（collapsed Gibbs）
- 明示的な分布関数（本手法）vs 暗黙的な依存構造（Pólya-urn）

## 引用すべき箇所

### 主要な貢献の説明

> "The key to the algorithm detailed in this paper, which also keeps the random distribution functions, is the introduction of a latent variable which allows a finite number, which is known, of objects to be sampled within each iteration of a Gibbs sampler." (p.1, Abstract)

→ 本手法の核心的アイデアの要約

> "We find a simple trick, based on the slice sampling schemes (Damien et al., 1999), which deals with the infiniteness. The introduction of latent variables makes finite the part of the random distribution function required to iterate through a Gibbs sampler." (p.1, Introduction)

→ スライスサンプリングの利用と有限化の原理

### 既存手法との比較

> "Ishwaran & Zarepour (2000) circumvent this [countable infiniteness] with an approximate method based on a truncation of the distributions. [...] Papaspiliopoulos & Roberts (2005) proposed an exact algorithm based on the notion of retrospective sampling. However, the algorithm itself becomes non-trivial when applied to the MDP model, and involves setting up a detailed balance criterion with connecting proposal moves." (p.1, Introduction)

→ 既存手法の限界と本手法の優位性

### 重要な性質

> "Note, it is quite clear that $A_w(u) = \{j: w_j > u\}$ is a finite set for all $u > 0$." (p.3, Section 2)

→ 有限化の数学的根拠

> "So we sample as many of the $w_k$'s until we are sure that we have all the $w_k > u_i$. [...] Hence, $k^*$ is not a loose approximation; it is an exact piece of information." (p.6, Section 3)

→ $k^*$が近似ではなく正確な値であることの強調

### アルゴリズムの利点

> "All of the conditional distributions are easy to sample and no accept/reject methods are needed." (p.1, Introduction)

→ 実装の容易さ

> "Hence, all the conditional densities are easy to sample and the Markov chain we have constructed is automatic. It requires no tuning nor retrospective steps." (p.7, Section 3)

→ チューニング不要の自動化

### 結論とまとめ

> "We have provided a simple and fast way to sample the MDP model. [...] It improves on current approaches in the following way: we know exactly how many of the $w_j$'s and $\theta_j$'s we need to sample at each iteration—it is $k^*$. This fundamental result eludes the alternative approaches." (p.9, Discussion)

→ 本手法の決定的な優位性

> "Retaining the random distribution function is useful as it removes the dependence between the $\theta_{k_i}$'s which exist in the Pólya-urn model." (p.9, Discussion)

→ ランダム分布関数を保持する動機

### 非共役への拡張

> "In the non-conjugate case, that is when $N(y|\theta)$ and $g_0(\theta)$ form a nonconjugate pair and perhaps difficult to sample, then a possible useful solution is again provided by the latent variable ideas presented in Damien et al. (1999, Sections 4 & 5)." (p.9, Discussion)

→ 非共役モデルへの拡張可能性

### スライスサンプリングの一般原理

> "The usefulness of the latent variable $u$ will become clear later on. A brief comment here is that the move from an infinite sum to a finite sum, given $u$, is going to make a lot of difference when sampling is involved." (p.3, Section 2)

→ 補助変数$u$の役割の本質的説明
