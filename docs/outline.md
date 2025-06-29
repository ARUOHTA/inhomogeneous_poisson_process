# 空間多項分布カウントデータに対するベイズノンパラメトリック手法

**著者**: 太田 歩
**所属**: 京都大学大学院情報学研究科

**共同研究者**: 原 尚幸（京都大学国際高等教育院）

## 1. はじめに

空間多項分布カウントデータは多様な科学分野で生じるデータであり、異なる空間位置での観測が位置依存パラメータを持つ多項分布に従う。このようなデータ構造は、生態学の種組成研究、考古学の人工物分布解析、疫学の疾病分類研究において広く見られる。従来のアプローチは多くの場合、パラメトリックな仮定や単純なカーネル平滑化に依存しており、複雑な空間依存性や根本的な確率過程における不確実性を適切に捉えることができない。

カーネル・スティック・ブレイキング過程（KSBP）は柔軟なノンパラメトリックベイズフレームワークを提供するが、空間多項分布データへの応用は限定的である。本研究では、KSBPフレームワークを拡張して空間多項分布カウントデータを扱い、位置依存多項分布パラメータに対する包括的な不確実性定量化を行う多項分布カーネル・スティック・ブレイキング過程（Multinomial-KSBP）を提案する。

## 2. 問題設定

空間領域 $\mathcal{D} \subset \mathbb{R}^2$ および観測位置 $s_1, \ldots, s_n \in \mathcal{D}$ を考える。各位置 $s_i$ において、総カウント数 $N_i = \sum_{k=1}^K y_{ik}$ を持つ多項分布カウントデータ $\mathbf{y}_i = (y_{i1}, \ldots, y_{iK})^\top \in \mathbb{N}^K$ を観測する。

目標は、任意の位置 $s \in \mathcal{D}$ における多項分布確率ベクトル $\boldsymbol{\pi}(s) = (\pi_1(s), \ldots, \pi_K(s))^\top$ を推定することであり、ここで $\sum_{k=1}^K \pi_k(s) = 1$ かつ $\pi_k(s) > 0$ である。

## 3. 多項分布カーネル・スティック・ブレイキング過程

位置依存多項分布パラメータを空間カーネルを持つ無限混合を通じてモデル化するMultinomial-KSBPを提案する：

1. 空間位置パラメータ: $\Gamma_h \sim H$ （$\mathcal{D}$ 上の一様分布）
2. スティック・ブレイキングパラメータ: $V_h \sim \text{Beta}(1, \lambda)$
3. 成分パラメータ: $\boldsymbol{\theta}_h \sim \text{Dirichlet}(\gamma_0/K \cdot \mathbf{1}_K)$

位置依存確率ベクトルは $\boldsymbol{\pi}(s) = \sum_{h=1}^{\infty} \pi_h(s) \boldsymbol{\theta}_h$ であり、空間適応的重みは：

$$
\pi_h(s) = V_h K(s, \Gamma_h) \prod_{l=1}^{h-1} [1 - V_l K(s, \Gamma_l)]
$$

ここで、$K(s, \Gamma_h)$ は空間依存性を制御するガウシアンカーネルである。$P(z_i = h \mid \mathbf{V}, \boldsymbol{\Gamma}) = \pi_h(s_i)$ を満たす補助変数 $z_i$ を導入すると：$z_i \mid \mathbf{V}, \boldsymbol{\Gamma}, s_i \sim \text{Categorical}(\pi_1(s_i), \pi_2(s_i), \ldots)$ および $\mathbf{y}_i \mid z_i, \{\boldsymbol{\theta}_h\} \sim \text{Multinomial}(N_i, \boldsymbol{\theta}_{z_i})$ となる。

事後推論については、無限混合構造を扱うためのスライスサンプリング技法を組み込んだ効率的なギブスサンプリングアルゴリズムを開発し、空間位置パラメータにはメトロポリス・ヘイスティングス更新を用いる。

## 4. モデルの性質と実証結果

Multinomial-KSBPは理論的に望ましい性質を示す：カーネル関数による空間局所性、無限混合による非パラメトリック柔軟性、ベイズフレームワーク内での自然な不確実性定量化である。$K(s, \Gamma_h) \equiv 1$ の場合、標準的なディリクレ過程混合に帰着し、その一般性を示している。

日本の考古学的黒曜石データにモデルを適用し、不確実性推定を伴う空間パターンの捕捉に成功した。データが疎な地域では不確実性が高くなり、データ不足の原理的な処理を実証している。結果は黒曜石交易ネットワークの時間的進化を明らかにし、先史時代の交換システムに関する新たな洞察を提供する。本アプローチは、厳密な不確実性定量化を伴う柔軟なノンパラメトリックベイズモデリングを通じて、空間多項分布解析における重要な進歩を表している。今後の研究では、セミパラメトリック拡張やガウス過程との関連を検討する。

## 参考文献

1. Dunson, D. B. and Park, J. H. (2008). *Kernel stick-breaking processes*. Biometrika, 95(2), 307--323.
2. Walker, S. G. (2007). *Sampling the Dirichlet mixture model with slices*. Communications in Statistics - Simulation and Computation, 36(1), 45--54.
3. Polson, N. G., Scott, J. G., and Windle, J. (2013). *Bayesian inference for logistic models using Pólya–Gamma latent variables*. Journal of the American Statistical Association, 108(504), 1339--1349.
