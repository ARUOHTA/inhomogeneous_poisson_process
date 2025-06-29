<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.min.js">
</script>
<script type="text/x-mathjax-config">
 MathJax.Hub.Config({
 tex2jax: {
 inlineMath: [['$', '$'] ],
 displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
 }
 });
</script>

# 黒曜石分布解析のための統計モデル

本文書では、考古学的な黒曜石分布データを分析するために開発された統計モデルの数学的定式化について説明する。

## 1. 問題設定

### 1.1 データの構造

- **対象領域**: $\mathcal{D} \subset \mathbb{R}^2$（関東地方：北緯34-37度、東経138-141度）
- **遺跡データ**:
  - 遺跡数: $n_X = 224$
  - 各遺跡の位置: $s_i \in \mathcal{D}$ $(i = 1, \ldots, n_X)$
  - 各遺跡における産地 $k$ の黒曜石出土数: $y_{ik}$ $(k = 1, \ldots, K)$
- **産地分類**: $K = 4$（神津島、信州、箱根、高原山）

### 1.2 研究目的

任意の地点 $s \in \mathcal{D}$ における黒曜石の産地構成比ベクトル

$$\boldsymbol{\pi}(s) = (\pi_1(s), \ldots, \pi_K(s)), \quad \sum_{k=1}^{K} \pi_k(s) = 1$$

を推定することが目的である。そのために以下の2つのモデルを構築する：

1. **遺跡の存在確率モデル**: 遺跡の空間的分布をモデル化
2. **産地構成比モデル**: 各地点における黒曜石の産地別比率をモデル化

## 2. 遺跡の存在確率モデル（非斉次ポアソン過程）

### 2.1 非斉次ポアソン過程の定義

$X$ を領域 $\mathcal{D}$ 上の計数過程とし、$\lambda(s)$ を $\mathcal{D}$ 上の連続な非負実数値可測関数とする。

**定義（非斉次ポアソン過程）**: 任意の可測集合 $D \subset \mathcal{D}$ に対して、

$$\Lambda(D) := \int_D \lambda(s) ds$$

とするとき、以下の条件を満たす計数過程 $X$ を強度 $\lambda(s)$ を持つ非斉次ポアソン過程という：

1. $X(D) \sim \text{Poisson}(\Lambda(D))$
2. 互いに素な集合 $D_1, D_2, \ldots, D_k$ に対して、$X(D_1), X(D_2), \ldots, X(D_k)$ は互いに独立

記法: $X \sim \text{IPP}(\lambda)$

### 2.2 尤度関数

観測データは領域 $\mathcal{D}$ 全体のイベント発生数 $n_X$ と、それぞれのイベント発生位置 $s_i$ $(i = 1, \ldots, n_X)$ の組である。非斉次ポアソン過程の尤度関数は：

$$P(X \mid \lambda) = \exp\left( -\int_{\mathcal{D}} \lambda(s) \, ds \right) \cdot \frac{1}{n_X!} \cdot \prod_{i=1}^{n_X} \lambda(s_i)$$

### 2.3 潜在変数アプローチ（Moreira and Gamerman, 2022）

積分項の計算困難性を回避するため、強度関数を以下のように分解：

$$\lambda(s) = \lambda^* \cdot q(s), \quad s \in \mathcal{D}$$

ここで：
- $\lambda^* > 0$: 強度関数の上限
- $q(s)$: 相対的な強度（0から1の値）

ロジスティック回帰モデルとして：

$$q(s) = \frac{\exp(\boldsymbol{W}(s)^\top \boldsymbol{\beta})}{1 + \exp(\boldsymbol{W}(s)^\top \boldsymbol{\beta})}$$

ここで $\boldsymbol{W}(s)$ は位置 $s$ での説明変数ベクトル（標高、傾斜角度、各産地からの距離など）。

### 2.4 潜在変数による拡張

偽不在の点過程 $U$ を潜在変数として導入：

$$U \sim \text{IPP}(\lambda^*(1-q))$$

同時尤度は：

$$P(X, U \mid \boldsymbol{\beta}, \lambda^*) = \exp(-\lambda^* |\mathcal{D}|) \cdot \frac{(\lambda^*)^{n}}{n_X! n_U!} \prod_{i=1}^{n} \frac{\{\exp(\boldsymbol{W}(s_i)^\top \boldsymbol{\beta})\}^{y_i}}{1 + \exp(\boldsymbol{W}(s_i)^\top \boldsymbol{\beta})}$$

ただし、$n = n_X + n_U$、$y_i$ は在/偽不在を表す二値変数。

### 2.5 事後分布とギブスサンプリング

事後分布：

- $\boldsymbol{\beta} \mid \boldsymbol{\omega}, X, U, \lambda^* \sim \mathcal{N}(\boldsymbol{m}, V)$
  - $V = (B_0^{-1} + W^\top \Omega W)^{-1}$
  - $\boldsymbol{m} = V(B_0^{-1} \boldsymbol{b}_0 + W^\top \Omega \boldsymbol{z})$
- $\lambda^* \mid \boldsymbol{\beta}, \boldsymbol{\omega}, X, U \sim \text{Ga}(m_0 + n, r_0 + |\mathcal{D}|)$
- $U \mid \boldsymbol{\beta}, \lambda^* \sim \text{IPP}(\lambda^*(1-q))$
- $\omega_i \mid \boldsymbol{\beta}, X, U \sim \text{PG}(1, \boldsymbol{W}(s_i)^\top \boldsymbol{\beta})$

ここで、$W$ は説明変数行列、$\Omega = \text{diag}(\omega_1, \ldots, \omega_n)$、$\text{PG}$ はPolya-Gamma分布。

## 3. 産地構成比モデル1：Nadaraya-Watson推定量

### 3.1 モデルの定式化

各遺跡における産地構成比を以下の回帰モデルで表現：

$$\frac{y_{ik}}{\sum_{k'} y_{ik'}} = \pi_k(s_i, \tilde{\boldsymbol{W}}(s_i)) + \epsilon$$

Nadaraya-Watson推定量による推定：

$$\hat{p}_k(s, \boldsymbol{w}) = \frac{\sum_{i=1}^{n_X} K_h(s-s_i) \prod_{l=1}^{p} K_h(w_l - \tilde{W}_l(s_i)) \cdot y_{ik}}{\sum_{i=1}^{n_X} K_h(s-s_i) \prod_{l=1}^{p} K_h(w_l - \tilde{W}_l(s_i)) \cdot \sum_{k'} y_{ik'}}$$

### 3.2 カーネル関数

ガウスカーネルを使用：

$$K_h(s-s_i) = \frac{1}{h^2} \exp\left(-\frac{d(s-s_i)^2}{2h^2}\right)$$

ここで、$h$ はバンド幅、$d(\cdot)$ はTobler's Hiking Functionに基づくコスト関数。

### 3.3 Tobler's Hiking Function

隣接する2地点間の移動速度（km/h）：

$$W = 6e^{-3.5|\frac{dh}{dx} + 0.05|}$$

ここで、$\frac{dh}{dx} = S = \tan\theta$ は勾配。

移動コストの設定：
- 陸上：上記の関数で計算
- 海上：一律4km/h
- 陸海境界：地形により異なるコスト（砂浜海岸は通常の50倍、岩石海岸は無限大）

### 3.4 最短距離の計算

1. 領域 $\mathcal{D}$ を5次メッシュコード（250m単位）でグリッド分割
2. 各グリッドの中心点を頂点とするグラフを構築
3. 隣接セル間の移動コストを辺の重みとして設定
4. 最短経路アルゴリズムにより全頂点間の最短コスト $C_u(v)$ を計算：

$$C_u(v) \leftarrow \min(\{C_u(v)\} \cup \{C_u(w) + t_{w \rightarrow v} \mid w \in \mathcal{N}_v\})$$

5. 任意の2点 $s, s'$ 間の距離：

$$d(s, s') := \frac{C_{v_s}(v_{s'}) + C_{v_{s'}}(v_s)}{2}$$

## 4. 産地構成比モデル2：Multinomial Kernel Stick-Breaking Process

### 4.1 Kernel Stick-Breaking Process (KSBP) の基礎

説明変数 $x \in \mathcal{X}$ に依存する確率測度の族 $\{G_x\}$ に対して、以下の独立なランダム要素の無限列を導入：

$$\{\Gamma_h, V_h, G_h^*\}_{h=1}^{\infty}$$

- $\Gamma_h \sim H$：位置パラメータ
- $V_h \sim \text{Beta}(a_h, b_h)$：重み生成用の確率変数
- $G_h^* \sim Q$：確率測度

KSBPは以下で定義される：

$$G_x = \sum_{h=1}^{\infty} \pi_h(x; V_h, \Gamma_h) G_h^*$$

ここで、

$$\pi_h(x; V_h, \Gamma_h) = W(x; V_h, \Gamma_h) \prod_{l < h} [1 - W(x; V_l, \Gamma_l)]$$

$$W(x; V_h, \Gamma_h) = V_h \cdot K(x, \Gamma_h)$$

### 4.2 Multinomial-KSBPモデル

#### 階層構造

観測空間 $\mathcal{X} \subset \mathbb{R}^d$ において、以下の階層モデルを構築：

$$\begin{aligned}
\Gamma_h &\sim H \\
V_h &\sim \text{Beta}(1, \lambda) \\
\theta_h &\sim \text{Dirichlet}\left(\frac{\gamma_0}{K}\mathbf{1}_K\right) \\
\pi_h(s) &= V_h K(s, \Gamma_h) \prod_{l < h} [1 - V_l K(s, \Gamma_l)] \\
\pi(s) &= \sum_h \pi_h(s) \theta_h \\
\mathbf{y}_i \mid \pi(s_i) &\sim \text{Multinomial}(N_i, \pi(s_i))
\end{aligned}$$

ここで、$\theta_h$ はクラスタ $h$ における産地構成比、$K(s, \Gamma)$ はガウスカーネル：

$$K(s, \Gamma) = \exp\left\{-\frac{\|s - \Gamma\|^2}{2h^2}\right\}$$

### 4.3 補助変数の導入と推論

無限混合を扱うため、補助変数 $z_i$ を導入：

$$z_i \mid \mathbf{V}, \Gamma, s_i \sim \text{Categorical}(\pi_1(s_i), \pi_2(s_i), \ldots)$$

$$\mathbf{y}_i \mid z_i, \{\theta_h\} \sim \text{Multinomial}(N_i, \theta_{z_i})$$

#### 完全条件付き分布

1. **$\theta_h$ の事後分布**：

$$\theta_h \mid \mathbf{z}, \mathbf{y} \sim \text{Dirichlet}\left(\frac{\gamma_0}{K} + S_{h1}, \ldots, \frac{\gamma_0}{K} + S_{hK}\right)$$

ここで、$S_{hk} = \sum_{i: z_i = h} y_{ik}$

2. **$z_i$ の事後分布**：

$$P(z_i = h \mid \text{rest}) \propto \pi_h(s_i) \prod_{k=1}^{K} \theta_{hk}^{y_{ik}}$$

3. **$V_h$ のサンプリング**（Walker, 2007のSlice Sampling）：

$$V_h \mid \text{rest} \sim \text{Beta}(1 + m_h, \lambda + r_h)$$

ここで、
- $m_h = \#\{i : z_i = h\}$
- $r_h = \#\{i : z_i > h, u_i < V_h K(s_i, \Gamma_h) \prod_{g<h}(1-V_g K(s_i, \Gamma_g))\}$

4. **$\Gamma_h$ の更新**：Metropolis-Hastings法を使用

### 4.4 実装上の設定

- ハイパーパラメータ：
  - $\lambda = 1$（スティック長のベータ事前分布の尺度パラメータ）
  - $\gamma_0 = 0.1$（ディリクレ事前分布の集中度パラメータ）
- 共変量：標高、傾斜角度、各産地からの距離、最寄りの川・湖からの距離
- MCMCイテレーション数：2000

## 5. まとめ

本研究では、考古学的な黒曜石分布データの分析のために以下の統計モデルを開発した：

1. **遺跡の存在確率モデル**：非斉次ポアソン過程に基づき、潜在変数アプローチによる効率的な推論を実現
2. **産地構成比モデル**：
   - 頻度主義的アプローチ：Nadaraya-Watson推定量とTobler距離を組み合わせた空間的平滑化
   - ベイズ的アプローチ：Multinomial-KSBPによる不確実性の定量化

これらのモデルにより、任意の地点における黒曜石の産地構成比を推定し、その不確実性を評価することが可能となった。
