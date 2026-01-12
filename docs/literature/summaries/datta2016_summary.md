# Datta et al. (2016) - Nearest Neighbor Gaussian Process

## 基本情報
- **タイトル**: Hierarchical Nearest-Neighbor Gaussian Process Models for Large Geostatistical Datasets
- **著者**: Abhirup Datta, Sudipto Banerjee, Andrew O. Finley, Alan E. Gelfand
- **ジャーナル**: Journal of the American Statistical Association (JASA)
- **年**: 2016
- **DOI**: 確認中
- **関連論点**: 論点1（空間的依存性）

## 主な貢献

1. **NNGP（Nearest Neighbor Gaussian Process）の理論的基盤**
   - 有向非巡回グラフ（DAG）に基づく疎な条件付き密度の構築
   - Kolmogorov整合性による正当な空間過程の保証
   - 疎精度行列を持つ正当なGaussian Processの構成

2. **大規模データへのスケーラビリティ**
   - 計算複雑度: $O(nm^3)$（$n$=観測点数、$m$=近傍数、$m \ll n$）
   - フルGP（$O(n^3)$）から劇的な計算量削減
   - 低ランク手法（$O(nr^2)$）を超える柔軟性とスケーラビリティ

3. **階層ベイズフレームワークへの統合**
   - NNGPを事前分布として使用可能（非退化な過程）
   - 空間変化係数モデル（spatially-varying coefficients）への自然な拡張
   - MCMCアルゴリズムの効率的実装

4. **低ランクモデルに対する優位性**
   - 正当な空間過程（低ランクは近似）
   - 任意位置での予測が自然に可能
   - 近傍観測が強く相関する場合でも性能維持

## 手法の概要

### NNGP の構築

**基本アイデア**: 多変量正規分布の同時密度を条件付き密度の積で表現し、各条件付けに用いる集合を小さな近傍集合に制限する。

**Step 1**: 参照集合$\mathcal{S} = \{\mathbf{s}_1, \ldots, \mathbf{s}_k\}$上の完全な同時密度：
$$
p(\mathbf{w}_\mathcal{S}) = p(\mathbf{w}(\mathbf{s}_1)) p(\mathbf{w}(\mathbf{s}_2) | \mathbf{w}(\mathbf{s}_1)) \cdots p(\mathbf{w}(\mathbf{s}_k) | \mathbf{w}(\mathbf{s}_{k-1}), \ldots, \mathbf{w}(\mathbf{s}_1))
$$

**Step 2**: 各点$\mathbf{s}_i$に対し、近傍集合$N(\mathbf{s}_i) \subset \mathcal{S} \setminus \{\mathbf{s}_i\}$（$|N(\mathbf{s}_i)| \leq m$）を定義：
$$
\tilde{p}(\mathbf{w}_\mathcal{S}) = \prod_{i=1}^k p(\mathbf{w}(\mathbf{s}_i) | \mathbf{w}_{N(\mathbf{s}_i)})
$$

**Step 3**: 親Gaussian Process $\mathbf{w}(\mathbf{s}) \sim GP(\mathbf{0}, \mathbf{C}(\cdot, \cdot | \boldsymbol{\theta}))$の場合：
$$
\tilde{p}(\mathbf{w}_\mathcal{S}) = \prod_{i=1}^k N(\mathbf{w}(\mathbf{s}_i) | \mathbf{B}_{\mathbf{s}_i} \mathbf{w}_{N(\mathbf{s}_i)}, \mathbf{F}_{\mathbf{s}_i})
$$

ここで：
- $\mathbf{B}_{\mathbf{s}_i} = \mathbf{C}_{\mathbf{s}_i, N(\mathbf{s}_i)} \mathbf{C}_{N(\mathbf{s}_i)}^{-1}$（回帰係数）
- $\mathbf{F}_{\mathbf{s}_i} = \mathbf{C}(\mathbf{s}_i, \mathbf{s}_i) - \mathbf{C}_{\mathbf{s}_i, N(\mathbf{s}_i)} \mathbf{C}_{N(\mathbf{s}_i)}^{-1} \mathbf{C}_{N(\mathbf{s}_i), \mathbf{s}_i}$（条件付き分散）

### 疎精度行列の構造

$\tilde{p}(\mathbf{w}_\mathcal{S})$は多変量正規分布$N(\mathbf{0}, \tilde{\mathbf{C}}_\mathcal{S})$で、その精度行列$\tilde{\mathbf{C}}_\mathcal{S}^{-1}$は疎：
- 非ゼロ要素数: 最大$km(m+1)q^2/2$（$q$=変量数）
- $m \ll k$の場合、劇的なスパース性

### 空間過程への拡張

任意の点$\mathbf{u} \notin \mathcal{S}$に対し、$N(\mathbf{u}) \subset \mathcal{S}$を$m$-最近傍集合として定義：
$$
\tilde{p}(\mathbf{w}_\mathcal{U} | \mathbf{w}_\mathcal{S}) = \prod_{i=1}^r N(\mathbf{w}(\mathbf{u}_i) | \mathbf{B}_{\mathbf{u}_i} \mathbf{w}_{N(\mathbf{u}_i)}, \mathbf{F}_{\mathbf{u}_i})
$$

任意の有限集合$\mathcal{V}$に対し、
$$
\tilde{p}(\mathbf{w}_\mathcal{V}) = \int \tilde{p}(\mathbf{w}_\mathcal{U} | \mathbf{w}_\mathcal{S}) \tilde{p}(\mathbf{w}_\mathcal{S}) \prod_{\{\mathbf{s}_i \in \mathcal{S} \setminus \mathcal{V}\}} d(\mathbf{w}(\mathbf{s}_i)), \quad \mathcal{U} = \mathcal{V} \setminus \mathcal{S}
$$

**Kolmogorov整合性**（Appendix D）により、これは正当な空間過程$NNGP(\mathbf{0}, \tilde{\mathbf{C}}(\cdot, \cdot | \boldsymbol{\theta}))$を定義する。

**共分散関数**:
$$
\tilde{\mathbf{C}}(\mathbf{v}_1, \mathbf{v}_2; \boldsymbol{\theta}) = \begin{cases}
\tilde{\mathbf{C}}_{\mathbf{s}_i, \mathbf{s}_j}, & \text{if } \mathbf{v}_1 = \mathbf{s}_i, \mathbf{v}_2 = \mathbf{s}_j \in \mathcal{S} \\
\mathbf{B}_{\mathbf{v}_1} \tilde{\mathbf{C}}_{N(\mathbf{v}_1), \mathbf{s}_j}, & \text{if } \mathbf{v}_1 \notin \mathcal{S}, \mathbf{v}_2 = \mathbf{s}_j \in \mathcal{S} \\
\mathbf{B}_{\mathbf{v}_1} \tilde{\mathbf{C}}_{N(\mathbf{v}_1), N(\mathbf{v}_2)} \mathbf{B}_{\mathbf{v}_2}' + \delta_{(\mathbf{v}_1=\mathbf{v}_2)} \mathbf{F}_{\mathbf{v}_1}, & \text{if } \mathbf{v}_1, \mathbf{v}_2 \notin \mathcal{S}
\end{cases}
$$

### 階層モデルとMCMC

**階層モデル**（空間変化係数モデル）:
$$
\begin{aligned}
\mathbf{y}(\mathbf{t}) &\sim N(\mathbf{X}(\mathbf{t})'\boldsymbol{\beta}(\mathbf{t}), \mathbf{\Psi}) \\
\boldsymbol{\beta}(\mathbf{t}) &= \mathbf{w}(\mathbf{t}) \\
\mathbf{w} &\sim NNGP(\mathbf{0}, \tilde{\mathbf{C}}(\cdot, \cdot | \boldsymbol{\theta}))
\end{aligned}
$$

**2つのMCMCアプローチ**:

1. **逐次Gibbsサンプラー**（$\mathcal{S} \neq \mathcal{T}$）:
   - $\mathbf{w}(\mathbf{s}_i)$を1つずつサンプリング
   - 各ステップ: $\mathbf{w}(\mathbf{s}_i) | \mathbf{w}_{N(\mathbf{s}_i)}, \mathbf{y}, \boldsymbol{\theta}$
   - 計算量: $O(nm^3)$（線形）
   - 収束が遅い場合がある

2. **ブロックMetropolis-Hastings**（$\mathcal{S} = \mathcal{T}$）:
   - $\mathbf{w}_\mathcal{S}$を一括提案
   - 提案分布: $\mathbf{w}_\mathcal{S} \sim N(\mathbf{0}, \tilde{\mathbf{C}}_\mathcal{S})$
   - 計算量: $O(nm^3)$
   - 収束が速い

### 近傍集合の選択

**Vecchiaの選択**（本論文の標準）: $N(\mathbf{s}_i)$を$\{\mathbf{s}_1, \ldots, \mathbf{s}_{i-1}\}$から$m$-最近傍（ユークリッド距離）

**位置の順序付け**:
- 座標の1つ（例：$x$座標）で順序付け
- 順序付けの選択は推論結果にほぼ影響しない（Appendix C）

**近傍数$m$の選択**:
- $m$が大きいほど親GPに近づく
- シミュレーション結果：$m = 10 \sim 20$で十分な性能

### 時空間動的モデルへの拡張

**動的GPモデル**（Gelfand et al. 2005）のNNGP版：
$$
\begin{aligned}
\mathbf{y}_t &\sim N(\mathbf{X}_t'\boldsymbol{\beta}_t, \mathbf{\Psi}) \\
\boldsymbol{\beta}_t &= \mathbf{G}\boldsymbol{\beta}_{t-1} + \mathbf{w}_t \\
\mathbf{w}_t &\sim NNGP(\mathbf{0}, \tilde{\mathbf{C}}(\cdot, \cdot | \boldsymbol{\theta}))
\end{aligned}
$$

フルGPを単にNNGPで置き換えるだけで正当な動的モデルを構成できる。

## 実証結果

### 1. シミュレーション研究（Section 5.1）

**設定**:
- $n = 1000$観測点
- Matérn共分散関数（$\nu = 1.5$）
- Out-of-sample予測性能の評価

**結果**（Table 1）:
- **NNGP（$m = 10$）**: RMSPE = 1.22、計算時間 = 10.3分
- **NNGP（$m = 20$）**: RMSPE = 1.22、計算時間 = 18.1分
- **GPP（低ランク、$r = 100$）**: RMSPE = 1.68、計算時間 = 25.6分
- **Full GP**: RMSPE = 1.20、計算時間 = 58.4分

**結論**:
- NNGPはフルGPとほぼ同等の予測精度
- 計算時間は5倍以上高速
- GPP（低ランク）より精度・速度ともに優れる

**$m$の選択**（Figure 1）:
- $m = 5 \sim 20$の範囲で性能評価
- RMSPE: $m = 10$以上でほぼ一定
- 予測区間幅: $m$が大きいほどわずかに広い
- $m = 10 \sim 15$が実用的な選択

### 2. USDA森林バイオマスデータ（Section 5.2）

**データ**:
- $n = 48{,}114$地点（米国中西部）
- 応答変数: 森林バイオマス（対数変換）
- 共変量: NDVI（正規化植生指数）

**モデル比較**（Table 2）:
- **Non-spatial**: DIC = 120,847、計算時間 = 7.3時間
- **NNGP Space-varying intercept**: DIC = 119,276、計算時間 = 8.5時間
- **NNGP Space-varying coefficients**: DIC = 118,901、計算時間 = 9.2時間

**結果**:
- 空間変化切片・係数モデルが非空間モデルより適合度が高い
- 約5万点のデータでも10時間以内に推論完了
- フルGPや低ランクモデルでは実装不可能な規模

**空間パターン**（Figure 3）:
- 空間的ランダム効果が明確な空間構造を示す
- NDVIの効果が空間的に変化
- 局所的なホットスポットを適切に捉える

### 3. 順序付けの頑健性（Appendix C）

**実験**: $x$座標、$y$座標、$x+y$座標で順序付けを変更

**結果**（Table 3）:
- パラメータ推定値は順序付けによらずほぼ一致
- 信頼区間も一貫
- 空間効果の推定表面も類似（Figure 5, 6）

## 限界・残された課題

### 著者が述べる限界

1. **近傍数$m$の選択**
   - 最適な$m$はデータ依存
   - 理論的なガイドラインは限定的
   - 経験的に$m = 10 \sim 20$が多くの場合に十分

2. **非ガウスデータへの拡張**
   - 本論文はGaussian responseに焦点
   - カウント、二値データへの拡張は今後の課題

3. **時空間モデルの計算効率**
   - 動的モデルでは時間方向の依存性で計算量増加
   - Kalmanフィルタ等との統合が必要

### 本研究の観点からの限界

1. **点過程への直接適用の困難**
   - 本論文は点参照データ（geostatistical data）
   - 点過程の強度関数への適用は別途考慮が必要

2. **Presence-only未対応**
   - 完全観測データを前提
   - Thinningモデルとの統合は扱っていない

3. **組成データの未対応**
   - 多変量NNGPはあるが、組成の制約（単体上）は考慮されていない

4. **空間変化係数の識別可能性**
   - 大規模データでは識別問題が起こりうる
   - 事前分布の慎重な設定が必要

## 本研究（MMCP）との関係

### 借用する要素

1. **NNGP の計算効率化手法**
   - 疎精度行列による$O(nm^3)$の計算量
   - 大規模空間データへの対応
   - $m$-最近傍による条件付けの枠組み

2. **階層ベイズモデリングへの統合**
   - NNGPを事前分布として使用
   - 空間過程の正当性（非退化）
   - MCMCによる推論

3. **空間変化係数モデルの概念**
   - 共変量効果が空間的に変化するモデリング
   - 多変量空間過程の扱い

4. **DAGに基づく条件付け構造**
   - 有向非巡回グラフによる依存構造
   - 近傍集合の定義と選択

### MMCPによる拡張

1. **点過程への適用**
   - Dattaは点参照データ（geostatistical）
   - MMCPはCox過程の強度関数にNNGPを適用
   - Moreira et al. (2024)で示された点過程+NNGPの統合

2. **Presence-onlyデータへの対応**
   - Dattaは完全観測データ
   - MMCPはthinningモデル+NNGPの統合
   - Moreira & Gamerman (2022)の枠組みとの統合

3. **組成マークの追加**
   - Dattaは実数値応答・ランダム効果
   - MMCPは組成値マーク（単体上）
   - Multinomial logit + Pólya-Gamma + NNGP

4. **Pólya-Gammaとの統合**
   - Dattaはガウス応答
   - MMCPは離散選択（組成）とNNGPの統合

### 位置づけ

Datta et al. (2016)は：
- 大規模空間データへの厳密なGP推論を可能にする理論的基盤
- 疎精度行列による計算効率化の実現
- 空間変化係数モデルの実装可能な形式
- 低ランクモデルに対する優位性の実証

MMCPは：
- この計算効率化手法を継承
- 点過程（強度関数）への適用
- Presence-onlyデータへの対応
- 組成マークとPólya-Gammaの統合

Moreira et al. (2024)が既にNNGP + Presence-only点過程を実装しており、MMCPはそこに組成マークを追加する形。

## 引用すべき箇所

### 計算スケーラビリティの課題（Section 1）
> "Spatial process models for analyzing geostatistical data entail computations that become prohibitive as the number of spatial locations become large. [...] However, model fitting usually involves the inverse and determinant of $\mathbf{C}(\boldsymbol{\theta})$, which typically require $\sim n^3$ floating point operations (flops) and storage of the order of $n^2$. These become prohibitive when $n$ is large and $\mathbf{C}(\boldsymbol{\theta})$ has no exploitable structure."

### 低ランクモデルの限界（Section 1）
> "Furthermore, low rank models perform poorly when neighboring observations are strongly correlated and the spatial signal dominates the noise (Stein 2014). Although bias-adjusted low-rank models tend to perform better (Finley et al. 2009; Banerjee et al. 2010; Sang and Huang 2012), they increase the computational burden."

### NNGPの新規性（Section 1）
> "Our intended inferential contribution is to offer substantial scalability for fully process-based inference on underlying, perhaps completely unobserved, spatial processes. [...] While sparsity has been effectively exploited by Vecchia (1988); Stein et al. (2004) [...] for approximating expensive likelihoods cheaply, a fully process-based modeling and inferential framework has, hitherto, proven elusive. The NNGP fills this gap."

### NNGPの正当性（Section 2.2）
> "These probability densities, defined on finite topologies, conform to Kolmogorov's consistency criteria and, hence, correspond to a valid spatial process over $\mathcal{D}$ (Appendix D)."

> "Unlike low rank processes, the NNGP is not a degenerate process. It is a proper, sparsity-inducing Gaussian process, immediately available as a prior in hierarchical modeling, and, as we show in the next section, delivers massive computational benefits."

### 計算複雑度の線形性（Section 3.3）
> "The floating point operations (flops) per iteration of this algorithm is linear in the number of spatial locations, thereby rendering substantial scalability."

### 順序付けの頑健性（Section 2.1）
> "However, the aforementioned authors have cogently demonstrated that the choice of the ordering has no discernible impact on the approximation of (1) by (3). Our own simulation experiments (see Appendix C) concur with these findings; inference based upon $\tilde{p}(\mathbf{w}_\mathcal{S})$ is extremely robust to the ordering of the locations."

### 時空間モデルへの柔軟な拡張（Section 4.3）
> "Thus, one retains exactly the same structure of process-based spatial dynamic models, e.g., as in Gelfand et al. (2005), and simply replaces the independent Gaussian process priors for $\mathbf{w}_t(\mathbf{s})$ with independent NNGP's to achieve computational tractability."

> "Being a well-defined process, the NNGP ensures a valid spatial dynamic model."

### 実証結果の示唆（Section 5）
> "Using a forestry example, we show how the NNGP delivers process-based inference for spatially-varying regression models at a scale where even low-rank processes, let alone full Gaussian processes, are unimplementable even in high-performance computing environments."
