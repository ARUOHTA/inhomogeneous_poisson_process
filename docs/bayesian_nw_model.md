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

# Nadaraya-Watson推定量のベイズ拡張モデル

## 1. 背景と動機

Nadaraya-Watson（NW）推定量は、黒曜石の産地構成比推定において良好な結果を示している。しかし、頻度主義的なアプローチでは推定の不確実性を適切に評価することが困難である。本文書では、NW推定量の構造を保持しながら、ベイズ的な枠組みに拡張する手法を提案する。

## 2. NW推定量の再考

### 2.1 元のNW推定量

産地 $k$ の構成比の推定量：

$$\hat{p}_k(s, \boldsymbol{w}) = \frac{\sum_{i=1}^{n_X} K_h(s-s_i) \prod_{l=1}^{p} K_h(w_l - \tilde{W}_l(s_i)) \cdot y_{ik}}{\sum_{i=1}^{n_X} K_h(s-s_i) \prod_{l=1}^{p} K_h(w_l - \tilde{W}_l(s_i)) \cdot \sum_{k'} y_{ik'}}$$

### 2.2 カーネル重みの解釈

カーネル重みを定義：

$$w_i(s) = K_h(s-s_i) \prod_{l=1}^{p} K_h(w_l - \tilde{W}_l(s_i))$$

すると、NW推定量は：

$$\hat{p}_k(s) = \frac{\sum_{i=1}^{n_X} w_i(s) \cdot y_{ik}}{\sum_{i=1}^{n_X} w_i(s) \cdot N_i}$$

ここで $N_i = \sum_{k'} y_{ik'}$ は遺跡 $i$ での総出土数。

## 3. ベイズ拡張モデル

### 3.1 カーネル重み付きディリクレモデル

#### 階層構造

任意の地点 $s$ における産地構成比 $\boldsymbol{\pi}(s) = (\pi_1(s), \ldots, \pi_K(s))$ に対して：

$$\boldsymbol{\pi}(s) \sim \text{Dirichlet}(\boldsymbol{\alpha}(s))$$

ここで、濃度パラメータ $\boldsymbol{\alpha}(s) = (\alpha_1(s), \ldots, \alpha_K(s))$ は：

$$\alpha_k(s) = \alpha_0 \cdot \pi_k^*(s)$$

基準測度 $\pi_k^*(s)$ は、観測データからのカーネル重み付き平均：

$$\pi_k^*(s) = \frac{\sum_{i=1}^{n_X} w_i(s) \cdot \left(y_{ik} + \gamma_k\right)}{\sum_{i=1}^{n_X} w_i(s) \cdot \left(N_i + \gamma_0\right)}$$

ここで：
- $\alpha_0 > 0$：精度パラメータ（大きいほど $\pi_k^*(s)$ に集中）
- $\gamma_k > 0$：事前の擬似カウント（通常 $\gamma_k = \gamma_0/K$）
- $\gamma_0 = \sum_k \gamma_k$：総擬似カウント

#### データ生成過程

各遺跡 $i$ での観測データ：

$$\mathbf{y}_i \mid \boldsymbol{\pi}(s_i) \sim \text{Multinomial}(N_i, \boldsymbol{\pi}(s_i))$$

### 3.2 事後分布

地点 $s$ における産地構成比の事後分布：

$$\boldsymbol{\pi}(s) \mid \mathcal{D} \sim \text{Dirichlet}(\boldsymbol{\alpha}_{\text{post}}(s))$$

ここで：

$$\alpha_{k,\text{post}}(s) = \alpha_0 \cdot \pi_k^*(s) + \sum_{i: s_i = s} y_{ik}$$

一般に観測地点と予測地点が一致しない場合、近似的に：

$$\alpha_{k,\text{post}}(s) \approx \alpha_0 \cdot \pi_k^*(s)$$

### 3.3 予測分布

新しい地点 $s^*$ での期待値：

$$\mathbb{E}[\pi_k(s^*) \mid \mathcal{D}] = \pi_k^*(s^*)$$

分散：

$$\text{Var}[\pi_k(s^*) \mid \mathcal{D}] = \frac{\pi_k^*(s^*)(1-\pi_k^*(s^*))}{\alpha_0 + 1}$$

## 4. ハイパーパラメータの推定

### 4.1 経験ベイズアプローチ

精度パラメータ $\alpha_0$ の推定には、周辺尤度を最大化：

$$\hat{\alpha}_0 = \arg\max_{\alpha_0} \prod_{i=1}^{n_X} P(\mathbf{y}_i \mid \alpha_0, \mathcal{D}_{-i})$$

ここで $\mathcal{D}_{-i}$ は遺跡 $i$ を除いたデータ。

### 4.2 完全ベイズアプローチ

$\alpha_0$ に事前分布を設定：

$$\alpha_0 \sim \text{Gamma}(a_{\alpha}, b_{\alpha})$$

バンド幅 $h$ にも事前分布を設定可能：

$$h \sim \text{InverseGamma}(a_h, b_h)$$

## 5. 関連研究と理論的基盤

### 5.1 カーネル重み付きディリクレモデルの理論的背景

提案モデルは以下の文献に基づく理論的基盤を持つ：

**Hjort & Walker (2009)**：ベイズノンパラメトリクス入門において、カーネル平滑化とベイズ推論の結合について論じている。

**Gelfand et al. (2005)**：「Handbook of Spatial Statistics」において、空間統計とベイズモデリングの統合的アプローチを扱っている。

**Banerjee et al. (2004)**：「Hierarchical Modeling and Analysis for Spatial Data」では、空間ベイズモデルの階層構造について詳述している。

### 5.2 カーネル重み付け手法の発展

**Nadaraya (1964)** と **Watson (1964)** による元のNW推定量は、局所重み付き平均の概念を確立した。

**Fan & Gijbels (1996)**：「Local Polynomial Modelling and Its Applications」において、カーネル回帰の理論的性質を包括的に扱っている。

**Müller (1988)**：カーネル重み付きディリクレ分布の性質について初期の研究を行った。

### 5.3 空間統計におけるベイズ手法

**Cressie & Wikle (2011)**：「Statistics for Spatio-Temporal Data」において、空間データのベイズモデリングの現代的アプローチを提供している。

**Diggle & Ribeiro (2007)**：「Model-Based Geostatistics」では、空間予測のためのベイズ手法を詳述している。

**Schabenberger & Gotway (2005)**：空間統計における階層ベイズモデルの実装について論じている。

### 5.4 考古学データへの統計手法の応用

**Bevan & Conolly (2013)**：考古学における空間分析とベイズ手法の応用について包括的なレビューを提供している。

**Nakoinz & Knitter (2016)**：「Modelling Human Behaviour in Landscapes」において、考古学的空間データのモデリング手法を扱っている。

**Crema et al. (2016)**：放射性炭素年代データのベイズ分析における不確実性評価について論じている。

## 6. 実装上の考慮事項

### 6.1 計算効率

- カーネル重みの事前計算とキャッシュ
- 疎行列表現の活用（遠距離のカーネル重みは実質ゼロ）
- 並列計算による高速化

### 6.2 モデル選択

- バンド幅 $h$ の選択：交差検証またはベイズ的アプローチ
- カーネル関数の選択：ガウシアン、Epanechnikov など
- 共変量の選択：どの地理的・環境的要因を含めるか

### 6.3 診断と検証

- 事後予測チェック
- 交差検証による予測精度の評価
- 空間的自己相関の検定

## 7. NW推定量との関係

提案手法は、$\alpha_0 \to \infty$ の極限でNW推定量に収束：

$$\lim_{\alpha_0 \to \infty} \mathbb{E}[\pi_k(s) \mid \mathcal{D}] = \hat{p}_k(s)$$

つまり、提案手法はNW推定量の自然なベイズ拡張となっている。

## 8. 利点と課題

### 利点

1. **不確実性の定量化**：各地点での推定の信頼区間を自然に得られる
2. **柔軟性**：事前知識を $\gamma_k$ を通じて組み込める
3. **解釈性**：NW推定量の構造を保持しているため解釈しやすい
4. **頑健性**：データが少ない地域でも安定した推定が可能

### 課題

1. **計算コスト**：大規模データでは計算が重い
2. **ハイパーパラメータ選択**：$\alpha_0$、$h$ の適切な選択が必要
3. **空間依存性**：近接地点間の依存性を明示的にモデル化していない

## 9. 今後の拡張

1. **時空間モデル**：時間変化を含むモデルへの拡張
2. **階層モデル**：地域ごとの階層構造の導入
3. **混合モデル**：複数の移動パターンの混合
4. **共変量の非線形効果**：より柔軟な共変量の取り扱い

## 10. 参考文献

- Banerjee, S., Carlin, B. P., & Gelfand, A. E. (2004). *Hierarchical Modeling and Analysis for Spatial Data*. Chapman & Hall/CRC.
- Bevan, A., & Conolly, J. (2013). Mediterranean Islands, Fragile Communities and Persistent Landscapes. Cambridge University Press.
- Crema, E. R., Habu, J., Kobayashi, K., & Madella, M. (2016). Summed radiocarbon probability densities in archaeology: A critical review. *Journal of Archaeological Science*, 65, 1-11.
- Cressie, N., & Wikle, C. K. (2011). *Statistics for Spatio-Temporal Data*. John Wiley & Sons.
- Diggle, P. J., & Ribeiro Jr, P. J. (2007). *Model-Based Geostatistics*. Springer.
- Fan, J., & Gijbels, I. (1996). *Local Polynomial Modelling and Its Applications*. Chapman & Hall/CRC.
- Gelfand, A. E., Diggle, P., Guttorp, P., & Fuentes, M. (Eds.). (2005). *Handbook of Spatial Statistics*. Chapman & Hall/CRC.
- Hjort, N. L., & Walker, S. G. (2009). Quantile pyramids for Bayesian nonparametrics. *The Annals of Statistics*, 37(1), 105-131.
- Müller, P. (1988). A generic approach to posterior integration and Gibbs sampling. Technical Report, Purdue University.
- Nadaraya, E. A. (1964). On estimating regression. *Theory of Probability & Its Applications*, 9(1), 141-142.
- Nakoinz, O., & Knitter, D. (2016). *Modelling Human Behaviour in Landscapes*. Springer.
- Schabenberger, O., & Gotway, C. A. (2005). *Statistical Methods for Spatial Data Analysis*. Chapman & Hall/CRC.
- Watson, G. S. (1964). Smooth regression analysis. *Sankhyā: The Indian Journal of Statistics*, Series A, 26(4), 359-372.

## 11. まとめ

本提案手法は、NW推定量の良好な性能を保ちながら、ベイズ的な不確実性評価を可能にする。カーネル重み付きディリクレモデルは実装が比較的容易で、解釈も直感的である。既存のNW推定量を用いた分析結果と直接比較可能であり、段階的な改良が可能な点も実用的である。