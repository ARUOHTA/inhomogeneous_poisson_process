# Moreira et al. (2024) - マーク付きPresence-only点過程とPreferential Sampling

## 基本情報
- **タイトル**: Presence-Only for Marked Point Process Under Preferential Sampling
- **著者**: Guido A. Moreira, Raquel Menezes, Laura Wise
- **ジャーナル**: Journal of Agricultural, Biological and Environmental Statistics (JABES)
- **年**: 2024（Online 2023）
- **DOI**: 確認中
- **関連論点**: 論点3（Presence-onlyデータ）、論点4（統合の困難さ）

## 主な貢献

1. **マーク付きPresence-only点過程モデルの提案**
   - Moreira & Gamerman (2022)の枠組みを連続マークに拡張
   - 各観測点に連続値（バイオマス等）が付随する場合のモデリング

2. **二重のPreferential Samplingへの対処**
   - 第1のバイアス: 点の観測自体がpreferential（presence-only）
   - 第2のバイアス: マーク値が高い場所が観測されやすい（漁師が魚が多い場所を選ぶ）

3. **NNGPによる空間依存構造の導入**
   - Gaussian Process $S(\cdot)$による空間的平滑化
   - Nearest Neighbor GP (Datta et al. 2016)で計算効率化

4. **pomppパッケージの開発**
   - RパッケージとしてCRANで公開

## 手法の概要

### モデル構造

Moreira & Gamerman (2022)を拡張した階層モデル：

$$
\begin{aligned}
X &\sim \text{IPP}(q(\cdot)p(\cdot)\lambda^*) \\
X' &\sim \text{IPP}(q(\cdot)(1-p(\cdot))\lambda^*) \\
U &\sim \text{IPP}((1-q(\cdot))\lambda^*) \\
Z(s) | s \in x \cup x' &\sim \text{logNormal}(W_z(s)\beta_z + S(s), \tau^2) \\
\text{logit } q(s) &= W_{\text{int}}(s)\beta_{\text{int}} \\
\text{logit } p(s) &= W_{\text{obs}}(s)\beta_{\text{obs}} + \gamma S(s) \\
S(\cdot) &\sim \text{NNGP}(0, \sigma^2\rho(\cdot))
\end{aligned}
$$

### 各成分の解釈

- **$X$**: 観測されたpresence locations
- **$X'$**: 観測されなかったoccurrence locations
- **$U$**: 種が存在しない場所
- **$Z(s)$**: マーク（バイオマス等）、対数正規分布
- **$S(\cdot)$**: 潜在空間過程（NNGP）
- **$\gamma$**: preferentiality parameter（マーク値と観測確率の関連）

### 二重Preferential Samplingの構造

1. **点のpreferential sampling**: $p(s)$がマークの空間構造$S(s)$に依存
2. **マークのpreferential sampling**: $\gamma > 0$なら高マーク領域が観測されやすい

### データ拡張とPólya-Gamma

Moreira & Gamerman (2022)と同様：
- $X'$と$U$の潜在過程で積分を消去
- Pólya-Gamma拡張でロジスティック回帰部分をGibbsサンプリング

### NNGPの導入

- $S(\cdot)$をNNGPでモデリング（計算量: $O(nm^3)$、$m$=近傍数）
- 条件付き分布の逐次計算で効率化
- 未観測点$x'$での$S(s)$と$Z(s)$の同時サンプリング

### 同時サンプリングの技術

未観測点での$(S(s), Z(s))$の同時サンプリング：

$$
\begin{aligned}
T &\sim \mathcal{N}\left(m\frac{\tau^2}{\tau^2+V}, \frac{\tau^2 V}{\tau^2+V}\right) \\
R &\sim \mathcal{N}\left(-m\frac{\tau^2}{\tau^2+2V}, (\tau^2+V)\frac{\tau^2}{\tau^2+2V}\right)
\end{aligned}
$$

変換: $S = T + \frac{VR}{\tau^2+V}$、$Z = \exp\{R + W_z(\cdot)\beta_z\}$

## 適用例

### ポルトガル南部イワシ漁業データ
- **データ**: 2011-2013年の1,211漁獲イベント（前処理後1,024点）
- **マーク**: イワシ漁獲量（kg）
- **共変量**: 水深（強度と観測確率の両方に使用）

### 結果
- 水深が強度に正の効果（浅い水域でイワシが多い）
- 観測確率への水深効果は弱い
- **$\gamma$が強く正**: 漁師が高バイオマス領域を選好
- 未観測occurrenceの予測: 22,500-30,500点
- 未観測バイオマス合計: 19,860-35,383トン

### シミュレーション研究
- 30データセットで正しく指定されたモデルとmisspecifiedモデルを検証
- パラメータは適切に推定
- $\gamma$（preferentiality）も正確に回復
- cloglogリンクでデータ生成してもロジットリンクで良好な結果

## 限界・残された課題

### 著者が述べる限界

1. **計算コストが高い**
   - NNGP条件付き分布の計算が各点で必要
   - MCMCの収束が遅い場合がある

2. **時間次元の未対応**
   - 漁業データでは規制変化など時間依存性があるが未考慮

3. **空間相関パラメータの固定**
   - $\theta$（$\sigma^2, \phi$）は事前に推定して固定
   - MCMCでのサンプリングは実装されていない

4. **マークの分布の制約**
   - 現在は連続正値（対数正規）のみ
   - カウント、多項などは将来の拡張

### 本研究の観点からの限界

1. **連続マークのみ対応**
   - 組成データ（$D$-部分単体上の値）は扱えない
   - multinomialリンクは提案していない

2. **単一種モデル**
   - 複数種の同時分布は範囲外

3. **計算スケーラビリティ**
   - NNGPでも大規模データでは困難

## 本研究（MMCP）との関係

### 借用する要素

1. **マーク付きPresence-only点過程の枠組み**
   - データ拡張による厳密ベイズ推論
   - thinning構造による観測バイアスの分離

2. **NNGPの活用**
   - 空間依存構造の効率的モデリング
   - 潜在空間過程$S(\cdot)$の導入

3. **Preferential Samplingの概念**
   - マーク値と観測確率の関連（$\gamma S(s)$）
   - ただしMMCPでは組成マークへの適用が必要

4. **Pólya-Gammaによる推論**
   - ロジスティック回帰部分のGibbsサンプリング

### MMCPによる拡張

1. **組成マークへの対応**
   - Moreira et al.は連続マーク（バイオマス）
   - MMCPは$D$-部分組成（黒曜石産地比率）
   - multinomial logitリンクとPólya-Gamma拡張

2. **複数「種」の同時モデリング**
   - 各黒曜石産地を「種」として扱う
   - 組成の各成分が空間的に相関

3. **空間変化係数**
   - 共変量効果が空間的に変化

### 位置づけ

Moreira et al. (2024)は、Presence-only点過程に連続マークと二重preferential samplingを導入し、NNGPで空間依存構造を効率的に扱う手法を提供した。MMCPは：

- **マークの種類**: 連続値 → 組成（単体上の値）
- **リンク関数**: logNormal → multinomial logit + Pólya-Gamma
- **preferential sampling**: $\gamma S(s)$ → 組成と観測確率の関連

という拡張を行う。特に組成マークのPólya-Gamma推論が本研究の技術的貢献となる。

## 引用すべき箇所

### 二重Preferential Samplingの定義（Section 1）
> "This work combines these concepts to handle opportunistically sampled marked point processes. Note that this implies preferentiality occurs in two ways: in the acquisition of presence-only data and in the biased collection of the marks."

### GPの追加の意義（Section 2.2）
> "The inclusion of a Gaussian Process in the intensity function extends Moreira and Gamerman (2022) presence-only model. This idea draws inspiration from the doubly stochastic process of Gonçalves and Gamerman (2018)."

### preferentialityパラメータの解釈（Section 2.2）
> "Parameter $\gamma$ measures the preferentiality of the biomass sampling procedure."

### 将来の拡張方向（Section 5）
> "Another important feature to consider is expanding the range of possible values for the marks. If the marks are not exclusively positive continuous values, a different distribution for $Z(\cdot)$ may be more appropriate. Real-valued, count, and binomial/multinomial responses are straightforward extensions."

### 計算コストの課題（Section 5）
> "The complexity of the proposed model is reflected in the computation time, which can be long due to the need to calculate the conditional distribution of each NNGP point."
