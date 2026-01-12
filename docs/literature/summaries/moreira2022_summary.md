# Moreira & Gamerman (2022) - Presence-onlyデータの厳密ベイズ推論

## 基本情報
- **タイトル**: Analysis of Presence-only Data via Exact Bayes, with Model and Effects Identification
- **著者**: Guido A. Moreira, Dani Gamerman
- **ジャーナル**: Annals of Applied Statistics (AOAS)
- **年**: 2022
- **DOI**: 10.1214/21-AOAS1569
- **関連論点**: 論点3（Presence-onlyデータ）

## 主な貢献

1. **厳密なベイズ推論の実現**
   - データ拡張（data augmentation）により、尤度関数の厄介な積分を回避
   - 従来手法の近似（格子離散化等）に頼らない厳密な推論

2. **識別可能性問題への対処**
   - 従来のthinned IPPモデルでは強度と観測確率の係数が非識別
   - 本手法では同じ共変量を両方に含めても識別可能

3. **潜在過程$X'$による予測**
   - 未観測のoccurrence過程$X'$の事後分布から直接予測
   - 観測バイアスを除いた真の種分布を推定

4. **bayesPOパッケージの開発**
   - Rパッケージとして実装を公開

## 手法の概要

### モデル構造

3つの条件付き独立な非均質ポアソン過程による拡張：

$$
\begin{aligned}
X &\sim \text{IPP}(q(\cdot)p(\cdot)\lambda^*) & \text{（観測されたpresence）}\\
X' &\sim \text{IPP}(q(\cdot)(1-p(\cdot))\lambda^*) & \text{（未観測のoccurrence）}\\
U &\sim \text{IPP}((1-q(\cdot))\lambda^*) & \text{（種が存在しない場所）}
\end{aligned}
$$

ここで：
- $q(s) = \frac{e^{Z(s)\beta}}{1+e^{Z(s)\beta}}$：強度関連関数（ロジスティック）
- $p(s) = \frac{e^{W(s)\delta}}{1+e^{W(s)\delta}}$：観測確率（thinning）
- $\lambda^*$：強度の上限

### 尤度の厳密計算

3過程の強度の和が$\lambda^*$（定数）となるため、積分が消える：

$$
q(s)p(s) + q(s)(1-p(s)) + (1-q(s)) = 1 \quad \forall s \in \mathcal{D}
$$

拡張尤度：

$$
L_x(q,p,\lambda^*,x',u) = \frac{e^{-\lambda^*|\mathcal{D}|}}{n_{x'}!n_u!} \lambda^{*n_x+n_{x'}+n_u} \prod_{s \in x} q(s)p(s) \prod_{s \in x'} q(s)(1-p(s)) \prod_{s \in u}(1-q(s))
$$

### Gibbsサンプラー

Pólya-Gamma augmentation（Polson et al. 2013）を用いた完全Gibbsサンプラー：

1. **$X'$と$U$のサンプリング**: HPP($\lambda^*|\mathcal{D}|$)から点を生成し、$q(s)$と$p(s)$でフィルタリング
2. **$\lambda^*$のサンプリング**: Gamma事後分布
3. **$\beta$のサンプリング**: Pólya-Gamma拡張によるロジスティック回帰（$x \cup x'$が成功、$u$が失敗）
4. **$\delta$のサンプリング**: 同様（$x$が成功、$x'$が失敗）

### 従来モデルとの関係

**命題1**: $\lambda^* \to \infty$かつ$\beta_0 \to -\infty$の極限で、提案モデルは従来の対数線形IPPに収束

$$
\lim_{\lambda^* \to \infty, \beta_0 \to -\infty} q(\cdot)\lambda^* = e^{Z(\cdot)\beta}
$$

## 適用例

### 1. ユーカリ（Eucalyptus Sparsifolia）データ
- **データ**: オーストラリア・ブルーマウンテン地域の230 presence locations
- **共変量**: 火災回数、気温、降水量、土壌タイプ（強度）; 道路・都市からの距離（観測確率）
- **結果**:
  - 従来モデル（AUC≈0.587）より提案モデル（AUC≈0.618）が優れた予測性能
  - 観測バイアス（道路・都市への近接性）を適切に分離

### 2. アンジコ（Anadenanthera Colubrina）データ
- **データ**: ブラジル半乾燥地域の140 presence locations
- **目的**: 同じ共変量を強度と観測確率に含めた場合の識別可能性検証
- **結果**: $\beta$と$\delta$の事後分布がほぼ独立（識別可能）

## 限界・残された課題

### 著者が述べる限界

1. **$\lambda^*$の識別可能性**
   - $\lambda^*$と$q(\cdot)$（特に切片$\beta_0$）の間に識別問題
   - 事前情報の導入で緩和可能

2. **計算時間の不確実性**
   - $X'$と$U$のサイズがランダムなため、各MCMC反復の計算量が変動
   - $\lambda^*$が大きいほど計算時間増加

3. **空間的残差相関の未対応**
   - 共変量で説明できない空間構造（Gaussian Process等）は本論文では扱わない
   - Renner et al. (2015)らが指摘する拡張の必要性

### 本研究の観点からの限界

1. **単一種モデル**
   - 複数種の同時モデリングは扱っていない

2. **組成データの未対応**
   - マーク（組成）を持つ点過程は範囲外

3. **空間変化係数の未導入**
   - 係数は空間的に一定と仮定

## 本研究（MMCP）との関係

### 借用する要素

1. **データ拡張による厳密ベイズ推論**
   - 潜在過程$X'$と$U$の導入による積分消去
   - thinningモデルの数学的構造

2. **Pólya-Gamma augmentationの適用**
   - ロジスティックリンクのGibbsサンプリング
   - 観測確率モデリングへの応用

3. **識別可能性の議論**
   - 強度と観測確率の分離可能性
   - パラメトリック形式の重要性

### MMCPによる拡張

1. **組成値マークの追加**
   - Moreira & Gamermanは点の位置のみ
   - MMCPは各点に組成（黒曜石産地比率）を付与

2. **空間的依存性の導入**
   - Gaussian Process / NNGPによる空間相関
   - 空間変化係数モデル

3. **multinomial logitとの統合**
   - 組成のモデリングにPólya-Gamma拡張を適用
   - 強度・観測確率・組成の同時推論

### 位置づけ

Moreira & Gamerman (2022)は、presence-onlyデータの厳密なベイズ推論フレームワークを確立し、データ拡張とPólya-Gammaによる効率的なMCMCを提供した。MMCPはこのフレームワークを基盤として：
- 組成値マークを追加
- 空間的依存構造を導入
- 複数種（黒曜石産地）の同時モデリング

を実現する。

## 引用すべき箇所

### データ拡張の動機（Section 3）
> "Our proposal uses the idea from Adams, Murray and MacKay (2009) and Gonçalves and Gamerman (2018). [...] In essence, the intensity function $\lambda(\cdot)$ is assumed to be upper bounded and is parameterized so that $\lambda(s) = q(s)\lambda^*$, where $q: \mathcal{D} \to (0,1)$."

### 厳密推論の利点（Section 1）
> "An additional contrast with the traditional model, as hinted above, is the use of exact inference in the sense that no approximation is applied to the model. This is achieved through a data augmentation technique with latent point processes."

### 識別可能性の改善（Section 4）
> "For our proposal, however, the set of covariates $Z(\cdot)$ relates to the intensity function through $q(\cdot)$ which must vary between 0 and 1. [...] Thus, the method usually fits an unthinned Poisson process, adding the observability covariates in the intensity set."

### 予測の解釈（Section 5.1）
> "The component that represents those very occurrences is $X'$, defined in equation (4). This process is sampled in the MCMC procedure, and these realizations can be seen as a sample from the predictive distribution of the infinite-dimensional process of unobserved occurrences."

### 将来の方向性（Section 8）
> "Currently, presence-only methodology development has somewhat stabilized. Instead, multiple-species modeling [...] and joining presence-only with presence-absence data [...] have gained attention."
