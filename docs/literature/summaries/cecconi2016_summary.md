# Cecconi et al. (2016) - Preferential Sampling下のベイズ地球統計

## 基本情報
- **タイトル**: Preferential sampling and Bayesian geostatistics: Statistical modeling and examples
- **著者**: Lorenzo Cecconi, Laura Grisotto, Dolores Catelan, Corrado Lagazio, Veronica Berrocal, Annibale Biggeri
- **ジャーナル**: Spatial Statistics (特集号か会議録と推測)
- **年**: 2016
- **DOI**: 記載なし（会議録の可能性）
- **関連論点**: 論点3（presence-onlyデータの扱い）、論点4（データ融合）

## 主な貢献
1. **Preferential samplingの一般的ベイズモデル**: 点過程（サンプリング位置）と測定過程が確率的に独立でない状況に対応する、共有空間成分モデルの提案
2. **3つのモデリング側面の明確化**:
   - 連続 vs 有限の空間サンプリングフレーム
   - 因果モデルと関連共変量
   - 推論目標（平均予測 vs 予測不確実性）
3. **2つの実データ応用**: 獣医寄生虫学（離散）と大気汚染曝露評価（連続）で異なる inferential goals を実証
4. **WinBUGS実装コード**: 再現可能な実装を提供

## 手法の概要

### Preferential Samplingの定義
**Diggle et al. (2010)**: 空間過程とサンプリング位置が確率的に独立でない状況

**重要な区別**:
- **確率的依存** (Stochastic dependence): Preferential samplingに該当
- **関数的依存** (Functional dependence): 共通共変量を通じた依存は厳密にはpreferential samplingでない
  - ただし、適切に共変量をモデルに含めれば区別は実質的に無意味

**Pati et al. (2011)の指摘**:
> "... accounting for informative sampling is only necessary when there is an association between the spatial surface of interest and the sampling density that cannot be explained by the [common] spatial covariates."

### 3つの過程によるモデル定式化
**Diggle et al. (2002)の枠組み**:
1. **Field process** $S$: 真の空間過程
2. **Measurement process** $Y$: 観測値
3. **Point process** $P$: サンプリング位置

**Preferential samplingの特徴**:
- $Y$ は $S$ に依存（当然）
- $P$ も $S$ に依存（通常の仮定では独立）

### 共有成分モデル (Shared Component Model)

**一般モデル** (Held et al., 2005):

$$
\begin{aligned}
X(x) &\sim PP(\lambda(x)) \\
\log(\lambda(x)) &= \alpha' + \beta' Z(x) + \eta'(x) + \delta^{-1} \xi(x) \\
Y(x) &\sim N(S(x), \tau^2) \\
S(x) &\sim GP[\mu(x), \sigma^2, \rho(d)] \\
\mu(x) &= \alpha'' + \beta'' Z(x) + \eta''(x) + \delta \xi(x)
\end{aligned}
$$

**記号**:
- $PP$: 点過程（空間的に変動する強度関数 $\lambda(x)$）
- $Z(x)$: 共変量
- $\xi(x)$: **共有潜在空間過程**（両過程の相関を誘導）
- $\eta'(x), \eta''(x)$: 過程特有の空間成分
- $\delta > 0$: チューニングパラメータ（共有成分の重要度を制御）
- $\tau^2$: ナゲット分散（測定誤差）
- $\rho(d)$: 距離 $d$ に依存する相関関数

**共有潜在過程の構造** (Besag et al., 1991):
$$
\xi(x) = u(x) + v(x)
$$
- $u(x)$: 空間的非構造ランダム項（heterogeneity）
- $v(x)$: 空間的構造ランダム項（clustering、CAR prior）

**チューニングパラメータ $\delta$ の解釈**:
- $\eta'(x) = 0$ の場合: 地球統計モデル with unmeasured covariate（強度 $\lambda(x)$ が代理変数）
- $\delta \to \infty$: Preferential sampling なし（点過程と場の過程が独立）
- $\delta$ 小: 強いpreferential sampling

### 2つの点過程の仕様

**1. 連続点過程**:
- 非均質ポアソン過程（Inhomogeneous Poisson Process）
- サンプリング位置が領域内の任意の点で発生可能
- 例: 大気汚染モニター配置

**2. 離散点過程**:
- 有限の候補地点集合が事前に定義
- 二項分布またはベルヌーイ分布でモデル化
- 例: 農場調査（農場数は有限）

### グリッドベース実装（簡略化モデル）
対象領域に細かいグリッドを重ね、グリッドセル $j$ $(j=1, \ldots, J)$ で：

$$
\begin{aligned}
X(x_j) &\sim PP(\lambda(x_j)) \\
\log(\lambda(x_j)) &= \alpha' + \beta' Z(x_j) + \eta'(x_j) + \delta^{-1}(u(x_j) + v(x_j)) \\
Y(x_{j(i)}) &\sim N(S(x_{j(i)}), \tau_Y) \\
S(x_{j(i)}) &\sim \text{SpatialExp}[\mu(x_{j(i)}); \psi(x_{j(i)}) = (f(\varphi, \gamma); \sigma)] \\
\mu(x_{j(i)}) &= \alpha'' + \beta'' Z(x_{j(i)}) + \eta''(x_{j(i)}) + \delta(u(x_{j(i)}) + v(x_{j(i)}))
\end{aligned}
$$

**SpatialExp**: 指数型共分散関数を持つガウス過程

**事前分布**:
- $\alpha', \alpha'' \sim \text{flat}()$（無情報一様分布）
- $\beta', \beta'' \sim N(0, \tau_\beta)$
- $u(x_j) \sim N(0, \tau_u)$（非構造ランダム効果）
- $v(x_j) \sim \text{CAR}(\bar{v}_{j \in i}, \tau_v)$（条件付き自己回帰）
- $\log(\delta) \sim TN(0, \tau_\delta)$（切断正規分布、$\delta$とその逆数が対称）

**予測分布**:
$$
[\tilde{\mathbf{y}} \mid \mathbf{y}; Z, u, v; \tilde{Z}, \tilde{u}, \tilde{v}] = \int [\tilde{\mathbf{y}} \mid \mathbf{y}, \Omega; ...] [\Omega \mid \mathbf{y}; ...] d\Omega
$$
$\Omega = (\alpha', \beta'', \delta; \varphi, \gamma, \sigma)$

## 実証結果

### Case Study 1: 獣医寄生虫学（カンパニア州羊農場）

**データ**:
- 総農場数: 8794
- サンプル農場数: 89（2013-2014年）
- アウトカム: 寄生虫感染（*Fasciola hepatica*）の有無
- Preferential sampling原因: 農場主の選択的報告（自主的サンプリング）

**モデル仕様**:
- グリッド: 10×10 km、184セル
- 点過程: 二項分布 $X(x_j) \sim \text{Binomial}(n(x_j), p(x_j))$
  - $n(x_j)$: セル $j$ の総農場数
  - $\text{logit}(p(x_j)) = \alpha' + \eta'(x_j) + \delta^{-1}(u(x_j) + v(x_j))$
- 共変量なし（リモートセンシング共変量は農場主の行動に影響しないと仮定）

**結果**:
- 観測感染率: 7.9%（90% CI: 3.7%–14.3%）
- 高リスク地域: 北西部、サレルノ県南部の湿地帯
- **共有成分 $\xi(x_j)$ の事後平均**（図3）:
  - "Flat"でなく、強い空間構造を示す → 選択的報告の存在を示唆
  - 感染率マップとの明確な関連

**示唆**: Preferential samplingを無視すると、感染確率の空間分布推定にバイアスが生じる

### Case Study 2: 大気汚染曝露評価（ロンバルディア州PM₁₀）

**データ**:
- グリッド: 4×4 km、1679セル
- モニター数: 58（2007年）
- アウトカム: PM₁₀年間平均濃度
- 共変量: Eulerianモデル（化学輸送モデル）の予測値
- Preferential sampling原因: モニター配置の政策的決定（都市部に集中）

**モデル仕様**:
- 点過程: 負の二項分布 $X(x_j) \sim \text{NegBin}(p(x_j), r)$
  - $p(x_j) = r/(r + \lambda(x_j))$
  - $\log(\lambda(x_j)) = \alpha' + \beta' Z(x_j) + \eta'(x_j) + \delta^{-1} v(x_j)$
- 測定過程: $\mu(x_{j(i)}) = \alpha'' + \beta'' Z(x_{j(i)}) + \eta''(x_{j(i)}) + \delta v(x_{j(i)})$
- 共変量 $Z$: Eulerianモデル出力（両過程に含む、Lee et al. 2015の"informative covariates"）

**結果**:
- **予測平均濃度**（図5）: Eulerianモデルとほぼ同一の空間分布
  - データ融合により平均表面は大きく変化せず
- **予測標準偏差**（図6, 8）: Preferential sampling調整の主要効果
  - 図8: 調整あり/なしの標準偏差比
  - モニターカバレッジが不十分な地域で**一貫して比 > 1**
  - → 調整により不確実性がより大きく評価される（正しい）
- **共有成分 $v(x_j)$ の事後平均**（図7）:
  - 共変量で説明できない残差空間変動
  - 社会政治的要因によるモニター配置の影響を反映

**主要知見**: Preferential samplingは**予測不確実性**に大きく影響
- 決定論的モデル（Eulerianモデル）では不確実性評価不可能
- 地球統計モデルで不確実性伝播を適切に評価
- モニターネットワークのカバレッジ不足地域での不確実性を過小評価しないために調整が重要

### シミュレーション研究

**設定**:
- 137グリッドセル、58観測地点（Case Study 2の設計を模倣）
- 3つのpreferential samplingシナリオ: High ($\delta = 0.4$)、Low ($\delta = 1.0$)、No ($\delta = 10000$)
- 各シナリオで100データセット生成
- MCMC: 100,000 iterations、burn-in 10,000

**比較モデル**:
- 共有モデル（提案手法）
- ベイズクリギング（preferential sampling無視）

**評価指標**（表1）:
- Bias（予測値と真値の平均絶対差）
- Variance（予測分散の平均）
- MSE（平均二乗誤差）

**結果**:

| モデル | Scenario 1 (High) | Scenario 2 (Low) | Scenario 3 (No) |
|--------|-------------------|------------------|-----------------|
| **Bayesian kriging** | | | |
| Bias | 0.129 | 0.073 | 0.052 |
| Variance | 0.0052 | 0.0047 | 0.0041 |
| MSE | 0.022 | 0.010 | 0.0068 |
| **Shared model** | | | |
| Bias | 0.030 | 0.046 | 0.047 |
| Variance | 0.0045 | 0.0042 | 0.0041 |
| MSE | 0.0054 | 0.0062 | 0.0062 |

**示唆**:
- Preferential sampling存在下で、ベイズクリギングは大きなバイアス
- 共有モデルは全シナリオでバイアス・MSEが大幅に小さい
- Scenario 3（no preferential sampling）でも僅かな改善 → 空間過程の強い構造が偶然の関連を誘導

## 限界・残された課題

### 著者が述べる限界
1. **識別可能性の問題**: 特異的空間成分 $\eta'(x), \eta''(x)$ と共有成分 $\xi(x)$ の同時推定
   - 一般モデルでは全潜在場の識別性に問題
   - 実例では特異的成分または共有成分のみを使用

2. **グリッドベース近似**: 細かいグリッドに依存
   - データが許さない状況では制限的
   - Diggle et al. (2013) のLog-Gaussian Cox過程による拡張が必要
   - 特定のパーティション（administrative units）に依存しない推論

3. **共変量選択の指針不足**:
   - どの共変量を両過程に含めるべきか、どの過程のみに含めるべきか
   - 効率性向上 vs バイアス増加のトレードオフ
   - Collider bias（傾向スコア文献のBrookhart et al. 2006参照）

4. **計算コスト**: MCMC with fine grid
   - 大規模問題では計算負荷が高い

### 追加の課題
1. **不確実性評価指標**: 標準偏差比（調整あり/なし）以外の評価方法
   - より詳細な不確実性表面の比較方法が必要

2. **時空間拡張**: Shaddick & Zidek (2015) のような時空間preferential sampling
   - 論文では空間のみ

3. **非ガウス応答**: 二項（Case Study 1）以外の非ガウス分布への拡張
   - 一般化線形混合モデルとの統合

4. **サンプリング設計最適化**: Preferential sampling下での最適設計
   - Diggle & Lophaven (2006) のベイズ地球統計設計との統合

## 本研究(MMCP)との関係

### 直接的な関係
1. **遺跡データのpreferential sampling**: 黒曜石遺跡の分布は無作為標本でない
   - 発掘調査の偏り（交通路沿い、特定地域への集中）
   - 遺跡の"発見されやすさ"と実際の分布の非独立性
   - Case Study 1の農場報告バイアスと類似構造

2. **点過程と共変量の同時モデリング**:
   - MMCPの点位置部分（IPP）は preferential samplingの枠組みと親和的
   - 黒曜石産地からの距離と遺跡分布の関係が、共有空間成分でモデル化可能

3. **データ融合と不確実性評価**:
   - トブラー距離（地形ベース）のような決定論的モデル出力との統合
   - Case Study 2の枠組みを考古学データに応用
   - 予測不確実性の定量化が重要

4. **有限vs連続サンプリングフレーム**:
   - 遺跡は離散的（有限の候補地点？）か連続的か
   - Case Study 1のように、グリッドセル内の潜在遺跡数をモデル化

### 本研究への示唆
- **IPPの拡張としてのpreferential sampling**:
  - 強度関数 $\lambda(s)$ 自体が共有空間過程に依存
  - 単純なIPPより現実的なモデル

- **共有成分による bias correction**:
  - 黒曜石組成と遺跡位置の確率的依存を明示的にモデル化
  - 「産地距離が近いところに遺跡が多い」だけでなく、
  - 「産地距離では説明できない共通の空間構造」を捉える

- **2段階モデリングへの統合**:
  - Stage 1（点位置）: Preferential samplingモデル
  - Stage 2（組成）: 点位置の共有成分を条件付けて組成をモデル化

- **チューニングパラメータ $\delta$ の推定**:
  - Preferential samplingの強さを定量化
  - $\delta$ の事後分布から、バイアスの程度を評価

### 借用すべき手法
1. **共有成分モデルの枠組み**:
   - BYM model (Besag et al., 1991): $\xi = u + v$（非構造 + 構造）
   - CAR prior for clustering component

2. **WinBUGS/JAGS実装**:
   - 論文のAppendixに完全なコード
   - MCMC実装の参考に

3. **グリッドベース近似**:
   - 細かいグリッドで空間過程を離散化
   - 計算効率とモデルの柔軟性のバランス

4. **標準偏差比による不確実性評価**:
   - Preferential sampling調整の効果を視覚化
   - 地図上で不確実性の空間分布を比較

### 拡張すべき点
1. **マーク付き点過程への統合**:
   - 点位置のpreferential sampling + 組成データ
   - 3段階階層: サンプリング強度 → 点位置 → 組成

2. **Log-Gaussian Cox過程との統合**:
   - Diggle et al. (2013) の枠組み
   - グリッド依存性を排除

3. **検出確率の明示的モデリング**:
   - Occupancy model (Dorazio 2014) との統合
   - 発見されていない遺跡の存在確率

4. **時間次元の追加**:
   - 時代ごとの遺跡分布（縄文・弥生など）
   - 時空間preferential sampling

## 引用すべき箇所

### Preferential samplingの定義（Introduction, p.1）
> "Preferential sampling refers to any situation in which the spatial process and the sampling locations are not stochastically independent."

**用途**: Preferential samplingの基本概念を説明

### 調整の必要性（Introduction, p.2, Pati et al.引用）
> "... accounting for informative sampling is only necessary when there is an association between the spatial surface of interest and the sampling density that cannot be explained by the [common] spatial covariates."

**用途**: 共変量で説明できない残差依存の重要性を強調

### バイアスの影響（Introduction, p.2）
> "Diggle et al., Pati et al. and Gelfand et al. showed that if preferential sampling is ignored, geostatistical inferences and predictions could be biased. Working in an air pollution context, Lee et al. illustrated that preferential sampling can also have a large impact on the validity and reliability of the health effect estimate."

**用途**: Preferential sampling無視のリスクを実証研究で示す

### 共有成分の役割（Methods, p.4）
> "To account for the dependence between the sampling process P and the field process S, we introduce in our model an underlying latent spatial process ξ(x) which induces a correlation between the two processes."

**用途**: 共有潜在過程の目的を明確化

### 2つの点過程仕様（Methods, p.5）
> "For the point process, P of the sampling location, we can distinguish two different possibilities: (1) P is an inhomogeneous point process whose intensity function depends on location x and whose 'events,' the sampling points, can occur at any point in the region; (2) P is a discrete point process with a finite set of eligible locations where 'events' can occur."

**用途**: 連続vs有限サンプリングフレームの使い分け

### 不確実性への影響（Case Study 2 Results, p.11）
> "However, differently from the deterministic numerical model, our Bayesian hierarchical model formulation allows us to obtain a measure of the prediction uncertainty through an appropriate propagation of uncertainty, which is not possible in the deterministic model."

**用途**: ベイズモデルによる不確実性定量化の利点

### 標準偏差への効果（Case Study 2 Results, p.12）
> "As the figure illustrates, we estimate a consistently greater standard deviation (ratio greater than 1) in the areas not properly covered by the air quality network when accounting for preferential sampling."

**用途**: Preferential sampling調整が予測不確実性に与える影響を数値で示す

### シミュレーション結果の解釈（Simulation, p.13）
> "As Table 1 indicates, Bayesian kriging behaves poorly in presence of preferential sampling, as expected... Notice that, as indicated by the results under scenario 3, our model is able to capture even small deviations from pure random sampling."

**用途**: 提案手法の優位性と感度を実証

### 残差空間成分の重要性（Discussion, p.15）
> "It is important to observe that if the systematic component of the model is correctly specified, the dependence between the sampling location process and the spatial process depends on the underlying latent residual shared spatial component."

**用途**: 共変量調整後も残る依存性の存在を強調

### 限界の認識（Discussion, p.16）
> "To approximate joint posterior densities and joint predictive densities, we used an MCMC approach working on a fine grid superimposed over our spatial domain of interest... Although advantageous and convenient in our applications, this can be limiting in situation where the data do not consent this reduction."

**用途**: グリッドベース近似の限界を正直に認める
