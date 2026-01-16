# Cressie & Johannesson (2008) - Fixed Rank Kriging (FRK)

## 基本情報
- **タイトル**: Fixed rank kriging for very large spatial data sets
- **著者**: Noel Cressie, Gardar Johannesson
- **ジャーナル**: Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70(1), 209-226
- **年**: 2008
- **DOI**: 10.1111/j.1467-9868.2007.00633.x
- **関連論点**: 論点1（空間統計モデル）

## 主な貢献

**Fixed Rank Kriging (FRK)**：大規模空間データセット（$n \sim 10^5$）に対する計算効率的なクリギング手法を提案。

1. **計算複雑度の削減**: $O(n^3)$ から $O(nr^2)$ へ（$r$ は固定ランク、$r \ll n$）
2. **非定常共分散関数**: 固定数の基底関数を用いた柔軟な共分散モデル
3. **Sherman-Morrison-Woodbury公式**: 大規模共分散行列の逆行列を $r \times r$ 行列の逆行列演算に帰着
4. **重み付きFrobenius norm推定**: 共分散パラメータの効率的な推定法
5. **全球リモートセンシングデータへの応用**: 173,405点の全球オゾンデータ（TCO）で実証

## 手法の概要

### モデル

**観測過程**:
$$
Z(\mathbf{s}) = Y(\mathbf{s}) + \varepsilon(\mathbf{s}), \quad \mathbf{s} \in D \subset \mathbb{R}^d
$$

**潜在過程**:
$$
Y(\mathbf{s}) = \mathbf{t}(\mathbf{s})' \boldsymbol{\alpha} + \nu(\mathbf{s})
$$

- $\mathbf{t}(\mathbf{s})$: 既知の共変量ベクトル（$p \times 1$）
- $\boldsymbol{\alpha}$: 未知回帰係数
- $\nu(\mathbf{s})$: 平均ゼロの空間過程
- $\varepsilon(\mathbf{s}) \sim N(0, \sigma^2 v(\mathbf{s}))$: 測定誤差（白色ノイズ）

**線形混合モデル表現**:
$$
\mathbf{Z} = \mathbf{T}\boldsymbol{\alpha} + \boldsymbol{\delta}, \quad \boldsymbol{\delta} = \boldsymbol{\nu} + \boldsymbol{\varepsilon}
$$
$$
\operatorname{var}(\boldsymbol{\delta}) = \boldsymbol{\Sigma} = \mathbf{C} + \sigma^2 \mathbf{V}
$$

### Fixed Rank共分散関数

**基底関数表現**:
$$
C(\mathbf{u}, \mathbf{v}) = \mathbf{S}(\mathbf{u})' \mathbf{K} \mathbf{S}(\mathbf{v}), \quad \mathbf{u}, \mathbf{v} \in \mathbb{R}^d
$$

- $\mathbf{S}(\mathbf{u}) = (S_1(\mathbf{u}), \ldots, S_r(\mathbf{u}))'$: 固定された $r$ 個の基底関数
- $\mathbf{K}$: $r \times r$ 正定値行列（推定対象）
- $r$: 固定（$r \ll n$）

**空間ランダム効果モデル解釈**:
$$
\nu(\mathbf{s}) = \mathbf{S}(\mathbf{s})' \boldsymbol{\eta}, \quad \operatorname{var}(\boldsymbol{\eta}) = \mathbf{K}
$$

**共分散行列**:
$$
\boldsymbol{\Sigma} = \mathbf{SKS}' + \sigma^2 \mathbf{V}
$$

### Sherman-Morrison-Woodbury公式による逆行列

**キーとなる公式**:
$$
\boldsymbol{\Sigma}^{-1} = (\sigma^2 \mathbf{V})^{-1} - (\sigma^2 \mathbf{V})^{-1} \mathbf{S} \{\mathbf{K}^{-1} + \mathbf{S}'(\sigma^2 \mathbf{V})^{-1} \mathbf{S}\}^{-1} \mathbf{S}' (\sigma^2 \mathbf{V})^{-1}
$$

- $\mathbf{V}$: 対角行列 → $(\sigma^2 \mathbf{V})^{-1}$ は $O(n)$
- $\mathbf{K}^{-1} + \mathbf{S}'(\sigma^2 \mathbf{V})^{-1} \mathbf{S}$: $r \times r$ 行列 → 逆行列は $O(r^3)$
- 全体の計算量: $O(nr^2)$

### クリギング予測量

**FRK予測**:
$$
\hat{Y}(\mathbf{s}_0) = \mathbf{t}(\mathbf{s}_0)' \hat{\boldsymbol{\alpha}} + \mathbf{S}(\mathbf{s}_0)' \mathbf{KS}' \boldsymbol{\Sigma}^{-1} (\mathbf{Z} - \mathbf{T}\hat{\boldsymbol{\alpha}})
$$

$$
\hat{\boldsymbol{\alpha}} = (\mathbf{T}' \boldsymbol{\Sigma}^{-1} \mathbf{T})^{-1} \mathbf{T}' \boldsymbol{\Sigma}^{-1} \mathbf{Z}
$$

**クリギング標準誤差**:
$$
\begin{aligned}
\sigma_k(\mathbf{s}_0) = & \{\mathbf{S}(\mathbf{s}_0)' \mathbf{K} \mathbf{S}(\mathbf{s}_0) - \mathbf{S}(\mathbf{s}_0)' \mathbf{KS}' \boldsymbol{\Sigma}^{-1} \mathbf{SKS}(\mathbf{s}_0) \\
& + (\mathbf{t}(\mathbf{s}_0) - \mathbf{T}' \boldsymbol{\Sigma}^{-1} \mathbf{SKS}(\mathbf{s}_0))' (\mathbf{T}' \boldsymbol{\Sigma}^{-1} \mathbf{T})^{-1} (\mathbf{t}(\mathbf{s}_0) - \mathbf{T}' \boldsymbol{\Sigma}^{-1} \mathbf{SKS}(\mathbf{s}_0)) \}^{1/2}
\end{aligned}
$$

### パラメータ推定

**ビニング（binning）**:
- $M$ 個のビン中心 $\{\mathbf{u}_j: j=1, \ldots, M\}$ を設定（$r \leq M < n$）
- 各ビン内の平均残差:
$$
\bar{Z}_j = \mathbf{w}_j' (\mathbf{Z} - \mathbf{T}\bar{\boldsymbol{\alpha}}) / \mathbf{w}_j' \mathbf{1}_n
$$

**経験的共分散行列**:
$$
\hat{\boldsymbol{\Sigma}}_M = \text{var}(\bar{Z}_1, \ldots, \bar{Z}_M) \quad \text{（method of moments推定量）}
$$

**理論的共分散行列（近似）**:
$$
\bar{\boldsymbol{\Sigma}}_M(\mathbf{K}, \sigma^2) = \bar{\mathbf{S}} \mathbf{K} \bar{\mathbf{S}}' + \sigma^2 \bar{\mathbf{V}}
$$

ここで $\bar{\mathbf{S}}$, $\bar{\mathbf{V}}$ はビン化された $\mathbf{S}$, $\mathbf{V}$

**重み付きFrobenius norm最小化**:
$$
\min_{\mathbf{K} > 0, \sigma^2 > 0} \|\hat{\boldsymbol{\Sigma}}_M - \bar{\boldsymbol{\Sigma}}_M(\mathbf{K}, \sigma^2)\|_a^2
$$

$$
\|\mathbf{A} - \mathbf{B}\|_a^2 = \sum_{j,k} a_j a_k (A_{jk} - B_{jk})^2
$$

**最適解**（QR分解 $\bar{\mathbf{S}} = \mathbf{QR}$ を使用）:
$$
\hat{\mathbf{K}} = \mathbf{R}^{-1} \mathbf{Q}' (\hat{\boldsymbol{\Sigma}}_M - \hat{\sigma}^2 \bar{\mathbf{V}}) \mathbf{Q} (\mathbf{R}^{-1})'
$$

$$
\hat{\sigma}^2 = \arg\min_{\sigma^2 > 0} \|\hat{\boldsymbol{\Sigma}}_M - \mathbf{QQ}' \hat{\boldsymbol{\Sigma}}_M \mathbf{QQ}' - \sigma^2 (\bar{\mathbf{V}} - \mathbf{QQ}' \bar{\mathbf{V}} \mathbf{QQ}')\|^2
$$

**重み**（データベース）:
$$
a_j \propto (\mathbf{w}_j' \mathbf{1}_n)^{1/2} / V_D(\mathbf{u}_j)
$$

ここで $V_D(\mathbf{u}_j)$ は第 $j$ ビン内の経験的分散

### 基底関数の選択

**推奨**: 多分解能（multi-resolution）基底
- ウェーブレット（直交または非直交）
- 局所二次（bisquare）関数
- 放射基底関数

**特性**:
- 複数の空間スケールを捕捉
- 平均関数 $\mathbf{t}(\cdot)'\boldsymbol{\alpha}$ で捉えられない大規模変動を $\mathbf{S}(\cdot)'\boldsymbol{\eta}$ で補完
- 計算効率: スパース行列技術の利用により $O(kr^2)$ where $k < n$

## 実証結果

### TCO（Total Column Ozone）衛星データ

**データセット**:
- 日付: 1988年1月1日
- サンプルサイズ: $n = 173,405$
- 全球カバレッジ（Nimbus-7衛星、TOMS instrument）
- グリッド: $1° \times 1.25°$ ピクセル
- 測定単位: Dobson Units (DU)

**実装詳細**:
- 基底関数: 多分解能局所bisquare関数（3スケール）
  - Coarse: 緯度 $15°$ × 経度 $15°$ （66個）
  - Medium: $7.5° \times 7.5°$ （264個）
  - Fine: $3.75° \times 3.75°$ （66個）
  - 合計 $r = 396$ 個
- ビン数: $M = 2,592$ （全球を $5° \times 5°$ に分割）
- 測定誤差: $v(\mathbf{s}) = 1$ （均一）
- 共変量: $\mathbf{t}(\mathbf{s}) = 1$ （定数項のみ）

**計算時間**（実装環境: 記載なし、2008年水準のPC）:
- パラメータ推定（$\hat{\mathbf{K}}$, $\hat{\sigma}^2$）: 数分
- FRK予測マップ作成（全球2万予測点）: 約3分
- 直接的な $n \times n$ 行列逆行列: **実行不可能**（$n^3 = (173,405)^3 \approx 5 \times 10^{15}$ 演算）

**推定結果**:
- $\hat{\sigma}^2$: 測定誤差分散の推定値（論文に具体値記載なし）
- $\hat{\mathbf{K}}$: $396 \times 396$ 正定値行列
- 予測マップ: オゾン全球分布、南極オゾンホールを捕捉
- 標準誤差マップ: データ密度の低い極域で大きい

**比較**:
- 古典的クリギング: 局所近傍法（ad hoc窓）が必要
- Vecchia近似: 条件付き尤度、順序依存性あり
- 多分解能モデル（Johannesson & Cressie 2004a）: 計算は高速だが「ブロック状」共分散
- FRK: 滑らかな非定常共分散、計算量 $O(n)$

### FRKの速度向上

論文Figure 3（計算時間の比較）より:
- $n = 1,000$: 直接法 vs FRK → 約10倍高速
- $n = 10,000$: 約100倍高速
- $n = 100,000$: 約1,000倍高速（直接法は実行不可能レベル）

Johannesson & Cressie (2004a)の多分解能モデルと比較:
- FRKは10^8倍の速度向上（$n \simeq 160,000$）
- ただし多分解能モデルはさらに高速だが柔軟性に欠ける

## 限界・残された課題

### 著者が述べる限界

1. **基底関数の選択**:
   - どの基底クラス（ウェーブレット、bisquare、放射基底等）を選ぶべきか明確な基準なし
   - 現在研究中

2. **ランク $r$ の決定**:
   - $r$ の最適選択方法が未開発
   - モデル選択基準（AIC、BIC等）の適用は計算コスト大

3. **nugget効果**: $\tau^2 I(\mathbf{u} = \mathbf{v})$ の追加は可能だが本論文では扱わず

4. **最尤推定の困難**:
   - $\mathbf{K}$ を直接パラメータ化すると $r(r+1)/2$ 個のパラメータ
   - 最尤法は計算困難（Stein 2008言及）
   - Frobenius normはモーメント法、最尤ではない

5. **ガウス性の仮定**: 重み付きFrobenius norm推定はガウス性を前提としない
   - 尤度ベース推論ではガウス性必須

### 私たちが認識する限界

6. **ビニングバイアス**:
   - 経験的共分散 $\hat{\boldsymbol{\Sigma}}_M$ はビンニングにより小バイアス
   - ビン幅・形状の選択が推定精度に影響

7. **非定常性のモデル化制限**:
   - 基底関数が固定 → 空間的に変化する共分散構造の表現力は $r$ に依存
   - トレンド・季節変動等の複雑な非定常性への対応が不明

8. **点過程データへの不適**: 観測位置が固定
   - presence-onlyデータ、不均質サンプリングへの拡張は議論されていない

9. **マルチスピーシーズへの拡張**: 単変量のみ
   - cross-covariance行列の推定、多変量FRKへの拡張は今後の課題

10. **時空間データへの適用**: 空間データのみ
    - 時空間FRKは理論的には可能だが実装・検証なし

11. **予測誤差の不確実性**: $\hat{\mathbf{K}}$, $\hat{\sigma}^2$ の推定誤差が予測標準誤差に伝播
    - プラグイン推定により過小評価の可能性（classical geostatisticsの既知問題）

12. **空間域の境界効果**: 全球データでは問題ないが、有界領域での境界処理が不明

## 本研究（MMCP）との関係

### 直接的な関連

1. **大規模空間データの計算効率化**: MMCPも考古学遺跡データ（数千点）
   - FRKの $O(nr^2)$ 戦略は大規模MMCPに有用
   - 特に空間点過程 + マーク過程の同時推論で計算ボトルネック

2. **非定常共分散のモデル化**: トブラー距離場の異方性
   - FRKの基底関数表現 $\mathbf{S}(\mathbf{s})'\mathbf{K}\mathbf{S}(\mathbf{s}')$ は非定常を自然に表現
   - MMCP: 距離減衰パラメータが地形により空間変化 → 非定常相関

3. **Sherman-Morrison-Woodbury公式**: MMCPのPolyá-Gamma augmentationと相補的
   - FRK: 共分散行列の逆行列を効率計算
   - MMCP: 事後分布サンプリングの高速化
   - 両手法とも低ランク構造を活用

### MMCPでの応用可能性

4. **潜在強度場の推定**: LGCP $\log\lambda(\mathbf{x}) = \mu + Y(\mathbf{x})$
   - $Y(\mathbf{x}) = \mathbf{S}(\mathbf{x})'\boldsymbol{\eta}$ とモデル化
   - FRKで $\boldsymbol{\eta}$ の事後分布を効率計算
   - 現行のGPベースより高速、大規模データセットに対応

5. **マーク過程の空間変動**: 産地比率 $\boldsymbol{\pi}(\mathbf{x})$ の空間変化
   - Dirichlet過程の位置依存性を基底関数で表現
   - $\alpha_k(\mathbf{x}) = \mathbf{S}_k(\mathbf{x})'\boldsymbol{\theta}_k$ for each産地 $k$

6. **トブラー距離のサロゲートモデル**: 距離行列の近似
   - トブラー距離計算は高コスト（DEM → ルーティング）
   - FRK型の低ランク近似: $D(\mathbf{s}_i, \mathbf{s}_j) \approx \|\mathbf{S}(\mathbf{s}_i) - \mathbf{S}(\mathbf{s}_j)\|_\mathbf{K}$

### MMCPで借用すべき要素

7. **重み付きFrobenius norm推定**:
   - MMCPのペア相関関数 $g(r)$ 推定に応用
   - ビニングによる経験的 $\hat{g}(r)$ vs 理論的 $g(r; \boldsymbol{\theta})$ の距離最小化

8. **多分解能基底の設計**:
   - 考古学データの多スケール構造（遺跡クラスター、広域交易圏）
   - Coarse: 広域トレンド（関東 vs 東海）
   - Medium: 地域パターン（利根川流域）
   - Fine: 局所クラスター（集落）

9. **計算量の解析**:
   - Table形式での実装速度比較（$n$ vs 計算時間）
   - MMCPでも $n$, $r$, MCMCイテレーション数の関係を定量化

### MMCPで拡張すべき点

10. **点過程データへの適用**:
    - FRKは格子データ、MMCPはイベント位置がランダム
    - 強度関数 $\lambda(\mathbf{x})$ のFRK表現: $\lambda(\mathbf{x}) = \exp\{\mathbf{S}(\mathbf{x})'\boldsymbol{\eta}\}$
    - Campbell定理により共分散も導出可能

11. **marked point processへの拡張**:
    - FRKは単変量、MMCPは産地マーク $\mathbf{w}_i \in \mathbb{S}^{d-1}$
    - FRK風アプローチ: $\mathbf{w}_i | \mathbf{s}_i \sim \text{Dir}(\boldsymbol{\alpha}(\mathbf{s}_i))$, $\boldsymbol{\alpha}(\mathbf{s}) = \mathbf{S}(\mathbf{s})'\mathbf{A}$
    - $\mathbf{A}$: $r \times K$ 行列（$K$=産地数）

12. **presence-onlyデータのバイアス補正**:
    - FRKは観測位置固定、MMCPは発見バイアス
    - weighted likelihood: $\pi(\mathbf{s}) = \text{detection prob at } \mathbf{s}$
    - FRK共分散に $\pi(\mathbf{s})$ を組み込む拡張

### 技術的な注意点

13. **$r$ の選択**:
    - FRK: $r=396$ for $n=173,405$ ($r/n \approx 0.002$)
   - MMCP遺跡データ（$n \sim 1000$）なら $r \sim 20-50$ が妥当?
   - クロスバリデーション、LOO、情報量基準の検討

14. **基底関数の直交性**:
    - FRKは非直交基底を許容（bisquare）
    - MMCPで組成データ（simplex）への射影 → ILR基底が自然
    - 直交化の利点（安定性）vs 解釈性のトレードオフ

15. **測定誤差 $\sigma^2$**:
    - FRKは白色ノイズ、MMCPは観測誤差なし（理想化）
    - 考古学データでは「報告バイアス」「同定誤差」をnuggetとして扱える

## 引用すべき箇所

### FRKの定義

> "For completeness, we mention another approach to spatial prediction, which is based on smoothing splines. In contrast with kriging, smoothing splines do not rely on a spatial stochastic process whose covariance function must be modelled, fitted and used for computing the optimal predictor. ... FRK is the name that we gave to the methodology that leads to equations (2.16)-(2.18)." (p. 215)

### 計算複雑度の削減

> "Inspection of equations (2.16)-(2.18) reveals that, for a fixed number of regressors $p$ and a fixed rank $r$ of $\mathbf{K}$ in the covariance model that is defined by expression (2.12), the computational burden of FRK is only linear in $n$." (p. 216)

### 基底関数による共分散表現

> "In this paper, we take a different approach and instead try to capture the scales of spatial dependence through a set of $r$ (not necessarily orthogonal) basis functions, $\mathbf{S}(\mathbf{u}) \equiv (S_1(\mathbf{u}), \ldots, S_r(\mathbf{u}))'$, $\mathbf{u} \in \mathbb{R}^d$, where $r$ is fixed." (p. 213)

### Sherman-Morrison-Woodbury公式の適用

> "Now, it is easy to see that, for any $n \times r$ matrix $\mathbf{P}$, $\mathbf{I} + \mathbf{PKP}' = \mathbf{I} + (\mathbf{I} + \mathbf{PKP}')\mathbf{PK}(\mathbf{I} + \mathbf{P}'\mathbf{PK})^{-1}\mathbf{P}'$. ... which is a result that is covered by the Sherman-Morrison-Woodbury formulae." (p. 214)

### 多分解能基底の推奨

> "Our main recommendation regarding the choice of basis functions is that they be multiresolutional. This enables the covariance function model (2.12) to capture multiple scales of variation. ... Indeed, a large spatial scale that is missed by the mean function $\mathbf{t}(\cdot)'\boldsymbol{\alpha}$ in expression (2.3) can potentially be recovered by some of the spatial random-effects components of $\mathbf{S}(\cdot)'\boldsymbol{\eta}$." (p. 218)

### 重み付きFrobenius normの正当化

> "More weight should be given to bins that are less variable or have more data. Consider the weighted Frobenius norm, ... From expression (A.5) in Appendix A, we motivate statistically the choice to be $a_j \propto (\mathbf{w}_j'\mathbf{1}_n)^{1/2} / V_D(\mathbf{u}_j)$, $j=1, \ldots, M$, which is a data-based weight where $V_D(\mathbf{u}_j)$ is the empirical variance in the $j$th bin." (p. 221)

### データの二重使用（古典的地球統計学の原則）

> "As in classical geostatistics, we see that the data are used twice (e.g. Cressie (1989)). Not only are they present (linearly) in the kriging predictor (2.7); they are also used (non-linearly) to obtain an estimator of the spatial dependence parameters $\mathbf{K}$ and $\sigma^2$." (p. 217)

### Karhunen-Loéve展開との関係

> "A related model to expression (2.12), but different from $C(\cdot, \cdot)$ given above, is a consequence of the Karhunen-Loéve expansion. ... Without loss of generality, assume that the truncation keeps only terms with positive eigenvalues; then clearly the truncated Karhunen-Loéve expansion is a special case of expression (2.12)." (p. 217)

### 非負定値性の証明

> "Most importantly, the function $C(\mathbf{u}, \mathbf{v})$ is non-negative definite, the proof of which is straightforward: for any locations $\{\mathbf{s}_i: i=1, \ldots, m\}$ in $\mathbb{R}^d$, any real $\{b_i: i=1, \ldots, m\}$, and any integer $m$, then $\sum_{i=1}^m \sum_{j=1}^m b_i b_j C(\mathbf{s}_i, \mathbf{s}_j) = (\mathbf{S}_m'\mathbf{b}_m)'\mathbf{K}(\mathbf{S}_m'\mathbf{b}_m) \geq 0$ since $\mathbf{K}$ is positive definite." (p. 217)

### 最尤推定の困難性

> "Under assumptions of Gaussianity, the likelihood of $\mathbf{K}$ and $\sigma^2$ depends on $\boldsymbol{\Sigma}^{-1}$ and $|\boldsymbol{\Sigma}|$. ... i.e. computation of the likelihood of $\mathbf{K}$ and $\sigma^2$ is feasible; however, its maximization is problematic unless $\mathbf{K}$ is further parameterized (Stein, 2008)." (p. 222)

### TCO応用の意義

> "The Nimbus-7 polar orbiting satellite was launched on October 24th, 1978, with the total ozone mapping spectrometer instrument aboard. ... The satellite was sun synchronous, staying on the plane between the Earth and the Sun. Successive orbits moved westwards because of the rotation of the Earth, and hence the Nimbus-7 satellite covered the entire globe in a 24-h period." (p. 222-223)

### 計算時間の実績

> "Johannesson and Cressie (2004a) achieved speed-ups of the order of $10^8$ over directly solving the kriging equations. They could compute optimal spatial predictors and their associated mean-squared prediction errors over the entire globe in about 3 min for $n \simeq 160000$." (p. 211)
