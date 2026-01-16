# Warton & Shepherd (2010) - Presence-onlyデータのポアソン点過程モデル

## 基本情報
- **タイトル**: Poisson point process models solve the "pseudo-absence problem" for presence-only data in ecology
- **著者**: David I. Warton, Leah C. Shepherd
- **ジャーナル**: The Annals of Applied Statistics, 4(3), 1383-1402
- **年**: 2010
- **DOI**: 10.1214/10-AOAS331
- **関連論点**: 論点3（presence-onlyデータとIPP）

## 主な貢献

Presence-onlyデータに対する**疑似欠損（pseudo-absence）問題**の解決策として**ポアソン点過程モデル**を提案。

1. **モデル仕様の改善**: 疑似欠損を人為的に生成せず、観測点 $\mathbf{y}$ のみをモデル化
2. **解釈可能性**: パラメータが疑似欠損の数・位置に依存しない
3. **実装の明確化**: 四分点（quadrature points）の選択理論を提供
4. **漸近的等価性の証明**: 疑似欠損ロジスティック回帰 → PPMに収束（$m \to \infty$）
5. **生態学的応用**: 樹種分布モデリング（Angophora costata）で実証

## 手法の概要

### Presence-onlyデータの定義

$$
\mathbf{y} = \{y_1, \ldots, y_n\}
$$

- $y_i \in \mathcal{A} \subset \mathbb{R}^2$: 種が報告された点位置
- $n$: 観測点数（ランダム）
- $\mathcal{A}$: 研究領域
- $\mathbf{x}_i = (x_{i1}, \ldots, x_{ik})$: 位置 $y_i$ での環境共変量

**特徴**:
- 欠損（absence）情報なし
- サンプリング位置・数が研究者の制御外
- 例: 博物館標本、アトラス記録、偶発的観測

### 非均質ポアソン点過程モデル

**仮定**:
1. **独立性**: 点イベント $y_1, \ldots, y_n$ は独立
2. **強度関数のlog-linear指定**:
$$
\log(\lambda_i) = \beta_0 + \sum_{j=1}^k x_{ij} \beta_j
$$

ここで $\lambda_i = \lambda(y_i)$ は単位面積あたりの期待イベント数

**対数尤度**:
$$
l(\boldsymbol{\beta}; \mathbf{y}) = \sum_{i=1}^n \log(\lambda_i) - \int_{y \in \mathcal{A}} \lambda(y) dy - \log(n!)
$$

### Berman-Turner四分近似

積分 $\int_{y \in \mathcal{A}} \lambda(y) dy$ を数値四分法で近似:
$$
\int_{y \in \mathcal{A}} \lambda(y) dy \approx \sum_{i=1}^m w_i \lambda_i
$$

**近似対数尤度（重み付きポアソン尤度）**:
$$
l_{\text{ppm}}(\boldsymbol{\beta}; \mathbf{y}, \mathbf{y}_0, \mathbf{w}) = \sum_{i=1}^m w_i \left( z_i \log(\lambda_i) - \lambda_i \right)
$$

ここで:
- $\mathbf{y}_0 = \{y_{n+1}, \ldots, y_m\}$: 四分点（quadrature points）
- $w_i$: 四分重み（点 $y_i$ 周辺の領域面積 $A_i$）
- $z_i = \frac{I(i \in \{1, \ldots, n\})}{w_i}$: 重み付き応答変数

**重要性**: GLM技術（一般化線形モデル、GAM等）が直接適用可能

### 四分点の選択

**推奨方法**: 規則的矩形グリッド
- 空間解像度を段階的に上げる
- 最大化対数尤度 $l_{\text{ppm}}(\hat{\boldsymbol{\beta}}; \mathbf{y}, \mathbf{y}_0, \mathbf{w})$ が収束するまで

**四分重みの計算**:
- タイリング法: 領域 $\mathcal{A}$ を矩形タイルに分割
- $w_i = |A_i|$ （点 $y_i$ の近傍領域の面積）
- Dirichlet tessellationも可能だが大サンプルでは実用的でない

**収束判定**:
- グリッド解像度を増やし、$l_{\text{ppm}}(\hat{\boldsymbol{\beta}})$ が安定するまで
- 実例（A. costata）: 約86,000四分点で収束

### 疑似欠損ロジスティック回帰

**従来の生態学的アプローチ**:
1. 疑似欠損 $\mathbf{y}_0 = \{y_{n+1}, \ldots, y_m\}$ をランダム生成
2. 二値応答変数 $I(i \in \{1, \ldots, n\})$ をロジスティック回帰:
$$
\log\left(\frac{p_i}{1 - p_i}\right) = \gamma_0 - \log(m - n) + \sum_{j=1}^k x_{ij} \gamma_j
$$

**問題点**:
- **モデル仕様**: 観測データ $\mathbf{y}$ でなく人工データ $\mathbf{y}_0$ も必要
- **解釈**: $p_i$ は疑似欠損数 $m$ に依存 → $m \to \infty$ で $p_i \to 0$
- **実装**: 疑似欠損の数・位置の選択が ad hoc

**対数尤度**:
$$
l_{\text{bin}}(\boldsymbol{\gamma}; \mathbf{y}, \mathbf{y}_0) = \sum_{i=1}^n \log(p_i) + \sum_{i=n+1}^m \log(1 - p_i)
$$

### 漸近的等価性の定理

**定理3.1**: $m \to \infty$ で、
$$
l_{\text{bin}}(\boldsymbol{\gamma}; \mathbf{y}, \mathbf{y}_0) = l_{\text{ppm}}(\boldsymbol{\gamma} - J_m; \mathbf{y}, \mathbf{y}_0, \mathbf{1}) + O(m^{-1})
$$

ここで $J_m = (\log m, 0, \ldots, 0)$, $\mathbf{1}$ は全要素1のベクトル

**意味**: 疑似欠損 = 四分点（重み無視）

**定理3.2**: 四分点が規則的グリッド（等重み $w_i = |\mathcal{A}|/m$）の場合、
$$
\hat{\boldsymbol{\gamma}} = \hat{\boldsymbol{\beta}} + J_{|\mathcal{A}|}
$$

Fisher情報量も等しい → **傾きパラメータと標準誤差は一致、切片のみ異なる**

**定理3.3**: 四分点がランダムサンプリングの場合、
$$
\hat{\boldsymbol{\gamma}} \xrightarrow{\mathcal{P}} \hat{\boldsymbol{\beta}} + J_{|\mathcal{A}|}
$$

確率収束（確率的に定理3.2が成立）

## 実証結果

### Angophora costata樹種分布モデリング

**データ**:
- 種: Angophora costata（オーストラリア固有樹種）
- 観測: 1972年以降の公園レンジャー報告721地点
- 領域: シドニー西方、Greater Blue Mountains World Heritage Area周辺86,000 km²
- 期間: 35年間

**環境共変量**（5変数）:
1. 最低気温（°C）
2. 最高気温（°C）
3. 年平均降水量（mm）
4. 1943年以降の火災回数
5. 「湿潤度」係数（局所水分指標）

**モデル仕様**:
- 線形項 + 二次項（全変数）
- 総パラメータ数: 16（切片 + 5線形 + 10二次 + 交互作用項）

**実装詳細**:
- 四分点: 規則的グリッド、段階的解像度上昇
- 収束時の四分点数: 約86,000点
- ソフトウェア: R package `spatstat` [Baddeley & Turner 2005]
- 四分重み計算: タイリング法（各タイルに1四分点）

**結果** (Table 2の一部):

| 変数 | PPM推定値 | PPM標準誤差 | ロジ回帰推定値 | ロジ回帰標準誤差 |
|------|----------|------------|--------------|----------------|
| 最低気温（線形） | 0.195 | 0.045 | 0.195 | 0.045 |
| 最低気温（二次） | -0.012 | 0.003 | -0.012 | 0.003 |
| 切片 | -5.234 | - | -16.891 | - |

**図1(c)**: 予測強度マップ（A. costata記録数/km²）
- 環境共変量4変数の二次関数
- 空間的に変動する分布パターンを捕捉
- 高地・低温域で強度低下

**図2(a)**: 収束診断
- PPM: $l_{\text{ppm}}(\hat{\boldsymbol{\beta}})$ は四分点数増加で収束
- ロジ回帰: $l_{\text{bin}}(\hat{\boldsymbol{\gamma}})$ は発散（$m \to \infty$ で $-\infty$）

**図2(b)**: パラメータ収束
- PPM・ロジ回帰とも傾きパラメータが同じ値に収束
- 標準誤差も一致
- 切片のみ約11.7の差（$\log|\mathcal{A}| = \log(86000) \approx 11.36$）

### 疑似欠損数の影響

**従来の生態学的推奨**（ad hoc）:
- 疑似欠損数 = 存在点数と同数（$m - n = n$）
- または $m - n = 10,000$ など固定値

**本研究の推奨**:
- 収束基準による選択（Figure 2a）
- 規則的グリッド（ランダムより効率的）
- 大規模データセット必要（$m \sim 10^5$）

**小サンプルでの問題**:
- $m$ が小さいと推定値が不安定（Figure 2b、低解像度領域）
- 疑似欠損の位置に結果が敏感 [Chefaoui & Lobo 2008引用]

## 限界・残された課題

### 著者が述べる限界

1. **独立性の仮定**:
   - 個体間相互作用でなく「報告間相互作用」が違反
   - 拡張可能（Baddeley & Turner 2005）だが実装は複雑

2. **環境共変量の完全性**:
   - 全ての関連変数がGISマップとして利用可能と仮定
   - 欠測・不完全な共変量への対応は議論していない

3. **計算コスト**:
   - 大規模四分点（$m \sim 10^5$）が必要
   - GISソフトウェアでの共変量抽出が律速
   - 空間的に不均質な四分重みの計算（Dirichlet tessellation）は大サンプルで非実用的

4. **観測過程のモデル化**:
   - 「報告されたイベント」と「実際の生物個体」の関係は暗黙
   - detection biasへの明示的対応なし

### 私たちが認識する限界

5. **空間相関の無視**:
   - 独立ポアソン点過程を仮定
   - クラスタリング、抑制パターンは扱えない
   - Cox過程（ランダム強度場）への拡張が必要

6. **時間的変動**:
   - 35年間のデータを集約（時間次元を無視）
   - 環境変動、分布域変化は考慮されていない

7. **検出確率の異質性**:
   - 全域で一様な検出確率を暗黙の仮定
   - 実際は道路近く、調査密度が高い地域で検出確率↑
   - thinned点過程としてのモデル化が理想

8. **マルチスピーシーズへの拡張**:
   - 単一種のみ
   - 複数種の共存・競争・相互作用は扱えない

9. **marked point processへの拡張**:
   - イベント位置のみ、マーク（属性情報）は考慮していない
   - 考古学データの産地マークのような構造には直接適用不可

10. **ベイズ推論の欠如**:
   - 最尤推定のみ
   - パラメータ不確実性の伝播、階層モデルへの拡張は今後の課題

## 本研究（MMCP）との関係

### 直接的な関連

1. **presence-onlyデータの理論的基礎**: MMCP考古学遺跡データも本質的にpresence-only
   - 遺跡の「欠損」は意味を持たない
   - Wartonの枠組みが直接適用可能

2. **非均質ポアソン点過程の推論**:
   - MMCP第1層: イベント位置 $\{\mathbf{s}_i\}$ ～ IPP with $\lambda(\mathbf{s})$
   - Wartonの四分近似 → GLM技術の利用
   - MCMCよりも高速な点推定が可能

3. **四分点選択の指針**:
   - MMCPの積分 $\int_\mathcal{A} \lambda(\mathbf{s}) d\mathbf{s}$ の数値計算
   - 規則的グリッド + 収束診断（Figure 2a型）
   - ad hocな選択を避ける理論的根拠

4. **疑似欠損 ≈ 背景点（background points）**:
   - MMCPの "case-control likelihood" [Fithian & Hastie 2013]
   - 背景点の数・配置が結果に影響 → Wartonの定理で正当化

### MMCPでの応用可能性

5. **トブラー距離場での強度推定**:
   - $\log\lambda(\mathbf{s}) = \mu + Y(\mathbf{s}) + \mathbf{x}(\mathbf{s})'\boldsymbol{\beta}$
   - $Y(\mathbf{s})$: 空間ランダム効果（GP）
   - $\mathbf{x}(\mathbf{s})$: 距離共変量（産地までのトブラー距離）
   - PPMで $\boldsymbol{\beta}$ を推定 → 産地選好性の定量化

6. **マーク付き点過程への拡張**:
   - Wartonは単純IPP、MMCPは産地マーク $\mathbf{w}_i$
   - 条件付き強度 $\lambda(\mathbf{s}, \mathbf{w}) = \lambda(\mathbf{s}) \cdot f(\mathbf{w} | \mathbf{s})$
   - $\lambda(\mathbf{s})$: WartonのPPMで推定
   - $f(\mathbf{w} | \mathbf{s})$: Dirichlet過程で推定

7. **計算効率化**:
   - 現行MMCP: MCMC（遅い）
   - Warton + Vecchia近似 + FRK → $O(n)$ 推論
   - 大規模遺跡データベース（数万点）への対応

### MMCPで借用すべき要素

8. **四分重み**:
   - MMCPの積分 $\int \lambda(\mathbf{s}) d\mathbf{s} = \sum w_i \lambda_i$
   - タイリング法、Voronoi tessellation
   - 重みの選択が推定精度に影響

9. **収束診断**:
   - Figure 2a型のプロット: 四分点数 vs 対数尤度
   - MMCPでも「四分点数は十分か？」の判定に有用
   - クロスバリデーションと組み合わせ

10. **GLM技術の活用**:
    - spatstat, mgcv等の既存ツール
    - MMCPのプロトタイピング・探索的分析で高速化
    - ベイズ階層モデルの初期値設定

### MMCPで拡張すべき点

11. **検出バイアスのモデル化**:
    - Wartonは検出確率を暗黙の仮定
    - MMCP: 遺跡発見バイアス $\pi(\mathbf{s})$（道路・開発密度）
    - thinned PPM: $\lambda_{\text{obs}}(\mathbf{s}) = \pi(\mathbf{s}) \lambda_{\text{true}}(\mathbf{s})$
    - $\pi(\mathbf{s})$ の推定が必須

12. **空間相関の導入**:
    - Wartonは独立IPP
    - MMCP: Log-Gaussian Cox Process
    - $\log\lambda(\mathbf{s}) = \mu + Y(\mathbf{s})$ where $Y \sim$ GP
    - Wartonの四分近似 + Cox過程 = 計算効率的LGCP

13. **時空間動態**:
    - 考古学データは時代情報あり（年代測定）
    - 時空間IPP: $\lambda(\mathbf{s}, t)$
    - Wartonの手法を時空間に拡張（文献にも言及あり）

14. **マルチスピーシーズ共存モデル**:
    - 複数産地の黒曜石が共存
    - multitype point process
    - 産地間の空間的競争・相補性

### 技術的な注意点

15. **切片の解釈**:
   - WartonのPPM切片 $\hat{\beta}_0$ ≠ ロジ回帰切片
   - $\hat{\gamma}_0 \approx \hat{\beta}_0 + \log|\mathcal{A}|$
   - MMCPで絶対強度（遺跡数/km²）を推定する際に注意

16. **モデル選択**:
   - Wartonは線形・二次項のみ
   - MMCP: トブラー距離の非線形効果、交互作用
   - AIC、BICによる選択（GLM枠組みで容易）

17. **大規模四分点のメモリ管理**:
   - $m = 86,000$ でもメモリ問題なし（2010年水準）
   - MMCPの3次元データ（標高含む）では注意
   - スパース行列、オンライン計算の活用

## 引用すべき箇所

### Presence-onlyデータの定義

> "Pearce and Boyce (2006) define presence-only data as 'consisting only of observations of the organism but with no reliable data on where the species was not found. Sources for these data include atlases, museum and herbarium records, species lists, incidental observation databases and radio-tracking studies.'" (p. 1383)

### 疑似欠損アプローチの問題点

> "We see three key weaknesses of the 'pseudo-absence' approach so widely used in ecology for analyzing presence-only data, which we describe concisely as problems of model specification, interpretation, and implementation. A sounder model specification would involve constructing a model for the observed data $\mathbf{y}$ only, rather than requiring us to generate new data $\mathbf{y}_0$ prior to constructing a model." (p. 1385)

### Berman-Turner四分近似

> "Berman and Turner (1992) showed that if the integral is estimated via numerical quadrature as $\int_{y \in \mathcal{A}} \lambda(y) dy \approx \sum_{i=1}^m w_i \lambda_i$, then the log-likelihood is (approximately) proportional to a weighted Poisson likelihood." (p. 1387)

### GLM技術の適用

> "Being able to write $l(\boldsymbol{\beta}; \mathbf{y})$ as a weighted Poisson likelihood has important practical significance because it implies that generalized linear modeling (GLM) techniques can be used for estimation and inference about $\boldsymbol{\beta}$. Further, adaptations of GLM techniques to other settings, such as generalized additive models [Hastie and Tibshirani (1990)], can then be readily applied to Poisson point process models also." (p. 1387)

### 四分点選択の推奨

> "We propose choosing quadrature points in a regular rectangular grid, and considering grids of increasing spatial resolution until the estimate of the maximized log-likelihood $l_{\text{ppm}}(\hat{\boldsymbol{\beta}}; \mathbf{y}, \mathbf{y}_0, \mathbf{w})$ has converged." (p. 1387)

### 定理3.1（漸近的等価性）

> "As $m \to \infty$, the logistic regression log-likelihood of equation (3.3) approaches the Poisson point process log-likelihood [equation (2.3)] but with all quadrature weights set to one: $l_{\text{bin}}(\boldsymbol{\gamma}; \mathbf{y}, \mathbf{y}_0) = l_{\text{ppm}}(\boldsymbol{\gamma} - J_m; \mathbf{y}, \mathbf{y}_0, \mathbf{1}) + O(m^{-1})$." (p. 1389)

### 疑似欠損 = 四分点

> "Theorem 3.1 has two interesting practical implications. First, it implies that the pseudo-absence points of presence-only logistic regression play the same role as quadrature points of a point process model, and so established guidelines on how to choose quadrature points (such as those of Section 2) can inform choice of pseudo-absences." (p. 1389)

### 定理3.2（等重みの場合）

> "That is, provided that quadrature points have been selected such that quadrature weights are equal, ignoring quadrature weights in a Poisson point process model does not change slope parameters nor their standard errors, although the intercept term will differ by $\log(|\mathcal{A}|/m)$." (p. 1390)

### 現行実践との乖離

> "Interestingly, current methods of selecting pseudo-absences in ecology [Pearce and Boyce (2006); Zarnetske, Edwards and Moisen (2007)] do not appear to be consistent with the best practice in low-dimensional numerical quadrature—points are usually selected at random rather than on a regular grid, and the number of pseudo-absences $(m - n)$ is more commonly chosen relative to the magnitude of the number of presences $(n)$ rather than based on some convergence criterion as in Figure 2(a)." (p. 1389)

### 切片の非収束

> "Note that together Theorems 3.1-3.2 suggest that when quadrature points (or, equivalently, pseudo-absences) are sampled in a regular grid at increasing resolution, the logistic regression parameter estimates and their standard errors will approach those of the point process model—with the exception of the intercept term, which diverges slowly to $-\infty$ as all $p_i \to 0$ at a rate inversely proportional to $m$. This nonconvergence of the intercept was also noticed by Owen (2007)." (p. 1390)

### 独立性の仮定

> "Note that the process being modeled here is locations where an organism has been reported rather than locations where individuals of the organism occur. Hence, the independence assumption would only be violated by interactions between records of sightings rather than by interactions between individual organisms per se." (p. 1387)

### 実践的意義

> "In particular, point process modeling offers a framework for choice of the number and location of pseudo-absences, both of which are currently chosen by ad hoc and sometimes ineffective methods in ecology, a point which we illustrate by example." (p. 1383, Abstract)
