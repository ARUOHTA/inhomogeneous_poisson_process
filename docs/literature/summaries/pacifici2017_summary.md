# Pacifici et al. (2017) - SDMsのための複数データソース統合フレームワーク

## 基本情報
- **タイトル**: Integrating multiple data sources in species distribution modeling: a framework for data fusion
- **著者**: Krishna Pacifici, Brian J. Reich, David A. W. Miller, Beth Gardner, Glenn Stauffer, Susheela Singh, Alexa McKerrow, Jaime A. Collazo
- **ジャーナル**: Ecology
- **年**: 2017
- **DOI**: 10.1002/ecy.1710
- **関連論点**: 論点4（Data fusion, Model validation）

## 主な貢献
1. **データ品質と量のトレードオフへの解決策**: 標準的な設計ベース調査データ（高品質・少数）と非標準的データ（低品質・大量、例：citizen science）を統合する柔軟なフレームワークを提示
2. **空間相関の明示的考慮**: Multivariate Conditional Autoregressive (MVCAR) modelを使用し、occurrenceとdetection errorに空間自己相関を導入（Dorazio 2014、Fithian et al. 2015が欠いていた点）
3. **3つのdata fusion approaches**: データ品質の不確実性に応じて、"Shared"（最も直接的）、"Correlation"（ロバスト）、"Covariates"（実装容易）の3つのアプローチを提案
4. **False positives/negativesの明示的モデル化**: Citizen scienceデータにおける誤同定（false positives）と不完全検出（false negatives）を同時に考慮

## 手法の概要

### 階層モデルの基本構造

**State model**（潜在的なoccupancy）:
$$\Pr(Z_i = 1 | \theta_{i0}) = q(\theta_{i0})$$

**Observation model**（データタイプ$j$ごと）:
$$Y_{ij} | Z_i, \theta_{ij} \sim f_j(y | Z_i, \theta_{ij})$$

**Random effects decomposition**:
$$\theta_{ij} = X_i^T \beta_j + \alpha_{ij}$$

### Multivariate Conditional Autoregressive (MVCAR) Model

空間相関と異なるrandom effectsの間の依存関係を同時にモデル化：
$$\alpha_i | \alpha_k \text{ for all } k \neq i \sim \text{Normal}(\rho \bar{\alpha}_i, \frac{1}{m_i} \Sigma)$$

ここで：
- $\rho \in (0, 1)$: 空間依存性の強さ
- $\bar{\alpha}_i$: 隣接領域（$m_i$個のneighbors）の平均
- $\Sigma$: $(J+1) \times (J+1)$ cross-covariance matrix

**重要な点**: $\Sigma_{jl}$は異なるrandom effects間の相関を表す。
- $\Sigma_{01} > 0$: OccupancyとBBS detection probabilityの正の相関
- $\Sigma_{02} > 0$: OccupancyとeBird abundanceの正の相関
- $\Sigma_{12} > 0$: BBS detection probabilityとeBird abundanceの正の相関

### 4つのモデルアプローチ

#### 1. Single model（ベースライン）

BBS dataのみを使用した空間occupancy model：
$$Y_{i1} | Z_i, \theta_{i1} \sim \text{Binomial}(N_i, Z_i p_i)$$

ここで$p_i = \Phi(\theta_{i1})$はdetection probability。

#### 2. Shared model（最も直接的）

両データソースがlatent stateを直接共有：

**BBS likelihood**:
$$Y_{i1} | Z_i, \theta_{i1} \sim \text{Binomial}(N_i, Z_i p_i)$$

**eBird likelihood**:
$$Y_{i2} | Z_i, \theta_{i2} \sim \text{Poisson}[E_i (Z_i \lambda_i + p_0)]$$

ここで：
- $p_i = \Phi(\theta_{i1})$: BBS detection probability
- $\lambda_i = \exp(\theta_{i2})$: eBird abundance（per unit effort）
- $E_i$: eBird effort（補助情報）
- $p_0 > 0$: False positive rate（誤同定による非存在セルでの観測確率）

**適用場面**: 両データソースが高品質で信頼できる場合。

**問題点**: 第2データソースがノイズの場合、それが第1データソースを圧倒してしまう。

#### 3. Correlation model（ロバスト）

Latent state $Z_i$をeBird likelihoodから除外：

**BBS likelihood**:
$$Y_{i1} | Z_i, \theta_{i1} \sim \text{Binomial}(N_i, Z_i p_i)$$

**eBird likelihood**（$Z_i$に条件付けない）:
$$Y_{i2} | \theta_{i2} \sim \text{Poisson}[E_i \lambda_i]$$

**情報共有メカニズム**: MVCAR cross-covariance $\Sigma$を通じて間接的に共有。
- $Y_{i2}$は$Z_i$を直接推定しない
- しかし$\theta_{i2}$（relative abundance）と$\theta_{i0}$（occupancy probability）の相関$\Sigma_{02}$を通じて間接的に影響

**適用場面**: 第2データソースの品質が不確実な場合。

**利点**: $\Sigma_{02} = \Sigma_{12} = 0$と設定すれば、データソース間のリンクを完全に切断可能（robustness）。

#### 4. Covariates model（実装容易）

第2データソースから共変量を構築し、それを$X_i$に含める：

**eBird dataから構築される共変量**:
- Log abundance: $X_{i1} = \log(\hat{A}_i + 0.1)$
- Naïve occupancy: $X_{i2} = \Phi^{-1}(0.01 + 0.98 \hat{O}_i)$
- Log effort: $X_{i3} = \log(\hat{E}_i)$

これらはGAM (Generalized Additive Model)で空間的にsmoothingされる。

**適用場面**:
- 実装が容易（既存のSDMフレームワークに組み込み可能）
- 計算負荷が低い（第2データソースの全地点を使わなくても良い）

**欠点**: 第2データソースの不確実性が伝播しない（uncertainty propagationなし）。

### Brown-headed Nuthatch (BHNU) 実証研究

**データ**:
- BBS data: 標準化された設計ベース調査（24.5 mile routeで50 stops）
- eBird data: Citizen scienceデータ（effort、観測者数などの補助情報あり）

**空間単位**: 0.25° × 0.25° lat/lon grid cells

**検証**: 2012年データで推定し、2007-2011年のBBS dataで検証。

## 実証結果

### Brown-headed Nuthatch 実データ分析

**Model comparison（Table 2）**:

| Model | MSE | Deviance |
|-------|-----|----------|
| Single | 9.28 | 3,705 |
| Covariate | 8.86 | 3,875 |
| Shared | **8.67** | **3,362** |
| Correlation | 8.64 | 3,388 |

- **全てのdata fusionアプローチがSingle modelより改善**
- SharedとCorrelationが最も低いMSEとdeviance
- SharedがdevianceでベストだがCorrelationとほぼ同等

**Spatial patterns（Fig. 2）**:
- Single model: 分布範囲を過大推定（Central Texas、Northern Arkansas、Tennesseeで非ゼロ確率）
- Data fusion models: 種の分布境界を明確に描写

**州別パフォーマンス（Table 3）**:
- 分布中心（AL、GA、LA、NC、SC）: 全モデルで類似
- 分布周縁部（AR、MD、TN）: Single modelとの乖離が大きい

**MVCAR parameters（Table 4）**:

| Model | Occupancy-Detection | Occupancy-Abundance | Spatial dependence $\rho$ |
|-------|---------------------|---------------------|---------------------------|
| Single | 0.82 (0.45, 0.93) | - | 1.000 (0.998, 1.000) |
| Shared | 0.08 (-0.36, 0.65) | 0.84 (0.49, 0.93) | 1.000 (0.999, 1.000) |
| Correlation | 0.50 (-0.02, 0.82) | **0.91 (0.75, 0.97)** | 1.000 (0.999, 1.000) |

- **強い空間依存性**（$\rho \approx 1$）
- **Occupancy-Abundance相関が高い**（SharedとCorrelationで$\Sigma_{02} > 0.8$）
- False positive rate $p_0$は極めて小さい（0.0002）

### シミュレーション研究

**設定**:
- $n = 400$ cells（20×20 grid）
- $N_i = 5$（第1データソースの観測数）、$E_i = 10$（第2データソースのeffort）
- Generating models: SharedまたはCorrelation
- 相関$\phi \in \{0.0, 0.8\}$（第2データソースの品質を表す）
- Training data: 25%または50%の第1データソース
- Detection probability: low（~0.3）またはhigh（~0.5）

**結果（Table 1 - Percent agreement）**:

**Generating model = Shared**:
| Scenario | Single | Covariate | Shared | Correlation |
|----------|--------|-----------|--------|-------------|
| $\phi=0.8$, 25% train, low det | 53.0% | 66.3% | **97.1%** | 93.2% |
| $\phi=0$, 25% train, low det | 52.8% | 63.9% | **92.6%** | 89.2% |

- Shared modelが最良（特に$\phi=0.8$で97.1%の精度）
- Correlation modelもほぼ同等
- Covariate modelはSingleより大幅に改善

**Generating model = Correlation**:
| Scenario | Single | Covariate | Shared | Correlation |
|----------|--------|-----------|--------|-------------|
| $\phi=0.8$, 25% train, low det | 52.7% | 58.5% | 54.6% | **67.9%** |
| $\phi=0$, 25% train, low det | 52.8% | 52.7% | 51.1% | 52.4% |

- **Correlation modelが最良**（$\phi=0.8$で67.9%）
- **Shared modelは$\phi=0$（情報なし）でSingleより悪化**（51.1% < 52.8%）
- Correlation modelは第2データソースが低品質でもrobust

**重要な発見**:
1. **両データソースが高品質の場合**: Shared > Correlation > Covariate > Single
2. **第2データソースが低品質の場合**: Correlation > Covariate > Single > Shared
3. **Shared modelは低品質データに脆弱**（誤った情報を過度に重視）
4. **Correlation modelはrobust**（品質不明時の安全な選択肢）

## 限界・残された課題

### 著者が指摘する限界
1. **eBird dataの補助情報への依存**: eBirdはeffort、調査時間、移動距離などの補助情報を提供するため、false positives/negativesを推定可能。しかし多くのcitizen scienceデータにはこのような情報がない
2. **パラメータ識別可能性**: 補助情報がない場合、Dorazio (2014)のようにFisher information matrixで識別可能性を確認する必要
3. **Covariates modelの不確実性伝播なし**: 第2データソースから構築した共変量の測定誤差が伝播しない
4. **Change of support**: 異なる空間スケールでのデータ統合にはchange of support priorが必要（言及のみで未実装）
5. **Species interactions**: 種間相互作用の明示的モデル化なし（共変量としての利用のみ提案）

### 本研究の視点からの限界
1. **MVCARとGPの違い**: MVCARはCAR prior（条件付き独立性）で、continuous Gaussian processとは異なる。本研究のようなcontinuous spatial domainではGPの方が自然
2. **Discretization requirement**: MVCARは空間を離散的なgrid cellsに分割する必要があり、point pattern dataへの直接適用が困難
3. **Computational efficiency**: 大規模データへのスケーラビリティの議論なし（NNGPなどの近似手法の言及なし）
4. **Marksの考慮なし**: Marked point processへの拡張なし
5. **Temporal dynamics**: 時空間動態の考慮なし（2012年のみの静的モデル）

## 本研究(MMCP)との関係

### 直接的な関連
- **データ統合の考え方**: 複数のデータソースを統合する原理（ただし本研究は単一データソース）
- **階層モデル構造**: State model（真の分布）+ Observation model（検出プロセス）
- **空間相関の重要性**: MVCARとGPの違いはあるが、空間相関を明示的にモデル化する必要性
- **False positives/negatives**: 考古学データにおける発見されていない遺跡（false negatives）と誤同定（false positives）の可能性

### 本研究で参考にできる点
1. **Robustness重視のアプローチ**: Correlation modelのように、データ品質が不確実な場合に間接的な情報共有を選択する考え方
2. **Cross-covariance structure**: 異なる過程（例：遺跡分布とobsidian産地からの距離）の相関構造を明示的にモデル化
3. **Spatial smoothing**: Covariates modelのように、GAMで空間的にsmoothingした共変量を使用する実践的アプローチ
4. **Out-of-sample validation**: 異なる時期のデータで検証する手法（本研究では時代間の検証に応用可能）

### 本研究との相違点
1. **データソース数**: 本研究は単一データソース（考古学的遺跡）のみで、標準化されたplanned surveysがない
2. **空間モデル**: MVCARではなくGaussian process（NNGP）を使用
3. **Point pattern vs grid**: 本研究はpoint pattern dataで、grid cellへの集約が不要
4. **Marks**: 本研究はobsidian sourceという離散marksを持つmarked point process
5. **Compositional data**: さらにEckardt et al. (2025)のように組成値marksへの拡張

### 本研究に適用するとしたら

**仮想的なdata fusion scenario**:
1. **Primary data source（高品質）**: 考古学的発掘調査データ（systematic surveys）
2. **Secondary data source（低品質）**: 表面採集データ（surface collections、観測バイアス大）

**Correlation modelの適用**:
- Primary data: 発掘調査での遺跡occurrenceと黒曜石組成
- Secondary data: 表面採集での黒曜石分布（検出バイアス大）
- Cross-covariance: 両者の間の空間相関を通じて情報を間接的に共有

ただし現実には本研究はsecondary data sourceを持たないため、Dorazio (2014)が指摘する「planned surveysなしでの推定」の限界を明示的に議論する必要がある。

### MVCARとGPの関係

Appendix S1で著者らはMVCARとIPPの関係を示唆している。本研究の文脈では：
- **MVCAR**: 離散的なgrid cells上のCAR prior（隣接構造ベース）
- **GP**: 連続空間上のcovariance function（距離ベース）
- **NNGP**: GPの近似で、MVCARのようにsparseな構造を持つ

本研究ではcontinuous spatial domainでpoint pattern dataを扱うため、GPの方が自然だが、MVCARの「隣接情報を借りる」考え方はNNGPの「近傍点から条件付ける」考え方と本質的に類似。

## 引用すべき箇所

### データ品質と量のトレードオフ
> "Efforts to parameterize SDMs often create a tension between the quality and quantity of data available to fit models." (Introduction)

SDMsにおけるデータ品質と量のトレードオフの問題を導入する際の引用。

### データ統合の必要性
> "Methods that allow for both data types to be used will maximize the useful information available for estimating species distributions." (Abstract)

データ統合の重要性を述べる引用。

### 空間相関の明示的考慮
> "We build on these other approaches by demonstrating a continuum of strategies that can be used to integrate data from standardized sampling protocols and non-standardized data types...we show how explicitly accounting for spatial autocorrelation in estimated parameters and estimation of detection bias (both false negatives and false positives) can be used to further deal with the challenges posed by the data types." (Introduction)

Dorazio (2014)、Fithian et al. (2015)に対する本論文の貢献（空間相関の明示的考慮）を述べる引用。

### Shared modelの脆弱性
> "This 'Shared' model is appropriate when all data sources are deemed reliable...However, this model may be problematic when the second data source is of poor quality. For example, consider the case where effort is high for the second data source but the resulting counts are purely noise due to unmodeled errors. In this extreme case the noisy second data source will overwhelm the first data source leading to poor estimates." (Methods - Shared model)

高品質データと低品質データを同等に扱うリスクを説明する引用（本研究の限界を議論する際に有用）。

### Correlation modelのrobustness
> "Unlike (3), the link between the two data sources can be severed by simply setting $\Sigma_{02} = \Sigma_{12} = 0$. Therefore, this approach should be preferred when the quality of the second data source is less certain." (Methods - Correlation model)

Correlation modelのrobustnessを説明する引用。

### シミュレーション結果のまとめ
> "In summary, the simulation study confirms that the Shared model is the most powerful when both data sources provide high-quality information about occupancy probability; the Correlation model is slightly less powerful, but more robust to contamination of the second data source (e.g., low quality of information); and the Covariates model is a simple and useful addition to the Single data source model." (Simulation results)

3つのアプローチの使い分けを簡潔にまとめた引用。

### 空間random effectsの役割
> "Modeling of spatial variation using spatial random effects was an important component of the modeling framework we present but not necessarily a requirement. However, spatial random-effects are likely to be especially useful to improve estimates by allowing information to be shared across data sets in cases where the locations of observations do not necessarily align, but are in close proximity or when improving upon standardized designs." (Discussion)

空間random effectsの役割を説明する引用（本研究でのGP/NNGPの使用を正当化する際に有用）。

### 本研究の限界として引用可能
本研究はprimary data sourceのみを使用しており、planned surveysに相当するデータがない点を限界として明記する際：

> "Ideally SDMs would be parameterized using high quality data collected under standardized design-based sampling protocols that include randomization, consistent sampling methods, and that control for observer effort and detection uncertainty. Data of this type are rare for many species and when available often lack the sampling extent necessary to produce range-wide estimates." (Introduction)

考古学データの性質（標準化されていない、観測バイアスが大きい）を議論する際の文脈引用。
