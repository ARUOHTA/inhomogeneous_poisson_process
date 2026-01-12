# Zens et al. (2023) - Ultimate Pólya Gamma Samplers: 不均衡データのための効率的MCMC

## 基本情報
- **タイトル**: Ultimate Pólya Gamma Samplers - Efficient MCMC for Possibly Imbalanced Binary and Categorical Data
- **著者**: Gregor Zens, Sylvia Frühwirth-Schnatter, Helga Wagner
- **ジャーナル**: Journal of the American Statistical Association
- **年**: 2023
- **DOI**: 10.1080/01621459.2023.2259030
- **関連論点**: 論点2（Pólya-Gamma augmentation, MCMC efficiency）

## 主な貢献
1. **新しいPólya-Gamma mixture representation**: Polson et al. (2013)のPólya-Gamma samplerを拡張し、latent variable representationを持つ新しい枠組みを提案
2. **Imbalanced Marginal Data Augmentation (iMDA)**: 不均衡データ（success probabilityが0または1に近い観測が多い）に対してロバストなMCMC boosting戦略を開発
3. **Location-based expansion**: 従来のscale-based expansion（Liu & Wu 1999）に加えて、location-based expansionを導入し、不均衡データでのmixing problemを解決
4. **汎用的な適用**: Binary、multinomial、binomial logit modelsに統一的に適用可能
5. **R package UPG**: 実装をCRANで公開

## 手法の概要

### Binary logit modelの基本

**Latent variable representation**:
$$y_i = I\{z_i > 0\}, \quad z_i = \log \lambda_i + \varepsilon_i, \quad \varepsilon_i \sim \mathcal{LO}$$

ここで$\mathcal{LO}$は logistic distribution: $f_\varepsilon(\varepsilon) = e^\varepsilon / (1 + e^\varepsilon)^2$

### 新しいPólya-Gamma mixture representation

**Key identity**（本論文の新しい貢献）:
$$f_\varepsilon(\varepsilon_i) = \frac{e^{\varepsilon_i}}{(1 + e^{\varepsilon_i})^2} = \frac{1}{4} \int e^{-\omega_i \varepsilon_i^2 / 2} p(\omega_i) d\omega_i$$

ここで$\omega_i \sim \mathcal{PG}(2, 0)$（Pólya-Gamma distribution）。

**Conditional posterior**:
$$\omega_i | \varepsilon_i \sim \mathcal{PG}(2, |\varepsilon_i|)$$

これはtilted Pólya-Gamma distributionで、効率的にサンプリング可能（Polson et al. 2013のアルゴリズムを使用）。

### Imbalanced Marginal Data Augmentation (iMDA)

**従来のMDAの限界**:
- Liu & Wu (1999)のscale-based expansion（working parameter $\delta$）は、バランスの取れたデータでは効果的
- しかし不均衡データ（例：10,000観測中2つのsuccess）では依然として mixing が遅い（Johndrow et al. 2019）

**本論文のiMDA戦略**:

#### Step 1: Location-based expansion（新規）

Working parameter $\gamma$を導入し、expanded modelを定義：
$$y_i = I\{\tilde{z}_i > \gamma\}, \quad \tilde{z}_i = \gamma + \log \lambda_i + \varepsilon_i$$

**Location-move proposal**:
1. $\tilde{\gamma} \sim \mathcal{N}(0, G_0)$をworking priorからサンプル
2. $\tilde{z}_i = z_i + \tilde{\gamma}$にshift
3. Posterior $\gamma | \omega, \tilde{z}, y \sim \mathcal{N}(g_N, G_N) I\{L(\tilde{\gamma}) \leq \gamma < U(\tilde{\gamma})\}$からサンプル   - $L(\tilde{\gamma}) = \max_{i: y_i=0} \tilde{z}_i$（$y_i=0$の最大utility）
   - $U(\tilde{\gamma}) = \min_{i: y_i=1} \tilde{z}_i$（$y_i=1$の最小utility）
4. $z_i^L = \tilde{z}_i - \gamma^{new} = z_i + \tilde{\gamma} - \gamma^{new}$に「補正」

**作用メカニズム**（Figure 2で実証）:
- Shift $\tilde{\gamma} - \gamma^{new}$がstep sizeを直接増加させる
- 特に posterior tailsでの step sizeが劇的に改善
- Systematic "push"により高posterior密度領域に誘導

#### Step 2: Scale-based expansion（既存手法の適用）

Working parameter $\delta$を導入し、expanded modelを定義：
$$y_i = I\{\tilde{z}_i > 0\}, \quad \tilde{z}_i = \sqrt{\delta} x_i \beta + \sqrt{\delta} \varepsilon_i$$

**Scale-move**:
1. $\tilde{\delta} \sim \mathcal{G}^{-1}(d_0, D_0)$をworking priorからサンプル
2. $\tilde{z}_i = \sqrt{\tilde{\delta}} z_i^L$にscale
3. Posterior $\delta | \tilde{\delta}, z^L, \omega \sim \mathcal{G}^{-1}(d_N, D_N(\tilde{\delta}))$からサンプル
4. Parameter posterior $\beta | \delta, \tilde{\delta}, z^L, \omega \sim \mathcal{N}(\sqrt{\tilde{\delta}/\delta} b_N, B_N)$からサンプル

### Ultimate Pólya-Gamma (UPG) Sampler（Algorithm 1）

**Scheme 3（完全版）**:
1. **(Z)** Latent utilities $z_i$をサンプル：$z_i | \beta, y_i \sim p(z_i | \lambda_i, y_i)$
2. **(PG)** Pólya-Gamma variables $\omega_i$をサンプル：$\omega_i | z_i, \beta \sim \mathcal{PG}(2, |z_i - x_i \beta|)$
3. **(B-L)** Location-based expansion：$z_i^L = z_i + \tilde{\gamma} - \gamma^{new}$
4. **(B-S)** Scale-based expansion and parameter update：$\beta$を$\delta, \tilde{\delta}, z^L, \omega$から条件付きGaussian posteriorでサンプル

**重要な性質**:
- 全てのステップがGibbs updateで実装可能（tuning不要）
- Conditionally Gaussianなので、state space models、random effects modelsへの拡張が容易
- 2層のlatent variablesを導入するにもかかわらず、iMDAにより autocorrelationが大幅に削減

### Multinomial logit modelへの拡張（Section 3）

**Differenced Random Utility Model (dRUM) representation**を使用：
$$y_{ij} = I\{z_{ij} > z_{ik} \text{ for all } k \neq j\}, \quad z_{ij} = \log \lambda_{ij} + \varepsilon_{ij}$$

Category $J$を reference categoryとし、$z_{ij}^* = z_{ij} - z_{iJ}$を定義すると：
$$y_{ij} = I\{z_{ij}^* > 0 \text{ for } j=1,\ldots,J-1\}$$

Logistic error differenceもlogistic distributionに従うため、同様のPólya-Gamma mixture representationとiMDAを適用。

### Binomial logit modelへの拡張（Section 4）

**新しいlatent variable representation**（Fussl et al. 2013を拡張）:
$$Y_i | N_i, \pi_i \sim \text{Binom}(N_i, \pi_i), \quad \text{logit}(\pi_i) = \log \lambda_i$$

Binomial outcomeを$N_i$個の latent binariesとして表現し、generalized logistic distributionの Pólya-Gamma mixture representationを導出。

## 実証結果

### Simulation study（Section 5、Figure 3）

**比較手法**:
- Plain DA (Scheme 1): データ拡張のみ
- PSW: Polson, Scott, Windle (2013)の original Pólya-Gamma sampler
- PX-DA (Scheme 2): Scale-based expansionのみ（Liu & Wu 1999）
- UPG (Scheme 3): Location + scale-based expansion（本論文）

**設定1: 不均衡度を変化**
- $N \in \{100, 500, 1000, 5000, 10000\}$
- Successes: 固定で2個（極めて不均衡）

**結果**（Figure 3 左上）:
- Plain DA, PSW, PX-DA: Inefficiency factorが$N$の増加とともに指数関数的に増大
- UPG: ほぼ flat（$N=10000$でも inefficiency factor < 10）

**設定2: Intercept $\beta_0$を変化**
- $N = 1000$固定
- $\beta_0 \in \{-5, -4, -3, -2, -1, 0\}$（success probabilityが変化）

**結果**（Figure 3 左下）:
- $\beta_0$が極端（-5や-4）の時、Plain DA, PSW, PX-DAは極めて非効率
- UPG: 全ての$\beta_0$でほぼ一定の効率（inefficiency factor < 50）

**Multinomial/Binomial logitでも同様の結果**（Figure 3 中央・右）

### Real data applications（Section 6）

**Application 1: Binary state space model**
- 1960-1985年のドイツの年次二値時系列データ
- Time-varyingなパラメータ$\beta_t$をstate space modelで推定
- UPG samplerにより効率的に推定可能

**Application 2: Mixture-of-experts model**
- Gating functionがmultinomial logit
- Expert functionsが異なるregression models
- iMDAなしでは convergence困難なケースでもUPGは効率的

## 限界・残された課題

### 著者が指摘する限界
1. **Pólya-Gamma samplingの計算コスト**: $\mathcal{PG}(b, c)$からのサンプリングは closed formでなく、infinite seriesの truncated sumを使用（Polson et al. 2013のアルゴリズム）。$b$や$c$が大きい場合にコスト増
2. **Tuning-free vs optimal tuning**: iMDAはtuning-freeだが、working priors $p(\gamma), p(\delta)$の選択には依然として裁量の余地。論文ではdefault values（$\mathcal{N}(0, 100)$, $\mathcal{G}^{-1}(0.01, 0.01)$）を推奨
3. **High-dimensional settings**: Covariate数が多い場合のscaling issuesについて十分な議論なし（Appendix A.1で alternative methodsとの比較で軽く言及）
4. **Spatial/temporal correlation**: 本論文はi.i.d.または conditional independenceを仮定。Spatial point process modelのようなより複雑な依存構造への適用は未検討

### 本研究の視点からの限界
1. **Point process modelへの直接適用なし**: 本論文はbinary/categorical regressionに焦点。IPPやCox processへの拡張は明示的に議論されていない
2. **NNGPとの統合**: 大規模空間データへの対応（NNGPとの組み合わせ）が未検討
3. **Marked point processへの拡張**: Marksの扱いについて議論なし
4. **Compositional marksへの適用**: Aitchison geometryとの統合は未検討

## 本研究(MMCP)との関係

### 直接的な関連
- **Pólya-Gamma augmentationの効率化**: 本研究もPólya-Gamma augmentationを使用する可能性があるため、UPG samplerの知見は直接適用可能
- **不均衡データへの対処**: 考古学データでは遺跡密度が地域によって大きく異なる（都市近郊は多数、山間部は少数）ため、iMDAの戦略が有用
- **Gibbs samplingの効率化**: MCMCの収束とmixingの改善は計算効率向上に直結

### 本研究で借用する要素
1. **Location-based expansion**: Multinomial logitでobsidian source選択をモデル化する際、特定の産地が支配的な地域（不均衡）でiMDAが有効
2. **Conditionally Gaussian framework**: State space modelや random effectsへの拡張が容易な枠組み
3. **R package UPG**: 実装の参照として有用（ただし本研究はPython実装）

### 本研究での具体的適用シナリオ

**Scenario 1: Binary presence/absence model with spatial covariates**
- $y_i = I\{\text{遺跡}i\text{が存在}\}$
- Spatial covariate: 黒曜石産地からの距離
- 不均衡データ: 産地近傍は遺跡多数、遠方は遺跡少数
- iMDAで効率的に推定

**Scenario 2: Multinomial logit for obsidian source selection**
- 遺跡$i$で発見された黒曜石のsource $j \in \{$信州, 神津島, 箱根, 高原山$\}$
- 不均衡: 特定地域では特定産地が支配的（例：関東では神津島が90%）
- iMDAなしではmixing困難、UPGで効率化

**Scenario 3: Marked IPP with Pólya-Gamma**
- 黒曜石遺跡の空間分布：IPP with intensity $\lambda(s)$
- Mark: obsidian source（multinomial）
- Pólya-Gamma augmentationでGibbs sampling
- iMDAで不均衡データに対応

### Polson et al. (2013)との比較

| 側面 | Polson et al. (2013) | Zens et al. (2023) |
|------|----------------------|--------------------|
| Representation | Marginal model（latent variable $z$なし） | Latent variable model |
| DA levels | 1層（$\omega$のみ） | 2層（$z$と$\omega$） |
| iMDA | 不可能（latent variableなし） | 可能 |
| 不均衡データ | 非効率 | 効率的 |
| Extensibility | 単純なregressionに限定 | State space、random effectsに拡張容易 |

本研究でPolson et al. (2013)を直接使う場合、不均衡データでのmixing problemが深刻。Zens et al. (2023)のUPG samplerを使うか、location/scale-based expansionの考え方を適用すべき。

### 本研究での限界として明記すべき点
第3章で引用する際、以下を明確に述べるべき：
- 本研究で不均衡データ（特定地域での遺跡の極端な多寡）が存在する場合、標準的なPólya-Gamma sampler（Polson et al. 2013）では収束が遅い可能性
- Zens et al. (2023)のiMDA戦略を適用することで改善可能だが、実装の複雑さとのトレードオフ
- 実装では、まず標準的なPólya-Gamma samplerを試し、mixing problemが観測された場合にUPGへの移行を検討

## 引用すべき箇所

### Pólya-Gamma augmentationの動機
> "While MCMC estimation based on DA is less straightforward for a logit model...Related latent variable representations with non-Gaussian errors exist for multinomial logit (MNL) models...A common solution relies on a scale-mixture representation of the non-Gaussian error distribution and introduces the corresponding scale parameters as a second level of DA." (Introduction)

Logit modelでのdata augmentationの必要性を説明する引用。

### Polson et al. (2013)の貢献
> "A seminal paper in this context is Polson, Scott, and Windle (2013) which avoids any explicit latent variable representation. They derive the Pólya-Gamma sampler that exploits a mixture representation of the non-Gaussian likelihood of the marginal model based on the Pólya-Gamma distribution and works with a single level of DA." (Introduction)

Original Pólya-Gamma samplerを紹介する引用。

### 不均衡データの問題
> "A commonly encountered challenge when working with MCMC methods based on DA is poor mixing. For binary and categorical regressions, this issue is especially pronounced for imbalanced data, where the success probability is either close to zero or one for the majority of the observations, see the excellent work of Johndrow et al. (2019). Neither the original Pólya-Gamma sampler of Polson, Scott, and Windle (2013) with a single level of DA, nor our new Pólya-Gamma sampler with two levels of DA, are an exception to this rule." (Introduction)

不均衡データでのmixing problemを説明する引用（本研究の考古学データにも該当する可能性）。

### iMDAの目的
> "To resolve this issue, we introduce imbalanced marginal data augmentation (iMDA) as a boosting strategy to make our new sampler as well as the original probit sampler of Albert and Chib (1993) robust to possibly imbalanced data." (Introduction)

iMDAの目的を簡潔に述べた引用。

### Location-based expansionのメカニズム（Figure 2の説明）
> "While step sizes increase everywhere, the improvement is particularly large in the tails of the posterior density in imbalanced datasets, where standard DA algorithms are usually highly inefficient. In addition, the shift-move evidently acts as a 'push into the right direction' that systematically leads the Markov chain back toward the highest posterior density region, effectively avoiding staying in the tails of the posterior distribution for too long." (Section 2.2)

Location-based expansionの作用メカニズムを説明する引用。

### Simulation resultのまとめ
> "The empirical inefficiency factors confirm that standard DA techniques exhibit extremely inefficient sampling behavior when confronted with imbalanced data, as shown theoretically and empirically in Johndrow et al. (2019). The MDA strategy we propose alleviates this issue and allows for rather efficient estimation also in highly imbalanced data settings." (Section 5)

Simulation結果の要約（本論文の主張を端的に述べた引用）。

### 新しいPólya-Gamma mixture representationの貢献
> "In this article, we propose a new sampling scheme involving the Pólya-Gamma distribution. Instead of working with the marginal model, we introduce a new mixture representation of the logistic distribution based on the Pólya-Gamma distribution in the latent variable representation of the logit model." (Introduction)

本論文の理論的貢献（Polson et al. 2013との違い）を明確にした引用。

### Conditionally Gaussian frameworkの利点
> "Our new Pólya-Gamma mixture representation has the advantage that the joint posterior distribution of all augmented variables is easy to sample from, as the Pólya-Gamma mixing variable follows a tilted Pólya-Gamma distribution conditional on the latent utilities. This allows to sample the unknown model parameters from a conditionally Gaussian model, facilitating posterior simulation in complex frameworks such as state space or random effects models." (Introduction)

State space modelsへの拡張可能性（本研究でも重要）を述べる引用。

### R package UPGの利用可能性
> "The underlying algorithms for probit regression and logistic regression models for binary, categorical and binomial outcomes have been made available in the R package UPG, which is available on CRAN (Zens, Frühwirth-Schnatter, and Wagner 2021)." (Introduction)

実装の参照可能性を述べる引用（Methods sectionで実装詳細を議論する際に有用）。
