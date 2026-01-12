# Phillips et al. (2006) - Maxent: 最大エントロピー種分布モデル

## 基本情報
- **タイトル**: Maximum entropy modeling of species geographic distributions
- **著者**: Steven J. Phillips, Robert P. Anderson, Robert E. Schapire
- **ジャーナル**: Ecological Modelling
- **年**: 2006
- **DOI**: 10.1016/j.ecolmodel.2005.03.026
- **関連論点**: 論点3（presence-onlyデータの扱い）

## 主な貢献
1. **Maxent手法の種分布モデリングへの導入**: 最大エントロピー原理をpresence-onlyデータからの種分布推定に初めて体系的に適用
2. **明確な数学的定式化**: 不完全情報からの確率分布推定という統計力学の原理を、生態学的ニッチモデリングへ翻訳
3. **GARPとの包括的比較**: 新熱帯区の2哺乳類種（低地種と高地種）で、Maxentが既存手法（GARP）を上回ることを実証
4. **オープンソース実装**: Javaで実装されたGUIツール（http://www.cs.princeton.edu/~schapire/maxent）を無料公開

## 手法の概要

### 最大エントロピー原理（Jaynes, 1957）
**基本原則**: 未知の確率分布を近似する際、既知の制約を満たしつつ、最もエントロピーの高い（＝最も広がった、一様分布に最も近い）確率分布を選ぶべき

**エントロピーの定義**:
$$
H(\hat{\pi}) = -\sum_{x \in X} \hat{\pi}(x) \ln \hat{\pi}(x)
$$

> "The maximum entropy principle can be interpreted as saying that no unfounded constraints should be placed on $\hat{\pi}$, or alternatively, [it] agrees with everything that is known, but carefully avoids assuming anything that is not known." (Jaynes, 1990)

### 機械学習的視点
**特徴（Features）**: 環境変数またはその関数 $f_1, \ldots, f_n$

**制約**: 各特徴の経験的平均（presence localitiesでの平均）を再現
$$
\hat{\pi}[f_j] = \tilde{\pi}[f_j], \quad \text{for each feature } f_j
$$

**Gibbs分布**: Maxent分布は以下の形式を取る
$$
q_\lambda(x) = \frac{e^{\lambda \cdot f(x)}}{Z_\lambda}
$$

- $\lambda$: 特徴重み（feature weights）のベクトル
- $Z_\lambda$: 正規化定数（partition function）
- Maxent分布 = サンプル点の尤度を最大化するGibbs分布

**ℓ1正則化**: 過学習を防ぐため、緩和された制約を使用
$$
|\hat{\pi}[f_j] - \tilde{\pi}[f_j]| \leq \beta_j
$$

この場合、Maxent分布は以下を最小化：
$$
\tilde{\pi}[-\ln(q_\lambda)] + \sum_j \beta_j |\lambda_j|
$$

- 第1項: Log loss（負の対数尤度）
- 第2項: パラメータの大きさへのペナルティ（Lasso in GLM context）

### 5つの特徴タイプ
1. **Linear features**: 連続変数 $f$ そのもの
   - 制約: $\hat{\pi}[f]$ が観測平均に近い

2. **Quadratic features**: $f^2$
   - Linear featureと併用で、分散 $\hat{\pi}[f^2] - \hat{\pi}[f]^2$ を制約
   - 種の許容範囲（tolerance）をモデル化

3. **Product features**: $f \cdot g$（2変数の積）
   - 共分散 $\hat{\pi}[fg] - \hat{\pi}[f]\hat{\pi}[g]$ を制約
   - 変数間相互作用を取り込む

4. **Threshold features**: $f$ が閾値以上なら1、それ以外0
   - 任意の応答曲線を閾値関数の線形結合で近似可能
   - 本研究では未使用

5. **Binary features**: カテゴリ変数用（k個のカテゴリに対してk個の2値特徴）
   - 各カテゴリの出現割合を制約

### 種分布モデリングへの応用
**データモデル**: Presence localities $x_1, \ldots, x_m$ を未知分布 $\pi$ からのサンプルと解釈

**理想的サンプリング**: ランダムにピクセルを選び、種が存在すれば記録
- この場合、$\pi = p(x \mid y=1)$（種が存在する条件付き分布）
- Bayesの定理より、$\pi \propto p(y=1 \mid x)$（出現確率に比例）

**現実的解釈**: サンプリングバイアスが存在するため、$\pi$ は**相対的な環境適合度指標**（relative index of environmental suitability）と解釈

### GLM/GAMとの関係
**類似性**:
- GLM (Gaussian logit model): $\text{logit}(p) = \alpha + \beta_1 f_1(x) + \gamma_1 f_1(x)^2 + \ldots$
- Maxent (linear + quadratic features): $\log(p) = \lambda_1 f_1(x) + \lambda_2 f_1(x)^2 + \ldots$
- GAM: $\text{logit}(p) = g_1(f_1(x)) + \ldots + g_n(f_n(x))$
- Maxent (threshold features): 同様の柔軟な応答曲線

**相違点**:
1. **データ要求**: GLM/GAMはabsenceデータ（または pseudo-absences）が必要
   - Maxent: Presenceのみで確率分布をモデル化、absenceの解釈不要
2. **Generative vs Discriminative**:
   - Maxent: Generative（$p(x \mid y=1)$ をモデル化）
   - GLM/GAM: Discriminative（$p(y \mid x)$ を直接モデル化）
   - 訓練データが少ない場合、generativeが有利（Ng & Jordan, 2001）

### 実装
**Sequential-update algorithm** (Dudík et al., 2004):
- 反復的に1つの重み $\lambda_j$ を選択し、正則化log lossを最小化するよう調整
- 決定論的、大域最適への収束保証

**Cumulative representation**: デフォルト出力
- ピクセル $x$ の値 = $x$ 以下の確率を持つ全ピクセルの確率和 × 100（%）
- 閾値 $t$ で二値化 → 約 $t\%$ の test localities omission

**ユーザーパラメータ**（デフォルト値）:
- 収束閾値: $10^{-5}$
- 最大反復数: 1000
- 正則化値 $\beta$: $10^{-4}$
- 使用特徴: Linear, Quadratic, Product, Binary

## 実証結果

### 対象種
1. **Bradypus variegatus**（ミユビナマケモノ）:
   - 低地種、ホンジュラス〜北アルゼンチン
   - 116 presence localities（博物館標本）

2. **Microryzomys minutus**（小型山地ネズミ）:
   - 高地種、アンデス1000-4000m
   - 88 presence localities

### 環境変数（0.05° × 0.05° 解像度）
- **気候**: IPCC 12変数（年間雲量、気温日較差、霜頻度、蒸気圧、降水量など）
- **標高**: USGS HYDRO1k
- **潜在植生**: 15主要生息地タイプ（Dinerstein et al., 1995）

### 評価方法
**10回のランダム分割**: 各種で10通りの train/test 分割

**閾値依存テスト**:
- Binomial omission test: 両手法とも全分割で有意（$p < 0.001$）
- Equalized predicted area test:
  - *B. variegatus*: 有意差なし
  - *M. minutus*: Maxentが有意に低いomission率（$p = 0.036, 0.014$）

**閾値非依存テスト**:
- ROC-AUC: ほぼ全分割でMaxent > GARP（有意差、$p < 0.05$）
  - *B. variegatus*: Maxent平均AUC 0.873, GARP 0.789
  - *M. minutus*: Maxent平均AUC 0.986, GARP 0.942

**潜在植生変数の追加**:
- Maxent: AUC向上（*B. variegatus*: 0.873 → 0.880, $p=0.093$）
  - 非森林地域（llanos, cerrado）を正しく除外
- GARP: ほとんど変化なし（0.789 → 0.780）

### 視覚的評価
**B. variegatus**: 両手法とも合理的な予測
- Maxent閾値 ≥ 1、GARP閾値 5-10 が適切
- Maxentは非森林地域（ベネズエラllanosなど）をより正確に除外

**M. minutus**: Maxentが顕著に優れた予測
- Maxent: アンデス山脈に限定（既知分布と一致）
- GARP: 過剰予測（メソアメリカ、ギアナ高地、ブラジル高原まで含む）

## 限界・残された課題

### 著者が述べる限界
1. **成熟度**: GLM/GAMほど成熟していない
   - 誤差推定の方法が少ない
   - Unconditional modelは機械学習で稀

2. **正則化の選択**: 最適な $\beta$ の決定方法に研究が必要
   - 変数選択手法との比較が不十分

3. **外挿時の注意**: 指数モデルは上限がない
   - 研究範囲外の環境条件への外挿時、予測値が極端に大きくなる
   - "Clamping"（範囲外の値を上下限にリセット）が必要

4. **専用ソフトウェア必要**: 標準統計パッケージに未実装

### 追加の課題
1. **サンプリングバイアス**: 道路沿い、特定地域への偏りの形式的対処は将来の課題
   - Zadrozny (2004) のアプローチを適用可能性あり

2. **Threshold選択**: 二値予測が必要な場合の閾値決定に明確なガイドラインが不足
   - Cumulative representation理論的基盤はあるが、実践的指針が不十分

3. **真のAbsenceとの比較**: Pseudo-absencesの質が評価できない
   - 検出確率を考慮したoccupancy modelとの統合が必要

4. **計算コスト**: 大規模データセットでの実行時間
   - Background pixelsのサンプリングで軽減可能

## 本研究(MMCP)との関係

### 直接的な関係
1. **IPPとの理論的同値性**: Renner et al. (2015) が後に示した、Maxent ≈ IPP with logistic link
   - MMCPの点過程部分でMaxentの数学的枠組みを援用可能
   - Background points = quadrature points in IPP

2. **Presence-onlyデータの原理的扱い**: 黒曜石遺跡データはpresence-only
   - Absenceを"扱わない"のではなく、"別の確率分布"として定式化
   - より原理的なアプローチ

3. **Feature engineeringの示唆**: 環境変数の関数（quadratic, product）を特徴として使用
   - MMCPでも、距離の2乗、距離と標高の積などを特徴に
   - 非線形応答曲線の柔軟なモデル化

4. **正則化の重要性**: ℓ1正則化で過学習を防止
   - MMCPでも、多数の空間基底関数を使う際に正則化が必須
   - Sparse solutionによる解釈性向上

### 本研究への示唆
- **最大エントロピー原理の適用**:
  - IPP強度関数推定で、既知情報（presence localities）以外の仮定を最小化
  - 自然な確率分布（exponential family）の選択

- **GLMとの関係明確化**:
  - Logistic regression with pseudo-absences ≈ Maxent
  - どちらを使うかは解釈と計算効率のトレードオフ

- **組成データへの拡張**:
  - Maxentは単変量応答だが、multivariate Maxent も理論的に可能
  - 各産地の組成を別々のMaxentでモデル化 → 組成の空間変動

- **ソフトウェア活用**:
  - Maxent software (現在はdismo Rパッケージに統合)
  - MMCPの点位置部分の初期モデルとして利用可能

### 借用すべき手法
1. **Sequential-update algorithm**: 凸最適化の効率的アルゴリズム
2. **Cumulative representation**: 確率値の解釈可能な表現
3. **Cross-validation with random splits**: 10回ランダム分割による頑健な評価
4. **ROC-AUC evaluation**: Presence-onlyデータでの標準評価指標

### 拡張すべき点
1. **Marked point process**: 点位置 + マーク（組成）の同時モデリング
   - Maxentは点位置のみ、MMCPは2段階構造

2. **空間相関の明示的モデル**: Maxentは独立ピクセル仮定
   - MMCPではGaussian processやCAR priorで空間相関をモデル化

3. **ベイズ推論**: Maxentは点推定（MAP）
   - MMCPでは事後分布全体を推定、不確実性定量化

4. **検出確率の組み込み**: Maxentは完全検出を仮定
   - Occupancy modelとの統合（Dorazio 2014）

## 引用すべき箇所

### Presence-onlyデータの価値（Introduction, p.2）
> "Modeling techniques that require only presence data are therefore extremely valuable, as vast stores of presence-only data exist (particularly in natural history museums and herbaria), whereas absence data are rarely available, especially for poorly sampled tropical regions where modeling potentially has the most value for conservation."

**用途**: 考古学データ（博物館所蔵遺物）のpresence-only性質を正当化

### 最大エントロピー原理（Methods, p.6）
> "The best approach is to ensure that the approximation satisfies any constraints on the unknown distribution that we are aware of, and that subject to those constraints, the distribution should have maximum entropy. This is known as the maximum-entropy principle."

**用途**: Maxent手法の哲学的基盤を説明

### GLMとの関係（Methods, p.11）
> "Despite these similarities, important differences exist between GLM/GAMs and Maxent, causing them to make different predictions. When GLM/GAMs are used to model probability of occurrence, absence data are required. When applied to presence-only data, background pixels must be used instead of true absences."

**用途**: Maxent vs GLM with pseudo-absences の違いを明確化

### Generative vs Discriminative（Methods, p.12）
> "Maxent is a generative approach, whereas GLM/GAMs are discriminative, and generative methods may give better predictions when the amount of training data is small (Ng & Jordan, 2001)."

**用途**: 少数サンプルでのMaxentの利点を示す

### 正則化の効果（Methods, p.8）
> "Regularization forces Maxent to focus on the most important features, and ℓ1-regularization tends to produce models with few nonzero λ_j values. Such models are less likely to overfit, because they have fewer parameters; as a general rule, the simplest explanation of a phenomenon is usually best (the principle of parsimony, Occam's Razor)."

**用途**: スパース解の利点とOccamの剃刀を議論

### 実証結果の要約（Results, p.13）
> "Both algorithms provided reasonable estimates of the species' range, far superior to the shaded outline maps available in field guides. All models were significantly better than random in both binomial tests of omission and receiver operating characteristic (ROC) analyses. The area under the ROC curve (AUC) was almost always higher for Maxent, indicating better discrimination of suitable versus unsuitable areas for the species."

**用途**: Maxentの有効性を実証データで示す

### 熱力学的根拠（Methods, p.10）
> "The applicability of the maximum entropy principle to species distributions is supported by thermodynamic theories of ecological processes. The second law of thermodynamics specifies that in systems without outside influences, processes move in a direction that maximizes entropy."

**用途**: 最大エントロピー原理の生態学的妥当性を議論

### Future work（Discussion, 未読部分からの推測）
Maxentの将来的展開として、サンプリングバイアスの形式的扱い、conditional modelによる presence/absence 統合などが議論されていると予想。

**用途**: MMCPでの拡張方向を示唆
