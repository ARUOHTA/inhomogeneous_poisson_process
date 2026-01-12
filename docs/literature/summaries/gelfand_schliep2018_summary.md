# Gelfand & Schliep (2018) - 空間点パターンのベイズ推論と計算

## 基本情報
- **タイトル**: Bayesian Inference and Computing for Spatial Point Patterns
- **著者**: Alan E. Gelfand, Erin M. Schliep
- **出版**: NSF-CBMS Regional Conference Series in Probability and Statistics, Vol. 10
- **年**: 2018
- **DOI**: 確認中
- **関連論点**: 論点3（点過程）、全般的な理論的基盤

## 主な貢献

1. **空間点パターンのベイズ推論の体系的教科書**
   - 点過程の理論からベイズ推論、計算手法までを統合的に解説
   - NSF-CBMS講義シリーズの教科書として出版

2. **階層モデリングの枠組みの提示**
   - `[data|process,parameters][process|parameters][parameters]` の3段階構造
   - 点過程を階層ベイズモデルとして扱う統一的フレームワーク

3. **モデル評価手法の議論**
   - 予測空間でのモデル比較（パラメータ空間ではなく）
   - CRPS, PMSE等のスコアリングルール
   - 事前/事後予測チェック

4. **Gaussian Processの詳細な解説**
   - 空間モデリングにおけるGPの役割
   - 共分散関数の選択と特性

## 手法の概要

### 空間データの3分類（Section 1.1）

1. **点参照データ（Geostatistical data）**
   - $Y(\mathbf{s})$が連続的な領域$D$上の選択された位置で観測
   - 例: 142地点でのツガの個体数

2. **地域データ（Areal data）**
   - 領域$D$が有限個の領域単位に分割
   - 例: グリッドセル単位の衛星観測データ

3. **点パターンデータ（Point pattern data）** ← 本書の焦点
   - 位置$\mathbf{s} \in D$自体がランダム
   - $Y(\mathbf{s}) = 1$（occurrence）、0が非可算無限個
   - 例: 森林内の樹木の位置

### マーク付き点過程の方向性

> "Introduction of marks forces us to think about whether we model the distribution of the marks and then provide a model for the points given the mark or whether we model the pattern of points and then assign a random mark to each point. The distributional specifications are quite different according to the direction of conditioning."

2つのアプローチ：
1. **マーク → 点**: $[マーク分布][点 | マーク]$
2. **点 → マーク**: $[点分布][マーク | 点]$

### ベイズ推論の原理（Section 1.2）

**基本形**:
$$
f(\boldsymbol{\theta} | \mathbf{Y}) = \frac{f(\mathbf{Y} | \boldsymbol{\theta})\pi(\boldsymbol{\theta})}{f(\mathbf{Y})} \propto f(\mathbf{Y} | \boldsymbol{\theta})\pi(\boldsymbol{\theta})
$$

**事後予測分布**:
$$
f(Y_0 | \mathbf{Y}) = \int f(Y_0 | \mathbf{Y}, \boldsymbol{\theta}) f(\boldsymbol{\theta} | \mathbf{Y}) d\boldsymbol{\theta}
$$

**ベイズ更新**: $Y_1$から$\pi(\theta | y_1)$を得て、これを$Y_2$の事前分布として使用

### 階層モデリング（Section 1.3）

**一般形**:
$$
f(\mathbf{y} | \boldsymbol{\theta})\pi(\boldsymbol{\theta} | \boldsymbol{\lambda})h(\boldsymbol{\lambda})
$$

- **第1段**: データ | プロセス, パラメータ
- **第2段**: プロセス | パラメータ
- **第3段**: パラメータの事前分布

**階層モデルの例**:
- 条件付き独立階層モデル（CIHM）
- ランダム効果モデル
- 欠測データと代入
- 潜在変数モデル
- 混合モデル

### モデル評価の原則

**予測空間での評価を推奨**（パラメータ空間ではなく）:

1. **事前 vs 事後予測チェック**
   - 事前予測: $f(\mathbf{Y}_{\text{rep}} | \text{model}) = \int f(\mathbf{Y}_{\text{rep}} | \boldsymbol{\theta})\pi(\boldsymbol{\theta})d\boldsymbol{\theta}$
   - 事後予測: $f(\mathbf{Y}_{\text{rep}} | \text{model}, \mathbf{Y}_{\text{obs}}) = \int f(\mathbf{Y}_{\text{rep}} | \boldsymbol{\theta})\pi(\boldsymbol{\theta} | \mathbf{Y}_{\text{obs}})d\boldsymbol{\theta}$

2. **Cross-validation**
   - 訓練データと検証データの分割
   - k-fold cross-validation

3. **評価基準**
   - PMSE（予測平均二乗誤差）
   - PMAE（予測平均絶対誤差）
   - CRPS（連続ランクプロバビリティスコア）
   - 予測区間の経験的カバレッジ vs 名目カバレッジ

**カバレッジの解釈**:
- 経験的 << 名目: 不確実性を過小評価（モデル不適切）
- 経験的 >> 名目: 不確実性を過大評価（過適合）

### Gaussian Process（Section 1.4）

**定義**: 任意の$n \geq 1$と位置集合$\{\mathbf{s}_1, \ldots, \mathbf{s}_n\}$に対し、$\mathbf{Y} = (Y(\mathbf{s}_1), \ldots, Y(\mathbf{s}_n))$が多変量正規分布に従う

**GPの利点**:
1. 平均関数と共分散関数で有限次元分布が決定
2. 周辺・条件付き分布が容易に計算可能
3. 階層モデルでのランダム効果として自然
4. 強定常性 = 弱定常性（Gaussianの場合）
5. 分布仮定の批判が困難（サンプルサイズ1）

**等方性共分散関数の例**（Table 1.1）:

| モデル | 共分散関数 $C(\|\mathbf{h}\|)$ |
|--------|-------------------------------|
| Spherical | $\sigma^2[1 - \frac{3}{2}\phi\|\mathbf{h}\| + \frac{1}{2}(\phi\|\mathbf{h}\|)^3]$ if $\|\mathbf{h}\| \leq 1/\phi$ |
| Exponential | $\sigma^2 \exp(-\phi\|\mathbf{h}\|)$ |
| Gaussian | $\sigma^2 \exp(-\phi^2\|\mathbf{h}\|^2)$ |
| Matérn ($\nu=3/2$) | $\sigma^2(1+\phi\|\mathbf{h}\|)\exp(-\phi\|\mathbf{h}\|)$ |

**Matérn共分散関数**:
$$
C(t) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(2\sqrt{\nu}\|\mathbf{h}\|\phi)^\nu K_\nu(2\sqrt{\nu}\|\mathbf{h}\|\phi)
$$

- $\nu$: 平滑性パラメータ
- $\nu = 1/2$: 指数関数
- $\nu \to \infty$: ガウス関数
- $\lfloor \nu \rfloor$: 平均二乗微分可能回数

## 本書の構成（推測）

Chapter 1で読めた内容から、後続章の構成は以下と推測：
- Chapter 2: 点過程の理論（Poisson過程、Cox過程等）
- Chapter 3: マーク付き点過程
- Chapter 4: ベイズ計算手法（MCMC, Gibbs sampling）
- Chapter 5以降: 応用例

## 限界・残された課題

### 教科書としての限界

1. **理論的扱い**
   - 教科書なので、最新の計算手法（NNGP等）は含まれない可能性
   - 2018年出版のため、それ以降の発展は含まれない

2. **応用範囲**
   - Presence-onlyデータへの言及はあるが、詳細は後続章次第

### 本研究の観点からの限界

1. **組成値マークの未対応**
   - マーク付き点過程の議論はあるが、組成データ特有の扱いは不明

2. **Pólya-Gamma拡張への言及なし**
   - Polson et al. (2013)の技法は本書出版時には利用可能だが、点過程への適用は含まれていない可能性

## 本研究（MMCP）との関係

### 借用する要素

1. **階層ベイズフレームワーク**
   - `[data|process,parameters][process|parameters][parameters]` の3段階構造
   - MMCPの数学的記述の基盤

2. **Gaussian Processの理論**
   - 空間依存構造のモデリング
   - 共分散関数の選択指針

3. **モデル評価の原則**
   - 予測空間での評価
   - CRPSやカバレッジチェック

4. **マーク付き点過程の概念**
   - マーク → 点 vs 点 → マークの条件付け方向

### 位置づけ

Gelfand & Schliep (2018)は、空間点パターンのベイズ推論における**教科書的基盤**を提供する。MMCPは：
- この階層ベイズフレームワークを継承
- 組成値マーク（単体上の値）という特殊ケースへの拡張
- Pólya-Gamma拡張による計算効率化
- Presence-onlyデータへの対応

を実現する。本書は理論的・哲学的基盤を、後続の技術的論文（Moreira & Gamerman 2022等）は具体的実装を提供する関係。

## 引用すべき箇所

### 階層モデリングの枠組み（Section 1.1）
> "We might conceptualize a process model specification of the form [data|process,parameters][process|parameters][parameters]. The first stage distribution brings in the data revealing how it is connected to the process."

### 点パターンデータの定義（Section 1.1）
> "Point pattern data, where now the set of locations in $D$ are themselves random; its index set gives the locations of random events that are the spatial point pattern."

### マーク付き点過程の条件付け（Section 1.1）
> "Introduction of marks forces us to think about whether we model the distribution of the marks and then provide a model for the points given the mark or whether we model the pattern of points and then assign a random mark to each point."

### ベイズ推論の優位性（Section 1.2）
> "One obtains an entire posterior distribution for an unknown rather than perhaps a point estimate and an asymptotic variance, as with usual classical inference."

### 予測空間での評価（Section 1.2）
> "Since the parameter space varies with the model [...] we avoid criteria which operate in the parameter space. We only consider model comparison (and model adequacy) in predictive space."

### Gaussian Processの役割（Section 1.4）
> "Gaussian processes play a crucial role in spatial modeling. They provide extremely flexible specifications for introducing spatial and spatio-temporal dependence."

### 階層モデリングの注意点（Section 1.3）
> "Hierarchical models must be handled with care. [...] We often specify models which are too large for the data to support, meaning we are overfitting the data."
