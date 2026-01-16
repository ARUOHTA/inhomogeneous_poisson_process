# Aitchison (1982) - 組成データの統計分析

## 基本情報
- **タイトル**: The Statistical Analysis of Compositional Data
- **著者**: J. Aitchison
- **所属**: University of Hong Kong
- **掲載**: Royal Statistical Society（Research Section、1982年1月13日講演）
- **年**: 1982
- **DOI**: 記載なし
- **関連論点**: 論点2（組成データの統計手法）、論点3（変換手法）

## 主な貢献
- **シンプレックスにおける統計分析の理論的枠組み**: 組成データの標本空間としてのシンプレックス$\mathbb{S}^d$の統計的扱い
- **Logistic-normal分布の導入**: 変換による多変量正規分布からシンプレックス上の分布を構成
- **独立性概念の体系化**: シンプレックス上の複数の独立性概念を定義し相互関係を明確化
- **Dirichlet分布の限界の指摘**: 既存の唯一のパラメトリッククラスが強すぎる独立性構造を持つことを示す
- **変換手法の理論的基盤**: additive logistic, multiplicative logistic, hybrid logistic変換の導入

## 手法の概要

### シンプレックスの定義と問題設定
**シンプレックス**（正単体）:

$$
\mathbb{S}^d = \{(x_1, \ldots, x_d): x_i > 0 \, (i=1, \ldots, d), \, x_1 + \ldots + x_d < 1\}
$$

**組成データ** (compositional data): シンプレックス$\mathbb{S}^d$上の点$\mathbf{x}$の集合

**問題の難しさ**:
1. 独立性の概念が不明確
2. 依存性の測度が不足
3. 満足なパラメトリック分布クラスの欠如
4. **閉じた配列の問題** (constant/bounded sum problem, closed arrays): Pearson (1897)のスプリアス相関以来の課題

### 基本操作

#### 1. 組成の基底 (Basis)
実際の量$\mathbf{w}^{(d+1)} \in \mathbb{P}^{d+1}$（正直交体）から組成への変換：

$$
\mathbf{x}^{(d+1)} = C(\mathbf{w}^{(d+1)}), \quad x_i = \frac{w_i}{T(\mathbf{w}^{(d+1)})}, \, i=1, \ldots, d+1
$$

ここで$T(\mathbf{w}^{(d+1)}) = w_1 + \ldots + w_{d+1}$

#### 2. 部分組成 (Subcomposition)
組成$\mathbf{x}^{(d+1)}$の部分ベクトル$\mathbf{x}^{(c)}$から構成される組成：

$$
C(\mathbf{x}^{(c)}) \in \mathbb{S}^{c-1}
$$

#### 3. 合併 (Amalgamation)
成分のグループ化による次元削減：

$$
t_j = x_{c_{j-1}+1} + \ldots + x_{c_j}, \quad j=1, \ldots, k+1
$$

ここで$0 = c_0 < c_1 < \ldots < c_k < c_{k+1} = d+1$

行列形式: $\mathbf{t}^{(k+1)} = \mathbf{A} \mathbf{x}^{(d+1)}$ （$\mathbf{A}$は0と1からなる行列）

#### 4. 分割 (Partition)
合併$\mathbf{t}$と各グループの部分組成$\mathbf{s}_j$の組：

$$
P(\mathbf{x}^{(d+1)}) = (\mathbf{t}; \mathbf{s}_1, \ldots, \mathbf{s}_{k+1})
$$

ここで部分組成$\mathbf{s}_j \in \mathbb{S}^{d_j}$ （$d_j = c_j - c_{j-1} - 1$）は：

$$
s_{jr} = \frac{x_{c_{j-1}+r}}{t_j}, \quad r=1, \ldots, d_j+1
$$

**重要な性質**: 変換$P: \mathbb{S}^d \to \mathbb{S}^k \times \prod_{j=1}^{k+1} \mathbb{S}^{d_j}$は一対一で、ヤコビアンは：

$$
\frac{D\mathbf{x}^{(d)}}{D(\mathbf{t}; \mathbf{s}_1, \ldots, \mathbf{s}_{k+1})} = t_1^{d_1} \cdots t_{k+1}^{d_{k+1}}
$$

### Dirichlet分布とその限界

**Dirichlet分布** $D^d(\boldsymbol{\alpha})$（パラメータ$\boldsymbol{\alpha} \in \mathbb{P}^{d+1}$）:

$$
f(\mathbf{x}^{(d)}) = \frac{\prod_{i=1}^{d+1} x_i^{\alpha_i - 1}}{\Delta(\boldsymbol{\alpha})}, \quad \mathbf{x}^{(d)} \in \mathbb{S}^d
$$

ここで$\Delta(\boldsymbol{\alpha}) = \Gamma(\alpha_1) \cdots \Gamma(\alpha_{d+1}) / \Gamma(\alpha_1 + \ldots + \alpha_{d+1})$

**Dirichlet分布の限界**:

1. **凸等確率輪郭の制約**: $\alpha_i > 1$ のとき、等確率輪郭が常に凸 → 凹パターン（多峰性など）を表現できない

2. **強すぎる独立性構造**:
   - **性質D1**: Dirichlet組成は、同じスケールパラメータを持つ独立なgamma分布の基底から構成可能
   - **性質D2**: 任意の分割$P(\mathbf{x}^{(d+1)}) = (\mathbf{t}; \mathbf{s}_1, \ldots, \mathbf{s}_{k+1})$に対して、$\mathbf{t} \Perp\!\!\!\Perp \mathbf{s}_1 \Perp\!\!\!\Perp \cdots \Perp\!\!\!\Perp \mathbf{s}_{k+1}$（相互独立）

→ Dirichlet分布は独立性仮説の極致であり、実データの柔軟なモデリングには不適

### 変換正規分布クラス (Transformed Normal Classes)

**基本アイデア**: 多変量正規分布$N^d(\boldsymbol{\mu}, \boldsymbol{\Sigma})$ on $\mathbb{R}^d$から、変換$f: \mathbb{R}^d \to \mathbb{S}^d$により$\mathbb{S}^d$上の分布$fN^d(\boldsymbol{\mu}, \boldsymbol{\Sigma})$を誘導

#### 基本的な変換（Table 1）

**1. Additive logistic変換** $a_d$:

$$
\begin{cases}
x_i = \frac{\exp(y_i)}{1 + \sum_{j=1}^d \exp(y_j)}, & i=1, \ldots, d \\
x_{d+1} = \frac{1}{1 + \sum_{j=1}^d \exp(y_j)}
\end{cases}
$$

逆変換：$y_i = \log(x_i / x_{d+1})$, $i=1, \ldots, d$

**2. Multiplicative logistic変換** $m_d$:

$$
y_i = \log \frac{x_i}{1 - \sum_{j=1}^{i} x_j}, \quad i=1, \ldots, d
$$

**3. Hybrid logistic変換** $h_d$:

$$
\begin{cases}
y_1 = \log \frac{x_1}{1 - x_1} \\
y_i = \log \frac{x_i}{(1 - \sum_{j=1}^{i-1} x_j)(1 - \sum_{j=1}^{i} x_j)}, & i=2, \ldots, d
\end{cases}
$$

**共通の性質**: 全ての変換でヤコビアン$D\mathbf{x}/D\mathbf{y} = x_1 x_2 \cdots x_{d+1}$

#### 複雑な変換の構成法

**A. 線形変換法**: $\mathbf{y}$を$\mathbf{Qy}$（$\mathbf{Q}$: 非特異行列）で置き換える

**例**: $q_{ii}=1$, $q_{i,i+1}=-1$ (差分)

$$
y_i = \log(x_i / x_{i+1}), \quad i=1, \ldots, d
$$

→ 隣接成分の比を用いた変換

**B. 分割変換法**: 分割$P(\mathbf{x}^{(d+1)}) = (\mathbf{t}; \mathbf{s}_1, \ldots, \mathbf{s}_{k+1})$に対して：

$$
f_0: \mathbb{R}^k \to \mathbb{S}^k, \quad f_j: \mathbb{R}^{d_j} \to \mathbb{S}^{d_j}, \, j=1, \ldots, k+1
$$

を選び、複合変換$(f_0; f_1, \ldots, f_{k+1})$を構成

### Logistic-normal分布

**定義**: $\mathbf{y} \sim N^d(\boldsymbol{\mu}, \boldsymbol{\Sigma})$のとき、$\mathbf{x} = a_d^{-1}(\mathbf{y})$の分布を**additive logistic-normal分布**と呼ぶ

**密度関数**:

$$
f(\mathbf{x}) = \frac{1}{(2\pi)^{d/2} |\boldsymbol{\Sigma}|^{1/2} x_1 \cdots x_{d+1}} \exp\left\{-\frac{1}{2}(\mathbf{y} - \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{y} - \boldsymbol{\mu})\right\}
$$

ここで$\mathbf{y} = (\log(x_1/x_{d+1}), \ldots, \log(x_d/x_{d+1}))^T$

**利点**:
- 多峰性、凹パターンを表現可能
- 正規分布理論の豊富な理論・手法が利用可能
- 共分散行列$\boldsymbol{\Sigma}$による柔軟な依存構造

## 独立性概念の体系化

### 分割$(x^{(c)}, x_{(c)})$または$(t; \mathbf{s}_1, \mathbf{s}_2)$（order 1）に対する独立性

**表記**（Dawid 1979）:
- $C(\mathbf{x}^{(c)}) \Perp\!\!\!\Perp C(\mathbf{x}_{(c)})$: 部分組成の独立
- $C(\mathbf{x}^{(c)}) \Perp\!\!\!\Perp C(\mathbf{x}_{(c)}) | T(\mathbf{x}^{(c)})$: 条件付き独立
- $\Perp\!\!\!\Perp \mathbf{w}^{(d+1)}$: $\mathbf{w}^{(d+1)}$が独立成分

#### 主要な独立性仮説

**1. Subcompositional independence** $\mathscr{S}$:
$$
\mathbf{s}_1 \Perp\!\!\!\Perp \mathbf{s}_2
$$

**2. Conditional subcompositional independence** $\mathscr{C}$:
$$
\mathbf{s}_1 \Perp\!\!\!\Perp \mathbf{s}_2 | t
$$

**3. Invariance** $\mathscr{I}_1, \mathscr{I}_2$:
- $\mathscr{I}_1$: $\mathbf{s}_1 \Perp\!\!\!\Perp t$ （部分組成$\mathbf{s}_1$が合併$t$と独立）
- $\mathscr{I}_2$: $\mathbf{s}_2 \Perp\!\!\!\Perp t$ （部分組成$\mathbf{s}_2$が合併$t$と独立）

**4. Neutrality** $\mathscr{N}_1, \mathscr{N}_2$:
- $\mathscr{N}_1 = \mathscr{C} \cap \mathscr{I}_1$: 条件付き独立 + $\mathbf{s}_1$が$t$と独立
- $\mathscr{N}_2 = \mathscr{C} \cap \mathscr{I}_2$: 条件付き独立 + $\mathbf{s}_2$が$t$と独立

**5. Partition independence** $\mathscr{P}$:
$$
\mathscr{P} = \mathscr{N}_1 \cap \mathscr{N}_2 = \mathbf{t} \Perp\!\!\!\Perp \mathbf{s}_1 \Perp\!\!\!\Perp \mathbf{s}_2
$$

（合併と2つの部分組成が相互独立）

#### 独立性仮説の階層構造（Figure 2）

$$
\begin{array}{ccc}
& \mathscr{C} & \\
/ & & \backslash \\
\mathscr{N}_1 & & \mathscr{N}_2 \\
\backslash & & / \\
& \mathscr{P} &
\end{array}
$$

**含意関係**:
- $\mathscr{P} \Rightarrow \mathscr{N}_1, \mathscr{N}_2$
- $\mathscr{N}_1, \mathscr{N}_2 \Rightarrow \mathscr{C}$
- $\mathscr{N}_1 \Rightarrow \mathscr{I}_1$, $\mathscr{N}_2 \Rightarrow \mathscr{I}_2$

### レベル$c$までの独立性

**定義**: 組成$\mathbf{x}^{(d+1)}$が独立性$H$を**レベル$c$まで**持つ ($H^c$) $\Leftrightarrow$ $H_k$ が$k=1, \ldots, c$で成立

**完全独立性** (complete independence property): $H^{d-1}$が成立

**例**: **Complete neutrality** $\mathscr{N}^{d-1}$ ← Connor & Mosimann (1969)

## パラメトリック検定

### 条件付きモデル

分割$(t; \mathbf{s}_1, \mathbf{s}_2)$に対して、条件付きモデル$(\mathbf{s}_1, \mathbf{s}_2 | t)$を構成：

$$
\mathbf{y}_1 = a_{c-1}^{-1}(\mathbf{s}_1) = \log(\mathbf{s}_1^{(c-1)} / s_{1c}), \quad \mathbf{y}_2 = a_{d-c}^{-1}(\mathbf{s}_2) = \log(\mathbf{s}_2^{(d-c)} / s_{2,d-c+1})
$$

$$
z = \log(t / (1-t))
$$

**モデル**$M$:

$$
(\mathbf{y}_1, \mathbf{y}_2 | z) \sim N^{d-1}\left(\begin{bmatrix} \boldsymbol{\alpha}_1 + \boldsymbol{\beta}_1 z \\ \boldsymbol{\alpha}_2 + \boldsymbol{\beta}_2 z \end{bmatrix}, \begin{bmatrix} \boldsymbol{\Sigma}_{11} & \boldsymbol{\Sigma}_{12} \\ \boldsymbol{\Sigma}_{21} & \boldsymbol{\Sigma}_{22} \end{bmatrix}\right)
$$

### パラメトリック対応（Table 2）

| 仮説 | パラメトリック制約 | 自由度 |
|------|-------------------|--------|
| $\mathscr{I}_1$ | $\boldsymbol{\beta}_1 = \mathbf{0}$ | $c-1$ |
| $\mathscr{I}_2$ | $\boldsymbol{\beta}_2 = \mathbf{0}$ | $d-c$ |
| $\mathscr{C}$ | $\boldsymbol{\Sigma}_{12} = \mathbf{0}$ | $(c-1)(d-c)$ |
| $\mathscr{N}_1$ | $\boldsymbol{\beta}_1 = \mathbf{0}, \boldsymbol{\Sigma}_{12} = \mathbf{0}$ | $(c-1)(d-c+1)$ |
| $\mathscr{N}_2$ | $\boldsymbol{\beta}_2 = \mathbf{0}, \boldsymbol{\Sigma}_{12} = \mathbf{0}$ | $c(d-c)$ |
| $\mathscr{P}$ | $\boldsymbol{\beta}_1 = \mathbf{0}, \boldsymbol{\beta}_2 = \mathbf{0}, \boldsymbol{\Sigma}_{12} = \mathbf{0}$ | $c(d-c) + c - 1$ |

### 尤度比検定

**検定統計量** (仮説$H$内でモデル$M$を検定):

$$
T = n \log \frac{|\hat{\boldsymbol{\Sigma}}_H|}{|\hat{\boldsymbol{\Sigma}}_M|} \sim \chi^2_{q_H}
$$

ここで$q_H$は制約の数

### 仮説の格子（Lattice）における検定戦略

Jeffreys (1961)の**簡潔性の公準** (simplicity postulate):
1. 最も単純な仮説$\mathscr{P}$から開始
2. 棄却されたら次のレベル（$\mathscr{N}_1, \mathscr{N}_2$）へ
3. 両方棄却されたら$\mathscr{C}$へ、以下同様

## 実証結果（応用例）

### Example 1: Skye lavas（スカイ島玄武岩）
- **データ**: 32標本、10主要酸化物組成
- **分割**: $(A, F, M | S, C, N, K, T, P)$ （AFM図の検証）
  - $A = \text{Al}_2\text{O}_3$, $F = \text{Fe}_2\text{O}_3$, $M = \text{MgO}$
  - $S = \text{SiO}_2$, $C = \text{CaO}$, $N = \text{Na}_2\text{O}$, $K = \text{K}_2\text{O}$, $T = \text{TiO}_2$, $P = \text{P}_2\text{O}_5$

**結果**（Figure 3の格子）:
- 全ての独立性仮説が棄却（有意水準 $p < 0.001$）
- $\mathscr{P}$, $\mathscr{N}_1$, $\mathscr{N}_2$, $\mathscr{C}$ すべて強く棄却
- **結論**: AFM部分組成は補完的な組成要素と独立ではない → AFM図単独の分析は疑わしい

**Neutrality検定**（Table 4）:
- 順序$(A|F|M|S|C|N|K|T|P)$でレベル$c=1, \ldots, 7$まで検定
- 全レベルで$\mathscr{N}^c$が強く棄却

### Example 2: Arctic lake sediments（北極湖堆積物）
- **データ**: 39標本、砂・シルト・粘土組成 + 水深$z$
- **モデル**: $\mathbf{x}^{(d+1)} | z \sim a_2N^2(\mathbf{g}(z), \boldsymbol{\Sigma})$
- **回帰関数**: $\mathbf{g}(z) = \boldsymbol{\alpha} + \boldsymbol{\beta} z + \boldsymbol{\gamma} z^2 + \boldsymbol{\delta} \log z + \boldsymbol{\epsilon} (\log z)^2$

**結果**:
- 線形回帰は棄却
- $\mathbf{g}(z) = \boldsymbol{\alpha} + \boldsymbol{\beta} \log z$（対数線形）または2次回帰が適合
- 残差が多変量正規性検定を通過

### Example 3: Glacial tills（氷河堆積物）
- **データ**: 93標本、4種類の礫組成 + 総礫数（サイズ情報）
- **分析**: 組成がサイズに依存するかを検証

### Example 4: Hong Kong household budgets（香港家計調査）
- **データ**: 2つの低所得住宅カテゴリーA, B（各41, 42世帯）、7支出グループ
- **回帰**: $\boldsymbol{\alpha} + \boldsymbol{\beta} \log(\text{total expenditure}) + \boldsymbol{\gamma} \log(\text{household size})$

**結果**:
- 仮説$\boldsymbol{\gamma} = \mathbf{0}$のみ棄却されない（世帯規模の影響なし）
- 総支出の影響は有意
- 適合残差が多変量正規性検定を通過

## 限界・残された課題

### 理論的限界

1. **Dirichlet分布を含む一般化の欠如**:
   - シンプレックス上の既存の独立性を持たない分布族の開発が未解決
   - Connor & Mosimann (1969), Darroch & James (1974), Mosimann (1975), James & Mosimann (1980)らの試みは限定的成功

2. **変換の選択**:
   - 指数関数以外の変換 $\mathbb{R}^1 \to \mathbb{P}^1$ の探索余地
   - 応用に応じた最適変換の選択基準が不明確

3. **$\mathscr{I}^c$と$\mathscr{N}^c$の関係**:
   - 予想: 変換正規モデルの枠組みでは$\mathscr{I}^c \equiv \mathscr{N}^c$（未証明）
   - $\mathscr{C}^c$と$\mathscr{N}^c$は明確に異なる

4. **ゼロ成分の問題** (Zero components):
   - 実データで頻繁に発生（支出項目なし、微量元素検出限界以下）
   - 対処: 小さな正の値（例: $0.05）で置換
   - 理論的に満足な解決策なし（Section 7.4で議論）

### 実装上の課題

5. **多変量正規性の検定**:
   - Andrews, Gnanadesikan & Warner (1973)の手法を利用
   - 大規模データでの計算負荷

6. **非正規性への対処**:
   - 変換しても正規性が得られない場合の代替手法
   - ノンパラメトリック手法への拡張

7. **高次元データ**:
   - 成分数$d$が大きい場合の推定精度低下
   - 次元削減・正則化の必要性

### 応用上の課題

8. **地質学的解釈**:
   - AFM図の独立性仮定の棄却 → 地球化学的解釈の再考
   - 既存の三角図法（ternary diagrams）の限界

9. **経済学的解釈**:
   - 家計支出の組成分析における因果推論
   - エンゲル則との関係

10. **サンプリング設計**:
    - 組成データに特化したサンプルサイズ決定
    - 情報的事前分布の構築

## 本研究（MMCP）との関係

### 借用すべき手法

1. **Additive logistic変換（alr変換）**:
   - **直接的な適用**: MMCPで黒曜石組成データ$\mathbf{x} \in \mathbb{S}^{d-1}$を$\mathbb{R}^{d-1}$に変換
   - $\mathbf{y} = a_{d-1}^{-1}(\mathbf{x}) = \log(\mathbf{x}^{(d-1)} / x_d)$
   - プロジェクトで使用しているilr変換の理論的基盤

2. **シンプレックスの基本操作**:
   - **Subcomposition**: 黒曜石産地ごとの組成比較（例: 信州 vs 神津島産の元素比）
   - **Amalgamation**: 微量元素の統合（例: 希土類元素の合計）
   - **Partition**: 主要元素 + 微量元素の階層構造

3. **Logistic-normal分布**:
   - 黒曜石組成の確率モデル
   - MMCPの尤度関数$f(\mathbf{y}_i | \mathbf{x}_i) = LN(\boldsymbol{\mu}(\mathbf{x}_i), \boldsymbol{\Sigma})$
   - 多峰性（複数産地の混合）を表現可能

4. **独立性検定の枠組み**:
   - 遺跡間での組成独立性の検定
   - 空間位置と組成の依存性検証
   - 条件付き独立性$\mathbf{s}_1 \Perp\!\!\!\Perp \mathbf{s}_2 | t$の応用

5. **回帰モデル**:
   - 組成 ← 距離、標高、時代などの共変量
   - $\mathbf{y} | \mathbf{z} \sim N^{d-1}(\boldsymbol{\alpha} + \mathbf{B}\mathbf{z}, \boldsymbol{\Sigma})$
   - トブラー距離を共変量とする組成回帰

### 拡張すべき点

1. **空間構造への統合**:
   - 本論文: 独立サンプルの組成分析
   - MMCP: 空間的に相関する遺跡の組成 → 空間Cox過程との統合が必要
   - **課題**: Logistic-normal + 空間相関構造の同時モデリング

2. **点過程への組み込み**:
   - Aitchison (1982): 組成データのみ
   - MMCP: 点過程（遺跡位置）+ マーク（組成）の同時モデル
   - **課題**: 強度関数$\lambda(\mathbf{s}, \mathbf{y})$での組成の扱い

3. **空間的に変化する組成分布**:
   - 産地からの距離で組成分布のパラメータ$\boldsymbol{\mu}(\mathbf{s}), \boldsymbol{\Sigma}(\mathbf{s})$が変化
   - **課題**: 空間的に変化する係数（SVC）との統合

4. **ベイズ推論への拡張**:
   - 本論文: 最尤推定、尤度比検定（頻度論的）
   - MMCP: ベイズ推論（事後分布、予測分布）
   - **課題**: Logistic-normal事前分布、MCMC実装

5. **Dirichlet分布との関係**:
   - Dirichlet $\subset$ Logistic-normal?（一般には含まない）
   - ベイズ推論でDirichlet事前分布を使う場合の整合性
   - **課題**: Dirichlet-multinomial vs Logistic-normal-multinomial

6. **閉じた配列問題の空間版**:
   - Pearson (1897)のスプリアス相関 → 空間的スプリアス相関?
   - 組成制約と空間相関の交互作用
   - **課題**: 組成データ + 空間統計の統合理論

7. **多峰性と産地混合**:
   - 黒曜石組成の多峰性 = 複数産地の混合?
   - Logistic-normal混合モデルの必要性
   - **課題**: 産地ラベルのクラスタリング vs 混合比の推定

### MMCPにおける位置づけ

**第3章での役割**:
- 組成データ分析の理論的基盤として第2.1節で引用
- alr/clr/ilr変換の原典
- Dirichlet分布の限界を指摘し、Logistic-normalの優位性を示す材料

**実装の基礎**:
- `bayesian_statistics/`モジュールでのilr変換実装の理論的根拠
- 組成データの可視化（ternary diagrams）の限界理解
- 変換後の正規性検定の必要性

**方法論的貢献**:
- 組成データを「閉じた配列」ではなく「シンプレックス上の確率分布」として扱う視点
- 独立性概念の体系化 → 空間的独立性への拡張
- 変換による統計的推論 → 空間過程への適用

## 引用すべき箇所

### シンプレックス上の統計分析の難しさ

> "The simplex, however, has proved to be an awkward space to handle statistically; the difficulties appear to lie in the scarcity of meaningful definitions of independence and of measures of dependence and in the absence of satisfactory parametric classes of distributions on $\mathbb{S}^d$." (p.1, Introduction)

→ 組成データ分析の根本的な難しさ

### Dirichlet分布の限界

> "A major obstacle to its use in the statistical analysis of compositional data is that it seldom, if ever, provides an adequate description of actual patterns of variability of compositions." (p.4, Section 2.2)

→ Dirichlet分布が実データに不適な理由

> "The Dirichlet class has so much independence structure built into its definition that it represents, not a convenient modelling class for compositional data but the ultimate in independence hypotheses." (p.4, Section 2.2)

→ Dirichlet分布 = 独立性仮説の極致

### 変換による解決策

> "In our view the way out of the impasse is simply to travel by a different route, escaping from the awkward constrictions of $\mathbb{S}^d$ into the wide open spaces of $\mathbb{R}^d$ through suitably selected transformations between $\mathbb{S}^d$ and $\mathbb{R}^d$." (p.5, Section 2.2)

→ 変換アプローチの動機づけ

### Logistic-normal分布の導入

> "The idea of inducing a tractable class of distributions over some awkward sample space from a proven and well-established class over some simpler space is at least a century old. McAlister (1879), faced with the 'awkward' sample space $\mathbb{P}^1$, saw that if he considered $y$ in $\mathbb{R}^1$ to be $N(\mu, \sigma^2)$ then the transformation $x = \exp(y)$ would induce a useful 'expnormal' distribution..." (p.5, Section 2.3)

→ Logistic-normalのアイデアの歴史的文脈（lognormalの類推）

### Additive logistic変換の定義

> "The additive logistic transformation $a_d: \mathbb{R}^d \to \mathbb{S}^d$ [is] already heavily exploited in other areas of statistical activity, such as logistic discriminant analysis (Cox, 1966; Day and Kerridge, 1967; Anderson, 1972) and in the analysis of binary data (Cox, 1970)." (p.5, Section 2.3)

→ alr変換の他分野での利用

### 独立性概念の重要性

> "Concepts of independence must often play an important role in any form of statistical analysis." (p.1, Introduction)

→ 統計分析における独立性の中心的役割

### AFM図の問題

> "Thus any analysis of AFM which subsumes that this subcomposition is independent of other aspects of the composition is surely suspect." (p.14, Section 6.4)

→ 実証分析の結論（地質学への含意）

### 閉じた配列の問題

> "An appreciation of the difficulty imposed by this confinement of data points [...] to a simplex is inherent in the comments of Pearson (1897) on spurious correlations, and in geological circles the difficulty has since become known as the constant or bounded sum problem and the problem of closed arrays." (p.2, Introduction)

→ Pearson (1897)以来の歴史的課題

### パーティションのヤコビアン

> "The transformation [...] is one-to-one, with Jacobian $D\mathbf{x}^{(d)} / D(\mathbf{t}; \mathbf{s}_1, \ldots, \mathbf{s}_{k+1}) = t_1^{d_1} \cdots t_{k+1}^{d_{k+1}}$." (p.3, Section 2.1)

→ 変数変換の数学的詳細（尤度計算に必要）

### Dirichlet分布の性質D2

> "If $\mathbf{x}^{(d+1)}$ is $D^d(\boldsymbol{\alpha})$ then, for partition (2.6), $\mathbf{t} \Perp\!\!\!\Perp \mathbf{s}_1 \Perp\!\!\!\Perp \cdots \Perp\!\!\!\Perp \mathbf{s}_{k+1}$." (p.4, Section 2.2)

→ Dirichlet分布の強い独立性構造の数学的表現

### 変換正規モデルの柔軟性

> "Rejection of all these hypotheses still leaves transformed normal models on the simplex as possible describers of patterns of variability of non-neutral compositional data." (p.16, Section 7.2)

→ 独立性が棄却されても、変換正規モデルは有効（重要な実用的コメント）
