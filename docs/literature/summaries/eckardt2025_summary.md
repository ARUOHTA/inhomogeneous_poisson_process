# Eckardt et al. (2025) - 組成値マーク付き空間点過程

## 基本情報
- **タイトル**: On spatial point processes with composition-valued marks
- **著者**: Matthias Eckardt, Mari Myllymäki, Sonja Greven
- **ジャーナル**: Spatial Statistics（推定）
- **年**: 2025
- **DOI**: 未確認
- **関連論点**: 論点2（組成データ）、論点3（点過程）、論点4（統合の困難さ）

## 主な貢献

1. **組成値マーク付き空間点過程の新クラス導入**
   - 空間点過程のマークが$D$-部分組成（単体$\mathbb{S}^D$上の値）である場合の理論的フレームワークを初めて体系化

2. **Aitchison幾何学の点過程への統合**
   - 組成データ分析の標準的な対数比変換（alr, clr, ilr）を空間点過程の文脈に統合
   - 単体上の内積・距離を用いた要約特性の定義

3. **2種類の要約特性の開発**
   - **Componentwise（成分別）**: 各組成成分間の空間的依存性を評価
   - **Compositional（組成全体）**: 組成全体としての空間的パターンを評価

4. **グローバルエンベロープテスト**
   - ランダムラベリング仮説（マークが独立同分布）の検定手法を組成値マークに拡張

## 手法の概要

### 組成データの変換

$D$-部分組成$\mathbf{c} = (c_1, \ldots, c_D) \in \mathbb{S}^D$に対して3種類の変換を定義：

1. **対数比変換（lr）**: $\text{lr}_{jl}(\mathbf{c}) = \log(c_j/c_l)$、$D^2$次元
2. **加法対数比（alr）**: $\text{alr}_j(\mathbf{c}) = \log(c_j/c_D)$、$D-1$次元
3. **中心化対数比（clr）**: $\text{clr}_j(\mathbf{c}) = \log(c_j/g(\mathbf{c}))$、$D$次元（ゼロ和制約）
4. **等長対数比（ilr）**: $\text{ilr}(\mathbf{c}) = \text{clr}(\mathbf{c}) \mathbf{H}_D^\top$、$D-1$次元（直交）

ここで$g(\mathbf{c}) = \left(\prod_{j=1}^D c_j\right)^{1/D}$は幾何平均。

### 要約特性

テスト関数$\mathfrak{t}_f$を用いた一般的な要約特性：

$$
\nabla_{\mathfrak{t}_f}^{\psi,jl}(r) = \mathbb{E}_{\circ,r}\left[\mathfrak{t}_f^{\psi,jl}(\psi_j(\mathbf{c}(\circ)), \psi_l(\mathbf{c}(\mathbf{r})))\right]
$$

主な特性：
- **マーク変差関数** $\gamma_{jl}^\psi(r)$: 距離$r$での組成の分散
- **条件付き期待値** $\tau_{jl}^\psi(r)$: 距離$r$でのマーク積の平均
- **Shimantaniの$\iota$関数**: 空間的自己相関の評価

### 組成全体の特性

clr変換を用いた組成全体のマーク変差関数：

$$
\gamma_{\mathbf{cc}}(r) = \sum_{j=1}^D \gamma_{j j}^{\text{clr}}(r) = \sum_{j=1}^{D-1} \gamma_{j j}^{\text{ilr}}(r)
$$

### 絶対情報との混合マーク

組成$\mathbf{c}$と総量$y$を組み合わせた混合マーク$\boldsymbol{\eta} = (y, \mathbf{c})$への拡張：

$$
\gamma_{\mathbf{cc},y}(r) = \mathbb{E}_{\circ,r}\left[\frac{1}{2}d_A(\mathbf{c}(\circ),\mathbf{c}(\mathbf{r}))^2 + \beta \cdot \frac{1}{2}d_+(y(\circ),y(\mathbf{r}))^2\right]
$$

### 推定

カーネル推定量を用いた二次積密度関数の比として推定：

$$
\widehat{\nabla_{\mathfrak{t}_f}^{\psi,jl}}(r) = \widehat{\varrho_{\mathfrak{t}_f}^{\psi,jl,(2)}}(r) / \widehat{\varrho^{(2)}}(r)
$$

### グローバルエンベロープテスト

1. マークを$s$回パーミュテーション
2. 各パーミュテーションでテスト統計量を計算
3. Extreme rank length (ERL)測度で順序付け
4. $100(1-\alpha)\%$グローバルエンベロープを構築

## 適用例

### 1. フィンランド森林データ
- **データ**: 349本の木（40m×40mプロット）
- **マーク**: 樹冠高/総高さ vs 基部高/総高さ（2部組成）+ 総高さ（絶対値）
- **結果**:
  - 中距離（約6m）で樹冠比率の類似性が期待より高い
  - 近距離（約1m）で対数比積が期待より小さい（競争効果）
  - 混合マーク分析では有意な依存性なし

### 2. スペイン自治体データ
- **データ**: 66自治体（人口1000人以上）
- **マーク**: 工業/建設/商業/サービスの4部門比率
- **結果**:
  - 近距離で組成の類似性
  - 中距離で分散増加
  - サービス部門が空間的異質性に最も寄与

## 限界・残された課題

### 著者が述べる限界
- 構造的ゼロ（一部成分が本質的にゼロ）の扱いは補遺で議論のみ
- $\alpha$変換など代替変換の理論的裏付けは限定的

### 本研究の観点からの限界

1. **記述統計的アプローチのみ**
   - ベイズ推論フレームワークを提供しない
   - パラメトリックモデルによる推定・予測は範囲外

2. **強度関数のモデリングなし**
   - 点過程の強度関数（どこに点が出現しやすいか）は扱っていない
   - 共変量との関係のモデリングなし

3. **空間的依存性の明示的モデルなし**
   - Gaussian Process等による空間相関の明示的モデリングなし
   - 距離に基づく要約特性の評価のみ

4. **計算効率への言及なし**
   - 大規模データへのスケーラビリティは議論されていない

## 本研究（MMCP）との関係

### 借用する要素

1. **組成値マークの基本概念**
   - $D$-部分組成を単体$\mathbb{S}^D$上で扱う枠組み
   - Aitchison幾何学（内積、距離、変換）の適用

2. **対数比変換の適用**
   - ilr/clr変換によるユークリッド空間への写像
   - 変換後の統計的処理の正当化

### MMCPによる拡張

1. **ベイズ推論フレームワークの追加**
   - Eckardtは要約特性の計算（記述統計）
   - MMCPはパラメータの事後分布推定（ベイズ推論）

2. **強度関数のモデリング**
   - Cox過程による強度関数の階層モデリング
   - 共変量効果の推定

3. **空間的依存性の明示的モデル**
   - Gaussian Process / NNGPによる空間相関
   - 空間変化係数モデル

4. **組成とPólya-Gammaの統合**
   - multinomial logitリンクによる組成のモデリング
   - Pólya-Gamma拡張によるGibbsサンプリング

5. **presence-onlyデータへの対応**
   - Eckardtは完全観測データを前提
   - MMCPはthinningによるpresence-onlyデータの扱い

### 位置づけ

Eckardt et al. (2025)は「組成値マーク付き空間点過程」という概念を確立し、記述的な要約特性を提供した。MMCPはこの概念を継承しつつ、以下を統合する推論フレームワークを構築：
- 点過程の強度モデリング（presence-only対応）
- 空間的依存性の明示的モデル（NNGP）
- 組成のベイズ推論（Pólya-Gamma）

## 引用すべき箇所

### 組成値マークの定義（Section 2）
> "We call a real-valued function $\mathbf{c}: \mathbb{R}^2 \to \mathbb{S}^D$ that assigns to each point $x_i$ a $D$-part composition $\mathbf{c}(x_i) = (c_1(x_i), \ldots, c_D(x_i))$ from the open unit $D$-simplex $\mathbb{S}^D$ a composition-valued mark."

### Aitchison幾何学の正当化（Section 2.1）
> "Our proposed set of different (functional) mark summary characteristics allows to decide on the mark independence assumption and investigate the pairwise dependencies of a new type of marked spatial point process. All developments are formalized through extended test functions, which generalise well-known interpretations to the present context."

### 変換の統一性（Section 2.2）
> "Transforming the composition-valued marks to the Euclidean space, the proposed tools can build on established methods for real-valued marks and can borrow strength from existing computational implementations."

### 組成全体 vs 成分別の分解（Section 3.5）
> "$\gamma_{\mathbf{cc}}(r) = \sum_{j=1}^D \gamma_{jj}^{\text{clr}}(r) = \sum_{j=1}^{D-1} \gamma_{jj}^{\text{ilr}}(r)$
> which allows to evaluate which individual components contribute to the overall mark characteristic."

### 結論部分（Section 6）
> "Combining methodological concepts for compositional data and spatial point processes, this paper introduces a novel class of composition-valued marked spatial point processes."
