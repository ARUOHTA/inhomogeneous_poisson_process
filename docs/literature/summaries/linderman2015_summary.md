# Linderman et al. (2015) - Stick-breaking Multinomial with Pólya-Gamma

## 基本情報
- **タイトル**: Dependent Multinomial Models Made Easy: Stick Breaking with the Pólya-Gamma Augmentation
- **著者**: Scott W. Linderman, Matthew J. Johnson, Ryan P. Adams（*共同第一著者）
- **所属**: Harvard University, Twitter (Adams)
- **年**: 2015
- **会議/ジャーナル**: NIPS 2015
- **関連論点**: 論点2（組成データ - 技術的基盤）

## 主な貢献

1. **Stick-breaking表現による多項分布のPólya-Gamma拡張**
   - $K$次元多項分布を$K-1$個の二項分布の再帰的な積で表現
   - Polson et al. (2013)の単一サイト更新から**ブロック更新**へ
   - 潜在変数$\boldsymbol{\omega}$を条件とすると$\boldsymbol{\psi}$全体がGaussian尤度

2. **Correlated Multinomial Modelsの効率的推論**
   - 多項パラメータ間の依存関係をGaussian事前分布でモデリング
   - Gaussian Process, 線形動的システムとの自然な統合
   - 既存のGaussianモデリングソフトウェアを活用可能

3. **3つの応用事例**
   - Correlated Topic Model (CTM)の効率化
   - 出生名データの時空間モデリング
   - 離散値時系列のための連続状態空間モデル

4. **ソフトウェア実装の公開**
   - github.com/HIPS/pgmult
   - モジュール式・拡張可能な設計

## 手法の概要

### Stick-breaking表現による多項分布の分解

**標準的な多項分布**:
$$
\text{Mult}(\mathbf{x} | N, \boldsymbol{\pi}), \quad \sum_{k=1}^K x_k = N, \quad \sum_{k=1}^K \pi_k = 1
$$

**Stick-breaking分解**（$K-1$個の二項分布の積）:
$$
\text{Mult}(\mathbf{x} | N, \boldsymbol{\pi}) = \prod_{k=1}^{K-1} \text{Bin}(x_k | N_k, \tilde{\pi}_k)
$$
ここで：
$$
\begin{aligned}
N_k &= N - \sum_{j<k} x_j, \quad N_1 = N \\
\tilde{\pi}_k &= \frac{\pi_k}{1 - \sum_{j<k} \pi_j}, \quad \tilde{\pi}_1 = \pi_1
\end{aligned}
$$

**解釈**: 各$\tilde{\pi}_k$は残りの確率質量（remaining probability mass）のうち、$k$番目の成分に割り当てられる割合。

### Logisticリンクとstick-breaking写像

**Logistic関数**: $\tilde{\pi}_k = \sigma(\psi_k)$、$\sigma(z) = 1/(1+e^{-z})$

**Stick-breaking写像**: $\boldsymbol{\pi}_{\text{SB}}: \mathbb{R}^{K-1} \to [0, 1]^K$

$$
\begin{aligned}
\pi_1 &= \sigma(\psi_1) \\
\pi_2 &= (1 - \sigma(\psi_1))\sigma(\psi_2) \\
&\vdots \\
\pi_k &= \left(\prod_{j=1}^{k-1}(1 - \sigma(\psi_j))\right)\sigma(\psi_k) \\
\pi_K &= \prod_{j=1}^{K-1}(1 - \sigma(\psi_j))
\end{aligned}
$$

### Pólya-Gamma augmentation

**ロジスティック形式での多項尤度**:
$$
\text{Mult}(\mathbf{x} | N, \boldsymbol{\psi}) = \prod_{k=1}^{K-1} \binom{N_k}{x_k} \frac{(e^{\psi_k})^{x_k}}{(1+e^{\psi_k})^{N_k}}
$$

**Polson et al. (2013)の積分恒等式** (1)を適用：
$$
\frac{(e^\psi)^a}{(1+e^\psi)^b} = 2^{-b} e^{\kappa\psi} \int_0^\infty e^{-\omega\psi^2/2} p(\omega | b, 0) d\omega, \quad \kappa = a - b/2
$$

**拡張モデル**: $a_k(x) = x_k$、$b_k(x) = N_k$として、各座標$\psi_k$に対応するPólya-Gamma潜在変数$\omega_k$を導入：
$$
p(\mathbf{x}, \boldsymbol{\omega} | \boldsymbol{\psi}) \propto \prod_{k=1}^{K-1} e^{(x_k - N_k/2)\psi_k - \omega_k\psi_k^2/2} \propto \mathcal{N}(\boldsymbol{\Omega}^{-1}\kappa(\mathbf{x}) | \boldsymbol{\psi}, \boldsymbol{\Omega}^{-1})
$$

ここで：
- $\boldsymbol{\Omega} = \text{diag}(\boldsymbol{\omega}) = \text{diag}(\omega_1, \ldots, \omega_{K-1})$
- $\kappa(\mathbf{x}) = \mathbf{x} - N(\mathbf{x})/2 = (x_1 - N_1/2, \ldots, x_{K-1} - N_{K-1}/2)$

**重要な性質**: 潜在変数$\boldsymbol{\omega}$を条件とすると、$\boldsymbol{\psi}$の尤度が**対角Gaussian**。

### Gaussian事前分布との共役性

**事前分布**: $\boldsymbol{\psi} \sim \mathcal{N}(\boldsymbol{\mu}_0, \boldsymbol{\Sigma}_0)$

**完全条件付き分布**:
$$
\boldsymbol{\psi} | \mathbf{x}, \boldsymbol{\omega} \sim \mathcal{N}(\boldsymbol{\mu}_n, \boldsymbol{\Sigma}_n)
$$
$$
\begin{aligned}
\boldsymbol{\Sigma}_n &= (\boldsymbol{\Sigma}_0^{-1} + \boldsymbol{\Omega})^{-1} \\
\boldsymbol{\mu}_n &= \boldsymbol{\Sigma}_n(\boldsymbol{\Sigma}_0^{-1}\boldsymbol{\mu}_0 + \kappa(\mathbf{x}))
\end{aligned}
$$

**Gibbsサンプラー**（2ステップ）:
1. $\boldsymbol{\omega} | \boldsymbol{\psi}, \mathbf{x} \sim \text{PG}(N(\mathbf{x}), \boldsymbol{\psi})$（要素ごとに独立）
2. $\boldsymbol{\psi} | \mathbf{x}, \boldsymbol{\omega} \sim \mathcal{N}(\boldsymbol{\mu}_n, \boldsymbol{\Sigma}_n)$（**ブロック更新**）

**従来手法（Polson et al. 2013のmulti-class logistic）との違い**:
- 従来: $\boldsymbol{\psi}$の各要素を**単一サイト更新**（1つずつサンプリング）
- 本手法: $\boldsymbol{\psi}$全体を**ブロック更新**（一括サンプリング）
- 収束速度の劇的な改善

### 単体上の分布の可視化（Figure 1）

**Gaussian事前分布が単体上に誘導する分布**:

- **相関Gaussian**（左図）: 確率質量が$\pi_1 = \pi_2$軸に集中
- **反相関Gaussian**（中央図）: 確率質量が単体の辺に集中（$\pi_1$が大きいとき$\pi_2$は小さい）
- **等方Gaussian**（右図）: 対称Dirichletに近似

**Appendix A**: Gaussian事前分布が誘導する$\boldsymbol{\pi}$上の密度の閉形式、およびDirichletをモーメントマッチングで近似するdiagonal Gaussianの式

## 応用事例

### 1. Correlated Topic Model (CTM)

**Latent Dirichlet Allocation (LDA)の拡張**:
- LDA [Blei et al. 2003]: トピック間は独立（Dirichlet事前分布）
- CTM [Blei & Lafferty 2007]: トピック間の相関をGaussianで表現
- 従来のCTM: Log-normal多項（LN-CTM）で推論困難

**Stick-breaking CTM (SB-CTM)**:
- トピック分布$\boldsymbol{\pi}_d$（文書$d$のトピック比率）に対してstick-breaking表現
- $\boldsymbol{\psi}_d \sim \mathcal{N}(\boldsymbol{\mu}, \boldsymbol{\Sigma})$（相関構造）
- 効率的なGibbsサンプリング

**結果**（Figure 2）:
- AP Newsコーパス、20 Newsgroupコーパスで評価
- LN-CTMより高いpredictive log-likelihood
- トピック間の相関を発見（例: 政治関連トピック間の正相関、経済と法律トピック間の負相関）

### 2. Gaussian Processes with Multinomial Observations

**モデル**（出生名データの時空間パターン）:
- 観測: 各位置$m$・時刻$t$での多項観測$\mathbf{x}_{mt}$
- 潜在過程: $\boldsymbol{\psi}_{m,k} \sim \text{GP}(\mathbf{0}, K_k(\cdot, \cdot))$（時空間GP）
- 各$k$に対して独立なGP

**推論**:
- Pólya-Gamma潜在変数$\omega_{m,k}$を導入
- $\boldsymbol{\omega}$を条件とすると、$\boldsymbol{\psi}_{:,k}$（$k$番目のトピックの全位置・時刻）の完全条件付き分布:
$$
\boldsymbol{\psi}_{:,k} | \mathbf{x}, \boldsymbol{\omega} \sim \mathcal{N}(\boldsymbol{\mu}_k, \boldsymbol{\Sigma}_k)
$$
- 既存のGP推論ソフトウェアを直接利用可能

### 3. Multinomial Linear Dynamical System (SBM-LDS)

**離散値時系列のための連続状態空間モデル**:

**モデル構造**:
$$
\begin{aligned}
\mathbf{z}_t &= \mathbf{A}\mathbf{z}_{t-1} + \boldsymbol{\epsilon}_t, \quad \boldsymbol{\epsilon}_t \sim \mathcal{N}(\mathbf{0}, \mathbf{Q}) \\
\boldsymbol{\psi}_t &= \mathbf{C}\mathbf{z}_t + \boldsymbol{\nu}_t, \quad \boldsymbol{\nu}_t \sim \mathcal{N}(\mathbf{0}, \mathbf{R}) \\
\mathbf{x}_t &\sim \text{Mult}(N_t, \boldsymbol{\pi}_{\text{SB}}(\boldsymbol{\psi}_t))
\end{aligned}$$

**推論**:
- Pólya-Gamma潜在変数: $\boldsymbol{\omega}_t | \mathbf{z}_t, \mathbf{C}, \mathbf{x}_t \sim \text{PG}(N(\mathbf{x}_t), \mathbf{C}\mathbf{z}_t)$
- 潜在状態$\{\mathbf{z}_t\}$: 既存のLDS推論アルゴリズム（Kalmanフィルタ・スムーザ）を直接適用
- システムパラメータ$(\mathbf{A}, \mathbf{C}, \mathbf{Q}, \mathbf{R})$: 標準的なEMまたはGibbsサンプリング

**適用例**:
- テキスト時系列（Figure 4で予測対数尤度を比較）
- DNA配列モデリング

## 実証結果

### CTMの効率性（Figure 2）

**AP Newsコーパス**:
- SB-CTM: 従来のLN-CTMより高いpredictive log-likelihood
- トピック相関の発見:
  - 正相関: (house, committee, congress, law) ↔ (Bush, Dukakis, president, campaign)
  - 負相関: (percent, year, billion, rate) ↔ (court, case, attorney, judge)

**20 Newsgroupコーパス**:
- SB-CTMが最も高い予測性能
- LDA（相関なし）より明確に優れる

### 時系列モデルの予測性能（Figure 4）

**比較対象**:
- SBM-LDS（提案手法）
- Hidden Markov Model (HMM)
- その他の時系列モデル

**結果**: SBM-LDSが最も高いpredictive log-likelihood

## 限界・残された課題

### 著者が述べる限界

1. **スケーラビリティ**
   - GP・LDSは大規模データで計算コスト増大
   - 近似推論手法（変分ベイズ等）への拡張が必要

2. **Stick-breaking順序の依存性**
   - 理論上は順序に依存（実際の影響は小さい）
   - 最適な順序付けの選択基準は未確立

3. **他のリンク関数との比較**
   - Softmax（log-normal）との理論的比較が不十分
   - 状況に応じた使い分けの指針が必要

### 本研究の観点からの限界

1. **点過程への直接適用**
   - 本論文は多項分布（カウントデータ）が対象
   - 点過程の強度関数への適用は別途考慮が必要

2. **空間過程の扱い**
   - GPは使用しているが、NNGPへの拡張は言及なし
   - 大規模空間データへの対応は今後の課題

3. **Presence-only未対応**
   - 完全観測多項データを前提
   - Thinningモデルとの統合は扱っていない

## 本研究（MMCP）との関係

### 借用する要素

1. **Stick-breaking表現の基本アイデア**
   - 多項分布を二項分布の積で表現
   - Logisticリンクによる確率simplex上への写像
   - $\boldsymbol{\pi}_{\text{SB}}(\boldsymbol{\psi})$の構成

2. **Pólya-Gammaブロック更新**
   - 潜在変数$\boldsymbol{\omega}$を条件とするGaussian尤度
   - $\boldsymbol{\psi}$全体の同時サンプリング
   - 収束速度の改善

3. **GPとの統合**
   - Gaussian事前分布の自然な拡張
   - 時空間依存性のモデリング

### MMCPによる拡張

1. **点過程の強度関数への適用**
   - Lindermanは多項観測（カウントデータ）
   - MMCPは点過程の強度関数を多項logitでモデリング
   - Marked point process with composition-valued marks

2. **NNGPとの統合**
   - LindermanはフルGP
   - MMCPはNNGP（計算効率化）
   - Moreira et al. (2024)が既に実装

3. **Presence-onlyデータへの対応**
   - Lindermanは完全観測
   - MMCPはthinningモデルとの統合
   - Moreira & Gamerman (2022)の枠組みと組み合わせ

4. **組成データへの特化**
   - Lindermanは一般的な多項分布
   - MMCPは組成データ（単体上の値）に特化
   - Aitchison幾何学との統合

### 位置づけ

Linderman et al. (2015)は：
- Polson et al. (2013)の多項分布への拡張
- Stick-breaking表現による効率的なブロック更新
- Gaussian models（GP, LDS）との自然な統合
- 相関多項モデルの実用化

MMCPは：
- このstick-breaking + Pólya-Gamma手法を継承
- 点過程の組成マークへの適用
- NNGPとの統合（計算効率化）
- Presence-onlyデータへの対応

Lindermanの手法は、MMCPにおける組成マークのモデリングと推論の技術的基盤を提供する。

## 引用すべき箇所

### 相関多項モデルの動機（Abstract）
> "Many practical modeling problems involve discrete data that are best represented as draws from multinomial or categorical distributions. [...] In all of these cases, we expect some form of dependency between the draws: [...] These dependencies are not naturally captured by the typical Dirichlet-multinomial formulation."

### 本手法の利点（Introduction）
> "Here, we leverage a logistic stick-breaking representation and recent innovations in Pólya-gamma augmentation to reformulate the multinomial distribution in terms of latent variables with jointly Gaussian likelihoods, enabling us to take advantage of a host of Bayesian inference techniques for Gaussian models with minimal overhead."

### Stick-breaking分解（Section 2.1）
> "First, rewrite the $K$-dimensional multinomial recursively in terms of $K-1$ binomial densities: [...] This decomposition of the multinomial density is a 'stick-breaking' representation where each $\tilde{\pi}_k$ represents the fraction of the remaining probability mass assigned to the $k$-th component."

### Gaussian尤度への変換（Section 2.1）
> "That is, conditioned on $\boldsymbol{\omega}$, the likelihood of $\psi$ under the augmented multinomial model is proportional to a diagonal Gaussian distribution."

### 従来手法との比較（Section 2.1）
> "The Pólya-gamma augmentation can be applied to such models [3, 5], but it only provides single-site Gibbs updating of $\psi$. This paper develops a joint augmentation in the sense that, given the auxiliary variables, the entire vector $\psi$ is resampled as a block in a single Gibbs update."

### CTMでの効率性（Section 3）
> "In this section we show that MCMC sampling in a correlated topic model based on the stick breaking construction (SB-CTM) can be significantly more efficient than sampling in the LN-CTM while maintaining the same integration advantage over EM."

### GPとの統合（Section 4）
> "To perform inference, introduce auxiliary Pólya-gamma variables, $\omega_{m,k}$ for each $\psi_{m,k}$. Conditioned on these variables, the conditional distribution of $\psi_{:,k}$ is [Gaussian], enabling us to use off-the-shelf GP inference algorithms."

### LDSへの拡張（Section 5）
> "The stick-breaking multinomial linear dynamical system (SBM-LDS) generates states via a linear Gaussian dynamical system but generates multinomial observations via the stick-breaking map. [...] the system parameters can similarly be updated using standard algorithms."

### Khan et al. (2017)との関係（Section 5）
> "The stick-breaking transformation used herein was applied to categorical models by Khan et al. [17], but they used local variational bound instead of the Pólya-gamma augmentation. Their promising results corroborate our findings of improved performance using this transformation."
