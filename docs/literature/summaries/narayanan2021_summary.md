# Narayanan et al. (2021) - 柔軟なマーク付き時空間点過程

## 基本情報
- **タイトル**: Flexible marked spatio-temporal point processes with applications to event sequences from association football
- **著者**: Santhosh Narayanan, Ioannis Kosmidis, Petros Dellaportas
- **ジャーナル**: ArXiv preprint (後にJournal of the Royal Statistical Society Series Cに掲載と推測)
- **年**: 2021 (ArXiv version: October 18, 2022)
- **DOI**: 未記載（ArXiv版）
- **関連論点**: 論点3（presence-onlyデータの扱い）、論点2（組成データの扱い）

## 主な貢献
1. **時間とマークの分解によるHawkes過程の一般化**: Marked Hawkes過程の特性（cross-excitation）をマーク空間にのみ適用し、発生時刻に別のモデルを指定する柔軟な枠組み
2. **共変量駆動のcross-excitation**: チーム能力などの共変量情報を直接モデルに組み込み、イベント間相互作用を説明
3. **Stanによるベイズ推論**: Hamiltonian Monte Carlo（No-U-Turn Sampler）による複雑な事後分布からの効率的サンプリング
4. **実データ応用**: サッカー（プレミアリーグ2013/14シーズン）の50万以上のイベントデータで、チーム能力推定、ホームアドバンテージ定量化、ゴール予測を実現

## 手法の概要

### 問題意識
**従来のMarked Hawkes過程の限界**:
- 発生時刻とマークを結合強度関数 $\lambda^*(t, m)$ で同時モデル化
- Self-excitation → 時間的クラスタリングを強制
- サッカーなど、イベントが時間的にクラスターしない現象には不適

**本研究のアイデア**:
> Restrict the excitation property exclusively to the space of marks, providing the freedom to specify a different model for the occurrence times.

### 尤度の因数分解（Cox 1975の多変量分布分解）

$$
\mathcal{L}(\mathcal{F}_{t_n} \mid \boldsymbol{\zeta}, \boldsymbol{\theta}) = \prod_{i=1}^{n} \{g(t_i \mid \mathcal{F}_{t_{i-1}}; \boldsymbol{\zeta}) f(m_i \mid t_i, \mathcal{F}_{t_{i-1}}; \boldsymbol{\theta})\} \{1 - G(T \mid \mathcal{F}_{t_n}; \boldsymbol{\zeta})\}
$$

**記号**:
- $g(\cdot)$: 発生時刻の条件付き密度関数
- $G(\cdot)$: 発生時刻の条件付き分布関数
- $f(\cdot)$: マークの条件付き確率質量関数
- $\boldsymbol{\zeta}$: 時刻のパラメータ
- $\boldsymbol{\theta}$: マークのパラメータ（$\boldsymbol{\zeta}$と独立）
- $\mathcal{F}_{t_i}$: 時刻 $t_i$ までのfiltration（過去の全履歴）

**最後の項**: $1 - G(T \mid \mathcal{F}_{t_n}; \boldsymbol{\zeta})$ は観測期間 $(0, T)$ 内に次のイベントが起こらない確率
- プロセスが $t_n$ で終了する場合（サッカーの前後半など）は不要

### マークの条件付き分布（Marked Hawkes過程から導出）

Marked Hawkes過程の結合強度関数から、マークの条件付き分布を導出：

$$
f(m_i \mid t_i, \mathcal{F}_{t_{i-1}}; \boldsymbol{\theta}) = \frac{\delta_{m_i} + \sum_{t_j < t_i} \alpha^* e^{-\beta(t_i - t_j)} \gamma_{m_j \to m_i}}{1 + \sum_{t_j < t_i} \alpha^* e^{-\beta(t_i - t_j)}}
$$

**パラメータの解釈**:
- $\delta_m \in (0, 1)$: **背景マーク確率**（background component）
  - イベントが背景過程のみで発生した場合のマーク $m$ の確率
- $\alpha^* \geq 0$: **興奮因子**（excitation factor）
  - 過去イベントの寄与の重み、大きいほど履歴への依存が強い
- $\beta > 0$: **減衰率**（decay rate）
  - 過去イベントの興奮が時間とともに減衰する指数率
- $\gamma_{m_j \to m_i} \in (0, 1)$: **変換率**（conversion rate）
  - マーク $m_j$ のイベントがマーク $m_i$ のイベントを引き起こす確率

**重要な洞察**: この式から、Marked Hawkes過程の $\mu$ と $\epsilon$ は一般的な $g(\cdot)$ の指定では識別不可能
- これらは時間次元の進化を特徴付けるが、マーク系列だけでは推定できない
- 時間とマークの分離により、この問題を回避

### 共変量駆動のcross-excitation

変換率 $\gamma_{m_j \to m}$ をBaseline-category logitでモデル化：

$$
\log\left(\frac{\gamma_{m_j \to m}}{\gamma_{m_j \to M}}\right) = \phi_{m_j \to m} + \boldsymbol{\omega}_m^T \boldsymbol{x} \quad (m = 1, \ldots, M-1)
$$

- $\boldsymbol{x}$: 共変量ベクトル（チーム情報、ホーム/アウェイなど）
- $\boldsymbol{\omega}_m$: 回帰パラメータ
- $\phi_{m_j \to m}$: ベースライン変換率（全共変量が0のとき）

### 時空間拡張

位置情報 $z_i$ を組み込んだ拡張尤度：

$$
\mathcal{L}(\mathcal{F}_{t_n} \mid \boldsymbol{\psi}) = \prod_{i=1}^{n} \{g(t_i \mid \mathcal{F}_{t_{i-1}}; \boldsymbol{\zeta}) h(z_i \mid t_i, \mathcal{F}_{t_{i-1}}; \boldsymbol{\eta}) f(m_i \mid t_i, z_i, \mathcal{F}_{t_{i-1}}; \boldsymbol{\theta})\}
$$

- $h(\cdot)$: 位置の条件付き確率質量/密度関数
- $\boldsymbol{\eta}$: 位置のパラメータ
- Filtration $\mathcal{F}_{t_i}$ は時刻・マーク・位置の全履歴を含む

## サッカーイベントデータへの応用

### データ
- **期間**: プレミアリーグ2013/14シーズン全380試合
- **イベント数**: 50万以上のtouch-ballイベント
- **イベントタイプ**: 22種類（Pass, Shot, Goal, Foul など）→ 30複合イベントに統合
- **属性**: 時刻、位置 $(x, y)$、チーム、選手ID、イベント結果

### モデル仕様

**1. 発生時刻**: ガンマ分布

$$
t_{si} - t_{si-1} \mid m_{si-1}, \boldsymbol{a}, \boldsymbol{b} \sim \text{Gamma}(a_{m_{si-1}}, b_{m_{si-1}})
$$

- 前イベントのマーク $m_{si-1}$ に応じた形状・率パラメータ
- Pass後は短い、Goal後は長いなど、イベントタイプ別の待ち時間を捉える

**2. 位置**: 離散1次マルコフ連鎖

$$
h(z_{si} \mid t_{si}, \mathcal{F}_{st_{si-1}}; \boldsymbol{\eta}) = \eta_{(z_{si-1}, m_{si-1}) \to z_{si}}
$$

- 3ゾーン（守備、ミッドフィールド、攻撃）への分割
- 状態空間: $\{1, \ldots, Z\} \times \{1, \ldots, M\}$（位置×マーク）
- 遷移確率行列 $\boldsymbol{\eta}$ でモデル化

**3. マーク**: 4つのモデル階層

**S$\beta$ (scalar $\beta$)**:

$$
f(m_{si} \mid t_{si}, \mathcal{F}_{st_{si-1}}; \boldsymbol{\theta}) = \frac{\delta_{m_{si}} + \sum_{t_{sj} < t_{si}} e^{\alpha - \beta(t_{si} - t_{sj})} \gamma_{m_{sj} \to m_{si}}}{1 + \sum_{t_{sj} < t_{si}} e^{\alpha - \beta(t_{si} - t_{sj})}}
$$

- 全イベントで共通の減衰率 $\beta$

**V$\beta$ (vector $\beta$)**:

$$
f(m_{si} \mid t_{si}, \mathcal{F}_{st_{si-1}}; \boldsymbol{\theta}) = \frac{\delta_{m_{si}} + \sum_{t_{sj} < t_{si}} e^{\alpha - \beta_{m_{sj}}(t_{si} - t_{sj})} \gamma_{m_{sj} \to m_{si}}}{1 + \sum_{t_{sj} < t_{si}} e^{\alpha - \beta_{m_{sj}}(t_{si} - t_{sj})}}
$$

- イベントタイプごとの減衰率 $\beta_m$

**M$\beta$ (matrix $\beta$)**:

$$
f(m_{si} \mid t_{si}, z_{si}, \mathcal{F}_{st_{si-1}}; \boldsymbol{\theta}) = \frac{\delta_{m_{si} \mid z_{si}} + \sum_{t_{sj} < t_{si}} e^{\alpha - \beta_{m_{sj} \to m_{si} \mid z_{si}}(t_{si} - t_{sj})} \gamma_{m_{sj} \to m_{si} \mid z_{si}}}{\sum_{m=1}^M [\delta_{m \mid z_{si}} + \sum_{t_{sj} < t_{si}} e^{\alpha - \beta_{m_{sj} \to m \mid z_{si}}(t_{si} - t_{sj})} \gamma_{m_{sj} \to m \mid z_{si}}]}
$$

- イベントペアごと・位置ごとの減衰率 $\beta_{m \to m' \mid z}$
- 位置依存の背景確率 $\delta_{m \mid z}$ と変換率 $\gamma_{m \to m' \mid z}$
- 例: Cornerは短期的にPass、長期的にShotを興奮 $\beta_{\text{Corner} \to \text{Pass} \mid 3} > \beta_{\text{Corner} \to \text{Shot} \mid 3}$

**M$\beta$A (matrix $\beta$ with abilities)**:

$$
\log\left(\frac{\gamma_{m_{sj} \to m \mid z}(c)}{\gamma_{m_{sj} \to M \mid z}(c)}\right) = \phi_{m_{sj} \to m \mid z} + \omega_{cm}
$$

- $c$: ボール保持チーム
- $\omega_{cm}$: チーム $c$ のマーク $m$ への変換能力
- チーム情報を直接組み込み

### 事前分布
- ガンマ分布のパラメータ $(a_m, b_m)$: 独立指数事前分布
- 位置遷移確率 $\boldsymbol{\eta}$: ディリクレ事前分布（共役）
- 背景マーク確率 $\boldsymbol{\delta}$: ディリクレ事前分布
- 興奮因子 $\alpha$: $N(0, \sigma_\alpha^2)$
- 減衰率 $\beta$（または $\beta_m, \beta_{m \to m' \mid z}$）: 指数事前分布
- 変換率パラメータ $\phi, \omega$: $N(0, \sigma_\gamma^2)$

### ベイズ推論
**事後分布の分解性**: $\boldsymbol{\zeta}, \boldsymbol{\eta}, \boldsymbol{\theta}$ はパラメータ共有なし
- 位置パラメータ $\boldsymbol{\eta}$: 共役事前分布 → ディリクレ事後分布（解析的）
- 時刻パラメータ $\boldsymbol{a}, \boldsymbol{b}$: Stan HMCでサンプリング
- マークパラメータ $\boldsymbol{\theta}$: Stan HMCでサンプリング

**Metropolis-within-Gibbsの失敗**:
- パラメータ間の強い相関と尤度の平坦性（Hawkes過程の典型的問題）
- 混合が極めて悪い → 計算不可能

**Stan + NUTS**: No-U-Turn Samplerによる自動チューニング
- Warm-up phaseで最適なステップサイズとmass matrixを自動調整
- 複雑な事後分布から効率的サンプリング

### Association Rule Learningによる複雑性削減
**問題**: M$\beta$とM$\beta$Aモデルは $M^2 Z$ 個の減衰率パラメータ + $M(M-1)Z$ 個のベースライン変換率
- 計算困難

**解決策**: データからスクリーニング
- Association rules (Agrawal et al., 1993) に基づく有意なイベント相互作用の抽出
- 有意でないパラメータを事前にゼロ固定または除外

## 実証結果

### モデル比較（out-of-sample log predictive density）
- Baseline models（興奮なし）< S$\beta$ < V$\beta$ < M$\beta$ < M$\beta$A
- 興奮ベースモデルが顕著に優れる
- 位置・チーム能力の組み込みで予測精度向上

### パラメータ解釈

**1. ホームアドバンテージの定量化**:
- ホームチームの変換率パラメータがアウェイより高い
- 例: Home_Pass_S → Home_Pass_Sの変換率 > Away_Pass_S → Away_Pass_Sの変換率

**2. チーム能力ランキング**:
- $\omega_{cm}$ からイベントタイプ別のチーム能力抽出
- Pass成功能力、Shot能力、Clear能力などで各チームをランク付け
- 伝統的な順位表とは異なる、プレースタイルを反映した評価

**3. プレースタイルの差異**:
- 攻撃的チーム: 高い $\gamma_{\text{Pass} \to \text{Shot} \mid 3}$
- ポゼッション重視チーム: 高い $\gamma_{\text{Pass} \to \text{Pass}}$

**4. Cross-excitationの可視化**:
- マーク間の全ペア変換率 $\gamma_{m \to m'}$ をヒートマップで表示
- Corner → Shot, Foul → Throw-in などの強い相互作用を特定

### 予測

**1. リアルタイムイベント確率**:
- 指定時間区間内にGoal、Corner、Foulなどが発生する確率
- シミュレーションベース: モデルからイベント系列を生成

**2. 試合結果予測**:
- 前半終了時点でのモデルから、最終スコアの事後予測分布
- 放送時の視聴体験向上に寄与

### Branching structure recovery
Hawkes過程と同様、隠れた分岐構造を復元可能：
- 各イベントが背景過程から発生したか、過去イベントに励起されたかの確率
- 自己興奮の定量化

## 限界・残された課題

### 著者が述べる限界
1. **計算コスト**: M$\beta$Aモデルは数千パラメータ
   - Association rulesでの削減が必須
   - さらに大規模なデータセット（複数シーズン）では課題

2. **時間不変の仮定**: チーム能力 $\omega_{cm}$ はシーズン通じて一定
   - 実際は試合ごと、シーズン内で変動
   - 動的モデルへの拡張が必要

3. **Separabilityとの違い**: 論文の分解は separable intensity ($\lambda^*(t, m) = \lambda_g^*(t) f^*(m)$) とは異なる
   - Separabilityは $f^*(m \mid t)$ が $t$ に依存しない
   - 本手法は $f(m \mid t, \mathcal{F}_{t_{i-1}})$ が $t$ と履歴に依存

4. **Event data quality**: データクリーニングが不可欠
   - 不可能なイベント系列（Goal後すぐにDribbleなど）
   - 手動アノテーションのエラー

### 追加の課題
1. **他のスポーツへの適用**: ラグビー、ホッケー、バスケなど
   - スポーツ固有の調整が必要

2. **選手レベルのモデリング**: チーム能力だけでなく、選手個人の能力
   - パラメータ数がさらに増大

3. **連続位置の扱い**: 論文では3ゾーンの離散化
   - 連続的な $(x, y)$ 座標のモデル化（空間ガウス過程など）

4. **非Poissonベースラインプロセス**: 背景過程 $\delta_m$ は単純な多項分布
   - より複雑な時間依存背景過程への拡張

5. **因果推論**: パラメータは相関を捉えるが、因果関係は不明
   - 反事実シミュレーション（もしこのPassがなかったら）が困難

## 本研究(MMCP)との関係

### 直接的な関係
1. **マーク付き点過程の枠組み**: 黒曜石遺跡データはマーク付き点過程
   - 点位置: 遺跡の空間座標
   - マーク: 黒曜石組成（産地割合）
   - 時間次元 → 時代（縄文・弥生など）に置き換え可能

2. **時間とマークの分解**: MMCPでも適用可能
   - Stage 1（点位置）: 空間IPP
   - Stage 2（組成）: 条件付き多項分布（Dirichlet回帰）
   - 論文と同様、パラメータ $\boldsymbol{\zeta}$ と $\boldsymbol{\theta}$ を独立にモデル化

3. **共変量駆動のマーク分布**: 黒曜石組成を距離・標高で説明
   - 論文のチーム能力 $\omega_{cm}$ ≈ 産地効果 $\omega_{産地 k}$
   - Baseline-category logitで組成をモデル化

4. **空間変動の取り込み**: M$\beta$モデルの位置依存パラメータ
   - MMCPでも、地域ごとに異なる $\gamma_{産地j \to 産地i \mid 地域z}$
   - 関東平野 vs 中部山岳で産地選好が異なるなど

### 本研究への示唆
- **Cross-excitationの考古学的解釈**:
  - 遺跡間の影響: 近くの遺跡の黒曜石利用が、当該遺跡の利用に影響
  - $\gamma_{産地j \to 産地i}$: 産地 $j$ を利用した遺跡の近くに、産地 $i$ を利用する遺跡が出現する確率
  - 交易ネットワークの再構築

- **背景マーク確率 $\delta_m$ vs 興奮成分**:
  - 背景 = 地理的・環境的要因（距離、標高）のみで決まる組成
  - 興奮 = 社会的・文化的ネットワークによる影響

- **減衰率 $\beta$ の解釈**:
  - 空間版: 距離による影響の減衰
  - 時間版（時代データがあれば）: 世代を超えた文化伝播の減衰

- **チーム能力 → 遺跡属性**:
  - 集落規模、生業形態、社会階層などを共変量として
  - 大規模集落は多様な産地を利用する能力が高い、など

### 借用すべき手法
1. **尤度の因数分解アプローチ**:
   - 点過程と組成過程を分離してモデル化
   - パラメータ推定の計算効率向上

2. **Stan + NUTS実装**:
   - Metropolis-Gibbsでは収束しない複雑モデル
   - HMCの自動チューニングで効率的サンプリング

3. **Association rulesによる複雑性削減**:
   - 黒曜石産地ペア全て（$K^2$個）をモデル化すると爆発
   - データから有意な産地間相互作用のみ抽出

4. **Out-of-sample predictive density**:
   - モデル比較の標準的指標
   - LOO-CVやWAICとの関係

5. **Branching structure**:
   - 各遺跡が「自然発生」か「近隣遺跡の影響」かを確率的に分類
   - 文化伝播のメカニズム解明

### 拡張すべき点
1. **Compositional dataへの適応**:
   - 論文のマークは離散カテゴリ（30種）
   - MMCPのマークは組成（連続、単体制約）
   - Aitchison幾何学との統合が必要

2. **空間的自己相関の明示的モデル**:
   - 論文は時間的依存のみ（Hawkes）
   - MMCPでは空間的依存（GP、CAR）を加える

3. **Preferential samplingとの統合**:
   - Cecconi et al. (2016) の共有成分モデル
   - 発掘調査の偏りを補正

4. **真の時空間モデル**:
   - 論文は時間次元のみ（サッカー試合の時刻）
   - 考古学では空間×時代の2次元
   - 時空間Hawkes過程（Reinhart 2018）

5. **Preferenceパラメータの解釈**:
   - 論文の変換率 $\gamma$ = 考古学の「産地選好」
   - 選好の空間・時間変動をモデル化

## 引用すべき箇所

### 因数分解の動機（Section 4.1, p.9）
> "The key insight in the current work is to derive the specification for the marks $f(\cdot \mid t_i, \mathcal{F}_{t_{i-1}}; \boldsymbol{\theta})$ from the joint conditional intensity function of a marked Hawkes process model, and then to specify a probability density function for the times $g(\cdot \mid \mathcal{F}_{t_{i-1}}; \boldsymbol{\zeta})$ best suited to our application. In this way, we can restrict the characteristic excitation property of marked Hawkes processes exclusively to the modelling of the marks."

**用途**: MMCPで点過程と組成を分離する理論的根拠

### 識別可能性の問題（Section 4.1, p.9）
> "Expression (5) makes it immediately apparent, that the parameters $\mu$ and $\epsilon$ of the marked Hawkes process... are not always identifiable for general specifications of $g(\cdot \mid \mathcal{F}_{t_{i-1}}; \boldsymbol{\zeta})$... Apart from a mathematical fact, this is also rather intuitive, because $\mu$ and $\epsilon$ in (2) characterise the evolution of the Hawkes process in the time dimension and the sequence of marks is not sufficient to identify them."

**用途**: 時間とマークを同時モデル化する際の識別性問題を議論

### Separabilityとの違い（Section 4.1, p.10）
> "We should highlight here that the marked point processes from the factorisation in (4) are generally different to the ones that result by assuming separability of the conditional intensity functions... A separable conditional intensity functions has the form $\lambda^*(t, m) = \lambda_g^*(t) f^*(m)$ and implies that the conditional distribution of the mark does not depend on the occurrence time $t$."

**用途**: 分離可能性と因数分解の違いを明確化

### パラメータ解釈（Section 4.2, p.10）
> "The background mark probability $\delta_m \in (0, 1)$ is the probability an event has a mark $m$ if the event is triggered solely by the background component. The excitation factor $\alpha^* \geq 0$ is a scaling factor applied to the contributions from the previous occurrences to the event mark probability. Large values of $\alpha^*$ indicate a stronger dependence of the process on its history."

**用途**: 背景成分と興奮成分の解釈を考古学文脈に翻訳

### 共変量駆動の柔軟性（Section 4.3, p.10）
> "The conditional distribution of marks with probability mass function (5), allows to drive the cross-excitation of the marks using covariates. The conversion rates $\gamma_{m_j \to m}$ can be linked to a covariate vector $\boldsymbol{x}$ observed at the current time through the baseline-category logit specification."

**用途**: 共変量（距離、標高など）で組成を説明する方法

### 時空間拡張（Section 4.4, p.11）
> "We can readily extend the factorisation of the likelihood in (4) to include conditional densities for the event locations, when the latter are observed."

**用途**: 空間位置を明示的にモデル化する拡張

### Branching structure（Section 4.1, p.10）
> "The proposed marked point process model also allows the recovery of the hidden branching structure of the process, a key feature of Hawkes Processes. In Section 6.8, we calculate the branching structure probabilities and quantify the relative contributions of the background process and previous occurrences to the triggering of a new event."

**用途**: 文化伝播vs環境決定論の定量化

### Metropolis-Gibbsの失敗（Section 5.4, p.14）
> "We have also implemented posterior sampling using a Metropolis-within-Gibbs procedure, which, though, proved to mix poorly in artificial data sets... rendering it computationally infeasible. As in the case of Hawkes process, the poor mixing stems from the presence of strong correlations between the model parameters as well as the flatness of the likelihood function."

**用途**: MCMC実装の困難性とHMCの必要性を説明

### Stan + NUTSの利点（Section 5.4, p.14）
> "Stan, on the other hand, implements the No-U-Turn Sampler that automatically calibrates tuning parameters in a warm-up phase and can efficiently sample from complex posterior distributions."

**用途**: Stanを選択する理由を説明

### 実データ応用の意義（Introduction, p.2）
> "The framework is used here for the modelling of in-game event sequences from association football, resulting not only in inferences about previously unquantified characteristics of game dynamics and extraction of event-specific team abilities, but also in predictions for the occurrence of events of interest."

**用途**: モデルの実践的有用性（説明＋予測の両立）を強調
