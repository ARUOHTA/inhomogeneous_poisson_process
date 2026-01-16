# Vecchia (1988) - Vecchia近似による大規模空間データ推定

## 基本情報
- **タイトル**: Estimation and Model Identification for Continuous Spatial Processes
- **著者**: A. V. Vecchia
- **ジャーナル**: Journal of the Royal Statistical Society: Series B (Statistical Methodology), 50(2), 297-312
- **年**: 1988
- **DOI**: 10.1111/j.2517-6161.1988.tb01729.x
- **関連論点**: 論点1（空間統計モデル）

## 主な貢献

**Vecchia近似**として知られる、大規模空間データセットに対する効率的な最尤推定法を提案。

1. **近似尤度関数 $L_m$**: 完全尤度 $L$ の計算複雑度 $O(n^3)$ を $O(nm^2)$ に削減
2. **反復推定手順**: $m=1, 2, \ldots, m'$ と段階的に増やしながらパラメータを推定
3. **モデル選択基準**: Akaike情報量基準 $A_m = \Lambda_m + 2d$ によるモデル識別
4. **有理スペクトル密度クラス**: 楕円異方性を含む柔軟なパラメトリックモデル
5. **非格子データへの応用**: 不規則配置データに対して最も効果的

## 手法の概要

### モデル

空間過程：
$$
Z(x, y) = \mathbf{f}^T(x, y) \boldsymbol{\beta} + \xi(x, y)
$$

観測：
$$
z_i = Z(x_i, y_i) + \eta_i, \quad i = 1, \ldots, n
$$

- $\boldsymbol{\beta}$: 回帰パラメータ
- $\xi(x, y)$: 2次定常ガウシアン過程、共分散 $\Gamma(u, v)$
- $\eta_i \sim N(0, \sigma_\eta^2)$: 測定誤差（i.i.d.）

### 有理スペクトル密度関数

$$
S(\kappa) = \sigma^2 \prod_{j=1}^q |\kappa^2 + \theta_j|^{2n_j} / \prod_{j=1}^p |\kappa^2 + \phi_j|^{2m_j}
$$

楕円異方性：
$$
\kappa^2 = [\lambda^{-1}(k_1 \cos\alpha - k_2 \sin\alpha)]^2 + [\lambda(k_1 \sin\alpha + k_2 \sin\alpha)]^2
$$

- $\lambda$: スケーリングパラメータ（異方性の程度）
- $\alpha$: 回転パラメータ（異方性の方向）
- $\theta_j \in \mathbb{R}$, $\phi_j \in \mathbb{R}^+$

対応する共分散関数：
$$
\Gamma(r) = \sigma^2 (2\pi)^{-1} (-1)^{M-1} \sum_{j=1}^p \frac{1}{(2m_j - 1)!} \frac{\partial^{2m_j-1}}{\partial \phi_j^{2m_j-1}} \{w_j K_0(r\sqrt{\phi_j})\}
$$

ここで $K_0(\cdot)$ は第2種変形Bessel関数、$r$ は異方性スケール後の距離。

### Vecchia近似

**完全尤度**:
$$
L(\mathbf{z}) = \prod_{i=1}^n p(z_i | z_j, 1 \leq j \leq i-1)
$$

**近似尤度（order $m$）**:
$$
L_m(\mathbf{z}) = \prod_{i=1}^n p(z_i | \mathbf{z}_{im})
$$

ここで $\mathbf{z}_{im}$ は $z_1, \ldots, z_{i-1}$ の中でユークリッド距離 $d_{ij} = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}$ が最も小さい $\min(i-1, m)$ 個の観測。

**明示的表現**:
$$
L_m(\mathbf{z}) = (2\pi)^{-n/2} (\sigma^2 \gamma_0)^{-n/2} \prod_{i=1}^n \omega_{im}^{-1/2} \exp\left[ -\frac{1}{2\sigma^2 \gamma_0} \sum_{i=1}^n \omega_{im}^{-1} e_{im}^2 \right]
$$

$$
\omega_{im} = 1 + v^2 - \mathbf{r}_{im}^T (R_{im} + v^2 I)^{-1} \mathbf{r}_{im}
$$
$$
e_{im} = \varepsilon_i - \mathbf{r}_{im}^T (R_{im} + v^2 I)^{-1} \boldsymbol{\varepsilon}_{im}
$$

- $v^2 = \sigma_\eta^2 / \text{var}\{\xi(x,y)\}$: 相対測定誤差分散
- $\mathbf{r}_{im} = \text{corr}(\xi_i, \boldsymbol{\xi}_{im})$
- $R_{im} = \text{corr}(\boldsymbol{\xi}_{im}, \boldsymbol{\xi}_{im})$

### 反復最尤推定

**Step 1**: $m=1$、$\boldsymbol{\beta}$ を最小二乗推定 $\widetilde{\boldsymbol{\beta}}_{\text{OLS}}$ に固定
- $\hat{\boldsymbol{\psi}}_1$ を数値最適化で求める

**Step 2**: $m=2, 3, \ldots$ と増加
- $\hat{\boldsymbol{\psi}}_{m-1}$ を初期値として $\hat{\boldsymbol{\psi}}_m$ を推定

**Step 3**: 収束判定
$$
\Lambda_m = -2 \log[L_m^*(\hat{\boldsymbol{\psi}}_m)]
$$

$\Lambda_m$ が安定する $m = m'$ を選択（通常 $m' \leq 10$）

**Step 4**: $\boldsymbol{\beta}$ の推定
- $m = m'$ で $\boldsymbol{\beta}$ を推定に含める:
$$
\hat{\boldsymbol{\beta}}_m = \left[ \sum_{i=1}^n \omega_{im}^{-1} \mathbf{g}_{im} \mathbf{g}_{im}^T \right]^{-1} \left[ \sum_{i=1}^n \omega_{im}^{-1} \mathbf{g}_{im} h_{im} \right]
$$

### モデル選択

Akaike情報量基準（AIC）:
$$
A_m = \Lambda_m + 2d
$$

ここで $d$ は $\hat{\boldsymbol{\psi}}_m$ の次元。

**手順**:
1. 等方性モデル（$\lambda=1, \alpha=0$）で複数のスペクトル密度を比較
2. 最小の $A_m$ を持つモデルを選択
3. 選択されたモデルに異方性パラメータ $\lambda, \alpha$ を追加してテスト

## 実証結果

### シミュレーション研究（6データセット）

**データセット1, 2**: $S(\kappa) = 100.0/(\kappa^2 + 0.1)^2$, $\lambda^2 = 2.00$, $\alpha = 30°$
- $n=100$, サンプリング領域 $30 \times 30$ および $45 \times 45$
- モデルA（正しい）vs モデルB（$(\kappa^2 + \phi)^{-4}$）
- $A_m$ の差: データセット1で約8、データセット2で約10 → モデルA選択
- 平滑性の違い: モデルAは平均二乗微分不可能、モデルBは2回微分可能

**データセット3, 4**: $S(\kappa) = 4.0(\kappa^2 + 0.0)^2/(\kappa^2 + 0.1)^6$, 等方性
- $n=100$, サンプリング領域 $40 \times 40$ および $60 \times 60$
- MA項 $\theta$ の有意性検定: わずかにモデルA（$\theta \neq \phi$）が優位
- 相関関数の形状に柔軟性を提供

**データセット5, 6**: $S(\kappa) = 1.0(\kappa^2 - 0.021)^2/(\kappa^2 + 0.015)^4$, $\lambda^2 = 2.64$, $\alpha = 60°$
- $n=100$, サンプリング領域 $150 \times 150$ および $225 \times 225$
- 強い異方性を正確に推定: $\hat{\lambda}^2 = 2.58, 5.73$, $\hat{\alpha} = 71°, 66°$

**計算時間**（Prime 9950ミニコンピュータ、データセット1）:
- $m=2$: 0.36秒、$m=8$: 2.11秒、完全尤度: 7.73秒
- 単純指数共分散の場合: $m=8$: 0.99秒、完全尤度: 2.19秒
- $m \leq 10$ で十分な近似、計算量は $n$ に線形

### 実データ: Saratoga Valley地下水データ

- **データ**: Wyoming州サラトガ盆地、93観測井の水位（1980年10-11月）
- **トレンド**: $Z(x,y) = \beta_1 + \beta_2 x + \beta_3 y + \xi(x,y)$
- **OLS推定**: $\widetilde{\beta}_1 + \widetilde{\beta}_2 x + \widetilde{\beta}_3 y = 2227.6 - 0.34x - 2.92y$
- **残差範囲**: 約160メートル、有意な非ランダム性

**モデル選択**:
- 4つの等方性モデルを比較（Table 6）
- 最良モデル: $S(\kappa) = \sigma^2(\kappa^2 + \theta)^2 / (\kappa^2 + \phi)^6$
- $m=6$ で安定、$A_6 = 652.32$

**最終推定**（$m=6$, $\boldsymbol{\beta}$ 含む）:
- $\hat{\beta}_1 = 2253.69$, $\hat{\beta}_2 = -0.67$, $\hat{\beta}_3 = -3.07$
- $\hat{\theta} = -0.212$, $\hat{\phi} = 0.182$
- $\hat{\sigma}\hat{\gamma}_0^{1/2} = 39.12$ m（空間過程の標準偏差）
- $\hat{v}\hat{\sigma}\hat{\gamma}_0^{1/2} = 4.39$ m（測定誤差の標準偏差）

**異方性テスト**: $\lambda, \alpha$ の追加は有意でない → 等方性モデルが適切

## 限界・残された課題

### 著者が述べる限界

1. **異方性の制限**: 回転・スケーリングのみの楕円異方性
   - より複雑な異方性形状には対応できない

2. **線形平均関数の仮定**: $\mathbf{f}^T(x,y)\boldsymbol{\beta}$ の特定が必要
   - ノンパラメトリックアプローチ（median polish等）との併用を示唆

3. **条件付け集合の選択**: ユークリッド距離による最近傍
   - 真の異方性が大きい場合（$\lambda \gg 1$）、共分散楕円内の点が最適だが実用上困難

4. **ガウス性の仮定**: $\xi(x,y)$ と $\eta_i$ が正規分布
   - ロバスト推定への拡張は今後の課題

### 私たちが認識する限界

5. **順序依存性**: $L_m$ は観測の順序に依存（小さい $m$ で顕著）
   - 推奨: $x$ または $y$ 座標の増加順に並べる

6. **点過程データには不適**: 固定サンプリング点での観測を前提
   - presence-onlyデータやポイントパターンには直接適用できない

7. **計算量の増加**: $m$ を大きくすると $O(nm^2)$ が支配的
   - 実用上 $m \leq 10$ が推奨されるが、強い空間相関では不十分な可能性

8. **スペクトル密度クラスの制限**: 実数パラメータのみ（複素数は排除）
   - より柔軟な相関構造への拡張は計算上困難

9. **収束判定の主観性**: $\Lambda_m$ の安定性を目視判断
   - 形式的な停止規則がない

10. **空間共変量の欠如**: 実データ例（Saratoga）は座標のみのトレンド
    - 環境共変量を含む拡張が自然だが議論されていない

## 本研究（MMCP）との関係

### 直接的な関連

1. **計算効率化の枠組み**: MMCPも大規模データへの対応が必須
   - Vecchia近似は $n \times n$ 共分散行列の逆行列・行列式計算を回避
   - 現代の空間統計ソフトウェア（GPyTorch、R-INLA等）でVecchia近似が標準実装

2. **近似戦略の共通性**: 条件付き独立性の活用
   - Vecchia: $z_i \perp z_{\{1,\ldots,i-1\} \setminus \text{nearest } m} \mid z_{\text{nearest } m}$
   - MMCP: 近似推論（Laplace近似、変分推論等）も条件付き構造を利用

3. **反復推定の哲学**: 粗い近似から始めて精緻化
   - $m=1 \to 2 \to \ldots \to m'$ の反復
   - MMCPのwarm startやhierarchical fittingに通じる

### MMCPでの応用可能性

4. **潜在ガウシアン過程の推論**: MMCP第1層（トブラー距離場）
   - MMCPの $\log\lambda(\mathbf{x}) = \mu + Y(\mathbf{x})$ で $Y$ がGP
   - Vecchia近似により $n$ 個の遺跡位置での $Y$ の事後分布を効率計算

5. **マーク過程の条件付きモデル**: 産地マーク $\mathbf{w}$ の推論
   - $p(\mathbf{w}_i | \mathbf{w}_{\text{neighbors}}, \mathbf{s}_i, \boldsymbol{\theta})$ の近似
   - ただし、マークは独立でないため直接適用は困難

6. **異方性の扱い**: トブラー距離の方向性
   - トブラー関数 $T(\theta) = a \exp(b \cdot |\tan\theta|)$ は方向依存
   - Vecchiaの楕円異方性パラメータ $(\lambda, \alpha)$ が類似の役割

### MMCPで借用すべき要素

7. **モデル選択基準 $A_m$**:
   - MMCPでも複数の相関関数（指数、Matérn、有理二次等）を比較
   - AICベースの選択は解釈可能で実装容易

8. **反復診断プロット**:
   - Figure 1: 真の相関関数 vs 推定相関関数
   - Figure 2: 2次元相関の等高線プロット
   - MMCPの事後チェックに有用（pair correlation function $g(r)$ の比較）

9. **計算時間のトレードオフ分析**:
   - Table 5のような $m$ vs 計算時間の定量化
   - MMCPでも近似精度とコストのバランスが重要

### MMCPで拡張すべき点

10. **点過程強度への適用**: Vecchiaは格子データ前提
    - MMCPはイベント位置 $\{\mathbf{s}_i\}_{i=1}^N$ がランダム
    - thinning algorithmやspatial birth-death processとの統合が必要

11. **マルチスピーシーズへの一般化**: Vecchiaは単変量
    - MMCPは産地マーク $\mathbf{w}_i \in \mathbb{S}^{d-1}$（組成データ）
    - 多変量ガウシアン過程への拡張（cross-covariance行列）

12. **非ガウシアン観測モデル**: Vecchiaは測定誤差のみ
    - MMCPはポアソン点過程（カウント分布）+ compositional likelihood
    - Diggle et al. (1998)の一般化線形予測とVecchia近似の統合が理想

### 技術的な注意点

13. **順序付けの重要性**: $d_{ij}$ と $i-j$ の相関
    - MMCPの空間データは通常 $(x, y)$ でソート可能
    - 3次元データ（標高含む）では最適順序が非自明

14. **小さい $m$ での不安定性**: Table 7の $m=1,2$ で推定が大きく変動
    - MMCPでも初期フェーズ（$m < 5$）は慎重に扱う

15. **完全尤度との一致**: $m \geq 10$ で $L_m \approx L$
    - データ密度・相関範囲による調整が必要

## 引用すべき箇所

### Vecchia近似の定義

> "Formally selecting which observations to include in $z_{im}$ according to some criterion such as maximising the multiple correlation coefficient between $z_i$ and $z_{im}$ is not feasible because $z_{im}$ would then depend on both the form of the model (4) and on specific parameter values. ... The only logical non-model-dependent method for choosing $z_{im}$ is to select those observations whose sample locations are closest to $z_i$ in some sense." (p. 303)

### 近似精度の経験的知見

> "Our experience from the analysis of numerous simulated and actual data sets is that there is generally a small value of $m$, say $m'$, such that fluctuations in $\Lambda_m$ become negligible for $m > m'$." (p. 306)

### 計算効率の利点

> "There should be no problems in computing the approximate likelihood functions for very large data sets as long as $m \leq 10$ is adequate, because the computation times and storage requirements involved in computing $-2\log \hat{L}_m$ increase in direct proportion to $n$." (p. 312)

### モデル識別の有効性（結論）

> "For most reasonable sampling schemes, virtually all the information necessary for estimation and model identification is contained in the approximate likelihood functions $L_m$ for $1 \leq m \leq 10$." (p. 319)

### 異方性推定のロバスト性

> "The model identification methods are robust to misspecification of anisotropy, thus allowing that the particular form of spectral density be selected assuming isotropy." (p. 319)

### 平滑性の識別

> "The iterative estimation statistics are effective for identifying the smoothness properties of a continuous spatial process based on a sparse set of observations at point locations." (p. 319)

### 有理スペクトル密度と微分可能性

> "It is straightforward to show that a process with spectral density function (4) is mean square differentiable up to order $(\Sigma 2m_j - \Sigma 2n_j - 2)$." (p. 310)

### 測定誤差の組み込み

> "With ground-water data of this type, there is generally a significant measurement error. Hence, the observations are assumed to be given by equation (2) and the measurement error parameter $v > 0$ is included in the estimation." (p. 315)

### 共分散関数の計算

> "To compute the process variance, $\Gamma(0)$, $K_0(r\sqrt{\phi_j})$ needs to be replaced by $-\log(\sqrt{\phi_j})$ (Vecchia (1985), proposition 3)." (p. 299)

### 反復推定の実用的アドバイス

> "$\beta$ is held fixed at $\widetilde{\beta}$ throughout the process of selecting $m'$, after which $\beta$ may be included in the estimation, if desired, to obtain $\widehat{\psi}_{m'}$ and $\widehat{\beta}_{m'}$. Furthermore, anisotropy need be included in the model only after selecting $m'$. This results in considerable savings in computation time while usually having a negligible effect on the selection process." (pp. 306-307)
