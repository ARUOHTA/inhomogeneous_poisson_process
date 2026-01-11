# マーク付き点過程モデル：理解と開発メモ

*作成日: 2026-01-05*

---

## 1. モデル概要

マーク付き点過程モデルは、以下の2つのコンポーネントを統合したベイズ階層モデルである。

### 1.1 点過程部分（Intensity Process）

遺跡の空間分布をモデル化する。

- **遺跡存在確率**: $q(s) = \text{sigmoid}(\eta_{\text{int}})$
- **強度関数**: $\lambda(s) = \lambda^* \cdot q(s)$
- **偽不在点**: $U \sim \text{IPP}(\lambda^*(1-q))$ via Poisson thinning

線形予測子:
$$\eta_{\text{int}}(s) = W_{\text{int}}(s)^\top \boldsymbol{\beta}_{\text{int}} + u_{\text{int}}(s)$$

### 1.2 マーク部分（Composition Process）

各遺跡における黒曜石の産地構成比をモデル化する。

- **産地構成比**: $\pi_k(s) = \text{softmax}(\eta_k)$
- **観測**: $y_i \sim \text{Multinomial}(N_i, \pi)$

線形予測子:
$$\eta_k(s) = W_z(s)^\top \boldsymbol{\beta}_k + u_k(s)$$

### 1.3 空間効果

両コンポーネントとも、NNGP（Nearest Neighbor Gaussian Process）による空間相関をモデル化。

---

## 2. ギブスサンプリングアルゴリズム

sec7.tex に基づく6ステップのギブスサンプラー:

### ステップ (a): 偽不在のサンプリング
$$U^{(\tau)} \sim \text{IPP}\bigl(\lambda^{*(\tau-1)}(1-q^{(\tau-1)})\bigr)$$

Poisson thinning アルゴリズム:
1. $N \sim \text{Poisson}(\lambda^* |\mathcal{D}|)$ をサンプル
2. $j=1,...,N$ について $(s_j, t_j) \sim \text{Uniform}(\mathcal{D})$
3. $u_j < 1-q(s_j,t_j)$ なら採用

### ステップ (b): 点過程側 PG 変数のサンプリング
$$\omega_i^{(\tau)} \sim \text{PG}(1, \eta_{\text{int},i}^{(\tau-1)})$$

**重要**: 点過程では $b=1$ を使用。

### ステップ (c): 点過程パラメータの更新
$$\theta_{\text{int}}^{(\tau)} = \begin{pmatrix} \boldsymbol{\beta}_{\text{int}}^{(\tau)} \\ \boldsymbol{u}_{\text{int}}^{(\tau)} \end{pmatrix} \sim \mathcal{N}(\boldsymbol{m}_{\text{int}}, V_{\text{int}})$$

完全条件付き分布:
- $V_{\text{int}}^{-1} = H^\top \Omega H + R_0$
- $\boldsymbol{m}_{\text{int}} = V_{\text{int}}(H^\top \Omega \boldsymbol{z} + R_0 \mu_0)$

ここで $\kappa_i = y_i - 1/2$, $z_i = \kappa_i / \omega_i$。

### ステップ (d): λ* の更新
$$\lambda^{*(\tau)} \sim \text{Gamma}(m_0 + n^{(\tau)}, r_0 + |\mathcal{D}|)$$

ここで $n^{(\tau)} = n_X + n_U^{(\tau)}$。

### ステップ (e): マーク側 PG 変数のサンプリング
$$\xi_{ik}^{(\tau)} \sim \text{PG}(N_i, \eta_{ik}^{(\tau-1)})$$

**重要**: マークでは $b = N_i$（総出土数）を使用。

### ステップ (f): マークパラメータの更新
$$\theta_k^{(\tau)} = \begin{pmatrix} \boldsymbol{\beta}_k^{(\tau)} \\ \boldsymbol{u}_k^{(\tau)} \end{pmatrix} \sim \mathcal{N}(\boldsymbol{m}_k, V_k)$$

完全条件付き分布（カテゴリ $k$ ごと）:
- $V_k^{-1} = H_z^\top \Omega_k H_z + R_{0,k}$
- $\boldsymbol{m}_k = V_k(H_z^\top \Omega_k \boldsymbol{z}_k + R_{0,k} \mu_{0,k})$

ここで $\tilde{\kappa}_{ik} = y_{ik} - N_i/2$, $z_{k,i} = \tilde{\kappa}_{ik} / \xi_{ik}$。

---

## 3. 実装構造

### 3.1 モジュール構成

```
bayesian_statistics/nngp/model/marked_point_process/
├── __init__.py          # 公開API
├── config.py            # MarkedPointProcessConfig
├── dataset.py           # MarkedPointProcessDataset, prepare_marked_point_process_dataset
├── intensity.py         # 点過程コンポーネント
├── composition.py       # マークコンポーネント
└── sampler.py           # MarkedPointProcessSampler, MarkedPointProcessResults
```

### 3.2 各モジュールの責務

| モジュール | 責務 |
|-----------|------|
| `config.py` | MCMC設定、カーネルパラメータ、距離事前分布パラメータ |
| `dataset.py` | データ準備、設計行列構築、距離特徴量計算 |
| `intensity.py` | `compute_q`, `sample_pseudo_absence`, `update_lambda_star`, `update_beta_intensity` |
| `composition.py` | `compute_kappa_tilde_mark`, `sample_xi_mark`, `softmax_with_baseline` |
| `sampler.py` | 6ステップのギブスサンプラー統合、結果クラス |

### 3.3 主要クラス

#### MarkedPointProcessConfig (`config.py:10`)
- MCMC設定: `n_iter`, `burn_in`, `thinning`, `seed`
- NNGP設定: `neighbor_count`
- カーネルパラメータ: `intensity_kernel_*`, `mark_kernel_*`
- λ* 事前分布: `lambda_prior_shape`, `lambda_prior_rate`
- 距離事前分布: `tau`, `alpha`, `source_weights`, `lambda_fixed`

#### MarkedPointProcessDataset (`dataset.py:17`)
- 遺跡データ: `site_coords`, `counts`, `total_counts`
- 設計行列: `design_matrix_intensity`, `design_matrix_marks`
- グリッドデータ: `grid_coords`, `valid_grids`
- 距離事前分布: `distance_features_*`, `prior_mean_intercept_*`

#### MarkedPointProcessSampler (`sampler.py:369`)
- 初期化: NNGP因子構築、PGサンプラー準備
- `_sample_intensity()`: ステップ (a)-(d)
- `_sample_marks()`: ステップ (e)-(f)
- `run()`: MCMCループ実行

#### MarkedPointProcessResults (`sampler.py:53`)
- 事後サンプル: `lambda_star_samples`, `beta_mark_samples`, `beta_int_samples`
- 予測メソッド: `predict_probabilities()`, `predict_intensity()`
- 効果分解: `decompose_effects()`

---

## 4. 距離事前分布

### 4.1 概念

データがない領域では産地からの距離に基づくベースライン期待値を使用し、データがある領域ではデータに引っ張られる。

### 4.2 数式

距離ベース log-ratio 特徴量:
$$g_k(s) = \log(p_{0k}(s)) - \log(p_{0K}(s))$$

ここで $p_{0k}$ は重み付き逆ソフトマックス:
$$p_{0k}(s) = \frac{w_k \cdot \exp(-Z_k(s)/\tau)}{\sum_{k'} w_{k'} \cdot \exp(-Z_{k'}(s)/\tau)}$$

事前平均:
$$\mu_k(s) = \lambda_k \cdot g_k(s)$$

### 4.3 ハイパーパラメータ

| パラメータ | 意味 | 典型値 |
|-----------|------|--------|
| `tau` | 温度（小さいほど距離による差が大きい） | 0.5 |
| `alpha` | 重要度の指数 | 1.0 |
| `source_weights` | 産地の重要度 | [2, 1, 0.05, 0.05] |
| `lambda_fixed` | λの固定値 | [1, 1, 1, 1] |

---

## 5. 点過程とマークの評価点の違い（重要）

```
┌─────────────────────────────────────────────────────────────────┐
│ 【点過程部分】                                                    │
│   評価点: X ∪ U = {(s_i, t_i) : i=1,...,n}                       │
│   - n = n_X + n_U （反復ごとにn_Uは変化）                         │
│   - 観測遺跡 X: y_i = 1                                          │
│   - 偽不在 U: y_i = 0                                            │
│   - NNGP因子は n 点で構築（毎反復で再構築）                       │
├─────────────────────────────────────────────────────────────────┤
│ 【マーク部分】                                                    │
│   評価点: X のみ = {(s_i, t_i) : i=1,...,n_X}                     │
│   - 偽不在 U にはマークデータ（産地出土数）がない                  │
│   - NNGP因子は n_X 点で構築（MCMC開始前に1回）                    │
└─────────────────────────────────────────────────────────────────┘
```

---

## 6. 現在の実装状況

### 6.1 実装済み

- [x] 設定クラス (`config.py`)
- [x] データセットクラス (`dataset.py`)
- [x] 点過程コンポーネント (`intensity.py`)
- [x] マークコンポーネント (`composition.py`)
- [x] 統合ギブスサンプラー (`sampler.py`)
- [x] 結果クラス・予測メソッド
- [x] 距離事前分布
- [x] 効果分解 (`decompose_effects`)

### 6.2 ノートブック

- `29_marked_point_process_demo.ipynb`: 基本デモ
- `29_marked_point_process_distance_prior.ipynb`: 距離事前分布付き（現在開発中）

---

## 7. 今後の作業予定

### 2026-01-05: バグ修正 - グリッド投影における事前平均調整

**問題**: `predict_probabilities(location="grid")` と `decompose_effects(location="grid")` の結果が異なっていた。

**原因**: `decompose_effects` には距離事前分布の事前平均調整があったが、`predict_probabilities` にはなかった。

**修正内容** (`sampler.py`):

1. **共通ヘルパーメソッド `_project_beta_to_grid()` を追加** (82-126行目)
   - サイトのβをグリッドにNNGP条件付き投影
   - 切片（j=0）の場合、事前平均調整を自動適用
   - 距離事前分布がない場合は調整をスキップ

2. **`predict_probabilities()` を修正** (157-200行目)
   - `sample_conditional=False`: ヘルパーメソッドを使用
   - `sample_conditional=True`: 各サンプルの投影時に事前平均調整を追加

3. **`decompose_effects()` を修正** (337-348行目)
   - 重複していたグリッド投影ロジックをヘルパーメソッドに置き換え

**効果**:
- DRY原則に従い、グリッド投影ロジックを一箇所に集約
- `predict_probabilities` と `decompose_effects["full"]` が一致するようになった
- 保守性・可読性の向上

**検証**: ノートブック29で動作確認済み ✅

---

### 次のステップ候補

#### A. ハイパーパラメータチューニング

| パラメータ | 現在値 | 影響 |
|-----------|-------|------|
| `mark_lengthscale` | 0.2 | 空間的な滑らかさ |
| `mark_variance` | 0.1 | 効果の大きさ |
| `tau` | 0.5 | 距離事前分布の温度 |
| `source_weights` | [2,1,0.05,0.05] | 産地の重要度 |
| `lambda_fixed` | [1,1,1,1] | 距離効果の強さ |

方法案:
- グリッドサーチ + LOOCV
- 事後予測チェック（posterior predictive check）
- 感度分析

#### B. モデルの考察・改善

1. **距離事前分布の妥当性**
   - 現在の逆ソフトマックス形式は適切か？
   - `tau`, `alpha` の解釈と設定根拠

2. **空間効果の構造**
   - 切片以外の共変量（標高など）の追加効果
   - カーネルパラメータの特徴量ごとの設定

3. **点過程部分との統合**
   - 遺跡存在確率と産地構成比の関係
   - 偽不在点の影響

4. **評価指標**
   - LOOCV以外の評価方法
   - 時期間の比較

---

## 8. データパイプライン

### 8.1 最終的に必要なファイル（`data/` ディレクトリ）

`ObsidianDataPreprocessor.load_data()` が読み込むファイル:

| ファイル名 | 内容 |
|-----------|------|
| `11_gdf_elevation.csv` | 標高データ（グリッド）+ 産地へのコスト距離 |
| `11_gdf_obsidian.csv` | 黒曜石出土データ（個票） |
| `11_gdf_sites.csv` | 遺跡データ（メタ情報） |

### 8.2 データフロー図

```
【原型データ】                           【中間ファイル】                    【最終ファイル】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

obsidian_data/obsidian.csv ─────┐
                                │
data/elevation_data/            │     09_obsidian_preprocess_1.ipynb
  G04-d-11_*-jgd_GML/ ──────────┼───────────┐
                                │           │
data/lake_data/                 │           ├──→ 9_gdf_elevation.csv
  W09-05_GML/ ──────────────────┤           │           │
                                │           │           │
data/river_node_data/           │           │           │ 10_1_obsidian_preprocess_2.ipynb
  W05-*_GML/ ───────────────────┤           │           ├────────────────────────────────┐
                                │           │           │                                │
data/river_stream_data/         │           │           ▼                                │
  W05-*_GML/ ───────────────────┴───────────┼──→ 9_gdf_river_stream.csv                  │
                                            │           │                                │
                                            │           │ 10_2_calculate_distance.py     │
                                            │           │           │                    │
                                            │           ▼           ▼                    │
                                            │    10_1_gdf_elevation_tobler.csv ◄─────────┘
                                            │                       │
                                            │                       ▼
                                            │    10_2_gdf_elevation_with_costs.csv
                                            │                       │
                                            │                       │ 11_obsidian_preprocess_3.ipynb
                                            ▼                       ▼
                                    9_obsidian_gdf.csv ─────────────┬──→ 11_gdf_elevation.csv
                                                                    ├──→ 11_gdf_obsidian.csv
                                                                    └──→ 11_gdf_sites.csv
```

### 8.3 原型データ（サーバーから手動で取得）

| データ | パス | 説明 |
|-------|------|------|
| 黒曜石データ | `obsidian_data/obsidian.csv` | 黒曜石出土記録の個票データ |
| 標高データ | `data/elevation_data/G04-d-11_*-jgd_GML/` | 国土地理院5次メッシュ標高 |
| 湖沼データ | `data/lake_data/W09-05_GML/` | 国土数値情報 湖沼データ |
| 河川ノード | `data/river_node_data/W05-*_GML/` | 国土数値情報 河川データ |
| 河川流路 | `data/river_stream_data/W05-*_GML/` | 国土数値情報 河川データ |

**必要なメッシュコード（標高データ）**:
```
5138, 5139, 5238, 5239, 5240, 5338, 5339, 5340, 5438, 5439, 5440, 5538, 5539, 5540, 5541
```

**必要な都道府県コード（河川データ）**:
```
07_07, 08_08, 08_09, 08_10, 08_11, 08_12, 08_13, 08_14, 07_15, 08_19, 08_20, 08_22
```

### 8.4 前処理パイプライン実行順序

```bash
# 1. 原型データを配置後、前処理を実行
cd notebooks/

# Step 1: 基本データ作成
jupyter nbconvert --execute 09_obsidian_preprocess_1.ipynb

# Step 2: Tobler歩行速度計算
jupyter nbconvert --execute 10_1_obsidian_preprocess_2.ipynb

# Step 3: 産地へのコスト距離計算（時間がかかる）
cd ../bayesian_statistics/
python 10_2_calculate_distance.py

# Step 4: 最終データ作成
cd ../notebooks/
jupyter nbconvert --execute 11_obsidian_preprocess_3.ipynb
```

---

## 参考文献

- 理論: `docs/251220_marked_markdown/` (sec2-sec9)
- 実装計画: `docs/implementation/MARKED_PP_IMPLEMENTATION_PLAN.md`
- コードベース分析: `docs/implementation/CODEBASE_ANALYSIS.md`
