# Experiments Module

修士論文 第5章の実験を実行するためのモジュール。

## クイックスタート

```bash
# 全実験を実行（モデル推定 → LOOCV → 図生成）
uv run python -m bayesian_statistics.experiments.run_all --all

# ヘルプを表示
uv run python -m bayesian_statistics.experiments.run_all --help
```

## コマンドラインオプション

| オプション | 説明 |
|-----------|------|
| `--run-models` | MMCP・NWモデルを全5期間で実行し、結果を保存 |
| `--mmcp-only` | MMCPモデルのみ実行（NWをスキップ） |
| `--run-loocv` | Leave-One-Out交差検証を実行 |
| `--generate-figures` | 保存済み結果から図を生成 |
| `--all` | 上記3つを順番に実行 |
| `--n-iter N` | MCMCイテレーション数（デフォルト: 1000） |
| `-q, --quiet` | 進捗バー・出力を抑制 |
| `-v, --verbose` | 詳細出力モード |

## 使用例

### モデル推定のみ

```bash
uv run python -m bayesian_statistics.experiments.run_all --run-models
```

5つの期間（早期・早々期、前期、中期、後期、晩期）について MMCP モデルと NW モデルを推定し、結果を `.npy` ファイルとして保存する。

### MMCPモデルのみ実行（NWをスキップ）

```bash
uv run python -m bayesian_statistics.experiments.run_all --mmcp-only
```

NWモデルの推定をスキップし、MMCPモデルのみを実行する。NW比較図は生成されないが、他の図は正常に生成される。

### 図生成のみ

```bash
uv run python -m bayesian_statistics.experiments.run_all --generate-figures
```

保存済みのモデル結果から論文用の図を生成する。事前に `--run-models` の実行が必要。

### LOOCV評価のみ

```bash
uv run python -m bayesian_statistics.experiments.run_all --run-loocv
```

両モデルのLOOCV評価を実行し、Aitchison距離による比較表（Table 5.2）を生成する。

### テスト実行（短いMCMC）

```bash
uv run python -m bayesian_statistics.experiments.run_all --run-models --n-iter 100
```

### 静かに全実行

```bash
uv run python -m bayesian_statistics.experiments.run_all --all -q
```

## 出力ディレクトリ構造

```
bayesian_statistics/experiments/output/
├── results/                              # モデル結果
│   ├── mmcp_period_0/
│   │   ├── grid_probs.npy               # グリッド上の事後確率 (K, n_grid)
│   │   ├── site_probs.npy               # サイト上の事後確率 (K, n_sites)
│   │   ├── effect_distance.npy          # 距離効果
│   │   ├── effect_intercept.npy         # 切片効果
│   │   ├── effect_intercept_adjustment.npy  # 事前分布からの調整量
│   │   ├── effect_full.npy              # 完全モデル効果
│   │   └── lambda_star.npy              # λ* MCMCサンプル
│   ├── mmcp_period_1/
│   ├── ...
│   ├── nw_period_0/
│   │   ├── grid_probs.npy
│   │   └── site_probs.npy
│   └── ...
├── figures/                              # 生成された図
│   ├── fig_5_1_all_origins_periods.png  # 全産地×全期間 (4×5)
│   ├── fig_5_2_effect_distance.png      # 距離効果
│   ├── fig_5_2b_effect_intercept_adjustment.png  # 調整効果
│   ├── fig_5_3_estimated_vs_observed.png # 推定 vs 観測
│   └── fig_5_4_lambda_diagnostics.png   # MCMC診断
└── tables/
    └── table_5_2_loocv.csv              # LOOCV結果表
```

## 生成される図

| ファイル | 内容 |
|---------|------|
| `fig_5_1_all_origins_periods.png` | 4産地×5期間の事後平均確率マップ |
| `fig_5_2_effect_distance.png` | 距離効果の空間分布 |
| `fig_5_2b_effect_intercept_adjustment.png` | データによる事前分布からの調整量 |
| `fig_5_3_estimated_vs_observed.png` | 推定確率 vs 観測比率の散布図 |
| `fig_5_4_lambda_diagnostics.png` | λ* のトレースプロットと事後分布 |
| `fig_5_5_distance_prior.png` | 距離ベース事前分布 p₀(s) |
| `fig_5_6_model_comparison.png` | NW vs MMCP モデル比較 |
| `fig_5_7_uncertainty.png` | 推定の不確実性（事後標準偏差） |
| `fig_5_8_loocv_comparison.png` | LOOCV結果の棒グラフ比較 |

## モジュール構成

```
experiments/
├── __init__.py
├── run_all.py           # メインスクリプト
├── config.py            # 実験設定（ExperimentConfig）
├── output.py            # 出力管理（ProgressManager）
├── models/              # モデルランナー
│   ├── mmcp_runner.py   # MMCPモデル
│   └── nw_runner.py     # Nadaraya-Watsonモデル
├── evaluation/          # 評価
│   └── loocv.py         # LOOCV評価器
└── visualization/       # 可視化
    ├── map_template.py  # MapPlotterクラス
    ├── period_origin_grid.py  # 期間×産地グリッド
    ├── scatter_plots.py # 散布図
    ├── diagnostics.py   # MCMC診断
    └── style.py         # スタイル設定
```

## 設定のカスタマイズ

`config.py` の `ExperimentConfig` クラスで以下を変更可能:

- `n_iter`: MCMCイテレーション数
- `burn_in`: バーンイン期間
- `grid_subsample_ratio`: グリッドサブサンプリング比率
- `origins`: 産地名リスト
- `time_periods`: 期間名マッピング
- `source_weights`: 産地重み
