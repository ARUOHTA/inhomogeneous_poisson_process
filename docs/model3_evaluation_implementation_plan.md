# Model 3 評価システム実装方針

## 1. 実装の大まかな方針

### 1.1 アーキテクチャ設計

```
実装の階層構造:
1. 基礎評価指標層: Compositional data評価指標
2. CV実行層: LOOCV実行システム
3. 統合評価層: 全時期・全産地評価システム
4. 可視化・レポート層: 結果表示とデバッグ
```

### 1.2 ファイル構成と役割

```
bayesian_statistics/
├── model3_compositional_metrics.py  # Compositional data評価指標
├── model3_loocv.py                  # LOOCV実行システム
├── model3_evaluation.py             # 全時期・全産地評価統合
└── model3_evaluation_report.py      # 結果可視化・レポート

scripts/
└── run_model3_evaluation.py         # CV実行用スクリプト

notebooks/
└── model3_evaluation_debug.ipynb    # 動作確認・デバッグ用
```

### 1.3 実装順序と依存関係

```
Phase 1: 基礎コンポーネント
├── model3_compositional_metrics.py
├── 基本的なLOOCV機能（model3_loocv.py）
└── 単体テスト・デバッグ

Phase 2: CV実行システム
├── 完全なLOOCV実装
├── NW推定量との統合
└── 小規模テスト実行

Phase 3: 統合・実用化
├── 全時期・全産地評価システム
├── 実行用スクリプト
└── 可視化・レポート機能
```

## 2. 技術的設計方針

### 2.1 データ構造

**入力データ形式**:
```python
# 各遺跡での産地構成比（compositional data）
composition_data = {
    site_id: {
        period: np.array([ratio_神津島, ratio_信州, ratio_箱根, ratio_高原山, ratio_その他])
        # 制約: sum(ratios) = 1.0
    }
}
```

**評価結果形式**:
```python
evaluation_results = {
    'period': int,
    'trial_results': [
        {
            'trial_id': int,
            'test_site_id': int,
            'observed_composition': np.array,
            'predicted_composition': np.array,
            'aitchison_distance': float,
            'total_variation': float
        }
    ],
    'summary_statistics': {
        'mean_aitchison_distance': float,
        'std_aitchison_distance': float,
        'mean_total_variation': float,
        'std_total_variation': float
    }
}
```

### 2.2 エラーハンドリング戦略

1. **Zero値の処理**: log-ratio計算時の0値に対する小さな正の値の追加
2. **Simplex制約の確認**: 構成比の合計が1になることの検証
3. **NaN/Inf値の検出**: 計算結果の妥当性チェック
4. **データ不足の処理**: 特定時期・遺跡でのデータ欠損への対応

### 2.3 計算効率化

1. **並列処理**: 複数のLOOCV試行の並列実行
2. **キャッシュ機能**: 重み行列の再利用
3. **バッチ処理**: 複数遺跡の同時予測
4. **メモリ管理**: 大規模データセットでのメモリ使用量最適化

## 3. 実装詳細仕様

### 3.1 Compositional Data評価指標

**Aitchison Distance**:
```python
def aitchison_distance(x: np.ndarray, y: np.ndarray, zero_replacement: float = 1e-6) -> float:
    """
    Aitchison距離の計算
    
    Parameters:
    - x, y: 産地構成比ベクトル（sum = 1）
    - zero_replacement: 0値の置換値
    
    Returns:
    - aitchison_distance: float
    """
```

**Total Variation**:
```python
def total_variation(x: np.ndarray, y: np.ndarray, zero_replacement: float = 1e-6) -> float:
    """
    Total Variationの計算
    
    Formula: 1/(2D) * sum_i sum_j (ln(x_i/x_j) - ln(y_i/y_j))^2
    """
```

### 3.2 LOOCV実行システム

**主要機能**:
```python
class LOOCVEvaluator:
    def __init__(self, preprocessor, config):
        self.preprocessor = preprocessor
        self.config = config
    
    def run_single_trial(self, period: int, test_site_id: int) -> dict:
        """単一のLOOCV試行を実行"""
    
    def run_period_evaluation(self, period: int, n_trials: int = 100) -> dict:
        """特定時期の全体評価を実行"""
    
    def run_all_periods_evaluation(self, n_trials: int = 100) -> dict:
        """全時期の評価を実行"""
```

### 3.3 NW推定量との統合

**再学習機能**:
```python
def retrain_nw_estimator(preprocessor, excluded_site_ids: List[int], config) -> NadarayaWatsonEstimator:
    """
    指定された遺跡を除外してNW推定量を再学習
    
    Steps:
    1. 除外遺跡のデータをマスク
    2. 重み行列の再計算
    3. NW推定量の再フィット
    """
```

## 4. テスト・デバッグ戦略

### 4.1 単体テスト

1. **評価指標の妥当性**: 既知の結果との比較
2. **データ前処理**: 構成比制約の確認
3. **数値安定性**: 極端なケースでの動作確認

### 4.2 統合テスト

1. **小規模データセット**: 2-3遺跡での動作確認
2. **計算時間測定**: 実用的な実行時間の確認
3. **結果の再現性**: ランダムシードでの一貫性確認

### 4.3 可視化デバッグ

1. **中間結果の表示**: 重み行列、予測値の可視化
2. **評価指標の分布**: ヒストグラム、箱ひげ図
3. **時期別・遺跡別の性能**: 詳細な比較分析

## 5. 実行・運用方針

### 5.1 実行環境

- **メモリ要件**: 16GB以上推奨（距離行列のサイズ考慮）
- **CPU**: マルチコア推奨（並列処理のため）
- **実行時間**: 100試行で数時間程度を想定

### 5.2 設定管理

```python
@dataclass
class EvaluationConfig:
    n_trials: int = 100
    zero_replacement: float = 1e-6
    random_seed: int = 42
    parallel_jobs: int = 4
    output_dir: str = "output/evaluation"
    save_intermediate: bool = True
```

### 5.3 結果出力

1. **CSV形式**: 詳細な試行結果
2. **JSON形式**: 統計サマリー
3. **図表**: 可視化結果（PNG/PDF）
4. **ログ**: 実行過程の記録

## 6. リスク管理

### 6.1 計算上のリスク

- **数値的不安定性**: Zero値、極端な比率への対処
- **計算時間**: 大規模データでの実行時間爆発
- **メモリ不足**: 距離行列の保持限界

### 6.2 データ品質リスク

- **データ欠損**: 特定時期・遺跡での不完全データ
- **異常値**: 予期しない構成比パターン
- **バージョン不整合**: 前処理データとの非互換

### 6.3 結果解釈リスク

- **統計的妥当性**: サンプルサイズ不足
- **考古学的妥当性**: 専門知識との齟齬
- **過度の一般化**: 限定的データからの過度な結論

---

この方針に基づいて、段階的に実装を進めます。まず基礎コンポーネントから開始し、動作確認を行いながら統合システムを構築します。