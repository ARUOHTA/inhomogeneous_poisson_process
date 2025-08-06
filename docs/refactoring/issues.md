# Model3リファクタリング問題・解決ログ

## 解決済み問題

### Phase 4: 統合デバッグフェーズ（2025-08-01）

#### I-001: KSBPModelのmodel_nameプロパティ不足
**現象**: 抽象メソッドmodel_nameが未実装でインスタンス化エラー
**原因**: BaseCompositionModelで必須のmodel_nameプロパティが未実装
**解決**: `@property def model_name(self) -> str: return "KSBP"`を追加
**ファイル**: `bayesian_statistics/models/composition/ksbp.py`

#### I-002: ModelComparisonでvariable_names引数不足
**現象**: `fit() missing 1 required positional argument: 'variable_names'`
**原因**: ModelComparisonでmodel.fit(preprocessor)のvariable_names引数不足
**解決**: ModelComparisonに__init__でvariable_names引数追加、fit()時に渡すよう修正
**ファイル**: `bayesian_statistics/models/evaluation/comparison.py`

#### I-003: インポートパス修正
**現象**: model_config.pyで古いインポートパスを使用
**解決**: 新しいディレクトリ構造に合わせてパスを更新
**ファイル**: `bayesian_statistics/models/config/model_config.py`

#### I-004: IPPモデル移動
**現象**: IPPモデルがintensity/ディレクトリに存在しない
**解決**: model3_ipp.pyをintensity/ディレクトリに移動
**ファイル**: `bayesian_statistics/models/intensity/ipp.py`

### 出力形式統一作業（2025-08-02）

#### I-005: predict_site_ratios()出力形式不統一
**現象**: AttributeError: 'DataFrame' object has no attribute 'items'
**原因**: 各モデルがDict[str, DataFrame]ではなくDataFrameを返していた
**解決**: 全モデルでDict[str, pl.DataFrame]形式に統一
**影響ファイル**:
- `nadaraya_watson.py`
- `bayesian_nw.py`

#### I-006: CompositionalMetricsメソッド不足
**現象**: AttributeError: 'CompositionalMetrics' object has no attribute 'bray_curtis_dissimilarity'
**原因**: 静的関数をインスタンスメソッドとして呼び出し
**解決**: `from .metrics import bray_curtis_dissimilarity`で直接インポート
**ファイル**: `bayesian_statistics/models/evaluation/unified_loocv.py`

#### I-007: Jensen-Shannon NaN問題
**現象**: Jensen-Shannon評価でNaN値が発生
**原因**: ゼロベクトルに対する処理が未実装
**解決**: ゼロベクトル判定と適切な距離値（0または log(2)）を返す処理追加
**ファイル**: `bayesian_statistics/models/evaluation/metrics.py`

## 発見された課題

### C-001: NadarayaWatsonモデルの予測精度
**現象**: 約80%の遺跡で全産地0予測
**原因**: データの希薄性と単一変数（average_elevation）のみ使用
**影響**: Jensen-Shannon、Bray-Curtisメトリクスで異常値
**推奨対策**:
- より多くの説明変数使用
- sigmaパラメータ調整（500→1000）
- zero_replacement値の調整（1e-6→1e-3）

## 注意事項

### 重要な制約
1. **新機能追加禁止**: 既存機能のみリファクタリング
2. **DRY原則厳守**: 重複コード排除最優先
3. **コード量削減**: 総行数を現在より減少
