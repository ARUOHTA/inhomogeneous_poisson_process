# Model 3rd 評価指標とクロスバリデーション戦略

## 概要
model3では2つの主要な推定問題があります：
1. **産地構成比推定（Nadaraya-Watson推定量）**: 各時期・産地における空間的な構成比分布の推定
2. **遺跡存在確率推定（IPP）**: 非斉次ポアソン過程による遺跡の存在確率推定

これらのモデルの性能評価と最適なハイパーパラメータ選択のための評価指標とクロスバリデーション手法を検討します。

## 1. 産地構成比推定（Nadaraya-Watson推定量）の評価

### 1.1 ハイパーパラメータ
- `sigma`: グリッド間のカーネルバンド幅（空間的スムージング）
- `sigma_for_sites`: 遺跡間のカーネルバンド幅（遺跡データのスムージング）

### 1.2 評価指標（Compositional Data Analysis）

産地構成比データは各遺跡での産地別出土比率の合計が1となる**compositional data**（組成データ）です。このため、simplex空間上のデータとして適切な評価指標を使用する必要があります。

#### 1.2.1 採用する主要評価指標（2つ）

1. **Aitchison Distance**（主要指標1）
   - Compositional dataの標準的距離測度
   - CLR変換後のユークリッド距離として定義
   - $d_{A}(x,y) = \sqrt{\sum_{i=1}^{D} \left(\ln\frac{x_i}{g(x)} - \ln\frac{y_i}{g(y)}\right)^2}$
   - ここで$g(x) = \sqrt[D]{\prod_{i=1}^{D} x_i}$は幾何平均
   - **主要性質**: Scale invariance, Perturbation invariance, Permutation invariance, Subcompositional coherence
   - **参考文献**: Aitchison, J. (1982). The statistical analysis of compositional data. Journal of the Royal Statistical Society Series B, 44(2), 139-177.

2. **Total Variation**（主要指標2）
   - $\text{totvar}(x,y) = \frac{1}{2D}\sum_{i=1}^{D}\sum_{j=1}^{D}(\ln(x_i/x_j) - \ln(y_i/y_j))^2$
   - 組成全体の変動性を定量化
   - 各成分ペア間のlog-ratio差の総合評価
   - **利点**: 全成分間の相対関係を同等に評価、解釈が直感的
   - **参考文献**: Aitchison, J. (1986). The Statistical Analysis of Compositional Data. Chapman & Hall.

#### 1.2.2 その他のCompositional Data指標（参考）

**Log-Ratio Transform Based Metrics**

1. **CLR (Centered Log-Ratio) Transform**
   - $\text{clr}(x) = \left(\ln\frac{x_1}{g(x)}, \ln\frac{x_2}{g(x)}, ..., \ln\frac{x_D}{g(x)}\right)$
   - CLR変換後のユークリッド距離がAitchison距離と等価
   - **利点**: 元の成分数を保持、全体に対する各部分の解釈が直感的

2. **ILR (Isometric Log-Ratio) Transform**
   - Simplex空間から実数空間への等距写像
   - 直交座標系を用いた変換により次元削減（D-1次元）
   - **利点**: 統計解析（回帰、PCA等）に適用可能
   - **参考文献**: Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G., Barceló-Vidal, C. (2003). Isometric logratio transformations for compositional data analysis. Mathematical Geology, 35(3), 279-300.

3. **ALR (Additive Log-Ratio) Transform**
   - $\text{alr}(x) = \left(\ln\frac{x_1}{x_D}, \ln\frac{x_2}{x_D}, ..., \ln\frac{x_{D-1}}{x_D}\right)$
   - 参照成分の選択が恣意的

4. **Variation Array**
   - $V_{ij} = \text{var}(\ln(x_i/x_j))$
   - 成分間の相対的変動を評価

5. **Subcompositional Coherence**
   - 部分組成での結果が全組成での結果と一致するかを確認
   - Compositional dataの基本原理の一つ

#### 1.2.3 二次的評価指標（従来手法との比較用）

1. **Bray-Curtis Dissimilarity**
   - $BC(x,y) = \frac{\sum_{i=1}^{D}|x_i - y_i|}{\sum_{i=1}^{D}(x_i + y_i)}$
   - 従来手法との比較のために併用

2. **Jensen-Shannon Divergence**
   - 情報理論的距離測度
   - Aitchison距離との性能比較に使用

**主要参考文献**:
- Aitchison, J. (1982). The statistical analysis of compositional data. Journal of the Royal Statistical Society Series B, 44(2), 139-177.
- Aitchison, J. (1986). The Statistical Analysis of Compositional Data. Chapman & Hall.
- Pawlowsky-Glahn, V., Egozcue, J.J., Tolosana-Delgado, R. (2015). Modeling and Analysis of Compositional Data. John Wiley & Sons.
- Pawlowsky-Glahn, V., Egozcue, J.J. (2001). Geometric approach to statistical analysis on the simplex. Stochastic Environmental Research and Risk Assessment, 15(5), 384-398.

### 1.3 クロスバリデーション手法

#### 1.3.1 Leave-One-Out Cross Validation (LOOCV)

**目的**: 個別遺跡での予測性能の客観的評価

**手法**:
1. 全遺跡から1つ（または少数）の遺跡をランダムに選択して除外
2. 残りの遺跡データでNadaraya-Watson推定量を学習
3. 除外した遺跡での産地構成比を予測
4. 実際の構成比と予測値をcompositional data評価指標で比較
5. この過程を異なる遺跡に対して繰り返し実行
6. 全ての試行結果の平均として最終的な評価指標を算出

**実装戦略**:
```python
def loocv_evaluation(sites_data, n_trials=100):
    """
    Leave-One-Out Cross Validationによる評価

    Parameters
    ----------
    sites_data : データ構造
        全遺跡のデータ
    n_trials : int
        試行回数（デフォルト100回）

    Returns
    -------
    evaluation_results : dict
        各評価指標の平均値と標準偏差
    """
    results = []

    for trial in range(n_trials):
        # ランダムに1つ（または数個）の遺跡を選択
        test_site_ids = random.sample(site_ids, k=1)  # または k=2-3
        train_site_ids = [id for id in site_ids if id not in test_site_ids]

        # 訓練データでNW推定量を学習
        nw_estimator = fit_nadaraya_watson(train_site_ids)

        # テスト遺跡で予測
        predicted_ratios = nw_estimator.predict(test_site_ids)
        observed_ratios = get_observed_ratios(test_site_ids)

        # Compositional data評価指標で評価
        aitchison_dist = calculate_aitchison_distance(predicted_ratios, observed_ratios)
        total_variation = calculate_total_variation(predicted_ratios, observed_ratios)

        results.append({
            'aitchison_distance': aitchison_dist,
            'total_variation': total_variation,
            'trial_id': trial,
            'test_sites': test_site_ids
        })

    return aggregate_results(results)
```

**評価対象**:
- 全時期: target_period = 0, 1, 2, 3, 4（早期・早々期、前期、中期、後期、晩期）
- 全産地: ["神津島", "信州", "箱根", "高原山", "その他"]の産地構成比ベクトル
- **注意**: 産地別個別評価ではなく、産地構成比ベクトル全体として評価

**メリット**:
- 各遺跡での個別予測性能を直接評価
- サンプルサイズが小さくても実行可能
- 計算コストが比較的低い
- 結果の解釈が直感的

### 1.4 評価実施方針

**LOOCV実行時の構成比ベクトル**:
各遺跡において、全産地の構成比を1つのベクトルとして扱います：
- $x = (x_{神津島}, x_{信州}, x_{箱根}, x_{高原山}, x_{その他})$
- ここで $\sum_{i} x_i = 1$（simplex制約）

**時期別評価**:
全時期について個別にLOOCVを実行し、時期ごとの予測性能を比較します。

## 2. 遺跡存在確率推定（IPP）の評価

### 2.1 ハイパーパラメータ
- MCMCパラメータ: `num_iterations`, `burn_in`
- 事前分布パラメータ: `prior_beta_mean`, `prior_beta_cov`, `prior_lambda_shape`, `prior_lambda_rate`

### 2.2 評価指標

#### 2.2.1 分類性能指標
1. **AUC-ROC**
   - 遺跡存在の確率予測の識別能力
   - 閾値に依存しない評価

2. **AUC-PR（Precision-Recall）**
   - 不均衡データ（遺跡は稀）に適した指標
   - 実用的な検出性能を評価

3. **キャリブレーション**
   - 予測確率と実際の頻度の一致度
   - Reliability diagram, Brier scoreで評価

#### 2.2.2 空間パターン評価
1. **Hotspot Detection Accuracy**
   - 高確率領域での遺跡検出率
   - 考古学的に意味のある領域の特定精度

2. **Cross-K Function**
   - 予測分布と実際の遺跡分布の空間的類似度

### 2.3 IPPモデルの評価戦略

IPPモデルについては、産地構成比推定とは異なる評価アプローチを採用します：

1. **全データ使用**: 遺跡存在確率の推定には全遺跡データを使用
2. **事後予測チェック**: MCMCサンプルを用いた予測分布の妥当性検証
3. **空間パターン比較**: 予測された確率分布と実際の遺跡分布の空間的類似度評価
4. **MCMC診断**: 収束性とサンプリング効率の確認

## 3. 統合評価戦略

### 3.1 評価の優先順位
1. **産地構成比推定（NW推定量）の評価**: compositional data評価指標を用いたLOOCVによる性能評価
2. **IPPモデルの妥当性検証**: MCMC診断と空間パターンの評価
3. **考古学的解釈可能性**: 専門知識との整合性確認

### 3.2 ドメイン知識との整合性
- 考古学的に知られているパターンとの一致度
- 専門家による定性的評価
- 既存研究との比較

### 3.3 計算効率性
- LOOCV試行回数と精度のトレードオフ
- メモリ使用量の最適化
- 全時期・全産地評価の効率化

## 4. 実装計画

### 4.1 評価パイプライン
1. **Compositional Data評価モジュール**: Aitchison距離等の専用指標実装
2. **LOOCV実行モジュール**: 効率的なクロスバリデーション
3. **可視化モジュール**: 評価結果の比較表示
4. **レポート生成モジュール**: 自動評価レポート作成

### 4.2 段階的実装
**Phase 1**: Compositional Data評価指標の実装
- Aitchison距離の計算
- Total Variationの算出
- その他参考指標（CLR, ILR, ALR変換等）

**Phase 2**: LOOCV評価システムの構築
- ランダムサンプリングによる遺跡除外
- NW推定量の再学習システム
- 全時期・全産地での評価指標の集計と統計処理

**Phase 3**: 性能評価と考古学的検証
- 時期別・産地構成比パターンの性能比較
- 結果の妥当性確認
- 考古学的解釈との整合性検証

### 4.3 ファイル構成
```
bayesian_statistics/
├── model3_compositional_metrics.py  # Compositional data評価指標
├── model3_loocv.py                  # LOOCV実装
├── model3_evaluation.py             # 全時期・全産地評価システム
└── model3_evaluation_report.py      # 評価結果のレポート生成
```

## 5. 期待される成果

1. **Compositional Data Analysis の適用**: 産地構成比データに適した評価指標の確立
2. **予測性能の定量化**: Aitchison距離・Total Variationによる予測精度の客観的評価
3. **時期別性能比較**: 全時期における産地構成比予測性能の比較分析
4. **考古学的洞察**: 各時期における産地構成比分布パターンの定量的理解

## 6. 注意事項

### 6.1 Compositional Dataの特性
- **Simplex制約**: 各遺跡での産地構成比の合計は1
- **相対情報**: 絶対量ではなく相対的比率が重要
- **Zero問題**: 特定産地が0の場合のlog-ratio変換への対処
- **Subcompositional coherence**: 部分組成の解析結果の一貫性

### 6.2 LOOCVの前提と限界
- **独立性の仮定**: 遺跡間の空間的相関の無視
- **サンプルサイズ**: 小規模データでの評価精度の限界
- **計算コスト**: 試行回数と計算時間のバランス
- **代表性**: 除外遺跡の選択バイアス

### 6.3 考古学的解釈の注意点
- **発掘バイアス**: 遺跡発見の地理的・時代的偏り
- **産地同定精度**: 化学分析の誤差と分類の不確実性
- **時代区分**: 早期・早々期の定義と遺物の帰属問題
- **流通網の複雑性**: 直接交換と間接交換の区別困難

### 6.4 統計的妥当性の確保
- **Multiple testing**: 複数ハイパーパラメータの同時検定
- **過適合の検出**: 訓練性能と汎化性能の乖離
- **結果の再現性**: ランダムサンプリングの影響評価

---

この戦略に基づいて、段階的に評価システムを構築し、model3の性能を包括的に評価していきます。
