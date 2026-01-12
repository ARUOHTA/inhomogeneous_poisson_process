# Barbet-Massin et al. (2012) - 擬似欠席データの選択ガイドライン

## 基本情報
- **タイトル**: Selecting pseudo-absences for species distribution models: how, where and how many?
- **著者**: Morgane Barbet-Massin, Frédéric Jiguet, Cécile Hélène Albert, Wilfried Thuiller
- **ジャーナル**: Methods in Ecology and Evolution
- **年**: 2012
- **DOI**: 10.1111/j.2041-210X.2011.00172.x
- **関連論点**: 論点3（presence-onlyデータの扱い）

## 主な貢献
1. **包括的なシミュレーション研究**: 仮想種を用いて、擬似欠席データ（pseudo-absences, PA）の選択方法（how, where, how many）が種分布モデル（SDM）の予測精度に与える影響を系統的に評価
2. **モデル手法別の最適戦略**: 7つのSDM（回帰、分類、機械学習）それぞれに最適なPA選択戦略を特定
3. **バイアスの影響評価**: 気候的・空間的にbiasedなpresenceデータでの最適戦略を明確化
4. **実践的ガイドライン**: 実務で直接使えるテーブル形式の推奨事項を提供

## 手法の概要

### 仮想種の作成
**2つの仮想種を欧州スケールで作成**:
- 環境変数: PCA第1・第2軸（温度・降水量から導出）
- ニッチ形状: 両変数に対してbell-shaped（ガウス型）応答曲線
- 確率閾値: 0.25で二値化した潜在分布
- 実現分布: 各ピクセルで確率に応じた二項分布から約2700個のpresenceを生成

**利点**: 真の分布を既知として、サンプリング設計やバイアスの影響を定量評価可能

### 3種類のサンプリングバイアス
1. **気候的バイアス (Climatically biased)**: ニッチの一部のみをサンプリング
   - ガウス応答曲線の平均をずらした確率面から約1000 presenceを抽出
   - 結果: fundamental climatic nicheの全範囲をカバーしない

2. **空間的バイアス1 (Countries bias)**: 分布の一部国を除外
   - 約1000 presenceを特定地域からのみ抽出

3. **空間的バイアス2 (Transportation bias)**: 道路・鉄道沿いのみサンプリング
   - 約1000 presenceを交通路付近からのみ抽出

### 4種類のPseudo-absence選択方法
1. **Random**: 全域からランダム選択（presenceを除く）
2. **SRE (Surface Range Envelope)**: presence-onlyモデルで推定された不適域からランダム選択
   - 環境的除外（climatic exclusion）
3. **1° far**: 全presenceから緯度経度1度以上離れた地点からランダム選択
4. **2° far**: 全presenceから緯度経度2度以上離れた地点からランダム選択
   - 地理的除外（geographical exclusion）

### 7つのSDM手法
**回帰手法**:
- GLM (Generalized Linear Model)
- GAM (Generalized Additive Model)
- MARS (Multiple Adaptive Regression Splines)

**分類手法**:
- MDA (Mixture Discriminant Analysis)
- CTA (Classification Tree Analysis)

**機械学習手法**:
- BRT (Boosted Regression Trees)
- RF (Random Forest)

### 実験デザイン
**検証する6つの質問**:
(a) どのprevalence（presences/absencesの比）が最も高精度か？
(b) PA選択の最適な反復回数は？
(c) 反復ごとの最適なPA数と重み付けスキームは？
(d) どのPA生成方法が最も高精度か？
(e) サンプリングバイアスは最適なPA使用にどう影響するか？
(f) どのパラメータ（PA数、選択方法、重み付け）が予測精度に最も影響するか？

**Presence数**: 30, 100, 300, 1000
**PA数**: 100, 300, 1000, 3000, 10000
**重み付け**: equal weight (presenceとabsenceの総重みを等しくする) vs. un-weighted
**反復**: 20回の異なるpresence選択 × 各条件で最大20回のPA選択

### 評価指標
- **AUC**: ROC曲線下面積
- **TSS (True Skill Statistics)**: sensitivity + specificity - 1
- **Sensitivity**: presenceを正しく予測した割合
- **Specificity**: absenceを正しく予測した割合
- **閾値**: TSSを最大化する確率閾値で二値化

## 実証結果

### (A) Prevalenceの影響
**3つのグループに分類**（図3）:

1. **GAM**: prevalenceの影響を受けない
2. **MARS, MDA**: prevalence（presencesの比率）が高いほど精度向上
3. **GLM, BRT, RF, CTA**: presenceが absenceの1/10に達するまで精度向上、その後plateau
   - CTA: presenceとabsenceが同数で最高精度

### (B) 反復回数の最適化
**PA数に依存** （図4）:
- **10000 PA**: 反復不要（1回で十分）
- **1000 PA**: GAM・CTAで5回、他は反復不要
- **100 PA**:
  - BRT: 4回
  - GLM, MARS, MDA, CTA, RF: 7回
  - GAM: 12回
  - ただし、100 PAでは変動が大きいため20回が安全

### (C) PA数と重み付けの最適化
**3つのグループに分類** （図5, 6）:

1. **GLM, GAM**:
   - PA数の影響小
   - Equal weightで精度向上
   - **推奨**: 10000 PA with equal weight

2. **CTA, BRT, RF**:
   - presenceと同数のPAで最高精度
   - PA数が異なる場合、equal weightで精度向上
   - **推奨**: same number as presences with equal weight

3. **MARS, MDA**:
   - 100 PA × 複数回の反復で最高精度
   - Equal weightで精度向上
   - **推奨**: 100 PA × 10 runs with equal weight

### (D) PA選択方法の影響
**回帰手法 (GLM, GAM, MARS)**:
- **Random**が最も高精度（一貫して）
- 例外: 気候的にbiasedなpresenceの場合、**2° far**が最適

**分類・機械学習手法 (MDA, BRT, CTA, RF)**:
- 選択方法の影響は小さいが、明確な傾向あり：
  - Few presences: **2° far**が最適
  - Many presences: **SRE**が最適
- Random: より高いspecificity
- 2° far/1° far: より高いsensitivity

### (E) バイアスの影響
**空間的バイアス**: Unbiasedと同様の結果
- 例外: MDAでは**Random**が最適

**気候的バイアス**:
- GLM, GAM, MARS: **Random**が不適
  - Few presences: **SRE**が最適
  - Many presences: **2° far**が最適
- MDA, CTA, BRT, RF: **2° far**が最適

**理由**: 気候的バイアスではfundamental nicheの全範囲がカバーされないため、環境的除外（SRE）や大きな地理的除外（2° far）がfalse absencesのリスクを低減

### (F) 最も影響の大きいパラメータ
**GLM, GAM**（図7）:
- 30 presences: presence選択のランダム性が最大要因
- 100+ presences: **PA選択方法**が最大要因

**MARS, MDA, CTA, BRT, RF**:
- 全presences数で**PA数**が最大要因
- PA選択方法の影響は副次的（presences数とともに増加）

**SDM手法間の変動** （図8）:
- 100+ presencesでは、PA選択の影響 < SDM手法間の変動
- モデル選択自体が最も重要な意思決定

### 評価指標間の相関
- AUC vs TSS: $r = 0.82 \pm 0.10$（高相関）
- 評価指標の選択は相対的な性能比較に影響しない

## 限界・残された課題

### 著者が述べる限界
1. **研究範囲の影響**: PA数の最適値は研究範囲の空間extent に依存
   - より広い範囲 → より多くのPAが必要（環境多様性が増すため）
   - 情報的なPAを十分選択するために、範囲に応じた調整が必要

2. **仮想種の単純化**: ガウス型応答曲線のみを検証
   - 他の応答曲線（非対称、複雑な形状）での検証が必要
   - ただし、相対的な性能比較には影響しないと著者は主張

3. **Ensemble forecastの課題**: 異なるSDMで最適なPA戦略が異なる
   - 全SDMで同じPAを使うと不公平な比較になる
   - 提案: SDMをグループ化（GLM+GAM, BRT+RF）して各グループで最適PA戦略を適用

4. **Sensitivity vs Specificity のトレードオフ**:
   - 閾値選択（TSS最大化）がsensitivity寄り
   - 別の閾値（sensitivity-specificity差最小化）ではspecificity寄り

### 追加の課題
1. **実データでの検証**: 仮想種では検証したが、複雑な実データでの性能は未確認
2. **時空間動態**: 分布の時間変化や動的過程は考慮されていない
3. **種間相互作用**: 生物学的相互作用（競争、捕食）は無視
4. **検出確率**: 不完全な検出（detection probability < 1）の影響は未評価
5. **計算コスト**: 多数の反復と大量のPAは計算負荷が高い

## 本研究(MMCP)との関係

### 直接的な関係
1. **Presence-onlyデータの扱い**: 黒曜石遺跡データはpresence-only
   - 遺跡が"ない"地点 ≠ 真のabsence（未発見・未調査の可能性）
   - Pseudo-absence戦略がモデル精度に直結

2. **IPPモデルへの応用**: MMCPの点過程部分は強度関数 $\lambda(s)$ をモデル化
   - Renner et al. (2015) が示したIPP-Maxent等価性から、PA選択戦略が重要
   - Background pointsの選択方法がIPPの強度推定に影響

3. **モデル選択の指針**: GLM系（logistic regression）を使うか、機械学習系を使うかでPA戦略が異なる
   - MMCPでGLMベースなら: 10000 PA with random selection
   - BRT/RFベースなら: presenceと同数のPA

4. **バイアス対応**: 遺跡データは空間的にbiasedの可能性（発掘調査の偏り）
   - 交通路沿い、特定地域への集中など
   - 空間的バイアス下では**Random**が頑健

### 本研究への示唆
- **IPPの2段階構造への適用**:
  - Stage 1（点位置モデリング）: presence-onlyから強度関数推定
  - PA選択がIPP強度 $\lambda(s)$ の推定精度に影響
  - ベースラインIPPモデルの構築で本論文のガイドライン適用

- **組成データへの拡張**:
  - 論文は二値（presence/absence）だが、MMCPはマーク付き点過程
  - 点位置のPA選択 → 組成データへの間接的影響
  - マーク（組成）モデルの条件付き分布推定でも、背景データの選択が重要

- **モデル診断**: TSS、AUC、Sensitivity、Specificityを複数評価
  - MMCPでも、点位置予測の評価でこれらの指標を利用可能
  - 考古学的には高sensitivityが望ましい（未発見遺跡の予測）

### 借用すべき手法
1. **仮想データシミュレーション**:
   - モデル検証のために、既知分布の仮想遺跡データを生成
   - 異なるPA戦略の性能比較

2. **反復平均（Ensemble）**:
   - 複数のPA選択を反復し、予測分布を平均化
   - 特に100〜1000 PAの場合、10回程度の反復が推奨

3. **Equal weighting scheme**:
   - Presenceとabsenceの総重みを等しくする重み付け
   - 多くのSDMで精度向上

4. **BIOMODパッケージ**: Rパッケージでの実装
   - `biomod` パッケージのPA選択関数を参考に

### 拡張すべき点
1. **IPP特有のPA**: 点過程強度のquadratureスキーム
   - Bermanの変換: $y_i \in \{0,1\}$ + integration weights
   - 本論文のPAより理論的に精緻なbackground points

2. **空間的自己相関の考慮**:
   - PAを空間的に均等配置（空間層別サンプリング）
   - 環境空間と地理空間の両方での層別化

3. **検出確率の明示的モデリング**:
   - Occupancy modelとの統合（Dorazio 2014）
   - Survey effortを共変量として組み込み

4. **マーク付きPA**:
   - 点位置だけでなく、組成データのabsenceも定義必要
   - 「その産地の黒曜石がない」≠ 「遺跡がない」

## 引用すべき箇所

### Presence-absenceモデルの優位性（Introduction, p.2）
> "Comparisons of various SDM show that presence-absence models tend to perform better than presence-only models (Elith et al. 2006). Thus, presence-absence models are increasingly used when only presence data is available, by creating artificial absence data (usually called pseudo-absences or background data)."

**用途**: Presence-onlyデータからpresence-absence モデルを構築する動機を説明

### 仮想実験の利点（Introduction, p.3）
> "Generalisation and application of the conclusions of these empirical studies are therefore of limited interest in general compared with conclusions from virtual experiments where results or patterns can be compared with the known truth."

**用途**: シミュレーション研究の重要性を正当化

### バイアスの種類（Introduction, p.3）
> "Geographically biased presence data could arise from sampling along main roads or railways, or within a subset of the countries where the species occurs. Climatically biased presence data can result either from a spatially biased sampling design, or from sampling that was not carried out over the whole environmental range of a given species."

**用途**: 考古学データのバイアス（交通路沿い、特定地域、標高範囲）を議論する際に引用

### 回帰手法でのPA数（Results, p.5）
> "The models could be separated into three groups according to the effect of prevalence on their predictive accuracy. GAM behaved differently from the others given this technique was not influenced by prevalence. The accuracy of MARS and MDA increased with prevalence, whereas the accuracy increased until an asymptote when the number of presences reached one tenth of the number of absences for GLM, BRT and RF or reached the same amount as the number of absences for CTA."

**用途**: GLM系モデルでのPA数の設定根拠を示す

### Random selectionの優位性（Results, p.7）
> "For GLM, GAM and MARS, randomly selected pseudo-absences produced the most accurate models."

**用途**: IPPモデル（GLMベース）でのrandom PA選択を正当化

### False absenceへの感度（Discussion, p.9）
> "Such differences regarding the best method of generating pseudo-absences indicate that regression techniques were less sensitive to false absences than classification and machine-learning techniques."

**用途**: GLM系がfalse absencesに頑健であることを強調

### 空間extentの影響（Discussion, p.9）
> "The optimal number of pseudo-absences to generate in each run is therefore likely to depend on the spatial extent of the study, which influences environmental variability. At a given spatial resolution, a higher number of pseudo-absences may be needed to optimise model performance for a larger spatial extent of the study, to ensure the selection of enough informative pseudo-absences."

**用途**: 研究範囲に応じたPA数の調整を議論する際に引用

### 実践的推奨事項（Conclusion, Table 1）
> "Overall, we recommend the use of a large number (e.g. 10000) of pseudo-absences with equal weighting for presences and absences when using GLM and GAM, averaging several runs with relatively fewer pseudo-absences (e.g. 100) with equal weighting for presences and absences with MARS and MDA, and using the same amount of pseudo-absences as the amount of available presences (averaging several runs if few pseudo-absences) for CTA, BRT and RF."

**用途**: MMCPのIPP部分で採用するPA戦略の具体的根拠として引用

### Ensemble forecastの課題（Discussion, p.9）
> "However, we have shown here that the optimal way of creating and using pseudo-absences information differs widely across SDM. The best way of using pseudo-absences through an ensemble forecast technique could therefore be to use pseudo-absences differently for each SDM."

**用途**: 複数モデルのアンサンブル時の注意点を議論

### Sensitivity最大化（Discussion, p.10）
> "In such studies, the 'SRE', '1 and 2° far' methods can be used as well as other methods for selecting pseudo-absences outside both spatially and climatically suitable areas. The selection of fewer pseudo-absences in each replicate also yielded better sensitivity."

**用途**: 考古学的応用（未発見遺跡の予測）でhigh sensitivityを重視する際の戦略
