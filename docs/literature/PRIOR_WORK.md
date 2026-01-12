# 先行研究の整理

修士論文 第3章（関連研究）執筆のための先行研究整理。

## 論点の再確認

| 論点 | 内容 | キーワード |
|------|------|-----------|
| **論点1** | 空間的依存性のモデリング | GP, NNGP, Vecchia, 空間変化係数 |
| **論点2** | 組成データのモデリング | Aitchison, multinomial logit, Pólya-Gamma, Dirichlet |
| **論点3** | Presence-onlyデータのモデリング | 点過程, Cox過程, pseudo-absence, thinning |
| **論点4** | なぜ統合されてこなかったか | 計算困難性, 識別可能性 |
| **論点5** | 本研究の位置づけ | MMCP, 統合モデル |

---

## 優先度リスト（未変換）

### 最優先 ★★★

| 論文 | 論点 | 理由 |
|------|------|------|
| **Eckardt et al. 2025** | 論点2・3・4 | Composition-valued marks（MMCPに最も近い）|
| **Gelfand & Schliep 2018** | 論点3 | 点パターンのベイズ推論（教科書的） |
| **Gelfand & Shirota 2019** | 論点3・4 | Presence-only + preferential sampling |
| **Gonçalves & Gamerman 2018** | 論点1・3 | Cox過程のベイズ推論 |

### 高優先 ★★

| 論文 | 論点 | 理由 |
|------|------|------|
| Fithian & Hastie 2013 | 論点3 | Presence-only理論的基礎 |
| Dorazio 2014 | 論点3 | Imperfect detection + presence-only |
| Pacifici et al. 2017 | 論点4 | データ統合フレームワーク |
| Zens et al. 2023 | 論点2 | Ultimate PG Samplers（計算効率） |
| Wehrhahn et al. 2022 | 論点2 | Compositional data + Bayesian nonparametric |
| Adams et al. 2009 | 論点1・3 | GP intensity Poisson process |

### 中優先 ★

| 論文 | 論点 | 理由 |
|------|------|------|
| Gelfand et al. 2003 | 論点1 | 空間変化係数モデル |
| Finley & Banerjee 2020 | 論点1 | spBayes実装 |
| Rue et al. 2009 | 論点1 | INLA（比較手法） |
| Barbet-Massin et al. 2012 | 論点3 | Pseudo-absence選択 |
| Phillips et al. 2006 | 論点3 | MaxEnt（比較手法） |
| Cecconi et al. 2016 | 論点3 | Preferential sampling |
| Narayanan et al. 2021 | 論点3 | Marked spatio-temporal point process |
| da Silva & Gamerman 2022 | 論点3 | Preferential sampling |

### 参考 ☆

| 論文 | 論点 | 理由 |
|------|------|------|
| Blei & Jordan 2006 | 論点2 | Dirichlet process（背景知識） |
| Walker 2007 | 論点2 | Slice sampling（背景知識） |
| Quintana et al. 2020 | 論点2 | Dependent Dirichlet process |

---

## 変換済み論文（35件）

### 最優先 ★★★（完了）

| ファイル | 論文 | 論点 |
|---------|------|------|
| `moreira2022.md` | Moreira & Gamerman (2022) | 論点3 |
| `moreira2024.md` | Moreira et al. (2024) | 論点3・4 |
| `polson2012.md` | Polson et al. (2013) | 論点2 |
| `linderman2015.md` | Linderman et al. (2015) | 論点2 |
| `renner2015.md` | Renner et al. (2015) | 論点3 |
| `datta2016.md` | Datta et al. (2016) | 論点1 |
| `eckardt2025.md` | Eckardt et al. (2025) | 論点2・3・4 |
| `gelfand_schliep2018_main.md` | Gelfand & Schliep (2018) | 論点3 |
| `gelfand_shirota2019.md` | Gelfand & Shirota (2019) | 論点3・4 |
| `goncalves_gamerman2018.md` | Gonçalves & Gamerman (2018) | 論点1・3 |

### 高優先 ★★（完了）

| ファイル | 論文 | 論点 |
|---------|------|------|
| `fithian2013.md` | Fithian & Hastie (2013) | 論点3 |
| `dorazio2014.md` | Dorazio (2014) | 論点3 |
| `pacifici2017.md` | Pacifici et al. (2017) | 論点4 |
| `zens2023.md` | Zens et al. (2023) | 論点2 |
| `wehrhahn2022.md` | Wehrhahn et al. (2022) | 論点2 |
| `adams2009.md` | Adams et al. (2009) | 論点1・3 |

### 中優先 ★（完了）

| ファイル | 論文 | 論点 |
|---------|------|------|
| `gelfand2003.md` | Gelfand et al. (2003) | 論点1 |
| `finley2020.md` | Finley & Banerjee (2020) | 論点1 |
| `rue2009.md` | Rue et al. (2009) | 論点1 |
| `barbet2012.md` | Barbet-Massin et al. (2012) | 論点3 |
| `phillips2006.md` | Phillips et al. (2006) | 論点3 |
| `cecconi2016.md` | Cecconi et al. (2016) | 論点3 |
| `narayanan2021.md` | Narayanan et al. (2021) | 論点3 |

### 参考 ☆（完了）

| ファイル | 論文 | 論点 |
|---------|------|------|
| `blei2006.md` | Blei & Jordan (2006) | 論点2 |
| `walker2007.md` | Walker (2007) | 論点2 |
| `quintana2020.md` | Quintana et al. (2020) | 論点2 |

### 基礎文献 ◆（追加）

文献ギャップ分析に基づき追加された基礎理論文献。

| ファイル | 論文 | 論点 | 重要性 |
|---------|------|------|--------|
| `aitchison1982.md` | Aitchison (1982) | 論点2 | **組成データ分析の創始論文**。対数比変換とAitchison幾何学の基礎。eckardt2025.mdで8回以上引用。 |
| `aitchison1983.md` | Aitchison (1983) | 論点2 | 組成データのPCA。主成分分析を単体上で行う手法。 |
| `pawlowsky2001.md` | Pawlowsky-Glahn & Egozcue (2001) | 論点2 | Aitchison幾何学の厳密な定式化。単体をヒルベルト空間として扱う理論的基盤。 |
| `moller1998.md` | Møller et al. (1998) | 論点1・3 | **Log Gaussian Cox過程の原論文**。GPで強度関数をモデル化する手法の理論的基礎。 |
| `diggle1998.md` | Diggle et al. (1998) | 論点1・3 | **「モデルベース地球統計学」の命名論文**。空間統計における尤度ベース推論の基礎。 |
| `vecchia1988.md` | Vecchia (1988) | 論点1 | **Vecchia近似の原論文**。NNGPの理論的基礎。大規模空間データへの対応手法。 |
| `cressie2008.md` | Cressie & Johannesson (2008) | 論点1 | **Fixed Rank Kriging**。NNGPの代替となる大規模空間データ近似手法。 |
| `warton2010.md` | Warton & Shepherd (2010) | 論点3 | **Pseudo-absence問題の理論的解決**。点過程モデリングがpseudo-absenceを不要にすることを証明。 |
| `gneiting2007.md` | Gneiting & Raftery (2007) | 論点4 | **Strictly Proper Scoring Rules**。予測性能評価の理論的基礎。モデル比較の正当化。 |

---

## 論点ごとの文献マップ

### 論点1: 空間的依存性

```
Gaussian Process (GP)
├── 基礎: Rasmussen & Williams (2006)
├── 空間統計
│   ├── Model-based Geostatistics: Diggle et al. (1998) ✓ ◆NEW
│   └── Banerjee et al. (2014)
├── 近似手法
│   ├── Vecchia近似: Vecchia (1988) ✓ ◆NEW
│   ├── NNGP: Datta et al. (2016) ✓
│   ├── Fixed Rank Kriging: Cressie & Johannesson (2008) ✓ ◆NEW
│   └── INLA: Rue et al. (2009) ✓
├── Cox過程
│   └── Log Gaussian Cox: Møller et al. (1998) ✓ ◆NEW
└── 空間変化係数: Gelfand et al. (2003) ✓
    └── spBayes実装: Finley & Banerjee (2020) ✓
```

### 論点2: 組成データ

```
組成データ分析
├── 古典的基礎 ◆NEW
│   ├── 創始論文: Aitchison (1982) ✓ ◆NEW
│   ├── 組成データPCA: Aitchison (1983) ✓ ◆NEW
│   └── Aitchison幾何学: Pawlowsky-Glahn & Egozcue (2001) ✓ ◆NEW
├── Multinomial Logit
│   ├── Pólya-Gamma: Polson et al. (2013) ✓
│   └── Stick-breaking + PG: Linderman et al. (2015) ✓
├── Ultimate PG: Zens et al. (2023) ✓
├── Bayesian nonparametric: Wehrhahn et al. (2022) ✓
└── Dirichlet Process
    ├── 変分推論: Blei & Jordan (2006) ✓
    ├── Slice sampling: Walker (2007) ✓
    └── Dependent DP: Quintana et al. (2020) ✓
```

### 論点3: Presence-only

```
Presence-only モデリング
├── 総説: Renner et al. (2015) ✓
├── 点過程理論
│   ├── Gelfand & Schliep (2018) ✓
│   ├── Gonçalves & Gamerman (2018) ✓
│   └── Log Gaussian Cox: Møller et al. (1998) ✓ ◆NEW
├── ベイズ推論
│   ├── Moreira & Gamerman (2022) ✓
│   └── Moreira et al. (2024) ✓
├── Preferential sampling
│   ├── Gelfand & Shirota (2019) ✓
│   └── Cecconi et al. (2016) ✓
├── 理論
│   ├── Fithian & Hastie (2013) ✓
│   └── Warton & Shepherd (2010) ✓ ◆NEW ← Pseudo-absence不要の理論的証明
├── Imperfect detection: Dorazio (2014) ✓
├── Pseudo-absence: Barbet-Massin et al. (2012) ✓
├── Marked point process: Narayanan et al. (2021) ✓
└── 比較手法: MaxEnt Phillips et al. (2006) ✓
```

### 論点4: 統合の困難さ

```
統合の課題
├── 計算困難性
│   ├── GP + 点過程 + multinomialの同時推論
│   └── 近似手法: Vecchia (1988) ✓, NNGP, Fixed Rank Kriging ◆NEW
├── 識別可能性
│   └── Moreira & Gamerman (2022)で議論 ✓
├── モデル評価
│   └── Scoring Rules: Gneiting & Raftery (2007) ✓ ◆NEW
├── データ統合: Pacifici et al. (2017) ✓
└── Composition marks: Eckardt et al. (2025) ✓
```

---

## 次のアクション

1. **各論文の要約作成**
   - 変換済み論文の要約をテンプレートに沿って作成

2. **論点ごとの議論整理**
   - 各論点について、文献から得られた知見をまとめる

3. **第3章の執筆**
   - 先行研究セクションのドラフト作成
