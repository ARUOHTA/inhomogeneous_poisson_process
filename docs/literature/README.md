# 文献レビュー

第3章（関連研究）執筆のための文献整理。

## 論文管理ワークフロー

### ツール構成

| ツール | 役割 | 場所 |
|--------|------|------|
| **Paperpile** | PDF論文の管理 | Google Drive (`My Drive/Paperpile/`) |
| **Mathpix CLI** | PDF → Markdown変換 | `mpx convert` (Node.js v20) |
| **変換スクリプト** | 簡易変換 | `scripts/pdf2md.sh` |
| **BibTeX同期** | 書誌情報の自動同期 | `Paperpile/paperpile.bib` |

### PDF変換手順

```bash
# 方法1: スクリプトを使用（推奨）
./scripts/pdf2md.sh "論文のパス.pdf"

# 方法2: Paperpileのショートカット（著者フォルダ/ファイル名）
./scripts/pdf2md.sh "P/Polson et al. 2012 - Bayesian inference....pdf"

# 方法3: 直接コマンド（Node.js v20 + script必須）
nvm use 20
script -q /dev/null mpx convert "input.pdf" "output.mmd"
```

### 注意事項

- Mathpix CLIはNode.js v20で動作（v22では互換性問題あり）
- `script -q /dev/null` はTTYエミュレーションに必要
- APIキーは `.env` に設定済み（`MATHPIX_APP_ID`, `MATHPIX_APP_KEY`）
- 変換後の `.md` ファイルはMathpix Markdown形式（LaTeX数式対応）

### BibTeX運用

- Paperpileの書誌情報は `Paperpile/paperpile.bib` に自動同期
- 必要な項目を `bibliography.bib` にコピーして使用
- Paperpileで論文を追加/編集すると `.bib` も自動更新

---

## 目的

- 各論文の貢献と限界を明確化
- 本研究（MMCP）との関係を整理
- 第3章の論点に対応する情報を抽出

## 論点との対応

| 論点 | 関連文献 |
|------|---------|
| 論点1: 空間的依存性 | Datta et al. (2016), ... |
| 論点2: 組成データ | Aitchison (1982), Polson et al. (2013), ... |
| 論点3: Presence-only | Moreira & Gamerman (2022), ... |
| 論点4: 統合の困難さ | （複数文献から総合） |
| 論点5: 本研究の位置づけ | （複数文献から総合） |

## 文献一覧（35件）

### 最優先 ★★★（10件）

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

### 高優先 ★★（6件）

| ファイル | 論文 | 論点 |
|---------|------|------|
| `fithian2013.md` | Fithian & Hastie (2013) | 論点3 |
| `dorazio2014.md` | Dorazio (2014) | 論点3 |
| `pacifici2017.md` | Pacifici et al. (2017) | 論点4 |
| `zens2023.md` | Zens et al. (2023) | 論点2 |
| `wehrhahn2022.md` | Wehrhahn et al. (2022) | 論点2 |
| `adams2009.md` | Adams et al. (2009) | 論点1・3 |

### 中優先 ★（7件）

| ファイル | 論文 | 論点 |
|---------|------|------|
| `gelfand2003.md` | Gelfand et al. (2003) | 論点1 |
| `finley2020.md` | Finley & Banerjee (2020) | 論点1 |
| `rue2009.md` | Rue et al. (2009) | 論点1 |
| `barbet2012.md` | Barbet-Massin et al. (2012) | 論点3 |
| `phillips2006.md` | Phillips et al. (2006) | 論点3 |
| `cecconi2016.md` | Cecconi et al. (2016) | 論点3 |
| `narayanan2021.md` | Narayanan et al. (2021) | 論点3 |

### 参考 ☆（3件）

| ファイル | 論文 | 論点 |
|---------|------|------|
| `blei2006.md` | Blei & Jordan (2006) | 論点2 |
| `walker2007.md` | Walker (2007) | 論点2 |
| `quintana2020.md` | Quintana et al. (2020) | 論点2 |

### 基礎文献 ◆（9件）

文献ギャップ分析により追加された基礎理論文献。

| ファイル | 論文 | 論点 |
|---------|------|------|
| `aitchison1982.md` | Aitchison (1982) - 組成データ分析の創始 | 論点2 |
| `aitchison1983.md` | Aitchison (1983) - 組成データPCA | 論点2 |
| `pawlowsky2001.md` | Pawlowsky-Glahn & Egozcue (2001) - Aitchison幾何学 | 論点2 |
| `moller1998.md` | Møller et al. (1998) - Log Gaussian Cox過程 | 論点1・3 |
| `diggle1998.md` | Diggle et al. (1998) - Model-based Geostatistics | 論点1・3 |
| `vecchia1988.md` | Vecchia (1988) - Vecchia近似 | 論点1 |
| `cressie2008.md` | Cressie & Johannesson (2008) - Fixed Rank Kriging | 論点1 |
| `warton2010.md` | Warton & Shepherd (2010) - Pseudo-absence問題の解決 | 論点3 |
| `gneiting2007.md` | Gneiting & Raftery (2007) - Scoring Rules | 論点4 |

詳細な優先度リストは `PRIOR_WORK.md` を参照。

## テンプレート

各論文のMarkdownは以下の構造で作成：

```markdown
# [著者名] ([年]) - [短いタイトル]

## 基本情報
- **タイトル**:
- **著者**:
- **ジャーナル**:
- **引用キー**: (bibliography.bibでのキー)

## 主な貢献
- （この論文が何を解決したか）

## 手法の概要
- （技術的な内容の要約）

## 限界・残された課題
- （著者自身が述べている限界、または私たちが認識する限界）

## 本研究との関係
- （MMCPのどの部分に関係するか、何を借りて何を拡張するか）

## 引用すべき箇所
- （論文で引用する際に使えそうな文・主張）
```
