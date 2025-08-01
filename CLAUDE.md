<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.min.js">
</script>
<script type="text/x-mathjax-config">
 MathJax.Hub.Config({
 tex2jax: {
 inlineMath: [['$', '$'] ],
 displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
 }
 });
</script>


# ベイズ統計プロジェクトガイド

## オイラーの公式
オイラーの公式は以下のように与えられる。

$$ e^{i x} = \cos{x} + i \sin{x} $$

##  $ \varepsilon - \delta $ 論法
任意の $ \varepsilon > 0 $ についてある $ \delta > 0 $ が存在して、任意の $ x \in \mathbb{R} $ に対して $ 0 < |x - a| < \delta $ ならば $ |f(x) - f(a)| < \varepsilon $ を満たすとき $ f(x) $ は $ a $ で連続であるという。


## プロジェクト概要
このプロジェクトは、考古学的な黒曜石分布データを分析するための、空間的に変化する係数を持つベイズスパース非均質ポアソン過程モデルを実装しています。

## プロジェクト構造
- `bayesian_statistics/` - ベイズ分析のためのコアPythonモジュール
- `notebooks/` - 実験と可視化のためのJupyterノートブック
- `data/` - 黒曜石遺跡と標高データのCSVファイルと距離行列
- `obsidian_data/` - 黒曜石考古学データの生データ
- `bayesPO_code/` - ベイズポアソン過程分析のためのRコード
- `docs/` - LaTeXドキュメントとプレゼンテーション
- `output/` - 生成された図と分析結果

## 主要な依存関係
- 科学計算: numpy, scipy, pandas, polars
- ベイズ分析: pypolyagamma, arviz, GPy
- 地理空間: geopandas, shapely, contextily, keplergl
- 可視化: matplotlib, seaborn, japanize-matplotlib

## 開発環境のセットアップ
プロジェクトは依存関係管理に`uv`を使用しています。依存関係をインストールするには：
```bash
uv sync
```

## コード品質
- リンティング: `ruff check .`
- フォーマット: `ruff format .`
- ノートブック用にnbqaを使用したpre-commitフックが設定済み

## 主要な分析コンポーネント
1. **1Dベイズ IPP**: シンプルな非均質ポアソン過程モデル
2. **2D空間モデル**: 2次元空間データへの拡張
3. **黒曜石データ処理**: 考古学遺跡データの前処理
4. **距離計算**: 地形ベースのルーティングのためのトブラー距離メトリクス
5. **可視化**: 様々なモデル出力と事後分布

## データファイル
- 標高データ: `*_gdf_elevation.csv`
- 遺跡位置: `*_gdf_sites.csv`
- 黒曜石産地: `*_gdf_obsidian.csv`
- 距離行列: `data/16_tobler_distance*/` ディレクトリ

## 主要なノートブック
- `12_model_1st.ipynb` - 最初の空間モデル
- `14_model_2nd.ipynb` - 改良を加えた第2版
- `16_model_3rd.ipynb` - トブラー距離を使用した第3版モデル
- `19_model_4.ipynb` - 第4版モデル
- `20_model_KSBP.ipynb` - カーネルスティックブレーキング過程モデル

## よくあるタスク
- 分析の実行: ノートブックを番号順に実行
- 可視化の生成: `13_`, `15_`, `17_` の可視化ノートブックを確認
- 新しいデータの処理: 前処理ノートブック (`09_`, `10_`, `11_`) から開始

## 黒曜石データについて
プロジェクトは日本の考古学的遺跡における黒曜石の分布を分析しています。主な黒曜石産地：
- 信州（長野県）
- 神津島
- 箱根
- 高原山

## モデルの進化
1. 第1版: 基本的な空間ポアソン過程
2. 第2版: 改良された空間構造
3. 第3版: 地形を考慮したトブラー距離の導入
4. 第4版: さらなる改良
5. KSBP版: より柔軟なノンパラメトリックベイズアプローチ

## Banerjee-Carlin-Gelfand 書籍翻訳プロジェクト
Banerjee, Carlin, and Gelfandの「Hierarchical Modeling and Analysis for Spatial Data」の各章を詳細な日本語要約として作成する際は、必ず以下の手順書を参照すること：

**翻訳手順書の場所**: `/home/ohta/dev/bayesian_statistics/docs/banerjee_carlin_gelfand/INSTRUCTION.md`

各章の翻訳を開始する前に、必ずINSTRUCTION.mdを読み、定められた手順と方針に従うこと。
