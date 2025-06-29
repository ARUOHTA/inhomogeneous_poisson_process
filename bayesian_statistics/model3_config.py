"""
Model 3rd の設定とパイプラインクラス
"""

import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import polars as pl

from .model3_ipp import InhomogeneousPoissonProcess
from .model3_nadaraya_watson import NadarayaWatsonEstimator
from .model3_preprocessing import ObsidianDataPreprocessor


@dataclass
class Model3Config:
    """Model 3rd の設定"""

    # データパス
    data_dir: str

    # グリッド設定
    x_min: float = 138
    x_max: float = 141
    y_min: float = 34
    y_max: float = 37

    # NW推定設定
    nw_sigma: float = 500
    nw_sigma_for_sites: float = 0.1

    # IPP設定
    mcmc_iterations: int = 30000
    burn_in: int = 5000

    # 説明変数（NW推定用）
    nw_variable_names: List[str] = field(
        default_factory=lambda: [
            "average_elevation",
            "average_slope_angle",
            "cost_kouzu",
            "cost_shinshu",
            "cost_hakone",
            "cost_takahara",
            "cost_river",
        ]
    )

    # 説明変数（IPP用）
    ipp_variable_names: List[str] = field(
        default_factory=lambda: [
            "average_elevation",
            "average_slope_angle",
            "cost_shinshu",
            "cost_river",
        ]
    )

    # 時期・産地設定
    time_periods: Dict[int, str] = field(
        default_factory=lambda: {
            0: "早期・早々期",
            1: "前期",
            2: "中期",
            3: "後期",
            4: "晩期",
        }
    )

    origins: List[str] = field(
        default_factory=lambda: ["神津島", "信州", "箱根", "高原山", "その他"]
    )


class Model3Pipeline:
    """Model 3rd のパイプライン"""

    def __init__(self, config: Model3Config):
        """
        Parameters
        ----------
        config : Model3Config
            設定オブジェクト
        """
        self.config = config
        self._preprocessor: Optional[ObsidianDataPreprocessor] = None
        self._nw_estimator: Optional[NadarayaWatsonEstimator] = None
        self._ipp_model: Optional[InhomogeneousPoissonProcess] = None
        self._results: Dict[str, Any] = {}

    def run_preprocessing(self) -> ObsidianDataPreprocessor:
        """
        前処理を実行

        Returns
        -------
        ObsidianDataPreprocessor
            前処理済みのデータ
        """
        print("=== 前処理を開始 ===")

        # 前処理クラスのインスタンス化
        self._preprocessor = ObsidianDataPreprocessor(self.config.data_dir)

        # データの読み込み
        print("データを読み込んでいます...")
        data = self._preprocessor.load_data()

        print(f"標高データ: {data['elevation'].shape}")
        print(f"黒曜石データ: {data['obsidian'].shape}")
        print(f"遺跡データ: {data['sites'].shape}")

        return self._preprocessor

    def run_nadaraya_watson(
        self, preprocessor: Optional[ObsidianDataPreprocessor] = None
    ) -> NadarayaWatsonEstimator:
        """
        Nadaraya-Watson推定を実行

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor, optional
            前処理済みのデータ（Noneの場合は内部のものを使用）

        Returns
        -------
        NadarayaWatsonEstimator
            学習済みの推定器
        """
        if preprocessor is None:
            preprocessor = self._preprocessor

        if preprocessor is None:
            raise ValueError(
                "前処理が実行されていません。run_preprocessing()を先に実行してください。"
            )

        print("\n=== Nadaraya-Watson推定を開始 ===")

        # 推定器のインスタンス化
        self._nw_estimator = NadarayaWatsonEstimator(
            sigma=self.config.nw_sigma,
            sigma_for_sites=self.config.nw_sigma_for_sites,
        )

        # モデルの学習
        print("重み行列を計算しています...")
        self._nw_estimator.fit(preprocessor, self.config.nw_variable_names)

        # 全時期・全産地の推定
        print("\n全時期・全産地の推定を実行しています...")
        ratio_df, ratio_sites_df = self._nw_estimator.predict_all_periods_origins(
            preprocessor,
            self.config.time_periods,
            self.config.origins,
        )

        # 結果の保存
        self._results["nw_ratio_df"] = ratio_df
        self._results["nw_ratio_sites_df"] = ratio_sites_df

        # データフレームの更新
        preprocessor._df_elevation = preprocessor.df_elevation.join(
            ratio_df, on=["x", "y"]
        )
        preprocessor._df_sites = preprocessor.df_sites.join(ratio_sites_df, on="遺跡ID")

        return self._nw_estimator

    def run_ipp(
        self, preprocessor: Optional[ObsidianDataPreprocessor] = None
    ) -> InhomogeneousPoissonProcess:
        """
        IPP推定を実行

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor, optional
            前処理済みのデータ（Noneの場合は内部のものを使用）

        Returns
        -------
        InhomogeneousPoissonProcess
            学習済みのモデル
        """
        if preprocessor is None:
            preprocessor = self._preprocessor

        if preprocessor is None:
            raise ValueError(
                "前処理が実行されていません。run_preprocessing()を先に実行してください。"
            )

        print("\n=== IPP推定を開始 ===")

        # モデルのインスタンス化
        region = [
            [self.config.x_min, self.config.x_max],
            [self.config.y_min, self.config.y_max],
        ]
        self._ipp_model = InhomogeneousPoissonProcess(
            variable_names=self.config.ipp_variable_names,
            region=region,
        )

        # モデルの学習
        print("MCMCサンプリングを実行しています...")
        mcmc_results = self._ipp_model.fit(
            preprocessor,
            num_iterations=self.config.mcmc_iterations,
            burn_in=self.config.burn_in,
        )

        # 遺跡存在確率の予測
        print("遺跡存在確率を計算しています...")
        grid_coords = preprocessor.create_grid_coords()
        site_probability = self._ipp_model.predict_site_probability(grid_coords)

        # 結果の保存
        self._results["ipp_beta_samples"] = mcmc_results["beta_samples"]
        self._results["ipp_lambda_star_samples"] = mcmc_results["lambda_star_samples"]
        self._results["site_probability"] = site_probability

        # データフレームの更新
        preprocessor._df_elevation = preprocessor.df_elevation.with_columns(
            pl.Series("site_probability", site_probability)
        )

        return self._ipp_model

    def run_full_pipeline(self) -> Dict[str, Any]:
        """
        フルパイプラインを実行

        Returns
        -------
        Dict[str, Any]
            実行結果の辞書
        """
        # 前処理
        preprocessor = self.run_preprocessing()

        # NW推定
        nw_estimator = self.run_nadaraya_watson(preprocessor)

        # IPP推定
        ipp_model = self.run_ipp(preprocessor)

        print("\n=== パイプライン実行完了 ===")

        return {
            "preprocessor": preprocessor,
            "nw_estimator": nw_estimator,
            "ipp_model": ipp_model,
            "results": self._results,
        }

    def save_results(self, output_dir: str):
        """
        結果を保存

        Parameters
        ----------
        output_dir : str
            出力ディレクトリ
        """
        if self._preprocessor is None:
            raise ValueError("パイプラインが実行されていません。")

        print(f"\n結果を保存しています: {output_dir}")

        # ディレクトリの作成
        os.makedirs(output_dir, exist_ok=True)

        # データフレームの保存
        self._preprocessor.df_elevation.write_csv(
            os.path.join(output_dir, "16_gdf_elevation_with_ratio_2nd.csv")
        )
        self._preprocessor.df_sites.write_csv(
            os.path.join(output_dir, "16_gdf_sites_with_ratio_2nd.csv")
        )

        print("保存が完了しました。")

    @property
    def preprocessor(self) -> ObsidianDataPreprocessor:
        """前処理オブジェクト"""
        if self._preprocessor is None:
            raise ValueError("前処理が実行されていません。")
        return self._preprocessor

    @property
    def nw_estimator(self) -> NadarayaWatsonEstimator:
        """NW推定器"""
        if self._nw_estimator is None:
            raise ValueError("NW推定が実行されていません。")
        return self._nw_estimator

    @property
    def ipp_model(self) -> InhomogeneousPoissonProcess:
        """IPPモデル"""
        if self._ipp_model is None:
            raise ValueError("IPP推定が実行されていません。")
        return self._ipp_model

    @property
    def results(self) -> Dict[str, Any]:
        """実行結果"""
        return self._results
