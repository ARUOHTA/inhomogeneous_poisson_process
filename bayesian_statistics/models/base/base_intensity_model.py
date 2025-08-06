"""
遺跡存在確率モデルの基底クラス

IPP（Inhomogeneous Poisson Process）等の空間点過程モデル用
統一インターフェース定義
"""

from abc import ABC, abstractmethod
from typing import Any, Dict

import numpy as np
import polars as pl


class BaseIntensityModel(ABC):
    """遺跡存在確率モデルの共通インターフェース"""

    @abstractmethod
    def fit(self, preprocessor, **kwargs) -> Dict[str, Any]:
        """
        モデルの学習

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ
        **kwargs
            モデル固有のパラメータ

        Returns
        -------
        Dict[str, Any]
            学習結果（MCMCサンプル等）
        """
        pass

    @abstractmethod
    def predict_probability(self, preprocessor) -> pl.DataFrame:
        """
        グリッド上での遺跡存在確率予測

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ

        Returns
        -------
        pl.DataFrame
            グリッド点の存在確率予測結果
            カラム: ['x', 'y', 'probability']
        """
        pass

    @abstractmethod
    def predict_intensity(self, coordinates: np.ndarray) -> np.ndarray:
        """
        指定座標での強度関数値予測

        Parameters
        ----------
        coordinates : np.ndarray
            予測したい座標 (N, 2) - [x, y]

        Returns
        -------
        np.ndarray
            各座標での強度関数値
        """
        pass

    @property
    @abstractmethod
    def model_name(self) -> str:
        """
        モデル名（評価結果の識別用）

        Returns
        -------
        str
            モデル名
        """
        pass

    @property
    @abstractmethod
    def results(self) -> Dict[str, Any]:
        """
        学習結果への統一アクセス

        Returns
        -------
        Dict[str, Any]
            学習結果（パラメータサンプル、予測結果等）
        """
        pass

    def create_inference_data(self):
        """
        ArviZ用診断データ作成（オプション）

        Returns
        -------
        arviz.InferenceData or None
            MCMC診断用データ
        """
        # デフォルト実装では何もしない
        # 必要に応じてサブクラスでオーバーライド
        return None
