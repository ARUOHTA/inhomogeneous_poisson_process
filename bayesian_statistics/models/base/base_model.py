"""
産地構成比モデルの基底クラス

DRY原則に基づく共通インターフェース定義
既存機能のリファクタリングのみ（新機能追加なし）
"""

from abc import ABC, abstractmethod
from typing import Any, Dict

import polars as pl


class BaseCompositionModel(ABC):
    """産地構成比モデルの共通インターフェース"""

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
            学習結果
        """
        pass

    @abstractmethod
    def predict_site_ratios(self, preprocessor) -> Dict[str, pl.DataFrame]:
        """
        遺跡での産地構成比予測

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ

        Returns
        -------
        Dict[str, pl.DataFrame]
            時期別の遺跡構成比予測結果
        """
        pass

    @abstractmethod
    def predict_grid_ratios(self, preprocessor) -> pl.DataFrame:
        """
        グリッド点での産地構成比予測

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ

        Returns
        -------
        pl.DataFrame
            グリッド点の構成比予測結果
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
        学習結果へのアクセス

        Returns
        -------
        Dict[str, Any]
            学習結果
        """
        pass
