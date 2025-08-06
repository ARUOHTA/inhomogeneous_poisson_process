"""
可視化の基底クラス

既存の可視化機能を整理・統合（新機能なし）
"""

from abc import ABC, abstractmethod
from typing import Any, Dict

import matplotlib.pyplot as plt
import polars as pl


class BaseVisualizer(ABC):
    """可視化の共通インターフェース"""

    @abstractmethod
    def plot_prediction_map(self, data: pl.DataFrame, **kwargs) -> plt.Figure:
        """
        予測結果の地図プロット

        Parameters
        ----------
        data : pl.DataFrame
            予測結果データ
        **kwargs
            プロット固有のパラメータ

        Returns
        -------
        plt.Figure
            作成された図
        """
        pass

    @abstractmethod
    def plot_model_diagnostics(self, model_results: Dict[str, Any]) -> plt.Figure:
        """
        モデル診断プロット

        Parameters
        ----------
        model_results : Dict[str, Any]
            モデルの学習結果

        Returns
        -------
        plt.Figure
            診断プロット
        """
        pass
