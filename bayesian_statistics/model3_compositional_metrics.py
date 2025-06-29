"""
Compositional Data Analysis用の評価指標
Aitchison Distance, Total Variation等の実装
"""

from typing import Optional, Tuple

import numpy as np


class CompositionalMetrics:
    """Compositional data分析用の評価指標クラス"""
    
    def __init__(self, zero_replacement: float = 1e-6):
        """
        Parameters
        ----------
        zero_replacement : float
            0値を置換する小さな正の値
        """
        self.zero_replacement = zero_replacement
    
    def _validate_composition(self, x: np.ndarray) -> np.ndarray:
        """
        構成比ベクトルの妥当性を確認し、必要に応じて修正
        
        Parameters
        ----------
        x : np.ndarray
            構成比ベクトル
            
        Returns
        -------
        np.ndarray
            修正された構成比ベクトル
        """
        # 負の値を0に置換
        x = np.maximum(x, 0.0)
        
        # 0値を小さな正の値に置換
        x = np.where(x == 0, self.zero_replacement, x)
        
        # 合計が1になるように正規化
        x = x / np.sum(x)
        
        return x
    
    def _geometric_mean(self, x: np.ndarray) -> float:
        """
        幾何平均を計算
        
        Parameters
        ----------
        x : np.ndarray
            正の値のベクトル
            
        Returns
        -------
        float
            幾何平均
        """
        # 対数の平均を取ってから指数を計算
        return np.exp(np.mean(np.log(x)))
    
    def clr_transform(self, x: np.ndarray) -> np.ndarray:
        """
        Centered Log-Ratio (CLR) 変換
        
        Parameters
        ----------
        x : np.ndarray
            構成比ベクトル
            
        Returns
        -------
        np.ndarray
            CLR変換されたベクトル
        """
        x = self._validate_composition(x)
        geom_mean = self._geometric_mean(x)
        return np.log(x / geom_mean)
    
    def aitchison_distance(self, x: np.ndarray, y: np.ndarray) -> float:
        """
        Aitchison距離を計算
        
        Aitchison距離は、CLR変換後のユークリッド距離として定義される
        
        Parameters
        ----------
        x : np.ndarray
            構成比ベクトル1 (sum = 1)
        y : np.ndarray
            構成比ベクトル2 (sum = 1)
            
        Returns
        -------
        float
            Aitchison距離
            
        References
        ----------
        Aitchison, J. (1982). The statistical analysis of compositional data. 
        Journal of the Royal Statistical Society Series B, 44(2), 139-177.
        """
        if len(x) != len(y):
            raise ValueError("構成比ベクトルの次元が一致しません")
        
        # CLR変換を適用
        clr_x = self.clr_transform(x)
        clr_y = self.clr_transform(y)
        
        # ユークリッド距離を計算
        return np.sqrt(np.sum((clr_x - clr_y) ** 2))
    
    def total_variation(self, x: np.ndarray, y: np.ndarray) -> float:
        """
        Total Variationを計算
        
        Total Variation = 1/(2D) * sum_i sum_j (ln(x_i/x_j) - ln(y_i/y_j))^2
        
        Parameters
        ----------
        x : np.ndarray
            構成比ベクトル1 (sum = 1)
        y : np.ndarray
            構成比ベクトル2 (sum = 1)
            
        Returns
        -------
        float
            Total Variation
            
        References
        ----------
        Aitchison, J. (1986). The Statistical Analysis of Compositional Data. 
        Chapman & Hall.
        """
        if len(x) != len(y):
            raise ValueError("構成比ベクトルの次元が一致しません")
        
        x = self._validate_composition(x)
        y = self._validate_composition(y)
        
        D = len(x)
        total_var = 0.0
        
        # 全ての成分ペアについて計算
        for i in range(D):
            for j in range(D):
                if i != j:  # i == j の場合は ln(1) = 0 なので省略可能
                    log_ratio_x = np.log(x[i] / x[j])
                    log_ratio_y = np.log(y[i] / y[j])
                    total_var += (log_ratio_x - log_ratio_y) ** 2
        
        return total_var / (2 * D)
    
    def variation_array(self, compositions: np.ndarray) -> np.ndarray:
        """
        Variation Arrayを計算
        
        V_ij = var(ln(x_i/x_j)) for all compositions
        
        Parameters
        ----------
        compositions : np.ndarray
            構成比行列 (n_samples, n_components)
            
        Returns
        -------
        np.ndarray
            Variation Array (n_components, n_components)
        """
        n_samples, n_components = compositions.shape
        V = np.zeros((n_components, n_components))
        
        # 各サンプルを妥当性チェック
        validated_compositions = np.array([
            self._validate_composition(comp) for comp in compositions
        ])
        
        for i in range(n_components):
            for j in range(n_components):
                if i != j:
                    log_ratios = np.log(validated_compositions[:, i] / validated_compositions[:, j])
                    V[i, j] = np.var(log_ratios)
                else:
                    V[i, j] = 0.0  # 対角成分は0
        
        return V
    
    def compute_all_metrics(
        self, 
        observed: np.ndarray, 
        predicted: np.ndarray
    ) -> dict:
        """
        全ての評価指標を一括計算
        
        Parameters
        ----------
        observed : np.ndarray
            観測された構成比ベクトル
        predicted : np.ndarray
            予測された構成比ベクトル
            
        Returns
        -------
        dict
            全ての評価指標を含む辞書
        """
        metrics = {}
        
        try:
            metrics['aitchison_distance'] = self.aitchison_distance(observed, predicted)
        except Exception as e:
            metrics['aitchison_distance'] = np.nan
            print(f"Aitchison距離の計算でエラー: {e}")
        
        try:
            metrics['total_variation'] = self.total_variation(observed, predicted)
        except Exception as e:
            metrics['total_variation'] = np.nan
            print(f"Total Variationの計算でエラー: {e}")
        
        # 追加の診断情報
        metrics['observed_sum'] = np.sum(observed)
        metrics['predicted_sum'] = np.sum(predicted)
        metrics['observed_validated'] = self._validate_composition(observed)
        metrics['predicted_validated'] = self._validate_composition(predicted)
        
        return metrics


def bray_curtis_dissimilarity(x: np.ndarray, y: np.ndarray) -> float:
    """
    Bray-Curtis Dissimilarity（比較参照用）
    
    BC(x,y) = sum|x_i - y_i| / sum(x_i + y_i)
    
    Parameters
    ----------
    x : np.ndarray
        構成比ベクトル1
    y : np.ndarray
        構成比ベクトル2
        
    Returns
    -------
    float
        Bray-Curtis Dissimilarity
    """
    numerator = np.sum(np.abs(x - y))
    denominator = np.sum(x + y)
    
    if denominator == 0:
        return np.nan
    
    return numerator / denominator


def jensen_shannon_divergence(x: np.ndarray, y: np.ndarray) -> float:
    """
    Jensen-Shannon Divergence（比較参照用）
    
    Parameters
    ----------
    x : np.ndarray
        確率分布1
    y : np.ndarray
        確率分布2
        
    Returns
    -------
    float
        Jensen-Shannon Divergence
    """
    # 正規化
    x = x / np.sum(x)
    y = y / np.sum(y)
    
    # 平均分布
    m = 0.5 * (x + y)
    
    # KL divergenceを計算（0*log(0) = 0と定義）
    def kl_divergence(p, q):
        # 0値の処理
        p = np.where(p == 0, 1e-10, p)
        q = np.where(q == 0, 1e-10, q)
        return np.sum(p * np.log(p / q))
    
    js_div = 0.5 * kl_divergence(x, m) + 0.5 * kl_divergence(y, m)
    
    return js_div