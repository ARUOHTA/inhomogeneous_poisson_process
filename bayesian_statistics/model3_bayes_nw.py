"""
ベイズ拡張Nadaraya-Watson推定量（簡素化版）

データがある遺跡のみを対象とする推定を標準動作とする
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl
from tqdm import tqdm

from .model3_nadaraya_watson import NadarayaWatsonEstimator
from .model3_preprocessing import ObsidianDataPreprocessor


class BayesianNadarayaWatson(NadarayaWatsonEstimator):
    """
    ベイズ拡張Nadaraya-Watson推定量

    データがある遺跡のみを対象とする推定を標準動作とし、
    0.25問題を根本的に解決する
    """

    def __init__(
        self,
        sigma: float,
        sigma_for_sites: float,
        alpha_0: float = 100.0,
        gamma_0: float = 1.0,
    ):
        """
        Parameters
        ----------
        sigma : float
            空間距離のバンド幅（km）
        sigma_for_sites : float
            共変量のバンド幅
        alpha_0 : float, default=100.0
            精度パラメータ（大きいほどNW推定量に近づく）
        gamma_0 : float, default=1.0
            事前の擬似カウント総数
        """
        super().__init__(sigma, sigma_for_sites)
        self.alpha_0 = alpha_0
        self.gamma_0 = gamma_0

        # ベイズ推定用の内部変数
        self._variable_names: Optional[List[str]] = None

    def fit(
        self, preprocessor: ObsidianDataPreprocessor, variable_names: List[str]
    ) -> Dict[str, np.ndarray]:
        """
        モデルを学習（既存のNW推定量と同じ重み計算）

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ
        variable_names : List[str]
            使用する説明変数名

        Returns
        -------
        Dict[str, np.ndarray]
            計算された重み行列
        """
        # 親クラスの重み計算をそのまま実行
        weights_dict = super().fit(preprocessor, variable_names)

        # 重み行列のNaNチェックと修正
        self._fix_nan_weights()

        # 変数名を保存（予測時に使用）
        self._variable_names = variable_names

        return weights_dict

    def _fix_nan_weights(self):
        """重み行列のNaN確認（置換はしない）"""
        if self._weights is not None:
            nan_mask = np.isnan(self._weights)
            nan_count = nan_mask.sum()

            if nan_count > 0:
                print(
                    f"情報: グリッド重み行列に{nan_count}個のNaNを検出（海上部分として保持）"
                )

        if self._weights_sites is not None:
            nan_mask = np.isnan(self._weights_sites)
            nan_count = nan_mask.sum()

            if nan_count > 0:
                print(f"警告: 遺跡重み行列に{nan_count}個のNaNを検出")
                # 遺跡重みのNaNは問題なので0で置換
                self._weights_sites = np.nan_to_num(self._weights_sites, nan=0.0)

    def _get_all_sources_data(
        self, preprocessor: ObsidianDataPreprocessor, target_period: int
    ) -> np.ndarray:
        """全産地の出土数を取得（正しい産地名マッピング使用）"""
        # 実際のデータに基づく産地名マッピング
        source_mapping = {
            "神津島": ["神津島"],
            "信州": ["諏訪", "和田峠", "蓼科", "男女倉"],
            "箱根": ["箱根", "天城"],
            "高原山": ["高原山"],
        }

        sources = ["神津島", "信州", "箱根", "高原山"]

        # 全遺跡IDのリストを取得
        max_site_id = preprocessor.df_obsidian["遺跡ID"].max()
        n_sites = max_site_id + 1

        # 各産地の出土数を格納する行列
        obsidian_counts = np.zeros((n_sites, len(sources)))

        # 対象時期のデータのみ抽出
        period_df = preprocessor.df_obsidian.filter(pl.col("時期") == target_period)

        for j, source in enumerate(sources):
            # 該当する実際の産地名のリストを取得
            actual_source_names = source_mapping[source]

            # 各実際の産地名でフィルタして合計
            total_counts = np.zeros(n_sites)

            for actual_source in actual_source_names:
                source_counts = (
                    period_df.filter(pl.col("産地") == actual_source)
                    .group_by("遺跡ID")
                    .agg([pl.len().alias("count")])
                    .join(
                        pl.DataFrame({"遺跡ID": np.arange(n_sites)}),
                        on="遺跡ID",
                        how="right",
                    )
                    .fill_null(0)
                    .sort("遺跡ID")["count"]
                    .to_numpy()
                )

                total_counts += source_counts

            obsidian_counts[:, j] = total_counts

        return obsidian_counts

    def calculate_base_measure(
        self,
        weights: np.ndarray,
        preprocessor: ObsidianDataPreprocessor,
        target_period: int,
    ) -> np.ndarray:
        """
        基準測度π*(s)を計算（元のNW推定量ベース、全産地同時計算）

        Parameters
        ----------
        weights : np.ndarray
            カーネル重み (n_points, n_sites) または (n_sites,)
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ
        target_period : int
            対象時期

        Returns
        -------
        np.ndarray
            正規化された基準測度 (n_points, n_sources) または (n_sources,)
            各行の合計は1.0
        """
        sources = ["神津島", "信州", "箱根", "高原山"]
        n_sources = len(sources)

        # 元のNW推定量の方法で各産地の比率を計算
        source_ratios = []

        for source in sources:
            # 対象産地のデータ取得（元のNW推定量の方法）
            counts, target_counts = preprocessor.preprocess_obsidian_data(
                target_period, source
            )

            if weights.ndim == 1:
                # 単一地点の場合
                weighted_total = np.sum(weights * counts)
                weighted_target = np.sum(weights * target_counts)

                # 比率計算（0除算を防ぐ）
                if weighted_total > 0:
                    ratio = weighted_target / weighted_total
                else:
                    ratio = 0.0

                source_ratios.append(ratio)

            else:
                # 複数地点の場合
                weighted_total = np.sum(weights * counts, axis=1)  # (n_points,)
                weighted_target = np.sum(weights * target_counts, axis=1)  # (n_points,)

                # 比率計算（0除算を防ぐ、NaN処理も含む）
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratios = np.where(
                        weighted_total > 0, weighted_target / weighted_total, 0
                    )
                    # NaNを0で置換
                    ratios = np.nan_to_num(ratios, nan=0.0)

                source_ratios.append(ratios)

        if weights.ndim == 1:
            # 単一地点の場合
            raw_ratios = np.array(source_ratios)  # (n_sources,)

            # 擬似カウントによる平滑化
            if self.gamma_0 > 0:
                pseudo_weight = self.gamma_0 / n_sources
                raw_ratios += pseudo_weight

            # 正規化
            total = raw_ratios.sum()
            if total > 0:
                result = raw_ratios / total
            else:
                result = np.ones(n_sources) / n_sources

            return result

        else:
            # 複数地点の場合
            raw_ratios = np.array(source_ratios).T  # (n_points, n_sources)

            # 擬似カウントによる平滑化
            if self.gamma_0 > 0:
                pseudo_weight = self.gamma_0 / n_sources
                raw_ratios += pseudo_weight

            # 各地点で正規化
            row_sums = raw_ratios.sum(axis=1, keepdims=True)

            # 有効な行のマスク
            valid_mask = row_sums.flatten() > 0

            result = np.zeros_like(raw_ratios)
            if valid_mask.any():
                result[valid_mask] = raw_ratios[valid_mask] / row_sums[valid_mask]

            # 合計がゼロの行は均等分布
            invalid_mask = ~valid_mask
            if invalid_mask.any():
                result[invalid_mask] = 1.0 / n_sources

            return result

    def predict_single(
        self,
        preprocessor: ObsidianDataPreprocessor,
        target_period: int,
        target_origin: str,
    ) -> Dict[str, np.ndarray]:
        """
        単一時期・産地を予測（データがある遺跡のみ）

        これが新しい標準動作です

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ
        target_period : int
            対象時期
        target_origin : str
            対象産地

        Returns
        -------
        Dict[str, np.ndarray]
            ratio_sites: 遺跡での比率（データがある遺跡のみ）
            site_ids: 対応する遺跡ID
            ratio_mesh: グリッドでの比率（通常通り全グリッド）
            n_sites_with_data: データがある遺跡数
            coverage: データカバレッジ率
        """
        if self._weights is None or self._weights_sites is None:
            raise ValueError(
                "モデルが学習されていません。fit()を先に実行してください。"
            )

        # 対象時期でデータがある遺跡を特定
        period_data = preprocessor.df_obsidian.filter(pl.col("時期") == target_period)
        sites_with_data = period_data["遺跡ID"].unique().sort().to_numpy()

        if len(sites_with_data) == 0:
            raise ValueError(f"時期{target_period}にはデータがありません")

        # 全産地の基準測度を計算
        pi_star_grid = self.calculate_base_measure(
            self._weights, preprocessor, target_period
        )
        pi_star_sites = self.calculate_base_measure(
            self._weights_sites, preprocessor, target_period
        )

        # 対象産地のインデックスを取得
        sources = ["神津島", "信州", "箱根", "高原山"]
        target_idx = sources.index(target_origin)

        # グリッドデータ（全グリッドポイント）
        lon_mesh, lat_mesh = preprocessor.create_meshgrid()
        pi_star_grid_target = pi_star_grid[:, target_idx]

        # データがある遺跡のみを抽出
        pi_star_sites_filtered = pi_star_sites[sites_with_data, target_idx]

        # グリッドデータを形状変換
        if pi_star_grid_target.shape[0] == lon_mesh.size:
            ratio_mesh = pi_star_grid_target.reshape(lon_mesh.shape)
        else:
            ratio_mesh = pi_star_grid_target

        # NaNを0で置換（推定不可能な場所）
        # これは元のNW推定量の動作と一致させるため
        ratio_mesh = np.nan_to_num(ratio_mesh, nan=0.0)

        return {
            "ratio_mesh": ratio_mesh,  # 全グリッド
            "ratio_sites": pi_star_sites_filtered,  # データがある遺跡のみ
            "site_ids": sites_with_data,  # 対応する遺跡ID
            "n_sites_with_data": len(sites_with_data),
            "coverage": len(sites_with_data) / len(pi_star_sites) * 100,
        }

    def predict_all_periods_origins(
        self,
        preprocessor: ObsidianDataPreprocessor,
        time_period_names: Dict[int, str],
        origin_order: List[str],
    ) -> Tuple[pl.DataFrame, Dict[int, pl.DataFrame]]:
        """
        全時期・全産地を予測（データがある遺跡のみ）

        これが新しい標準動作です

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ
        time_period_names : Dict[int, str]
            時期名の辞書
        origin_order : List[str]
            産地のリスト（最後の要素\"その他\"は除外）

        Returns
        -------
        ratio_df : pl.DataFrame
            グリッド上の比率（通常通り）
        ratio_sites_dict : Dict[int, pl.DataFrame]
            各時期のデータがある遺跡での比率
        """
        # グリッドデータは通常通り
        lon_mesh, lat_mesh = preprocessor.create_meshgrid()
        ratio_df = pl.DataFrame({"x": lon_mesh.ravel(), "y": lat_mesh.ravel()})

        # 時期ごとの遺跡データ辞書
        ratio_sites_dict = {}

        for target_period in tqdm(time_period_names.keys(), desc="時期"):
            # この時期でデータがある遺跡を特定
            period_data = preprocessor.df_obsidian.filter(
                pl.col("時期") == target_period
            )
            sites_with_data = period_data["遺跡ID"].unique().sort()

            if len(sites_with_data) == 0:
                print(f"警告: 時期{target_period}にはデータがありません")
                continue

            # この時期の遺跡データフレームを初期化
            period_sites_df = pl.DataFrame({"遺跡ID": sites_with_data})

            for target_origin in origin_order[:-1]:  # "その他"を除外
                print(
                    f"時期{target_period}, 産地{target_origin}: {len(sites_with_data)}遺跡"
                )

                # 予測実行
                results = self.predict_single(
                    preprocessor, target_period, target_origin
                )

                # グリッドデータの追加（通常通り）
                ratio_df = ratio_df.join(
                    pl.DataFrame(
                        {
                            "x": lon_mesh.ravel(),
                            "y": lat_mesh.ravel(),
                            f"ratio_{target_period}_{target_origin}": results[
                                "ratio_mesh"
                            ].ravel(),
                        }
                    ),
                    on=["x", "y"],
                )

                # 遺跡データの追加（データがある遺跡のみ）
                period_sites_df = period_sites_df.join(
                    pl.DataFrame(
                        {
                            "遺跡ID": results["site_ids"],
                            f"比率_{target_period}_{target_origin}": results[
                                "ratio_sites"
                            ],
                        }
                    ),
                    on="遺跡ID",
                )

            # この時期の結果を保存
            ratio_sites_dict[target_period] = period_sites_df

        return ratio_df, ratio_sites_dict

    def sample_posterior(
        self,
        preprocessor: ObsidianDataPreprocessor,
        target_period: int,
        target_origin: str,
        n_samples: int = 1000,
    ) -> Dict[str, np.ndarray]:
        """
        事後分布からサンプリング

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ
        target_period : int
            対象時期
        target_origin : str
            対象産地
        n_samples : int, default=1000
            サンプル数

        Returns
        -------
        Dict[str, np.ndarray]
            事後分布からのサンプル
        """
        from scipy.stats import dirichlet
        from tqdm import tqdm

        if self._weights is None or self._weights_sites is None:
            raise ValueError(
                "モデルが学習されていません。fit()を先に実行してください。"
            )

        # 全産地の基準測度を計算
        pi_star_grid = self.calculate_base_measure(
            self._weights, preprocessor, target_period
        )
        pi_star_sites = self.calculate_base_measure(
            self._weights_sites, preprocessor, target_period
        )

        # ディリクレ分布のパラメータ（全産地）
        alpha_grid = self.alpha_0 * pi_star_grid
        alpha_sites = self.alpha_0 * pi_star_sites

        # グリッド上でのサンプリング
        n_grid_points = alpha_grid.shape[0]
        n_sources = alpha_grid.shape[1]
        samples_grid = np.zeros((n_samples, n_grid_points, n_sources))

        for i in tqdm(range(n_grid_points), desc="グリッドサンプリング"):
            if np.sum(alpha_grid[i]) > 0 and not np.any(np.isnan(alpha_grid[i])):
                samples_grid[:, i, :] = dirichlet.rvs(alpha_grid[i], size=n_samples)
            else:
                # 無効な場合は均等分布
                samples_grid[:, i, :] = 1.0 / n_sources

        # 遺跡でのサンプリング
        n_sites = alpha_sites.shape[0]
        samples_sites = np.zeros((n_samples, n_sites, n_sources))

        for i in tqdm(range(n_sites), desc="遺跡サンプリング"):
            if np.sum(alpha_sites[i]) > 0 and not np.any(np.isnan(alpha_sites[i])):
                samples_sites[:, i, :] = dirichlet.rvs(alpha_sites[i], size=n_samples)
            else:
                # 無効な場合は均等分布
                samples_sites[:, i, :] = 1.0 / n_sources

        # 対象産地のインデックスを取得
        sources = ["神津島", "信州", "箱根", "高原山"]
        target_idx = sources.index(target_origin)

        return {
            "samples_grid": samples_grid,  # 全産地のサンプル (n_samples, n_grid_points, n_sources)
            "samples_sites": samples_sites,  # 全産地のサンプル (n_samples, n_sites, n_sources)
            "target_idx": target_idx,  # 対象産地のインデックス
            "sources": sources,  # 産地名リスト
        }
