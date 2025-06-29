"""
Inhomogeneous Poisson Process (IPP) クラス
遺跡の存在確率を推定する
"""

from typing import Dict, List, Optional, Tuple, Union

import arviz as az
import numpy as np
import polars as pl
from pypolyagamma import PyPolyaGamma
from tqdm import tqdm

from .model3_preprocessing import ObsidianDataPreprocessor


class IntensityFunction:
    """強度関数を管理するクラス"""
    
    def __init__(self, design_matrix_func, beta: np.ndarray, lambda_star: float):
        self.design_matrix_func = design_matrix_func
        self.beta = beta
        self.lambda_star = lambda_star
    
    def q(self, x: np.ndarray) -> np.ndarray:
        """q(x) を計算"""
        X = self.design_matrix_func(x)
        eta = X @ self.beta
        return self.link_function(eta)
    
    @staticmethod
    def link_function(eta: np.ndarray) -> np.ndarray:
        """ロジスティックリンク関数"""
        return 1 / (1 + np.exp(-eta))
    
    def update_beta(self, beta: np.ndarray):
        """β パラメータを更新"""
        self.beta = beta
    
    def update_lambda_star(self, lambda_star: float):
        """λ* パラメータを更新"""
        self.lambda_star = lambda_star
    
    def copy(self):
        """コピーを作成"""
        return IntensityFunction(
            self.design_matrix_func, self.beta.copy(), self.lambda_star
        )


class InhomogeneousPoissonProcess:
    """非斉次ポアソン過程による遺跡存在確率の推定"""

    def __init__(self, variable_names: List[str], region: List[List[float]]):
        """
        Parameters
        ----------
        variable_names : List[str]
            使用する説明変数名
        region : List[List[float]]
            領域 [[x_min, x_max], [y_min, y_max]]
        """
        self.variable_names = variable_names
        self.region = region
        self._beta_samples: Optional[np.ndarray] = None
        self._lambda_star_samples: Optional[np.ndarray] = None
        self._intensity_func: Optional[IntensityFunction] = None
        self._W_grids: Optional[np.ndarray] = None
        self._preprocessor: Optional[ObsidianDataPreprocessor] = None

    def create_design_matrix(self, xy: np.ndarray) -> np.ndarray:
        """
        デザイン行列を作成

        Parameters
        ----------
        xy : np.ndarray
            座標データ（(N, 2)）

        Returns
        -------
        np.ndarray
            デザイン行列（(N, p+1)）
        """
        if self._W_grids is None or self._preprocessor is None:
            raise ValueError("モデルが初期化されていません。fit()を先に実行してください。")
        
        index = self._convert_to_grid_indices(xy)
        
        # 切片を追加してデザイン行列を作成
        return np.column_stack([np.ones(len(index)), self._W_grids[index]])

    def _convert_to_grid_indices(self, points: np.ndarray) -> np.ndarray:
        """座標をグリッドのインデックスに変換"""
        if self._preprocessor is None:
            raise ValueError("モデルが初期化されていません。fit()を先に実行してください。")
        
        grid_info = self._preprocessor.create_grid_info()
        x, y = points[:, 0], points[:, 1]
        
        index_list = self._get_grid_xy(x, y, grid_info)
        
        x_index = index_list[:, 0]
        y_index = index_list[:, 1]
        
        return y_index * grid_info["n_grid_x"] + x_index

    @staticmethod
    def _get_grid_xy(
        x: Union[float, np.ndarray],
        y: Union[float, np.ndarray],
        grid_info: dict,
    ) -> Union[Optional[tuple[int, int]], np.ndarray]:
        """(x, y) が格子のどこに含まれるかを O(1) で求める"""
        x_min_c = grid_info["x_min_center"]
        x_max_c = grid_info["x_max_center"]
        y_min_c = grid_info["y_min_center"]
        y_max_c = grid_info["y_max_center"]
        dx = grid_info["delta_x"]
        dy = grid_info["delta_y"]

        # メッシュ全体の外枠
        x_left_bound = x_min_c - dx / 2
        x_right_bound = x_max_c + dx / 2
        y_bottom_bound = y_min_c - dy / 2
        y_top_bound = y_max_c + dy / 2

        # x, y を ndarray 化
        x_arr = np.atleast_1d(np.asarray(x))
        y_arr = np.atleast_1d(np.asarray(y))

        if x_arr.shape != y_arr.shape:
            raise ValueError("x と y の shape が一致していません。")

        # 範囲判定
        in_bounds_mask = (
            (x_left_bound <= x_arr)
            & (x_arr <= x_right_bound)
            & (y_bottom_bound <= y_arr)
            & (y_arr <= y_top_bound)
        )

        # grid_x, grid_y を計算
        gx_arr = np.floor((x_arr - x_left_bound) / dx).astype(int)
        gy_arr = np.floor((y_arr - y_bottom_bound) / dy).astype(int)

        # 範囲外は -1 にしておく
        gx_arr[~in_bounds_mask] = -1
        gy_arr[~in_bounds_mask] = -1

        # スカラーだった場合
        if x_arr.ndim == 0:
            if not in_bounds_mask:
                return None
            return (gx_arr.item(), gy_arr.item())

        # 配列の場合
        return np.stack([gx_arr, gy_arr], axis=-1)

    def poisson_thinning_U(
        self,
        intensity_func: IntensityFunction,
        valid_grids: np.ndarray,
        volume: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Poisson thinning を用いて U をサンプリング"""
        # 領域の計算
        mins = np.array([self.region[i][0] for i in range(len(self.region))])
        maxs = np.array([self.region[i][1] for i in range(len(self.region))])

        # λ* × 領域の体積
        lambda_total = intensity_func.lambda_star * volume
        N = np.random.poisson(lambda_total)

        # validなグリッド内からN個のポイントを生成
        candidate_points = np.empty((0, len(self.region)))
        while len(candidate_points) < N:
            # バッチサイズは残りの必要数
            batch_size = N - len(candidate_points)
            # 一様サンプリング
            candidates = np.random.uniform(mins, maxs, size=(batch_size, len(self.region)))

            # グリッドインデックスに変換してvalidかチェック
            grid_indices = self._convert_to_grid_indices(candidates)
            valid_mask = (valid_grids[grid_indices].T)[0]
            # validなものだけ追加
            if len(valid_mask) > 0:
                candidates = candidates[valid_mask]
            candidate_points = np.vstack([candidate_points, candidates])

        # 正確にN個に調整
        candidate_points = np.array(candidate_points[:N])

        # q(x) を計算
        q_candidates = intensity_func.q(candidate_points)
        u = np.random.uniform(size=N)
        # (1 - q(x)) でフィルタリング
        U = candidate_points[u < (1 - q_candidates)]

        return U, candidate_points

    @staticmethod
    def _check_nan_inf(variable, name):
        """NaNやInfのチェック"""
        if np.isnan(variable).any():
            print(f"{name} has NaN")
        if np.isinf(variable).any():
            print(f"{name} has Inf")

    def mcmc_sampler(
        self,
        intensity_func: IntensityFunction,
        y_obs: np.ndarray,
        X_obs: np.ndarray,
        valid_grids: np.ndarray,
        volume: float,
        num_iterations: int,
        prior_beta_mean: np.ndarray,
        prior_beta_cov: np.ndarray,
        prior_lambda_shape: float,
        prior_lambda_rate: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """MCMC サンプラー"""
        # 初期化
        beta_samples = []
        lambda_star_samples = []
        pg = PyPolyaGamma()
        beta = intensity_func.beta
        lambda_star = intensity_func.lambda_star

        for iteration in tqdm(range(num_iterations)):
            # 1. U のサンプリング
            U_events, _ = self.poisson_thinning_U(intensity_func, valid_grids, volume)
            n_U = len(U_events)

            # 2. 結合データの作成
            X_combined = np.concatenate([X_obs, U_events]) if n_U > 0 else X_obs
            self._check_nan_inf(X_combined, "X_combined")
            W_combined = intensity_func.design_matrix_func(X_combined)
            self._check_nan_inf(W_combined, "W_combined")
            self._check_nan_inf(intensity_func.design_matrix_func(X_obs), "W_obs")
            y_combined = np.concatenate([y_obs, np.zeros(n_U)]) if n_U > 0 else y_obs
            n_combined = len(y_combined)

            # 3. ω のサンプリング
            psi = W_combined @ beta
            omega = self._multi_pgdraw_vectorized(pg, np.ones(n_combined), psi)

            # 4. β のサンプリング
            Omega = np.diag(omega)
            z_tilde = (y_combined - 0.5) / omega
            V_inv = np.linalg.inv(prior_beta_cov) + W_combined.T @ Omega @ W_combined
            V = np.linalg.inv(V_inv)

            m = V @ (
                np.linalg.inv(prior_beta_cov) @ prior_beta_mean
                + W_combined.T @ Omega @ z_tilde
            )
            beta = np.random.multivariate_normal(m, V)
            intensity_func.update_beta(beta)

            # 5. λ* のサンプリング
            n_total = n_combined
            lambda_shape = prior_lambda_shape + n_total
            lambda_rate = prior_lambda_rate + np.prod(
                [self.region[i][1] - self.region[i][0] for i in range(len(self.region))]
            )
            lambda_star = np.random.gamma(shape=lambda_shape, scale=1 / lambda_rate)
            intensity_func.update_lambda_star(lambda_star)

            # サンプルの保存
            beta_samples.append(beta)
            lambda_star_samples.append(lambda_star)

        return beta_samples, lambda_star_samples

    @staticmethod
    def _multi_pgdraw_vectorized(pg, a, c):
        """ベクトル化されたPolyaGamma分布からのサンプリング"""
        # polyagammaはベクトル化をサポートしていないので、ループで実装
        n = len(a)
        samples = np.zeros(n)
        for i in range(n):
            samples[i] = pg.pgdraw(a[i], c[i])
        return samples

    def fit(
        self,
        preprocessor: ObsidianDataPreprocessor,
        num_iterations: int = 30000,
        burn_in: int = 5000,
    ) -> Dict[str, np.ndarray]:
        """
        モデルを学習

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ
        num_iterations : int
            MCMCのイテレーション数
        burn_in : int
            バーンイン期間

        Returns
        -------
        Dict[str, np.ndarray]
            MCMCサンプルの辞書
        """
        self._preprocessor = preprocessor
        
        # 説明変数の取得
        self._W_grids, _ = preprocessor.create_explanatory_variables(self.variable_names)
        
        # 有効グリッドの取得
        valid_grids = (
            preprocessor.df_elevation.sort(["y", "x"])
            .select(["is_valid"])
            .to_numpy()
            .astype(bool)
        )
        
        # 面積を計算
        volume = np.prod(
            [self.region[i][1] - self.region[i][0] for i in range(len(self.region))]
        ) * np.mean(valid_grids)
        
        # 観測データの準備
        n = preprocessor.df_sites.shape[0]
        X_obs = preprocessor.df_sites.select(["経度", "緯度"]).to_numpy().astype(np.float64)
        y_obs = np.ones(n)
        
        # 事前分布の設定
        p = len(self.variable_names) + 1
        prior_beta_mean = np.zeros(p)
        prior_beta_cov = np.eye(p) * 10
        prior_lambda_shape = 2.0
        prior_lambda_rate = 1.0
        
        # 初期値の設定
        beta_init = np.zeros(p)
        lambda_star_init = 30
        
        self._intensity_func = IntensityFunction(
            self.create_design_matrix, beta_init, lambda_star_init
        )
        
        # MCMC の実行
        beta_samples, lambda_star_samples = self.mcmc_sampler(
            intensity_func=self._intensity_func,
            y_obs=y_obs,
            X_obs=X_obs,
            valid_grids=valid_grids,
            volume=volume,
            num_iterations=num_iterations,
            prior_beta_mean=prior_beta_mean,
            prior_beta_cov=prior_beta_cov,
            prior_lambda_shape=prior_lambda_shape,
            prior_lambda_rate=prior_lambda_rate,
        )
        
        # バーンイン
        self._beta_samples = np.array(beta_samples)[burn_in:, :]
        self._lambda_star_samples = np.array(lambda_star_samples)[burn_in:]
        
        # 平均値で強度関数を更新
        self._intensity_func.update_beta(self._beta_samples.mean(axis=0))
        self._intensity_func.update_lambda_star(self._lambda_star_samples.mean())
        
        return {
            "beta_samples": self._beta_samples,
            "lambda_star_samples": self._lambda_star_samples,
        }

    def predict_site_probability(self, grid_coords: np.ndarray) -> np.ndarray:
        """
        グリッド座標での遺跡存在確率を予測

        Parameters
        ----------
        grid_coords : np.ndarray
            グリッド座標

        Returns
        -------
        np.ndarray
            存在確率
        """
        if self._intensity_func is None:
            raise ValueError("モデルが学習されていません。fit()を先に実行してください。")
        
        # grid_coords は (N, 2) でラジアン単位なので、度に変換
        X_grids = grid_coords / np.pi * 180
        X_grids = X_grids[:, [1, 0]]  # (lat, lon) -> (lon, lat)
        
        return self._intensity_func.q(X_grids)

    def create_inference_data(self) -> az.InferenceData:
        """
        ArviZ用のInferenceDataオブジェクトを作成

        Returns
        -------
        az.InferenceData
            推論データ
        """
        if self._beta_samples is None or self._lambda_star_samples is None:
            raise ValueError("モデルが学習されていません。fit()を先に実行してください。")
        
        # サンプル数とパラメータ数の取得
        num_samples, num_beta = self._beta_samples.shape
        
        # パラメータ名とサンプルの辞書を作成
        posterior_samples = {}
        
        # βの各成分を辞書に追加
        for i in range(num_beta):
            param_name = f"beta_{i}"
            posterior_samples[param_name] = self._beta_samples[:, i]
        
        # λ*を辞書に追加
        posterior_samples["lambda_star"] = self._lambda_star_samples
        
        # InferenceDataオブジェクトを作成
        idata = az.from_dict(posterior=posterior_samples)
        
        return idata

    @property
    def beta_samples(self) -> np.ndarray:
        """βのMCMCサンプル"""
        if self._beta_samples is None:
            raise ValueError("モデルが学習されていません。fit()を先に実行してください。")
        return self._beta_samples

    @property
    def lambda_star_samples(self) -> np.ndarray:
        """λ*のMCMCサンプル"""
        if self._lambda_star_samples is None:
            raise ValueError("モデルが学習されていません。fit()を先に実行してください。")
        return self._lambda_star_samples

    @property
    def intensity_func(self) -> IntensityFunction:
        """強度関数"""
        if self._intensity_func is None:
            raise ValueError("モデルが学習されていません。fit()を先に実行してください。")
        return self._intensity_func