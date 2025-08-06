"""
KSBP (Kernel Stick-Breaking Process) モデルの実装
ノートブック20_model_KSBP.ipynbからの移植
"""

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import polars as pl
from numpy.random import default_rng
from scipy.stats import beta, dirichlet
from tqdm import tqdm

from ..base.base_model import BaseCompositionModel
from ..preprocessing.data_preprocessor import ObsidianDataPreprocessor


class KSBPModel(BaseCompositionModel):
    """KSBP (Kernel Stick-Breaking Process) モデル"""

    def __init__(
        self,
        variable_names: List[str],
        kappa_coords: float = 0.002,
        kappa_costs: float = 2.0,
        kappa_elevation: float = 10000.0,
        kappa_angle: float = 3.0,
        kappa_river: float = 2.0,
        gamma: float = 0.1,
        alpha0: float = 0.1,
        j_max: int = 80,
        n_iter: int = 2000,
        burnin: int = 100,
        seed: int = 0,
    ):
        """
        Parameters
        ----------
        variable_names : List[str]
            使用する説明変数名のリスト
        kappa_coords : float
            座標用の長さ尺度パラメータ
        kappa_costs : float
            コスト系変数用の長さ尺度パラメータ
        kappa_elevation : float
            標高用の長さ尺度パラメータ
        kappa_angle : float
            角度用の長さ尺度パラメータ
        kappa_river : float
            河川用の長さ尺度パラメータ
        gamma : float
            ベータ分布のパラメータ
        alpha0 : float
            ディリクレ分布のパラメータ
        j_max : int
            最大クラスタ数
        n_iter : int
            MCMC反復数
        burnin : int
            バーンイン期間
        seed : int
            乱数シード
        """
        self.variable_names = variable_names
        self.kappa_coords = kappa_coords
        self.kappa_costs = kappa_costs
        self.kappa_elevation = kappa_elevation
        self.kappa_angle = kappa_angle
        self.kappa_river = kappa_river
        self.gamma = gamma
        self.alpha0 = alpha0
        self.j_max = j_max
        self.n_iter = n_iter
        self.burnin = burnin
        self.seed = seed

        self._fit_results: Optional[Dict[str, Any]] = None
        self._mu_z: Optional[np.ndarray] = None
        self._sd_z: Optional[np.ndarray] = None

    def _build_design_matrices(
        self,
        site_coords: np.ndarray,
        W_sites: np.ndarray,
        grid_coords: np.ndarray,
        W_grids: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """設計行列を構築"""
        mu_z = W_sites.mean(0)
        sd_z = W_sites.std(0, ddof=0)
        W_sites_z = (W_sites - mu_z) / sd_z
        W_grids_z = (W_grids - mu_z) / sd_z
        return (
            np.hstack([site_coords, W_sites_z]),
            np.hstack([grid_coords, W_grids_z]),
            mu_z,
            sd_z,
        )

    def _ard_kernel(
        self, X1: np.ndarray, X2: np.ndarray, ell: np.ndarray
    ) -> np.ndarray:
        """ARD Gaussianカーネル"""
        diff2 = (X1[:, None, :] - X2[None, :, :]) ** 2 / (ell**2)
        return np.exp(-0.5 * diff2.sum(-1))

    def _log_multinom(
        self, counts: np.ndarray, probs: np.ndarray, eps: float = 1e-20
    ) -> np.ndarray:
        """多項分布の対数尤度"""
        return np.sum(counts * np.log(probs + eps), axis=-1)

    def _compute_pi(self, psi: np.ndarray) -> np.ndarray:
        """スティックブレーキング過程からπを計算"""
        psi = np.clip(psi, 1e-300, 1 - 1e-16)
        cum = np.cumprod(1 - psi, axis=1)
        pi = psi.copy()
        if pi.shape[1] > 1:
            pi[:, 1:] *= cum[:, :-1]
        return pi

    def _create_length_scales(self) -> np.ndarray:
        """変数名に基づいて長さ尺度を作成"""
        ell = np.array([self.kappa_coords, self.kappa_coords])  # lat, lon

        for var_name in self.variable_names:
            if "elevation" in var_name:
                ell = np.append(ell, self.kappa_elevation)
            elif "angle" in var_name:
                ell = np.append(ell, self.kappa_angle)
            elif "cost" in var_name:
                ell = np.append(ell, self.kappa_costs)
            elif "river" in var_name:
                ell = np.append(ell, self.kappa_river)
            else:
                ell = np.append(ell, 1.0)

        return ell

    def _train_ddp(
        self,
        X_sites: np.ndarray,
        Y: np.ndarray,
        X_grids: np.ndarray,
        ell: np.ndarray,
    ) -> Dict[str, Any]:
        """DDPトレーニング"""
        rng = default_rng(self.seed)
        Nsites, dim = X_sites.shape
        p_cov = dim - 2
        K = Y.shape[1]
        N = Y.sum(1)

        if X_grids is None:
            X_grids = X_sites
        mins, maxs = X_grids.min(0), X_grids.max(0)

        # 初期φをグリッドからランダム抽出
        grid_coords = X_grids[:, :2]
        idx_sp = rng.choice(len(grid_coords), self.j_max, replace=True)
        phi = np.empty((self.j_max, dim))
        phi[:, :2] = grid_coords[idx_sp]

        cov_z = X_grids[:, 2:]
        low = np.nanmin(cov_z, 0)
        high = np.nanmax(cov_z, 0)
        phi[:, 2:] = rng.uniform(low, high, size=(self.j_max, p_cov))

        v = beta.rvs(1, self.gamma, size=self.j_max, random_state=rng)
        theta = dirichlet.rvs([self.alpha0 / K] * K, size=self.j_max, random_state=rng)
        z = rng.integers(self.j_max, size=Nsites)
        u = rng.uniform(0, 0.5, size=Nsites)

        n_store = (self.n_iter - self.burnin) // 1
        phi_tr = np.empty((n_store, self.j_max, dim))
        v_tr = np.empty((n_store, self.j_max))
        theta_tr = np.empty((n_store, self.j_max, K))
        trace_k = []
        s_ix = 0

        for it in tqdm(range(self.n_iter), desc="KSBP MCMC"):
            Kmat = self._ard_kernel(X_sites, phi, ell)
            psi = v * Kmat
            pi = self._compute_pi(psi)

            # Z Gibbs
            lp = np.log(pi) + self._log_multinom(Y[:, None, :], theta[None, :, :])
            pz = np.exp(lp - lp.max(1, keepdims=True))
            pz /= pz.sum(1, keepdims=True)
            z = np.array([rng.choice(self.j_max, p=pz[i]) for i in range(Nsites)])

            # slice
            u = rng.uniform(0, pi[np.arange(Nsites), z])
            L = int((np.cumsum(psi, 1) < (1 - u.min())).any(0).sum()) + 1
            L = min(max(L, 3), self.j_max)
            active = slice(0, L)

            # θ Dirichlet
            for h in range(L):
                y_sum = Y[z == h].sum(0)
                theta[h] = rng.dirichlet(
                    self.alpha0 / K + y_sum if y_sum.any() else [self.alpha0 / K] * K
                )

            # v 共役 Beta
            m = np.array([(z == h).sum() for h in range(L)])
            bmask = u[:, None] < psi[:, active]
            r = np.array([(bmask[:, h] & (z > h)).sum() for h in range(L)])
            v[active] = beta.rvs(1 + m, self.gamma + r, random_state=rng)

            # φ MH
            for h in range(L):
                prop = phi[h] + rng.normal(scale=0.05, size=dim)
                if ((prop < mins).any()) or ((prop > maxs).any()):
                    continue
                psi_new = psi.copy()
                psi_new[:, h] = (
                    v[h] * self._ard_kernel(X_sites, prop[None, :], ell)[:, 0]
                )
                pi_new = self._compute_pi(psi_new)
                idx = z >= h
                if np.log(rng.uniform()) < np.sum(
                    np.log(pi_new[np.arange(Nsites)[idx], z[idx]] + 1e-20)
                    - np.log(pi[np.arange(Nsites)[idx], z[idx]] + 1e-20)
                ):
                    phi[h] = prop
                    psi = psi_new
                    pi = pi_new

            # v MH
            eps = 1e-12
            for h in range(L):
                v_h = np.clip(v[h], eps, 1.0 - eps)
                logit = np.log(v_h) - np.log1p(-v_h)
                logit_p = logit + rng.normal(scale=0.3)
                v_p = 1 / (1 + np.exp(-logit_p))
                v_p = np.clip(v_p, eps, 1.0 - eps)

                J_old = np.log(v_h) + np.log1p(-v_h)
                J_new = np.log(v_p) + np.log1p(-v_p)

                log_prior = (
                    beta.logpdf(v_p, 1, self.gamma)
                    - beta.logpdf(v_h, 1, self.gamma)
                    + (J_new - J_old)
                )

                psi_new = psi.copy()
                psi_new[:, h] = v_p * Kmat[:, h]
                pi_new = self._compute_pi(psi_new)
                idx = z >= h
                log_like = (
                    np.log(pi_new[np.arange(Nsites)[idx], z[idx]])
                    - np.log(pi[np.arange(Nsites)[idx], z[idx]])
                ).sum()

                log_alpha = log_prior + log_like
                if np.log(rng.uniform()) < log_alpha:
                    v[h] = v_p
                    psi, pi = psi_new, pi_new

            trace_k.append(np.unique(z).size)
            if it >= self.burnin and (it - self.burnin) % 1 == 0:
                phi_tr[s_ix], v_tr[s_ix], theta_tr[s_ix] = phi, v, theta
                s_ix += 1

        return dict(
            phi_trace=phi_tr,
            v_trace=v_tr,
            theta_trace=theta_tr,
            phi=phi,
            v=v,
            theta=theta,
            ell=ell,
            trace_clusters=np.array(trace_k),
        )

    def fit(self, preprocessor: ObsidianDataPreprocessor) -> Dict[str, Any]:
        """
        KSBPモデルを学習

        Parameters
        ----------
        preprocessor : ObsidianDataPreprocessor
            前処理済みのデータ

        Returns
        -------
        Dict[str, Any]
            学習結果
        """
        print("=== KSBP学習を開始 ===")

        # データの準備
        data = preprocessor.create_X_Y(
            self.variable_names,
            {0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"},
            ["神津島", "信州", "箱根", "高原山", "その他"],
        )

        site_coords = preprocessor.create_site_coords()
        grid_coords = preprocessor.create_grid_coords()
        W_grids, W_sites = preprocessor.create_explanatory_variables(
            self.variable_names
        )

        # NaN値を0で埋める
        W_grids = np.nan_to_num(W_grids, nan=0.0)

        # 設計行列の構築
        X_sites, X_grids, mu_z, sd_z = self._build_design_matrices(
            site_coords, W_sites, grid_coords, W_grids
        )

        self._mu_z = mu_z
        self._sd_z = sd_z

        # 長さ尺度の設定
        ell = self._create_length_scales()

        print(f"設計行列のサイズ: X_sites={X_sites.shape}, X_grids={X_grids.shape}")
        print(f"長さ尺度: {ell}")

        # すべての時期を統合したデータで学習
        Y_all = data["Y"]  # (遺跡数, 時期数×産地数)

        # 時期×産地の形に変換（最初の時期のみ使用）
        n_origins = 4  # "その他"を除く
        Y_period0 = Y_all[:, :n_origins]  # 最初の時期のデータ

        # 学習実行
        results = self._train_ddp(X_sites, Y_period0, X_grids, ell)
        results["X_sites"] = X_sites
        results["X_grids"] = X_grids
        results["mu_z"] = mu_z
        results["sd_z"] = sd_z

        self._fit_results = results

        print("=== KSBP学習完了 ===")

        return results

    def _pi_from_psi_batch(self, psi: np.ndarray, eps: float = 1e-12) -> np.ndarray:
        """バッチ版π計算"""
        psi = np.clip(psi, eps, 1.0 - eps)
        log_psi = np.log(psi)
        log1m = np.log1p(-psi)
        log_cum = np.cumsum(log1m, axis=-1)
        log_cum = np.concatenate(
            [np.zeros_like(log_cum[..., :1]), log_cum[..., :-1]], axis=-1
        )
        return np.exp(log_psi + log_cum)

    def predict(
        self,
        X_new: np.ndarray,
        aggregate: str = "mean",
        batch_size: int = 1,
        S_reduced: int = 20,
        return_sd: bool = False,
        eps: float = 1e-20,
    ) -> np.ndarray:
        """予測実行"""
        if self._fit_results is None:
            raise ValueError(
                "モデルが学習されていません。fit()を先に実行してください。"
            )

        phi_tr, v_tr, theta_tr = (
            self._fit_results["phi_trace"],
            self._fit_results["v_trace"],
            self._fit_results["theta_trace"],
        )
        ell = self._fit_results["ell"]
        dim = X_new.shape[1]

        # サンプル間引き
        if phi_tr.shape[0] > S_reduced:
            idx = np.random.choice(phi_tr.shape[0], S_reduced, replace=False)
            phi_tr, v_tr, theta_tr = phi_tr[idx], v_tr[idx], theta_tr[idx]

        S, J, _ = phi_tr.shape
        m, K = X_new.shape[0], theta_tr.shape[2]

        if aggregate == "mean":
            sum1 = np.zeros((m, K))
            sum2 = np.zeros((m, K))
        else:
            raise NotImplementedError(f"aggregate='{aggregate}' は未実装です")

        for s0 in tqdm(range(0, S, batch_size), desc="KSBP予測"):
            s1 = min(S, s0 + batch_size)
            phi_b, v_b, th_b = phi_tr[s0:s1], v_tr[s0:s1], theta_tr[s0:s1]

            K_b = self._ard_kernel(X_new, phi_b.reshape(-1, dim), ell)
            K_b = K_b.reshape(s1 - s0, m, J)
            psi_b = v_b[:, None, :] * K_b
            pi_b = self._pi_from_psi_batch(psi_b)
            probs_b = np.einsum("bmj,bjk->bmk", pi_b, th_b)
            denom = probs_b.sum(2, keepdims=True)
            probs_b /= np.where(denom < eps, 1.0, denom)

            sum1 += probs_b.sum(0)
            sum2 += (probs_b**2).sum(0)

        mean = sum1 / S
        if return_sd:
            var = (sum2 - (sum1**2) / S) / (S - 1 + eps)
            sd = np.sqrt(np.maximum(var, 0.0))
            return mean, sd
        return mean

    def predict_site_ratios(
        self, preprocessor: ObsidianDataPreprocessor
    ) -> Dict[str, pl.DataFrame]:
        """遺跡の産地構成比を予測（BaseCompositionModel準拠・標準化済み）"""
        if self._fit_results is None:
            raise ValueError("モデルが学習されていません。")

        X_sites = self._fit_results["X_sites"]
        pi_sites = self.predict(X_sites)

        # 遺跡IDを取得
        unique_sites = preprocessor.df_obsidian.unique(subset=["遺跡ID"]).sort("遺跡ID")

        # 結果をDataFrameに変換（主要4産地のみ）
        ratio_sites_df = pl.DataFrame({"遺跡ID": unique_sites["遺跡ID"]})

        origins = ["神津島", "信州", "箱根", "高原山"]
        for i, origin in enumerate(origins):
            ratio_sites_df = ratio_sites_df.with_columns(
                pl.Series(f"ratio_0_{origin}", pi_sites[:, i])
            )

        # 現在は時期0のデータのみ対応（将来拡張可能）
        return {"0": ratio_sites_df}

    def predict_grid_ratios(
        self, preprocessor: ObsidianDataPreprocessor
    ) -> pl.DataFrame:
        """グリッドの産地構成比を予測"""
        if self._fit_results is None:
            raise ValueError("モデルが学習されていません。")

        X_grids = self._fit_results["X_grids"]
        pi_grids = self.predict(X_grids, S_reduced=5)

        # グリッド座標を取得
        lon_mesh, lat_mesh = preprocessor.create_meshgrid()

        # 結果をDataFrameに変換
        ratio_df = pl.DataFrame({"x": lon_mesh.ravel(), "y": lat_mesh.ravel()})

        origins = ["神津島", "信州", "箱根", "高原山"]
        for i, origin in enumerate(origins):
            ratio_df = ratio_df.with_columns(
                pl.Series(f"ratio_0_{origin}", pi_grids[:, i])
            )

        return ratio_df

    @property
    def model_name(self) -> str:
        """モデル名"""
        return "KSBP"

    @property
    def results(self) -> Dict[str, Any]:
        """学習結果"""
        if self._fit_results is None:
            raise ValueError("モデルが学習されていません。")
        return self._fit_results
