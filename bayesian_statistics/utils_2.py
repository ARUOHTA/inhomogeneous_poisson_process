from typing import Callable, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from pypolyagamma import PyPolyaGamma
from tqdm import tqdm


def multi_pgdraw_vectorized(
    pg: PyPolyaGamma, B: np.ndarray, C: np.ndarray
) -> np.ndarray:
    """multi_pgdraw のベクトル化バージョン"""
    return np.array([pg.pgdraw(b, c) for b, c in zip(B, C)])


class IntensityFunction:
    def __init__(
        self,
        design_matrix_func: Callable[[np.ndarray], np.ndarray],
        beta: np.ndarray = None,
        lambda_star: float = None,
        link_function: Callable[[np.ndarray], np.ndarray] = None,
    ):
        self.design_matrix_func = design_matrix_func
        self.beta = beta
        self.lambda_star = lambda_star
        if link_function is None:
            # デフォルトでロジスティックリンク関数を使用
            self.link_function = lambda eta: 1 / (1 + np.exp(-eta))
        else:
            self.link_function = link_function

    def copy(self):
        return IntensityFunction(
            design_matrix_func=self.design_matrix_func,
            beta=self.beta.copy(),
            lambda_star=self.lambda_star,
            link_function=self.link_function,
        )

    def q(self, x: np.ndarray) -> np.ndarray:
        """q(x) を計算"""
        X = self.design_matrix_func(x)
        eta = X @ self.beta
        return self.link_function(eta)

    def lambda_func(self, x: np.ndarray) -> np.ndarray:
        """λ(x) = q(x) * λ* を計算"""
        return self.q(x) * self.lambda_star

    def update_beta(self, beta: np.ndarray):
        self.beta = beta

    def update_lambda_star(self, lambda_star: float):
        self.lambda_star = lambda_star


def create_design_matrix(x: np.ndarray) -> np.ndarray:
    """設計行列を生成する関数。x の形状は (N, D)"""
    if x.ndim == 1:
        # 1次元の場合
        x = x[:, np.newaxis]
    # x の形状は (N, D)
    # 例として、多項式特徴量と相互作用項を作成
    # x1, x2 の場合
    poly_terms = [np.ones(len(x))]  # バイアス項
    for i in range(x.shape[1]):
        poly_terms.append(x[:, i])  # 一次項
        poly_terms.append(x[:, i] ** 2)  # 二次項
    # 相互作用項（2次元以上の場合）
    if x.shape[1] > 1:
        poly_terms.append(x[:, 0] * x[:, 1])
    return np.column_stack(poly_terms)


def generate_events(
    region: Tuple[np.ndarray, np.ndarray], intensity_func: IntensityFunction
) -> np.ndarray:
    """非斉次ポアソン過程によるイベントの生成（任意の次元）"""
    events = []

    # 領域の計算
    mins = np.array([region[i][0] for i in range(len(region))])
    maxs = np.array([region[i][1] for i in range(len(region))])
    volume = np.prod(maxs - mins)

    # λ(x) の最大値を計算
    grid_size = 100
    grids = [np.linspace(mins[i], maxs[i], grid_size) for i in range(len(region))]
    mesh = np.meshgrid(*grids)
    grid_points = np.column_stack([m.flatten() for m in mesh])
    lambdas = intensity_func.lambda_func(grid_points)
    max_lambda = np.max(lambdas)

    # 全体のポアソン過程から候補点を生成
    lambda_total = max_lambda * volume
    N = np.random.poisson(lambda_total)
    candidate_points = np.random.uniform(mins, maxs, size=(N, len(region)))

    # Poisson thinning
    lambda_candidates = intensity_func.lambda_func(candidate_points)
    acceptance_probs = lambda_candidates / max_lambda
    u = np.random.uniform(size=N)
    accepted = u < acceptance_probs
    events = candidate_points[accepted]

    return events


def generate_U(
    region: Tuple[np.ndarray, np.ndarray], intensity_func: IntensityFunction
) -> np.ndarray:
    """非斉次ポアソン過程による潜在変数 U の生成（任意の次元）"""
    U = []

    # 領域の計算
    mins = np.array([region[i][0] for i in range(len(region))])
    maxs = np.array([region[i][1] for i in range(len(region))])
    volume = np.prod(maxs - mins)

    # (1 - q(x)) * λ* の最大値を計算
    grid_size = 100
    grids = [np.linspace(mins[i], maxs[i], grid_size) for i in range(len(region))]
    mesh = np.meshgrid(*grids)
    grid_points = np.column_stack([m.flatten() for m in mesh])
    lambdas = (1 - intensity_func.q(grid_points)) * intensity_func.lambda_star
    max_lambda = np.max(lambdas)

    # 全体のポアソン過程から候補点を生成
    lambda_total = max_lambda * volume
    N = np.random.poisson(lambda_total)
    candidate_points = np.random.uniform(mins, maxs, size=(N, len(region)))

    # Poisson thinning
    q_candidates = intensity_func.q(candidate_points)
    lambda_candidates = (1 - q_candidates) * intensity_func.lambda_star
    acceptance_probs = lambda_candidates / max_lambda
    u = np.random.uniform(size=N)
    accepted = u < acceptance_probs
    U = candidate_points[accepted]

    return U


def poisson_thinning_U(
    region: Tuple[np.ndarray, np.ndarray], intensity_func: IntensityFunction
) -> Tuple[np.ndarray, np.ndarray]:
    """Poisson thinning を用いて U をサンプリング（任意の次元）"""
    # 領域の計算
    mins = np.array([region[i][0] for i in range(len(region))])
    maxs = np.array([region[i][1] for i in range(len(region))])
    volume = np.prod(maxs - mins)

    # λ* × 領域の体積
    lambda_total = intensity_func.lambda_star * volume
    N = np.random.poisson(lambda_total)
    candidate_points = np.random.uniform(mins, maxs, size=(N, len(region)))

    # q(x) を計算
    q_candidates = intensity_func.q(candidate_points)
    u = np.random.uniform(size=N)
    # (1 - q(x)) でフィルタリング
    U = candidate_points[u < (1 - q_candidates)]

    return U, candidate_points


def generate_IPP(
    region: Tuple[np.ndarray, np.ndarray],
    intensity_func: IntensityFunction,
    generate_U_flag: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """非斉次ポアソン過程モデルに基づいてデータを生成する（任意の次元）"""
    # 観測データ X の生成
    X_events = generate_events(region, intensity_func)

    # 潜在変数 U の生成
    U_events = np.array([])
    if generate_U_flag:
        U_events, _ = poisson_thinning_U(region, intensity_func)

    # データの結合
    X = np.concatenate([X_events, U_events]) if U_events.size > 0 else X_events
    n_X = len(X_events)
    n_U = len(U_events)
    n = n_X + n_U

    # 設計行列の作成
    W = intensity_func.design_matrix_func(X)

    # 応答変数 y の作成
    y = (
        np.concatenate([np.ones(n_X), np.zeros(n_U)]).astype(int)
        if n_U > 0
        else np.ones(n_X)
    )

    return W, y, X


def mcmc_sampler(
    region: Tuple[np.ndarray, np.ndarray],
    intensity_func: IntensityFunction,
    y_obs: np.ndarray,
    X_obs: np.ndarray,
    num_iterations: int,
    prior_beta_mean: np.ndarray,
    prior_beta_cov: np.ndarray,
    prior_lambda_shape: float,
    prior_lambda_rate: float,
):
    """MCMC サンプラー"""
    # 初期化
    beta_samples = []
    lambda_star_samples = []
    pg = PyPolyaGamma()
    beta = intensity_func.beta
    lambda_star = intensity_func.lambda_star

    for iteration in tqdm(range(num_iterations)):
        # 1. U のサンプリング
        U_events, candidate_points = poisson_thinning_U(region, intensity_func)
        n_U = len(U_events)

        # 2. 結合データの作成
        X_combined = np.concatenate([X_obs, U_events]) if n_U > 0 else X_obs
        W_combined = intensity_func.design_matrix_func(X_combined)
        y_combined = np.concatenate([y_obs, np.zeros(n_U)]) if n_U > 0 else y_obs
        n_combined = len(y_combined)

        # 3. ω のサンプリング
        psi = W_combined @ beta
        omega = multi_pgdraw_vectorized(pg, np.ones(n_combined), psi)

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
            [region[i][1] - region[i][0] for i in range(len(region))]
        )
        lambda_star = np.random.gamma(shape=lambda_shape, scale=1 / lambda_rate)
        intensity_func.update_lambda_star(lambda_star)

        # サンプルの保存
        beta_samples.append(beta)
        lambda_star_samples.append(lambda_star)

    return beta_samples, lambda_star_samples


def plot_intensity(T: float, intensity_func: IntensityFunction):
    """強度関数のプロット"""
    x_lin = np.linspace(0, T, 1000)
    intensity = intensity_func.lambda_func(x_lin)
    plt.figure(figsize=(8, 6))
    plt.plot(x_lin, intensity, c="black", label="Intensity λ(t)")
    plt.xlabel("Time")
    plt.ylabel("Intensity")
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.legend(loc="upper right")
    plt.show()


def plot_events(
    t0: np.ndarray, t1: np.ndarray, T: float, intensity_func: IntensityFunction
):
    """イベントと背景データのプロット"""
    plt.figure(figsize=(8, 6))
    x_lin = np.linspace(0, T, 1000)
    intensity = intensity_func.lambda_func(x_lin)
    if len(t0) > 0:
        plt.scatter(
            t0,
            np.zeros(len(t0)),
            c="grey",
            alpha=0.4,
            label=f"Pseudo-absence U: {len(t0)}",
        )
    plt.scatter(
        t1, np.zeros(len(t1)), c="green", alpha=0.3, label=f"Presence X: {len(t1)}"
    )
    plt.plot(x_lin, intensity, c="black", label="Intensity λ(t)")
    plt.xlabel("Time")
    plt.ylabel("Intensity")
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.legend(loc="upper right")
    plt.show()


def plot_intensity_2d(
    region: Tuple[np.ndarray, np.ndarray],
    intensity_func: IntensityFunction,
    num_points: int = 100,
):
    """2次元の強度関数をプロット"""
    x_min, x_max = region[0]
    y_min, y_max = region[1]
    x = np.linspace(x_min, x_max, num_points)
    y = np.linspace(y_min, y_max, num_points)
    X, Y = np.meshgrid(x, y)
    grid_points = np.column_stack([X.ravel(), Y.ravel()])
    intensity = intensity_func.lambda_func(grid_points).reshape(num_points, num_points)

    plt.figure(figsize=(8, 6))
    plt.contourf(X, Y, intensity, levels=50, cmap="viridis")
    plt.colorbar(label="Intensity λ(x, y)")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Intensity Function")
    plt.show()


def plot_events_2d(
    region: Tuple[np.ndarray, np.ndarray],
    events: np.ndarray,
    intensity_func: IntensityFunction,
    num_points: int = 100,
    plot_intensity: bool = True,
):
    """2次元のイベントと強度関数をプロット"""
    if plot_intensity:
        x_min, x_max = region[0]
        y_min, y_max = region[1]
        x = np.linspace(x_min, x_max, num_points)
        y = np.linspace(y_min, y_max, num_points)
        X, Y = np.meshgrid(x, y)
        grid_points = np.column_stack([X.ravel(), Y.ravel()])
        intensity = intensity_func.lambda_func(grid_points).reshape(
            num_points, num_points
        )

    if plot_intensity:
        plt.figure(figsize=(8, 6))
        plt.contourf(X, Y, intensity, levels=50, cmap="viridis", alpha=0.6)
        plt.colorbar(label="Intensity λ(x, y)")
    else:
        plt.figure(figsize=(6.45, 6))
    plt.scatter(
        events[:, 0],
        events[:, 1],
        c="red",
        s=10,
        label=f"Observed Events (N={len(events)})",
        marker="x",
    )
    plt.xlabel("X")
    plt.ylabel("Y")
    if plot_intensity:
        plt.title("Events and Intensity Function")
    else:
        plt.title("Events")
    plt.legend()
    plt.show()


def plot_intensity_posterior_2d(
    region: Tuple[np.ndarray, np.ndarray],
    beta_samples: np.ndarray,
    lambda_star_samples: np.ndarray,
    intensity_func: IntensityFunction,
    num_points: int = 100,
    credible_interval: float = 0.95,
    true_intensity_func: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    events: Optional[np.ndarray] = None,
    cmap: str = "viridis",
):
    x_min, x_max = region[0]
    y_min, y_max = region[1]
    x = np.linspace(x_min, x_max, num_points)
    y = np.linspace(y_min, y_max, num_points)
    X_grid, Y_grid = np.meshgrid(x, y)
    grid_points = np.column_stack([X_grid.ravel(), Y_grid.ravel()])
    num_grid_points = grid_points.shape[0]
    num_samples = beta_samples.shape[0]

    # 強度関数のサンプルを格納する配列
    lambda_samples = np.zeros((num_samples, num_grid_points))

    # 各サンプルについて強度関数を計算
    for i in range(num_samples):
        beta = beta_samples[i]
        lambda_star = lambda_star_samples[i]
        # IntensityFunctionを更新
        intensity_func.update_beta(beta)
        intensity_func.update_lambda_star(lambda_star)
        # 強度関数を計算
        lambda_samples[i, :] = intensity_func.lambda_func(grid_points)

    # 統計量の計算
    lambda_mean = np.mean(lambda_samples, axis=0).reshape(num_points, num_points)

    # プロット
    plt.figure(figsize=(12, 8))
    plt.contourf(X_grid, Y_grid, lambda_mean, levels=50, cmap=cmap, alpha=0.6)
    plt.colorbar(label="Posterior Mean of λ(x, y)")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Posterior Mean of Intensity Function")

    # イベントのプロット（オプション）
    if events is not None:
        plt.scatter(
            events[:, 0],
            events[:, 1],
            c="red",
            marker="x",
            s=10,
            label="Observed Events",
        )
        plt.legend()

    plt.show()
