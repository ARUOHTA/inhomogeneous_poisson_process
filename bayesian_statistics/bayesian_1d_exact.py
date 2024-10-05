import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Callable, Tuple, Optional
from pypolyagamma import PyPolyaGamma
from tqdm import tqdm

class IntensityFunction:
    def __init__(
        self, 
        design_matrix_func: Callable[[np.ndarray], np.ndarray], 
        beta: np.ndarray,
        lambda_star: float,
        link_function: Callable[[np.ndarray], np.ndarray] = None
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
            link_function=self.link_function
        )

    def q(self, x: np.ndarray) -> np.ndarray:
        """q(t) を計算"""
        X = self.design_matrix_func(x)
        eta = X @ self.beta
        return self.link_function(eta)

    def lambda_func(self, x: np.ndarray) -> np.ndarray:
        """λ(t) = q(t) * λ* を計算"""
        return self.q(x) * self.lambda_star

    def update_beta(self, beta: np.ndarray):
        self.beta = beta

    def update_lambda_star(self, lambda_star: float):
        self.lambda_star = lambda_star

def create_design_matrix(x: np.ndarray) -> np.ndarray:
    """設計行列を生成する関数"""
    return np.vstack([np.ones_like(x), x, x**2]).T

def generate_events(T: float, intensity_func: IntensityFunction) -> np.ndarray:
    """非斉次ポアソン過程 IPP(q(t) * λ*) によるイベントの生成"""
    events = []
    t = 0

    # λ(t) の最大値を計算
    x_lin = np.linspace(0, T, 1000)
    lambdas = intensity_func.lambda_func(x_lin)
    max_lambda = np.max(lambdas)
    
    while t < T:
        u = np.random.uniform()
        t += -np.log(u) / max_lambda
        if t < T:
            lambda_t = intensity_func.lambda_func(np.array([t]))[0]
            if np.random.uniform() < lambda_t / max_lambda:
                events.append(t)

    return np.array(events)

def generate_U(T: float, intensity_func: IntensityFunction) -> np.ndarray:
    """非斉次ポアソン過程 IPP((1 - q(t)) * λ*) による潜在変数 U の生成"""
    U = []
    t = 0

    # (1 - q(t)) * λ* の最大値を計算
    x_lin = np.linspace(0, T, 1000)
    lambdas = (1 - intensity_func.q(x_lin)) * intensity_func.lambda_star
    max_lambda = np.max(lambdas)
    
    while t < T:
        u = np.random.uniform()
        t += -np.log(u) / max_lambda
        if t < T:
            q_t = intensity_func.q(np.array([t]))[0]
            lambda_t = (1 - q_t) * intensity_func.lambda_star
            if np.random.uniform() < lambda_t / max_lambda:
                U.append(t)

    return np.array(U)

def poisson_thinning_U(T: float, intensity_func: IntensityFunction) -> Tuple[np.ndarray, np.ndarray]:
    """Poisson thinning を用いて U をサンプリング"""
    # 全体の Poisson 点過程から候補点を生成
    lambda_total = intensity_func.lambda_star * T  # λ* × 領域の長さ
    N = np.random.poisson(lambda_total)
    t_candidates = np.random.uniform(0, T, N)

    # q(t) を計算
    q_t = intensity_func.q(t_candidates)
    u = np.random.uniform(size=N)
    # (1 - q(t)) でフィルタリング
    U = t_candidates[u < (1 - q_t)]

    return U, t_candidates

def plot_intensity(T: float, intensity_func: IntensityFunction):
    """強度関数のプロット"""
    x_lin = np.linspace(0, T, 1000)
    intensity = intensity_func.lambda_func(x_lin)
    plt.figure(figsize=(8, 6))
    plt.plot(x_lin, intensity, c='black', label='Intensity λ(t)')
    plt.xlabel('Time')
    plt.ylabel('Intensity')
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.legend(loc='upper right')
    plt.show()

def plot_events(t0: np.ndarray, t1: np.ndarray, T: float, intensity_func: IntensityFunction):
    """イベントと背景データのプロット"""
    plt.figure(figsize=(8, 6))
    x_lin = np.linspace(0, T, 1000)
    intensity = intensity_func.lambda_func(x_lin)
    if len(t0) > 0:
        plt.scatter(t0, np.zeros(len(t0)), c='grey', alpha=0.4, label=f"Pseudo-absence U: {len(t0)}")
    plt.scatter(t1, np.zeros(len(t1)), c='green', alpha=0.3, label=f"Presence X: {len(t1)}")
    plt.plot(x_lin, intensity, c='black', label='Intensity λ(t)')
    plt.xlabel('Time')
    plt.ylabel('Intensity')
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.legend(loc='upper right')
    plt.show()

def generate_IPP(
    T: float, 
    intensity_func: IntensityFunction,  
    generate_U_flag: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """1次元の非斉次ポアソン過程モデルに基づいてデータを生成する"""
    # 観測データ X の生成
    t_X = generate_events(T, intensity_func)

    # 潜在変数 U の生成
    t_U = np.array([])
    if generate_U_flag:
        t_U, _ = poisson_thinning_U(T, intensity_func)

    # 時間の結合
    t = np.concatenate([t_X, t_U])
    n_X = len(t_X)
    n_U = len(t_U)
    n = n_X + n_U

    # 設計行列の作成
    X = intensity_func.design_matrix_func(t)

    # 応答変数 y の作成
    y = np.concatenate([np.ones(n_X), np.zeros(n_U)]).astype(int)

    return X, y, t

def multi_pgdraw_vectorized(pg: PyPolyaGamma, B: np.ndarray, C: np.ndarray) -> np.ndarray:
    """multi_pgdraw のベクトル化バージョン"""
    return np.array([pg.pgdraw(b, c) for b, c in zip(B, C)])

# MCMC のサンプル実装
def mcmc_sampler(
    T: float,
    intensity_func: IntensityFunction,
    W_obs: np.ndarray,
    y_obs: np.ndarray,
    num_iterations: int,
    prior_beta_mean: np.ndarray,
    prior_beta_cov: np.ndarray,
    prior_lambda_shape: float,
    prior_lambda_rate: float
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
        t_U, t_candidates = poisson_thinning_U(T, intensity_func)
        n_U = len(t_U)

        # 2. 結合データの作成
        t_combined = np.concatenate([W_obs[:, 1], t_U])  # X_obs[:, 1] は時間 t を示すと仮定
        W_combined = intensity_func.design_matrix_func(t_combined)
        y_combined = np.concatenate([y_obs, np.zeros(n_U)])
        n_combined = len(y_combined)

        # 3. ω のサンプリング
        psi = W_combined @ beta
        omega = multi_pgdraw_vectorized(pg, np.ones(n_combined), psi)

        # 4. β のサンプリング
        Omega = np.diag(omega)
        z_tilde = (y_combined - 0.5) / omega
        V_inv = np.linalg.inv(prior_beta_cov) + W_combined.T @ Omega @ W_combined
        V = np.linalg.inv(V_inv)
        m = V @ (np.linalg.inv(prior_beta_cov) @ prior_beta_mean + W_combined.T @ Omega @ z_tilde)
        beta = np.random.multivariate_normal(m, V)
        intensity_func.update_beta(beta)

        # 5. λ* のサンプリング
        n_total = n_combined
        lambda_shape = prior_lambda_shape + n_total
        lambda_rate = prior_lambda_rate + T
        lambda_star = np.random.gamma(shape=lambda_shape, scale=1/lambda_rate)
        intensity_func.update_lambda_star(lambda_star)

        # サンプルの保存
        beta_samples.append(beta)
        lambda_star_samples.append(lambda_star)

    return beta_samples, lambda_star_samples




def plot_posterior_distributions(beta_samples, lambda_star_samples, true_beta=None, true_lambda_star=None):
    """
    パラメータの事後分布をプロットする関数。

    Parameters:
    - beta_samples: MCMC で得られた β のサンプル（二次元配列）
    - lambda_star_samples: MCMC で得られた λ* のサンプル（一次元配列）
    - true_beta: 真の β の値（既知の場合、プロットに赤線で表示されます）
    - true_lambda_star: 真の λ* の値（既知の場合、プロットに赤線で表示されます）
    """
    num_params = beta_samples.shape[1]
    
    plt.figure(figsize=(15, 4))
    
    # β の事後分布のプロット
    for i in range(num_params):
        plt.subplot(1, num_params + 1, i + 1)
        sns.kdeplot(beta_samples[:, i], shade=True, color='blue')
        plt.title(f'Posterior of β{i}')
        if true_beta is not None:
            plt.axvline(true_beta[i], color='red', linestyle='--', label='True value')
            plt.legend()
        plt.xlabel(f'β{i}')
        plt.ylabel('Density')
    
    # λ* の事後分布のプロット
    plt.subplot(1, num_params + 1, num_params + 1)
    sns.kdeplot(lambda_star_samples, shade=True, color='green')
    plt.title('Posterior of λ*')
    if true_lambda_star is not None:
        plt.axvline(true_lambda_star, color='red', linestyle='--', label='True value')
        plt.legend(loc='upper right')
    plt.xlabel('λ*')
    plt.ylabel('Density')
    
    plt.tight_layout()
    plt.show()

def plot_trace_plots(beta_samples, lambda_star_samples, true_beta=None, true_lambda_star=None):
    """
    パラメータのトレースプロットを作成する関数。
    """
    num_params = beta_samples.shape[1]
    
    plt.figure(figsize=(15, 4))
    
    # β のトレースプロット
    for i in range(num_params):
        plt.subplot(1, num_params + 1, i + 1)
        plt.plot(beta_samples[:, i], color='blue')
        if true_beta is not None:
            plt.axhline(true_beta[i], color='red', linestyle='--', label='True value')
            plt.legend(loc='upper right')
        plt.title(f'Trace of β{i}')
        plt.xlabel('Iteration')
        plt.ylabel(f'β{i}')
    
    # λ* のトレースプロット
    plt.subplot(1, num_params + 1, num_params + 1)
    plt.plot(lambda_star_samples, color='green')
    if true_lambda_star is not None:
        plt.axhline(true_lambda_star, color='red', linestyle='--', label='True value')
        plt.legend(loc='upper right')
    plt.title('Trace of λ*')
    plt.xlabel('Iteration')
    plt.ylabel('λ*')
    
    plt.tight_layout()
    plt.show()


def plot_intensity_posterior_with_events(
    beta_samples: np.ndarray,
    lambda_star_samples: np.ndarray,
    intensity_func: IntensityFunction,
    T: float,
    events: np.ndarray,
    num_points: int = 100,
    credible_interval: float = 0.95,
    true_intensity_func: Optional[Callable[[np.ndarray], np.ndarray]] = None
):
    """
    強度関数の事後推定値と信用区間、観測イベントをプロットする関数。

    Parameters:
    - beta_samples: MCMCで得られたβのサンプル（二次元配列）
    - lambda_star_samples: MCMCで得られたλ*のサンプル（一次元配列）
    - intensity_func: IntensityFunctionクラスのインスタンス（betaとlambda_starは任意）
    - T: 時間範囲の上限値
    - events: 観測されたイベントのタイムスタンプの配列
    - num_points: 時間範囲を分割するポイント数
    - credible_interval: 信用区間の幅（0から1の値、デフォルトは0.95）
    - true_intensity_func: 真の強度関数（Callable）、既知の場合はプロットに追加可能
    """

    # タイムポイントの設定
    t_values = np.linspace(0, T, num_points)
    num_samples = beta_samples.shape[0]

    # 強度関数のサンプルを格納する配列
    lambda_samples = np.zeros((num_samples, num_points))

    # 各サンプルについて強度関数を計算
    for i in range(num_samples):
        beta = beta_samples[i]
        lambda_star = lambda_star_samples[i]
        # IntensityFunctionを更新
        intensity_func.update_beta(beta)
        intensity_func.update_lambda_star(lambda_star)
        # 強度関数を計算
        lambda_samples[i, :] = intensity_func.lambda_func(t_values)

    # 統計量の計算
    lambda_mean = np.mean(lambda_samples, axis=0)
    lambda_median = np.median(lambda_samples, axis=0)
    lower_percentile = (1 - credible_interval) / 2 * 100
    upper_percentile = (1 + credible_interval) / 2 * 100
    lambda_lower = np.percentile(lambda_samples, lower_percentile, axis=0)
    lambda_upper = np.percentile(lambda_samples, upper_percentile, axis=0)

    # プロット
    plt.figure(figsize=(8, 6))

    # 強度関数の信用区間の塗りつぶし
    plt.fill_between(
        t_values, lambda_lower, lambda_upper,
        color='skyblue', alpha=0.5,
        label=f'{int(credible_interval*100)}% Credible Interval'
    )

    # 強度関数の事後平均
    plt.plot(
        t_values, lambda_mean,
        color='skyblue', linewidth=2,
        label='Posterior Mean'
    )

    # 強度関数の事後中央値
    #plt.plot(
    #    t_values, lambda_median,
    #    color='navy', linestyle='--', linewidth=2,
    #    label='Posterior Median'
    #)

    # 真の強度関数をプロット（もし既知の場合）
    if true_intensity_func is not None:
        lambda_true = true_intensity_func(t_values)
        plt.plot(
            t_values, lambda_true,
            color='red', linestyle='-', linewidth=2,
            label='True Intensity'
        )

    # 観測イベントのプロット
    plt.scatter(
        events, np.zeros_like(events),
        marker='o', color='gray', s=20,
        label='Observed Events', zorder=5
    )

    # プロットの装飾
    plt.xlabel('Time')
    plt.ylabel('Intensity λ(t)')
    plt.title('Posterior of Intensity Function with Observed Events')
    plt.legend(loc='upper right')
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()
