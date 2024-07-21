import matplotlib.pyplot as plt
import numpy as np
import scipy

def lamda(x, beta):
    """対数リンク関数を用いた強度関数

    Args:
        x (ndarray): 説明変数 (n, 1)
        beta (ndarray): パラメータ (3, 1)

    Returns:
        float: xにおける強度関数の値
    """
    X = np.vstack([np.ones(x.shape), x, x**2]).T
    return np.exp(X@beta)
    
def generate_IPP(T, beta_true, n0, Z, verbose=False):
    """１次元の非斉次ポアソン過程モデルに基づいてデータを生成する

    Args:
        T (int, float): 生成するデータの最大時刻
        beta_true (ndarray): 真のパラメータ (3, 1)
        n0 (int): 偽不在データの数
        Z (int): 重み付き尤度における偽不在に対する重み
        verbose (bool, optional): 途中の計算結果を出力するかどうか. Defaults to False.

    Returns:
        X (ndarray): 説明変数 (n, 3)
        y (ndarray): presenceかpseudo-absenceを表す二値変数 (n, )
        t (ndarray): イベントの発生時刻 (n, )
        z (ndarray): 重み付き尤度の重み (n, )
    """
    # 1次元の非斉次ポアソン過程モデル

    

    x_lin = np.linspace(0, T, 1000)
    intensity = lamda(x_lin, beta_true)
    
    if verbose:
        plt.plot(x_lin, intensity, c='black')
        plt.xlabel('Time')
        plt.ylabel('Intensity')
        plt.grid(True, linestyle='--', linewidth=0.5)
        plt.show()

    # サンプリングするイベントを格納するリスト
    events = []
    t = 0
    
    max_lambda = np.max(intensity)

    while t < T:
        # 次の候補イベントの間隔を指数分布からサンプリング
        u = np.random.uniform()
        t += -np.log(u) / max_lambda
        
        # 時刻 t での強度関数の値
        if t < T:
            lambda_t = lamda(np.array([t]), beta_true)
            
            # 強度関数の値に基づいてイベントを受け入れるか拒否するか決定
            if np.random.uniform() < lambda_t / max_lambda:
                events.append(t)

    t1 = np.array(events)    

    n1 = len(t1)
    
    # background data: 一様分布からのサンプリング
    t0 = np.random.uniform(0, T, n0)

    if verbose:
        plt.scatter(t0, np.zeros(n0), c='grey', alpha=0.4, label='pseudo-absence')
        plt.scatter(t1, np.zeros(len(t1)), c='green', alpha=0.3, label='presence')
        plt.plot(x_lin, intensity, c='black')
        plt.xlabel('Time')
        plt.ylabel('Intensity')
        plt.grid(True, linestyle='--', linewidth=0.5)
        plt.legend()
        plt.show()

    n = n1 + n0
    t = np.concatenate([t1, t0])
    X = np.vstack([np.ones(n), t, t**2]).T
    y = np.concatenate([np.ones(n1), np.zeros(n0)]).astype(int)

    z = np.concatenate([np.ones(len(t1)), np.ones(len(t0))*Z]).astype(int)
    
    return X, y, t, z
