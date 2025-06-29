import h5py
import numpy as np
import numpy.random as npr
from numpy.linalg import inv
from pypolyagamma import PyPolyaGamma
from tqdm import tqdm

from bayesian_statistics.utils import generate_IPP, multi_pgdraw_vectorized


def gibbs_sampling(
    X: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    b: np.ndarray,
    B: np.ndarray,
    SIM: int,
    batch_size: int,
    output_file: str,
) -> None:
    """ギブスサンプリングを実行し、結果をHDF5ファイルに保存する"""
    pg = PyPolyaGamma()
    beta_hat = npr.multivariate_normal(b, B)
    k = y - z / 2

    inv_B = inv(B)
    X_T = X.T

    with h5py.File(output_file, "w", libver="latest") as f:
        dset = f.create_dataset("beta_samples", (SIM, 3), dtype="float64", chunks=True)
        f.swmr_mode = True  # データセット作成後にSWMRモードを有効化

        for sim in tqdm(range(0, SIM, batch_size)):
            batch_size_actual = min(batch_size, SIM - sim)
            batch = np.zeros((batch_size_actual, 3))

            for i in tqdm(range(batch_size_actual)):
                Omega_diag = multi_pgdraw_vectorized(pg, z, X @ beta_hat)
                V = inv(X_T @ np.diag(Omega_diag) @ X + inv_B)
                m = V @ (X_T @ k + inv_B @ b)
                beta_hat = npr.multivariate_normal(m, V)
                batch[i] = beta_hat

            dset[sim : sim + batch_size_actual] = batch
            f.flush()  # データをファイルに書き込む

    print("サンプリングと保存が完了しました！")


# 使用例
SIM: int = 10000000
batch_size: int = 10000
output_file: str = "output/beta_samples.h5"

# 変数の定義
T: int = 5
beta_true: np.ndarray = np.array([4, 7, -4])
n0: int = 1000
Z: int = 10000

X, y, t, z = generate_IPP(T, beta_true, Z, background=True, n0=n0)

print(f"存在の数: {len(y[y == 1])}")
print(f"不在の数: {n0}")

# 平均ベクトルと共分散行列の設定
b: np.ndarray = np.zeros(3)
B: np.ndarray = np.diag(np.ones(3))

gibbs_sampling(X, y, z, b, B, SIM, batch_size, output_file)

# サンプルの読み込みと統計量の計算
with h5py.File(output_file, "r", libver="latest", swmr=True) as f:
    beta_samples = f["beta_samples"][:]

mean: np.ndarray = np.mean(beta_samples, axis=0)
var: np.ndarray = np.var(beta_samples, axis=0)

print("最終的な平均:", mean)
print("最終的な分散:", var)
