import numpy as np
from numpy.linalg import inv
import numpy.random as npr
from pypolyagamma import PyPolyaGamma
from tqdm import tqdm
import h5py

from utils import lamda, generate_IPP

def multi_pgdraw_vectorized(pg, B, C):
    """Vectorized version of multi_pgdraw"""
    return np.array([pg.pgdraw(b, c) for b, c in zip(B, C)])

def gibbs_sampling(X, y, z, b, B, SIM, batch_size, output_file):
    pg = PyPolyaGamma()
    beta_hat = npr.multivariate_normal(b, B)
    k = y - z/2
    
    inv_B = inv(B)
    X_T = X.T
    
    with h5py.File(output_file, 'w', libver='latest') as f:
        dset = f.create_dataset("beta_samples", (SIM, 3), dtype='float64', chunks=True)
        f.swmr_mode = True  # Enable SWMR mode after creating datasets
        
        for sim in tqdm(range(0, SIM, batch_size)):
            batch_size_actual = min(batch_size, SIM - sim)
            batch = np.zeros((batch_size_actual, 3))
            
            for i in tqdm(range(batch_size_actual)):
                Omega_diag = multi_pgdraw_vectorized(pg, z, X @ beta_hat)
                V = inv(X_T @ np.diag(Omega_diag) @ X + inv_B)
                m = V @ (X_T @ k + inv_B @ b)
                beta_hat = npr.multivariate_normal(m, V)
                batch[i] = beta_hat
            
            dset[sim:sim+batch_size_actual] = batch
            f.flush()  # Ensure data is written to file
    
    print("Finished sampling and saving!")

# 使用例
SIM = 10000000
batch_size = 10000
output_file = "output/beta_samples.h5"

# 変数の定義
T = 5
beta_true = np.array([4, 7, -4])
n0 = 1000
Z = 10000

X, y, t, z = generate_IPP(T, beta_true, n0, Z, verbose=True)

print(f"Number of presence: {len(y[y==1])}")
print(f"Number of absence: {n0}")

# 平均ベクトルと共分散行列の設定
b = np.zeros(3)
B = np.diag(np.ones(3))

gibbs_sampling(X, y, z, b, B, SIM, batch_size, output_file)

# サンプルの読み込みと統計量の計算
with h5py.File(output_file, 'r', libver='latest', swmr=True) as f:
    beta_samples = f['beta_samples'][:]

mean = np.mean(beta_samples, axis=0)
var = np.var(beta_samples, axis=0)

print("Final mean:", mean)
print("Final variance:", var)