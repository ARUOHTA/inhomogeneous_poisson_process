# %%
import matplotlib.pyplot as plt
import numpy as np
from   numpy.linalg import inv
import numpy.random as npr
from   pypolyagamma import PyPolyaGamma
from scipy.stats import multivariate_normal

import matplotlib.cm as cm
from tqdm import tqdm

import pandas as pd
import statsmodels.api as sm
import scipy
import scipy.optimize as opt
import seaborn as sns


# %%
from utils import lamda, generate_IPP

# %% [markdown]
# ## データの生成

# %%
T = 5
beta_true = np.array([4, 7, -4])
n0 = 1000
Z = 10000

X, y, t, z = generate_IPP(T, beta_true, n0, Z, verbose=True)

print(f"Number of presence: {len(y[y==1])}")
print(f"Number of absence: {n0}")


# %% [markdown]
# ## ベイズ推定

# %%
from   numpy.linalg import inv
import numpy.random as npr
from   pypolyagamma import PyPolyaGamma

def multi_pgdraw(pg, B, C):
    """Utility function for calling `pgdraw` on every pair in vectors B, C.
    """
    return np.array([pg.pgdraw(b, c) for b, c in zip(B, C)])

# %%
# 平均ベクトルと共分散行列の設定
b = np.zeros(3)
B = np.diag(np.ones(3))

# Peform Gibbs sampling for SIM iterations.
pg         = PyPolyaGamma()
SIM          = 10000000
beta_hat   = npr.multivariate_normal(b, B)
k          = y - z/2

beta_list = []
for _ in tqdm(range(SIM)):
    # ω ~ PG(1, x*β).
    Omega_diag = multi_pgdraw(pg, z, X @ beta_hat)
    # β ~ N(m, V).
    V         = inv(X.T @ np.diag(Omega_diag) @ X + inv(B))
    m         = np.dot(V, X.T @ k + inv(B) @ b)
    beta_hat  = npr.multivariate_normal(m, V)
    beta_list.append(beta_hat)
    
beta_list = np.array(beta_list)

# %%
# beta_listをcsvに保存
pd.DataFrame(beta_list).to_csv(f"beta_list_{SIM}.csv")
print("finished!")