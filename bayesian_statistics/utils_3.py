from abc import ABC, abstractmethod
from typing import Callable

import numpy as np
from pypolyagamma import PyPolyaGamma


def multi_pgdraw_vectorized(
    pg: PyPolyaGamma, B: np.ndarray, C: np.ndarray
) -> np.ndarray:
    """multi_pgdraw のベクトル化バージョン"""
    return np.array([pg.pgdraw(b, c) for b, c in zip(B, C)])


# 抽象クラスの定義
class DesignMatrix(ABC):
    @abstractmethod
    def create(self, x: np.ndarray) -> np.ndarray:
        pass


# 多項式設計行列クラス
class PolynomialDesignMatrix(DesignMatrix):
    def __init__(self, degree: int = 2, include_interaction: bool = True):
        self.degree = degree
        self.include_interaction = include_interaction

    def create(self, x: np.ndarray) -> np.ndarray:
        if x.ndim == 1:
            x = x[:, np.newaxis]
        N, D = x.shape
        terms = [np.ones(N)]  # バイアス項

        # 各次元の多項式項を追加
        for d in range(D):
            for deg in range(1, self.degree + 1):
                terms.append(x[:, d] ** deg)

        # 相互作用項を追加
        if self.include_interaction and D > 1:
            for i in range(D):
                for j in range(i + 1, D):
                    terms.append(x[:, i] * x[:, j])

        return np.column_stack(terms)


# サイン・コサイン設計行列クラス
class TrigonometricDesignMatrix(DesignMatrix):
    def __init__(self, omega: float = np.pi / 5):
        self.omega = omega

    def create(self, x: np.ndarray) -> np.ndarray:
        if x.ndim == 1:
            x = x[:, np.newaxis]
        N, D = x.shape
        terms = [np.ones(N)]  # バイアス項

        # 各次元のサイン・コサイン項を追加
        for d in range(D):
            terms.append(np.sin(self.omega * x[:, d]))
            terms.append(np.cos(self.omega * x[:, d]))

        return np.column_stack(terms)


# 複合設計行列クラス
class CompositeDesignMatrix(DesignMatrix):
    def __init__(self, design_matrices: list):
        self.design_matrices = design_matrices

    def create(self, x: np.ndarray) -> np.ndarray:
        terms = []
        for dm in self.design_matrices:
            dm_terms = dm.create(x)
            # バイアス項の重複を避ける
            if len(terms) > 0 and np.array_equal(dm_terms[:, 0], np.ones(len(x))):
                dm_terms = dm_terms[:, 1:]
            terms.append(dm_terms)
        return np.column_stack(terms)


# IntensityFunction クラスの修正
class IntensityFunction:
    def __init__(
        self,
        design_matrix: DesignMatrix,
        beta: np.ndarray,
        lambda_star: float,
        link_function: Callable[[np.ndarray], np.ndarray] = None,
    ):
        self.design_matrix = design_matrix
        self.beta = beta
        self.lambda_star = lambda_star
        if link_function is None:
            # デフォルトでロジスティックリンク関数を使用
            self.link_function = lambda eta: 1 / (1 + np.exp(-eta))
        else:
            self.link_function = link_function

    def copy(self):
        return IntensityFunction(
            design_matrix=self.design_matrix,
            beta=self.beta.copy(),
            lambda_star=self.lambda_star,
            link_function=self.link_function,
        )

    def q(self, x: np.ndarray) -> np.ndarray:
        """q(x) を計算"""
        X = self.design_matrix.create(x)
        eta = X @ self.beta
        return self.link_function(eta)

    def lambda_func(self, x: np.ndarray) -> np.ndarray:
        """λ(x) = q(x) * λ* を計算"""
        return self.q(x) * self.lambda_star

    def update_beta(self, beta: np.ndarray):
        self.beta = beta

    def update_lambda_star(self, lambda_star: float):
        self.lambda_star = lambda_star


# イベント生成関数などは変更なし


# テストコード
def test_design_matrix():
    # テストデータの作成
    x_1d = np.linspace(0, 10, 5)
    x_2d = np.array([[i, j] for i in range(5) for j in range(5)])

    # 多項式設計行列のテスト
    poly_dm = PolynomialDesignMatrix(degree=2, include_interaction=True)
    W_1d_poly = poly_dm.create(x_1d)
    W_2d_poly = poly_dm.create(x_2d)

    assert W_1d_poly.shape == (5, 3)  # バイアス + 1次項 + 2次項
    assert W_2d_poly.shape == (25, 7)  # バイアス + 各次元の1次・2次項 + 相互作用項

    # サイン・コサイン設計行列のテスト
    trig_dm = TrigonometricDesignMatrix(omega=np.pi / 5)
    W_1d_trig = trig_dm.create(x_1d)
    W_2d_trig = trig_dm.create(x_2d)

    assert W_1d_trig.shape == (5, 3)  # バイアス + sin + cos
    assert W_2d_trig.shape == (25, 5)  # バイアス + 各次元の sin, cos

    # 複合設計行列のテスト
    composite_dm = CompositeDesignMatrix([poly_dm, trig_dm])
    W_1d_composite = composite_dm.create(x_1d)
    W_2d_composite = composite_dm.create(x_2d)

    # 1次元の場合の特徴量数を計算
    num_features_1d = (
        W_1d_poly.shape[1] + W_1d_trig.shape[1] - 1
    )  # バイアス項の重複を避ける
    assert W_1d_composite.shape == (5, num_features_1d)

    # 2次元の場合の特徴量数を計算
    num_features_2d = (
        W_2d_poly.shape[1] + W_2d_trig.shape[1] - 1
    )  # バイアス項の重複を避ける
    assert W_2d_composite.shape == (25, num_features_2d)

    print("All tests passed for design matrices.")


# テストの実行
test_design_matrix()
