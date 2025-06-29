from typing import Optional

import numpy as np
import polars as pl

##############################################
# 1) ダミーの df_elevation を作る
##############################################
# ここでは「x, y が 5次メッシュの中心」という想定で等間隔に並んでいる例を簡単に作成
#   x: 130.0000, 130.0100, 130.0200, ...
#   y:  30.0000,  30.0100,  30.0200, ...
# 合計 5×5 = 25 グリッド程度の小さな例

x_vals = [
    130.0000 + 0.01 * i for i in range(5)
]  # 130.0000, 130.0100, 130.0200, 130.0300, 130.0400
y_vals = [
    30.0000 + 0.01 * j for j in range(5)
]  #  30.0000,  30.0100,  30.0200,  30.0300,  30.0400

# polars DataFrame で (x, y) の組み合わせを作成
df_elevation = pl.DataFrame({"x": x_vals}).join(
    pl.DataFrame({"y": y_vals}), how="cross"
)


##############################################
# 2) 先に示した「中心座標からステップを求める」関数
##############################################
def prepare_grid_info_center(df_elevation: pl.DataFrame):
    """
    df_elevation の x, y が「メッシュの中心座標」。
    最小中心, 最大中心, ステップ幅 (delta_x, delta_y) を算出して返す。
    """
    x_unique = df_elevation.select("x").unique().to_series().to_list()
    y_unique = df_elevation.select("y").unique().to_series().to_list()
    x_unique.sort()
    y_unique.sort()

    # ステップ幅 = 隣り合う中心座標の差
    delta_x = x_unique[1] - x_unique[0]
    delta_y = y_unique[1] - y_unique[0]

    return {
        "x_min_center": x_unique[0],
        "x_max_center": x_unique[-1],
        "y_min_center": y_unique[0],
        "y_max_center": y_unique[-1],
        "delta_x": delta_x,
        "delta_y": delta_y,
    }


def get_grid_xy_for_center_mesh(
    lon: float, lat: float, grid_info: dict
) -> Optional[tuple[int, int]]:
    """
    (lon, lat) がどのメッシュ (grid_x, grid_y) に含まれるかを O(1) で求める。
    範囲外なら None を返す。
    """
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

    # 範囲外なら None
    if not (x_left_bound <= lon <= x_right_bound):
        return None
    if not (y_bottom_bound <= lat <= y_top_bound):
        return None

    # (lon, lat) が何番目の格子か = floor( (lon - 左境界) / dx ), floor( (lat - 下境界) / dy )
    gx = int(np.floor((lon - x_left_bound) / dx))
    gy = int(np.floor((lat - y_bottom_bound) / dy))

    return gx, gy


##############################################
# 3) テスト用データを用意
##############################################
# 「入力 (lon, lat)」「期待される (grid_x, grid_y)」のペアで定義
#   grid_x, grid_y の定義:
#   - 左下のメッシュ(130.0000, 30.0000) を (0, 0) とする
#   - x方向(経度)に大きくなるほど grid_x が増える
#   - y方向(緯度)に大きくなるほど grid_y が増える
#
# メッシュの幅は dx = 0.01, dy = 0.01
# 中心(130.0000, 30.0000) の左境界は 129.995, 下境界は 29.995
#
# メッシュ一覧 (中心):
#   grid_x=0 → x=130.0000
#   grid_x=1 → x=130.0100
#   ...
#   grid_x=4 → x=130.0400
#   grid_y=0 → y=30.0000
#   grid_y=4 → y=30.0400
#
# つまり、例えば (x=130.0000, y=30.0000) は (grid_x=0, grid_y=0)
# そのメッシュの境界 → x ∈ [129.995, 130.005), y ∈ [29.995, 30.005)
#
# テストポイント例:
#   1) メッシュ中心ピッタリ
#   2) メッシュ境界ギリギリ
#   3) 範囲外
#   4) 連続的にいくつか
#   5) ランダム
##############################################

# 1) メッシュ中心ピッタリに該当する事例
test_points_center = [
    # (lon, lat, expected (gx, gy))
    (130.0000, 30.0000, (0, 0)),  # 左下
    (130.0100, 30.0000, (1, 0)),  # 1つ右
    (130.0200, 30.0100, (2, 1)),  # さらに右、上
    (130.0400, 30.0400, (4, 4)),  # 右上端
]

# 2) メッシュ境界付近
#   例： (130.0050, 30.0049) は、
#        x=130.0050 → dx=0.01 なら、130.0050 - 129.995 = 0.0100 なので、ちょうど右メッシュの左境界。
#        ただし浮動小数点誤差には要注意。
test_points_boundary = [
    # このあたりは人間で「どちらに入るか」計算してみる
    (129.9950, 29.9950, (0, 0)),  # ちょうど (grid_x=0, grid_y=0) メッシュの左下角
    (130.0049999, 30.0049999, (0, 0)),  # ほぼギリギリ(境界内)→(0,0)
    (
        130.0050,
        30.0049999,
        (1, 0),
    ),  # x=130.0050 は(0, 0)メッシュの右境界・左境界スレスレ
    (
        130.0049999,
        30.0050,
        (0, 1),
    ),  # y=30.0050 は(0,0)メッシュの上境界と(0,1)の下境界ギリギリ
    (130.0050, 30.0050, (1, 1)),  # 境界上、誤差のさじ加減で(1,1)に入る例
]

# 3) 範囲外
test_points_out_of_range = [
    (129.9949, 30.0000, None),  # 左にはみ出す
    (130.0000, 29.9949, None),  # 下にはみ出す
    (130.0401, 30.0000, None),  # 右にはみ出す (右境界は 130.0450)
    (130.0000, 30.0451, None),  # 上にはみ出す (上境界は 30.0450)
]

# 4) 連続的にチェック
test_points_sequential = []
for i in range(5):  # grid_x: 0 ~ 4
    for j in range(5):  # grid_y: 0 ~ 4
        center_x = 130.0000 + 0.01 * i
        center_y = 30.0000 + 0.01 * j
        # (中心からちょっとズラした座標)→確実に同じメッシュに入るはず
        test_points_sequential.append((center_x + 0.002, center_y - 0.001, (i, j)))

# 5) ランダムに少し
rng = np.random.default_rng(seed=42)
test_points_random = []
for _ in range(5):
    # メッシュ全体の “左下境界” = 129.995, “右上境界” = 130.045
    # そこからランダムに
    rx = rng.uniform(129.995, 130.045)
    ry = rng.uniform(29.995, 30.045)
    # 期待値は計算してもよいが、ここでは実行時に確認
    test_points_random.append((rx, ry, None))


##############################################
# 4) テスト実行
##############################################
if __name__ == "__main__":
    grid_info = prepare_grid_info_center(df_elevation)

    print("=== Grid Info ===")
    for k, v in grid_info.items():
        print(f"{k}: {v}")
    print()

    def run_test_case(lon, lat, expected):
        """(lon, lat) -> get_grid_xy_for_center_mesh() して結果を表示。"""
        result = get_grid_xy_for_center_mesh(lon, lat, grid_info)
        if expected is not None:
            # 期待値がある場合はアサート
            assert result == expected, (
                f"Test fail: (lon, lat)=({lon}, {lat}), expected={expected}, got={result}"
            )
        print(f"OK: (lon={lon}, lat={lat}) -> {result} (expected={expected})")

    # セットごとに実行
    print("=== 1) メッシュ中心ピッタリ ===")
    for lon, lat, exp in test_points_center:
        run_test_case(lon, lat, exp)
    print()

    print("=== 2) メッシュ境界付近 ===")
    for lon, lat, exp in test_points_boundary:
        run_test_case(lon, lat, exp)
    print()

    print("=== 3) 範囲外 ===")
    for lon, lat, exp in test_points_out_of_range:
        run_test_case(lon, lat, exp)
    print()

    print("=== 4) 連続チェック (各メッシュの近傍座標) ===")
    for lon, lat, exp in test_points_sequential:
        run_test_case(lon, lat, exp)
    print()

    print("=== 5) ランダム座標 ===")
    for lon, lat, _ in test_points_random:
        # ランダムは事前に期待値を計算していないので、結果だけ表示
        result = get_grid_xy_for_center_mesh(lon, lat, grid_info)
        print(f"(lon={lon:.6f}, lat={lat:.6f}) -> result={result}")
