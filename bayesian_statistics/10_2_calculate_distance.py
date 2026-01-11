import os
from heapq import heappop, heappush

import geopandas as gpd
import numpy as np
import polars as pl
from shapely import wkt
from shapely.geometry import Point
from tqdm import tqdm


def prepare_data_for_dijkstra(
    df: pl.DataFrame, target_geometry, use_average_angle=False
):
    """
    dfは以下の列を持つと仮定:
    - grid_x: int
    - grid_y: int
    - geometry: shapely Polygon（WGS84など位置参照不要な単純な多角形）
    - travel_time_east, travel_time_west, travel_time_north, travel_time_south: float or NaN

    前処理として:
    - (grid_x, grid_y)から1次元インデックスへの対応付け
    - travel_timeを2D配列化
    - target_geometryとの交差セルを特定
    """

    # grid_x, grid_yが0-basedの連続整数と仮定
    max_x = df["grid_x"].max()
    max_y = df["grid_y"].max()

    width = max_x + 1
    height = max_y + 1

    # 1次元インデックスへのマッピング関数
    def idx(x, y):
        return y * width + x

    n_cells = width * height

    # min_costs用に∞で初期化
    # float64で良い
    min_costs = np.full(n_cells, np.inf, dtype=np.float64)

    # travel_timeをまとめる
    # 方向は0:east,1:west,2:north,3:southとする
    directions = [("east", 1, 0), ("west", -1, 0), ("north", 0, 1), ("south", 0, -1)]
    dir_map = {"east": 0, "west": 1, "north": 2, "south": 3}

    travel_times = np.full((4, height, width), np.nan, dtype=np.float64)

    # dfをイテレートしてtravel_timeを格納
    # geometryやtravel_time_*列を参照するため一度to_dictやiterrowsを使う
    for row in tqdm(df.iter_rows(named=True), total=len(df)):
        gx = row["grid_x"]
        gy = row["grid_y"]
        # 各方向のtravel_timeをセット
        for d in directions:
            dname = d[0]
            val = (
                row[f"travel_time_{dname}"]
                if not use_average_angle
                else row["travel_time"]
            )
            if not np.isnan(val):
                travel_times[dir_map[dname], gy, gx] = val

    # target_geometryと交差するセルを特定
    print("converting dataframe to WKT...")
    gdf = df.to_pandas()
    gdf.geometry = gdf.geometry.apply(wkt.loads)
    gdf = gpd.GeoDataFrame(gdf, geometry="geometry")
    gdf.crs = "EPSG:4326"

    # intersectionの計算
    print("calculating intersection...")
    df = df.with_columns(pl.Series(gdf.intersects(target_geometry)).alias("intersect"))

    # intersect列がTrueのセルのmin_costsを0に設定
    starting_indices = []
    df_intersect = df.filter(pl.col("intersect") == True)
    for row in tqdm(df_intersect.iter_rows(named=True), total=len(df_intersect)):
        gx = row["grid_x"]
        gy = row["grid_y"]
        cell_idx = idx(gx, gy)
        min_costs[cell_idx] = 0.0
        starting_indices.append(cell_idx)

    if len(starting_indices) == 0:
        raise ValueError("対象ジオメトリと交差するセルが見つかりませんでした。")

    return min_costs, travel_times, width, height, starting_indices, directions


def run_dijkstra(min_costs, travel_times, width, height, directions):
    def idx(x, y):
        return y * width + x

    visited = np.full_like(min_costs, False, dtype=bool)
    queue = []
    total_cells = width * height  # グリッド内の全セル数

    # 進捗バー用のpbarを作成
    with tqdm(total=total_cells, desc="Processing cells") as pbar:
        # 開始セルをキューに登録
        visited_count = 0
        for i in range(len(min_costs)):
            if min_costs[i] == 0.0:
                y, x = divmod(i, width)
                heappush(queue, (0.0, x, y))
                visited[idx(x, y)] = True
                visited_count += 1
        pbar.update(visited_count)

        while queue:
            current_cost, cx, cy = heappop(queue)
            cidx = idx(cx, cy)

            # 隣接セルへ
            for d_i, (dname, dx, dy) in enumerate(directions):
                nx, ny = cx + dx, cy + dy
                # 範囲外ならスキップ
                if nx < 0 or nx >= width or ny < 0 or ny >= height:
                    continue

                # travel_timeがNaNなら移動不可
                t = travel_times[d_i, cy, cx]
                if np.isnan(t):
                    continue

                nidx = idx(nx, ny)
                if not visited[nidx]:  # 未訪問の場合のみカウント
                    visited[nidx] = True
                    pbar.update(1)

                new_cost = current_cost + t
                if new_cost < min_costs[nidx]:
                    min_costs[nidx] = new_cost
                    heappush(queue, (new_cost, nx, ny))

    return min_costs


# 2つをまとめて実行する関数
def calculate_minimum_cost_to_geometry(
    df, target_geometry, geometry_name="", use_average_angle=False
):
    min_costs, travel_times, width, height, starting_indices, directions = (
        prepare_data_for_dijkstra(
            df, target_geometry, use_average_angle=use_average_angle
        )
    )

    print(f"number of starting cells: {len(starting_indices)}")

    min_costs = run_dijkstra(min_costs, travel_times, width, height, directions)

    min_costs = np.where(min_costs == np.inf, np.nan, min_costs)

    # grid_x, grid_y, min_cost_minutsという列を持つpl.DataFrameを作成
    min_costs_df = pl.DataFrame(
        {
            "grid_x": np.tile(np.arange(width), height),
            "grid_y": np.repeat(np.arange(height), width),
            f"cost_{geometry_name}": min_costs,
        }
    )

    return min_costs_df


if __name__ == "__main__":
    data_dir = "data"

    # データの読み込み
    df_elevation = pl.read_csv(os.path.join(data_dir, "10_1_gdf_elevation_tobler.csv"))
    # 川のポリゴンを読み込む
    df_river_stream = pl.read_csv(os.path.join(data_dir, "9_gdf_river_stream.csv"))

    river_polygon = gpd.GeoSeries.from_wkt(
        df_river_stream.get_column("geometry").to_numpy()
    ).union_all()

    polygon_dict = {
        "kouzu": Point(139.1517940530124, 34.214991203764754),
        "shinshu": Point(138.1431305514158, 36.14658805493071),
        "hakone": Point(139.0446125901514, 35.221867157762105),
        "takahara": Point(139.7766240432928, 36.900342242149065),
        "river": river_polygon,
    }

    for key, value in polygon_dict.items():
        print(f"calculating distance to {key}...")
        min_costs_df = calculate_minimum_cost_to_geometry(
            df_elevation, value, geometry_name=key, use_average_angle=True
        )
        df_elevation = df_elevation.join(
            min_costs_df, on=["grid_x", "grid_y"], how="left"
        )
        print("done.")
        print("")

    # 保存
    df_elevation.write_csv(
        os.path.join(data_dir, "10_2_gdf_elevation_with_costs.csv")
    )
