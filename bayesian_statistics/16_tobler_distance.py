# %%
import os
import pickle

import geopandas as gpd
import numpy as np
import polars as pl
from shapely import wkt
from tqdm import tqdm


# %%
def create_site_coords(df: pl.DataFrame) -> np.ndarray:
    """
    遺跡の座標をラジアンに変換

    Parameters
    ----------
    df : pl.DataFrame
        入力データフレーム

    Returns
    -------
    np.ndarray
        遺跡の座標（ラジアン）
    """
    # 座標をラジアンに変換
    coords = (
        df.select(
            [
                (pl.col("遺跡ID")),
                (pl.col("緯度") * np.pi / 180).alias("lat_rad"),
                (pl.col("経度") * np.pi / 180).alias("lon_rad"),
            ]
        )
        .unique(subset=["遺跡ID"])
        .sort("遺跡ID")
    )

    # 座標と出土数を numpy 配列に変換
    site_coords = np.column_stack(
        [coords["lat_rad"].to_numpy(), coords["lon_rad"].to_numpy()]
    )

    return site_coords


# %%
from pathlib import Path

data_dir = str(Path(__file__).parent.parent / "data") + "/"

SCHEMA = pl.Schema(
    {
        "grid_x": pl.Int64,
        "x": pl.Float64,
        "grid_y": pl.Int64,
        "y": pl.Float64,
        "mesh_code_5th": pl.Int64,
        "average_elevation": pl.Float64,
        "maximum_elevation": pl.Float64,
        "minimum_elevation": pl.Float64,
        "minimum_elevation_code": pl.Int64,
        "maximum_slope_angle": pl.Float64,
        "maximum_slope_direction": pl.Int64,
        "minimum_slope_angle": pl.Float64,
        "minimum_slope_direction": pl.Int64,
        "average_slope_angle": pl.Float64,
        "geometry": pl.Utf8,
        "is_sea": pl.Boolean,
        "walking_velocity": pl.Float64,
        "travel_time": pl.Float64,
        "elevation_diff_east": pl.Float64,
        "angle_east": pl.Float64,
        "walking_velocity_east": pl.Float64,
        "travel_time_east": pl.Float64,
        "elevation_diff_west": pl.Float64,
        "angle_west": pl.Float64,
        "walking_velocity_west": pl.Float64,
        "travel_time_west": pl.Float64,
        "elevation_diff_north": pl.Float64,
        "angle_north": pl.Float64,
        "walking_velocity_north": pl.Float64,
        "travel_time_north": pl.Float64,
        "elevation_diff_south": pl.Float64,
        "angle_south": pl.Float64,
        "walking_velocity_south": pl.Float64,
        "travel_time_south": pl.Float64,
        "cost_kouzu": pl.Float64,
        "cost_shinshu": pl.Float64,
        "cost_hakone": pl.Float64,
        "cost_takahara": pl.Float64,
        "cost_river": pl.Float64,
        "x_meter": pl.Float64,
        "y_meter": pl.Float64,
    }
)

df_elevation = pl.read_csv(
    os.path.join(data_dir, "11_gdf_elevation.csv"),
    null_values=["nan"],  # "nan"をnullとして扱う
)
df_elevation = (
    # データ型を変換
    df_elevation.cast(SCHEMA)
    # 全く欠損していない行をis_validとする
    .with_columns(
        [
            pl.all_horizontal(
                [~pl.col(col).is_null() for col in df_elevation.columns]
            ).alias("is_valid")
        ]
    )
)

df_obsidian = pl.read_csv(os.path.join(data_dir, "11_gdf_obsidian.csv"))
df_sites = pl.read_csv(os.path.join(data_dir, "11_gdf_sites.csv"))

time_period_name = {0: "早期・早々期", 1: "前期", 2: "中期", 3: "後期", 4: "晩期"}

origin_order = ["神津島", "信州", "箱根", "高原山", "その他"]

sigma = 500
sigma_for_sites = 0.1

# =======================================================================================

# メッシュグリッドを作成
lon_mesh, lat_mesh = np.meshgrid(
    df_elevation["x"].unique().sort(),  # 経度の一意な値
    df_elevation["y"].unique().sort(),  # 緯度の一意な値
)

print(lon_mesh.shape, lat_mesh.shape)

# 遺跡の位置: (遺跡数, 2)
site_coords = create_site_coords(df_obsidian)

# グリッドの位置座標: (グリッド数, 2)
grid_coords = np.column_stack(
    [lat_mesh.ravel() * np.pi / 180, lon_mesh.ravel() * np.pi / 180]
)


# 海岸線のfrictionを高くする
# 海岸線の定義: is_seaがFalseかつ、elevation_diff_{direction}の4つのうち1~3つがNullである場合

df_elevation = df_elevation.with_columns(
    [
        (
            (
                (
                    pl.col("elevation_diff_north").is_null().cast(pl.Int64)
                    + pl.col("elevation_diff_east").is_null().cast(pl.Int64)
                    + pl.col("elevation_diff_south").is_null().cast(pl.Int64)
                    + pl.col("elevation_diff_west").is_null().cast(pl.Int64)
                ).is_in([1, 2, 3])
            )
            & (pl.col("is_sea") == False)
        ).alias("is_coast")
    ]
)

# 海岸沿いのfrictionを高くする

# 極端に高くする場所のポリゴン一覧

from shapely.geometry import Polygon

izu = {
    "type": "伊豆半島",
    "coordinates": [
        [
            [138.7518662555113, 35.042617144097115],
            [138.6707917599528, 34.56053181508489],
            [139.11491388354187, 34.52700264837797],
            [139.2211370720799, 35.07275597779107],
            [139.3281342140192, 35.305497704601535],
            [139.31655815998965, 35.34834874650015],
            [138.9909080638121, 35.204024799263756],
            [138.86887475645986, 35.080895934583104],
            [138.85793636642387, 35.07985935793645],
            [138.7518662555113, 35.042617144097115],
        ]
    ],
}
bousou = {
    "type": "房総半島",
    "coordinates": [
        [
            [139.82365476349196, 35.2922877248991],
            [139.6952824547859, 34.91167841077183],
            [139.97838837375795, 34.86098321580854],
            [140.16386293220444, 35.07131446513346],
            [140.4122282892051, 35.136203527893834],
            [140.4358783809528, 35.34544112887992],
            [139.82365476349196, 35.2922877248991],
        ]
    ],
}

izu_polygon = Polygon([[x, y] for x, y in izu["coordinates"][0]])
bousou_polygon = Polygon([[x, y] for x, y in bousou["coordinates"][0]])

# 海岸沿いは、travel_timeを、travel_time * multiple_normal倍にする
# 伊豆半島と房総半島の海岸沿いは、travel_timeを、travel_time * multiple_izu_bousou倍にする

print("converting dataframe to WKT...")
gdf = df_elevation.to_pandas()
gdf.geometry = gdf.geometry.apply(wkt.loads)
gdf = gpd.GeoDataFrame(gdf, geometry="geometry")
gdf.crs = "EPSG:4326"

# intersectionの計算
print("calculating intersection...")
df_elevation = df_elevation.with_columns(
    pl.Series(gdf.intersects(izu_polygon)).alias("intersect_izu")
)
df_elevation = df_elevation.with_columns(
    pl.Series(gdf.intersects(bousou_polygon)).alias("intersect_bousou")
)

multiple_normal = 50
multiple_izu_bousou = 1000

df_elevation = df_elevation.with_columns(
    [
        # average
        pl.when(
            # 海岸沿いかつ、伊豆半島または房総半島と交差している場合
            pl.col("is_coast") & (pl.col("intersect_izu") | pl.col("intersect_bousou"))
        )
        .then(pl.col("travel_time") * multiple_izu_bousou)
        .when(
            # 海岸沿いかつ、伊豆半島または房総半島と交差していない場合
            pl.col("is_coast")
            & ((~pl.col("intersect_izu")) & (~pl.col("intersect_bousou")))
        )
        .then(pl.col("travel_time") * multiple_normal)
        .otherwise(pl.col("travel_time"))
        .alias("travel_time"),
        # north
        pl.when(
            # 海岸沿いかつ、伊豆半島または房総半島と交差している場合
            pl.col("is_coast") & (pl.col("intersect_izu") | pl.col("intersect_bousou"))
        )
        .then(
            pl.when(pl.col("elevation_diff_north").is_null())
            .then(pl.col("travel_time_north") * multiple_izu_bousou)
            .otherwise(pl.col("travel_time_north"))
        )
        .when(
            # 海岸沿いかつ、伊豆半島または房総半島と交差していない場合
            pl.col("is_coast")
            & ((~pl.col("intersect_izu")) & (~pl.col("intersect_bousou")))
        )
        .then(
            pl.when(pl.col("elevation_diff_north").is_null())
            .then(pl.col("travel_time_north") * multiple_normal)
            .otherwise(pl.col("travel_time_north"))
        )
        .otherwise(pl.col("travel_time_north"))
        .alias("travel_time_north"),
        # south
        pl.when(
            pl.col("is_coast") & (pl.col("intersect_izu") | pl.col("intersect_bousou"))
        )
        .then(
            pl.when(pl.col("elevation_diff_south").is_null())
            .then(pl.col("travel_time_south") * multiple_izu_bousou)
            .otherwise(pl.col("travel_time_south"))
        )
        .when(
            pl.col("is_coast")
            & ((~pl.col("intersect_izu")) & (~pl.col("intersect_bousou")))
        )
        .then(
            pl.when(pl.col("elevation_diff_south").is_null())
            .then(pl.col("travel_time_south") * multiple_normal)
            .otherwise(pl.col("travel_time_south"))
        )
        .otherwise(pl.col("travel_time_south"))
        .alias("travel_time_south"),
        # east
        pl.when(
            pl.col("is_coast") & (pl.col("intersect_izu") | pl.col("intersect_bousou"))
        )
        .then(
            pl.when(pl.col("elevation_diff_east").is_null())
            .then(pl.col("travel_time_east") * multiple_izu_bousou)
            .otherwise(pl.col("travel_time_east"))
        )
        .when(
            pl.col("is_coast")
            & ((~pl.col("intersect_izu")) & (~pl.col("intersect_bousou")))
        )
        .then(
            pl.when(pl.col("elevation_diff_east").is_null())
            .then(pl.col("travel_time_east") * multiple_normal)
            .otherwise(pl.col("travel_time_east"))
        )
        .otherwise(pl.col("travel_time_east"))
        .alias("travel_time_east"),
        # west
        pl.when(
            pl.col("is_coast") & (pl.col("intersect_izu") | pl.col("intersect_bousou"))
        )
        .then(
            pl.when(pl.col("elevation_diff_west").is_null())
            .then(pl.col("travel_time_west") * multiple_izu_bousou)
            .otherwise(pl.col("travel_time_west"))
        )
        .when(
            pl.col("is_coast")
            & ((~pl.col("intersect_izu")) & (~pl.col("intersect_bousou")))
        )
        .then(
            pl.when(pl.col("elevation_diff_west").is_null())
            .then(pl.col("travel_time_west") * multiple_normal)
            .otherwise(pl.col("travel_time_west"))
        )
        .otherwise(pl.col("travel_time_west"))
        .alias("travel_time_west"),
    ]
)


# %%
def get_site_grid_coords(
    df_obsidian: pl.DataFrame, df_elevation: pl.DataFrame
) -> np.ndarray:
    """
    遺跡の座標をグリッド座標に変換

    Parameters
    ----------
    df_obsidian : pl.DataFrame
        遺跡のデータフレーム
    df_elevation : pl.DataFrame
        標高のデータフレーム

    Returns
    -------
    np.ndarray
        遺跡のグリッド座標
    """

    x_spacing = df_elevation["x"].unique().sort().diff().mode()[0]
    y_spacing = df_elevation["y"].unique().sort().diff().mode()[0]
    x_first_center = df_elevation["x"].min()
    y_first_center = df_elevation["y"].min()

    coords = (
        df_obsidian.select([(pl.col("遺跡ID")), (pl.col("緯度")), (pl.col("経度"))])
        .unique(subset=["遺跡ID"])
        .sort("遺跡ID")
        .with_columns(
            [
                ((pl.col("経度") - (x_first_center - (x_spacing / 2))) / x_spacing)
                .cast(pl.Int64)
                .alias("grid_x"),
                ((pl.col("緯度") - (y_first_center - (y_spacing / 2))) / y_spacing)
                .cast(pl.Int64)
                .alias("grid_y"),
            ]
        )
    )

    site_grid_coords = np.column_stack(
        [coords["grid_x"].to_numpy(), coords["grid_y"].to_numpy()]
    )

    return site_grid_coords


# %%
# ====================
# polygonを保存
# ====================

# polygon_list = df_elevation.filter(pl.col("geometry").is_not_null()).select(pl.col("geometry")).to_numpy().squeeze().tolist()
# polygon_list = [wkt.loads(p) for p in polygon_list]
# polygon = unary_union(polygon_list)
# gpd.GeoSeries([polygon]).to_file(os.path.join(data_dir, "polygon.shp"))

# %%
# ====================
# 以前の関数を再利用して試してみる
# ====================

import os
from heapq import heappop, heappush

import geopandas as gpd
import numpy as np
import polars as pl
from shapely import wkt


def process_travel_times(df, height, width, use_average_angle):
    # 方向の定義
    directions = ["east", "west", "north", "south"]

    if use_average_angle:
        # 平均角度を使用する場合
        # 各方向に同じtravel_timeを使用
        result = df.select(
            [
                pl.col("grid_x"),
                pl.col("grid_y"),
                pl.col("travel_time").alias("value"),
                pl.lit(list(range(4))).alias("direction"),
            ]
        ).explode("direction")
    else:
        # 方向ごとの個別のtravel_timeを使用
        result = df.melt(
            id_vars=["grid_x", "grid_y"],
            value_vars=[f"travel_time_{d}" for d in directions],
            value_name="value",
        ).with_columns(
            [
                pl.col("variable")
                .str.replace("travel_time_", "")
                .replace(dict(zip(directions, range(4))))
                .cast(pl.Int64)
                .alias("direction")
            ]
        )

    # nanを除外し、3次元配列に変換
    travel_times = np.full((4, height, width), np.nan, dtype=np.float64)

    # DataFrameから値を抽出して配列に格納
    valid_data = result.filter(pl.col("value").is_not_nan()).select(
        [pl.col("direction"), pl.col("grid_y"), pl.col("grid_x"), pl.col("value")]
    )

    # convert to numpy array for efficient indexing
    indices = valid_data.select(["direction", "grid_y", "grid_x"]).to_numpy()
    values = valid_data.select("value").to_numpy().flatten()

    # efficient array assignment
    travel_times[indices[:, 0], indices[:, 1], indices[:, 2]] = values

    return travel_times


def prepare_data_for_dijkstra(df: pl.DataFrame, target_grids, use_average_angle=False):
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

    travel_times = process_travel_times(
        df, height, width, use_average_angle=use_average_angle
    )

    # intersect列がTrueのセルのmin_costsを0に設定
    starting_indices = []
    cell_idx = idx(target_grids[0], target_grids[1])
    min_costs[cell_idx] = 0.0
    starting_indices.append(cell_idx)

    return min_costs, travel_times, width, height, starting_indices, directions


def run_dijkstra(min_costs, travel_times, width, height, directions):
    def idx(x, y):
        return y * width + x

    visited = np.full_like(min_costs, False, dtype=bool)
    queue = []

    # 開始セルをキューに登録
    visited_count = 0
    for i in range(len(min_costs)):
        if min_costs[i] == 0.0:
            y, x = divmod(i, width)
            heappush(queue, (0.0, x, y))
            visited[idx(x, y)] = True
            visited_count += 1

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

            new_cost = current_cost + t
            if new_cost < min_costs[nidx]:
                min_costs[nidx] = new_cost
                heappush(queue, (new_cost, nx, ny))

    return min_costs


# 2つをまとめて実行する関数
def calculate_minimum_cost_to_geometry(
    df, target_grids, geometry_name="", use_average_angle=False
):
    min_costs, travel_times, width, height, starting_indices, directions = (
        prepare_data_for_dijkstra(df, target_grids, use_average_angle=use_average_angle)
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


# %%
site_grid_coords = get_site_grid_coords(df_obsidian, df_elevation)

min_costs_matrix = np.zeros((len(grid_coords), len(site_grid_coords)))
for i, target_grids in enumerate(tqdm(site_grid_coords)):
    min_costs, travel_times, width, height, starting_indices, directions = (
        prepare_data_for_dijkstra(df_elevation, target_grids, use_average_angle=True)
    )
    min_costs = run_dijkstra(min_costs, travel_times, width, height, directions)

    min_costs = np.where(min_costs == np.inf, np.nan, min_costs)

    with open(
        os.path.join(
            data_dir, "16_tobler_distance_with_coast_50_average", f"distance_siteID_{i}"
        ),
        mode="wb",
    ) as fo:
        pickle.dump(min_costs, fo)
