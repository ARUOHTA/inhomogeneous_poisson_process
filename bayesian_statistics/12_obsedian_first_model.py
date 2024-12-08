import os
from typing import Tuple

import numpy as np
import polars as pl

# %% [markdown]
# 地点ごとの産地構成比を計算

# %%
def create_base_grid(
    df: pl.DataFrame,
    grid_size: int = 100
) -> tuple[np.ndarray, np.ndarray]:
    """
    解析領域全体をカバーする等間隔グリッドを作成
    
    Parameters
    ----------
    df : pl.DataFrame
        入力データフレーム（緯度・経度を含む）
    grid_size : int
        グリッドの分割数
        
    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        緯度・経度のメッシュグリッド (lat_mesh, lon_mesh)
    """
    lon_min, lon_max = df['経度'].min(), df['経度'].max()
    lat_min, lat_max = df['緯度'].min(), df['緯度'].max()
    
    lon_grid = np.linspace(lon_min, lon_max, grid_size)
    lat_grid = np.linspace(lat_min, lat_max, grid_size)
    lon_mesh, lat_mesh = np.meshgrid(lon_grid, lat_grid)
    
    return lon_mesh, lat_mesh

def create_elevation_grid(df_elevation: pl.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """
    df_elevationのすべての地点を含むグリッドを作成
    
    Parameters
    ----------
    df_elevation : pl.DataFrame
        地形データフレーム（x: 経度, y: 緯度を含む）
        
    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        緯度・経度のメッシュグリッド (lat_mesh, lon_mesh)
    """
    # ユニークな経度・緯度の値を取得（ソートされた形で）
    unique_lons = df_elevation['x'].unique().sort()
    unique_lats = df_elevation['y'].unique().sort()
    
    # メッシュグリッドを作成
    lon_mesh, lat_mesh = np.meshgrid(
        unique_lons.to_numpy(),  # 経度の一意な値
        unique_lats.to_numpy()   # 緯度の一意な値
    )
    
    return lon_mesh, lat_mesh

def preprocess_data(
    df: pl.DataFrame,
    target_period: int,
    target_origin: str
) -> Tuple[np.ndarray, np.ndarray]:
    """
    解析用のデータを前処理
    
    Parameters
    ----------
    df : pl.DataFrame
        入力データフレーム
    target_period : int
        対象時期
    target_origin : str
        対象産地カテゴリ
        
    Returns
    -------
    counts : np.ndarray
        各遺跡での出土数（インデックスが遺跡ID）
    target_counts : np.ndarray
        対象産地の出土数（インデックスが遺跡ID）
    """
    # 全遺跡IDのリストを取得
    max_site_id = df['遺跡ID'].max()
    
    # 対象時期のデータのみ抽出
    period_df = df.filter(pl.col('時期') == target_period)
    
    # 全体のカウント
    counts = (
        period_df
        .group_by('遺跡ID')
        .agg([pl.len().alias('count')])
        .join(
            pl.DataFrame({
                '遺跡ID': np.arange(max_site_id + 1)
            }),
            on='遺跡ID',
            how='right'
        )
        .fill_null(0)
        .sort('遺跡ID')['count']
        .to_numpy()
    )
    
    # 対象産地のカウント
    target_counts = (
        period_df
        .filter(pl.col('産地カテゴリ') == target_origin)
        .group_by('遺跡ID')
        .agg([pl.len().alias('count')])
        .join(
            pl.DataFrame({
                '遺跡ID': np.arange(max_site_id + 1)
            }),
            on='遺跡ID',
            how='right'
        )
        .fill_null(0)
        .sort('遺跡ID')['count']
        .to_numpy()
    )
    
    return counts, target_counts

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
        df_result
        .select([
            (pl.col("遺跡ID")), 
            (pl.col('緯度') * np.pi / 180).alias('lat_rad'),
            (pl.col('経度') * np.pi / 180).alias('lon_rad')
        ])
        .unique(subset=["遺跡ID"])
        .sort("遺跡ID")
    )
    
    # 座標と出土数を numpy 配列に変換
    site_coords = np.column_stack([
        coords['lat_rad'].to_numpy(),
        coords['lon_rad'].to_numpy()
    ])

    return site_coords

# %%
def calculate_weights_matrix(
    grid_coords: np.ndarray, #(N, 2)
    site_coords: np.ndarray, #(M, 2)
    sigma: float,
) -> np.ndarray:
    """
    重み行列を計算
    """
    R = 6371  # 地球の半径(km)
    
    # 通常の距離計算
    dlat = grid_coords[:, np.newaxis, 0] - site_coords[np.newaxis, :, 0]
    dlon = grid_coords[:, np.newaxis, 1] - site_coords[np.newaxis, :, 1]
    
    a = (np.sin(dlat/2)**2 + 
         np.cos(grid_coords[:, np.newaxis, 0]) * 
         np.cos(site_coords[np.newaxis, :, 0]) * 
         np.sin(dlon/2)**2)
    
    distances = 2 * R * np.arcsin(np.sqrt(a))
    
    # 重みの初期計算
    weights = np.exp(-0.5 * (distances**2) / (sigma**2)) / (2 * np.pi * sigma**2)

    return weights


# %%
def update_weights_matrix(weights, grid_coords, df_elevation, lon_mesh, lat_mesh):
    # 地形マスクの作成
    land_points = df_elevation.select([
        'x',
        'y',
        pl.col('is_sea').cast(pl.Boolean)
    ]).to_numpy()
    
    lons_1d = lon_mesh[0, :]
    lats_1d = lat_mesh[:, 0]
    land_mask = np.full(lon_mesh.shape, False)
    
    x_indices = np.searchsorted(lons_1d, land_points[:, 0])
    y_indices = np.searchsorted(lats_1d, land_points[:, 1])
    valid_points = (
        (x_indices > 0) & 
        (x_indices < len(lons_1d)) & 
        (y_indices > 0) & 
        (y_indices < len(lats_1d))
    )
    is_sea = land_points[valid_points, 2].astype(bool)
    land_mask[y_indices[valid_points], x_indices[valid_points]] = ~is_sea
    
    # grid_coordsの各点について、対応するland_maskの値を取得
    grid_lons = grid_coords[:, 1] * 180/np.pi  # ラジアンから度に変換
    grid_lats = grid_coords[:, 0] * 180/np.pi
    
    grid_x_indices = np.searchsorted(lons_1d, grid_lons)
    grid_y_indices = np.searchsorted(lats_1d, grid_lats)
    
    # インデックスが有効範囲内にあることを確認
    valid_grid_points = (
        (grid_x_indices > 0) & 
        (grid_x_indices < len(lons_1d)) & 
        (grid_y_indices > 0) & 
        (grid_y_indices < len(lats_1d))
    )
    
    # 海上の点の重みを0に設定
    grid_is_land = np.zeros(len(grid_coords), dtype=bool)
    grid_is_land[valid_grid_points] = land_mask[
        grid_y_indices[valid_grid_points],
        grid_x_indices[valid_grid_points]
    ]
    
    # 海上の点からの重みをすべて0に
    weights = weights * grid_is_land[:, np.newaxis]
    
    return weights

def calculate_ratios(
    weights: np.ndarray,
    counts: np.ndarray,
    target_counts: np.ndarray
) -> np.ndarray:
    """
    重み付き比率を計算
    
    Parameters
    ----------
    weights : np.ndarray
        重み行列
    counts : np.ndarray
        各遺跡での出土数
    target_counts : np.ndarray
        対象産地の出土数
        
    Returns
    -------
    ratios : np.ndarray
        各グリッド点での重み付き比率
    """
    # 重み付き合計を計算
    weighted_total = np.sum(weights * counts, axis=1)
    weighted_target = np.sum(weights * target_counts, axis=1)
    
    # 比率計算（0除算を防ぐ）
    ratios = np.where(
        weighted_total > 0,
        weighted_target / weighted_total,
        0
    )
    
    return ratios


# %%
def calculate_obsedian_weight(
    df: pl.DataFrame,
    df_elevation: pl.DataFrame,
    sigma: int
) -> None:
    """
    黒曜石の分布をプロットする関数
    
    Parameters
    ----------
    df : pl.DataFrame
        入力データフレーム
    df_elevation : pl.DataFrame
        標高データフレーム
    sigma : int
        ガウス分布の標準偏差
    target_period : int
        対象時期
    target_origin : str
        対象産地カテゴリ
    grid_size : int
        グリッドのサイズ
    figsize : tuple
        プロットのサイズ
    """
    # グリッドの作成
    # lon_mesh, lat_mesh = create_base_grid(df, grid_size)

    lon_mesh, lat_mesh = create_elevation_grid(df_elevation)

    print(lon_mesh.shape, lat_mesh.shape)

    # データの前処理
    site_coords = create_site_coords(df)

    # グリッド座標の準備
    grid_coords = np.column_stack([
        lat_mesh.ravel() * np.pi / 180,
        lon_mesh.ravel() * np.pi / 180
    ])

    print("creating weights matrix...")
    print(f"grid_coords: {grid_coords.shape}, site_coords: {site_coords.shape}")
    # 陸地のみの重み行列の計算
    weights = calculate_weights_matrix(
        grid_coords, site_coords, sigma
    )

    print(f"weights: {weights.shape}")

    print("updating weights matrix...")
    # 重みの更新
    weights_updated = update_weights_matrix(weights, grid_coords, df_elevation, lon_mesh, lat_mesh)

    return lon_mesh, lat_mesh, weights_updated

def calculate_obsedian_ratio(
        df: pl.DataFrame,
        target_period: int,
        target_origin: str,
        weights: np.ndarray,
        lon_mesh: np.ndarray,
        lat_mesh: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:

    # ここからtarget_period, target_originに依存する処理
    counts, target_counts = preprocess_data(
        df, target_period, target_origin
    )

    # 重み付き比率の計算
    ratios = calculate_ratios(weights, counts, target_counts)

    ratio_mesh = ratios.reshape(lon_mesh.shape)

    return lon_mesh, lat_mesh, ratio_mesh
    

def calculate_site_ratios_fast(
    df: pl.DataFrame,
    sigma: float,
    target_period: int,
    target_origin: str
) -> pl.DataFrame:
    """
    各遺跡地点での比率を高速計算
    
    Parameters
    ----------
    df : pl.DataFrame
        入力データフレーム
    radius_km : float
        計算する円の半径(km)
    target_period : int
        対象時期(0-4)
    target_origin : str
        対象産地カテゴリ
        
    Returns
    -------
    pl.DataFrame
        緯度、経度、比率を含むDataFrame
    """
    # 遺跡の一意な地点を取得
    unique_sites = df.unique(subset=['遺跡ID']).sort('遺跡ID')
    
    # データの前処理（全遺跡の出土データ）
    site_coords = create_site_coords(df)
    
    # 計算対象の遺跡の座標をラジアンに変換
    target_coords = np.column_stack([
        unique_sites['緯度'].to_numpy() * np.pi / 180,
        unique_sites['経度'].to_numpy() * np.pi / 180
    ])
    
    # 距離行列の計算
    # 陸地のみの重み行列の計算
    weights = calculate_weights_matrix(
            target_coords, site_coords, sigma
    )

    # 重みの更新
    weights_updated = update_weights_matrix(weights, target_coords, df_elevation, lon_mesh, lat_mesh)

    # ここからtarget_period, target_originに依存する処理
    counts, target_counts = preprocess_data(
        df, target_period, target_origin
    )
    
    # 比率の計算
    ratios = calculate_ratios(weights_updated, counts, target_counts)

    # いったんすべての結果をDataFrameに変換
    result_df = pl.DataFrame({
        '遺跡ID': unique_sites['遺跡ID'],
        f"比率_{target_period}_{target_origin}": ratios
    })
    
    return result_df

# %%

data_dir = "/home/ohta/dev/bayesian_statistics/data/"

df_elevation = pl.read_csv(os.path.join(data_dir, "gdf_elevation_with_costs.csv"))
df_result = pl.read_csv(os.path.join(data_dir, "11_gdf_obsedian.csv"))
df_sites = pl.read_csv(os.path.join(data_dir, "11_gdf_sites.csv"))

time_period_name = {
    0: "早期・早々期",
    1: "前期",
    2: "中期",
    3: "後期",
    4: "晩期"
}

origin_order = ["神津島", "信州", "箱根", "高原山", "その他"]

sigma = 14
sigma_for_sites = 0.1

# %%
#calculate_obsedian_weight
lon_mesh, lat_mesh, weights = calculate_obsedian_weight(df_result, df_elevation, sigma)

# %%
# list(time_period_name.keys()), origin_order[:-1]のすべての組み合わせを作成

ratio_df = pl.DataFrame({
    'x': lon_mesh.ravel(),
    'y': lat_mesh.ravel()
})

ratio_sites_df = pl.DataFrame({
    '遺跡ID': df_sites['遺跡ID']
})
    

for target_period in time_period_name.keys():
    for target_origin in origin_order[:-1]:

        print(f"target_period: {target_period}, target_origin: {target_origin}")

        lon_mesh, lat_mesh, ratio_mesh = calculate_obsedian_ratio(
            df = df_result, 
            target_period = target_period,
            target_origin = target_origin,
            weights = weights,
            lon_mesh = lon_mesh,
            lat_mesh = lat_mesh
        )

        
        ratio_df = ratio_df.join(
            pl.DataFrame({
                'x': lon_mesh.ravel(),
                'y': lat_mesh.ravel(),
                f"ratio_{target_period}_{target_origin}": ratio_mesh.ravel()
            }),
            on=["x", "y"]
        )

        ratio_sites_df = ratio_sites_df.join(
            calculate_site_ratios_fast(
                df = df_result,
                sigma = sigma_for_sites,
                target_period = target_period,
                target_origin = target_origin
            ),
            on="遺跡ID"
        )

df_elevation = df_elevation.join(
    ratio_df,
    on=["x", "y"]
)
df_sites = df_sites.join(
    ratio_sites_df,
    on="遺跡ID"
)

# %%
df_elevation.write_csv(os.path.join(data_dir, "12_gdf_elevation_with_ratio.csv"))
df_sites.write_csv(os.path.join(data_dir, "12_gdf_sites_with_ratio.csv"))
