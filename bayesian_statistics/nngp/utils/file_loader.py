from typing import Optional

import geopandas as gpd
import polars as pl
from shapely import wkt

from .s3_manager import S3Manager
from .s3uri import S3URI


def load_geometry_csv(
    file_uri: str | S3URI,
    s3_manager: Optional[S3Manager] = None,
    geometry_col: str = "geometry",
    source_crs: str = "EPSG:4326",
    target_crs: Optional[str] = None,
    rename: Optional[dict] = None,
    usecols: Optional[list] = None,
    dissolve_col: Optional[str] = None,
) -> gpd.GeoDataFrame:
    """
    汎用的なGeoDataFrame読み込み関数（CSV形式）
    target_crsを指定すると、CRS変換を行う。

    Args:
        file_uri (str | S3URI): ファイルパスまたはS3 URI
        s3_manager (S3Manager, optional): S3経由で取得する場合に使用
        geometry_col (str): WKT形式が格納されているカラム名
        source_crs (str): 入力データのCRS（デフォルト: EPSG:4326）
        target_crs (str, optional): 出力データのCRS（例: EPSG:6677）
        rename (dict, optional): カラム名の変換に使用（例: {"type": "fclass"}）
        usecols (list, optional): 読み込むカラムを制限
        dissolve_col (str, optional): dissolveしたい場合のキー（例: "fclass"）

    Returns:
        gpd.GeoDataFrame
    """
    if s3_manager is None:
        file_uri = S3URI(file_uri)
        s3_manager = S3Manager(file_uri.bucket)

    file_path = s3_manager.download(file_uri)
    pl_df = pl.read_csv(file_path, columns=usecols)

    if rename:
        pl_df = pl_df.rename(rename)

    pd_df = pl_df.to_pandas()
    pd_df[geometry_col] = pd_df[geometry_col].map(wkt.loads)

    gdf = gpd.GeoDataFrame(pd_df, geometry=geometry_col, crs=source_crs)

    if target_crs:
        gdf = gdf.to_crs(target_crs)

    if dissolve_col:
        gdf = gdf.dissolve(by=dissolve_col).reset_index()

    return gdf
