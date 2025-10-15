import geopandas as gpd


def get_buffer_meshes(
    gdf: gpd.GeoDataFrame,
    buffer_meter: float = 500,
    plane_epsg: int = 4326,
    latlon_epsg: int = 4326,
) -> gpd.GeoDataFrame:
    """
    指定のGeoDataFrameにバッファを追加し、メッシュコードを付与する
    Args:
        gdf (gpd.GeoDataFrame): 対象のGeoDataFrame
        buffer_meter (float): バッファの距離（メートル）
        plane_epsg (int): 平面座標系のEPSGコード
        latlon_epsg (int): 緯度経度座標系のEPSGコード
    Returns:
        gpd.GeoDataFrame: バッファを追加したGeoDataFrame
    """
    buffer_mesh = gdf.copy()
    buffer_mesh = buffer_mesh.to_crs(epsg=plane_epsg)
    buffer_mesh["geometry"] = buffer_mesh.geometry.buffer(buffer_meter)
    buffer_mesh = buffer_mesh.to_crs(epsg=latlon_epsg)

    return buffer_mesh
