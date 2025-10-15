import geopandas as gpd
from shapely import wkt


def map_mesh_to_keycode(
    gdf_mesh: gpd.GeoDataFrame, areas: gpd.GeoDataFrame, plane_epsg: str
):
    """
    Map mesh codes to key codes based on spatial intersection with areas.
    Args:
        gdf_mesh (gpd.GeoDataFrame): GeoDataFrame containing mesh codes and geometries.
        areas (gpd.GeoDataFrame): GeoDataFrame containing areas with key codes.
        plane_epsg (str): EPSG code for the plane coordinate system.
    Returns:
        gdf_mesh2keycode (gpd.GeoDataFrame): GeoDataFrame with mesh codes and corresponding key codes.
        dict_mesh2keycode (dict): Dictionary mapping mesh codes to key codes.
    """
    gdf_mesh = gdf_mesh.to_crs(plane_epsg)
    areas = areas.to_crs(plane_epsg)

    gdf_mesh2keycode = gdf_mesh.sjoin(areas, how="left", predicate="intersects")

    gdf_mesh2keycode["geometry_string"] = gdf_mesh2keycode.geometry.astype(str).map(
        wkt.loads
    )

    dict_mesh2keycode = dict(zip(gdf_mesh2keycode.mesh_code, gdf_mesh2keycode.KEY_CODE))

    return gdf_mesh2keycode, dict_mesh2keycode
