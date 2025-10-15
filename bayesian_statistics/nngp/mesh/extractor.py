# ---------- extractor.py ----------
# warningの非表示
import warnings
from multiprocessing import Pool

import geopandas as gpd
from shapely.ops import unary_union
from shapely.strtree import STRtree
from tqdm import tqdm

from .constants import DEFAULT_MESH_LEVEL
from .converter import convert_mesh_to_gdf
from .generator import generate_candidate_mesh

# 緯度経度座標のままbufferを実行すると警告が出るが、意図的なので無視する
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r"Geometry is in a geographic CRS. Results from 'buffer' are likely incorrect.*",
)


def check_intersects(mesh_geom, boundary_union, str_tree):
    possible_intersects = str_tree.query(mesh_geom)
    return mesh_geom.intersects(boundary_union) if any(possible_intersects) else False


def extract_mesh(gdf, candidate_mesh):
    df_mesh = convert_mesh_to_gdf(candidate_mesh)
    str_tree = STRtree(gdf.geometry.buffer(0.001))
    boundary_union = unary_union(gdf.geometry.to_list())
    args = [(geom, boundary_union, str_tree) for geom in df_mesh.geometry]
    with Pool(8) as p:
        intersects = list(
            tqdm(
                p.starmap(check_intersects, args),
                total=len(df_mesh),
                desc="Extracting mesh",
            )
        )
    return df_mesh[intersects].copy()


def area_to_gdf_mesh(areas: gpd.GeoDataFrame, mesh_level: int = DEFAULT_MESH_LEVEL):
    """
    エリア(gdf)から交差するメッシュを抽出し、GeoDataFrameと対応辞書を返す。

    Args:
        areas (gpd.GeoDataFrame): 対象エリア
        mesh_level (int): メッシュレベル（例: 6）

    Returns:
        Tuple[gpd.GeoDataFrame, Dict[int, Polygon]]:
            - メッシュのGeoDataFrame
            - メッシュコードとジオメトリの辞書
    """
    bounds = areas.total_bounds
    print("Generate candidate mesh")
    candidate_mesh = generate_candidate_mesh(bounds, mesh_level)
    gdf_tgt_area_mesh = extract_mesh(areas, candidate_mesh)
    dict_mesh2geo = dict(
        zip(gdf_tgt_area_mesh.mesh_code.astype(int), gdf_tgt_area_mesh.geometry)
    )
    return gdf_tgt_area_mesh, dict_mesh2geo
