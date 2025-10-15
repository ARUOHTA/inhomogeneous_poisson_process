# ---------- converter.py ----------
from typing import List

import geopandas as gpd
import pandas as pd

from .constants import EPSG_CODE
from .geometry import create_mesh_geometry


def convert_mesh_to_gdf(area_mesh: List[str]) -> gpd.GeoDataFrame:
    df_mesh = pd.DataFrame(area_mesh, columns=["mesh_code"]).astype({"mesh_code": int})
    df_mesh["geometry"] = create_mesh_geometry(df_mesh["mesh_code"])
    return gpd.GeoDataFrame(df_mesh, geometry="geometry", crs=EPSG_CODE)
