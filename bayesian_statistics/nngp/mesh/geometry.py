# ---------- geometry.py ----------

import jismesh.utils as ju
from shapely.geometry import Polygon
from tqdm import tqdm


def create_mesh_geometry(code: list[int] | int) -> list[Polygon]:
    if isinstance(code, int):
        code = [code, code]

    if len(code) < 2:
        code = [code[0], code[0]]

    lat1, lng1 = ju.to_meshpoint(code, 0, 0)
    lat2, lng2 = ju.to_meshpoint(code, 1.000000001, 1)

    return [
        Polygon(
            [
                (lng1[i], lat1[i]),
                (lng1[i], lat2[i]),
                (lng2[i], lat2[i]),
                (lng2[i], lat1[i]),
                (lng1[i], lat1[i]),
            ]
        )
        for i in tqdm(range(len(code)), desc="Creating Mesh Geometry")
    ]


def get_latbin(mesh_level: int) -> float:
    return {
        1: 40 / 60,
        2: 5 / 60,
        3: 30 / 3600,
        4: 15 / 3600,
        5: 7.5 / 3600,
        6: 3.75 / 3600,
    }[int(mesh_level)]


def get_lngbin(mesh_level: int) -> float:
    return {
        1: 1,
        2: 7 / 60 + 30 / 3600,
        3: 45 / 3600,
        4: 22.5 / 3600,
        5: 11.25 / 3600,
        6: 5.625 / 3600,
    }[int(mesh_level)]
