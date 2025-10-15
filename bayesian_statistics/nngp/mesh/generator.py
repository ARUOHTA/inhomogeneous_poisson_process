# ---------- generator.py ----------
import jismesh.utils as ju

from .geometry import get_latbin, get_lngbin


def generate_candidate_mesh(bounds, mesh_level):
    minx, miny, maxx, maxy = bounds
    lat_bin = get_latbin(mesh_level)
    lng_bin = get_lngbin(mesh_level)
    minx, miny = minx - lng_bin / 2, miny - lat_bin / 2
    maxx, maxy = maxx + lng_bin / 2, maxy + lat_bin / 2
    lat_range = int((maxy - miny) / lat_bin) + 2
    lng_range = int((maxx - minx) / lng_bin) + 2
    coords = list(
        {
            (miny + lat_bin * i, minx + lng_bin * j)
            for i in range(-2, lat_range + 1)
            for j in range(-2, lng_range + 1)
        }
    )
    lats, lngs = zip(*coords)
    return ju.to_meshcode(lats, lngs, mesh_level)
