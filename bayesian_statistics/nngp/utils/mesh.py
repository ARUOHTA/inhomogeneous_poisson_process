import functools

import jismesh.utils as ju
import numpy as np
from shapely.geometry import Polygon


class MeshCalculator:
    def __init__(self):
        self.lat_factors = np.array(
            [
                2 / 3,
                1 / 12,
                1 / 120,
                1 / 240,
                1 / 480,
                1 / 960,
                1 / 1920,
                1 / 3840,
                1 / 11520,
            ]
        )
        self.lon_factors = np.array(
            [1, 1 / 8, 1 / 80, 1 / 160, 1 / 320, 1 / 640, 1 / 1280, 1 / 2560, 1 / 7680]
        )
        self.mesh_value_lookup = {
            1: np.array([0, 0]),
            2: np.array([0, 1]),
            3: np.array([1, 0]),
            4: np.array([1, 1]),
        }
        self.mesh_value_9digit_lookup = {
            1: np.array([0, 0]),
            2: np.array([0, 1]),
            3: np.array([0, 2]),
            4: np.array([1, 0]),
            5: np.array([1, 1]),
            6: np.array([1, 2]),
            7: np.array([2, 0]),
            8: np.array([2, 1]),
            9: np.array([2, 2]),
        }

    @functools.lru_cache(maxsize=1)
    def _unit_lat(self, level):
        return self.lat_factors[level - 1]

    @functools.lru_cache(maxsize=1)
    def _unit_lon(self, level):
        return self.lon_factors[level - 1]

    def coordinates_to_mesh_code(self, latitudes, longitudes, target_level=6):
        lat = latitudes
        lon = longitudes - 100

        def calc_mesh_part(value, unit, offset=0):
            return np.floor((value - offset) / unit).astype(int)

        mesh_codes = []
        offsets_lat, offsets_lon = 0, 0

        for level in range(
            1, target_level + 1
        ):  # Use self.level to limit the mesh level
            unit_lat = self._unit_lat(level)
            unit_lon = self._unit_lon(level)
            mesh_lat = calc_mesh_part(lat, unit_lat, offset=offsets_lat)
            mesh_lon = calc_mesh_part(lon, unit_lon, offset=offsets_lon)

            if level < 4:
                offsets_lat += mesh_lat * unit_lat
                offsets_lon += mesh_lon * unit_lon
                mesh_codes.extend([mesh_lat, mesh_lon])
            elif level < 9:
                mesh_value = mesh_lat * 2 + mesh_lon + 1
                offsets_lat += (mesh_value // 3) * unit_lat
                offsets_lon += ((mesh_value + 1) % 2) * unit_lon
                mesh_codes.append(mesh_value)
            else:
                mesh_value = mesh_lat * 3 + mesh_lon + 1
                mesh_codes.append(mesh_value)

        return ["".join(map(str, row)) for row in np.column_stack(mesh_codes)]

    @functools.lru_cache(maxsize=1024)
    def mesh_code_to_coordinates(self, mesh_code, lat_=1.0, lon_=1.0):
        mesh_code = str(mesh_code)
        target_level = len(mesh_code) - 5
        mesh_values = np.array([int(x) for x in mesh_code])

        lat = mesh_values[0:2].dot([10, 1]) * self.lat_factors[0]
        lon = mesh_values[2:4].dot([10, 1]) + 100

        lat += np.sum(mesh_values[4:7:2] * self.lat_factors[1:3])
        lon += np.sum(mesh_values[5:8:2] * self.lon_factors[1:3])

        for i, value in enumerate(mesh_values[8:-1]):
            lat_lon_offset = self.mesh_value_lookup[value]
            lat += lat_lon_offset[0] * self.lat_factors[3 + i]
            lon += lat_lon_offset[1] * self.lon_factors[3 + i]

        if target_level == 9:
            lat_lon_offset = self.mesh_value_9digit_lookup[mesh_values[-1]]
        else:
            lat_lon_offset = self.mesh_value_lookup[mesh_values[-1]]

        lat += lat_lon_offset[0] * self.lat_factors[target_level - 1]
        lon += lat_lon_offset[1] * self.lon_factors[target_level - 1]

        lat += lat_ * self.lat_factors[target_level - 1]
        lon += lon_ * self.lon_factors[target_level - 1]

        return lat, lon

    def mesh_code_to_polygon(self, mesh_code):
        lower_left = self.mesh_code_to_coordinates(mesh_code, lat_=0, lon_=0)
        upper_right = self.mesh_code_to_coordinates(mesh_code, lat_=1, lon_=1)

        return Polygon(
            [
                [lower_left[1], lower_left[0]],
                [upper_right[1], lower_left[0]],
                [upper_right[1], upper_right[0]],
                [lower_left[1], upper_right[0]],
                [lower_left[1], lower_left[0]],
            ]
        )

    def split_mesh_code_to_level(self, mesh_code, target_level):
        assert len(mesh_code) >= 5, (
            "Input mesh_code must be at least 5 characters long."
        )
        assert 1 <= target_level <= 9, "target_level must be between 1 and 9."

        def recursive_split(current_mesh_code, current_level):
            if current_level == target_level:
                return [current_mesh_code]

            next_level_mesh_codes = []
            if current_level < 8:
                # For levels 1 to 8, split into 4 parts
                for i in range(1, 5):
                    next_level_mesh_codes.extend(
                        recursive_split(current_mesh_code + str(i), current_level + 1)
                    )
            else:
                # For level 8 to 9, split into 9 parts
                for i in range(1, 10):
                    next_level_mesh_codes.extend(
                        recursive_split(current_mesh_code + str(i), current_level + 1)
                    )

            return next_level_mesh_codes

        return recursive_split(mesh_code, len(mesh_code) - 5)
