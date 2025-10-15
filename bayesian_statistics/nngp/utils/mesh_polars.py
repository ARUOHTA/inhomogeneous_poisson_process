from typing import Tuple

import polars as pl
import polars_st as pst


def mesh_code_to_coordinates(
    mesh_code_col: str, lat_offset: float = 1.0, lon_offset: float = 1.0
) -> Tuple[float, float]:
    """
    メッシュコードのカラム名を受け取り、座標を返す。
    """
    lat_factors = [
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
    lon_factors = [
        1,
        1 / 8,
        1 / 80,
        1 / 160,
        1 / 320,
        1 / 640,
        1 / 1280,
        1 / 2560,
        1 / 7680,
    ]

    mesh_code_clean = (
        pl.when(
            pl.col(mesh_code_col).is_null()
            | (pl.col(mesh_code_col).cast(pl.Utf8) == "")
            | (pl.col(mesh_code_col).cast(pl.Utf8).str.len_chars() < 4)
        )
        .then(pl.lit("0000"))
        .otherwise(pl.col(mesh_code_col).cast(pl.Utf8))
    )

    mesh_level = mesh_code_clean.str.len_chars() - 5

    def get_digit(pos):
        return (
            pl.when(mesh_code_clean.str.len_chars() > pos)
            .then(
                mesh_code_clean.str.slice(pos, 1)
                .str.to_integer(strict=False)
                .fill_null(0)
            )
            .otherwise(0)
        )

    # Calculate base coordinates (level 1)
    lat = (get_digit(0) * 10 + get_digit(1)) * lat_factors[0]
    lon = (get_digit(2) * 10 + get_digit(3)) + 100

    # Add level 2
    lat = lat + pl.when(mesh_level >= 2).then(get_digit(4) * lat_factors[1]).otherwise(
        0.0
    )
    lon = lon + pl.when(mesh_level >= 2).then(get_digit(5) * lon_factors[1]).otherwise(
        0.0
    )

    # Add level 3
    lat = lat + pl.when(mesh_level >= 3).then(get_digit(6) * lat_factors[2]).otherwise(
        0.0
    )
    lon = lon + pl.when(mesh_level >= 3).then(get_digit(7) * lon_factors[2]).otherwise(
        0.0
    )

    # Add levels 4-8
    for level in range(4, 9):
        digit_idx = 4 + level
        digit = get_digit(digit_idx)

        lat_add = pl.when(digit.is_in([3, 4])).then(1.0).otherwise(0.0)
        lon_add = pl.when(digit.is_in([2, 4])).then(1.0).otherwise(0.0)

        lat = lat + pl.when(mesh_level >= level).then(
            lat_add * lat_factors[level - 1]
        ).otherwise(0.0)
        lon = lon + pl.when(mesh_level >= level).then(
            lon_add * lon_factors[level - 1]
        ).otherwise(0.0)

    # Add level 9
    digit_9 = get_digit(13)
    lat_add_9 = (
        pl.when(digit_9.is_in([1, 2, 3]))
        .then(0.0)
        .when(digit_9.is_in([4, 5, 6]))
        .then(1.0)
        .when(digit_9.is_in([7, 8, 9]))
        .then(2.0)
        .otherwise(0.0)
    )
    lon_add_9 = (
        pl.when(digit_9.is_in([1, 4, 7]))
        .then(0.0)
        .when(digit_9.is_in([2, 5, 8]))
        .then(1.0)
        .when(digit_9.is_in([3, 6, 9]))
        .then(2.0)
        .otherwise(0.0)
    )

    lat = lat + pl.when(mesh_level == 9).then(lat_add_9 * lat_factors[8]).otherwise(0.0)
    lon = lon + pl.when(mesh_level == 9).then(lon_add_9 * lon_factors[8]).otherwise(0.0)

    # Apply offset based on mesh level
    offset_factor_lat = (
        pl.when(mesh_level == 1)
        .then(lat_factors[0])
        .when(mesh_level == 2)
        .then(lat_factors[1])
        .when(mesh_level == 3)
        .then(lat_factors[2])
        .when(mesh_level == 4)
        .then(lat_factors[3])
        .when(mesh_level == 5)
        .then(lat_factors[4])
        .when(mesh_level == 6)
        .then(lat_factors[5])
        .when(mesh_level == 7)
        .then(lat_factors[6])
        .when(mesh_level == 8)
        .then(lat_factors[7])
        .when(mesh_level == 9)
        .then(lat_factors[8])
        .otherwise(0.0)
    )

    offset_factor_lon = (
        pl.when(mesh_level == 1)
        .then(lon_factors[0])
        .when(mesh_level == 2)
        .then(lon_factors[1])
        .when(mesh_level == 3)
        .then(lon_factors[2])
        .when(mesh_level == 4)
        .then(lon_factors[3])
        .when(mesh_level == 5)
        .then(lon_factors[4])
        .when(mesh_level == 6)
        .then(lon_factors[5])
        .when(mesh_level == 7)
        .then(lon_factors[6])
        .when(mesh_level == 8)
        .then(lon_factors[7])
        .when(mesh_level == 9)
        .then(lon_factors[8])
        .otherwise(0.0)
    )

    final_lat = lat + lat_offset * offset_factor_lat
    final_lon = lon + lon_offset * offset_factor_lon

    return final_lat, final_lon


def mesh_code_to_polygon(mesh_code_col: str) -> pst.GeoExpr:
    """
    メッシュコードが含まれているカラム名から、ポリゴンを作成するpolars_st.GeoExprを返す。
    polars-stというライブラリが必要
    """
    try:
        import polars_st as pst
    except ImportError:
        raise ImportError("polars-st is required for polygon creation.")

    # メッシュの左下角と右上角の座標を取得
    ll_lat, ll_lon = mesh_code_to_coordinates(mesh_code_col, 0.0, 0.0)
    ur_lat, ur_lon = mesh_code_to_coordinates(mesh_code_col, 1.0, 1.0)

    # polars-stのrectangle関数を使用してGeoExprを作成
    return pst.rectangle(ll_lon, ll_lat, ur_lon, ur_lat)
