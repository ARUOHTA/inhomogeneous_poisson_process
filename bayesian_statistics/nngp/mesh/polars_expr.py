# ---------- polars_expr.py ----------
import polars as pl


def create_meshcode_expr(lat_col, lon_col, level=6):
    def split_decimal(col):
        int_part = col.cast(pl.Int64)
        decimal_part = col - int_part
        return int_part, decimal_part

    lat = pl.col(lat_col)
    lon = pl.col(lon_col)

    lat_80, lat_d = split_decimal(lat * 1.5)
    lat_10, lat_d = split_decimal(lat_d * 8)
    lat_1, lat_d = split_decimal(lat_d * 10)
    lat_500, lat_d = split_decimal(lat_d * 2)
    lat_250, lat_d = split_decimal(lat_d * 2)
    lat_125, _ = split_decimal(lat_d * 2)

    lon_80, lon_d = split_decimal(lon - 100)
    lon_10, lon_d = split_decimal(lon_d * 8)
    lon_1, lon_d = split_decimal(lon_d * 10)
    lon_500, lon_d = split_decimal(lon_d * 2)
    lon_250, lon_d = split_decimal(lon_d * 2)
    lon_125, _ = split_decimal(lon_d * 2)

    mesh_code = (lat_80 * 100 + lon_80).cast(pl.Int64).cast(pl.Utf8)
    if level >= 2:
        mesh_code = pl.concat_str(
            [
                mesh_code,
                lat_10.abs().cast(pl.Int64).cast(pl.Utf8),
                lon_10.abs().cast(pl.Int64).cast(pl.Utf8),
            ]
        )
    if level >= 3:
        mesh_code = pl.concat_str(
            [
                mesh_code,
                lat_1.abs().cast(pl.Int64).cast(pl.Utf8),
                lon_1.abs().cast(pl.Int64).cast(pl.Utf8),
            ]
        )
    if level >= 4:
        mesh_code = pl.concat_str(
            [mesh_code, (lat_500 * 2 + lon_500 + 1).cast(pl.Int64).cast(pl.Utf8)]
        )
    if level >= 5:
        mesh_code = pl.concat_str(
            [mesh_code, (lat_250 * 2 + lon_250 + 1).cast(pl.Int64).cast(pl.Utf8)]
        )
    if level >= 6:
        mesh_code = pl.concat_str(
            [mesh_code, (lat_125 * 2 + lon_125 + 1).cast(pl.Int64).cast(pl.Utf8)]
        )
    return mesh_code.cast(pl.Int64)


def add_meshcode(
    df, lat_col="latitude", lon_col="longitude", level=6, output_col="meshcode"
):
    return df.with_columns(
        [create_meshcode_expr(lat_col, lon_col, level).alias(output_col)]
    )


def add_meshcode_lazy(
    lf, lat_col="latitude", lon_col="longitude", level=6, output_col="meshcode"
):
    return lf.with_columns(
        [create_meshcode_expr(lat_col, lon_col, level).alias(output_col)]
    )
