import polars as pl
from src.poi.expression import (
    expr_category_full,
    expr_meshcode,
)

usecols = [
    "latitude",
    "longitude",
    "fsq_category_labels",
    "locality",
]

dtype = {
    "latitude": "object",
    "longitude": "object",
    "fsq_category_labels": "str",
    "locality": "str",
}


def load_and_filter_poi_data(
    poi_path: str,
    city_names: list[str],
    lat_min: float,
    lat_max: float,
    lon_min: float,
    lon_max: float,
) -> pl.LazyFrame:
    """
    POIデータを読み込み、都市名に基づいてフィルタする
    """
    df = pl.scan_csv(
        poi_path,
        separator="\t",
        has_header=True,
        ignore_errors=True,
        quote_char="",  # データの中の"が原因でエラーになることがあるため、quote_charを空にする
    )
    df = df.filter(
        pl.col("locality").is_in(city_names)
        & pl.col("latitude").is_not_null()
        & pl.col("longitude").is_not_null()
    )

    if lat_min and lat_max and lon_min and lon_max:
        df = df.filter(
            (pl.col("latitude") >= lat_min)
            & (pl.col("latitude") <= lat_max)
            & (pl.col("longitude") >= lon_min)
            & (pl.col("longitude") <= lon_max)
        )

    return df


def assign_poi_features(
    df: pl.LazyFrame, level: int, mesh_level: int, mesh_set: set[str]
) -> pl.LazyFrame:
    """
    POIデータにカテゴリ・メッシュコードを付与し、メッシュでフィルタする
    """
    return df.with_columns(
        [
            expr_category_full(level=level),
            expr_meshcode("latitude", "longitude", mesh_level),
        ]
    ).filter(pl.col("mesh_code").is_in(mesh_set))


def agg_poi(df) -> pl.DataFrame:
    """
    メッシュ×カテゴリでPOI数を集計する（pivot形式）

    Parameters
    ----------
    df : pl.LazyFrame または pl.DataFrame または pandas.DataFrame
        POIデータフレーム

    Returns
    -------
    pl.DataFrame
        集計結果（pivot形式）
    """
    # pandas.DataFrameをpolars.DataFrameに変換
    if hasattr(df, "to_polars"):
        df = df.to_polars()

    # polars.DataFrameをpolars.LazyFrameに変換
    if isinstance(df, pl.DataFrame):
        df = df.lazy()

    # 集計処理
    return (
        df.group_by(["mesh_code", "category"])
        .len()
        .collect()
        .pivot("category", index="mesh_code", values="len")
        .fill_null(0)
    )
