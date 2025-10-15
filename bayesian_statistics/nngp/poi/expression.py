import ast

import polars as pl
from src.mesh.polars_expr import create_meshcode_expr


def parse_fsq_category_labels(x: str, level: int) -> str:
    """
    fsq_category_labels カラムから指定レベルのカテゴリを抽出するパーサ。

    Args:
        x (str): fsq_category_labels の文字列（例: "[[\"Food\", \"Bakery\"]]"）
        level (int): 抽出する階層レベル

    Returns:
        str: 指定レベルのカテゴリ、または 'Unknown'
    """
    try:
        labels = ast.literal_eval(x)
        # " -> ' に変換
        return labels[0][min(level, len(labels[0]) - 1)]
    except Exception:
        return "Unknown"


def expr_parse_fsq_category(level: int) -> pl.Expr:
    """
    fsq_category_labels から指定レベルのカテゴリを抽出する Polars Expression

    Args:
        level (int): 抽出するカテゴリ階層（例: 0: top-level, 1: sub-level）

    Returns:
        pl.Expr: カテゴリ名（pl.Utf8）
    """
    return (
        pl.col("fsq_category_labels")
        .fill_null("[['Unknown']]")
        .map_elements(
            lambda x: parse_fsq_category_labels(x, level), return_dtype=pl.Utf8
        )
    )


def expr_category_full(level=1) -> pl.Expr:
    """
    category_first と category_second を結合したカテゴリ名を生成する expression
    Returns:
        pl.Expr: 例: "Food_Bakery"
    """
    if level == 1:
        return expr_parse_fsq_category(0).alias("category")
    elif level == 2:
        return (
            expr_parse_fsq_category(0).alias("category_first")
            + "_"
            + expr_parse_fsq_category(1).alias("category_second")
        ).alias("category")
    else:
        raise ValueError("level must be 1 or 2")


def expr_meshcode(
    lat_col: str = "latitude", lon_col: str = "longitude", mesh_level: int = 3
) -> pl.Expr:
    """
    緯度経度からメッシュコードを生成する expression

    Args:
        lat_col (str): 緯度カラム名
        lon_col (str): 経度カラム名
        mesh_level (int): メッシュコードレベル（例: 1〜5）

    Returns:
        pl.Expr: メッシュコード列
    """
    return (
        create_meshcode_expr(lat_col, lon_col, mesh_level)
        .cast(pl.Int64)
        .alias("mesh_code")
    )
