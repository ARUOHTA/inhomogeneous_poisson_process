import os
import re
from pathlib import Path
from typing import List, Tuple, Union

import boto3
import pandas as pd
import polars as pl
from botocore.exceptions import ClientError

from .s3uri import S3URI

# 型エイリアス: S3 URI の入力として許容する形式
S3INPUT = Union[str, Path, S3URI, Tuple[str, str]]


class S3Manager:
    def __init__(self, bucket: str, cache_dir: str = "/tmp") -> None:
        """
        S3Manager は S3 バケットへのアクセスを管理するクラス。

        Args:
            bucket (str): 操作対象の S3 バケット名
            cache_dir (str): ダウンロードファイルを保存するローカルキャッシュディレクトリ
        """
        self.bucket = bucket
        self.cache_dir = cache_dir
        self.client = boto3.client("s3")

    def _to_s3uri(self, s3_input: S3INPUT) -> S3URI:
        """
        入力（文字列・パス・タプル）を S3URI に変換する内部ユーティリティ。

        Args:
            s3_input (S3INPUT): S3のキーやURI情報

        Returns:
            S3URI: 標準化された S3 URI オブジェクト
        """
        if isinstance(s3_input, S3URI):
            return s3_input
        elif isinstance(s3_input, (str, Path)):
            s3_input = str(s3_input)
            if s3_input.startswith("s3://"):
                return S3URI(s3_input)
            else:
                return S3URI(f"s3://{self.bucket}/{s3_input.lstrip('/')}")
        elif isinstance(s3_input, tuple) and len(s3_input) == 2:
            bucket, key = s3_input
            return S3URI(f"s3://{bucket}/{key.lstrip('/')}")
        else:
            raise TypeError("Unsupported input type for S3 URI")

    def exists(self, s3_input: S3INPUT) -> bool:
        """
        指定したS3 URIが存在するか確認。

        Args:
            s3_input (S3INPUT): 対象のS3 URIやキー

        Returns:
            bool: 存在する場合はTrue
        """
        uri = self._to_s3uri(s3_input)
        try:
            self.client.head_object(Bucket=uri.bucket, Key=uri.key)
            return True
        except ClientError as error:
            if error.response["Error"]["Code"] == "404":
                return False
            raise error

    def download(self, s3_input: S3INPUT) -> str:
        """
        S3からローカルにファイルをダウンロード（キャッシュ機能付き）。

        Args:
            s3_input (S3INPUT): 対象のS3 URI

        Returns:
            str: ローカルに保存されたファイルパス
        """
        uri = self._to_s3uri(s3_input)
        local_path = Path(self.cache_dir) / uri.key
        local_path.parent.mkdir(parents=True, exist_ok=True)

        if local_path.exists():
            s3_last_modified = self.client.head_object(Bucket=uri.bucket, Key=uri.key)[
                "LastModified"
            ].timestamp()
            local_last_modified = local_path.stat().st_mtime
            if local_last_modified >= s3_last_modified:
                return str(local_path)

        self.client.download_file(uri.bucket, uri.key, str(local_path))
        return str(local_path)

    def upload(self, local_path: str, s3_input: S3INPUT) -> None:
        """
        ローカルファイルをS3にアップロード。

        Args:
            local_path (str): ローカルファイルのパス
            s3_input (S3INPUT): アップロード先のS3 URI
        """
        uri = self._to_s3uri(s3_input)
        self.client.upload_file(local_path, uri.bucket, uri.key)

    def upload_dir(self, local_dir: str, prefix: Union[str, Path, S3URI] = "") -> int:
        """
        指定ディレクトリ以下のすべてのファイルをS3にアップロード。

        Args:
            local_dir (str): アップロードするローカルディレクトリ
            prefix (Union[str, Path, S3URI]): S3側でのプレフィックス（キーの先頭）

        Returns:
            int: アップロードしたファイル数
        """
        base_path = Path(local_dir).resolve()
        prefix_uri = self._to_s3uri(prefix)

        file_count = 0

        # ここでは、Path.rglob()をモックできるように、1回の呼び出しでファイル一覧を取得
        for file_path in base_path.rglob("*"):
            # is_file()がモックされるため、直接使用
            if file_path.is_file():
                # relative_to()がモックされるため、直接使用
                relative_path = file_path.relative_to(base_path)

                # S3キーを構築
                if prefix_uri.key:
                    s3_key = f"{prefix_uri.key}/{str(relative_path)}"
                else:
                    s3_key = str(relative_path)

                # パスのセパレータを統一（Windowsの場合に必要）
                s3_key = s3_key.replace(os.sep, "/")

                # ファイルパス文字列を取得（__str__メソッドはモック可能）
                local_file_path = str(file_path)

                # アップロードを実行
                self.client.upload_file(local_file_path, prefix_uri.bucket, s3_key)
                file_count += 1

        return file_count

    def list_keys(self, prefix: str | S3URI = "", pattern: str = r".*") -> List[str]:
        """
        S3内の指定prefix以下のオブジェクト一覧を正規表現でフィルタして返す。

        Args:
            prefix (str): キーのプレフィックス
            pattern (str): 正規表現でフィルタするパターン

        Returns:
            List[str]: 条件に一致するキー一覧
        """
        if isinstance(prefix, S3URI):
            prefix = prefix.key

        paginator = self.client.get_paginator("list_objects_v2")
        matched = []
        for page in paginator.paginate(Bucket=self.bucket, Prefix=prefix):
            for obj in page.get("Contents", []):
                key = obj["Key"]
                if re.search(pattern, key):
                    matched.append(key)
        return matched

    def to_csv(
        self,
        data: Union[pd.DataFrame, pl.DataFrame],
        s3_input: S3INPUT,
        index: bool = False,
        encoding: str = "utf-8",
        is_delete: bool = True,
    ) -> None:
        """
        任意のDataFrameをローカルに保存し、S3にCSV形式でアップロード。

        Args:
            data (Union[pd.DataFrame, pl.DataFrame]): アップロードするデータ
            s3_input (S3INPUT): アップロード先のS3 URI
            index (bool): pandas の場合、indexを含めるか
            encoding (str): 出力ファイルのエンコーディング
            is_delete (bool): アップロード後にローカルファイルを削除するかどうか
        """
        uri = self._to_s3uri(s3_input)
        local_path = Path(self.cache_dir) / uri.key
        local_path.parent.mkdir(parents=True, exist_ok=True)

        if isinstance(data, pd.DataFrame):
            data.to_csv(local_path, index=index, encoding=encoding)
        elif isinstance(data, pl.DataFrame):
            data.write_csv(str(local_path))
        else:
            raise TypeError(
                "Unsupported data type. Must be pandas or polars DataFrame."
            )

        self.upload(str(local_path), uri)

        if is_delete:
            os.remove(local_path)
