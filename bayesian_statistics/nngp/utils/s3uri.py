from pathlib import PurePosixPath
from typing import Union
from urllib.parse import urlparse


class S3URI:
    _uri: str  # Type annotation for the _uri attribute

    def __init__(self, uri: Union[str, "S3URI"]):
        """
        S3 URIを表すクラス
        Parameters:
            uri (str or S3URI): S3 URI (例: "s3://bucket/key")
        """
        if isinstance(uri, S3URI):
            self._uri = uri._uri
        else:
            if not self.is_valid_s3_uri(uri):
                raise ValueError(f"Invalid S3 URI: {uri}")
            self._uri = uri

    def __truediv__(self, other: Union[str, "S3URI"]) -> "S3URI":
        """
        S3 URIの結合を行う
        Parameters:
            other (str or S3URI): 結合するS3 URI
        Returns:
            S3URI: 結合されたS3 URI
        """
        base = self._uri.rstrip("/")
        tail = str(other).lstrip("/")
        return S3URI(f"{base}/{tail}")

    def __str__(self):
        """
        S3 URIを文字列として返す
        """
        return self._uri

    def __repr__(self):
        """
        S3 URIの文字列表現を返す
        """
        return f"S3URI({self._uri!r})"

    def __contains__(self, item: str) -> bool:
        """
        Enables usage like: 'foo/bar' in s3uri

        Args:
            item (str): Substring to search for in the full S3 URI

        Returns:
            bool: True if item is found in the S3 URI string representation
        """
        return item in str(self)

    def __eq__(self, other) -> bool:
        """
        S3 URIの等価性を比較する
        Args:
            other (str or S3URI): 比較するS3 URI
        Returns:
            bool:
        """
        if isinstance(other, S3URI):
            return (self.bucket, self.key) == (other.bucket, other.key)
        if isinstance(other, str):
            return str(self) == other
        return NotImplemented

    def __hash__(self) -> int:
        """
        S3 URIのハッシュ値を計算する
        Returns:
            int: S3 URIのハッシュ値
        """
        return hash((self.bucket, self.key))

    @staticmethod
    def is_valid_s3_uri(uri: str) -> bool:
        """
        S3 URIが有効かどうかを確認する
        Parameters:
            uri (str): S3 URI
        Returns:
            bool: 有効な場合はTrue
        """
        parsed = urlparse(uri)
        # バケットだけのS3 URI (例: s3://bucket) も有効と見なす
        return parsed.scheme == "s3" and bool(parsed.netloc)

    @property
    def bucket(self) -> str:
        """
        S3 URIからバケット名を取得する
        Returns:
            str: バケット名
        """
        return self._uri[5:].split("/")[0]

    @property
    def key(self) -> str:
        """
        S3 URIからキーを取得する
        """
        parts = self._uri[5:].split("/")
        return "/".join(parts[1:]) if len(parts) > 1 else ""

    @property
    def path(self) -> PurePosixPath:
        """
        S3 URIをPurePosixPathとして返す
        """
        if not self.key:
            # キーが空の場合は空のPathを作成し、文字列化時に空文字列となるようカスタムクラスを返す
            class EmptyPath(PurePosixPath):
                def __str__(self):
                    return ""

            return EmptyPath(".")
        return PurePosixPath(self.key)

    @property
    def name(self) -> str:
        """
        S3 URIからファイル名を取得する
        """
        if not self.key or self.key.endswith("/"):
            return ""
        return self.path.name

    @property
    def stem(self) -> str:
        """
        S3 URIから拡張子を除いたファイル名を取得する
        """
        if not self.key or self.key.endswith("/"):
            return ""
        return self.path.stem

    @property
    def suffix(self) -> str:
        """
        S3 URIから拡張子を取得する
        """
        if not self.key or self.key.endswith("/"):
            return ""
        return self.path.suffix

    @property
    def parent(self) -> "S3URI":
        """
        S3 URIの親ディレクトリを取得する
        Returns:
            S3URI: 親ディレクトリのS3 URI
        """
        if not self.key:
            # バケットのみの場合は自分自身を返す
            return self

        parent_path = self.path.parent
        if str(parent_path) == ".":
            return S3URI(f"s3://{self.bucket}")
        else:
            return S3URI(f"s3://{self.bucket}/{parent_path}")
