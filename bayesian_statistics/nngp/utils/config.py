from collections.abc import Mapping
from pathlib import Path

import yaml


class _Namespace:
    """dict を attribute でアクセスできるように薄く包むヘルパ."""

    def __init__(self, mapping: Mapping[str, object]):
        for key, value in mapping.items():
            # ネストしていればさらに _Namespace で包む
            if isinstance(value, Mapping):
                value = _Namespace(value)
            # Python の識別子でないキーはエラーにしたい場合はここで検証
            self.__dict__[key] = value

    # 欠けている属性アクセスでわかりやすいエラーにする
    def __getattr__(self, item):
        raise AttributeError(f"{item!r} は設定に存在しません")

    # dict 互換を最低限サポート（必要なら拡張）
    def __iter__(self):
        return iter(self.__dict__)

    def __getitem__(self, item):
        return self.__dict__[item]

    def items(self):
        return self.__dict__.items()

    def __repr__(self):
        # print() したとき見やすいよう簡易表現
        return f"{self.__class__.__name__}({self.__dict__})"


class Config(_Namespace):
    """YAML から設定を読み込み、任意の深さをドットでアクセスできるクラス."""

    def __init__(self, yaml_path: str | Path):
        yaml_path = Path(yaml_path)
        with yaml_path.open("r", encoding="utf-8") as f:
            raw = yaml.safe_load(f) or {}
        if not isinstance(raw, Mapping):
            raise ValueError("ルートはマッピング型である必要があります")
        super().__init__(raw)


# ------- 使用例 ------- #
# config_path = 'config/hogehoge.yaml'
# config = Config(config_path)
