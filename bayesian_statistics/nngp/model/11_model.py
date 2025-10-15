"""Backward-compatibility shim for code that imports ``11_model``.

The actual implementation now lives in :mod:`multinomial_model`.
"""

from .multinomial_model import *  # noqa: F401,F403

__all__ = [name for name in globals() if not name.startswith("_")]
