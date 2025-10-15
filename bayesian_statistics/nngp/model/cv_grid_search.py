"""Cross-validation is not yet implemented for the multinomial NNGP model."""

from __future__ import annotations

from typing import Dict


def compute_cv_scores(*args, **kwargs) -> Dict[str, float]:
    raise NotImplementedError(
        "Cross-validation for the multinomial NNGP model has not been implemented."
    )


__all__ = ["compute_cv_scores"]
