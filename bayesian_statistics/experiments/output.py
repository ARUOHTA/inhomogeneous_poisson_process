"""Unified output and progress management for experiments.

This module provides centralized control over console output and progress bars
to reduce noise and provide consistent formatting across all experiment scripts.
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Iterator, TypeVar

from tqdm import tqdm

T = TypeVar("T")


class ProgressManager:
    """Unified progress and output management.

    Controls verbosity levels and provides consistent output formatting.
    Uses tqdm.write() to ensure compatibility with progress bars.

    Verbosity levels:
    - QUIET (0): No output except critical errors
    - NORMAL (1): Section headers and main progress bars only
    - VERBOSE (2): All output including nested progress bars

    Attributes
    ----------
    verbosity : int
        Current verbosity level (0-2).

    Examples
    --------
    >>> pm = ProgressManager(verbosity=1)
    >>> pm.section("Running MMCP model")
    >>> for item in pm.progress(items, desc="Processing"):
    ...     process(item)
    """

    QUIET = 0
    NORMAL = 1
    VERBOSE = 2

    def __init__(self, verbosity: int = 1):
        """Initialize ProgressManager.

        Parameters
        ----------
        verbosity : int
            Verbosity level: 0=quiet, 1=normal, 2=verbose.
        """
        self.verbosity = verbosity
        self._depth = 0

    def section(self, title: str) -> None:
        """Print main section header.

        Uses '=' separators for major sections.

        Parameters
        ----------
        title : str
            Section title to display.
        """
        if self.verbosity >= self.NORMAL:
            tqdm.write("")
            tqdm.write("=" * 70)
            tqdm.write(title)
            tqdm.write("=" * 70)

    def subsection(self, title: str) -> None:
        """Print subsection header.

        Uses '---' formatting for subsections within a major section.

        Parameters
        ----------
        title : str
            Subsection title to display.
        """
        if self.verbosity >= self.NORMAL:
            tqdm.write(f"\n--- {title} ---")

    def info(self, message: str) -> None:
        """Print informational message.

        Displayed at NORMAL verbosity and above.

        Parameters
        ----------
        message : str
            Message to display.
        """
        if self.verbosity >= self.NORMAL:
            tqdm.write(message)

    def detail(self, message: str) -> None:
        """Print detailed message.

        Only displayed at VERBOSE verbosity. Indented with two spaces.

        Parameters
        ----------
        message : str
            Detailed message to display.
        """
        if self.verbosity >= self.VERBOSE:
            tqdm.write(f"  {message}")

    def progress(
        self,
        iterable: Iterator[T],
        desc: str | None = None,
        total: int | None = None,
        **kwargs,
    ) -> Iterator[T]:
        """Create progress bar with nesting control.

        At QUIET verbosity, returns the iterable unchanged (no progress bar).
        At NORMAL verbosity, nested progress bars (depth > 0) are suppressed.
        At VERBOSE verbosity, all progress bars are shown.

        Parameters
        ----------
        iterable : Iterator[T]
            Iterable to wrap with progress bar.
        desc : str, optional
            Description to show in progress bar.
        total : int, optional
            Total number of items (for iterators without __len__).
        **kwargs
            Additional arguments passed to tqdm.

        Returns
        -------
        Iterator[T]
            The iterable, optionally wrapped with tqdm.
        """
        if self.verbosity == self.QUIET:
            return iterable
        if self._depth > 0 and self.verbosity < self.VERBOSE:
            return iterable
        return tqdm(iterable, desc=desc, total=total, **kwargs)

    @contextmanager
    def nested(self):
        """Context manager for nested loops.

        Increments depth counter to suppress inner progress bars
        when verbosity is NORMAL.

        Yields
        ------
        None
        """
        self._depth += 1
        try:
            yield
        finally:
            self._depth -= 1

    def should_show_progress(self) -> bool:
        """Check if progress bars should be shown at current depth.

        Returns
        -------
        bool
            True if progress bars should be displayed.
        """
        if self.verbosity == self.QUIET:
            return False
        if self._depth > 0 and self.verbosity < self.VERBOSE:
            return False
        return True


# Global default instance (can be overridden)
_default_progress: ProgressManager | None = None


def get_progress_manager() -> ProgressManager:
    """Get the global ProgressManager instance.

    Returns
    -------
    ProgressManager
        Global progress manager, created with default settings if needed.
    """
    global _default_progress
    if _default_progress is None:
        _default_progress = ProgressManager()
    return _default_progress


def set_progress_manager(pm: ProgressManager) -> None:
    """Set the global ProgressManager instance.

    Parameters
    ----------
    pm : ProgressManager
        Progress manager to use globally.
    """
    global _default_progress
    _default_progress = pm
