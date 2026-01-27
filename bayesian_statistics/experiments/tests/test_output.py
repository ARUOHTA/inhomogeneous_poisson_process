"""Tests for output and progress management."""

import io
from unittest.mock import patch

import pytest


def test_progress_manager_verbosity_levels():
    """verbosity levels が正しく設定されることを確認"""
    from bayesian_statistics.experiments.output import ProgressManager

    pm_quiet = ProgressManager(verbosity=0)
    assert pm_quiet.verbosity == ProgressManager.QUIET

    pm_normal = ProgressManager(verbosity=1)
    assert pm_normal.verbosity == ProgressManager.NORMAL

    pm_verbose = ProgressManager(verbosity=2)
    assert pm_verbose.verbosity == ProgressManager.VERBOSE


def test_progress_manager_section_output():
    """section出力がverbosityに応じて制御されることを確認"""
    from bayesian_statistics.experiments.output import ProgressManager

    # QUIET: 出力なし
    pm_quiet = ProgressManager(verbosity=0)
    with patch("tqdm.tqdm.write") as mock_write:
        pm_quiet.section("Test Section")
        mock_write.assert_not_called()

    # NORMAL: 出力あり
    pm_normal = ProgressManager(verbosity=1)
    with patch("tqdm.tqdm.write") as mock_write:
        pm_normal.section("Test Section")
        assert mock_write.call_count >= 3  # 空行 + '=' + タイトル + '='


def test_progress_manager_nested_depth():
    """nested context managerがdepthを正しく管理することを確認"""
    from bayesian_statistics.experiments.output import ProgressManager

    pm = ProgressManager(verbosity=1)
    assert pm._depth == 0

    with pm.nested():
        assert pm._depth == 1
        with pm.nested():
            assert pm._depth == 2
        assert pm._depth == 1
    assert pm._depth == 0


def test_progress_manager_progress_suppression():
    """nested時にNORMAL verbosityでprogressが抑制されることを確認"""
    from bayesian_statistics.experiments.output import ProgressManager

    pm = ProgressManager(verbosity=1)

    # トップレベル: プログレスバー表示
    assert pm.should_show_progress() is True

    # ネストされたレベル: 抑制
    with pm.nested():
        assert pm.should_show_progress() is False

    # VERBOSE: ネストでも表示
    pm_verbose = ProgressManager(verbosity=2)
    with pm_verbose.nested():
        assert pm_verbose.should_show_progress() is True


def test_get_progress_manager_default():
    """get_progress_managerがデフォルトインスタンスを返すことを確認"""
    from bayesian_statistics.experiments.output import (
        ProgressManager,
        get_progress_manager,
        set_progress_manager,
    )

    # 既存の状態をリセット
    import bayesian_statistics.experiments.output as output_module

    output_module._default_progress = None

    pm = get_progress_manager()
    assert isinstance(pm, ProgressManager)
    assert pm.verbosity == ProgressManager.NORMAL


def test_set_progress_manager():
    """set_progress_managerがグローバルインスタンスを設定することを確認"""
    from bayesian_statistics.experiments.output import (
        ProgressManager,
        get_progress_manager,
        set_progress_manager,
    )

    custom_pm = ProgressManager(verbosity=2)
    set_progress_manager(custom_pm)

    retrieved = get_progress_manager()
    assert retrieved is custom_pm
    assert retrieved.verbosity == 2

    # クリーンアップ: デフォルトに戻す
    import bayesian_statistics.experiments.output as output_module

    output_module._default_progress = None


def test_progress_manager_info_detail():
    """info と detail が verbosity に応じて制御されることを確認"""
    from bayesian_statistics.experiments.output import ProgressManager

    # QUIET: info も detail も出力なし
    pm_quiet = ProgressManager(verbosity=0)
    with patch("tqdm.tqdm.write") as mock_write:
        pm_quiet.info("Info message")
        pm_quiet.detail("Detail message")
        mock_write.assert_not_called()

    # NORMAL: info は出力、detail は出力なし
    pm_normal = ProgressManager(verbosity=1)
    with patch("tqdm.tqdm.write") as mock_write:
        pm_normal.info("Info message")
        assert mock_write.call_count == 1
        mock_write.reset_mock()
        pm_normal.detail("Detail message")
        mock_write.assert_not_called()

    # VERBOSE: 両方出力
    pm_verbose = ProgressManager(verbosity=2)
    with patch("tqdm.tqdm.write") as mock_write:
        pm_verbose.info("Info message")
        pm_verbose.detail("Detail message")
        assert mock_write.call_count == 2
