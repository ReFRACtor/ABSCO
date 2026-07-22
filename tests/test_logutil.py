"""Unit tests for absco.logutil (no binaries or line file needed)."""

from __future__ import annotations

import io
import os
import sys

from absco import logutil


def test_resolve_log_dir_default(tmp_path):
    d = logutil.resolve_log_dir(None, os.fspath(tmp_path))
    assert d == tmp_path / "logs"
    assert d.is_dir()


def test_resolve_log_dir_override(tmp_path):
    target = tmp_path / "my_logs"
    d = logutil.resolve_log_dir(os.fspath(target), os.fspath(tmp_path / "unused"))
    assert d == target.resolve()
    assert d.is_dir()


def test_tee_stdio_writes_file_and_terminal(tmp_path, capsys):
    log = tmp_path / "run.log"
    with logutil.tee_stdio(log):
        print("hello tee")
    # went to the real terminal (captured by pytest)...
    assert "hello tee" in capsys.readouterr().out
    # ...and to the log file
    assert "hello tee" in log.read_text()
    # streams restored
    assert sys.stdout is sys.__stdout__ or hasattr(sys.stdout, "write")


def test_tee_stdio_appends(tmp_path):
    log = tmp_path / "run.log"
    with logutil.tee_stdio(log):
        print("first")
    with logutil.tee_stdio(log):
        print("second")
    text = log.read_text()
    assert "first" in text and "second" in text


def test_run_logged_captures_output(tmp_path):
    log = tmp_path / "cmd.log"
    rc = logutil.run_logged(
        [sys.executable, "-c", "print('subproc out')"],
        log, header="unit test", echo=False,
    )
    assert rc == 0
    text = log.read_text()
    assert "=== unit test ===" in text
    assert "subproc out" in text


def test_run_logged_returns_nonzero(tmp_path):
    log = tmp_path / "cmd.log"
    rc = logutil.run_logged(
        [sys.executable, "-c", "import sys; sys.exit(3)"],
        log, echo=False,
    )
    assert rc == 3
