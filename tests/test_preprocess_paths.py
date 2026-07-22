"""Unit tests for preprocess config parsing that need no binaries or line file.

These exercise ``configure.readConfig`` in isolation (via ``__new__`` so the
artifact ``file_check`` in ``__init__`` is bypassed).
"""

from __future__ import annotations

import os

from absco import config_assistant as ca
from absco.preprocess import configure


def _read_config(path):
    """Run only readConfig on a config file, returning the populated object."""
    obj = configure.__new__(configure)
    obj.configFile = os.fspath(path)
    obj.prompt_user = False
    obj.readConfig()
    return obj


def _write(path, intdir):
    cfg = ca.build_config([820.0], [910.0], [4e-4], ["h2o"])
    cfg["output"]["intdir"] = intdir
    with open(path, "w") as fh:
        ca.write_config(cfg, fh)


def test_intdir_dot_becomes_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = tmp_path / "c.ini"
    _write(cfg, ".")
    ini = _read_config(cfg)
    assert ini.intdir == os.fspath(tmp_path)


def test_relative_intdir_resolved_absolute(tmp_path, monkeypatch):
    # split_config writes intdir like "./01"; it must resolve to an absolute path
    # (compute.py chdir's into run dirs and then uses intdir-derived paths)
    monkeypatch.chdir(tmp_path)
    cfg = tmp_path / "c.ini"
    _write(cfg, "./01")
    ini = _read_config(cfg)
    assert os.path.isabs(ini.intdir)
    assert ini.intdir == os.fspath(tmp_path / "01")


def test_intdir_expanduser(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = tmp_path / "c.ini"
    _write(cfg, "~/some_absco_intdir_xyz")
    ini = _read_config(cfg)
    assert ini.intdir == os.path.expanduser("~/some_absco_intdir_xyz")
    # clean up the directory readConfig created
    try:
        os.rmdir(ini.intdir)
    except OSError:
        pass
