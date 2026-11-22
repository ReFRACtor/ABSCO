"""A band whose two edges are equal has no width and produces a table of NaN, so it is refused."""

from __future__ import annotations

import pytest

from absco import config_assistant as ca
from absco import preprocess as preproc

pytestmark = pytest.mark.integration


def _config(tmp_path, wn1, wn2):
    cfg = ca.build_config(wn1=[wn1], wn2=[wn2], outres=[1.2e-3], molnames=["ch4"])
    path = tmp_path / "ABSCO_config.ini"
    with open(path, "w") as fh:
        cfg.write(fh)
    return str(path)


def test_a_band_with_equal_edges_is_refused(absco_data_dir, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit, match="wn1 equal to wn2"):
        preproc.configure(_config(tmp_path, 9208.9872, 9208.9872), prompt_user=False)


def test_a_band_with_a_width_is_accepted(absco_data_dir, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    ini = preproc.configure(_config(tmp_path, 9208.9872, 9209.9872), prompt_user=False)
    assert ini.channels["wn2"][0] > ini.channels["wn1"][0]
