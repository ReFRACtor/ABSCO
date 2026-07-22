"""Fast unit tests for the config assistant (no binaries or line file needed)."""

from __future__ import annotations

import configparser
import io

import numpy as np
import pytest

from absco import config_assistant as ca


def test_suggest_lblres_matches_notebook():
    # outres 0.01 cm-1 -> lblres 7.8125e-05 (create_absco_config-tropomi.ipynb)
    lblres = ca.suggest_lblres(0.01)
    assert np.isclose(lblres[0], 7.8125e-05)


@pytest.mark.parametrize("outres", [0.01, 1.2e-3, 1.0e-3, 5.0e-2])
def test_suggest_lblres_power_of_two_and_bounded(outres):
    lblres = ca.suggest_lblres(outres)[0]
    ratio = outres / lblres
    # exact power of 2
    assert np.isclose(np.log2(ratio), round(np.log2(ratio)))
    # at or below the LBLRTM target resolution
    assert lblres <= ca.TARGET_LBLRES + 1e-15


def test_suggest_lblres_fine_outres_not_coarsened():
    # when outres is already finer than the target, lblres must not exceed it
    outres = 1.0e-4
    lblres = ca.suggest_lblres(outres)[0]
    assert lblres <= outres


def test_suggest_lblres_multiband():
    outres = [0.01, 1.2e-3, 1.5e-4]
    lblres = ca.suggest_lblres(outres)
    assert lblres.size == 3
    for o, l in zip(outres, lblres):
        assert np.isclose(np.log2(o / l), round(np.log2(o / l)))


def test_build_config_fills_channels_and_leaves_paths_blank():
    cfg = ca.build_config(
        wn1=[4166.0], wn2=[4358.0], outres=[0.01], molnames=["ch4", "co", "H2O"]
    )
    assert cfg["channels"]["wn1"].strip() == "4166.0000"
    assert cfg["channels"]["wn2"].strip() == "4358.0000"
    assert cfg["channels"]["lblres"].strip().startswith("7.8125")
    # molecules lowercased
    assert cfg["molecules"]["molnames"] == "ch4 co h2o"
    # path fields stay blank for runtime resolution
    assert cfg["runLNFL"]["lnfl_path"].strip() == ""
    assert cfg["data files"]["pfile"].strip() == ""


def test_build_config_rejects_non_power_of_two_lblres():
    with pytest.raises(ValueError):
        ca.build_config(
            wn1=[750.0], wn2=[850.0], outres=[1.0e-3],
            molnames=["h2o"], lblres=[3.0e-4],  # 1e-3/3e-4 is not 2**n
        )


def test_build_config_band_count_mismatch():
    with pytest.raises(ValueError):
        ca.build_config(
            wn1=[750.0, 900.0], wn2=[850.0], outres=[1.0e-3], molnames=["h2o"]
        )


def test_estimate_ram_positive():
    ram = ca.estimate_ram_gb([4166], [4358], [0.01], ["ch4", "co", "h2o"])
    assert ram is not None and ram > 0


def test_write_config_adds_comments_and_defaults():
    cfg = ca.build_config(
        wn1=[4166.0], wn2=[4358.0], outres=[0.01], molnames=["h2o"]
    )
    buf = io.StringIO()
    ca.write_config(cfg, buf)
    text = buf.getvalue()

    # a per-field description comment is present
    assert "; HITRAN molecule names" in text
    assert "; LBLRTM computation resolution" in text
    # blank data-file fields note the packaged default that will be used
    assert "; (blank -> using packaged/default:" in text
    assert "AIRS_P_air.txt" in text

    # still parses as a valid config with the expected values
    parsed = configparser.ConfigParser()
    parsed.read_string(text)
    assert parsed["molecules"]["molnames"] == "h2o"
    assert parsed["channels"]["lblres"].strip().startswith("7.8125")
    assert parsed["data files"]["pfile"].strip() == ""
