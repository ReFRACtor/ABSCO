"""A band that is skipped (the molecule is not active there) must not shift the settings of the bands
after it, and a molecule with nothing to compute must not crash the run.

calcABSCO loops over the TAPE3 files that exist. Inside the loop it used the position in that list
to pick the spectral-degradation kernel and the cross-section flag, which are indexed by the band's
position in the configuration. After a skipped band every later band got the settings of the one
before it, without an error.
"""

from __future__ import annotations

import glob

import numpy as np
import pytest

from absco import compute as absco
from absco import config_assistant as ca
from absco import preprocess as preproc

pytestmark = pytest.mark.integration


def _ini(tmp_path, wn1, wn2, outres, mol):
    cfg = ca.build_config(wn1=wn1, wn2=wn2, outres=outres, molnames=[mol])
    with open(tmp_path / "ABSCO_config.ini", "w") as fh:
        cfg.write(fh)
    return preproc.configure(str(tmp_path / "ABSCO_config.ini"), prompt_user=False)


def test_later_bands_keep_their_own_kernel_after_a_skipped_band(absco_data_dir, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # CO is active at 2000-2250 and 4150-4350 cm-1 but not at 3000; the three bands have three resolutions
    ini = _ini(tmp_path, [2100.0, 3000.0, 4200.0], [2104.0, 3004.0, 4204.0], [1.2e-3, 4.8e-3, 2.4e-3], "co")
    assert ini.doBand["CO"] == [True, False, True]
    k = absco.makeABSCO(ini, "CO", vmrWV=np.nan, vmrO2=np.nan)
    k.lnflT5("CO")
    k.runLNFL()
    k = absco.makeABSCO(ini, "CO", debug=True)
    k.lblT5("CO")
    assert k.calcABSCO("CO") is not False
    k.arrABSCO()
    step = [float(np.median(np.diff(k.freq[a:b + 1]))) for a, b in k.indBands]
    assert len(step) == 2
    assert step[0] == pytest.approx(1.2e-3, rel=1e-3)
    assert step[1] == pytest.approx(2.4e-3, rel=1e-3), "the third band was degraded with the second band's kernel"


def test_a_molecule_with_no_active_band_gives_no_table_and_no_crash(absco_data_dir, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    ini = _ini(tmp_path, [3000.0], [3004.0], [1.2e-3], "co")
    assert ini.doBand["CO"] == [False]
    k = absco.makeABSCO(ini, "CO", vmrWV=np.nan, vmrO2=np.nan)
    k.lnflT5("CO")
    k.runLNFL()
    k = absco.makeABSCO(ini, "CO", debug=True)
    k.lblT5("CO")
    assert k.calcABSCO("CO") is False
    assert not glob.glob(str(tmp_path / "nc_ABSCO" / "*.nc"))
