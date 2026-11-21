"""The band edges written to the LBLRTM input must be the configured ones.

Record 1.3 of the TAPE5 holds the band limits in a 10-character field. They were written with
'%10.3e', which keeps four significant digits, so an edge at or above 10,000 cm-1 moved by up to 5 cm-1
(10649 became 10650) and the table covered a different range from the one configured.
"""

from __future__ import annotations

import glob

import numpy as np
import pytest

from absco import compute as absco
from absco import config_assistant as ca
from absco import preprocess as preproc

pytestmark = pytest.mark.integration


# an active molecule for each edge, so lblT5 writes a TAPE5
@pytest.mark.parametrize("mol, wn1, wn2", [("H2O", 10649.0, 10654.0), ("O2", 13089.9996, 13094.9996), ("O2", 29861.0, 29866.0), ("CH4", 4288.0, 4293.0)])
def test_band_edges_are_written_in_full(absco_data_dir, tmp_path, monkeypatch, mol, wn1, wn2):
    monkeypatch.chdir(tmp_path)
    cfg = ca.build_config(wn1=[wn1], wn2=[wn2], outres=[1.2e-3], molnames=[mol.lower()])
    with open(tmp_path / "ABSCO_config.ini", "w") as fh:
        cfg.write(fh)
    ini = preproc.configure(str(tmp_path / "ABSCO_config.ini"), prompt_user=False)
    assert ini.doBand[mol] == [True], "pick an active molecule for this band"
    kobj = absco.makeABSCO(ini, mol, debug=True, vmrWV=1e3, vmrO2=2e5)
    kobj.lblT5(mol)
    f = sorted(glob.glob(str(tmp_path / "LBL_Runs" / "TAPE5_dir" / mol / "*")))[0]
    rec13 = open(f).read().split("\n")[3]                       # record 1.3: V1, V2 in two 10-character fields
    v1, v2 = float(rec13[:10]), float(rec13[10:20])
    assert v1 == pytest.approx(wn1, abs=1e-4)
    assert v2 == pytest.approx(wn2, abs=1e-4)
