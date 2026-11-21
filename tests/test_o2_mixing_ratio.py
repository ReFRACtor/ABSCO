"""The O2 mixing ratios of an O2 table must reach the LBLRTM input.

absco-generate runs O2 at two fixed mixing ratios (ini.vmrO2) so the table has an O2 axis. Before this
was tested, the value only entered the TAPE5 file name: the two TAPE5 files for the two ratios were
byte-identical and the two halves of the table were the same.
"""

from __future__ import annotations

import glob
import os

import numpy as np
import pytest

from absco import compute as absco
from absco import config_assistant as ca
from absco import preprocess as preproc

pytestmark = pytest.mark.integration


def _tape5_o2_vmrs(path):
    """O2 mixing ratio [ppmv] of every level in a TAPE5 (molecule 7 of record 3.6)."""
    lines = open(path).read().split("\n")
    vmr = []
    for i, l in enumerate(lines):
        if l.rstrip().endswith("AA"):
            nums = []
            j = i + 1
            while len(nums) < 47:                              # 47 molecules, 8 per line
                nums += [float(x) for x in lines[j].split()]
                j += 1
            vmr.append(nums[6])
    return np.array(vmr), lines


def test_o2_tape5_carries_the_o2_mixing_ratio(absco_data_dir, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = ca.build_config(wn1=[13100.0], wn2=[13110.0], outres=[1.2e-3], molnames=["o2"])
    with open(tmp_path / "ABSCO_config.ini", "w") as fh:
        cfg.write(fh)
    ini = preproc.configure(str(tmp_path / "ABSCO_config.ini"), prompt_user=False)
    made = {}
    for o2 in ini.vmrO2:
        kobj = absco.makeABSCO(ini, "O2", debug=True, vmrWV=ini.wv_vmr[0], vmrO2=o2)
        kobj.lblT5("O2")
        f = sorted(glob.glob(str(tmp_path / "LBL_Runs" / "TAPE5_dir" / "O2" / ("*vmrO2%06.0f" % o2))))
        assert f, "no TAPE5 written for O2 = %g ppmv" % o2
        made[o2] = f[0]
    lo, hi = sorted(made)
    v_lo, lines_lo = _tape5_o2_vmrs(made[lo])
    v_hi, lines_hi = _tape5_o2_vmrs(made[hi])
    assert np.allclose(v_lo, lo) and np.allclose(v_hi, hi), (v_lo[:3], v_hi[:3])
    assert lines_lo != lines_hi                                # the files differ, not only in name
    # water vapour is still set from its own axis
    assert os.path.basename(made[lo]).count("_vmrWV") == 1
