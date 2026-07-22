"""End-to-end integration test for the ABSCO generation pipeline.

Exercises the real LNFL -> LBLRTM -> netCDF flow through the refactored package on
a narrow H2O band, using debug mode (3 pressure levels) to keep it fast. Requires
the compiled binaries and AER line file; see the ``absco_data_dir`` fixture, which
skips this test when they are absent (unless ABSCO_TEST_BUILD=1).
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


def _write_config(path):
    """Write a minimal single-band H2O config via the config assistant."""
    cfg = ca.build_config(
        wn1=[850.0], wn2=[855.0], outres=[1.2e-3], molnames=["h2o"]
    )
    with open(path, "w") as fh:
        cfg.write(fh)


def test_generate_pipeline_produces_readable_table(absco_data_dir, tmp_path, monkeypatch):
    # run entirely inside a scratch working directory
    monkeypatch.chdir(tmp_path)

    config_file = tmp_path / "ABSCO_config.ini"
    _write_config(config_file)

    ini = preproc.configure(str(config_file), prompt_user=False)
    assert ini.molnames == ["H2O"]

    # --- LNFL: build TAPE3 (mirrors absco-generate -lnfl) ---
    for mol in ini.molnames:
        kobj = absco.makeABSCO(ini, mol, vmrWV=np.nan, vmrO2=np.nan)
        kobj.lnflT5(mol)
        kobj.runLNFL()

    tape3 = glob.glob(str(tmp_path / "TAPE3_dir" / "H2O" / "TAPE3_*"))
    assert tape3, "LNFL did not produce a TAPE3 file"

    # --- LBLRTM + netCDF for the water-vapor-affected molecule (debug=3 P levels) ---
    mol = "H2O"
    vmr_objs = []
    for ppm in ini.wv_vmr:
        kobj = absco.makeABSCO(ini, mol, debug=True, vmrWV=ppm)
        kobj.lblT5(mol)
        kobj.calcABSCO(mol)
        kobj.arrABSCO()
        vmr_objs.append(kobj)
    combined = absco.combineVMR(vmr_objs)
    combined.makeNC(mol)

    # --- temp/intermediate dirs must be under the working directory ---
    for sub in ("LNFL_Runs", "LBL_Runs", "TAPE3_dir", "nc_ABSCO"):
        assert (tmp_path / sub).is_dir(), f"{sub} not created under cwd"

    # --- output netCDF exists and has the expected structure ---
    nc_files = glob.glob(str(tmp_path / "nc_ABSCO" / "*.nc"))
    assert len(nc_files) == 1, nc_files
    out_nc = nc_files[0]

    import xarray as xr

    with xr.open_dataset(out_nc) as ds:
        for var in ("Cross_Section", "Spectral_Grid", "T_level", "P_level", "H2O_VMR"):
            assert var in ds.variables, f"missing variable {var}"
        # Cross_Section dims: nfreq x ntemp x nlay x nvmr
        assert ds["Cross_Section"].ndim == 4
        nfreq = ds.sizes["nfreq"]
        assert ds.sizes["nvmr"] == len(ini.wv_vmr)
        assert ds["Spectral_Grid"].size == nfreq
        # the coefficients should contain finite, positive values somewhere
        xsec = np.asarray(ds["Cross_Section"])
        assert np.isfinite(xsec).any()
        assert np.nanmax(xsec) > 0


def test_read_tables_returns_coefficient(absco_data_dir, tmp_path, monkeypatch):
    """absco-read locates a coefficient in a freshly generated table."""
    monkeypatch.chdir(tmp_path)
    config_file = tmp_path / "ABSCO_config.ini"
    _write_config(config_file)

    ini = preproc.configure(str(config_file), prompt_user=False)
    mol = "H2O"
    for m in ini.molnames:
        kobj = absco.makeABSCO(ini, m, vmrWV=np.nan, vmrO2=np.nan)
        kobj.lnflT5(m)
        kobj.runLNFL()
    vmr_objs = []
    for ppm in ini.wv_vmr:
        kobj = absco.makeABSCO(ini, mol, debug=True, vmrWV=ppm)
        kobj.lblT5(mol)
        kobj.calcABSCO(mol)
        kobj.arrABSCO()
        vmr_objs.append(kobj)
    absco.combineVMR(vmr_objs).makeNC(mol)

    out_nc = glob.glob(str(tmp_path / "nc_ABSCO" / "*.nc"))[0]

    from absco.cli.read_tables import testABSCO

    # In debug mode only a few pressure levels are populated; pick a level whose
    # temperatures are finite so the reader's nanargmin has a valid slice.
    import xarray as xr
    with xr.open_dataset(out_nc) as ds:
        p_levels = np.asarray(ds["P_level"])
        t_levels = np.asarray(ds["T_level"])
        finite_levels = [
            i for i in range(p_levels.size - 1)  # last level is TOA (disallowed)
            if np.isfinite(p_levels[i]) and np.isfinite(t_levels[i]).any()
        ]
        assert finite_levels, "no populated pressure level found"
        iP = finite_levels[0]
        p_level = float(p_levels[iP])
        t_level = float(t_levels[iP][np.isfinite(t_levels[iP])][0])
        wv = float(np.asarray(ds["H2O_VMR"])[0])
        wn = float(np.asarray(ds["Spectral_Grid"])[ds.sizes["nfreq"] // 2])

    reader = testABSCO({
        "ncFile": out_nc,
        "in_pressure": p_level,
        "in_temp": t_level,
        "in_spectral": [wn, "cm-1"],
        "in_h2o": wv,
        "in_o2": 190000.0,
        "tolerance": 0.05,
    })
    reader.valueLocate()
    # readABSCO prints the value; here we just assert indices were resolved
    assert reader.iP >= 0 and reader.iWN >= 0
