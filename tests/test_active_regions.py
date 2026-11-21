"""The file of active regions (data/active_regions.csv) and the tool that checks it against the line list.

A band is processed for a molecule only if it overlaps one of the molecule's regions. The regions were a
hand-compiled table in preprocess.py; when they were narrower than the line list a band was skipped and the
generator stopped with an AttributeError. The coverage test below fails if the line list has strong lines
that no region covers, so the same gap cannot come back, for example with a new line file.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from absco import active_regions, line_summary, paths
from absco import config_assistant as ca
from absco import preprocess as preproc
from absco.cli import active_regions as cli

REGION_FILE = paths.data_file(active_regions.FILE)
SUMMARY_FILE = paths.data_file(line_summary.FILE)


def test_the_packaged_file_loads_and_is_consistent():
    regions = active_regions.load()
    assert regions
    by = active_regions.by_molecule()                      # raises if regions of one molecule overlap
    assert set(by) == {r.molecule for r in regions}
    assert all(r.basis and r.source in active_regions.SOURCES for r in regions)
    # a continuum entry says where the continuum comes from
    assert all(r.reference for r in regions if "continuum" in r.basis)


def test_every_allowed_molecule_has_a_region():
    allowed = ["H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2", "NO", "SO2", "NO2", "NH3", "HNO3", "OCS", "H2CO", "N2",
               "HCN", "C2H2", "HCOOH", "C2H4", "CH3OH", "CCL4", "CF4", "F11", "F12", "F22", "ISOP", "PAN", "HDO", "BRO"]
    assert set(allowed) <= set(active_regions.by_molecule())


def test_the_file_is_packaged_in_the_wheel():
    text = (Path(__file__).resolve().parent.parent / "pyproject.toml").read_text()
    assert '"**/*.csv"' in text                            # package-data picks up the region file and the summary


@pytest.mark.parametrize("bad, message", [
    ("X,100,50,lines,2018,,\n", "not below"),
    ("X,1,2,magic,2018,,\n", "unknown basis"),
    ("X,1,2,lines,someone,,\n", "unknown source"),
])
def test_a_malformed_file_is_refused(tmp_path, bad, message):
    f = tmp_path / "r.csv"
    f.write_text("molecule,wn1,wn2,basis,source,reference,note\n" + bad)
    with pytest.raises(ValueError, match=message):
        active_regions.load(str(f))


def test_overlapping_regions_of_one_molecule_are_refused(tmp_path):
    f = tmp_path / "r.csv"
    f.write_text("molecule,wn1,wn2,basis,source,reference,note\nX,1,10,lines,2018,,\nX,5,20,lines,2018,,\n")
    with pytest.raises(ValueError, match="overlap"):
        active_regions.by_molecule(str(f))


# bands the line list has lines in but the 2018 table did not cover
@pytest.mark.parametrize("mol, wn1, wn2", [
    ("nh3", 4288.0, 4488.0), ("nh3", 6537.0, 6592.0), ("nh3", 10000.0, 10300.0),
    ("ch4", 9659.0, 9859.0), ("ch4", 11000.0, 11200.0),
    ("h2o", 25100.0, 25300.0), ("h2o", 41000.0, 41900.0),
    ("o2", 13300.0, 13500.0), ("o2", 14189.0, 14389.0), ("o2", 34100.0, 34300.0),
    ("h2co", 1700.0, 1780.0),
])
def test_bands_the_old_table_skipped_are_processed(tmp_path, monkeypatch, mol, wn1, wn2):
    monkeypatch.chdir(tmp_path)
    cfg = ca.build_config(wn1=[wn1], wn2=[wn2], outres=[3.84e-2], molnames=[mol])
    with open(tmp_path / "ABSCO_config.ini", "w") as fh:
        cfg.write(fh)
    ini = preproc.configure(str(tmp_path / "ABSCO_config.ini"), prompt_user=False)
    assert ini.doBand[mol.upper()][0] is True


@pytest.mark.parametrize("mol, wn1, wn2", [("nh3", 10500.0, 10600.0), ("ch4", 11600.0, 11800.0), ("o2", 2000.0, 3000.0)])
def test_bands_beyond_the_line_list_are_still_skipped(tmp_path, monkeypatch, mol, wn1, wn2):
    monkeypatch.chdir(tmp_path)
    cfg = ca.build_config(wn1=[wn1], wn2=[wn2], outres=[3.84e-2], molnames=[mol])
    with open(tmp_path / "ABSCO_config.ini", "w") as fh:
        cfg.write(fh)
    ini = preproc.configure(str(tmp_path / "ABSCO_config.ini"), prompt_user=False)
    assert ini.doBand[mol.upper()][0] is False


def test_every_strong_bin_of_the_line_list_is_covered():
    """Peak optical depth of the standard column above 1e-3 between 400 and 50000 cm-1 needs a region."""
    summary = line_summary.read(SUMMARY_FILE)
    column = line_summary.columns(paths.data_file("VMR/USS_AIRS_profile.csv"))
    missing = cli.missing_bins(summary, column, active_regions.load())
    text = "; ".join("%s %s" % (m, ", ".join("%d-%d" % tuple(r) for r in cli.ranges(b))) for m, b in sorted(missing.items()))
    assert not missing, "bins with peak tau above %g that no region covers: %s" % (cli.THRESHOLD, text)


def test_the_check_finds_the_gaps_of_the_old_table(tmp_path):
    """The same check on the 2018 table (the NH3, CH4, H2O and O2 regions it had) must report those gaps."""
    old = {"H2O": [(100, 25000)], "CH4": [(100, 9100)], "NH3": [(750, 1200), (1450, 1800), (3200, 3600)],
           "O2": [(0, 150), (1350, 1850), (6200, 6500), (7500, 8500), (9100, 11000), (12990, 13224), (15000, 29870), (36000, 50000)]}
    regions = [r for r in active_regions.load() if r.molecule not in old]
    regions += [active_regions.Region(m, a, b, ("lines",), "2018") for m, v in old.items() for a, b in v]
    summary = line_summary.read(SUMMARY_FILE)
    column = line_summary.columns(paths.data_file("VMR/USS_AIRS_profile.csv"))
    missing = cli.missing_bins(summary, column, regions, threshold=1e-4)
    assert {"CH4", "H2O", "NH3", "O2"} <= set(missing)


def test_the_check_uses_a_peak_not_a_band_mean():
    """One strong narrow line in a 10 cm-1 bin must count even though the bin mean is far below the threshold."""
    summary = {"X": [(100.0, 1, 1e-3 * np.pi * 0.07)]}     # S/(pi gamma) with gamma 0.07, column 1 -> peak tau 1e-3
    assert line_summary.peak_tau(summary, {"X": 1.0})["X"][0][1] == pytest.approx(1e-3 * np.pi * 0.07)
    lines = {"X": (np.array([100.5]), np.array([1e-20]), np.array([0.07]))}
    (b, n, v), = line_summary.summarize(lines)["X"]
    assert (b, n) == (100.0, 1) and v == pytest.approx(1e-20 / (np.pi * 0.07))


def test_a_new_line_file_changes_the_summary_only_if_the_lines_changed():
    """When the line file is present the packaged summary must equal one rebuilt from it."""
    tape1 = paths.line_file_paths()["tape1_path"]
    if not os.path.isfile(tape1):
        pytest.skip("AER line file not in the data directory")
    rebuilt = line_summary.summarize(line_summary.read_lines(tape1))
    packaged = line_summary.read(SUMMARY_FILE)
    assert set(rebuilt) == set(packaged)
    for mol in rebuilt:
        a = {b: (n, v) for b, n, v in rebuilt[mol]}
        p = {b: (n, v) for b, n, v in packaged[mol]}
        assert a.keys() == p.keys(), mol
        assert all(a[b][0] == p[b][0] and a[b][1] == pytest.approx(p[b][1], rel=1e-3) for b in a), mol
