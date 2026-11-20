"""Tests of absco.validate: TAPE3 and LBLRTM output checks (no binaries or data directory needed)."""

from __future__ import annotations

import struct

import numpy as np
import pytest

from absco import validate as v

N = v.N_LINES


def _panel(lines):
    """Pack ``lines`` into one 250-line TAPE3 record. Each line is a dict with ``vnu``,
    optional ``mol``, and optional ``wv_co2`` = (width, temperature exponent, shift)."""
    vnu = np.zeros(N, dtype="<f8")
    w = np.zeros(37 * N, dtype="<f4")
    mol = np.zeros(N, dtype="<i4")
    flag = np.zeros((N, 7), dtype="<i4")
    dat = np.zeros((N, 21), dtype="<f4")
    for k, ln in enumerate(lines):
        vnu[k] = ln["vnu"]
        mol[k] = ln.get("mol", 101)
        w[0 * N + k] = ln.get("sp", 1e-25)
        w[1 * N + k] = 0.07
        w[2 * N + k] = 1500.0
        w[4 * N + k] = 0.07
        w[5 * N + k] = 0.7
        if "wv_co2" in ln:
            flag[k, v.WV_CO2_FLAG] = 1
            dat[k, v.WV_CO2_DATA] = ln["wv_co2"]
    w[3 * N:4 * N] = mol.view("<f4")
    w[8 * N:15 * N] = flag.reshape(-1).view("<f4")
    w[15 * N:36 * N] = dat.reshape(-1)
    payload = vnu.tobytes() + w.tobytes()
    assert len(payload) == v.PANEL_BYTES
    return struct.pack("<i", len(payload)) + payload + struct.pack("<i", len(payload))


def _header():
    payload = b" LNFL header".ljust(1664, b" ")
    return struct.pack("<i", len(payload)) + payload + struct.pack("<i", len(payload))


def _write_tape3(path, panels):
    path.write_bytes(_header() + b"".join(panels))


def _write_wv_co2(dirpath, rows):
    lines = ["> test file", "> second header line"]
    for vnu, hw, temp, shift in rows:
        lines.append("%13.6f %10.3e %6.4f %10.4f %5.2f %9.6f   1 1 1" % (vnu, 1e-25, hw, 1500.0, temp, shift))
    (dirpath / v.WV_CO2_FILE).write_text("\n".join(lines) + "\n")


GOOD = [{"vnu": 835.249667, "wv_co2": (0.1189, 0.31, 0.0)},
        {"vnu": 835.5532, "wv_co2": (0.1305, 0.28, 0.01005)},
        {"vnu": 836.0087}]            # an H2O line with no CO2 broadening
REF = [(835.249667, 0.1189, 0.31, 0.0), (835.5532, 0.1305, 0.28, 0.01005), (836.0087, 0.12, 0.3, 0.0)]


def test_good_tape3_passes(tmp_path):
    t3 = tmp_path / "TAPE3_H2O_a"
    _write_tape3(t3, [_panel(GOOD)])
    _write_wv_co2(tmp_path, REF)
    v.validate_tape3(t3, brd_dir=tmp_path)


def test_nonfinite_width_is_an_error_and_names_the_line(tmp_path):
    lines = [dict(GOOD[0], wv_co2=(np.nan, 0.31, 0.0)), GOOD[1]]
    t3 = tmp_path / "TAPE3_H2O_a"
    _write_tape3(t3, [_panel(lines)])
    _write_wv_co2(tmp_path, REF)
    with pytest.raises(v.Tape3Error) as err:
        v.validate_tape3(t3, brd_dir=tmp_path)
    msg = str(err.value)
    assert "non-finite" in msg and "835.249667" in msg and "PR 17" in msg


def test_nonfinite_line_strength_is_an_error(tmp_path):
    t3 = tmp_path / "TAPE3_CO2_a"
    _write_tape3(t3, [_panel([{"vnu": 700.0, "mol": 201, "sp": np.inf}])])
    with pytest.raises(v.Tape3Error, match="non-finite"):
        v.validate_tape3(t3)


def test_wrong_but_finite_broadening_is_an_error(tmp_path):
    # the failure of the COMMON block mismatch: finite values that are not the file's
    lines = [dict(GOOD[0], wv_co2=(6.817, 0.1622, 0.71)), GOOD[1]]
    t3 = tmp_path / "TAPE3_H2O_a"
    _write_tape3(t3, [_panel(lines)])
    _write_wv_co2(tmp_path, REF)
    with pytest.raises(v.Tape3Error) as err:
        v.validate_tape3(t3, brd_dir=tmp_path)
    assert "differ from wv_co2_brd_param" in str(err.value) and "835.249667" in str(err.value)


def test_line_missing_from_the_broadening_file_is_an_error(tmp_path):
    t3 = tmp_path / "TAPE3_H2O_a"
    _write_tape3(t3, [_panel([{"vnu": 900.0, "wv_co2": (0.1, 0.3, 0.0)}])])
    _write_wv_co2(tmp_path, REF)
    with pytest.raises(v.Tape3Error, match="not in the file"):
        v.validate_tape3(t3, brd_dir=tmp_path)


def test_co2_broadened_lines_without_the_broadening_file_is_an_error(tmp_path):
    t3 = tmp_path / "TAPE3_H2O_a"
    _write_tape3(t3, [_panel(GOOD)])
    with pytest.raises(v.Tape3Error, match="cannot be checked"):
        v.validate_tape3(t3, brd_dir=tmp_path)


def test_molecule_without_co2_broadening_needs_no_broadening_file(tmp_path):
    t3 = tmp_path / "TAPE3_O3_a"
    _write_tape3(t3, [_panel([{"vnu": 1000.0 + k, "mol": 301} for k in range(5)])])
    v.validate_tape3(t3)                                  # no brd_dir at all


def test_padding_slots_are_ignored(tmp_path):
    # unused slots (vnu <= 0) may hold anything, including NaN
    panel = _panel(GOOD)
    t3 = tmp_path / "TAPE3_H2O_a"
    raw = bytearray(_header() + panel)
    # overwrite the line-strength word of the unused slot 200 with NaN
    struct.pack_into("<f", raw, len(_header()) + 4 + 8 * N + 4 * 200, float("nan"))
    t3.write_bytes(bytes(raw))
    _write_wv_co2(tmp_path, REF)
    v.validate_tape3(t3, brd_dir=tmp_path)


def test_not_a_tape3_is_an_error(tmp_path):
    t3 = tmp_path / "TAPE3_H2O_a"
    t3.write_bytes(b"not a tape3")
    with pytest.raises(v.Tape3Error, match="no line-data records"):
        v.validate_tape3(t3)


def test_od_with_nan_is_an_error_and_gives_the_range():
    wn = np.linspace(850.0, 855.0, 11)
    od = np.ones(11)
    v.check_od_finite(wn, od, "H2O P=986 T=290")
    od[3:6] = np.nan
    with pytest.raises(v.Tape3Error, match="3 of 11 points, from 851.5000 to 852.5000"):
        v.check_od_finite(wn, od, "H2O P=986 T=290")
