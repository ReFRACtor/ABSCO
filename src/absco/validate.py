"""Checks of the binary line file (TAPE3) written by LNFL and of the optical depths
written by LBLRTM.

A bad value in TAPE3 does not stop LNFL or LBLRTM: LBLRTM writes NaN optical depths over
a block of the spectrum around the line, and the NaN ends up in the absorption table.
The checks here turn that into an error at the point it is first visible.

TAPE3 layout (single-precision LNFL build, little endian).  The file is a sequence of
Fortran unformatted sequential records, each ``[4-byte length][payload][4-byte length]``.
Line data are in records of ``PANEL_BYTES`` bytes, ``N_LINES`` lines each, laid out as
the unnamed COMMON block in ``LNFL/src/lnfl.f``::

    VNU3(250)            real*8
    STR3, ALF3, EPP3     real*4 x 250 each
    MOL3                 int*4 x 250   (a real*4 for line coupling lines, IFLG < 0)
    HWHMS, TMPALF, PSHIFT real*4 x 250 each
    IFLG                 int*4 x 250
    ADDFLAG(7,250)       int*4   (Fortran order: the 7 flags of line 1, then line 2, ...)
    ADDDATA(21,250)      real*4  (same order)
    SDEP_DATA(250)       real*4

Slots with ``VNU3 <= 0`` are unused padding.  For H2O lines ``ADDFLAG(2) == 1`` marks a line
with CO2 broadening, and ``ADDDATA(4:6)`` then holds its width, temperature exponent and
shift, copied by ``BRDMATCH`` from ``wv_co2_brd_param``.
"""

import os
import struct

import numpy as np

N_LINES = 250
PANEL_WORDS = 39 * N_LINES
PANEL_BYTES = 4 * PANEL_WORDS

#: ADDFLAG index (0-based) and ADDDATA slice (0-based) of the CO2 broadening of H2O lines
WV_CO2_FLAG = 1
WV_CO2_DATA = slice(3, 6)

#: file name of the CO2 broadening of H2O lines, under the extra broadening directory
WV_CO2_FILE = "wv_co2_brd_param"

#: tolerances for comparing the single-precision TAPE3 values with the text file
RTOL = 1e-5
ATOL = 1e-7

_HINT = ("LNFL (AER-RC/LNFL e79c278, v3.2-30) writes these values when COMMON /WV_CO2_BRD/ is "
         "declared with different sizes; use an LNFL with the fix (AER-RC/LNFL PR 17), rebuild it "
         "with absco-build, and delete the TAPE3 files so they are regenerated.")


class Tape3Error(RuntimeError):
    """A TAPE3 or an LBLRTM output file contains bad values."""


def read_panels(path):
    """Return the line-data panels of a TAPE3 as a list of dicts of numpy arrays."""
    with open(path, "rb") as f:
        raw = f.read()
    panels = []
    i = 0
    while i + 8 <= len(raw):
        n = struct.unpack("<i", raw[i:i + 4])[0]
        if n <= 0 or i + 8 + n > len(raw):
            break
        if struct.unpack("<i", raw[i + 4 + n:i + 8 + n])[0] != n:
            break
        if n == PANEL_BYTES:
            panels.append(_unpack(raw[i + 4:i + 4 + n]))
        i += 8 + n
    if not panels:
        raise Tape3Error("%s: no line-data records of %d bytes found; the file is not a TAPE3 "
                         "from a single-precision LNFL build, or it is empty" % (path, PANEL_BYTES))
    return panels


def _unpack(buf):
    n = N_LINES
    vnu = np.frombuffer(buf, dtype="<f8", count=n)
    w = np.frombuffer(buf, dtype="<f4", offset=8 * n)
    return {
        "vnu": vnu,
        "sp": w[0:n], "alfa": w[n:2 * n], "epp": w[2 * n:3 * n],
        "mol": w[3 * n:4 * n].view("<i4"),
        "hwhm": w[4 * n:5 * n], "tmpalf": w[5 * n:6 * n], "pshift": w[6 * n:7 * n],
        "iflg": w[7 * n:8 * n].view("<i4"),
        "brd_flag": w[8 * n:15 * n].view("<i4").reshape(n, 7),
        "brd_data": w[15 * n:36 * n].reshape(n, 21),
        "speed_dep": w[36 * n:37 * n],
    }


_FLOAT_FIELDS = ("sp", "alfa", "epp", "hwhm", "tmpalf", "pshift", "speed_dep")


def find_nonfinite(panels):
    """List ``(wavenumber, field)`` for every non-finite value in the line data."""
    bad = []
    for p in panels:
        used = p["vnu"] > 0
        for name in _FLOAT_FIELDS:
            for k in np.flatnonzero(used & ~np.isfinite(p[name])):
                bad.append((float(p["vnu"][k]), name))
        for k, j in zip(*np.nonzero(used[:, None] & ~np.isfinite(p["brd_data"]))):
            bad.append((float(p["vnu"][k]), "broadening data %d" % (j + 1)))
        for k in np.flatnonzero(~np.isfinite(p["vnu"])):
            bad.append((float("nan"), "vnu"))
    return bad


def read_wv_co2(path):
    """Read ``wv_co2_brd_param``: arrays of wavenumber, width, temperature exponent and shift."""
    vnu, hw, temp, shift = [], [], [], []
    with open(path, errors="replace") as f:
        for line in f:
            if line[:1] in ">%":
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                v, _, h, _, t, s = (float(x) for x in parts[:6])
            except ValueError:
                continue
            vnu.append(v)
            hw.append(h)
            temp.append(t)
            shift.append(s)
    order = np.argsort(vnu)
    return tuple(np.asarray(a, dtype=float)[order] for a in (vnu, hw, temp, shift))


def find_wv_co2_mismatch(panels, wv_co2_path):
    """List ``(wavenumber, tape3 values, file values)`` for H2O lines whose CO2 broadening in
    TAPE3 differs from ``wv_co2_brd_param``.  The file is read only if the TAPE3 has such lines."""
    flagged = []
    for p in panels:
        sel = (p["vnu"] > 0) & (p["iflg"] >= 0) & (p["mol"] % 100 == 1) & (p["brd_flag"][:, WV_CO2_FLAG] == 1)
        for k in np.flatnonzero(sel):
            with np.errstate(invalid="ignore"):          # a NaN pattern can warn when widened
                flagged.append((float(p["vnu"][k]), p["brd_data"][k, WV_CO2_DATA].astype(float)))
    if not flagged:
        return []
    if not wv_co2_path or not os.path.isfile(wv_co2_path):
        raise Tape3Error("TAPE3 has %d CO2-broadened H2O lines but %s was not found at %r, so their "
                         "broadening cannot be checked" % (len(flagged), WV_CO2_FILE, wv_co2_path))
    v_ref, hw, temp, shift = read_wv_co2(wv_co2_path)
    out = []
    for v, got in flagged:
        i = int(np.searchsorted(v_ref, v - 5e-7))
        if i >= len(v_ref) or abs(v_ref[i] - v) > 5e-7:
            out.append((v, got, None))
            continue
        want = np.array([hw[i], temp[i], shift[i]])
        if not np.all(np.isfinite(got)) or not np.allclose(got, want, rtol=RTOL, atol=ATOL):
            out.append((v, got, want))
    return out


def validate_tape3(path, brd_dir=None):
    """Raise :class:`Tape3Error` if a TAPE3 holds non-finite line data or CO2 broadening values
    that differ from ``wv_co2_brd_param`` (looked for in ``brd_dir``)."""
    panels = read_panels(path)
    problems = []

    bad = find_nonfinite(panels)
    if bad:
        shown = "; ".join("%.6f cm-1 (%s)" % b for b in bad[:5])
        problems.append("%d non-finite value(s) in the line data, first: %s" % (len(bad), shown))

    wv = os.path.join(brd_dir, WV_CO2_FILE) if brd_dir else None
    wrong = find_wv_co2_mismatch(panels, wv)
    if wrong:
        def fmt(w):
            return "not in the file" if w is None else "file (%s)" % ", ".join("%.6g" % x for x in w)
        shown = "; ".join("%.6f cm-1: TAPE3 (%s) vs %s" % (v, ", ".join("%.6g" % x for x in g), fmt(w))
                          for v, g, w in wrong[:3])
        problems.append("%d CO2-broadened H2O line(s) differ from %s, first: %s"
                        % (len(wrong), WV_CO2_FILE, shown))

    if problems:
        raise Tape3Error("%s: %s. %s" % (path, " | ".join(problems), _HINT))


def check_od_finite(wavenumber, od, label):
    """Raise :class:`Tape3Error` if an LBLRTM optical-depth spectrum has non-finite values."""
    od = np.asarray(od)
    bad = ~np.isfinite(od)
    if not bad.any():
        return
    idx = np.flatnonzero(bad)
    wn = np.asarray(wavenumber)
    raise Tape3Error("LBLRTM optical depths for %s are not finite: %d of %d points, from %.4f to %.4f "
                     "cm-1. LBLRTM does not stop on this, but the table would contain NaN. Check the "
                     "TAPE3 (absco.validate.validate_tape3) and the LNFL build."
                     % (label, idx.size, od.size, wn[idx[0]], wn[idx[-1]]))
