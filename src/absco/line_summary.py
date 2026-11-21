"""Where the AER line list has absorption, per molecule, as a table small enough to ship.

The measure is the peak optical depth of the strongest line in a bin of the standard column:

    tau_peak = S * N / (pi * gamma)

with ``S`` the line strength (cm-1 / (molecule cm-2) at 296 K), ``N`` the molecule's column in the
standard profile (molecules cm-2) and ``gamma`` the line's own air-broadened half width at the
surface (cm-1).  A bin-mean optical depth would average a narrow strong line away, by a factor of
about 45 at the surface and about 500 where the width is Doppler limited.  The summary keeps, per
molecule and per bin of ``BIN`` cm-1, the largest ``S / (pi gamma)`` and the line count; the column
is applied when the summary is read, so a different profile does not need the line file.

The AER line file (``aer_v_3.9``, the ``tape1_path`` of the configuration) is a HITRAN-format
file: ``I2 I1 F12.6 D10.3 E10.3 F5.4 F5.4 ...``, that is molecule, isotopologue, wavenumber,
strength, Einstein A, air width, self width.  Lines are read as LNFL reads them.
"""

from __future__ import annotations

import csv
import os

import numpy as np

BIN = 10.0                      # cm-1
FILE = "active_regions_summary.csv"
G = 9.80665                     # m s-2
M_AIR = 4.8097e-26              # kg, mean dry air molecule

#: HITRAN molecule numbers (the first two columns of a line) to the names used here
MOLECULE_NUMBER = {
    1: "H2O", 2: "CO2", 3: "O3", 4: "N2O", 5: "CO", 6: "CH4", 7: "O2", 8: "NO", 9: "SO2", 10: "NO2",
    11: "NH3", 12: "HNO3", 19: "OCS", 20: "H2CO", 22: "N2", 23: "HCN", 26: "C2H2", 32: "HCOOH",
    38: "C2H4", 39: "CH3OH", 42: "CF4",
}

#: profile columns that hold each molecule (the column names in the packaged VMR profile)
PROFILE_COLUMN = {
    "SO2": "SO2_HI", "NO2": "NO2_HI", "HNO3": "HNO3_HI", "CF4": "CF4_HI",
}


def parse_line(text: str):
    """Molecule number, wavenumber [cm-1], strength and air width [cm-1] of one line, or None."""
    if len(text) < 45 or text[0] in ">%":
        return None
    try:
        mol = int(text[0:2])
        wn = float(text[3:15])
        s = float(text[15:25].replace("D", "E"))
        gamma = float(text[35:40])
    except ValueError:
        return None
    return mol, wn, s, gamma


def read_lines(path: str, molecules=None):
    """``{molecule: (wn, S, gamma)}`` arrays from an AER HITRAN-format file."""
    keep = None if molecules is None else {n for n, m in MOLECULE_NUMBER.items() if m in molecules}
    cols: dict[int, list] = {}
    with open(path, errors="replace") as fh:
        for text in fh:
            r = parse_line(text)
            if r is None or (keep is not None and r[0] not in keep) or r[0] not in MOLECULE_NUMBER:
                continue
            cols.setdefault(r[0], []).append(r[1:])
    return {MOLECULE_NUMBER[n]: tuple(np.array(c).T) for n, c in cols.items()}


def summarize(lines: dict, bin_: float = BIN):
    """Per molecule: ``(bin_start, n_lines, max S/(pi gamma))`` for every bin that has lines."""
    out = {}
    for mol, (wn, s, gamma) in lines.items():
        gamma = np.where(gamma > 0, gamma, 0.07)               # a missing width: a typical surface value
        peak = s / (np.pi * gamma)
        idx = np.floor(wn / bin_).astype(int)
        order = np.argsort(idx, kind="stable")
        idx, peak = idx[order], peak[order]
        starts = np.concatenate([[0], np.flatnonzero(np.diff(idx)) + 1])
        counts = np.diff(np.concatenate([starts, [len(idx)]]))
        best = np.maximum.reduceat(peak, starts)
        out[mol] = [(float(idx[i] * bin_), int(c), float(b)) for i, c, b in zip(starts, counts, best)]
    return out


def write(summary: dict, path: str) -> None:
    with open(path, "w", newline="") as fh:
        fh.write("# Per molecule and %g cm-1 bin of the AER line list: line count and the largest S/(pi*gamma)\n" % BIN)
        fh.write("# (S: cm-1/(molecule cm-2) at 296 K, gamma: air half width, cm-1). Made by absco-active-regions.\n")
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["molecule", "bin_start", "n_lines", "max_s_over_pi_gamma"])
        for mol in sorted(summary):
            for b, n, v in summary[mol]:
                w.writerow([mol, "%g" % b, n, "%.4e" % v])


def read(path: str) -> dict:
    out: dict[str, list] = {}
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh) if r and not r[0].startswith("#")]
    for mol, b, n, v in rows[1:]:
        out.setdefault(mol, []).append((float(b), int(n), float(v)))
    return out


def columns(profile_csv: str) -> dict:
    """Molecules cm-2 of each molecule in the packaged standard profile (level amounts of the layers)."""
    with open(profile_csv, newline="") as fh:
        rows = list(csv.DictReader(fh))
    p = np.array([float(r["P"]) for r in rows]) * 100.0                      # Pa
    n_air = np.abs(np.diff(p)) / (G * M_AIR) / 1e4                           # molecules cm-2 per layer
    out = {}
    for mol in set(MOLECULE_NUMBER.values()):
        key = PROFILE_COLUMN.get(mol, mol)
        if key not in rows[0]:
            continue
        v = np.array([float(r[key]) for r in rows])
        out[mol] = float(np.sum(0.5 * (v[1:] + v[:-1]) * n_air))
    return out


def peak_tau(summary: dict, column: dict) -> dict:
    """``{molecule: [(bin_start, tau_peak), ...]}`` for the molecules with a column."""
    return {m: [(b, v * column[m]) for b, _, v in bins] for m, bins in summary.items() if m in column}
