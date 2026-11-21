#!/usr/bin/env python
"""absco-active-regions: compare the active-region file with the AER line list.

Reads the line list (the ``tape1_path`` of a configuration, by default the one in the data
directory), reduces it to the largest S/(pi gamma) per molecule and 10 cm-1 bin, and with the
column of the standard profile gives a peak optical depth per bin.  Two outputs:

* ``--write-summary FILE`` writes the reduced table (the packaged ``active_regions_summary.csv``);
* the default is a report: for each molecule, the bins above ``--threshold`` that no region
  covers ("missing"), and the ``lines`` regions with no bin above it ("no signal").  ``--check``
  makes the exit status non-zero if bins are missing that are not in the file's accepted gaps.
"""

import argparse
import os
import sys

import numpy as np

from absco import active_regions, line_summary, paths

#: tau threshold of the coverage check and the wavenumber range it applies to
THRESHOLD = 1e-3
CHECK_RANGE = (400.0, 50000.0)


def covered(regions, bin_start, bin_):
    mid = bin_start + bin_ / 2.0
    return any(r.wn1 <= mid <= r.wn2 for r in regions)


def missing_bins(summary, column, regions, threshold=THRESHOLD, lo=CHECK_RANGE[0], hi=CHECK_RANGE[1], accepted=()):
    """``{molecule: [bin_start, ...]}`` of bins with peak tau above ``threshold`` that no region
    of the molecule covers, inside [lo, hi], apart from the bins listed in ``accepted``."""
    by_mol = {}
    for r in regions:
        by_mol.setdefault(r.molecule, []).append(r)
    skip = {(m, float(b)) for m, b in accepted}
    out = {}
    for mol, bins in line_summary.peak_tau(summary, column).items():
        for b, tau in bins:
            if tau > threshold and lo <= b < hi and (mol, b) not in skip and not covered(by_mol.get(mol, []), b, line_summary.BIN):
                out.setdefault(mol, []).append(b)
    return out


def quiet_regions(summary, column, regions, threshold=THRESHOLD):
    """``lines`` regions in which no bin of the molecule reaches ``threshold`` (a stale entry, or a weak band)."""
    tau = {m: dict(v) for m, v in line_summary.peak_tau(summary, column).items()}
    out = []
    for r in regions:
        if "lines" not in r.basis or r.molecule not in tau:
            continue
        bins = [b for b in tau[r.molecule] if r.wn1 <= b + line_summary.BIN / 2 <= r.wn2]
        if not any(tau[r.molecule][b] > threshold for b in bins):
            out.append(r)
    return out


def ranges(starts, bin_=line_summary.BIN):
    out = []
    for b in sorted(starts):
        if out and abs(b - out[-1][1]) < 1e-9:
            out[-1][1] = b + bin_
        else:
            out.append([b, b + bin_])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Compare the active-region file with the AER line list.")
    ap.add_argument("--line-file", help="AER line file (default: the one in the data directory)")
    ap.add_argument("--summary", help="reduced line table to read instead of the line file (default: the packaged one)")
    ap.add_argument("--write-summary", metavar="FILE", help="write the reduced line table made from the line file and stop")
    ap.add_argument("--regions", help="region file to check (default: the packaged one)")
    ap.add_argument("--profile", default=paths.data_file("VMR/USS_AIRS_profile.csv"))
    ap.add_argument("--threshold", type=float, default=THRESHOLD, help="peak optical depth (default %(default)g)")
    ap.add_argument("--check", action="store_true", help="exit with status 1 if bins above the threshold are not covered")
    a = ap.parse_args(argv)

    if a.write_summary or a.line_file:
        path = a.line_file or paths.line_file_paths()["tape1_path"]
        print("reading %s" % path, flush=True)
        summary = line_summary.summarize(line_summary.read_lines(path))
        if a.write_summary:
            line_summary.write(summary, a.write_summary)
            print("wrote %s (%d molecules)" % (a.write_summary, len(summary)))
            return 0
    else:
        summary = line_summary.read(a.summary or paths.data_file(line_summary.FILE))

    column = line_summary.columns(a.profile)
    regions = active_regions.load(a.regions)
    miss = missing_bins(summary, column, regions, a.threshold)
    quiet = quiet_regions(summary, column, regions, a.threshold)
    print("threshold: peak optical depth %g of the standard column, %g-%g cm-1" % ((a.threshold,) + CHECK_RANGE))
    if not miss:
        print("every bin above the threshold is covered by a region")
    for mol in sorted(miss):
        print("  %-6s not covered: %s" % (mol, ", ".join("%d-%d" % tuple(r) for r in ranges(miss[mol]))))
    for r in quiet:
        print("  %-6s region %g-%g (%s, source %s): no bin above the threshold" % (r.molecule, r.wn1, r.wn2, "+".join(r.basis), r.source))
    return 1 if (a.check and miss) else 0


if __name__ == "__main__":
    sys.exit(main())
