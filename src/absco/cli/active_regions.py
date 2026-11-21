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


def print_regions(regions):
    """One block per molecule: its regions in cm-1 with basis and source, and the total width."""
    by_mol = {}
    for r in regions:
        by_mol.setdefault(r.molecule, []).append(r)
    print("ACTIVE REGIONS PER MOLECULE")
    print("%-6s %10s %10s  %-22s %s" % ("mol", "wn1 [cm-1]", "wn2 [cm-1]", "basis", "source"))
    for mol in sorted(by_mol):
        rs = sorted(by_mol[mol], key=lambda r: r.wn1)
        for i, r in enumerate(rs):
            print("%-6s %10g %10g  %-22s %s" % (mol if i == 0 else "", r.wn1, r.wn2, "+".join(r.basis), r.source))
        if len(rs) > 1:
            print("%-6s %10s %10g  (total width, %d regions)" % ("", "", sum(r.wn2 - r.wn1 for r in rs), len(rs)))


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Report the spectral regions in which each molecule is active, and compare them with the AER line list.",
        epilog="""\
With no arguments, using the packaged region file, line-list summary and profile, it prints:
  1. ACTIVE REGIONS PER MOLECULE: a table with one row per region (wn1 and wn2 in cm-1, basis,
     source), and the total width for a molecule that has several regions.
  2. DEVIATIONS OF THE REGION FILE FROM THE LINE LIST:
       - bins above the threshold that no region covers (these make --check fail);
       - regions with lines but no bin above the threshold (reported only).
Exit status: 0, or 1 if --check is given and some bin above the threshold is not covered.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line-file", metavar="FILE",
                    help="AER line list (HITRAN-format text file, about 800 MB) to read instead of the packaged "
                         "summary of it; reading takes about 25 s (default: use the packaged summary)")
    ap.add_argument("--summary", metavar="CSV",
                    help="CSV of the largest S/(pi*gamma) per molecule and 10 cm-1 bin, as written by "
                         "--write-summary (default: the one packaged with absco)")
    ap.add_argument("--write-summary", metavar="CSV",
                    help="read the line list (--line-file, or the one named in the configuration), write this "
                         "summary CSV and stop")
    ap.add_argument("--regions", metavar="CSV",
                    help="active-region CSV to report and check, with columns molecule,wn1,wn2,basis,source,"
                         "reference,note (default: the one packaged with absco, src/absco/data/active_regions.csv)")
    ap.add_argument("--profile", metavar="CSV", default=paths.data_file("VMR/USS_AIRS_profile.csv"),
                    help="standard atmosphere profile CSV (columns ALT,P,T and one volume mixing ratio column "
                         "per molecule) that sets each molecule's column (default: %(default)s)")
    ap.add_argument("--threshold", type=float, default=THRESHOLD, help="peak optical depth of one standard column above which a 10 cm-1 bin counts as active (default %(default)g)")
    ap.add_argument("--check", action="store_true", help="exit with status 1 if bins above the threshold are not covered")
    only = ap.add_mutually_exclusive_group()
    only.add_argument("--summary-only", action="store_true", help="print only the regions per molecule")
    only.add_argument("--deviations-only", action="store_true", help="print only the deviations from the line list")
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
    if not a.deviations_only:
        print_regions(regions)
    if not a.summary_only:
        if not a.deviations_only:
            print()
        print("DEVIATIONS OF THE REGION FILE FROM THE LINE LIST")
        print("peak optical depth threshold %g of the standard column, checked over %g-%g cm-1" % ((a.threshold,) + CHECK_RANGE))
        print("lines above the threshold that no region covers (these fail --check):")
        if not miss:
            print("  none")
        for mol in sorted(miss):
            print("  %-6s %s" % (mol, ", ".join("%d-%d" % tuple(r) for r in ranges(miss[mol]))))
        print("regions with lines but no bin above the threshold (reported, do not fail --check):")
        if not quiet:
            print("  none")
        for r in quiet:
            print("  %-6s %g-%g (%s, source %s)" % (r.molecule, r.wn1, r.wn2, "+".join(r.basis), r.source))
    return 1 if (a.check and miss) else 0


if __name__ == "__main__":
    sys.exit(main())
