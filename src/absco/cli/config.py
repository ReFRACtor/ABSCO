"""``absco-config`` -- assistant that writes an ABSCO_config.ini for absco-generate.

The user typically supplies only the spectral range, the output resolution, and the
molecules; the tool suggests a matching LBLRTM resolution (``lblres``) and writes a
config from the bundled template, leaving data/executable/line-file paths blank so
they resolve from the installed package / data dir at run time.

Inputs come from flags, with an interactive prompt fallback for any required value
that was omitted.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

from absco import config_assistant as ca


def build_parser():
    parser = argparse.ArgumentParser(
        prog="absco-config",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Generate an ABSCO_config.ini. Specify the output resolution, "
        "spectral range, and molecules; lblres is suggested automatically.",
    )
    parser.add_argument(
        "--outres",
        type=float,
        nargs="+",
        help="Output (degraded) resolution in --units. One value applied to all "
        "bands, or one per band. For um/nm it is a wavelength spacing and is "
        "converted to a cm-1 spacing at each band center. Required (prompted "
        "if omitted).",
    )
    parser.add_argument(
        "--begin",
        type=float,
        nargs="+",
        help="Band start value(s) in --units (wavenumber or wavelength). Pair with --end.",
    )
    parser.add_argument(
        "--end",
        type=float,
        nargs="+",
        help="Band end value(s) in --units (wavenumber or wavelength). Pair with --begin.",
    )
    parser.add_argument(
        "--range",
        type=float,
        nargs=2,
        metavar=("BEGIN", "END"),
        help="Convenience for a single band: --range BEGIN END (in --units).",
    )
    parser.add_argument(
        "--molnames",
        nargs="+",
        help="HITRAN molecule names (space separated, case-insensitive). "
        "Required (prompted if omitted).",
    )
    parser.add_argument(
        "--units",
        default="cm-1",
        choices=["cm-1", "um", "nm"],
        help="Units of wn1/wn2/outres. um/nm inputs are converted to cm-1 and the "
        "config is always written in cm-1 (LBLRTM's native wavenumber grid).",
    )
    parser.add_argument(
        "--wv-vmr",
        type=float,
        nargs="+",
        help="Water vapor VMR values [ppmv] for the [vmr] wv_vmr field "
        "(default: template value).",
    )
    parser.add_argument(
        "--lblres",
        type=float,
        nargs="+",
        help="Override the suggested LBLRTM resolution(s) [cm-1]. Must keep "
        "outres/lblres an exact power of 2.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="ABSCO_config.ini",
        help="Path to write the config file.",
    )
    parser.add_argument(
        "-y",
        "--no-prompt",
        action="store_false",
        dest="interactive",
        default=True,
        help="Never prompt; error out if a required value is missing.",
    )
    return parser


def _prompt(message):
    """Prompt the user for a line of input (returns stripped string)."""
    return input(message).strip()


def _resolve_bands(args):
    """Determine begin/end arrays from --range or --begin/--end, prompting if needed."""
    begin, end = args.begin, args.end

    if args.range is not None:
        if begin or end:
            sys.exit("Use either --range or --begin/--end, not both.")
        begin, end = [args.range[0]], [args.range[1]]

    if not begin or not end:
        if not args.interactive:
            sys.exit("Spectral range required: use --range or --begin/--end.")
        raw = _prompt("Spectral range as 'BEGIN END' (in %s): " % args.units)
        parts = raw.split()
        if len(parts) != 2:
            sys.exit("Please provide exactly two numbers for the range.")
        begin, end = [float(parts[0])], [float(parts[1])]

    begin = np.atleast_1d(np.asarray(begin, dtype=float))
    end = np.atleast_1d(np.asarray(end, dtype=float))
    if begin.size != end.size:
        sys.exit("--begin and --end must have the same number of values.")
    return begin, end


def _resolve_outres(args, n_bands):
    """Determine the outres array (one per band), prompting if needed."""
    outres = args.outres
    if not outres:
        if not args.interactive:
            sys.exit("Output resolution required: use --outres.")
        raw = _prompt("Output resolution outres [cm-1] (one value, or one per band): ")
        parts = raw.split()
        if not parts:
            sys.exit("Please provide at least one outres value.")
        outres = [float(p) for p in parts]

    outres = np.atleast_1d(np.asarray(outres, dtype=float))
    if outres.size == 1 and n_bands > 1:
        outres = np.repeat(outres, n_bands)
    if outres.size != n_bands:
        sys.exit(
            "--outres must have 1 value or one per band (%d bands)." % n_bands
        )
    return outres


def _resolve_molecules(args):
    """Determine the molecule list, prompting if needed."""
    mols = args.molnames
    if not mols:
        if not args.interactive:
            sys.exit("Molecules required: use --molnames.")
        raw = _prompt("Molecules (space separated, e.g. 'ch4 co h2o'): ")
        mols = raw.split()
    if not mols:
        sys.exit("At least one molecule is required.")
    return mols


def main():
    args = build_parser().parse_args()

    wn1, wn2 = _resolve_bands(args)
    outres = _resolve_outres(args, wn1.size)
    molnames = _resolve_molecules(args)

    # LBLRTM works on a wavenumber grid: convert wavelength (um/nm) inputs to cm-1
    # so lblres, the RAM estimate, and the written config are all in cm-1.
    try:
        wn1_cm1, wn2_cm1, outres_cm1 = ca.convert_bands_to_cm1(
            wn1, wn2, outres, args.units)
    except ValueError as exc:
        sys.exit("Error: %s" % exc)

    # an explicit --lblres is always a cm-1 spacing; a suggestion derives from the
    # (converted) cm-1 output resolution
    lblres = ca.suggest_lblres(outres_cm1) if args.lblres is None else \
        np.atleast_1d(np.asarray(args.lblres, dtype=float))

    header_note = None
    if args.units.lower() != "cm-1":
        header_note = [
            "Generated by absco-config from %s inputs, converted to cm-1." % args.units,
            "wn1/wn2/lblres/outres below are wavenumbers [cm-1] (LBLRTM's native grid).",
        ]

    try:
        config = ca.build_config(
            wn1, wn2, outres, molnames,
            units=args.units, wv_vmr=args.wv_vmr, lblres=lblres,
        )
    except ValueError as exc:
        sys.exit("Error: %s" % exc)

    with open(args.output, "w") as fh:
        ca.write_config(config, fh, header_note=header_note)

    # report (all in cm-1, as written)
    print("Wrote %s" % args.output)
    if args.units.lower() != "cm-1":
        print("  (converted from %s to cm-1)" % args.units)
    print("  units   = cm-1")
    print("  wn1     = %s" % config["channels"]["wn1"])
    print("  wn2     = %s" % config["channels"]["wn2"])
    print("  outres  = %s" % config["channels"]["outres"])
    print("  lblres  = %s  (suggested; outres/lblres is a power of 2)"
          % config["channels"]["lblres"])
    print("  molecules = %s" % config["molecules"]["molnames"])

    ram = ca.estimate_ram_gb(wn1_cm1, wn2_cm1, outres_cm1, molnames)
    if ram is not None:
        print("\nEstimated peak RAM for generation: up to %.3f GB" % ram)


if __name__ == "__main__":
    main()
