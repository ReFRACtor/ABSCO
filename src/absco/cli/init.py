"""``absco-init`` -- initialize the runtime artifacts for an installed ABSCO tool.

Intended for users who installed a (prebuilt-binary) wheel: it stages the Fortran
executables into the user data directory and downloads/extracts the AER line file
from Zenodo.  Idempotent -- artifacts already present are skipped unless ``--force``
is given.  Path resolution at run time is handled by :mod:`absco.paths`, so this
command never edits ``ABSCO_config.ini``.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from absco import artifacts, paths


def build_parser():
    parser = argparse.ArgumentParser(
        prog="absco-init",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Download the AER line file and stage executables into the "
        "ABSCO data directory so the tool can run from any working directory.",
    )
    parser.add_argument(
        "--data-dir",
        help="Directory for runtime artifacts. Overrides $ABSCO_DATA_DIR and the "
        "platform default for this run.",
    )
    parser.add_argument(
        "--record",
        type=int,
        default=artifacts.DEFAULT_ZENODO_RECORD,
        help="Zenodo record ID for the AER line file.",
    )
    parser.add_argument(
        "--lines-only",
        action="store_true",
        help="Only fetch/stage the line file; do not stage executables.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download/re-stage artifacts even if they already exist.",
    )
    return parser


def _stage_bundled_executables(force: bool) -> None:
    """Stage wheel-bundled executables (absco/_bin) into the data dir bin.

    If the wheel did not bundle binaries, report that the user should build them
    with ``absco-build`` instead; this is not fatal so the line-file download can
    still proceed.
    """
    wheel_bin = Path(__file__).resolve().parent.parent / "_bin"
    exes = sorted(p for p in wheel_bin.glob("*") if p.is_file()) if wheel_bin.is_dir() else []

    if not exes:
        print(
            "No wheel-bundled executables found in {}.\n"
            "  If you installed a source/pure-python build, compile them with "
            "`absco-build` instead.".format(wheel_bin)
        )
        return

    for exe in exes:
        artifacts.stage_executable(exe, force=force)


def main():
    args = build_parser().parse_args()

    if args.data_dir:
        os.environ[paths.DATA_DIR_ENV] = os.path.abspath(
            os.path.expanduser(args.data_dir)
        )

    data_dir = paths.data_dir(create=True)
    print(f"Using ABSCO data directory: {data_dir}")

    if not args.lines_only:
        _stage_bundled_executables(force=args.force)

    extracted = artifacts.fetch_line_file(record=args.record, force=args.force)
    artifacts.stage_line_file(extracted, force=args.force)

    # report final resolution status
    lnfl, lbl = paths.lnfl_exe(), paths.lblrtm_exe()
    lf = paths.line_file_paths()
    xs_ok = os.path.isdir(lf["xs_path"]) and os.path.isfile(lf["fscdxs"])
    print("\nInitialization complete.")
    print(f"  line file : {paths.line_file_root()}")
    print(f"  lnfl exe  : {lnfl or '(not found -- run absco-build)'}")
    print(f"  lblrtm exe: {lbl or '(not found -- run absco-build)'}")
    if not xs_ok:
        print("  cross sections: (not found -- AER line file v3.9 no longer ships "
              "them; run `absco-build` to stage the LBLRTM cross-sections)")


if __name__ == "__main__":
    main()
