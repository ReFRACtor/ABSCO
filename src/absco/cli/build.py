"""``absco-build`` -- developer command to compile LNFL/LBLRTM from source.

For contributors working from a source checkout (rather than a prebuilt binary
wheel).  It compiles the Fortran executables from the ``LNFL`` and ``LBLRTM``
submodule directories and stages them into the ABSCO data directory, giving the
same layout an ``absco-init`` user would have.  Optionally also fetches the AER line
file so a from-source developer is fully set up in one command.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from absco import artifacts, paths


def _source_root():
    """Best guess at the repository checkout containing the LNFL/LBLRTM submodules.

    For an editable install this package lives at ``<repo>/src/absco/cli/build.py``,
    so the checkout is three parents up. Used only to default ``--lnfl-path`` /
    ``--lblrtm-path`` so ``absco-build`` works from any working directory; an
    explicit flag always overrides.
    """
    return Path(__file__).resolve().parents[3]


def build_parser():
    src_root = _source_root()
    parser = argparse.ArgumentParser(
        prog="absco-build",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Compile the LNFL and LBLRTM Fortran executables from a source "
        "checkout and stage them into the ABSCO data directory.",
    )
    parser.add_argument(
        "-c",
        "--compiler",
        default="gfortran",
        help="Fortran compiler to build with (gfortran, ifort, pgf90). "
        "Case-insensitive.",
    )
    parser.add_argument(
        "--lnfl-path",
        default=os.fspath(src_root / "LNFL"),
        help="Top-level LNFL source directory (contains build/).",
    )
    parser.add_argument(
        "--lblrtm-path",
        default=os.fspath(src_root / "LBLRTM"),
        help="Top-level LBLRTM source directory (contains build/).",
    )
    parser.add_argument(
        "--data-dir",
        help="Directory for staged artifacts. Overrides $ABSCO_DATA_DIR and the "
        "platform default for this run.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite already-staged executables.",
    )
    parser.add_argument(
        "--lines",
        action="store_true",
        help="Also download and stage the AER line file from Zenodo.",
    )
    parser.add_argument(
        "--record",
        type=int,
        default=artifacts.DEFAULT_ZENODO_RECORD,
        help="Zenodo record ID for the AER line file (used with --lines).",
    )
    return parser


def main():
    args = build_parser().parse_args()

    if args.data_dir:
        os.environ[paths.DATA_DIR_ENV] = os.path.abspath(
            os.path.expanduser(args.data_dir)
        )

    data_dir = paths.data_dir(create=True)
    print(f"Using ABSCO data directory: {data_dir}")

    lnfl_exe = artifacts.compile_model("lnfl", args.lnfl_path, compiler=args.compiler)
    artifacts.stage_executable(lnfl_exe, force=args.force)

    lbl_exe = artifacts.compile_model("lblrtm", args.lblrtm_path, compiler=args.compiler)
    artifacts.stage_executable(lbl_exe, force=args.force)
    # LBLRTM v12.17 reads MT_CKD continuum data (netCDF) from its run dir
    artifacts.stage_lblrtm_data_files(args.lblrtm_path, force=args.force)

    if args.lines:
        extracted = artifacts.fetch_line_file(record=args.record, force=args.force)
        artifacts.stage_line_file(extracted, force=args.force)

    print("\nBuild complete.")
    print(f"  lnfl exe  : {paths.lnfl_exe()}")
    print(f"  lblrtm exe: {paths.lblrtm_exe()}")
    if args.lines:
        print(f"  line file : {paths.line_file_root()}")


if __name__ == "__main__":
    main()
