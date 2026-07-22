"""Resolution of bundled data files and runtime artifacts for the ABSCO tool.

Two categories of files are located here:

**Bundled data** shipped inside the wheel (:mod:`absco.data`): pressure/temperature
grids, VMR profiles, the XS/line lookup CSV, and the config template.  These are
resolved with :func:`data_file` via :mod:`importlib.resources`, so they work from an
installed package without the user specifying any path.

**Runtime artifacts** that are too large or platform-specific to ship on PyPI: the
Fortran executables (LNFL, LBLRTM) and the AER line file database.  These live in a
user data directory resolved by :func:`data_dir` (``$ABSCO_DATA_DIR`` if set,
otherwise a :mod:`platformdirs` location).  ``absco-init`` populates that directory;
:func:`lnfl_exe`, :func:`lblrtm_exe`, and :func:`line_file_paths` locate the pieces
within it.

Every resolver returns an absolute path.  A user overrides any default simply by
setting the corresponding value in ``ABSCO_config.ini`` to an explicit path.
"""

from __future__ import annotations

import glob
import os
from importlib.resources import as_file, files
from pathlib import Path

import platformdirs

__all__ = [
    "data_file",
    "config_template",
    "data_dir",
    "bin_dir",
    "line_file_root",
    "lnfl_exe",
    "lblrtm_exe",
    "line_file_paths",
    "lblrtm_data_files",
]

# Data files LBLRTM (>= v12.11) reads from its run directory at run time, e.g. the
# MT_CKD water-vapor continuum coefficients. Staged into <data_dir>/bin alongside
# the executables and symlinked into the LBL run dir by absco.compute.
LBLRTM_RUNTIME_DATA = ["absco-ref_wv-mt-ckd.nc"]

# Environment variable that overrides the runtime data directory.
DATA_DIR_ENV = "ABSCO_DATA_DIR"

# Default relative names for the bundled data files (under absco/data).
DEFAULT_DATA_FILES = {
    "pfile": "PT_grid/AIRS_P_air.txt",
    "ptfile": "PT_grid/build_temp_array.txt",
    "vmrfile": "VMR/USS_AIRS_profile.csv",
    "hdofile": "VMR/HDO_example_profile.csv",
    "xs_lines": "FSCDXS_line_params.csv",
}

# Layout of the AER line file database beneath the data dir.  Mirrors the structure
# staged by absco.artifacts (line_file/<aer_ver>, line_file/lncpl_lines, ...).
LINE_FILE_SUBDIR = "AER_Line_File"


def data_file(relpath: str) -> str:
    """Return an absolute path to a bundled data file under ``absco.data``.

    ``relpath`` is a POSIX-style path relative to the package data root, e.g.
    ``"PT_grid/AIRS_P_air.txt"``.  Works for a normal filesystem install; for the
    (rare) zipped-install case the resource is materialized to a temp file.
    """
    resource = files("absco.data").joinpath(relpath)
    try:
        # Fast path: the resource is a real file on disk (the normal wheel install).
        return os.fspath(Path(str(resource)).resolve())
    except (TypeError, ValueError):
        # Fallback for non-filesystem loaders (zipimport): extract to a temp path.
        with as_file(resource) as p:
            return os.fspath(Path(p).resolve())


def default_data_file(key: str) -> str:
    """Absolute path to the packaged default for a config data-file key."""
    return data_file(DEFAULT_DATA_FILES[key])


def config_template() -> str:
    """Absolute path to the bundled ``ABSCO_config.template.ini``."""
    return data_file("ABSCO_config.template.ini")


def data_dir(create: bool = False) -> Path:
    """Return the runtime data directory for artifacts (line file, executables).

    Resolution order: ``$ABSCO_DATA_DIR`` if set, else
    ``platformdirs.user_data_dir("absco")``.  With ``create=True`` the directory is
    created if missing.
    """
    env = os.environ.get(DATA_DIR_ENV)
    root = Path(env).expanduser() if env else Path(platformdirs.user_data_dir("absco"))
    root = root.resolve()
    if create:
        root.mkdir(parents=True, exist_ok=True)
    return root


def bin_dir() -> Path:
    """Directory that holds the Fortran executables in the data dir (``<data_dir>/bin``)."""
    return data_dir() / "bin"


def line_file_root() -> Path:
    """Top-level directory of the AER line file database in the data dir."""
    return data_dir() / LINE_FILE_SUBDIR


def _wheel_bin_dir() -> Path:
    """Directory of executables bundled inside the wheel (``absco/_bin``), if any."""
    return Path(__file__).resolve().parent / "_bin"


def _find_exe(stem: str) -> str | None:
    """Locate an executable named/prefixed ``stem``.

    Searches, in order: the wheel-bundled ``absco/_bin`` directory, then
    ``<data_dir>/bin``.  In each directory an exact-name match wins; otherwise a
    versioned build (e.g. ``lnfl_v3.2_linux_gnu_sgl``) is matched by glob.
    """
    for d in (_wheel_bin_dir(), bin_dir()):
        exact = d / stem
        if exact.is_file():
            return os.fspath(exact.resolve())
        matches = sorted(glob.glob(os.fspath(d / f"{stem}_*")))
        if matches:
            return os.fspath(Path(matches[0]).resolve())
    return None


def lnfl_exe() -> str | None:
    """Absolute path to the LNFL executable, or ``None`` if not initialized."""
    return _find_exe("lnfl")


def lblrtm_exe() -> str | None:
    """Absolute path to the LBLRTM executable, or ``None`` if not initialized."""
    return _find_exe("lblrtm")


def line_file_paths() -> dict:
    """Return default paths for the AER line-file components in the data dir.

    Keys mirror the ``ABSCO_config.ini`` fields: ``tape1_path``, ``tape2_path``,
    ``extra_params``, ``xs_path``, ``fscdxs``.  ``tape1_path`` points at the AER
    line-parameter directory (``line_file/aer_v_*``), discovered by glob so the exact
    version-stamped name does not need to be hard-coded; if none is present the
    conventional ``line_file/aer_v_3.6`` path is returned so error messages are useful.
    """
    root = line_file_root()

    aer_matches = sorted(glob.glob(os.fspath(root / "line_file" / "aer_v_*")))
    tape1 = Path(aer_matches[0]) if aer_matches else root / "line_file" / "aer_v_3.6"

    return {
        "tape1_path": os.fspath(tape1),
        "tape2_path": os.fspath(root / "line_file" / "lncpl_lines"),
        "extra_params": os.fspath(root / "extra_brd_params"),
        "xs_path": os.fspath(root / "xs_files" / "xs"),
        "fscdxs": os.fspath(root / "xs_files" / "FSCDXS"),
    }


def lblrtm_data_files() -> list:
    """Absolute paths to LBLRTM runtime data files staged in the data dir.

    Looks in the wheel-bundled ``absco/_bin`` first, then ``<data_dir>/bin`` (where
    ``absco-build`` / ``absco-init`` stage them). Returns only the files that exist,
    so callers can symlink whatever is present into the LBL run directory.
    """
    found = []
    for name in LBLRTM_RUNTIME_DATA:
        for d in (_wheel_bin_dir(), bin_dir()):
            candidate = d / name
            if candidate.is_file():
                found.append(os.fspath(candidate.resolve()))
                break
    return found
