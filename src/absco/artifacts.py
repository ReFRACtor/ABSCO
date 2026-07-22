"""Build, fetch, and stage the runtime artifacts ABSCO needs.

This module holds the three separable concerns that used to live in
``build_models.py``:

* :func:`compile_model` -- compile the LNFL or LBLRTM Fortran executable from a
  source checkout (``make -f make_lnfl`` / ``make_lblrtm``).
* :func:`fetch_line_file` -- download and extract the AER line file database from
  Zenodo.
* :func:`stage_executable` / :func:`stage_line_file` -- place compiled/bundled
  executables and the extracted line file into the user data directory layout that
  :mod:`absco.paths` resolves (``<data_dir>/bin`` and ``<data_dir>/AER_Line_File``).

The commands ``absco-build`` (developer compile) and ``absco-init`` (end-user
artifact setup) are thin wrappers around these functions.  Unlike the old
``build_models.py`` this module never rewrites ``ABSCO_config.ini``; path resolution
happens at runtime in :mod:`absco.paths`.
"""

from __future__ import annotations

import contextlib
import glob
import os
import shutil
import stat
import subprocess
import sys
import sysconfig
import tarfile
import tempfile
from pathlib import Path

from absco import paths

__all__ = [
    "DEFAULT_ZENODO_RECORD",
    "COMPILER_MAP",
    "platform_os",
    "compile_model",
    "fetch_line_file",
    "stage_executable",
    "stage_line_file",
]

# Default Zenodo record for the AER line file (from the original build_models.py).
DEFAULT_ZENODO_RECORD = 3837550

# ABSCO compiler name -> (make FC_TYPE label, ini/exe token).
COMPILER_MAP = {
    "gfortran": "GNU",
    "ifort": "INTEL",
    "pgf90": "PGI",
}

# Per-model make configuration: (product name, makefile, precision token).
_MODEL_MAKE = {
    "lnfl": ("lnfl", "make_lnfl", "sgl"),
    "lblrtm": ("lblrtm", "make_lblrtm", "dbl"),
}


def platform_os() -> str:
    """Return the LBLRTM/LNFL platform token for the current OS."""
    if sys.platform in ("linux", "linux2"):
        return "linux"
    if sys.platform == "darwin":
        return "osx"
    if sys.platform == "win32":
        return "mingw"
    raise RuntimeError("Could not determine OS platform for the Fortran build")


@contextlib.contextmanager
def _gfortran_compat_path(compiler: str):
    """Yield an env whose PATH shims ``gfortran`` to allow argument mismatches.

    Only active for gfortran builds and only when the real gfortran accepts
    ``-fallow-argument-mismatch`` (gfortran >= 10).  Otherwise the current
    environment is yielded unchanged.
    """
    env = dict(os.environ)
    if compiler.lower() != "gfortran":
        yield env
        return

    real = shutil.which("gfortran")
    # -std=legacy relaxes obsolete-feature errors; -fallow-argument-mismatch
    # degrades external-call argument mismatches (type/rank/element-count) to
    # warnings. Both are needed to build the legacy AER code on gfortran >= 10.
    candidate_flags = ["-std=legacy", "-fallow-argument-mismatch"]
    flags = [f for f in candidate_flags if real and _accepts_flag(real, f)]
    if real is None or not flags:
        yield env
        return

    tmpdir = tempfile.mkdtemp(prefix="absco-fc-")
    try:
        shim = Path(tmpdir) / "gfortran"
        shim.write_text(f'#!/bin/sh\nexec "{real}" {" ".join(flags)} "$@"\n')
        shim.chmod(shim.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
        env["PATH"] = tmpdir + os.pathsep + env.get("PATH", "")
        yield env
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _accepts_flag(compiler_exe: str, flag: str) -> bool:
    """Return True if ``compiler_exe`` accepts ``flag`` (probed via --help)."""
    try:
        out = subprocess.run(
            [compiler_exe, flag, "--version"],
            capture_output=True,
            text=True,
        )
        return out.returncode == 0
    except OSError:
        return False


def compile_model(model: str, source_dir, compiler: str = "gfortran") -> str:
    """Compile LNFL or LBLRTM from a source checkout and return the exe path.

    ``model`` is ``"lnfl"`` or ``"lblrtm"``; ``source_dir`` is the top-level
    submodule directory (containing ``build/``).  Mirrors the make invocation of the
    original ``build_models.py``: ``make -f <makefile> <os><COMPILER><precision>`` run
    from ``<source_dir>/build``.  Returns the absolute path to the produced executable.
    """
    model = model.lower()
    if model not in _MODEL_MAKE:
        raise ValueError(f"Unknown model {model!r}; expected 'lnfl' or 'lblrtm'")

    compiler = compiler.lower()
    if compiler not in COMPILER_MAP:
        raise ValueError(
            f"{compiler!r} is not a supported compiler "
            f"({', '.join(sorted(COMPILER_MAP))})"
        )

    product, makefile, precision = _MODEL_MAKE[model]
    source_dir = Path(source_dir).resolve()
    build_dir = source_dir / "build"
    if not build_dir.is_dir():
        raise FileNotFoundError(f"No build directory found at {build_dir}")

    target = f"{platform_os()}{COMPILER_MAP[compiler]}{precision}"
    print(f"Building {product} (target {target}) in {build_dir}")

    # The LNFL/LBLRTM makefiles pin FC=gfortran and the FCFLAG list, so there is no
    # clean make-level hook to add flags. gfortran >= 10 rejects the legacy code's
    # argument type mismatches by default; shim `gfortran` on PATH to append
    # -fallow-argument-mismatch (a no-op on older gfortran that lack the flag).
    with _gfortran_compat_path(compiler) as env:
        status = subprocess.call(
            ["make", "-f", makefile, target], cwd=os.fspath(build_dir), env=env
        )
    if status != 0:
        raise RuntimeError(f"{product} build failed (make returned {status})")

    # make writes the exe one level up as <product>_<version>_<os>_<comp>_<prec>
    pattern = os.fspath(
        source_dir / f"{product}_*_{platform_os()}_{compiler}_{precision}"
    )
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"Build reported success but no executable matched {pattern}"
        )
    return os.fspath(Path(matches[0]).resolve())


def stage_executable(exe_path, dest_dir=None, force: bool = False) -> str:
    """Copy a compiled/bundled executable into the data dir ``bin`` directory.

    Preserves the versioned filename so :func:`absco.paths._find_exe` can match it.
    Returns the staged path.  Skips the copy when an identically named file already
    exists unless ``force`` is set.
    """
    exe_path = Path(exe_path).resolve()
    if dest_dir is None:
        dest_dir = paths.bin_dir()
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    dest = dest_dir / exe_path.name
    if dest.exists() and not force:
        print(f"{dest} already exists, skipping")
        return os.fspath(dest)

    shutil.copy2(exe_path, dest)
    dest.chmod(dest.stat().st_mode | 0o111)
    print(f"Staged {exe_path.name} -> {dest}")
    return os.fspath(dest)


def _zenodo_cli():
    """Return the argv prefix for invoking the zenodo_get console script.

    Prefer the installed ``zenodo_get`` entry point; fall back to
    ``python -m zenodo_get`` so the build works even if the script dir is not on
    PATH.  Using the CLI (rather than importing zenodo_get internals) keeps this
    stable across zenodo_get versions.
    """
    exe = shutil.which("zenodo_get")
    if exe:
        return [exe]
    scripts = sysconfig.get_path("scripts")
    candidate = Path(scripts) / "zenodo_get"
    if candidate.exists():
        return [os.fspath(candidate)]
    return [sys.executable, "-m", "zenodo_get"]


def fetch_line_file(
    record=DEFAULT_ZENODO_RECORD, download_dir=None, force: bool = False
):
    """Download and extract the AER line file tarball from Zenodo.

    Returns the path to the extracted top-level directory (the ``.tar.gz`` name
    minus the extension), which follows the AER convention of being the line-file
    directory itself.  Downloads into ``download_dir`` (default: a temp area under
    the data dir) and skips re-download/re-extract when the artifacts already exist
    unless ``force`` is set.
    """
    if download_dir is None:
        download_dir = paths.data_dir(create=True) / "_download"
    download_dir = Path(download_dir)
    download_dir.mkdir(parents=True, exist_ok=True)

    # First list the record's files (wget mode does not download) so we can learn
    # the tarball name without hard-coding the AER version.
    wget_list = download_dir / "line_file_list.txt"
    subprocess.check_call(
        _zenodo_cli() + ["-r", str(record), "-w", os.fspath(wget_list),
                         "-o", os.fspath(download_dir)]
    )
    urls = wget_list.read_text().splitlines()
    tarballs = [u for u in urls if u.strip().endswith(".tar.gz")]
    if not tarballs:
        raise RuntimeError(f"No .tar.gz found in Zenodo record {record}")

    tar_name = os.path.basename(tarballs[0].strip())
    tar_path = download_dir / tar_name
    extract_dir = download_dir / tar_name[:-len(".tar.gz")]

    if extract_dir.exists() and not force:
        print(f"{extract_dir} already extracted, reusing")
        return os.fspath(extract_dir)

    if not tar_path.exists() or force:
        print(f"Downloading Zenodo record {record} ({tar_name})")
        subprocess.check_call(
            _zenodo_cli() + ["-r", str(record), "-g", "*.tar.gz",
                             "-o", os.fspath(download_dir)]
        )

    print(f"Extracting {tar_path}")
    with tarfile.open(tar_path) as tar:
        tar.extractall(download_dir)
    return os.fspath(extract_dir)


def stage_line_file(extracted_dir, dest_root=None, force: bool = False) -> str:
    """Move the extracted line-file directory into place under the data dir.

    Places the contents at ``<data_dir>/AER_Line_File`` (the layout
    :func:`absco.paths.line_file_paths` expects).  Skips when the destination
    already exists unless ``force`` is set.
    """
    extracted_dir = Path(extracted_dir)
    if dest_root is None:
        dest_root = paths.line_file_root()
    dest_root = Path(dest_root)

    if dest_root.exists():
        if not force:
            print(f"{dest_root} already exists, using existing line file")
            return os.fspath(dest_root)
        shutil.rmtree(dest_root)

    dest_root.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(os.fspath(extracted_dir), os.fspath(dest_root))
    print(f"Staged line file -> {dest_root}")
    return os.fspath(dest_root)
