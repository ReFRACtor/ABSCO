"""Shared pytest fixtures for the ABSCO test suite.

The end-to-end pipeline test needs the compiled LNFL/LBLRTM executables and the
AER line file staged in an ABSCO data directory.  Those are large/expensive to
produce, so the ``absco_data_dir`` fixture:

* uses ``$ABSCO_DATA_DIR`` (or the platformdirs default) if it already contains a
  working set of artifacts, and
* otherwise skips the integration test -- unless ``ABSCO_TEST_BUILD=1`` is set, in
  which case it builds the executables from the submodules and downloads the line
  file (the ~385 MB Zenodo archive) into a persistent cache under ``.pytest_cache``.

Fast unit tests do not use this fixture and always run.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from absco import artifacts, paths

# repo root (…/absco-refractor); tests/ lives directly under it
REPO_ROOT = Path(__file__).resolve().parent.parent


def _artifacts_present() -> bool:
    """True if both executables and the line file resolve in the current data dir."""
    return (
        paths.lnfl_exe() is not None
        and paths.lblrtm_exe() is not None
        and paths.line_file_root().exists()
    )


@pytest.fixture(scope="session")
def absco_data_dir():
    """Return a data dir populated with LNFL/LBLRTM binaries + AER line file.

    Skips the test if artifacts are absent and ABSCO_TEST_BUILD is not set.
    """
    # If the ambient data dir already has everything, just use it.
    if _artifacts_present():
        return paths.data_dir()

    if os.environ.get("ABSCO_TEST_BUILD") != "1":
        pytest.skip(
            "LNFL/LBLRTM binaries and AER line file not found. "
            "Set ABSCO_TEST_BUILD=1 to build+download them (slow, ~385 MB), "
            "or run `pixi run build-fortran` and `absco-init` first."
        )

    # Build/fetch into a persistent cache so repeated runs are cheap.
    cache = REPO_ROOT / ".pytest_cache" / "absco_data"
    os.environ[paths.DATA_DIR_ENV] = str(cache)
    cache.mkdir(parents=True, exist_ok=True)

    if paths.lnfl_exe() is None:
        artifacts.stage_executable(
            artifacts.compile_model("lnfl", REPO_ROOT / "LNFL")
        )
    if paths.lblrtm_exe() is None:
        artifacts.stage_executable(
            artifacts.compile_model("lblrtm", REPO_ROOT / "LBLRTM")
        )
        artifacts.stage_lblrtm_data_files(REPO_ROOT / "LBLRTM")

    if not paths.line_file_root().exists():
        extracted = artifacts.fetch_line_file()
        artifacts.stage_line_file(extracted)

    assert _artifacts_present(), "artifact setup did not produce a usable data dir"
    return paths.data_dir()
