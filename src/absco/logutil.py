"""Logging helpers: tee the driver's console output and capture subprocess logs.

``absco-generate`` uses these to record both its own stdout/stderr and the output
of the LNFL/LBLRTM subprocess runs into a log directory (default ``<intdir>/logs``),
while still showing everything on the terminal.
"""

from __future__ import annotations

import contextlib
import os
import subprocess
import sys
from pathlib import Path

__all__ = ["resolve_log_dir", "tee_stdio", "run_logged"]


def resolve_log_dir(log_dir, intdir, create=True):
    """Return the absolute log directory, defaulting to ``<intdir>/logs``.

    ``log_dir`` overrides the default when given (relative paths are resolved
    against the current working directory). Created if ``create`` is set.
    """
    if log_dir:
        path = Path(os.path.expanduser(str(log_dir))).resolve()
    else:
        path = Path(intdir) / "logs"
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path


class _Tee:
    """File-like object that writes to several streams at once."""

    def __init__(self, *streams):
        self._streams = streams

    def write(self, data):
        for s in self._streams:
            s.write(data)
        return len(data)

    def flush(self):
        for s in self._streams:
            s.flush()


@contextlib.contextmanager
def tee_stdio(log_path):
    """Duplicate ``sys.stdout``/``sys.stderr`` to ``log_path`` and the terminal.

    Restores the original streams on exit. The log file is opened in append mode so
    repeated runs accumulate rather than clobber.
    """
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    orig_out, orig_err = sys.stdout, sys.stderr
    with open(log_path, "a", buffering=1) as fh:
        try:
            sys.stdout = _Tee(orig_out, fh)
            sys.stderr = _Tee(orig_err, fh)
            yield log_path
        finally:
            sys.stdout, sys.stderr = orig_out, orig_err


def run_logged(cmd, log_path, cwd=None, env=None, header=None, echo=True):
    """Run ``cmd`` (list) appending its stdout+stderr to ``log_path``.

    Writes an optional ``header`` line before the output. When ``echo`` is set the
    combined output is also written to the current ``sys.stdout`` (so it still shows
    on the terminal / is captured by an outer :func:`tee_stdio`). Returns the
    subprocess return code.
    """
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with open(log_path, "a", buffering=1) as fh:
        if header:
            fh.write("\n=== %s ===\n" % header)
            if echo:
                sys.stdout.write("%s\n" % header)

        proc = subprocess.Popen(
            cmd, cwd=cwd, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True,
        )
        for line in proc.stdout:
            fh.write(line)
            if echo:
                sys.stdout.write(line)
        proc.stdout.close()
        return proc.wait()
