"""ABSCO: absorption coefficient table generation using LBLRTM and LNFL.

This package wraps the ABSCO pipeline (preprocessing, LNFL/LBLRTM execution, and
netCDF output) as an installable tool that can be run from an arbitrary working
directory.  See :mod:`absco.paths` for how bundled data and runtime artifacts
(line file, Fortran executables) are located.
"""

__version__ = "2.0"
