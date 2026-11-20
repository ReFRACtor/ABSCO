"""Estimate of the peak memory used while generating an ABSCO table.

Shared by :meth:`absco.preprocess.configure.calcRAM` (the prompt printed before a
run) and :func:`absco.config_assistant.estimate_ram_gb` (the estimate printed by
``absco-config``) so the two cannot drift apart.

What is being estimated.  ``absco-generate`` builds the whole table in memory
before writing it: one table per water-vapour (and, for O2, per O2) mixing ratio,
then a stacked copy of all of them (``compute.combineVMR``), then a transposed,
masked copy that is handed to netCDF4 (``compute.makeNC``).  Several full copies of
the final table therefore exist at once.  The table itself has

    n_out_wn * n_temp * n_lay * n_vmr  float64 values,

where ``n_lay`` is two fewer than the number of pressure levels in the pressure
file (the surface level has no layer below it and the TOA level is skipped), and
``n_vmr`` is 1 for most molecules, 2 for those whose continua depend on water vapour
(two H2O values) and 4 for O2 (two H2O values for each of two O2 values).  The
wavenumber axis is stored once, not per cell, and is negligible.

Measured peaks (O2, the case with the most mixing-ratio tables; see
``plans/absco-refractor-streaming-output_plan.md`` in the research repository):

* 12,900-13,250 cm-1 at outres 1.2e-3 cm-1 (a 13.4 GB table): peak at least 53.3 GB
  (3.9 tables), from a memory sampler that stopped before the end of the run, so the
  true peak is higher;
* 12,900-13,100 cm-1 with 3 pressure levels (a 0.16 GB table): peak 0.78 GB,
  4.4 tables once the 0.08 GB of the interpreter and libraries is subtracted.

``PEAK_COPIES`` is set to 5, above both, so the estimate is not below what was
measured.  It was measured for O2 only; for the other molecules it rests on the same
sequence of operations in ``combineVMR`` and ``makeNC``, and has not been measured.
Molecules are processed one after another, but the estimate for several molecules is
the sum of the per-molecule estimates, which is conservative; whether the memory of
one molecule is returned before the next starts has not been measured.
"""

import numpy as np

#: temperatures per pressure layer in the packaged pressure/temperature grid
N_TEMP = 15

#: molecules whose continua depend on water vapour (table gets an H2O mixing-ratio axis)
MOL_H2O = ("CO2", "N2", "O2", "H2O", "HDO")

#: number of full copies of the final table held at the peak (measured, rounded up)
PEAK_COPIES = 5


def n_layers(n_pressure):
    """Layers in the output table for ``n_pressure`` levels in the pressure file."""
    return max(int(n_pressure) - 2, 1)


def n_mixing_ratios(mol):
    """Number of mixing-ratio combinations stored along the table's extra axes."""
    mol = str(mol).upper()
    if mol == "O2":
        return 4
    return 2 if mol in MOL_H2O else 1


def table_gb(n_out_wn, n_pressure, mol):
    """Size [GB] of the final float64 table of one molecule."""
    return 8.0 * float(n_out_wn) * N_TEMP * n_layers(n_pressure) * n_mixing_ratios(mol) / 1e9


def peak_gb(n_out_wn, n_pressure, mol):
    """Estimated peak memory [GB] while generating the table of one molecule."""
    return PEAK_COPIES * table_gb(n_out_wn, n_pressure, mol)


def n_output_points(wn1, wn2, outres):
    """Output spectral points summed over all bands (cm-1 inputs)."""
    wn1 = np.atleast_1d(np.asarray(wn1, dtype=float))
    wn2 = np.atleast_1d(np.asarray(wn2, dtype=float))
    outres = np.atleast_1d(np.asarray(outres, dtype=float))
    return float(np.sum((wn2 - wn1) / outres))


def estimate_gb(wn1, wn2, outres, n_pressure, molnames):
    """Return ``(per_molecule, total)`` peak estimates [GB]; ``per_molecule`` maps name -> GB."""
    n_wn = n_output_points(wn1, wn2, outres)
    per = {m: peak_gb(n_wn, n_pressure, m) for m in molnames}
    return per, float(sum(per.values()))
