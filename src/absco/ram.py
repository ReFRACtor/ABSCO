"""Estimate of the peak memory used while generating an ABSCO table.

Shared by :meth:`absco.preprocess.configure.calcRAM` (the prompt printed before a
run) and :func:`absco.config_assistant.estimate_ram_gb` (the estimate printed by
``absco-config``) so the two cannot drift apart.

What is being estimated.  ``absco-generate`` writes each LBLRTM spectrum to a scratch file
as soon as it is computed (``absco.scratch``) and ``compute.makeNC`` compresses the table
into the netCDF one block of ``freq_chunk`` frequencies at a time, so memory does not grow
with the table.  The scratch files hold the uncompressed float64 table, so the run
directory needs that much disk (``table_gb``), in addition to the netCDF.

Memory is a fixed part plus a few copies of one block:

    block = freq_chunk * n_temp * n_lay * n_vmr * 8 bytes
    peak  = BASE_GB + BLOCK_COPIES * block

``n_lay`` is two fewer than the number of pressure levels in the pressure file (the
surface level has no layer below it and the TOA level is skipped), and ``n_vmr`` is 1 for
most molecules, 2 for those whose continua depend on water vapour (two H2O values) and 4
for O2 (two H2O values for each of two O2 values).

Measured (O2, 12,900-13,250 cm-1 at outres 1.2e-3 cm-1, 98 levels, freq_chunk 171, one
core; see ``plans/absco-refractor-streaming-output_plan.md`` in the research repository):
the peak was 0.44 GB for a 13.4 GB table, sampled every 10 s over the whole run.  The
same case before the streaming change peaked at 53.3 GB or more.  ``BASE_GB`` covers the
interpreter, libraries and per-run bookkeeping; ``BLOCK_COPIES`` allows for the block as
read, its masked copy and mask, and the compressor's buffer.  Only the fixed part has been
measured; the block part rests on the sequence of operations in ``makeNC`` and has not
been tested with a large ``freq_chunk``.  Molecules are processed one after another, and
the scratch files of one are deleted before the next starts, so the largest single value is
what matters; the total shown is the sum, which is conservative.
"""

import numpy as np

#: temperatures per pressure layer in the packaged pressure/temperature grid
N_TEMP = 15

#: molecules whose continua depend on water vapour (table gets an H2O mixing-ratio axis)
MOL_H2O = ("CO2", "N2", "O2", "H2O", "HDO")

#: memory of the interpreter, libraries and bookkeeping [GB] (measured peak 0.44 GB)
BASE_GB = 0.6

#: copies of one block of frequencies held while it is written
BLOCK_COPIES = 4

#: frequency chunk of the template configuration, used when the caller does not know it
DEFAULT_FREQ_CHUNK = 5000


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


def peak_gb(n_out_wn, n_pressure, mol, freq_chunk=DEFAULT_FREQ_CHUNK):
    """Estimated peak memory [GB] while generating the table of one molecule."""
    n_f = min(float(freq_chunk), float(n_out_wn))
    block = table_gb(n_f, n_pressure, mol)
    return BASE_GB + BLOCK_COPIES * block


def n_output_points(wn1, wn2, outres):
    """Output spectral points summed over all bands (cm-1 inputs)."""
    wn1 = np.atleast_1d(np.asarray(wn1, dtype=float))
    wn2 = np.atleast_1d(np.asarray(wn2, dtype=float))
    outres = np.atleast_1d(np.asarray(outres, dtype=float))
    return float(np.sum((wn2 - wn1) / outres))


def estimate_gb(wn1, wn2, outres, n_pressure, molnames, freq_chunk=DEFAULT_FREQ_CHUNK):
    """Return ``(per_molecule, total)`` peak memory estimates [GB]; ``per_molecule`` maps name -> GB."""
    n_wn = n_output_points(wn1, wn2, outres)
    per = {m: peak_gb(n_wn, n_pressure, m, freq_chunk) for m in molnames}
    return per, float(sum(per.values()))


def disk_gb(wn1, wn2, outres, n_pressure, molnames):
    """Scratch disk [GB] needed in the run directory: the largest single table (the scratch
    files of one molecule are deleted before the next)."""
    n_wn = n_output_points(wn1, wn2, outres)
    return float(max(table_gb(n_wn, n_pressure, m) for m in molnames))
