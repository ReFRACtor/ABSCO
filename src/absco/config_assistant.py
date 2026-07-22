"""Build an ``ABSCO_config.ini`` for the generation step.

The user typically specifies only the spectral range, the desired output
resolution, and the molecules; this module suggests a matching LBLRTM computation
resolution (``lblres``) and writes a config from the bundled template, leaving all
data/executable/line-file paths blank so :mod:`absco.paths` resolves them at run
time.

The ``lblres`` suggestion mirrors ``notebooks/create_absco_config-tropomi.ipynb``:
target ``lblres <= TARGET_LBLRES`` while keeping ``outres/lblres`` an exact power of
2 (a hard requirement enforced in :mod:`absco.preprocess`).

Inputs may be given in wavenumbers (cm-1) or wavelengths (um/nm). Because LBLRTM and
the degradation kernel work on a wavenumber grid, wavelength inputs are converted to
cm-1 (see :func:`convert_bands_to_cm1`) and the config is always written in cm-1.
"""

from __future__ import annotations

import configparser
from configparser import ConfigParser

import numpy as np

from absco import paths

__all__ = [
    "TARGET_LBLRES",
    "suggest_lblres",
    "convert_bands_to_cm1",
    "build_config",
    "write_config",
    "estimate_ram_gb",
    "WN_FORMAT",
    "RES_FORMAT",
]

# Multiplicative constants to convert a wavelength [um or nm] to wavenumber [cm-1]
# via wavenumber = CONST / wavelength.
_WL_TO_CM1 = {"um": 1e4, "nm": 1e7}

# One-line description per config field, drawn from the README configuration-field
# table, emitted as a comment above each key by write_config().
FIELD_DOC = {
    ("data files", "header"): "80-character header written into each final OD file (coarse grid)",
    ("data files", "pfile"): "text file with one pressure level [mbar] per line; the pressures on which ABSCOs are calculated",
    ("data files", "ptfile"): "text file mapping each pressure level to its allowed layer-average temperatures [K]",
    ("data files", "vmrfile"): "CSV of interpolated/extrapolated volume mixing ratios (VMRs) for the full profile",
    ("data files", "hdofile"): "CSV HDO profile (used when HDO is among the molecules)",
    ("data files", "xs_lines"): "CSV of species/bands where HITRAN recommends line parameters over cross sections",
    ("channels", "wn1"): "starting spectral point(s) for each band (in 'units'); space-delimited",
    ("channels", "wn2"): "ending spectral point(s) for each band (in 'units'); space-delimited",
    ("channels", "lblres"): "LBLRTM computation resolution(s) [cm-1]; outres/lblres must be a power of 2",
    ("channels", "outres"): "output (degraded) resolution(s) [cm-1] after spectral degradation",
    ("channels", "units"): "spectral units for wn1/wn2: cm-1, um, or nm",
    ("vmr", "wv_vmr"): "two water-vapor VMR values [ppmv] for H2O/CO2/O2/N2 (continua depend on water vapor)",
    ("molecules", "molnames"): "HITRAN molecule names, space-delimited, case-insensitive",
    ("makeTAPE5", "scale"): "continuum/extinction scaling factor(s)",
    ("makeTAPE5", "tape5_dir"): "subdirectory (under the LNFL/LBL run dirs) for generated TAPE5 files",
    ("runLNFL", "lnfl_run_dir"): "directory (under intdir) where LNFL runs occur",
    ("runLNFL", "tape1_path"): "TAPE1 ASCII line file used in LNFL runs",
    ("runLNFL", "tape2_path"): "TAPE2 ASCII line-coupling file used in LNFL runs (O2, CO2, CH4)",
    ("runLNFL", "lnfl_path"): "LNFL executable",
    ("runLNFL", "extra_params"): "directory of broadening / speed-dependence parameter files",
    ("runLNFL", "tape3_dir"): "directory (under intdir) for LNFL output binary line files (TAPE3)",
    ("runLBL", "lbl_path"): "LBLRTM executable",
    ("runLBL", "xs_path"): "LBLRTM cross-section file directory",
    ("runLBL", "fscdxs"): "cross-section lookup file used with xs_path",
    ("runLBL", "lbl_run_dir"): "directory (under intdir) where LBLRTM runs occur",
    ("output", "intdir"): "top-level directory for intermediate and output files ('.' = current working directory)",
    ("output", "outdir"): "directory (under intdir) where output netCDFs are written",
    ("output", "sw_ver"): "software version string recorded in the output netCDF",
    ("output", "out_file_desc"): "run description used in the output file name (no spaces)",
    ("output", "nc_compress"): "netCDF compression level (0-9)",
    ("output", "freq_chunk"): "chunk size for the frequency dimension of the cross-section dataset",
}

# Upper bound for the LBLRTM computation resolution [cm-1]. LBLRTM wants a fine
# grid (~1.5e-4 or smaller); see the notebook and ABSCO docs.
TARGET_LBLRES = 1.5e-4

# String formats used when writing the .ini, matching the notebook so the values
# round-trip cleanly (and the power-of-2 ratio survives the text representation).
WN_FORMAT = "{:.4f}"
RES_FORMAT = "{:0.10e}"

# Molecules whose continua depend on water vapor -> table gets an extra VMR
# dimension (doubles the RAM estimate). Mirrors preprocess.configure.molH2O.
_MOL_H2O = ["CO2", "N2", "O2", "H2O", "HDO"]


def suggest_lblres(outres):
    """Suggest an LBLRTM resolution for each output resolution.

    For each ``outres`` value, returns ``lblres = outres / 2**power`` where
    ``power = ceil(log2(outres / TARGET_LBLRES))`` (clamped so ``power >= 0``, i.e.
    ``lblres`` is never coarser than ``outres``). This guarantees ``outres/lblres``
    is an exact power of 2 and ``lblres <= TARGET_LBLRES`` whenever ``outres`` is
    coarser than the target.

    Accepts a scalar or array-like; always returns a 1-D numpy array.
    """
    outres = np.atleast_1d(np.asarray(outres, dtype=float))
    if np.any(outres <= 0):
        raise ValueError("outres values must be positive")

    power = np.ceil(np.log2(outres / TARGET_LBLRES)).astype(int)
    power = np.maximum(power, 0)
    return outres / 2.0 ** power


def _verify_power_of_two(outres, lblres):
    """Assert outres/lblres is an exact power of 2 after .ini string formatting.

    Mirrors the notebook's round-trip check: format each value the way it will be
    written to the file, parse it back, and confirm log2(outres/lblres) is integral.
    Raises ValueError if not (which would make preprocess.configure reject the file).
    """
    for o, l in zip(np.atleast_1d(outres), np.atleast_1d(lblres)):
        o_rt = float(RES_FORMAT.format(o))
        l_rt = float(RES_FORMAT.format(l))
        ratio = o_rt / l_rt
        power = np.log2(ratio)
        if not np.isclose(power, round(power), atol=1e-9):
            raise ValueError(
                "outres/lblres = %g/%g = %g is not a power of 2 after formatting"
                % (o_rt, l_rt, ratio)
            )


def estimate_ram_gb(wn1, wn2, outres, molnames):
    """Best-effort RAM estimate [GB] for the run, mirroring preprocess.calcRAM.

    Uses the packaged pressure grid for the layer count. Returns the total across
    all molecules (H2O-affected molecules counted twice for their extra VMR
    dimension). Returns None if the pressure grid cannot be read.
    """
    try:
        pressures = np.loadtxt(paths.default_data_file("pfile"))
        n_p = int(np.atleast_1d(pressures).size)
    except Exception:
        return None

    wn1 = np.atleast_1d(np.asarray(wn1, dtype=float))
    wn2 = np.atleast_1d(np.asarray(wn2, dtype=float))
    outres = np.atleast_1d(np.asarray(outres, dtype=float))

    # number of output spectral points across all bands
    n_out_wn = float(np.sum((wn2 - wn1) / outres))

    # 16 bytes (wavenumber + absco, float64) x nOutWN x 15 temperatures x nLayers
    per_mol_gb = 16 * 2 * n_out_wn * 15 * n_p / 1e9

    total = 0.0
    for mol in molnames:
        scale = 2 if mol.upper() in _MOL_H2O else 1
        total += scale * per_mol_gb
    return total


def convert_bands_to_cm1(wn1, wn2, outres, units):
    """Convert band bounds and output resolution from ``units`` to cm-1.

    ``units`` is ``cm-1`` (returned unchanged), ``um``, or ``nm``. For wavelength
    units the bounds are converted with ``wavenumber = CONST / wavelength`` (and
    reordered so wn1 < wn2), and each band's ``outres`` -- a *spacing* expressed in
    the same wavelength units -- is converted to a cm-1 spacing at the band-center
    wavelength using ``|d(wavenumber)| = CONST / wavelength**2 * d(wavelength)``.

    Returns ``(wn1_cm1, wn2_cm1, outres_cm1)`` as 1-D numpy arrays.
    """
    wn1 = np.atleast_1d(np.asarray(wn1, dtype=float))
    wn2 = np.atleast_1d(np.asarray(wn2, dtype=float))
    outres = np.atleast_1d(np.asarray(outres, dtype=float))

    units = str(units).lower()
    if units == "cm-1":
        return wn1, wn2, outres
    if units not in _WL_TO_CM1:
        raise ValueError("units must be 'cm-1', 'um', or 'nm' (got %r)" % units)

    const = _WL_TO_CM1[units]
    wl_lo = np.minimum(wn1, wn2)
    wl_hi = np.maximum(wn1, wn2)

    # bounds: wavenumber = const / wavelength (larger wavelength -> smaller wavenumber)
    v1 = const / wl_hi
    v2 = const / wl_lo

    # spacing at band center: |dv| = const / wl_center**2 * d(wavelength)
    wl_center = 0.5 * (wl_lo + wl_hi)
    outres_cm1 = const / wl_center ** 2 * outres

    return v1, v2, outres_cm1


def build_config(wn1, wn2, outres, molnames, units="cm-1", wv_vmr=None,
                 lblres=None, overrides=None):
    """Return a ConfigParser for ABSCO_config.ini built from the bundled template.

    Fills ``[channels]`` (wn1/wn2/lblres/outres/units) and ``[molecules] molnames``;
    optionally sets ``[vmr] wv_vmr`` and any ``overrides`` (a dict of
    ``(section, key) -> value`` for custom data-file paths etc.). All path fields
    left in the template stay blank so absco.paths resolves them at run time.

    ``units`` is the units of the *input* ``wn1``/``wn2``/``outres`` (``cm-1``,
    ``um``, or ``nm``). Wavelength inputs are converted to cm-1 via
    :func:`convert_bands_to_cm1`, and the written config is ALWAYS in cm-1 with
    ``units = cm-1`` -- LBLRTM and the degradation kernel operate on a wavenumber
    grid, so ``lblres``/``outres`` must be cm-1 spacings.

    ``lblres`` defaults to :func:`suggest_lblres` of the (converted) ``outres`` when
    not given; if provided explicitly it is interpreted as a cm-1 spacing. Raises
    ValueError unless every band's outres/lblres is an exact power of 2.
    """
    wn1 = np.atleast_1d(np.asarray(wn1, dtype=float))
    wn2 = np.atleast_1d(np.asarray(wn2, dtype=float))
    outres = np.atleast_1d(np.asarray(outres, dtype=float))

    if not (wn1.size == wn2.size == outres.size):
        raise ValueError("wn1, wn2, and outres must have the same number of bands")

    orig_units = str(units).lower()
    # convert bounds + output spacing to cm-1; the config is always written in cm-1
    wn1, wn2, outres = convert_bands_to_cm1(wn1, wn2, outres, orig_units)

    if lblres is None:
        lblres = suggest_lblres(outres)
    else:
        lblres = np.atleast_1d(np.asarray(lblres, dtype=float))
        if lblres.size != outres.size:
            raise ValueError("lblres must have the same number of bands as outres")

    _verify_power_of_two(outres, lblres)

    config = ConfigParser()
    config.read(paths.config_template())

    config["channels"]["wn1"] = " ".join(WN_FORMAT.format(v) for v in wn1)
    config["channels"]["wn2"] = " ".join(WN_FORMAT.format(v) for v in wn2)
    config["channels"]["lblres"] = " ".join(RES_FORMAT.format(v) for v in lblres)
    config["channels"]["outres"] = " ".join(RES_FORMAT.format(v) for v in outres)
    config["channels"]["units"] = "cm-1"

    config["molecules"]["molnames"] = " ".join(m.lower() for m in molnames)

    # record the installed package version so the output netCDF's sw_ver matches
    from absco import __version__ as absco_version
    config["output"]["sw_ver"] = str(absco_version)

    if wv_vmr is not None:
        config["vmr"]["wv_vmr"] = " ".join(str(v) for v in np.atleast_1d(wv_vmr))

    if overrides:
        for (section, key), value in overrides.items():
            if not config.has_section(section):
                config.add_section(section)
            config[section][key] = str(value)

    return config


def _resolved_default(section, key):
    """Return the path a blank config field will resolve to at run time, or None.

    Mirrors the resolution in absco.paths so the generated config can annotate each
    blank path field with the bundled/data-dir file that will actually be used.
    """
    data_keys = {
        "pfile": "pfile",
        "ptfile": "ptfile",
        "vmrfile": "vmrfile",
        "hdofile": "hdofile",
        "xs_lines": "xs_lines",
    }
    if section == "data files" and key in data_keys:
        try:
            return paths.default_data_file(data_keys[key])
        except Exception:
            return None
    if section == "runLNFL" and key == "lnfl_path":
        return paths.lnfl_exe()
    if section == "runLBL" and key == "lbl_path":
        return paths.lblrtm_exe()
    line_defaults = None
    if (section, key) in {
        ("runLNFL", "tape1_path"), ("runLNFL", "tape2_path"),
        ("runLNFL", "extra_params"), ("runLBL", "xs_path"), ("runLBL", "fscdxs"),
    }:
        line_defaults = paths.line_file_paths()
        return line_defaults.get(key)
    return None


def write_config(config, fh, header_note=None):
    """Write ``config`` (a ConfigParser) to file object ``fh`` with comments.

    Each field is preceded by a one-line description (from :data:`FIELD_DOC`), and
    any blank path field gets an extra comment naming the packaged/data-dir file
    that will be used at run time (or a note that the artifact is not yet available).

    ``header_note`` is an optional string (or list of lines) written as a comment
    block at the top of the file -- used to record e.g. that wavelength inputs were
    converted to cm-1.
    """
    if header_note:
        lines = [header_note] if isinstance(header_note, str) else list(header_note)
        for line in lines:
            fh.write("; %s\n" % line)
        fh.write("\n")

    for section in config.sections():
        fh.write("[%s]\n" % section)
        for key, value in config.items(section):
            doc = FIELD_DOC.get((section, key))
            if doc:
                fh.write("; %s\n" % doc)
            if value.strip() == "":
                default = _resolved_default(section, key)
                if default is not None:
                    fh.write("; (blank -> using packaged/default: %s)\n" % default)
                elif (section, key) in {("runLNFL", "lnfl_path"), ("runLBL", "lbl_path")}:
                    fh.write("; (blank -> resolved from the data dir; run absco-build or absco-init)\n")
            fh.write("%s = %s\n" % (key, value))
            # blank line after each parameter for readability
            fh.write("\n")
