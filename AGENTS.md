# AGENTS.md

This file provides guidance to AI agents when working with code in this repository.

## Overview

ABSCO (ABSoprtion COefficient) is a table generation system for atmospheric radiative transfer applications. It computes absorption coefficients across thermal IR to UV spectral ranges using LBLRTM (Line-By-Line Radiative Transfer Model) and LNFL (LiNe FiLe) Fortran codes with the AER line parameter database.

The software generates netCDF tables of absorption coefficients indexed by wavenumber, pressure, temperature, and water vapor VMR (for H₂O-affected molecules). These tables support remote sensing instruments and atmospheric retrieval algorithms.

## Architecture

### Core Pipeline
1. **Configuration** (`ABSCO_preprocess.py`) - Parse `ABSCO_config.ini`, validate inputs, check spectral ranges, verify molecule activity, compute spectral degradation kernels
2. **Line File Generation** - LNFL converts ASCII line parameters to binary TAPE3 files required by LBLRTM
3. **Optical Depth Calculation** (`ABSCO_compute.py`) - LBLRTM computes optical depths at specified P/T/spectral points
4. **Absorption Coefficient Calculation** - Divide optical depths by layer amounts to produce cross sections
5. **Output** - Write multi-dimensional netCDF with dimensions: [nfreq × ntemp × nlay × nvmr]

### Key Design Patterns
- **H₂O-affected molecules** (H₂O, CO₂, O₂, N₂) require dual VMR calculations because their continua depend on water vapor amount
- **TAPE3 caching** - Binary line files are generated once per molecule/band and reused across P/T calculations
- **Spectral degradation** - LBLRTM runs at high resolution (e.g., 1.5e-4 cm⁻¹) then degrades output to reduce storage (e.g., 1.2e-3 cm⁻¹)
- **XS vs line parameters** - Some molecules have both cross sections and line parameters; `FSCDXS_line_params.csv` encodes HITRAN recommendations

### Directory Structure
- `src/absco/` - Installable Python package (src layout)
  - `src/absco/_common/` - Shared utilities vendored from AER-RC/common: `lblTools.py` (LBLRTM I/O), `RC_utils.py`, `utils.py`, `FortranFile.py`. These are vendored (copied in), not a submodule, so the installed package has no runtime dependence on the `common` repository.
  - `src/absco/data/` - Bundled data files resolved from the installed package (see `src/absco/paths.py`): `PT_grid/` (default AIRS instrument grid), `VMR/` profiles, `FSCDXS_line_params.csv`, and `ABSCO_config.template.ini`
- `LBLRTM/` - Fortran source code submodule (AER-RC/LBLRTM v12.17). Has nested submodules `aer_rt_utils` and `cross-sections` — requires `--recursive` init.
- `LNFL/` - Fortran line file converter submodule (AER-RC/LNFL master, v3.2-30). Has nested submodule `aer_rt_utils` — requires `--recursive` init.
- `VMR/` - Volume mixing ratio profile generator scripts
- `AER_Line_File/` - Line parameter database (v3.7, downloaded via Zenodo by `absco-init`)

## Common Commands

### Initial Setup
```bash
# Clone with submodules
git clone --recursive git@github.com:ReFRACtor/ABSCO.git

# Build Fortran executables and download line file
./build_models.py -c gfortran -i ABSCO_config.ini

# Setup environment
conda env create -n absco -f environment.yml
conda activate absco
```

### Running ABSCO Generation

```bash
# Full end-to-end run (LNFL + LBLRTM + netCDF output)
./run_LBLRTM_ABSCO.py -e2e

# Generate binary line files only (TAPE3)
./run_LBLRTM_ABSCO.py -lnfl

# Run LBLRTM and generate netCDF (assumes TAPE3 exists)
./run_LBLRTM_ABSCO.py -lbl

# Specify custom config file
./run_LBLRTM_ABSCO.py -e2e -i custom_config.ini
```

### Utilities

```bash
# Read absorption coefficient from output netCDF
./read_ABSCO_tables.py nc_ABSCO/output.nc -p 500 -T 250 -s 800 cm-1 -wv 10000

# Split config for parallel processing
./split_config.py -i ABSCO_config.ini -n 4

# Join multiple ABSCO tables
./join_tables.py -i table1.nc table2.nc -o combined.nc
```

## Configuration (`ABSCO_config.ini`)

Key parameters to modify for new runs:

- **Spectral range**: `wn1`, `wn2` (minimum 200 cm⁻¹ extent recommended for accuracy)
- **Resolution**: `lblres` (LBLRTM resolution), `outres` (degraded output resolution)
- **Molecules**: `molnames` (case-insensitive HITRAN names; leave empty for auto-detection)
- **Water vapor VMR**: `wv_vmr` (space-delimited ppmv values, e.g., `1.0e1 4.0e4`)
- **Pressure grid**: `pfile` (default: `PT_grid/AIRS_P_air.txt`)
- **Output directory**: `outdir` (relative to `intdir`)

Critical: Run `./build_models.py -c <compiler> -i ABSCO_config.ini` after cloning to auto-populate paths for executables and line files.

## Allowed Molecules

H₂O, CO₂, O₃, N₂O, CO, CH₄, O₂, NO, SO₂, NO₂, NH₃, HNO₃, OCS, H₂CO, N₂, HCN, C₂H₂, HCOOH, C₂H₄, CH₃OH, CCL4, CF4, F11, F12, F22, ISOP, PAN, HDO, BRO, O₂-O₂

For HDO: set `tape1_path` to `01_h2o_162_only` subdirectory and use `hdofile` for VMR profile.

## Important Notes

- **Python 3 only** - Python 2 not tested or supported
- **Compiler options**: `gfortran`, `ifort`, `pgf90` (gfortran pinned to 11.2.0 in pixi for building the legacy Fortran; requires `netcdf-fortran` in the env for LBLRTM v12.17)
- **Model versions**: LNFL master (v3.2-30), LBLRTM v12.17, AER LPD v3.7
- **H₂O-affected molecules output extra dimension**: netCDF includes `H2O_VMR` dimension for H₂O, CO₂, O₂, N₂
- **O₂ special case**: Uses two O₂ VMRs (1.9e5, 2.3e5 ppmv) in addition to two H₂O VMRs
- **TAPE3 existence check**: LNFL skips regeneration if TAPE3 already exists for molecule/band
- **Memory/disk usage**: Large files possible; software estimates size before run
- **Fill values**: Invalid P/T combinations populated with NaN in netCDF

## Submodule Updates

```bash
# Update all submodules to latest commits
git submodule update --init --recursive

# Update specific submodule
cd LBLRTM && git pull origin master
```

## VMR Profile Generation

To use non-default atmospheric profiles:

```bash
cd VMR
python standard_atm_profiles.py  # generates CSV with interpolated VMRs
# Set output filename as `vmrfile` in ABSCO_config.ini
```

For HDO profiles:
```bash
python standard_atm_HDO_profile.py
# Set output as `hdofile` in config
```
