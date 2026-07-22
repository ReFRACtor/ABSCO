# ABSoprtion COefficient (ABSCO) Table Generation

---
**Contents**

1. [Quickstart](#quickstart)
2. [Introduction](#intro)
3. [Dependencies](#dependencies)
4. [Setup](#setup)
5. [HITRAN Cross Sections](#xs)
6. [Driver Script](#driver)
7. [Parallel Processing](#parallel)
8. [Output](#output)

**Tables**

[AER-maintained Models](#Table1)<br>
[Configuration File Fields](#Table2)<br>
[Output netCDF Dimensions](#Table3)<br>
[Output netCDF Variables](#Table4)<br>

# Quickstart <a name="quickstart"></a>

ABSCO is packaged as an installable Python tool (`absco`) with console commands:

| Command | Purpose |
| :--- | :--- |
| `absco-build` | (developer) compile the LNFL/LBLRTM executables from the submodules and stage them into the data directory |
| `absco-init` | (end user) download the AER line file and stage prebuilt binaries into the data directory |
| `absco-config` | write an `ABSCO_config.ini` from a spectral range, output resolution, and molecule list (suggests `lblres` automatically) |
| `absco-generate` | run the LNFL → LBLRTM → netCDF pipeline |
| `absco-read` | read a coefficient out of an output netCDF |

The tool is designed to be installed once and then run from any working directory;
all intermediate/temporary files are written under the current directory, and the
bundled data files and runtime artifacts (executables, line file) are located
automatically (see [Setup](#setup)).

## With pixi (recommended for development)

One-time setup from the repository checkout:

```
git clone --recursive git@github.com:ReFRACtor/ABSCO.git
cd ABSCO
pixi install                 # creates the environment and installs absco (editable)
pixi run build-fortran       # init submodules + compile LNFL/LBLRTM into the data dir
pixi run init                # download + stage the AER line file
```

Then activate the environment with `pixi shell` and run the `absco-*` commands from
any working directory (they are on `PATH` inside the shell):

```
pixi shell                   # activate the environment (run from the ABSCO checkout)
mkdir -p ~/absco_work && cd ~/absco_work
absco-config --outres 0.01 --begin 4166 --end 4358 --molnames ch4 co h2o
absco-generate -e2e -y
absco-read nc_ABSCO/*.nc -p 500 -T 250 -s 4200 cm-1 -wv 10 -tol 0.05
```

## With pip

```
pip install absco            # (prebuilt binary wheel; or `pip install -e .` from a checkout)
absco-init                   # download the AER line file (+ stage bundled binaries)
cd my_work_dir
absco-config --outres 0.01 --begin 4166 --end 4358 --molnames ch4 co h2o
absco-generate -e2e -y
absco-read nc_ABSCO/*.nc -p 500 -T 250 -s 4200 cm-1 -wv 10 -tol 0.05
```

> **Note on the Fortran build:** LNFL and LBLRTM are legacy AER Fortran. `absco-build`
> shims the compiler with `-std=legacy -fallow-argument-mismatch` so a current
> `gfortran` builds them, and the environment must provide `netcdf-fortran` (LBLRTM
> v12.17 reads an MT_CKD continuum netCDF at run time). The pixi environment sets this
> up for you.

# Introduction <a name="intro"></a>

Software that can generate a set of ABSCO tables that span the thermal IR to UV spectral range.

Assuming the user has set the `user.name` and `user.email` Git configuration variables and has setup their Github account SSH keys, the repository can be cloned by using:

```
git clone --recursive git@github.com:ReFRACtor/ABSCO.git
```

Note the `--recursive` keyword -- it will force the clone to copy the necessary subroutines. Otherwise, the user will have to do (assuming the repository is just checked out with the default `ABSCO` name in the working directory):

```
cd ABSCO
git submodule update --init --recursive
```

# Dependencies <a name="dependencies"></a>

This library depends on a number of standard Python libraries (`os`, `sys`, `configparser`, `subprocess`, `glob`, `argparse`), some widely distributed third-party libraries (see "Python Packages" subsection), and some ad hoc subroutines that originate from the [AER-RC/common](<https://github.com/AER-RC/common>) repository (`utils`, `RC_utils`, `lblTools`, `FortranFile`). These modules are vendored into the package as `absco._common`, so the installed tool has no runtime dependence on the `common` repository. Additionally, some AER-maintained models and line parameters are required to run the ABSCO software.

## Python Packages

Runtime dependencies are declared in `pyproject.toml` (and mirrored in `environment.yml`
and `pixi.toml`) and installed automatically with the package:

  - Python ≥ 3.10
  - numpy ≥ 1.24
  - scipy ≥ 1.10
  - pandas ≥ 2.0
  - netCDF4 ≥ 1.6
  - xarray ≥ 2023.1
  - platformdirs ≥ 3
  - [`zenodo_get`](https://gitlab.com/dvolgyes/zenodo_get) ≥ 1.3 (fetches the line file from Zenodo)

Building the Fortran executables additionally requires a Fortran compiler (`gfortran`),
`make`, and `netcdf-fortran`; the pixi environment provides these. **Python 3 only --
Python 2 is not supported.**

## LNFL, LBLRTM, and the AER Line File

LNFL (LiNe FiLe) FORTRAN code that converts ASCII text line parameter files to the binary files that LBLRTM expects (TAPE3) is located in [its own repository](<https://github.com/AER-RC/LNFL>) and is a submodule under the `LNFL` subdirectory. LBLRTM (Line-By-Line Radiative Transfer Model) FORTRAN code also has [its own Git repository](<https://github.com/AER-RC/LBLRTM>) and is a submodule under `LBLRTM`; in the context of this software it calculates optical depth at specified pressures, temperatures, and spectral ranges. Both of these submodules themselves contain nested submodules (`aer_rt_utils`, and for LBLRTM `cross-sections`), so they must be cloned **recursively** (`git clone --recursive`, or `git submodule update --init --recursive`).

The AER line parameter database (LPD) is distributed as a set of ASCII text files archived as a [dataset in Zenodo](https://zenodo.org/record/18881607) (record `18881607`, `aer_v_3.9.tgz`) rather than in Git; `absco-init` (and `absco-build --lines`) fetch and extract it with the `zenodo_get` library and stage it under the data directory. A different record can be selected with `--record`.

Periodically, the models and LPD are updated to reflect new line parameters, a new continuum, or bug fixes. These revisions can have significant effects on the model output. The versions currently pinned as submodules are listed in [Table 1](#Table1).

**AER-maintained models pinned as ABSCO submodules** <a id="Table1"></a>

| Model | Version (GitHub Tag/Release) |
| :---: | :---: |
| LNFL | master (v3.2-30) |
| LBLRTM | v12.17 |
| LPD | v3.9 |

The executables are built and staged with `absco-build` (developer, compiles from the
submodules) or downloaded/staged with `absco-init` (end user, prebuilt binary wheel):

```
% pixi run build-fortran      # = git submodule update --init --recursive; absco-build
```

`absco-build` compiles LNFL (single precision) and LBLRTM (double precision) for the
current OS/compiler, then copies the resulting executables into the data directory as
`<data_dir>/bin/model_version_os_compiler_precision` (e.g. `lblrtm_v12.17_linux_gnu_dbl`).
It also stages the LBLRTM MT_CKD continuum netCDF that v12.17 needs at run time, and
the LBLRTM cross-section files (`xs/` and `FSCDXS`) -- as of AER line file v3.9 the
line-file archive no longer ships an `xs_files/` directory, so the cross sections come
from LBLRTM's `cross-sections` submodule instead. The executables, line file, and
cross sections are located at run time (see [Setup](#setup)), so **no paths need to be
written into `ABSCO_config.ini`**. Supported compilers are `gfortran` (default),
`ifort`, and `pgf90`; use `-h` for options.

The runtime artifacts live in a data directory resolved in this order: the
`ABSCO_DATA_DIR` environment variable, otherwise a per-user location from
[`platformdirs`](https://pypi.org/project/platformdirs/) (e.g. `~/.local/share/absco`).
Existing artifacts are not overwritten unless `--force` is given.

# Setup (Configuration File) <a name="setup"></a>

`ABSCO_config.ini` contains the inputs expected from the user. **The easiest way to
create one is `absco-config`**, which needs only the spectral range, output resolution,
and molecules and fills in the rest (including a suggested `lblres`) from the bundled
template:

```
absco-config --outres 0.01 --begin 4166 --end 4358 --molnames ch4 co h2o
```

`absco-config` writes the data/executable/line-file path fields **blank**; ABSCO then
resolves them at run time -- bundled data (`pfile`, `ptfile`, `vmrfile`, `hdofile`,
`xs_lines`) from the installed package, and the executables and line-file components
(`lnfl_path`, `lbl_path`, `tape1_path`, `tape2_path`, `extra_params`, `xs_path`,
`fscdxs`) from the data directory populated by `absco-init`/`absco-build`. You only need
to set one of these fields if you want to override a default with a custom file. The
`--lblres` computation and the required power-of-2 `outres/lblres` ratio are handled
automatically.

### Wavelength units

LBLRTM and the spectral-degradation kernel operate on a wavenumber (cm<sup>-1</sup>)
grid, so `lblres`/`outres` must be cm<sup>-1</sup> spacings. If you specify `--units um`
or `--units nm`, `absco-config` treats `--begin`/`--end` as wavelengths and `--outres` as a
wavelength spacing (in those same units), then **converts everything to cm<sup>-1</sup>
before writing the file**:

- band bounds convert via `wavenumber = C / wavelength` (C = 10<sup>4</sup> for um,
  10<sup>7</sup> for nm) and are reordered so `wn1 < wn2`;
- `outres` is converted to a cm<sup>-1</sup> spacing at each band's center wavelength
  (`|dν| = C / λ² · dλ`);
- `lblres` is then suggested from the converted cm<sup>-1</sup> `outres`.

The generated config always has `units = cm-1` (with a header comment noting the
original units), so it is consistent with what LBLRTM expects. An explicit `--lblres`
is always interpreted as a cm<sup>-1</sup> spacing.

For reference, the full set of parameters recognized in the file is given in
[Table 2](#Table2). Any field left blank falls back to the resolved default described
above.

**ABSCO configuration file (`ABSCO_config.ini` by default) items** <a id="Table2"></a>

| Field | Parent Directory | Notes |
| :---: | :---: | :--- |
| pfile | packaged / `PT_grid` | text file with 1 pressure level in millibars per line. these will be the pressures on which the ABSCOs are calculated. leave blank to use the packaged default, or set an absolute path to a custom file |
| ptfile | packaged / `PT_grid` | for every pressure level, there are different allowed temperatures. this file contains a set of pressures and their permitted temperatures. leave blank to use the packaged default, or set an absolute path |
| vmrfile | packaged / `VMR` | CSV file (see `VMR/standard_atm_profiles.py`) providing interpolated/extrapolated volume mixing ratios (VMRs) for the entire profile. leave blank to use the packaged default, or set an absolute path |
| xs_lines | packaged | this file contains the species names for when XS and line parameters exist and line parameter usage is recommended by HITRAN. leave blank to use the packaged default, or set an absolute path |
| wn1, wn2 | N/A | starting and ending spectral points for every desired band. can be in wavenumbers, microns, or nanometers. 200 cm<sup>-1</sup> should be the minimum extent of the window -- broader windows increase the accuracy of the calculations |
| lblres, outres | N/A | spectral resolution at which LBLRTM is run and spectral resolution of the output (after spectral degradation). for now, this should be in wavenumbers |
| units | N/A | spectral units ("cm<sup>-1</sup>", "um", and "nm" are allowed) |
| wv_vmr | N/A | water vapor VMR (in **ppmv**) that will be used for all levels in a given profile for H<sub>2</sub>O, CO<sub>2</sub>, O<sub>2</sub>, and N<sub>2</sub>, since their continua are dependent on water vapor amounts |
| molnames | N/A | HITRAN molecule names, space delimited, case-insensitive. really should only be one molecule per run |
| scale | N/A | multiplicative factors used for continuum or extinction scaling (separate factors for H<sub>2</sub>O self continuum, H<sub>2</sub>O foreign continuum, CO<sub>2</sub> continuum, O<sub>3</sub> continuum, O<sub>2</sub> continuum, N<sub>2</sub> continuum, and Rayleigh extinction) |
| tape5_dir | LNFL_Runs and LBL_Runs | Directory underneath both lnfl_run_dir and lbl_run_dir to which TAPE5 files will be written (should just be a single string) |
| lnfl_run_dir | `intdir` | Directory where LNFL runs will occur (created under `intdir`, one subdirectory per molecule). |
| tape1_path | data dir | Full path to the TAPE1 ASCII line file used in LNFL runs. leave blank to resolve from the data dir (`AER_Line_File/line_file/aer_v_*`) |
| tape2_path | data dir | Full path to the TAPE2 ASCII line coupling file used in LNFL runs with O<sub>2</sub>, CO<sub>2</sub>, and CH<sub>4</sub>. leave blank to resolve from the data dir |
| lnfl_path | data dir | Full path to the LNFL executable. leave blank to resolve from the data dir (`bin/`) |
| extra_params | data dir | directory with CO<sub>2</sub>-CO<sub>2</sub>, CO<sub>2</sub>-H<sub>2</sub>O, O<sub>2</sub>-H<sub>2</sub>O, O<sub>2</sub>-O<sub>2</sub>, and H<sub>2</sub>O-CO<sub>2</sub> broadening parameter specifications. leave blank to resolve from the data dir |
| tape3_dir | `intdir` | directory (under `intdir`) to which LNFL output (binary line files) will be written |
| lbl_path | data dir | Full path to the LBLRTM executable. leave blank to resolve from the data dir (`bin/`) |
| xs_path | data dir | Full path to the LBLRTM cross section file directory. leave blank to resolve from the data dir |
| fscdxs | data dir | Full path to the cross section "lookup" file used with xs_path. leave blank to resolve from the data dir |
| lbl_run_dir | `intdir` | Directory where LBL runs will occur (created under `intdir`, one subdirectory per molecule). |
| intdir | N/A | Directory where intermediate and output directories (`lnfl_run_dir`, `lbl_run_dir`, `tape3_dir`, and `outdir`) will be written.|
| outdir | `intdir` | directory relative to `intdir` where the netCDFs will be written. |
| sw_ver | N/A | used in the output netCDF file to specify the software version number |
| out_file_desc | N/A | used in the output netCDF file. allows use to document run more specifically |
| nc_compress | N/A | compression level for netCDF output |

`ABSCO_config.ini` can be named something else, but that will have to be specified at the command line (otherwise it's assumed that `ABSCO_config.ini` in the working directory is the configuration file to use):

```
absco-generate -i your_ini_file -e2e
```

# HITRAN Cross Section Usage <a name="xs"></a>

Some molecules have both line parameters and XS parameters.  HITRAN makes recommendations on the preferred parameters given the species and the band, and these are taken into account in the error checking that the `absco.preprocess` module does.  Molecules where line parameters are recommended, the associated bands, and the flag (0: only XS exist, 1: both exist, use line params, 2: both exist, use XS) are recorded in `FSCDXS_line_params.csv`, which was generated with a separate script not in version control. For now, if there is any overlap of the user-specified region and a HITRAN-recommended XS region, the XS parameters are used.

# Running the `absco-generate` Driver  <a name="driver"></a>

The generation pipeline is built from two modules -- `absco.preprocess` and `absco.compute`. Calculating cross sections for a given layer amounts to an optical depth (OD) calculation followed by dividing the OD by a layer amount (molecule number density) -- i.e. an LNFL and LBLRTM run, both driven by `absco.compute`. `absco-generate` selects which stage(s) to run via `-lnfl`, `-lbl`, or `-e2e` (end to end), and `-y` skips the interactive warning prompts. Before the line file generation and OD computation, `absco.preprocess` determines a number of things:

  - Configuration file read
  - Verify existence of necessary paths
  - Channels check
    - are channel inputs consistent (limits, input and output resolution have the same number of elements)?
    - are units defined and valid (wavenumber, micron, or nanometer)?
    - is the expected ratio of the output resolution to the LBLRTM resolution an integer exponent of 2?
    - are band limits in the correct order?
    - conversion to wavenumber fo use in LBLRTM
    - is the band no larger than 2000 cm<sup>-1</sup>? Break up into separate bins if it is
  - Verify that input molecules are allowed
  - Ensure that only two H<sub>2</sub>O values are provided
  - Confirm that each of the remaining configuration attributes has only one value assigned to it
  - Verify all necessary attributes have been found in the configuration file
  - Check that the input pressures are monotonic, then force them to be in descending order (surface to TOA)
  - Ensure that the same number of pressures are in the arrays from the input VMR (`vmrfile`) and P (`pfile`) files (but here is no check that the pressures are equal)
  - Find the molecules that are radiatively active in the given spectral regions and ask user if they want to include omitted or extra species (this provision allows the `molnames` field in the configuration file to be empty, but if it is, all active molecules will be processed)
  - Determine when to use cross section parameters as opposed to the line parameter database (as specified by `FSCDXS_line_params.csv`)
  - Calculate the kernel and weights for spectral degradation
  - Ascertain the source (AER LPD v3.9 or HITRAN 2012) for each molecule and each band
  - Compute the amount of memory needed for the calculation and output file

These items are part of a `configure` object that is required input for the `makeABSCO` class.

Note that the objects contain the information for all user specified species. In the driver, we loop over each molecule so that a single netCDF output file is made for each molecule. This allows for parallelization of the processing (e.g., using a single core per molecule), but the code to facilitate it has not yet been written (partially because of concerns of the amount of HD space and RAM is needed for a single molecule).

Additionally, this driver handles H<sub>2</sub>O-affected molecules (CO<sub>2</sub>, N<sub>2</sub>, O<sub>2</sub>, H<sub>2</sub>O) and O<sub>2</sub> by looping over VMRs (in ppmv) for both species, and then generating ABSCO objects at these VMRs. The water vapor VMRs are provided in the configuration file, while the oxygen VMRs are hard set at 1.9 * 10<sup>5</sup> and 2.3 * 10<sup>5</sup>. For this processing, the H<sub>2</sub> and O<sub>2</sub> profiles are held constant at a given VMR while others (N2 and CO2) follow the profiles given by `vmrfile` (see [Table 2](#Table2)).

## Defaults

The user must provide a spectral range in the configuration file. If multiple bands are provided, each band must have assigned to it an associated starting value, ending value, LBLRTM resolution, and output (netCDF) resolution.

Species specification is optional, but nothing is done if not molecule is provided. If `molnames` is not assigned anything, the code will check to see what molecules are radiatively active in the given band, then provide them to the user in a standard output prompt.

### Pressure levels

In the `VMR` subdirectory, `standard_atm_profiles.py` should be run if the user ever wants to use a different profile (rather than the default 1976 United States standard atmopshere provided in the repository). This module utilizes a standard atmosphere (the different kinds of standard atmospheres computed by LBLRTM are documented in the constructor of the vmrProfiles class) and the pressures expected by the user and performs an interpolation of the associated volume mixing ratios onto the user-specified grid. The interpolated profile is then used as a user-input atmosphere in the TAPE5 files that are generated for LBLRTM. Whatever the user chooses to be the output file name of `standard_atm_profiles.py` should be entered into the `vmrfile` field in `ABSCO_config.ini` (see [Table 2](#Table2)).

As documented in [Table 2](#Table2), pressures are extracted from `pfile`, and each pressure has an associated range of acceptable temperature values listed in `ptfile`. Both of these files are underneath the `PT_Grid`. By default, the pressures that are used are based on the AIRS instrument (`PT_Grid/AIRS_P_air.txt`). `build_temp_array.txt` is the associated `ptfile`.

### Allowed Molecules

While LBLRTM handles a number of molecules via line parameters and cross sections, the allowed molecules for ABSCO processing is more limited. From `absco.preprocess`, where we verify that the user-provided molecule can be processed:

```
allowed = ['H2O', 'CO2', 'O3', 'N2O', 'CO', 'CH4', 'O2', \
  'NO', 'SO2', 'NO2', 'NH3', 'HNO3', 'OCS', 'H2CO', 'N2', \
  'HCN', 'C2H2', 'HCOOH', 'C2H4', 'CH3OH', 'CCL4', 'CF4', \
  'F11', 'F12', 'F22', 'ISOP', 'PAN', 'HDO', 'BRO', 'O2-O2']
```

The molecule names match those in the HITRAN database and are case-insensitive.

## Binary Line Files (TAPE3)

Line files need to be generated for every molecule and spectral range. Depending on the range and the number of lines for a given species in the range, this can be a time consuming process. However, the TAPE3 files likely only need to be generated once and can be saved to disk, which can be done by setting the `-lnfl` keyword:

```
absco-generate -lnfl
```

In the call, we assume `ABSCO_config.ini` to be the configuration file, which contains the molecule name, spectral range, and output TAPE3 subdirectory, all of which are incorporated into the file naming convention of the TAPE3 files: `working_dir/TAPE3_dir/molecule/TAPE3_molname_wn1-wn2`. In the examples that follow, We will assume `ABSCO_config.ini` is populated with its default values (i.e., is unchanged from its original checkout).

LNFL runs are performed inside an `LNFL_Runs` directory (also defined in `ABSCO_config.ini`). Links to the LNFL executable and necessary input files (line coupling TAPE2, broadening parameters, full ASCII TAPE1 line file), and TAPE5 files that direct LNFL on what to do are also automatically generated and saved in the `LNFL_Runs/TAPE5_dir` subdirectory by default.

## Optical Depth and Absorption Coefficient Calculation

For the absorption coefficient calculation, LBLRTM must be run to compute a) optical depths, and b) layer amounts (done with the LBLATM subroutine). Once TAPE3 files are generated for the specified molecule and bands, LBLRTM TAPE5s (LBLRTM input file with state specifications and radiative transfer parameters) are written for every specified pressure level, temperature, and band combination. In each iteration, the optical depth spectrum is converted to absorption coefficients (*k* values) by dividing them by the layer amount for the given molecule. This *k* spectrum is then degraded to lessen the amount of RAM and hard drive space needed for the output.

The LBLRTM runs are followed by array manipulation (i.e., tranposing arrays to the dimensions that were agreed upon, adding fill values, etc.) and writing the necessary arrays to an output netCDF for the given species. LBLRTM is run inside the directory specified by `lbl_run_dir` in [Table 2](#Table2).

To run this LBLRTM process, use the driver script with the `-lbl` keyword.

```
absco-generate -lbl
```

The process is repeated over both water vapor VMRs for species whose continua are affected by water vapor (CO<sub>2</sub>, O<sub>2</sub>, and N<sub>2</sub>). Separate objects for each VMR are instantiated, then their ABSCO arrays are combined.

## End-to-end Run

Initially, it may be best to just run LNFL and LBLRTM in series, i.e., the end-to-end run. This can be done in the driver with the `e2e` keyord:

```
absco-generate -e2e
```

Separating the LNFL and LBL runs may be useful after the user has generated all of the TAPE3 files that they need, but it will not be detrimental to run the end-to-end mode everytime. The only bottlenecks are the LNFL runs and the loop over all LBL cases. The latter will happen whenever ABSCOs are computed, and the former will not take a noticeable amount of time because LNFL will not be run if the expected TAPE3 exists.

## Logging

`absco-generate` automatically records its output to a log directory (default:
`<intdir>/logs`, e.g. `logs/` under the working directory). Three files are written:

- `absco-generate.log` -- the full driver output (everything printed to the console,
  which is still shown on the terminal as well);
- `lnfl.log` -- the LNFL subprocess output, one section per TAPE5 run;
- `lblrtm.log` -- the LBLRTM subprocess output, one section per P/T/VMR run.

Use `--log_dir DIR` (or `-log DIR`) to write the logs elsewhere. Log files are opened
in append mode, so re-running accumulates history rather than overwriting it.

```
absco-generate -e2e --log_dir /path/to/logs
```

# Parallel Processing (split and join) <a name="parallel"></a>

A single `absco-generate` run processes all molecules and the full spectral range in
series, which can be slow and memory-hungry. Two helper commands let you break a run
into independent pieces, run them concurrently (e.g. one per core, or across a cluster),
and stitch the resulting tables back together.

## `absco-split-config`

Splits one `ABSCO_config.ini` into several smaller configuration files -- by molecule
and, optionally, by spectral chunk -- so each can be generated independently.

```
absco-split-config CONFIG_FILE [-m MOLECULES_PER_FILE] [-c CHUNK_SIZE]
```

| Option | Default | Description |
| :--- | :---: | :--- |
| `config_file` (positional) | -- | Source configuration file to split |
| `-m`, `--molecules_per_file` | `1` | Number of molecules per output config. `0` puts all molecules in one file |
| `-c`, `--chunk_size` | none | Width of each spectral sub-window, in the config's `units`. Omit (or `0`) to keep each band whole |

For every (spectral chunk × molecule group) combination it writes a new config named
`<base>-<NNN>-<wn1>_<wn2>-<molecules><ext>` (e.g. `ABSCO_config-001-4166_4358-ch4.ini`).
Each output config also gets a numbered subdirectory appended to its `intdir`
(`01`, `02`, ...) so concurrent runs write their intermediate/output files to separate
locations and do not clobber one another.

```
# one config per molecule, splitting the band into 100 cm-1 chunks
absco-split-config ABSCO_config.ini -m 1 -c 100
# then run each generated config (in parallel, on separate cores/nodes)
for cfg in ABSCO_config-*.ini; do absco-generate -e2e -y -i "$cfg"; done
```

## `absco-join-tables`

Concatenates ABSCO tables for the **same molecule** along their spectral axis. All inputs
must have been generated with the same configuration parameters (pressures, temperatures,
resolution, VMRs); only the spectral range differs. Input order does not matter -- the
tables are sorted by spectral extent before joining.

```
absco-join-tables FILENAME [FILENAME ...] [-o OUT_FILENAME]
```

| Option | Default | Description |
| :--- | :---: | :--- |
| `FILENAME ...` (positional) | -- | Two or more input `.nc` tables to combine |
| `-o`, `--out_filename` | derived | Output file. If omitted, a name is derived from the molecule and the combined extent, e.g. `H2O_4166-4358_joined.nc` |

```
absco-join-tables nc_ABSCO/01/H2O_*.nc nc_ABSCO/02/H2O_*.nc -o H2O_full.nc
```

The `scripts/run_multiple_configs.sh` and `scripts/join_multiple.sh` helper scripts drive
this workflow: the former runs a set of split config files concurrently across the
available cores (with `absco-generate`), and the latter joins the per-chunk outputs back
into one table per molecule (with `absco-join-tables`). Both expect the environment to be
active (e.g. `pixi shell`) so the `absco-*` commands are on `PATH`.

# Output netCDF <a name="output"></a>

Output from this software is relatively simple but potentially very large (when the software is run, a prompt will warn the user of the potential size of the file and memory footprint of the code). This is primarily because of the many dimensions of the cross section arrays (`[nP x nT x nSpec]`), with the `nSpec` dimension being the biggest contributor. The number of spectral points can be manipulated by the input bandwidth, resolution, and spectral degradation. The output are stored in a netCDF in the `outdir` directory listed in [Table 2](#Table2) (`nc_ABSCO` by default). A single netCDF is generated for each molecule.

Cross section array dimensions are also dependent on the molecule, because species whose continuum is dependent on water vapor amount contain an extra dimension for the water vapor VMR. Absorption coefficients are calculated at two different H<sub>2</sub>O VMRs at all applicable pressures and temperatures, so another dimension is necessary to store all of the calculations. All possible dimensions are given in [Table 3](#Table3).

**netCDF Dimensions** <a id="Table3"></a>

| Dimension Name | Description |
| :---: | :---: |
| nfreq | Number of spectral points (wavenumbers, microns, nanometers) |
| nlev | Number of pressure levels |
| nlay | Number of pressure layers (nlev-1) |
| ntemp | Number of temperatures. This is always 15, regardless of the input pressures; and any invalid temperature (for a given pressure) will just be populated with fill values for *k* |
| nranges | Number of spectral ranges (AKA bands, channels)
| nranges_lims | Number of limits for spectral ranges. This is always two because  there is a single minimum and single maximum spectral point associated with each range |
| nvmr | Number of H<sub>2</sub>O VMR values. This will always be 2 |

Since their continua are dependent on water vapor content, netCDF variables for H<sub>2</sub>O, CO<sub>2</sub>, O<sub>2</sub>, and N<sub>2</sub> will include all of the dimensions listed in [Table 3](#Table3). [Table 4](#Table4) contains the names, dimensions, units, valid ranges, and descriptions of variables for these four molecules.

**Output Arrays for Water Vapor-affected Species** <a id="Table4"></a>

| Array Name | Dimensions | Units | Range | Description |
| :---: | :---: | :---: | :---: | :--- |
| P_level | (nlev) | mbar | [0, 1050] | User-provided layer <br> boundary pressures |
| P_layer | (nlay x ntemp) | mbar | [0, 1050] | LBLRTM-calculated <br> layer pressures |
| Cross_Section | (nfreq x ntemp x nlay x nvmr) | cm<sup>-2</sup> | [0, 10<sup>-20</sup>] | Absorption coefficients *k* <br> calculated from <br> LBLRTM optical depths <br> and layer amounts |
| Spectral_Grid | (nfreq) | cm<sup>-1</sup> | [0, 50000] | Spectral points corresponding to <br> ABSCOs in a single layer <br> for a single temperature <br> and in a given spectral range |
| Temperature | (nlev x ntemp) | K | [180, 320] | Applicable temperatures <br> associated with each <br> layer boundary pressure |
| Extent_Ranges | (nranges x nranges_lims) | cm<sup>-1</sup> | [0, 50000] | Starting and ending spectral points <br> for each input channel |
| Extent_Indices | (nranges x nranges_lims) | unitless | [0, "infinity"] | Pairs of indices defining the <br> start and end indices of the <br> Cross_Section frequency dimension <br> for non-continuous calculation regions |
| H2O_VMR | (nvmr) | ppmv | [0, 50000] | Water vapor volume mixing ratio |
| O2_VMR | (nvmr) | ppmv | [0, 250000] | Oxygen volume mixing ratio. <br>Only added to oxygen *k* array|

All other allowed molecules will conform to a similar convention as [Table 4](#Table4), only without the H<sub>2</sub>O VMR dimension. There is thus no `H2O_VMR` array for these species.

In the global attributes, there is a "source" field. There are only two sources -- HITRAN 2012 and AER LPD v3.9 -- and they can vary by band. The code accounts for the by-band differences in source. All fill values in the netCDF files are NaN (not-a-number).

## Reading the Output

The `absco-read` command reads a file like the ones in the [Output netCDF](#output) section and prints the cross section value at a given coordinate (along with the zero-offset array coordinates of the absorption coefficient):

```
% absco-read -h
usage: absco-read [-h] [-p IN_PRESSURE] [-T IN_TEMP]
                  [-s IN_SPECTRAL IN_SPECTRAL] [-wv IN_H2O] [-o2 IN_O2]
                  [-tol TOLERANCE]
                  ncFile

Read netCDF generated with absco-generate and print out absorption coefficient
(k) at a given pressure, temperature, spectral point, and water vapor amount
(if molecule continuum is affected by water vapor).

positional arguments:
  ncFile                Output netCDF generated by absco-generate.

options:
  -h, --help            show this help message and exit
  -p, --in_pressure IN_PRESSURE
                        Reference pressure level [mbar] for which k is
                        retrieved. There are two pressure boundaries for a
                        given layer, and this is the lower bound (ie, closest
                        to surface). (default: 1050.0)
  -T, --in_temp IN_TEMP
                        Reference temperature [K] for which k is retrieved.
                        (default: 230.0)
  -s, --in_spectral IN_SPECTRAL IN_SPECTRAL
                        Reference spectral point AND units [cm-1, um, or nm]
                        for which k is retrieved. (default: [500, 'cm-1'])
  -wv, --in_h2o IN_H2O  Reference water vapor VMR (ppmv) for which k is
                        retrieved IF the specified molecule is H2O, CO2, O2,
                        or N2. (default: 10.0)
  -o2, --in_o2 IN_O2    Reference oxygen VMR (ppmv) for which k is retrieved
                        IF the specified molecule is O2. (default: 190000.0)
  -tol, --tolerance TOLERANCE
                        Tolerance used when searching for floating point
                        matches in *each* of the dimensions. This should be a
                        relative tolerance (e.g. 0.01 would mean P from netCDF
                        is within 1% of in_pressure). (default: None)
```
