# ABSoprtion COefficient (ABSCO) Table Generation

NOTE: this README will be updated more thoroughly and accurately once the software is completed.

Software that can generate a set of ABSCO tables that span the thermal IR to UV spectral range.

Assuming the user has set the `user.name` and `user.email` Git configuration variables and has setup their Github account SSH keys, the repository can be cloned by using:

```
git clone --recursive git@github.com:pernak18/ABSCO.git
```

Note the `--recursive` keyword -- it will force the clone to copy the necessary subroutines. Otherwise, the user will have to do (assuming the repository is just checked out with the default `ABSCO` name in the working directory):

```
cd ABSCO
git submodule init
git submodule update
```

# Dependencies

This library depends on a number of standard Python libraries (`os`, `sys`, `configparser`, `subprocess`, `glob`, `argparse`), some widely distributed third-party libraries (see "Python Packages" subsection), and some ad hoc subroutines that are available as a GitHub repository in the pernak18 account (`utils`, `RC_utils`, `lblTools`). The latter are located in the `common` subdirectory, which is its own repository (<https://github.com/pernak18/common>) but also a submodule that is cloned along with the ABSCO repository if submodules are updated (e.g., `--recursive` clone). Additionally, some AER-maintained models and line parameters are required to run the ABSCO software.

## Python Packages

The following libraries were installed with miniconda (<https://conda.io/docs/user-guide/install/index.html>) for Python 3.6.3:

  - numpy (v1.13.3)
  - scipy (v1.0.0)
  - pandas (v0.23.0)
  - netCDF4 (v1.3.1)

All are used in this ABSCO library. The software is optimized for Python 3 usage -- Python 2 has not been tested and is not recommended.

## LNFL, LBLRTM, and the AER Line File

LNFL (LiNe FiLe) FORTRAN code that converts ASCII text line parameter files to the binary files that LBLRTM expects (TAPE3) is located in its own repository (<https://github.com/pernak18/LNFL>). Because LNFL has been declared a submodule of the ABSCO library, using the `--recursive` keyword in the clone of this ABSCO repository will also clone the LNFL source code the is necessary. The source code is fetched and saved under the `LNFL` subdirectory.

LBLRTM (Line-By-Line Radiative Transfer Model) FORTRAN code also has its own Git repository (<https://github.com/pernak18/LBLRTM>) and is a declared ABSCO submodule. It is stored under the `LBLRTM` subdirectory. LBLRTM in the context of this software simply calculates optical depth at specified pressures, temperatures, and spectral ranges.

The AER line parameter database (LPD) is distributed as a set of ASCII text files in the `AER_Line_File` directory.

Periodically, the models and LPD will be updated to reflect new line parameters, a new continuum, or bug fixes. These revisions can have significant effects on the model output. For reference, the model and LPD version numbers associated with the initial release of the ABSCO software are:

| Model | Version |
| :---: | :---: |
| LNFL | v3.1 |
| LBLRTM | v12.9 |
| LPD | v3.6 |
---

Currently, there are no plans on updating these three repositories. In the future, we may set up a separate AER account that will contain model code for the public rather than hosting in my personal account.

LNFL and LBLRTM can be built with the `build_models.py` script:

```
% ./build_models.py -c gfortran -i ABSCO_config.ini
```

This script call specifies a `gfortran` compiler (`-c`) and replaces the paths to the executables in `ABSCO_config.ini` with the paths of the newly-built executables. Other valid compilers are `ifort` and `pgf90`. Use the `-h` option with the script for more options. Path replacement also occurs with the line file-associated paths (`tape1_path`, `tape2_path`, `extra_params`, `xs_path`, and `fscdxs`). If `-i` is not set, no path replacement occurs even though the executable are compiled. The two executables follow the naming convention `model_version_os_compiler_precision`, e.g., `lblrtm_v12.9_linux_intel_dbl`.

# Setup (Configuration File)

With the exception of the `--run_lnfl`, `--run_lbl`, and `--end_to_end` (alternatively `-lnfl`, `-lbl`, or `-e2e`) keywords that dictate which executable will be utilized, `ABSCO_config.ini` contains all of the inputs expected from the user. All of the following parameters are expected in the file:

| Field | Notes |
| :---: | :--- |
| header | 80-character header that is written to LBL optical depth files but is otherwise not used|
| pfile | text file with 1 pressure layer \[mbar\] per line. these will be the pressures on which the ABSCOs are calculated. this is a *relative* path with respect to the working directory |
---

`ABSCO_config.ini` can be named something else, but that will have to be specified at the command line (otherwise it's assumed that `ABSCO_config.ini` is the configuration file to use):

```
ABSCO_tables.py -i your_ini_file
```

## FSCDXS

Some molecules have both line parameters and XS parameters.  HITRAN makes recommendations on the preferred parameters given the species and the band, and these are taken into account in the error checking that the `ABSCO_preprocess.py` module does.  Molecules where line parameters are recommended, the associated bands, and the flag (0: only XS exist, 1: both exist, use line params, 2: both exist, use XS) are recorded in `FSCDXS_line_params.csv`, which was generated with a separate, non-version controlled script.

# Binary Line Files (TAPE3 files)

Line file need to be generated for every molecule and spectral range. Depending on the range and the number of lines for a given species in the range, this can be a time consuming process. However, the TAPE3 files likely only need to be generated once and can be saved to disk, which can be done by setting the `-lnfl` keyword:

```
ABSCO_tables.py -lnfl
```

In the call, we assume `ABSCO_config.ini` to be the configuration file, which contains the molecule name, spectral range, and output TAPE3 subdirectory, all of which are incorporated into the file naming convention of the TAPE3 files: `working_dir/TAPE3_dir/molecule/TAPE3_wn1-wn2`.

LNFL runs are performed inside an `LNFL_Runs` directory (also in `ABSCO_config.ini`). Links to the LNFL executable and necessary input files (line coupling TAPE2, broadening parameters, full ASCII TAPE1 line file), and TAPE5 files that direct LNFL on what to do are also automatically generated and saved in the `LNFL_Runs/TAPE5_dir` subdirectory by default.

# Optical Depth Files (ODint_001)

```
ABSCO_tables.py -lbl
```

# End-to-end Run

```
ABSCO_tables.py -e2e
```

# Defaults

User must provide a spectral range. Species specification is optional -- if `molnames` is not assigned anything, the code will check to see what molecules are radiatively active in the given band.

# User Options

can be set in .ini file that follows ABSCO_config.ini convention.

## Pressure levels

In the `VMR` subdirectory, `standard_atm_profiles.py` should be run if the user ever wants to use a different profile (rather than the default 1976 United States standard atmopshere provided in the repository). This module utilizes a standard atmosphere (the different kinds of standard atmospheres computed by LBLRTM are documented in the constructor of the vmrProfiles class) and the pressures expected by the user and performs an interpolation of the associated volume mixing ratios onto the user-specified grid. The interpolated profile is then used as a user-input atmosphere in the TAPE5 files that are generated for LBLRTM. Whatever the user chooses to be the output file name of `standard_atm_profiles.py` should be entered into the `vmrfile` field in `ABSCO_config.ini`.

For whatever atmosphere is used, one must also calculate the broadening density at each pressure layer for each molecule (i.e., when all of the other molecules are zeroed out in the TAPE5). This can also be done with the `standard_atm_profiles.py` module using the `--broad` keyword. This should be done whenever a new profile is calculated. The corresponding output file should be entered into the `broadfile` field in `ABSCO_config.ini` as well.

## Allowed Molecules

While LBLRTM handles a number of molecules via line parameters and cross sections, the allowed molecules for ABSCO processing is more limited.

```
allowed = ['H2O', 'CO2', 'O3', 'N2O', 'CO', 'CH4', 'O2', \
  'NO', 'SO2', 'NO2', 'NH3', 'HNO3', 'OCS', 'H2CO', 'N2', \
  'HCN', 'C2H2', 'HCOOH', 'C2H4', 'CH3OH', 'CCL4', 'CF4', \
  'F11', 'F12', 'F22', 'ISOP', 'PAN', 'HDO', 'BRO', 'O2-O2']
```

water vapor is treated as a number of particles, depending on the number of PWV values provided. in `ABSCO_tables.py`, LBLRTM TAPE5 files are generated in a loop over PWV, and separate TAPE5s are generated for the self and foreign continua.

all HITRAN molecule names (these are the molecules for which we have line parameters)

```
# might be useful later...list the HITRAN molecule names
lfMolDir = '/nas/project/rc_static/models/' + \
  'aer_line_parameters/AER_line_files/aer_v_3.6/' + \
  'line_files_By_Molecule/*'
molDirs = sorted(glob.glob(lfMolDir))

# the upper() takes care of the Br problem
htMols = [os.path.basename(md).split('_')[1].upper() for \
  md in molDirs]
print(htMols)
```

# Reading the Output

```
% read_ABSCO_tables.py -h
usage: read_ABSCO_tables.py [-h] [-p IN_PRESSURE] [-T IN_TEMP]
                            [-s IN_SPECTRAL IN_SPECTRAL] [-wv IN_H2O]
                            [-tol TOLERANCE]
                            ncFile

Read in netCDF generated with ABSCO_tables.py module and print out absorption
coefficient (k) at a given pressure, temperature, spectral point, and water
vapor amount (if molecule continuum is affected by water vapor).

positional arguments:
  ncFile                Output netCDF generated by ABSCO_tables.py.

optional arguments:
  -h, --help            show this help message and exit
  -p IN_PRESSURE, --in_pressure IN_PRESSURE
                        Reference pressure level [mbar] for which k is
                        retrieved. There are two pressure boundaries for a
                        given layer, and this is the lower bound (ie, closest
                        to surface). (default: 1050.0)
  -T IN_TEMP, --in_temp IN_TEMP
                        Reference temperature [K] for which k is retrieved.
                        (default: 230.0)
  -s IN_SPECTRAL IN_SPECTRAL, --in_spectral IN_SPECTRAL IN_SPECTRAL
                        Reference spectral point AND units [cm-1, um, or nm]
                        for which k is retrieved. (default: [500, 'cm-1'])
  -wv IN_H2O, --in_h2o IN_H2O
                        Reference water vapor VMR (ppmv) for which k is
                        retrieved IF the specified molecule is H2O, CO2, O2,
                        or N2. (default: 10000.0)
  -tol TOLERANCE, --tolerance TOLERANCE
                        Tolerance used when searching for floating point
                        matches in *each* of the dimensions. This should be a
                        relative tolerance (e.g. 0.01 would mean P from netCDF
                        is within 1% of in_pressure). (default: None)
```

