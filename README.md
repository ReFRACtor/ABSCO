# ABSoprtion COefficient (ABSCO) Table Generation

NOTE: this README will be updated more thoroughly and accurately once the software is completed.

Software that can generate a set of ABSCO tables that span the thermal IR to UV spectral range.

Assuming the user has set the `user.name` and `user.email` Git configuration variables and has setup their Github account SSH keys, the repository can be cloned by using:

```
git clone --recursive git@github.com:pernak18/ABSCO.git
```

Note the `--recursive` keyword -- it will force the clone to copy the necessary subroutines.

# Dependencies

This library depends on a number of standard Python libraries (`os`, `sys`, `configparser`, `subprocess`, `glob`, `argparse`), some widely distributed third-party libraries (see "Python Packages" subsection), and some ad hoc subroutines that are available as a GitHub repository in the pernak18 account (`utils`, `RC_utils`, `lblTools`). The latter are located in the `common` subdirectory.

In addition to Python dependencies, source code for a couple of AER models (LNFL and LBLRTM) as well as the line parameter database are required. LNFL and LBLRTM code is available underneath the `LNFL` and `LBLRTM` subdirectories, respectively, and the AER line parameter database is distributed as a set of ASCII text files in the `AER\_Line\_File` directory. Both models can be built with the `build_models.py` script:

```
% ./build_models.py -c gfortran -i ABSCO_config.ini
```

This script call specifies a `gfortran` compiler (`-c`) and replaces the paths to the executables in ABSCO\_config.ini with the paths of the newly-built executables. Other valid compilers are `ifort` and `pgf90`. Use the `-h` option with the script for more options.

## Python Packages

The following libraries were installed with miniconda for Python 3 (the software is also optimized for Python 3 usage -- Python 2 has not been tested):

  - numpy
  - scipy
  - pandas
  - netCDF4

All are used in this ABSCO library.

## LNFL

## LBLRTM

# Setup (Configuration File)

`ABSCO_config.ini` contains all of the inputs expected from the user.

## FSCDXS

Some molecules have both line parameters and XS parameters.  HITRAN makes recommendations on the preferred parameters given the species and the band, and these are taken into account in the error checking that the ABSCO\_preprocess.py module does.  Molecules where line parameters are recommended, the associated bands, and the flag (0: only XS exist, 1: both exist, use line params, 2: both exist, use XS) are recorded in `FSCDXS_line_params.csv`, which was generated with a separate, non-version controlled script.  Currently, there is only one molecule (in one band) where the line parameters are recommended over the XS coefficients (CH<sub>3</sub>CN), and it is not a molecule that ABSCO_tables.py processes.

# Binary Line Files

Running LNFL

# Optical Depth Files

Running LBLRTM

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

