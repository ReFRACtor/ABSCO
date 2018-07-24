# ABSoprtion COefficient (ABSCO) Table Generation

Software that can generate a set of ABSCO tables that span the thermal IR to UV spectral range.

Assuming the user has set the `user.name` and `user.email` Git configuration variables and has setup their Github account SSH keys, the repository can be cloned by using:

```
git clone --recursive git@github.com:pernak18/ABSCO.git
```

Note the `--recursive` keyword -- it will force the clone to copy the necessary subroutines.

# Dependencies

In `common` subroutine directory. Some dependencies are Python modules, others are FORTRAN source code that will need to be built by the user (note to Rick: do a quick Makefile). NOTE: only Linux executables built with the PGI compiler have been tested with the software in this repository.

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

## Pressure levels

In the `VMR` subdirectory, `standard_atm_profiles.py` should be run if the user ever wants to use a different profile (rather than the default 1976 United States standard atmopshere provided in the repository). This module utilizes a standard atmosphere (the different kinds of standard atmospheres computed by LBLRTM are documented in the constructor of the vmrProfiles class) and the pressures expected by the user and performs an interpolation of the associated volume mixing ratios onto the user-specified grid. The interpolated profile is then used as a user-input atmosphere in the TAPE5 files that are generated for LBLRTM. Whatever the user chooses to be the output file name of `standard_atm_profiles.py` should be entered into the `vmrfile` field in `ABSCO_config.ini`.

For whatever atmosphere is used, one must also calculate the broadening density at each pressure layer for each molecule (i.e., when all of the other molecules are zeroed out in the TAPE5). This can also be done with the `standard_atm_profiles.py` module using the `--broad` keyword. This should be done whenever a new profile is calculated. The corresponding output file should be entered into the `broadfile` field in `ABSCO_config.ini` as well.

