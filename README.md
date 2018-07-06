# ABSoprtion COefficient (ABSCO) Table Generation

Software that can generate a set of ABSCO tables that span the thermal IR to UV spectral range.

Assuming the user has set the `user.name` and `user.email` Git configuration variables and has setup their Github account SSH keys, the repository can be cloned by using:

```
git clone --recursive git@github.com:pernak18/ABSCO.git
```

Note the `--recursive` keyword -- it will force the clone to copy the necessary subroutines.

# Dependencies

In `common` subroutine directory. Some dependencies are Python modules, others are FORTRAN source code that will need to be built by the user (note to Rick: do a quick Makefile).

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

# Binary Line Files

Running LNFL

# Optical Depth Files

Running LBLRTM

# Defaults

User must provide a spectral range. Species specification is optional -- if `molnames` is not assigned anything, the code will check to see what molecules are radiatively active in the given band.

# User Options

## Pressure levels

In `VMR` subdirectory, run standard_atm_profiles.py.

