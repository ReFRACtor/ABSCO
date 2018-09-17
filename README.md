# ABSoprtion COefficient (ABSCO) Table Generation

---
**Contents**

1. [Introduction](#intro)
2. [Dependencies](#dependencies)
3. [Setup](#setup)
4. [HITRAN Cross Sections](#xs)
5. [Driver Script](#driver)
6. [Output](#output)

# Introduction <a name="intro"></a>

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

# Dependencies <a name="dependencies"></a>

This library depends on a number of standard Python libraries (`os`, `sys`, `configparser`, `subprocess`, `glob`, `argparse`), some widely distributed third-party libraries (see "Python Packages" subsection), and some ad hoc subroutines that are available as a GitHub repository in the pernak18 account (`utils`, `RC_utils`, `lblTools`). The latter are located in the `common` subdirectory, which is [its own repository] (<https://github.com/pernak18/common>) but also a submodule that is cloned along with the ABSCO repository if submodules are updated (e.g., `--recursive` clone). Additionally, some AER-maintained models and line parameters are required to run the ABSCO software.

## Python Packages

The following libraries were [installed with miniconda] (<https://conda.io/docs/user-guide/install/index.html>) for Python 3.6.3:

  - numpy (v1.13.3)
  - scipy (v1.0.0)
  - pandas (v0.23.0)
  - netCDF4 (v1.3.1)
  - xarray (v0.10.8)

All are used in this ABSCO library. **The software is optimized for Python 3 usage -- Python 2 has not been tested and is not recommended.**

## LNFL, LBLRTM, and the AER Line File

LNFL (LiNe FiLe) FORTRAN code that converts ASCII text line parameter files to the binary files that LBLRTM expects (TAPE3) is located in [its own repository] (<https://github.com/pernak18/LNFL>). Because LNFL has been declared a submodule of the ABSCO library, using the `--recursive` keyword in the clone of this ABSCO repository will also clone the LNFL source code the is necessary. The source code is fetched and saved under the `LNFL` subdirectory.

LBLRTM (Line-By-Line Radiative Transfer Model) FORTRAN code also has [its own Git repository] (<https://github.com/pernak18/LBLRTM>) and is a declared ABSCO submodule. It is stored under the `LBLRTM` subdirectory. LBLRTM in the context of this software simply calculates optical depth at specified pressures, temperatures, and spectral ranges.

The AER line parameter database (LPD) is distributed as a set of ASCII text files in the `AER_Line_File` directory. Currently, it is available on the AER external Git server and is linked as a submodule of this ABSCO repository (so a `--recursive` clone will take care of this dependency as well). End users will need to contact Rick Pernak or Karen Cady-Pereira so that they can be granted access to the LPD repo.

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

# Setup (Configuration File) <a name="setup"></a>

With the exception of the `--run_lnfl`, `--run_lbl`, and `--end_to_end` (alternatively `-lnfl`, `-lbl`, or `-e2e`) keywords that dictate which executable will be run, `ABSCO_config.ini` contains all of the inputs expected from the user. All of the following parameters are expected in the file:

| Field | Notes |
| :---: | :--- |
| header | 80-character header that is written to LBL optical depth files but is otherwise not used|
| pfile | text file with 1 pressure level in millibars per line. these will be the pressures on which the ABSCOs are calculated. this can be a *relative* path with respect to the working directory or an absolute path|
| ptfile | for every pressure level, there are different allowed temperatures. this file contains a set of pressures and their permitted temperatures|
| vmrfile | CSV files generated with VMR/standard_atm_profiles.py that provide interpolated/extrapolated volume mixing ratios (VMRs) for entire user-specified profile (can be relative to working directory) |
| xs_lines | this file contains the species names for when XS and line parameters exist and line parameter usage is recommended by HITRAN |
| wn1, wn2 | starting and ending spectral points for every desired band. can be in wavenumbers, microns, or nanometers |
| lblres, outres | spectral resolution at which LBLRTM is run and spectral resolution of the output (after spectral degradation). for now, this should be in wavenumbers |
| units | spectral units ("cm<sup>-1</sup>", "um", and "nm" are allowed) |
| wv_vmr | water vapor VMR (in **ppmv**) that will be used for all levels in a given profile for H<sub>2</sub>O, CO<sub>2</sub>, O<sub>2</sub>, and N<sub>2</sub>, since their continua are dependent on water vapor amounts |
| molnames | HITRAN molecule names, space delimited, case-insensitive. really should only be one molecule per run |
| scale | multiplicative factors used for continuum or extinction scaling (separate factors for H<sub>2</sub>O self continuum, H<sub>2</sub>O foreign continuum, CO<sub>2</sub> continuum, O<sub>3</sub> continuum, O<sub>2</sub> continuum, N<sub>2</sub> continuum, and Rayleigh extinction) |
| tape5_dir | Directory underneath both lnfl_run_dir and lbl_run_dir to which TAPE5 files will be written (should just be a single string) |
| lnfl_run_dir | Path to directory where LNFL runs will occur. Additional subdirectories (one for each molecule) will be created underneath this directory. can be full or relative path. assignment can be automated with build_models.py |
| tape1_path | Full path to the full TAPE1 ASCII line file that should be used in LNFL runs. assignment can be automated with build_models.py |
| tape2_path | Full path to the full TAPE2 ASCII line coupling file that should be used in LNFL runs with O<sub>2</sub>, CO<sub>2</sub>, and CH<sub>4</sub>. assignment can be automated with build_models.py |
| lnfl_path | Full path to LNFL executable to be run. assignment can be automated with build_models.py |
| extra_params | directory with CO<sub>2</sub>-CO<sub>2</sub>, CO<sub>2</sub>-H<sub>2</sub>O, O<sub>2</sub>-H<sub>2</sub>O, O<sub>2</sub>-O<sub>2</sub>, and H<sub>2</sub>O-CO<sub>2</sub> broadening parameter specifications. assignment can be automated with build_models.py |
| tape3_dir | directory relative to the working directory to which LNFL output (binary line files)  will be written|
| lbl_path | Full path to LBLRTM executable to be run. assignment can be automated with build_models.py |
| xs_path | Full path to LBLRTM cross section file directories. assignment can be automated with build_models.py |
| fscdxs | Full path to cross section "lookup" file used with xs_path. assignment can be automated with build_models.py |
| lbl_run_dir | Path to directory where LBL runs will occur. Additional subdirectories (one for each molecule) will be created underneath this directory. Can be a relative path with respect to working directory |
| outdir | directory where the netCDFs will be written. can be relative to working directory |
| sw_ver | used in the output netCDF file to specify the software version number |
| out_file_desc | used in the output netCDF file. allows use to document run more specifically |
| nc_compress | compression level for netCDF output |
---

`ABSCO_config.ini` can be named something else, but that will have to be specified at the command line (otherwise it's assumed that `ABSCO_config.ini` is the configuration file to use):

```
run_LBLRTM_ABSCO.py -i your_ini_file
```

# HITRAN Cross Section Usage <a name="xs"></a>

Some molecules have both line parameters and XS parameters.  HITRAN makes recommendations on the preferred parameters given the species and the band, and these are taken into account in the error checking that the `ABSCO_preprocess.py` module does.  Molecules where line parameters are recommended, the associated bands, and the flag (0: only XS exist, 1: both exist, use line params, 2: both exist, use XS) are recorded in `FSCDXS_line_params.csv`, which was generated with a separate, script not in version control. For now, if there is any overlap of the user-specified region and a HITRAN-recommended XS region, the XS parameters are used.

# Running the `run_LBLRTM_ABSCO.py` Driver Script  <a name="driver"></a>

## Defaults

The user must provide a spectral range in the configuration file. If multiple bands are provided, each band must have assigned to it an associated starting value, ending value, LBLRTM resolution, and output (netCDF) resolution.

Species specification is optional, but nothing is done if not molecule is provided. If `molnames` is not assigned anything, the code will check to see what molecules are radiatively active in the given band, then provide them to the user in a standard output prompt.

### Pressure levels

In the `VMR` subdirectory, `standard_atm_profiles.py` should be run if the user ever wants to use a different profile (rather than the default 1976 United States standard atmopshere provided in the repository). This module utilizes a standard atmosphere (the different kinds of standard atmospheres computed by LBLRTM are documented in the constructor of the vmrProfiles class) and the pressures expected by the user and performs an interpolation of the associated volume mixing ratios onto the user-specified grid. The interpolated profile is then used as a user-input atmosphere in the TAPE5 files that are generated for LBLRTM. Whatever the user chooses to be the output file name of `standard_atm_profiles.py` should be entered into the `vmrfile` field in `ABSCO_config.ini`.

### Allowed Molecules

While LBLRTM handles a number of molecules via line parameters and cross sections, the allowed molecules for ABSCO processing is more limited. From `ABSCO_preprocess.py`, where we verify that the user-provided molecule can be processed:

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
run_LBLRTM_ABSCO.py -lnfl
```

In the call, we assume `ABSCO_config.ini` to be the configuration file, which contains the molecule name, spectral range, and output TAPE3 subdirectory, all of which are incorporated into the file naming convention of the TAPE3 files: `working_dir/TAPE3_dir/molecule/TAPE3_wn1-wn2`. In the examples that follow, We will assume `ABSCO_config.ini` is populated with its default values (i.e., is unchanged from its original checkout).

LNFL runs are performed inside an `LNFL_Runs` directory (also defined in `ABSCO_config.ini`). Links to the LNFL executable and necessary input files (line coupling TAPE2, broadening parameters, full ASCII TAPE1 line file), and TAPE5 files that direct LNFL on what to do are also automatically generated and saved in the `LNFL_Runs/TAPE5_dir` subdirectory by default.

## Optical Depth and Absorption Coefficient Calculation

For the absorption coefficient calculation, LBLRTM must be run to compute a) optical depths, and b) layer amounts (done with the LBLATM subroutine). Once TAPE3 files are generated for the specified molecule and bands, LBLRTM TAPE5s (LBLRTM input file with state specifications and radiative transfer parameters) are written for every specified pressure level, temperature, and band combination. In each iteration, the optical depth spectrum is converted to absorption coefficients (*k* values) by dividing them by the layer amount for the given molecule. This *k* spectrum is then degraded to lessen the amount of RAM and hard drive space needed for the output.

The LBLRTM runs are followed by array manipulation (i.e., tranposing arrays to the dimensions that were agreed upon, adding fill values, etc.) and writing the necessary arrays to an output netCDF for the given species.

To run this LBLRTM process, use the driver script with the `-lbl` keyword.

```
run_LBLRTM_ABSCO.py -lbl
```

The process is repeated over both water vapor VMRs for species whose continua are affected by water vapor (CO<sub>2</sub>, O<sub>2</sub>, and N<sub>2</sub>). Separate objects for each VMR are instantiated, then their ABSCO arrays are combined.

## End-to-end Run

Initially, it may be best to just run LNFL and LBLRTM in series, i.e., the end-to-end run. This can be done in the driver with the `e2e` keyord:

```
run_LBLRTM_ABSCO.py -e2e
```

Separating the LNFL and LBL runs may be useful after the user has generated all of the TAPE3 files that they need, but it will not be detrimental to run the end-to-end mode everytime. The only bottlenecks are the LNFL runs and the loop over all LBL cases. The latter will happen whenever ABSCOs are computed, and the former will not take a noticeable amount of time because LNFL will not be run if the expected TAPE3 exists.

# Output netCDF <a name="output"></a>

Output from this software is relatively simple but potentially very large (when the software is run, a prompt will warn the user of the potential size of the file and memory footprint of the code). This is primarily because of the many dimensions of the cross section arrays (`[nP x nT x nSpec]`), with the `nSpec` dimension being the biggest contributor. The number of spectral points can be manipulated by the input bandwidth, resolution, and spectral degradation.

Cross section array dimensions are also dependent on the molecule, because species whose continuum is dependent on water vapor amount contain an extra dimension for the water vapor VMR. Absorption coefficients are calculated at two different H<sub>2</sub>O VMRs at all applicable pressures and temperatures, so another dimension is necessary to store all of the calculations. netCDF headers for H<sub>2</sub>O, CO<sub>2</sub>, O<sub>2</sub>, and N<sub>2</sub> will follow the convention of this example:

```
% ncdump -h nc_ABSCO/O2_00000-00150_v0.0_init.nc 
netcdf O2_00000-00150_v0.0_init {
dimensions:
	nfreq = 375001 ;
	nlev = 3 ;
	nlay = 2 ;
	ntemp = 15 ;
	nranges = 1 ;
	nranges_lims = 2 ;
	nvmr = 2 ;
variables:
	double P_level(nlev) ;
		P_level:_FillValue = NaN ;
		P_level:units = "mbar" ;
		P_level:long_name = "Pressure Levels" ;
		P_level:valid_range = 0., 1050. ;
		P_level:description = "User-provided layer boundary pressures" ;
	double P_layer(nlay, ntemp) ;
		P_layer:_FillValue = NaN ;
		P_layer:units = "mbar" ;
		P_layer:long_name = "Layer Pressures" ;
		P_layer:valid_range = 0., 1050. ;
		P_layer:description = "LBLRTM-calculated layer pressures" ;
	double Cross_Section(nfreq, ntemp, nlay, nvmr) ;
		Cross_Section:_FillValue = NaN ;
		Cross_Section:units = "cm2/molecule" ;
		Cross_Section:long_name = "Absorption Coefficients" ;
		Cross_Section:valid_range = 0., 1.e-20 ;
		Cross_Section:description = "Absorption coefficients (K) calculated from LBLRTM optical depths and layer amounts." ;
	double Spectral_Grid(nfreq) ;
		Spectral_Grid:_FillValue = NaN ;
		Spectral_Grid:units = "cm-1" ;
		Spectral_Grid:long_name = "Spectral Points" ;
		Spectral_Grid:valid_range = 0., 50000. ;
		Spectral_Grid:description = "Spectral points corresponding to ABSCOs in a single layer for a single temperature and in a given spectral range (wavenumbers, microns, or nanometers)." ;
	double Temperature(nlev, ntemp) ;
		Temperature:_FillValue = NaN ;
		Temperature:units = "K" ;
		Temperature:long_name = "Temperature Levels" ;
		Temperature:valid_range = 180., 320. ;
		Temperature:description = "Applicable temperatures associated with each layer boundary pressure" ;
	double Extent_Ranges(nranges, nranges_lims) ;
		Extent_Ranges:_FillValue = NaN ;
		Extent_Ranges:units = "cm-1" ;
		Extent_Ranges:long_name = "Spectral Ranges" ;
		Extent_Ranges:valid_range = 0., 50000. ;
		Extent_Ranges:description = "Starting and ending spectral points for each input channel." ;
	int64 Extent_Indices(nranges, nranges_lims) ;
		Extent_Indices:_FillValue = 9223372036854775807L ;
		Extent_Indices:units = "N/A" ;
		Extent_Indices:long_name = "Spectral Ranges Reference Indices" ;
		Extent_Indices:valid_range = 0L, 9223372036854775807L ;
		Extent_Indices:description = "Pairs of indices defining the start and end indices of the Cross_Secion frequency dimension for non-continuous calculation regions." ;
	double H2O_VMR(nvmr) ;
		H2O_VMR:_FillValue = NaN ;
		H2O_VMR:units = "ppmv" ;
		H2O_VMR:long_name = "Water Vapor Mixing Ratio" ;
		H2O_VMR:valid_range = 0., 50000. ;
		H2O_VMR:description = "Water vapor amount that influences the continua of [CO2 N2 O2] molecules." ;

// global attributes:
		:_NCProperties = "version=1|netcdflibversion=4.4.1.1|hdf5libversion=1.10.1" ;
		:description = "Absorption coefficients for O2 as a function of pressure, temperature, H2O VMR, wavenumber, and band" ;
		:source = "Band 1: AER v3.6" ;
}
```

All other allowed molecules will conform to a similar convention, only without the H<sub>2</sub>O VMR dimension:

```
% ncdump -h nc_ABSCO/O3_00500-00600_v0.0_init.nc 
netcdf O3_00500-00600_v0.0_init {
dimensions:
	nfreq = 250001 ;
	nlev = 3 ;
	nlay = 2 ;
	ntemp = 15 ;
	nranges = 1 ;
	nranges_lims = 2 ;
variables:
	double P_level(nlev) ;
		P_level:_FillValue = NaN ;
		P_level:units = "mbar" ;
		P_level:long_name = "Pressure Levels" ;
		P_level:valid_range = 0., 1050. ;
		P_level:description = "User-provided layer boundary pressures" ;
	double P_layer(nlay, ntemp) ;
		P_layer:_FillValue = NaN ;
		P_layer:units = "mbar" ;
		P_layer:long_name = "Layer Pressures" ;
		P_layer:valid_range = 0., 1050. ;
		P_layer:description = "LBLRTM-calculated layer pressures" ;
	double Cross_Section(nfreq, ntemp, nlay) ;
		Cross_Section:_FillValue = NaN ;
		Cross_Section:units = "cm2/molecule" ;
		Cross_Section:long_name = "Absorption Coefficients" ;
		Cross_Section:valid_range = 0., 1.e-20 ;
		Cross_Section:description = "Absorption coefficients (K) calculated from LBLRTM optical depths and layer amounts." ;
	double Spectral_Grid(nfreq) ;
		Spectral_Grid:_FillValue = NaN ;
		Spectral_Grid:units = "cm-1" ;
		Spectral_Grid:long_name = "Spectral Points" ;
		Spectral_Grid:valid_range = 0., 50000. ;
		Spectral_Grid:description = "Spectral points corresponding to ABSCOs in a single layer for a single temperature and in a given spectral range (wavenumbers, microns, or nanometers)." ;
	double Temperature(nlev, ntemp) ;
		Temperature:_FillValue = NaN ;
		Temperature:units = "K" ;
		Temperature:long_name = "Temperature Levels" ;
		Temperature:valid_range = 180., 320. ;
		Temperature:description = "Applicable temperatures associated with each layer boundary pressure" ;
	double Extent_Ranges(nranges, nranges_lims) ;
		Extent_Ranges:_FillValue = NaN ;
		Extent_Ranges:units = "cm-1" ;
		Extent_Ranges:long_name = "Spectral Ranges" ;
		Extent_Ranges:valid_range = 0., 50000. ;
		Extent_Ranges:description = "Starting and ending spectral points for each input channel." ;
	int64 Extent_Indices(nranges, nranges_lims) ;
		Extent_Indices:_FillValue = 9223372036854775807L ;
		Extent_Indices:units = "N/A" ;
		Extent_Indices:long_name = "Spectral Ranges Reference Indices" ;
		Extent_Indices:valid_range = 0L, 9223372036854775807L ;
		Extent_Indices:description = "Pairs of indices defining the start and end indices of the Cross_Secion frequency dimension for non-continuous calculation regions." ;

// global attributes:
		:_NCProperties = "version=1|netcdflibversion=4.4.1.1|hdf5libversion=1.10.1" ;
		:description = "Absorption coefficients for O3 as a function of pressure, temperature, wavenumber, and band" ;
		:source = "Band 1: AER v3.6" ;
}
```

In the global attributes, there is a "source" field. There are only two sources -- HITRAN 2012 and AER LPD v3.6 -- and they can vary by band. The code accounts for the by-band differences in source.

## Reading the Output

We also designed a module to read the `ABSCO_compute.py` output. The following script reads in a file like the ones in the [Output](#output) section and prints out the cross section value at a given coordinate (it also prints out the array coordinates -- i.e., zero-offset coordinates -- of the absorption coefficient):

```
% read_ABSCO_tables.py -h
usage: read_ABSCO_tables.py [-h] [-p IN_PRESSURE] [-T IN_TEMP]
                            [-s IN_SPECTRAL IN_SPECTRAL] [-wv IN_H2O]
                            [-tol TOLERANCE]
                            ncFile

Read in netCDF generated with run_LBLRTM_ABSCO.py module and print out absorption
coefficient (k) at a given pressure, temperature, spectral point, and water
vapor amount (if molecule continuum is affected by water vapor).

positional arguments:
  ncFile                Output netCDF generated by run_LBLRTM_ABSCO.py.

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

