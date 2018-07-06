#!/usr/bin/env python

# standard Python libraries
# for Python 3 compatibility
from __future__ import print_function

import os, sys, argparse

sys.path.append('../common')
import utils

class vmrProfiles():
  def __init__(self, hiFile, xsFile, pressureFile, \
    outFile='vmr_profiles.csv', stanAtm=6):
    """
    Read in CSV data for both the HITRAN and XS molecules and perform
    a linear interpolation onto the user-specified pressure grid.  
    The output file (CSV) will contain all molecules together

    Standard atmosphere index (int) can also be specified.
    From LBLATM:
      1 Tropical
      2 Mid-lat summer
      3 Mid-lat winter
      4 Subarctic summer
      5 Subarctic winter
      6 US Standard, 1976
    """

    self.HITRAN = str(hiFile)
    self.XS = str(xsFile)
    self.pressureFile = str(pressureFile)
    self.outFile = str(outFile)
    self.stanAtm = int(stanAtm)
  # end constructor

  def calcVMR(self):
    """
    Calculate VMR profile on user-specified pressure grid
    """

    # miniconda installs
    import pandas as pd
    from pandas import DataFrame as DF
    import numpy as np
    from scipy import interpolate

    # user-specified pressures
    inP = np.loadtxt(self.pressureFile)

    # make sure they are going surface to TOA
    if np.diff(inP)[0] > 0: inP = inP[::-1]

    # import CSV as data frames
    hiDF = pd.read_csv(self.HITRAN)
    xsDF = pd.read_csv(self.XS)

    # only the first 7 HITRAN molecules have VMR profiles that 
    # are not constant as a function of standard atmosphere, so only 
    # grab the columns corresponding to the requested stanAtm
    # also keep the corresponding pressure grid
    first7 = ['H2O', 'CO2', 'O3', 'CO', 'N2O', 'CH4', 'O2']
    first7 = ['%s_%1d' % (f7, self.stanAtm) for f7 in first7]
    hiNames = first7 + list(hiDF.keys().values[61:])
    xsNames = xsDF.keys().values
    stanAtmP = hiDF['P%1d' % self.stanAtm].values

    # now loop over each variable assemble interpolation function
    outDict = {}
    for mol in hiNames:
      # interpolate linearly onto input P grid, 
      # extrapolating where required
      vmrInterp = interpolate.interp1d(stanAtmP, hiDF[mol].values, \
        kind='linear', fill_value='extrapolate')
      vmrLBL = vmrInterp(inP)

      # for output, we can remove the standard atmosphere index substr
      if '_%1d' % self.stanAtm in mol: mol = mol.split('_')[0]

      # handle molecules that have line parameters and XS values
      if mol in xsNames: mol += '_HI'

      outDict[mol] = vmrLBL
    # end HITRAN mol loop

    # now do the same thing for the XS molecules
    for mol in xsNames:
      if mol == 'ALTX': continue

      vmrInterp = interpolate.interp1d(stanAtmP, xsDF[mol].values, \
        kind='linear', fill_value='extrapolate')
      vmrLBL = vmrInterp(inP)

      # handle molecules that have line parameters and XS values
      if mol in hiNames: mol += '_XS'

      outDict[mol] = vmrLBL
    # end XS mol loop

    # convert dictionary to data frame and write the CSV
    DF.from_dict(outDict).to_csv(\
      self.outFile, float_format='%10.3E', index=False, na_rep='nan')
    print('Wrote %s' % self.outFile)
  # end calcVMR()
  
# end vmrProfiles()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(\
    description='Utilize the data blocks from LBLATM for ' + \
    'standard atmosphere VMR profiles, then ' + \
    'interpolate/extrapolate to the desired pressure grid.')
  parser.add_argument('-hi', '--csvHITRAN', type=str, \
    default='LBLATM_Standard_Profiles.csv', \
    help='CSV file that contains LBLATM VMR profile blocks ' + \
    'for all 6 standard atmospheres.  The VMRs exist for ' + \
    'all HITRAN molecules.')
  parser.add_argument('-xs', '--csvXS', type=str, \
    default='XS_LBLATM_Standard_Profiles.csv', \
    help='CSV file that contains LBLATM VMR profile blocks ' + \
    'for all 6 standard atmospheres.  The VMRs exist for ' + \
    'all XS molecules (i.e., those where no HITRAN line ' + \
    'parameters exist.')
  parser.add_argument('-p', '--pressures', type=str, \
    default='../PT_grid/AIRS_P_air.txt', \
    help='List of pressures to use (1 pressure per line).')
  parser.add_argument('--standard_atm', type=int, default=6, \
    help='Index of LBLRTM standard atmosphere (TRP, USS, SAW, etc.)')
  parser.add_argument('-o', '--outfile', type=str, \
    default='vmr_profiles.csv', \
    help='Name of output CSV file.')
  args = parser.parse_args()

  hiCSV = args.csvHITRAN; xsCSV = args.csvXS; pFile = args.pressures
  for f in [hiCSV, xsCSV, pFile]: utils.file_check(f)

  vmrProf = vmrProfiles(hiCSV, xsCSV, pFile, \
    outFile=args.outfile, stanAtm=args.standard_atm)
  vmrProf.calcVMR()
# endif main()

