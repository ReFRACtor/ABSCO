#!/usr/bin/env python

# standard Python libraries
# for Python 3 compatibility
from __future__ import print_function

import os, sys, argparse

# miniconda installs
import numpy as np
import pandas as pd
from pandas import DataFrame as DF

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
    from scipy import interpolate

    # user-specified pressures
    inP = np.loadtxt(self.pressureFile)

    # make sure they are going surface to TOA
    if np.diff(inP)[0] > 0: inP = inP[::-1]

    # import CSV as data frames
    hiDF = pd.read_csv(self.HITRAN)
    xsDF = pd.read_csv(self.XS)

    # only the first 7 HITRAN molecules, the broadening density, P, 
    # and T have VMR profiles that are not constant as a function of 
    # standard atmosphere, so only grab the columns corresponding to 
    # the requested stanAtm also keep the corresponding pressure grid
    stanAtmList = ['H2O', 'CO2', 'O3', 'CO', 'N2O', \
      'CH4', 'O2', 'BRD']
    stanAtmList = ['%s_%1d' % (param, self.stanAtm) for param in \
      stanAtmList]

    # P and T have different LBLATM naming conventions...
    ptList = ['%1s%1d' % (pt, self.stanAtm) for pt in ['P', 'T']]
    hiNames = ptList + stanAtmList + list(hiDF.keys().values[61:])
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
      if mol in ptList: mol = mol[0]

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

LBLDEFAULT = '/nas/project/rc_static/models/aer_lblrtm/' + \
  'lblrtm_v12.9/lblrtm_v12.9_linux_pgi_dbl'

class broadener():
  def __init__(self, inObj, lblPath=LBLDEFAULT):
    """
    After vmrProfiles() object is constructed, run LBLATM (subroutine 
    of LBLRTM that calculates density profiles) to compute the correct
    broadening densities when all only one specified molecule is "on"

    This class is dependent on vmrProfiles attributes to remind the 
    user that the broadening densities need to be recalculated every 
    time a new profile is provided

    inObj -- vmrProfiles object
    lblPath -- string, path to LBLRTM executable
    """

    utils.file_check(inObj.outFile)
    utils.file_check(lblPath)

    self.vmrObj = inObj
    self.workDir = os.getcwd()
    self.dirT5 = '%s/LBLATM/TAPE5_dir' % self.workDir
    self.dirT7 = '%s/LBLATM/TAPE7_dir' % self.workDir
    self.lblPath = str(lblPath)

    # check if output directories exist
    for outDir in ['%s/LBLATM' % self.workDir, self.dirT5, self.dirT7]:
      if not os.path.exists(outDir): os.mkdir(outDir)
    # end outDir loop

    self.outCSV = '%s/%s_broadener.csv' % \
      (self.workDir, inObj.outFile[:-4])

    csvDat = pd.read_csv(inObj.outFile)
    molNames = list(csvDat.keys().values)

    # "molecules" to skip
    rm = ['P', 'T', 'BRD', '??????']
    for mol in rm:
      if mol in molNames: molNames.remove(mol)
    # end mol loop

    self.molNames = list(molNames)

    # weird...record 3.5 only allows 39 molecules, but HITRAN has 47 
    # and that is the max allowed in record 2.1 (NMOL)
    self.molMaxLBL = 39

    # the band for LBLATM is arbitrary because no radiative transfer 
    # calculation is done (same with resolution)
    self.startWN = 4700.0
    self.endWN = 4800.00

    # pressures should be monotonically decreasing since we grabbed 
    # them from LBLATM
    self.pLev = csvDat['P'].values
    self.nP = self.pLev.size
    self.tLev = csvDat['T'].values
    self.vmr = csvDat
  # end constructor

  def makeT5(self):
    """
    Generate a small LBLRTM TAPE5 with the parameters necessary to 
    run LBLATM and generate a TAPE7 file. One for each molecule.
    """

    molHITRAN = ['H2O', 'CO2', 'O3', 'N2O', 'CO', 'CH4', 'O2', \
      'NO', 'SO2', 'NO2', 'NH3', 'HNO3', 'OH', 'HF', 'HCL', 'HBR', \
      'HI', 'CLO', 'OCS', 'H2CO', 'HOCL', 'N2', 'HCN', 'CH3CL', \
      'H2O2', 'C2H2', 'C2H6', 'PH3', 'COF2', 'SF6', 'H2S', 'HCOOH', \
      'HO2', 'O', 'CLONO2', 'NO+', 'HOBR', 'C2H4', 'CH3OH', 'CH3BR', \
      'CH3CN', 'CF4', 'C4H2', 'HC3N', 'H2', 'CS', 'SO3']

    # first generate all of the records needed for the TAPE5
    # records 1.2, 1.3, 3.1, and 3.2 are independent of molecule

    # record1.2: HI=0: HIRAC not activated, no LBL calculation used
    # CN=0: no continuum
    # OD=0, MG=0: no LBL optical depth computation, final layer
    record12 = ' HI=0 F4=0 CN=0 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 ' + \
      'AM=1 MG=0 LA=0 OD=0 XS=0'

    # record 1.3 is kinda long...first, band limits
    record13 = '%10.3e%10.3e' % (self.startWN, self.endWN)

    # concatenate (NOT append) 6 zeros in scientific notation
    # using defaults for SAMPLE, DVSET, ALFAL0, AVMASS, 
    # DPTMIN, and DPTFAC params
    record13 += ''.join(['%10.3e' % 0 for i in range(6)])

    # line rejection not recorded and 1e-4 output OD spectral
    # resolution
    record13 += '%4s%1d%5s%10.3e' % (' ', 0, ' ', 0)

    # records required with IATM=1 (provide some doc on this rec)
    # US Standard atmosphere, path type 2 (slant from H1 to H2), 2 
    # pressure levels, no zero-filling, full printout, 7 molecules, 
    # print to TAPE7
    record31 = '%5d%5d%5d%5d%5d%5d%5d' % \
      (0, 2, -self.nP, 0, 0, self.molMaxLBL, 1)

    # record 3.2: pressure limits, nadir SZA
    record32 = '%10.3f%10.3f%10.3f' % (self.pLev[0], self.pLev[-1], 0)

    # record 3.3b: pressure boundaries at each level
    record33 = ''
    for iP, p in enumerate(self.pLev):
      record33 += '%10.3f' % p

      # eight pressures per line (and no need for new line at end)
      if ((iP+1) % 8) == 0 and iP < self.nP: record33 += '\n'
    # end record36 loop

    # for record 3.5
    strForm = ['A', 'A', ' ', ' '] + ['A'] * self.molMaxLBL
    strForm = ''.join(strForm)
    altSfc = '%10.3E' % 0

    for mol in self.molNames:
      # records 1.1 and 3.4 all depend on molecule
      record11 = '$ LBLATM run for %s, broadener calc' % mol

      # record 3.4: user profile header for given molecule
      record34 = '%5d%24s' % (-self.nP, 'User profile for %s' % mol)

      # some molecules have line parameters and XS, use the former
      # and don't do any XS
      if '_XS' in mol: continue
      hiMol = mol.split('_')[0] if '_HI' in mol else str(mol)
      if mol == 'NOPLUS': hiMol = 'NO+'
      if hiMol not in molHITRAN: continue

      recs = [record11, record12, record13, record31, record32, \
        record33, record34]

      for iP, p, t in zip(np.arange(self.nP), self.pLev, self.tLev):
        # record 3.5: level and unit info for record 3.6
        # using a fill space for "ZM" because whatever i would 
        # provide for that field would be ignored
        # really all we're doing is P and T units (in mbar and K)
        # and using the default (blank) format and units for profile
        # info (E10.3 VMR)
        alt = str(altSfc) if p == self.pLev[0] else ''
        record35 = '%10s%10.3E%10.3E%5s%s' % (alt, p, t, '', strForm)
        recs.append(record35)

        # record 3.6: provide profile info at a given level, but 
        # only for the broadener (density) and given mol (VMR)
        lblAll = np.repeat(0.0, self.molMaxLBL)

        # fill in the VMR for the given mol (needs to be in ppmv)
        iMatch = molHITRAN.index(hiMol)
        if iMatch >= self.molMaxLBL: continue
        lblAll[iMatch] = self.vmr[mol][iP] * 1e6

        # insert fill broadening density -- the eighth "molecule"
        lblAll = np.insert(lblAll, 7, 0)

        # start building the string for record 3.6
        record36 = ''
        for iVMR, vmr in enumerate(lblAll):
          record36 += '%10.3E' % vmr

          # eight molecules per line (but only 48 molecules, and 
          # no need for new line at end)
          if ((iVMR+1) % 8) == 0 and iVMR < self.molMaxLBL:
            record36 += '\n'
        # end record36 loop

        recs.append(record36)

      # end level loop

      # write the TAPE5 for a given molecule
      outFile = '%s/%s_TAPE5' % (self.dirT5, mol)
      outFP = open(outFile, 'w')
      for rec in recs: outFP.write('%s\n' % rec)
      outFP.write('%%%%%%\n')
      outFP.close()
      print('Wrote %s' % outFile)
    # end mol loop
  # end makeT5()

  def runLBLATM(self):
    """
    Run LBLATM with the TAPE5s generated in makeT5(), save TAPE7s in 
    their own directory
    """
  # end runLBLATM()

  def writeCSV(self):
    """
    Write a CSV file that consolidates all of the LBLATM broadener 
    information in it (density as a function of molecule and level)
    """
  # end runLBLATM
# end broadener

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
  parser.add_argument('-b', '--broad', action='store_true', \
    help='After calculating the VMR profiles, calculate the ' + \
    'broadening parameter for every layer and every molecule.  ' + \
    'This calculation is done by turning "on" one molecule at a ' + \
    'time.  Output is written to a variant of outfile.')
  parser.add_argument('-lbl', '--lbl_path', default=LBLDEFAULT, \
    help='Full path to LBLRTM executable that will be used if ' + \
    '--broad is used.')
  args = parser.parse_args()

  hiCSV = args.csvHITRAN; xsCSV = args.csvXS; pFile = args.pressures
  for f in [hiCSV, xsCSV, pFile]: utils.file_check(f)

  vmrProf = vmrProfiles(hiCSV, xsCSV, pFile, \
    outFile=args.outfile, stanAtm=args.standard_atm)
  vmrProf.calcVMR()

  if args.broad:
    broadObj = broadener(vmrProf, lblPath=args.lbl_path)
    broadObj.makeT5()
    broadObj.runLBLATM()
    broadObj.writeCSV()
  # endif broad
# endif main()

