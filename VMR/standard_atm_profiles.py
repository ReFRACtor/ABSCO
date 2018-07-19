#!/usr/bin/env python

# standard Python libraries
# for Python 3 compatibility
from __future__ import print_function

import os, sys, glob, argparse

# miniconda installs
import numpy as np
import pandas as pd
from pandas import DataFrame as DF

sys.path.append('../common')
import utils

LBLDEFAULT = '/nas/project/rc_static/models/aer_lblrtm/' + \
  'lblrtm_v12.9/lblrtm_v12.9_linux_pgi_dbl'
XSDEFAULT = '/nas/project/rc_static/models/aer_line_parameters/' + \
  'AER_line_files/aer_v_3.6/xs_files_v3.6'

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

class broadener():
  def __init__(self, inObj, lblPath=LBLDEFAULT, \
    doXS=False, xsPath=XSDEFAULT, debug=False):
    """
    After vmrProfiles() object is constructed, run LBLATM (subroutine 
    of LBLRTM that calculates density profiles) to compute the correct
    broadening densities when all only one specified molecule is "on"

    This class is dependent on vmrProfiles attributes to remind the 
    user that the broadening densities need to be recalculated every 
    time a new profile is provided

    inObj -- vmrProfiles object
    lblPath -- string, path to LBLRTM executable
    doXS -- boolean, instead of calculating broadener for each 
      molecule in inObj.outFile, just do one profile that will be used
      for all XS species
    xsPath -- string, path to directory with FSCDXS file and xs/ dir
    debug -- boolean, does not run LBLRTM
    """

    utils.file_check(inObj.outFile)
    utils.file_check(lblPath)

    if doXS:
      utils.file_check('%s/FSCDXS' % xsPath)
      utils.file_check('%s/xs' % xsPath)
    # endif XS

    self.vmrObj = inObj
    self.topDir = os.getcwd()
    self.workDir = '%s/LBLATM' % os.getcwd()
    self.dirT5 = '%s/TAPE5_dir' % self.workDir
    self.dirT7 = '%s/TAPE7_dir' % self.workDir
    self.lblPath = str(lblPath)
    self.lblExe = 'lblrtm'
    self.xsPath = str(xsPath)
    self.debug = bool(debug)

    # check if output directories exist
    for outDir in [self.workDir, self.dirT5, self.dirT7]:
      if not os.path.exists(outDir): os.mkdir(outDir)
    # end outDir loop

    self.outCSV = '%s/%s_broadener.csv' % \
      (self.topDir, inObj.outFile[:-4])

    csvDat = pd.read_csv(inObj.outFile)
    molNames = list(csvDat.keys().values)

    # "molecules" to skip
    rm = ['P', 'T', 'BRD', '??????']
    for mol in rm:
      if mol in molNames: molNames.remove(mol)
    # end mol loop

    self.molNames = ['XS'] if doXS else list(molNames)
    self.molMaxLBL = 39 if doXS else 47
    self.doXS = bool(doXS)

    # the band for LBLATM is arbitrary because no radiative transfer 
    # calculation is done (same with resolution)
    self.startWN = 4700.00
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
    xsStr = 1 if self.doXS else 0
    record12 = ' HI=0 F4=0 CN=0 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 ' + \
      'AM=1 MG=0 LA=0 OD=0 XS=%1d' % xsStr

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
    nMol = 7 if self.doXS else self.molMaxLBL
    record31 = '%5d%5d%5d%5d%5d%5d%5d' % \
      (0, 2, -self.nP, 0, 0, nMol, 1)

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
    strForm = ['A', 'A', ' ', ' '] + ['A'] * nMol
    strForm = ''.join(strForm)
    altSfc = '%10.3E' % 0

    allT5 = []
    for mol in self.molNames:
      # records 1.1 and 3.4 all depend on molecule
      record11 = '$ LBLATM run for %s, broadener calc' % mol

      # record 3.4: user profile header for given molecule
      record34 = '%5d%24s' % (-self.nP, 'User profile for %s' % mol)

      recs = [record11, record12, record13, record31, record32, \
        record33, record34]

      if not self.doXS:
        # some molecules have line parameters and XS, use the former
        # and don't do the XS ('XS' will be separate)
        if '_XS' in mol: continue
        hiMol = mol.split('_')[0] if '_HI' in mol else str(mol)
        if mol == 'NOPLUS': hiMol = 'NO+'
        if hiMol not in molHITRAN: continue

        # for record 3.6
        iMatch = molHITRAN.index(hiMol)
        if iMatch >= self.molMaxLBL: continue
      # endif XS

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
        # fill in the VMR for the given mol (needs to be in ppmv)
        lblAll = np.repeat(0.0, nMol)
        if not self.doXS: lblAll[iMatch] = self.vmr[mol][iP] * 1e6

        # insert fill broadening density -- the eighth "molecule"
        lblAll = np.insert(lblAll, 7, 0)

        # start building the string for record 3.6
        record36 = ''
        for iVMR, vmr in enumerate(lblAll):
          record36 += '%10.3E' % vmr

          # eight molecules per line, new line every 8 
          # (but only 48 molecules, and no need for new line at end)
          if ((iVMR+1) % 8) == 0 and iVMR < nMol: record36 += '\n'
        # end record36 loop

        recs.append(record36)

      # end level loop

      if self.doXS:
        # let's just use a single molecule from an LBLATM-determined
        # standard  atmosphere rather than providing a user profile 
        # (like the line parameter molecules). the choice of standard
        # atmosphere is not crucial because we're just interested in
        # the broadener, and the XS densities are relatively small

        # record 3.7
        record37 = '%5d%5d%5d' % (1, 1, 0)
        recs.append(record37)

        # record 3.7.1
        record371 = 'F11'
        recs.append(record371)
      # endif XS

      # write the TAPE5 for a given molecule
      if '_HI' in mol: mol = mol.replace('_HI', '')
      outFile = '%s/%s_TAPE5' % (self.dirT5, mol)

      if not self.debug:
        outFP = open(outFile, 'w')
        for rec in recs: outFP.write('%s\n' % rec)
        outFP.write('%%%%%%\n')
        outFP.close()
        print('Wrote %s' % outFile)
      # endif debug

      allT5.append(outFile)
    # end mol loop

    self.allT5 = list(allT5)

    return self
  # end makeT5()

  def runLBLATM(self):
    """
    Run LBLATM with the TAPE5s generated in makeT5(), save TAPE7s in 
    their own directory
    """

    # standard Python library
    import subprocess as sub

    os.chdir(self.workDir)
    lblExe = self.lblExe
    if not os.path.islink(lblExe): os.symlink(self.lblPath, lblExe)

    if self.doXS:
      targets = ['FSCDXS', 'xs']
      sources = ['%s/%s' % (self.xsPath, tar) for tar in targets]
      for src, tar in zip(sources, targets):
        if not os.path.islink(tar): os.symlink(src, tar)
    # endif XS

    allT7 = []
    for t5 in self.allT5:
      if os.path.islink('TAPE5'): os.unlink('TAPE5')
      os.symlink(t5, 'TAPE5')

      # for the output file, basically replace "TAPE5" with "TAPE7"
      outFile = '%s/%s' % \
        (self.dirT7, os.path.basename(t5).replace('TAPE5', 'TAPE7'))

      if not self.debug:
        status = sub.run(['./%s' % lblExe])
        if status.returncode:
          print('LBLRTM did not work for %s' % os.path.basename(t5))
          continue
        # endif returncode

        os.rename('TAPE7', outFile)
        print('Wrote %s' % outFile)
      # endif debug

      allT7.append(outFile)
    # end TAPE5 loop

    os.chdir(self.topDir)

    self.allT7 = list(allT7)

    return self
  # end runLBLATM()

  def writeCSV(self, supplementT7=None):
    """
    Write a CSV file that consolidates all of the LBLATM broadener 
    information in it (density as a function of molecule and level)

    supplement -- list, additional TAPE7 files to include in the CSV.
      this is helpful since XS and HITRAN molecules are processed in 
      separate broadener objects
    """

    # ABSCO submodule (in ../common)
    import RC_utils as RC

    allT7 = list(self.allT7) if supplementT7 is None else \
      self.allT7 + supplementT7

    brdDict = {}
    for t7 in allT7:
      mol = os.path.basename(t7).split('_')[0]
      brdDict[mol] = RC.readTAPE7(t7)['broadener']
    # end TAPE7 loop

    if not self.debug:
      DF.from_dict(brdDict).to_csv(\
        self.outCSV, float_format='%10.3E', index=False, na_rep='nan')
      print('Wrote %s' % self.outCSV)
    # endif debug
  # end writeCSV()
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
  parser.add_argument('-x', '--csvXS', type=str, \
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
    '--broad is specified.')
  parser.add_argument('-xs', '--xs_path', default=XSDEFAULT, \
    help='Full path to directory with FSCDX file and xs ' + \
    'subdirectory that will be used if --broad is specified.')
  parser.add_argument('-t', '--test', action='store_true', \
    help='Assumes script has already been run once and thus ' + \
    'the necessary TAPE5s and TAPE7s have been generated.  ' + \
    'Useful for debugging.')
  args = parser.parse_args()

  hiCSV = args.csvHITRAN; xsCSV = args.csvXS; pFile = args.pressures
  for f in [hiCSV, xsCSV, pFile]: utils.file_check(f)

  vmrProf = vmrProfiles(hiCSV, xsCSV, pFile, \
    outFile=args.outfile, stanAtm=args.standard_atm)
  vmrProf.calcVMR()

  if args.broad:
    print('Calculating broadening density profiles')
    # HITRAN molecules object
    broadObj = broadener(vmrProf, lblPath=args.lbl_path, \
      debug=args.test)
    broadObj.makeT5()
    broadObj.runLBLATM()

    # XS object
    xsBroadObj = broadener(vmrProf, lblPath=args.lbl_path, \
      debug=args.test, doXS=True, xsPath=args.xs_path)
    xsBroadObj.makeT5()
    xsBroadObj.runLBLATM()

    # there is probably a better way to combine the TAPE7s than what 
    # i am currently doing, but for now the it gets the job done
    broadObj.writeCSV(supplementT7=xsBroadObj.allT7)
  # endif broad
# endif main()

