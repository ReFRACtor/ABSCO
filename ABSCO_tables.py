#!/usr/bin/env python

# standard Python libraries
# for Python 3 compatibility
from __future__ import print_function

import os, sys, glob, argparse
import subprocess as sub

# miniconda-installed libs
import numpy as np
import pandas as pd

# local module (part of the ABSCO library)
import ABSCO_preprocess as preproc

# path to GIT common submodules
sys.path.append('common')
import utils
import RC_utils as RC
import lblTools

def makeSymLinks(sources, targets):
  """
  Loop over input files and make symbolic links for them
  """

  for source, target in zip(sources, targets):
    # UNCOMMENT THIS WHEN NOT DEBUGGING
    utils.file_check(source)
    if os.path.exists(target): os.unlink(target)
    os.symlink(source, target)
  # end link loop

  return None
# end makeSymLinks()

class makeABSCO():
  """
  - Build TAPE5s (LNFL and LBLRTM) for each molecule of interest
  - Build TAPE3 (binary line file) for specified bands
  - Run LBLRTM to generate ODInt files
  - Use the ODInt files to calculate absorption coefficients and store
    in netCDF file

  Generate absorption coefficient tables (ABSCO as a function of 
  wavenumber, pressure, temperature, and band) for specified molecule
  """

  def __init__(self, inObj):

    """
    Inputs
      inObj -- preproc.configure instance

    Keywords
    """

    # UNCOMMENT THIS WHEN NOT DEBUGGING
    # make sure paths exist before proceeding
    for path in inObj.paths: utils.file_check(path)

    # make output directories
    for outDir in inObj.outDirs:
      if outDir == 'TAPE5_dir':
        # make a T5 subdir for LNFL and LBL
        lnflDirT5 = '%s/%s' % (inObj.lnfl_run_dir, outDir)
        if not os.path.exists(lnflDirT5): os.mkdir(lnflDirT5)

        lblDirT5 = '%s/%s' % (inObj.lbl_run_dir, outDir)
        if not os.path.exists(lblDirT5): os.mkdir(lblDirT5)
      else:
        if not os.path.exists(outDir): os.mkdir(outDir)
      # endif outDir
    # end outDir loop

    inP = np.array(inObj.pressures)

    # determine the temperature "levels" at each pressure boundary
    # [temperature] = K, [pressure] = mbar
    gridP, gridT1, gridT2 = np.loadtxt(inObj.ptfile, unpack=True)
    tLevList = []
    allT = []
    for pLev in inP:
      # find closest value to pLev in PT grid
      iLocP = np.argmin(np.abs(gridP-pLev))

      # generate an array of temperatures for this pressure in 10 K
      # increments
      gridT = np.arange(gridT1[iLocP], gridT2[iLocP]+10, 10)
      tLevList.append(gridT)

      # this will be used for determination of T dimension
      allT += list(gridT)
    # end P loop

    # grab the user provided VMR profile for the species of interest
    # TO DO: _HI and _XS species (so right now, this will not work 
    # for SO2, NO_2, HNO3, CLONO2, CH3CN, or CF4)
    userProfAll = pd.read_csv(inObj.vmrfile)
    userProf = {}
    userProf['P'] = userProfAll['P'].values
    userProf['T'] = userProfAll['T'].values
    userProf['BRD'] = userProfAll['BRD'].values
    for mol in inObj.molnames:
      if mol in userProfAll.keys().values:
        userProf[mol] = userProfAll[mol].values
      # endif mol
    # end mol loop

    # set class attributes
    # state, etc. atts
    self.headerOD = inObj.header
    self.cntnmScale = float(inObj.scale)
    self.pLev = np.array(inP)
    self.nP = inP.size
    self.tLev = list(tLevList)
    self.nT = np.unique(np.array(allT)).size
    self.pressures = np.array(inP)
    self.bands = dict(inObj.channels)
    self.nBands = len(inObj.channels['res'])
    self.molNames = list(inObj.molnames)
    self.doBand = dict(inObj.doBand)

    # grab the profiles
    self.vmrProf = dict(userProf)

    # LNFL atts
    self.runDirLNFL = str(inObj.lnfl_run_dir)
    self.pathLNFL = str(inObj.lnfl_path)
    self.pathT1 = str(inObj.tape1_path)
    self.pathT2 = str(inObj.tape2_path)
    self.dirExtras = str(inObj.extra_params)
    self.dirT3 = str(inObj.tape3_dir)

    # LBL atts
    self.runDirLBL = str(inObj.lbl_run_dir)
    self.pathLBL = str(inObj.lbl_path)
    self.pathXSDB = str(inObj.xs_path)
    self.pathListXS = str(inObj.fscdxs)
    self.dirT5 = str(inObj.tape5_dir)
    self.fineOD = str(inObj.od_dir)
    self.coarseOD = str(inObj.absco_dir)

    # all HITRAN molecule names (these are the molecules for which we
    # have line parameters)
    """
    # might be useful later...list the HITRAN molecule names
    lfMolDir = '/nas/project/rc_static/models/' + \
      'aer_line_parameters/AER_line_files/aer_v_3.6/' + \
      'line_files_By_Molecule/*'
    molDirs = sorted(glob.glob(lfMolDir))

    # the upper() takes care of the Br problem
    htMols = [os.path.basename(md).split('_')[1].upper() for \
      md in molDirs]
    print(htMols)
    """
    self.HITRAN = ['H2O', 'CO2', 'O3', 'N2O', 'CO', 'CH4', 'O2', \
      'NO', 'SO2', 'NO2', 'NH3', 'HNO3', 'OH', 'HF', 'HCL', 'HBR', \
      'HI', 'CLO', 'OCS', 'H2CO', 'HOCL', 'N2', 'HCN', 'CH3CL', \
      'H2O2', 'C2H2', 'C2H6', 'PH3', 'COF2', 'SF6', 'H2S', 'HCOOH', \
      'HO2', 'O', 'CLONO2', 'NO+', 'HOBR', 'C2H4', 'CH3OH', 'CH3BR', \
      'CH3CN', 'CF4', 'C4H2', 'HC3N', 'H2', 'CS', 'SO3']

    # for cd'ing back into the cwd
    self.topDir = os.getcwd()
  # end constructor()

  def lnflT5(self):
    """
    Make a TAPE5 for every band and molecule that can be used as 
    input into an LNFL run. We will only be using the extra broadening
    option (so we are not writing to an ASCII TAPE7 and are including 
    line coupling.
    """

    # this is part of record 3 (LNFL instructions), all molecules off
    molIndInit = np.repeat('0', 47)

    # the other part of record 3 -- let's always keep line coupling on
    # so we get the O2, CO2, and CH4 coupling params; suppress any 
    # output to an ASCII TAPE7; and always use extra broadening params
    holInd = 'EXBRD'
    for mol in self.molNames:
      outDirT5 = '%s/%s/%s' % (self.runDirLNFL, self.dirT5, mol)
      if not os.path.exists(outDirT5): os.mkdir(outDirT5)

      try:
        iMol = self.HITRAN.index(mol)
      except:
        print('Could not find %s in HITRAN names' % mol)
        continue
      # end exception

      for iBand in range(self.nBands):
        if self.doBand[mol][iBand] is False: continue
        # band specification	
        wvn1 = self.bands['wn1'][iBand]
        wvn2 = self.bands['wn2'][iBand]
        record1 = 'TAPE5 for %s, %.f-%.f' % (mol, wvn1, wvn2)
        record2 = '%10.3f%10.3f' % (wvn1, wvn2)

        # switch molecule "on" for LNFL
        molInd = np.array(molIndInit)
        molInd[iMol] = '1'
        record3 = '%47s%4s%40s' % (''.join(list(molInd)), ' ', holInd)

        # write LNFL TAPE5 for this band and molecule
        recs = [record1, record2, record3]

        # making WN1 and WN2 ints just to keep "." out of name
        outFile = '%s/TAPE5_%s_%05d-%05d' % \
          (outDirT5, mol, \
           self.bands['wn1'][iBand], self.bands['wn2'][iBand])
        outFP = open(outFile, 'w')
        for rec in recs: outFP.write('%s\n' % rec)
        outFP.close()
        print('Wrote %s' % outFile)
      # end band loop
    # end mol loop
  # end lnflT5()

  def runLNFL(self):
    """
    Run executable to generate a binary line file (TAPE3) for usage
    in LBLRTM.  Do this for each TAPE5 input available for a given 
    molecule.  Each TAPE3 contains a subset of the full TAPE1 line 
    file (i.e., only the lines for the given band and molecule).
    """

    # link to extra broadening and speed dependence parameters
    # equivalent to `ln -s full_path`
    os.chdir(self.runDirLNFL)
    extras = glob.glob('%s/*_param' % self.dirExtras)
    # UNCOMMENT THIS WHEN NOT DEBUGGING
    """
    if len(extras) == 0:
      print('No broadening or speed dependence parameters found')
      print('Returning')
      sys.exit(1)
    # endif extras
    """

    slExtras = [os.path.basename(extra) for extra in extras]
    makeSymLinks(extras, slExtras)

    # other LNFL links that will not change
    srcLNFL = [self.pathLNFL, self.pathT1]
    tarLNFL = ['lnfl', 'TAPE1']

    makeSymLinks(srcLNFL, tarLNFL)

    for mol in self.molNames:
      # line coupling molecules have other targets
      if mol in ['CO2', 'CH4', 'O2']:
        makeSymLinks([self.pathT2], ['TAPE2'])
      # endif TAPE2

      inDirT5 = '%s/%s' % (self.dirT5, mol)
      inT5 = sorted(glob.glob('%s/TAPE5_*' % inDirT5))
      outDirT3 = '%s/%s/%s' % (self.topDir, self.dirT3, mol)
      if not os.path.exists(outDirT3): os.mkdir(outDirT3)

      if len(inT5) == 0:
        print('Found no TAPE5s for %s' % mol)
        continue
      # endif nT5

      for t5 in inT5:
        print('Running LNFL for %s' % os.path.basename(t5))
        if os.path.islink('TAPE5'): os.unlink('TAPE5')
        # making some assumptions about file naming convention...
        split = os.path.basename(t5).split('_')
        band = split[-1]
        os.symlink(t5, 'TAPE5')
        sub.call(['./lnfl'])
        os.rename('TAPE3', '%s/TAPE3_%s_%s' % (outDirT3, mol, band))
      # end TAPE5 loop

    # end mol loop

    os.chdir(self.topDir)
    print(os.getcwd())
  # end runLNFL()

  def lblT5(self):
    """
    Make a TAPE5 for every temperature, band, and molecule that can 
    be used as input into an LBLRTM run

    see lblrtm_instructions.html for doc on each TAPE5 record
    """

    # multiplicative continuum factors (because CN=6 in record12)
    scales = np.repeat(0.0, 7)

    # records required with IATM=1 (provide some doc on this rec)
    # US Standard atmosphere, path type 2 (slant from H1 to H2), 2 
    # pressure levels, no zero-filling, full printout, 7 molecules,
    record31 = '%5d%5d%5d%5d%5d%5d' % (0, 2, -2, 0, 0, 7)

    for mol in self.molNames:
      doXS = 0 if mol in self.HITRAN else 1

      # TO DO: WORK ON THIS EVENTUALLY
      if doXS == 1: continue

      # record1.2: HI=9: central line contribution omitted
      # CN=6: continuum scale factor for given molecules used
      # OD=1, MG=1: optical depth computation, layer-by-layer
      record12 = ' HI=9 F4=0 CN=6 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 ' + \
        'AM=1 MG=1 LA=0 OD=1 XS=%1d' % doXS

      # record 3.4: user profile header for given molecule
      record34 = '%5d%24s' % (self.nP, 'User profile for %s' % mol)

      # this is for record 3.6
      n36 = len(self.HITRAN) if doXS == 0 else 7

      for iBand in range(self.nBands):

        if self.doBand[mol][iBand] is False: continue

        # continuum scale factors
        if mol == 'H2O':
          scales[:2] = self.cntnmScale
        elif mol == 'CO2':
          scales[2] = self.cntnmScale
        elif mol == 'O3':
          scales[3] = self.cntnmScale
        elif mol == 'O2':
          scales[4] = self.cntnmScale
        elif mol == 'N2':
          scales[5] = self.cntnmScale
        # endif CNTNM scales

        # generate free-format record 1.2
        record12a = ' '.join(scales.astype(str))

        outDirT5 = '%s/%s/%s' % (self.runDirLBL, self.dirT5, mol)
        if not os.path.exists(outDirT5): os.mkdir(outDirT5)

        for iP, pLev in enumerate(self.pLev):
          # need 2 P bounds for LBLATM
          if iP == 0: continue
          pArr = [self.pLev[iP-1], self.pLev[iP]]

          # record 3.2: pressure limits for all levels, nadir SZA
          record32 = '%10.3f%10.3f%10.3f' % (pArr[0], pArr[1], 0)

          # record 3.3b: pressure levels
          record33b = '%10.3f%10.3f' % (pArr[0], pArr[1])

          # record 3.5: level and unit info for record 3.6
          # using a fill space for "ZM" because whatever i would 
          # provide for that field would be ignored
          # really all we're doing is P and T units (in mbar and K)
          # and using the default (blank) format and units for profile
          # info (E10.3 VMR)
          record35 = '%10s%10.3E%10.3E%5sAA' % \
            ('', pLev, self.vmrProf['T'][iP], '')

          # record 3.6: provide profile info at a given level, but 
          # only for the broadener (density) and given mol (VMR)
          if doXS == 0:
            hitranAll = np.repeat(0.0, n36)

            # fill in the VMR for the given mol
            iMatch = self.HITRAN.index(mol)
            hitranAll[iMatch] = self.vmrProf[mol][iP]

            # insert the broadening density -- the eighth "molecule"
            hitranAll = np.insert(\
              hitranAll, 7, self.vmrProf['BRD'][iP])

            # start building the string for record 3.6
            record36 = ''
            for iDen, den in enumerate(hitranAll):
              record36 += '%10.3E' % den

              # eight molecules per line (but only 48 molecules, and 
              # no need for new line at end)
              if ((iDen+1) % 8) == 0 and iDen < 47: record36 += '\n'
            # end record36 loop
          else:
            # TO DO: XS record 3.7.1
            continue
          # end record36

          for tLev in self.tLev[iP]:
            outFile = '%s/TAPE5_%s_P%09.4fmb_T%05.1fK_%05d-%05d' % \
              (outDirT5, mol, pLev, tLev, self.bands['wn1'][iBand], \
               self.bands['wn2'][iBand])

            # TAPE5 records
            # record 1.3 is kinda long...first, band limits
            record13 = '%10.3e%10.3e' % \
              (self.bands['wn1'][0], self.bands['wn2'][0])

            # concatenate (NOT append) 6 zeros in scientific notation
            # using defaults for SAMPLE, DVSET, ALFAL0, AVMASS, 
            # DPTMIN, and DPTFAC params
            record13 += ''.join(['%10.3e' % 0 for i in range(6)])

            # line rejection not recorded and 1e-4 output OD spectral
            # resolution
            record13 += '%4s%1d%5s%10.3e' % \
              (' ', 0, ' ', self.bands['res'][0])

            # write the TAPE5 for this set of params
            recs = [record12, record12a, record13, record31, \
              record32, record33b, record34, record35, record36]

            if os.path.exists(outFile):
              print('WARNING: Overwriting %s' % outFile)

            outFP = open(outFile, 'w')
            outFP.write('$ %s\n' % self.headerOD)
            for rec in recs: outFP.write('%s\n' % rec)
            outFP.write('%%%%')
            outFP.close()
          # end temperature loop
        # end pressure loop
      # end band loop
    # end molecule loop

    return True
  # end lblT5()

  def runLBL(self):
    """
    Run LBLRTM for each TAPE5 made in lblT5

    This can be run in parallel for each molecule, but that means 
    each molecule should have its own LBL_Run directory (as specified
    in the input configuration file)
    """

    workSubDir = '%s/%s' % (self.topDir, self.runDirLBL)
    os.chdir(workSubDir)

    # aliases for symlinks; should correspond to lblFiles
    # these are identical for all LBL runs for this task
    targets = ['lblrtm', 'xs', 'FSCDXS']
    sources = [self.pathLBL, self.pathXSDB, self.pathListXS]
    makeSymLinks(sources, targets)

    for mol in self.molNames:
      outDirOD = '%s/%s/%s' % (self.topDir, self.fineOD, mol)
      if not os.path.exists(outDirOD): os.mkdir(outDirOD)

      # set up working subdirectory
      # find TAPE3s and their associated TAPE5s (for every TAPE3, 
      # there is a TAPE5 for every pressure and every temperature)
      searchStr = '%s/%s/%s/TAPE3*' % (self.topDir, self.dirT3, mol)
      molT3 = sorted(glob.glob(searchStr))
      if len(molT3) == 0:
        print('Found no TAPE3 files for %s' % mol)
        continue
      # endif nT3

      molT5 = []
      for t3 in molT3:
        base = os.path.basename(t3)
        band = base.split('_')[-1]
        searchStr = '%s/%s/%s/%s/TAPE5_*_%s' % \
          (self.topDir, self.runDirLBL, self.dirT5, mol, band)
        molT5 = sorted(glob.glob(searchStr))

        if len(molT5) == 0:
          print('Found no TAPE5s for %s' % mol)
          continue
        # endif t5

        if os.path.islink('TAPE3'): os.unlink('TAPE3')
        os.symlink(t3, 'TAPE3')

        for t5 in molT5:
          base = os.path.basename(t5)
          print(base)
          if os.path.islink('TAPE5'): os.unlink('TAPE5')
          os.symlink(t5, 'TAPE5')

          # grab extension for use in renaming the ODint LBL output 
          # file
          ext = base.replace('TAPE5_', '')
          sub.call(['./lblrtm'])
          odStr = 'ODint_001'

          # if all ODs are zero, remove the file and continue to next 
          # iteration (this saves HD space)
          freq, od = lblTools.readOD(odStr, double=True)
          if od.min() == 0 and od.max() == 0:
            os.remove(odStr)
            continue
          # endif zero OD

          os.rename(odStr, '%s/%s' % \
            (outDirOD, odStr.replace('001', ext)))
          break
        # end T5 loop
      # end t3 loop
    # end molecule loop

    return True
  # end runLBL()

# end makeABSCO()

if __name__ == '__main__':

  parser = argparse.ArgumentParser(\
    description='Generate ABSCO tables for user-specified molecules.')
  parser.add_argument('--config_file', type=str, \
    default='ABSCO_config.ini', \
    help='Configuration file that contains the file and ' + \
    'directory names necessary for the makeABSCO class.')

  # argument switches (booleans, no assignments)
  parser.add_argument('-lft5', '--lnfl_tape5', action='store_true', \
    help='Make the LBLRTM TAPE5 files for each species, channel, ' + \
    'pressure, and temperature.')
  parser.add_argument('-lnfl', '--run_lnfl', action='store_true', \
    help='Run LBLRTM for each of the TAPE5s generated and saved ' + \
    'in lblT5().')
  parser.add_argument('-lblt5', '--lbl_tape5', action='store_true', \
    help='Make the LBLRTM TAPE5 files for each species, channel, ' + \
    'pressure, and temperature.')
  parser.add_argument('-lbl', '--run_lbl', action='store_true', \
    help='Run LBLRTM for each of the TAPE5s generated and saved ' + \
    'in lblT5().')
  parser.add_argument('-e2e', '--end_to_end', action='store_true', \
    help='Runs the entire process from tape 5 generation to ' + \
    'post-processing (rather than entering all of the keywords ' + \
    'separately).')
  args = parser.parse_args()

  iniFile = args.config_file; utils.file_check(iniFile)

  # configuration object instantiation
  ini = preproc.configure(iniFile)

  # verification of molecules and associated bands to process
  # this should be called to limit unnecessary runs of LNFL and LBL
  ini.molProcess()

  for iniName in ini.molnames:
    if iniName in ini.dunno: sys.exit('Cannot do %s yet' % iniName)

  # ABSCO object instantiation
  absco = makeABSCO(ini)

  if args.run_lnfl: 
    absco.lnflT5()
    absco.runLNFL()
  # end LNFL

  if args.run_lbl:
    absco.lblT5()
    absco.runLBL()
  # end LBL

  # haven't tested this yet, but no reason it won't work...right?
  if args.end_to_end:
    sys.exit('No e2e yet')
    absco.lblT5()
    absco.runLBL()
  # endif e2e
# end main()

