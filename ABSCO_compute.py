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

  def __init__(self, inObj, inMol, debug=False, \
    vmrWV=None, vmrO2=None):

    """
    Inputs
      mol -- string, name of molecule to process
      inObj -- preproc.configure instance

    Keywords
      debug -- boolean, only for testing purposes (probably a bit 
        obsolete as well, except that it does not do a loop over all 
        TAPE5s)
      vmrWV -- float, water vapor mixing ratio [ppmv] to be used in 
        H2O, CO2, and N2 profiles 
      vmrO2 -- float, oxygen mixing ratio [ppmv] to be used in 
        O2 profiles 
    """

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

    # for cd'ing back into the directories with the Git repo
    self.gitDir = os.getcwd()

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
    inUserProf = pd.read_csv(inObj.vmrfile)
    userProf = {}
    userProf['P'] = inUserProf['P'].values
    userProf['T'] = inUserProf['T'].values
    userProf['ALT'] = inUserProf['ALT'].values
    allMolCSV = inUserProf.keys().values

    for mol in inObj.molnames:
      # extract VMRs for each specified molecule
      if mol in allMolCSV:
        userProf[mol] = inUserProf[mol].values
      else:
        # This should be identical to "if mol in self.xsLines" 
        # NO2 and SO2 have identical XS and HITRAN densities
        # HNO3 and CF4 do not. whether we use the XS or HITRAN density
        # is determined in lblT5()
        hiMol = '%s_HI' % mol
        xsMol = '%s_XS' % mol

        if hiMol in allMolCSV:
          userProf[hiMol] = inUserProf[hiMol].values

        if xsMol in allMolCSV:
          userProf[xsMol] = inUserProf[xsMol].values
      # endif mol

      # this is the only special case so far: allowed molecule that 
      # has an LBLRTM alias that is not the allowed string
      if mol == 'F22': userProf[mol] = inUserProf['CHCLF2'].values
    # end mol loop

    # set class attributes
    # state, etc. atts
    self.headerOD = str(inObj.header)
    self.cntnmScale = float(inObj.scale)
    self.pLev = np.array(inP)
    self.nP = inP.size
    self.tLev = list(tLevList)
    self.allT = np.unique(np.array(allT))
    self.nT = self.allT.size
    self.pressures = np.array(inP)
    self.bands = dict(inObj.channels)
    self.nBands = len(inObj.channels['lblres'])
    self.molNames = list(inObj.molnames)
    self.doBand = dict(inObj.doBand)
    self.degradeKern = list(inObj.kernel)
    self.spectralUnits = str(inObj.spectral_units)
    self.dataSource = dict(inObj.sources)

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
    self.doXS = dict(inObj.doXS)
    self.molH2O = list(inObj.molH2O)
    self.molMaxLBL = 47
    self.vmrWV = float(vmrWV) if inMol in self.molH2O else None
    self.vmrArrWV = np.array(inObj.wv_vmr)

    # with O2, the VMR array is used with the netCDF; single O2 values
    # are used with each ABSCO run
    self.vmrArrO2 = np.array(inObj.vmrO2)
    self.vmrO2 = float(vmrO2) if inMol == 'O2' else None

    self.HITRAN = ['H2O', 'CO2', 'O3', 'N2O', 'CO', 'CH4', 'O2', \
      'NO', 'SO2', 'NO2', 'NH3', 'HNO3', 'OH', 'HF', 'HCL', 'HBR', \
      'HI', 'CLO', 'OCS', 'H2CO', 'HOCL', 'N2', 'HCN', 'CH3CL', \
      'H2O2', 'C2H2', 'C2H6', 'PH3', 'COF2', 'SF6', 'H2S', 'HCOOH', \
      'HO2', 'O', 'CLONO2', 'NO+', 'HOBR', 'C2H4', 'CH3OH', 'CH3BR', \
      'CH3CN', 'CF4', 'C4H2', 'HC3N', 'H2', 'CS', 'SO3']

    try:
      self.iMol = self.HITRAN.index(mol)
    except:
      print('Could not find %s in HITRAN names' % mol)
      self.iMol = None
    # end exception

    # XS species and molecules with XS and line parameters
    self.xsNames = list(inObj.xsNames)
    self.xsLines = list(inObj.xsLines)

    # for final output netCDF
    self.version = str(inObj.sw_ver)
    self.runDesc = str(inObj.out_file_desc)
    self.outDir = str(inObj.outdir)
    self.compress = int(inObj.nc_compress)
    self.freq_chunk = int(inObj.freq_chunk)

    self.debug = bool(debug)
  # end constructor()

  def lnflT5(self, mol):
    """
    For a given molecule, make a TAPE5 that can be used as input into
    an LNFL run for every band. We will only be using the extra 
    broadening option (so we are not writing to an ASCII TAPE7 and 
    are including line coupling).
    """

    # this is part of record 3 (LNFL instructions), all molecules off
    molIndInit = np.repeat('0', 47)

    # the other part of record 3 -- let's always keep line coupling on
    # so we get the O2, CO2, and CH4 coupling params; suppress any 
    # output to an ASCII TAPE7; and always use extra broadening params
    holInd = 'EXBRD'

    outDirT5 = '%s/%s/%s' % (self.runDirLNFL, self.dirT5, mol)
    if not os.path.exists(outDirT5): os.mkdir(outDirT5)

    for iBand in range(self.nBands):
      if self.doBand[mol][iBand] is False: continue
      # band specification	
      wvn1 = self.bands['wn1'][iBand]
      wvn2 = self.bands['wn2'][iBand]

      # spectral range should go +/- 25 cm-1 to incorporate 
      # non-negligible contributions from all "nearby" lines
      # and LNFL is smart enough to work with wvn1 < 0
      record1 = 'TAPE5 for %s, %.f-%.f' % (mol, wvn1, wvn2)
      record2 = '%10.3f%10.3f' % (wvn1-25, wvn2+25)

      # switch molecule "on" for LNFL; we still need a TAPE3 for LBL
      # with XS molecules even though they have no molecule number, 
      # so just don't turn on any molecules
      molInd = np.array(molIndInit)
      if self.iMol is not None: molInd[self.iMol] = '1'
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
    if len(extras) == 0:
      print('No broadening or speed dependence parameters found')
      print('Returning')
      sys.exit(1)
    # endif extras

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
      outDirT3 = '%s/%s' % (self.dirT3, mol)
      if not os.path.exists(outDirT3): os.mkdir(outDirT3)

      if len(inT5) == 0:
        print('Found no LBL TAPE5s for %s' % mol)
        continue
      # endif nT5

      for t5 in inT5:
        # making some assumptions about file naming convention...
        split = os.path.basename(t5).split('_')
        band = split[-1]
        outT3 = '%s/TAPE3_%s_%s' % (outDirT3, mol, band)
        if os.path.exists(outT3): continue

        print('Running LNFL for %s' % os.path.basename(t5))
        if os.path.islink('TAPE5'): os.unlink('TAPE5')
        os.symlink(t5, 'TAPE5')
        sub.call(['./lnfl'])
        os.rename('TAPE3', outT3)
      # end TAPE5 loop

    # end mol loop

    os.chdir(self.gitDir)
  # end runLNFL()

  def lblT5(self, mol):
    """
    For a given molecule, make a TAPE5 for every temperature and band 
    that can be used as input into an LBLRTM run

    see lblrtm_instructions.html for doc on each TAPE5 record

    mol -- string, molecule name from preproc.readConfig.allowed
    """

    pLevArr = np.array(self.pLev)

    # even though we are inputing an entire user profile, we will only
    # process 1 layer (2 levels) at a time
    nLevT5 = 2

    # multiplicative continuum factors (because CN=6 in record12)
    scales = np.repeat(0.0, 7)

    # records required with IATM=1 (we're using LBLATM to calculate
    # layer amounts for species -- we only have level amounts)
    # US Standard atmosphere, path type 2 (slant from H1 to H2), 2 
    # pressure levels, no zero-filling, full printout, 7 molecules, 
    # write to TAPE7
    record31 = '%5d%5d%5d%5d%5d%5d%5d' % \
      (0, 2, -nLevT5, 0, 0, self.molMaxLBL, 1)

    # record 3.4: user profile header for given molecule
    record34 = '%5d%24s' % (-self.nP, 'User profile for %s' % mol)

    for iBand in range(self.nBands):

      if self.doBand[mol][iBand] is False: continue

      if mol in self.xsNames:
        # omit lines
        doXS, optHI, optF4 = 1, 9, 0
      elif mol in self.xsLines and self.doXS[mol][iBand]:
        doXS, optHI, optF4 = 1, 9, 0
      else:
        # use lines +/- 25 cm-1
        doXS, optHI, optF4 = 0, 1, 1
      # endif doXS

      # record1.2: HI, F4: spectral line application
      # CN=6: continuum scale factor for given molecules used
      # OD=1, MG=1: optical depth computation, layer-by-layer
      record12 = ' HI=%1d F4=%1d CN=6 AE=0 EM=0 SC=0 FI=0 PL=0 ' % \
        (optHI, optF4) + 'TS=0 AM=1 MG=1 LA=0 OD=1 XS=%1d' % doXS
      record12 += '%20s' % '1'

      # continuum scale factors
      if mol == 'H2O':
        scales[0] = self.cntnmScale
        scales[1] = self.cntnmScale
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

      # record 1.3 is kinda long...first, band limits
      record13 = '%10.3e%10.3e' % \
        (self.bands['wn1'][iBand], self.bands['wn2'][iBand])

      # concatenate (NOT append) 6 zeros in scientific notation
      # using defaults for SAMPLE, DVSET, ALFAL0, AVMASS, 
      # DPTMIN, and DPTFAC params
      record13 += ''.join(['%10.3e' % 0 for i in range(6)])

      # line rejection not recorded and output OD spectral resolution
      record13 += '%4s%1d%5s%10.3e' % \
        (' ', 0, ' ', self.bands['lblres'][0])

      outDirT5 = '%s/%s/%s' % (self.runDirLBL, self.dirT5, mol)
      if not os.path.exists(outDirT5):
        os.makedirs(outDirT5, exist_ok=True)

      for constantT in self.allT:

        # generate entire user profile for given T and write to TAPE5
        records35_36, records381_382 = [], []
        for iP, pLev in enumerate(pLevArr):
          # records 3.5 and 3.6 should be repeated for every level
          # record 3.5: level and unit info for record 3.6
          # really all we're doing is P and T units (in mbar and K)
          # and using the default (blank) format and units for profile
          # info (E10.3 VMR)
          record35 = '%10.3E%10.3E%10.3E%5sAA\n' % \
            (self.vmrProf['ALT'][iP], pLev, constantT, '')
          records35_36.append(record35)

          # now determine the density to use for the level
          if mol in self.xsLines:
            # handle the "double agents" -- HITRAN and XS params are 
            # available, and density profiles for each are stored
            if doXS:
              levVMR = self.vmrProf['%s_XS' % mol][iP]
            else:
              levVMR = self.vmrProf['%s_HI' % mol][iP]
            # endif doXS
          else:
            levVMR = self.vmrProf[mol][iP]
          # endif doXS

          # eventually decided to use VMR in ppmv
          levVMR *= 1e6

          # record 3.6: provide VMR at a given level, but only for the 
          # given mol (VMR)
          lblAll = np.repeat(0.0, self.molMaxLBL)

          # fill in the VMR for the given mol
          if not doXS:
            lblAll[self.iMol] = float(levVMR)
            if mol in self.molH2O: lblAll[0] = self.vmrWV
          # endif no XS

          # start building the string for record 3.6
          record36 = ''
          for iDen, den in enumerate(lblAll):
            record36 += '%10.3E' % den

            # eight molecules per line
            if ((iDen+1) % 8) == 0: record36 += '\n'
          # end record36 loop

          # end of record36, if nMol is not divisible by 8
          if record36[-1] != '\n' and iP != self.nP-1:record36 += '\n'
          records35_36.append(record36)

          # record 3.8.1: boundary pressure (really just rec 3.3b)
          # record 3.8.2: layer molecule VMR, which should have been
          # defined in when constructing records 3.5 and 3.6
          # can combine these two since we're only doing 1 XS 
          # molecule per layer
          if doXS: 
            records381_382.append('%10.3f%5s1\n%10.3E' % \
              (pLev, '', float(levVMR)) )
            if iP != self.nP-1: records381_382.append('\n')
          # end doXS
        # end pressure loop

        # combine the record lists into single strings
        records35_36 = ''.join(records35_36)

        if doXS:
          # record 3.7: 1 molecule, user-provided profile
          record37 = '%5d%5d' % (1, 0)

          # record 3.7.1: XS molecule name
          record371 = '%10s' % mol

          # record 3.8: n pressure levels, pressure used for "height"
          record38 = '%5d%5d %s User Profile' % (self.nP, 1, mol)

          # records 3.8.1 and 3.8.2: user profile for XS mol
          records381_382 = ''.join(records381_382)

          xsRecs = [record37, record371, record38, records381_382]
        # end record36

        # write the TAPE5 for this set of params
        recs = [record12, record12a, record13, record31, \
          record34, records35_36]

        if doXS: recs += xsRecs

        # start writing the TAPE5s, 1 layer (2 levels) at a time
        for iP, pLev in enumerate(pLevArr):
          # skip surface since we need two levels
          if iP == 0: continue

          # pressure levels might not span all temperatures
          if constantT not in self.tLev[iP]: continue

          # keep recs constant for this pLev loop for reusability
          finalRecs = list(recs)

          # record 3.2: observer pressure limits, nadir SZA
          record32 = '%10.3f%10.3f%10.3f' % (pLevArr[iP-1], pLev, 0)
          finalRecs.insert(4, record32)

          record33b = '%10.3f%10.3f' % (pLevArr[iP-1], pLev)

          outFile = '%s/TAPE5_%s_P%09.4f_T%05.1fK_%05d-%05d' % \
            (outDirT5, mol, pLev, constantT, \
             self.bands['wn1'][iBand], self.bands['wn2'][iBand])

          finalRecs.insert(5, record33b)

          # throw the water vapor records into the appropriate spot
          if mol in self.molH2O: outFile = '%s_vmrWV%06.0f' % \
            (outFile, self.vmrWV)

          if mol == 'O2': outFile = '%s_vmrO2%06.0f' % \
            (outFile, self.vmrO2)

          if os.path.exists(outFile):
            print('WARNING: Overwriting %s' % outFile)

          outFP = open(outFile, 'w')
          outFP.write('$ %s\n' % self.headerOD)
          for rec in finalRecs: outFP.write('%s\n' % rec)
          outFP.write('%%%%')
          outFP.close()
        # end pressure level loop
      # end temperature loop
    # end band loop

    return True
  # end lblT5()

  def calcABSCO(self, mol):
    """
    Run LBLRTM for each TAPE5 made in lblT5(). Extract spectrum and 
    pressure layers from run results. Then to save time and space, 
    we will just calculate ABSCOs directory instead of transporting 
    the data for their computation into another method

    This can be run in parallel for each molecule, but that means 
    each molecule should have its own configuration file (each of 
    which should specify a different lbl_run_dir)
    """

    # find TAPE3s and their associated TAPE5s (for every TAPE3, 
    # there is a TAPE5 for every pressure and every temperature)
    molT3 = []
    for wn1, wn2 in zip(self.bands['wn1'], self.bands['wn2']):
      searchStr = '%s/%s/TAPE3_%s_%05d-%05d' % (\
        self.dirT3, mol, mol, wn1, wn2)

      # there should be only 1 match
      pathT3 = sorted(glob.glob(searchStr))
      if len(pathT3) == 1: molT3.append(pathT3[0])
    # end band loop

    if len(molT3) != len(self.bands['wn1']):
        print('Did not find TAPE3 files for all bands, ' + \
          'found %d files for %d bands' % \
          (len(molT3), len(self.bands['wn1'])))

    if len(molT3) == 0:
      print('Found no TAPE3 files for %s' % mol)
      return False
    # endif nT3

    # set up working subdirectory
    os.chdir(self.runDirLBL)

    # aliases for symlinks; should correspond to lblFiles
    # these are identical for all LBL runs for this task
    targets = ['lblrtm', 'xs', 'FSCDXS']
    sources = [self.pathLBL, self.pathXSDB, self.pathListXS]
    makeSymLinks(sources, targets)

    # outList is going to be an nBand-element list of dictionaries 
    # that contain nLay x nT x nWN arrays of ABSCOs ('ABSCO' field)
    # and nLay x nT arrays of layer pressures ('layerP' field)
    # the idea with layerP is "this layer pressure is associated with
    # the layer bounded by this lower level pressure and temperature"
    # and we will only need to store it once since it remains constant
    # over all bands
    outList = []
    for iBand, t3 in enumerate(molT3):
      # should be one TAPE3 per band, and we are assuming they're 
      # sorted by band
      base = os.path.basename(t3)
      band = base.split('_')[-1]

      bandDict = {}

      # wavenumber array of spectrum will remain the same for a 
      # given band
      bandWN = None

      # chose the first 2 levels for testing because that's when the
      # number of temperatures started to change
      pLevs = self.pLev[:3] if self.debug else list(self.pLev)

      # initialize output for given pressure
      # (which are dim [nT x nWN] and [nT])
      # the dictionary fields will be level pressures -- since each 
      # P has a different number of corresponding T values, we 
      # cannot simply make an nP x nT x nWN array
      # "pLayP": layer pressure (P) associated with level pressure (p)
      # "TLayP": layer temperature (T) associated with p
      levP = []
      pABSCO, pLayP, pLayT = {}, {}, {}
      for iP, pLev in enumerate(pLevs):
        # skip the TOA level because (from Karen):
        # LBL is calculating the altitude of the observer level and 
        # the profiles slightly differently, and coming up with an 
        # observer height above the top of the profile
        if iP == self.nP-1: continue

        # this is a bit different than self.pLev -- it represents 
        # what levels we processed rather than what we expected to use
        # skipping the last level is the difference here
        # but we count the first level because it is the lower bound
        # of the first layer
        levP.append(pLev)

        # do not expect anything for surface level since we skip it
        # in lblT5() (we need 2 levels for 1 layer)
        if iP == 0: continue

        tempABSCO, tempLayP, tempLayT = [], [], []
        for iT, tLev in enumerate(self.tLev[iP]):

          # find LBL TAPE5s corresponding to band
          searchStr = '%s/%s/%s/TAPE5_%s_P%09.4f_T%05.1fK_%s' % \
            (self.runDirLBL, self.dirT5, mol, mol, pLev, tLev, band)

          if mol in self.molH2O:
            searchStr += '_vmrWV%06.0f' % (self.vmrWV)

          if mol == 'O2': searchStr += '_vmrO2%06.0f' % (self.vmrO2)

          molT5 = sorted(glob.glob(searchStr))

          # should be one LBL TAPE5 per allowed T on a given P level
          if len(molT5) == 0:
            errMsg = 'Found no LBL TAPE5 for ' + \
              '%s, P=%-9.4f mbar, T=%-5.1f K' % \
              (os.path.basename(t3), pLev, tLev)
            print(errMsg)

            # still should append to the lists associated with each
            # temperature
            tempLayP.append(np.nan)
            tempLayT.append(np.nan)
            tempABSCO.append([np.nan])
            continue
          else:
            t5 = molT5[0]
          # endif t5

          # setup the LBL run
          if os.path.islink('TAPE3'): os.unlink('TAPE3')
          os.symlink(t3, 'TAPE3')

          base = os.path.basename(t5)
          print(base)
          if os.path.islink('TAPE5'): os.unlink('TAPE5')
          os.symlink(t5, 'TAPE5')

          # files that should be generated and have required info
          t7 = 'TAPE7'
          odFile = 'ODint_001'

          # first delete them if they already exist
          for lblFile in [odFile, t7]:
            if os.path.exists(lblFile): os.remove(lblFile)
          # end lblFile loop

          # run the model
          ext = base.replace('TAPE5_', '')
          status = sub.call(['./lblrtm'])

          # if LBL does not finish, no OD file is generated
          if not os.path.exists(odFile):
            tempLayP.append(np.nan)
            tempABSCO.append([np.nan])
            continue
          # endif odFile check

          # grab necessary parameters from TAPE7
          # store pressure layer calculated in LBLATM
          if self.doXS[mol][iBand]:
            t7Dict = RC.readXS('TAPE7', mol)
          else:
            t7Dict = RC.readTAPE7('TAPE7')
          # endif XS

          tempLayP.append(t7Dict['p_lay'][0])
          tempLayT.append(t7Dict['T_lay'][0])
          molDen = t7Dict['densities'][self.iMol][0]

          # extract the spectrum
          wnFine, odFine = lblTools.readOD(odFile, double=True)
          abscoFine = odFine / molDen

          # calculate absorption coefficients then degrade the 
          # spectrum. convolve ABSCO with weighting associated with 
          # kernel then resample; most of this taken directly from 
          # IDL code resample_grid.pro
          kernel = self.degradeKern[iBand]
          nDegrade = kernel.size - 1
          coarseRes = np.arange(0, wnFine.size, nDegrade)
          abscoCoarse = np.convolve(abscoFine, kernel)
          abscoCoarse[0], abscoCoarse[-1] = \
            abscoFine[0], abscoFine[-1]
          tempABSCO.append(abscoFine[coarseRes])

          # have to do 2nd conditional in case the first couple of 
          # runs did not extend the entire spectral range
          bandWN = wnFine[coarseRes]
        # end temperature loop

        # dimensions in tempABSCO are not consistent right now because
        # of the NaNs. we will fix this in arrABSCO()
        pABSCO[str(pLev)] = np.array(tempABSCO)
        pLayP[str(pLev)] = np.array(tempLayP)
        pLayT[str(pLev)] = np.array(tempLayT)
      # end pressure loop

      # save the important parameters in outList
      # number of frequencies will be the same for every LBL run
      bandDict['wavenum'] = np.array(bandWN)
      bandDict['nWN'] = bandWN.size
      bandDict['ABSCO'] = dict(pABSCO)
      bandDict['layerP'] = dict(pLayP)
      bandDict['levelP'] = np.array(levP)
      bandDict['layerT'] = dict(pLayT)
      outList.append(bandDict)
    # end t3 (band) loop

    os.chdir(self.gitDir)

    self.ABSCO = list(outList)

    return self
  # end calcABSCO()

  def arrABSCO(self):
    """
    Organize an array from the complicated dictionary returned by 
    calcABSCO(). This is where the (nWN x nP x nT) arrays 
    are generated.

    Not every LBL run produces spectra of the same size for a given 
    band, and we know that different temperatures apply depending on 
    the pressure, so this is where we reconcile differences in 
    dimensions
    """

    # construct the array, filling in NaNs with pABSCO
    wnAll = []
    for iBand, bandABSCO in enumerate(self.ABSCO):
      # each ABSCO dictionary has a different key for each P run
      pKeys = bandABSCO['ABSCO'].keys()
      numP = len(pKeys)

      nBandWN = int(bandABSCO['nWN'])
      bandArr = np.ones((numP, self.nT, nBandWN)) * np.nan
      wnAll.append(bandABSCO['wavenum'])

      # only need to do this once because it's the same for all bands
      if iBand == 0:
        layPArr = np.ones((numP, self.nT)) * np.nan
        layTArr = np.ones((numP, self.nT)) * np.nan
        levelT = np.ones((numP+1, self.nT)) * np.nan
        iMatchT = np.where(np.in1d(self.allT, self.tLev[0]))[0]
        if iMatchT.size != 0: levelT[0, iMatchT] = self.allT[iMatchT]
      # endif band

      for iP, pLev in enumerate(pKeys):
        # we skipped over the surface level P in the LBL runs, so 
        # indexing has to be changed accordingly
        offsetP = iP + 1

        # our T array spans from 180-320, but not every temperature 
        # is used. whatever range is used, it is in increments of 10 K
        # so let's find the indices of the T values that are used
        # for this pressure that correspond to the allT array
        iMatchT = np.where(np.in1d(self.allT, self.tLev[offsetP]))[0]

        # NaN spectrum for this P over all T if no match
        if iMatchT.size == 0: continue

        # layer P and level T population
        if iBand == 0:
          layPArr[iP, iMatchT] = bandABSCO['layerP'][pLev]
          layTArr[iP, iMatchT] = bandABSCO['layerT'][pLev]
          levelT[offsetP, iMatchT] = self.allT[iMatchT]
        # endif band 0

        # now build the spectrum array
        inArr = bandABSCO['ABSCO'][pLev]
        for iT in range(len(inArr)):
          # FILL IN EMPTY SPECTRA
          # here, we assume that if the LBL run did not extend the 
          # entire spectral range, then it stopped early and we 
          # append fill values to the end of the array
          nWN = len(inArr[iT])
          if nWN < nBandWN:
            nMiss = nBandWN - nWN
            inArr[iT] = np.hstack( \
              (inArr[iT], np.repeat(np.nan, nMiss) ))
          # endif nWN
          bandArr[iP, iMatchT[iT], :] = np.array(inArr[iT])
        # end T loop
      # end pLev loop

      arrABSCO = np.array(bandArr) if iBand == 0 else \
        np.append(arrABSCO, bandArr, axis=2)

    # end band loop

    # stack the frequencies and keep track of the start/end indices
    # for each band
    iStart, iEnd = [], []
    start, end = 0, 0
    for iBand, band in enumerate(wnAll):
      if iBand == 0: 
        wnStack = list(band)
      else:
        start = int(end)
        wnStack += list(band)
      # endif first band

      end += len(band)
      iStart.append(start)
      iEnd.append(end-1)
    # end band loop

    self.freq = np.array(wnStack)
    self.indBands = np.array([iStart, iEnd]).T

    # replace the ABSCO dictionary with the array we'll use in output
    # ABSCO axes need to be (nfreq, ntemp, npress)
    outAx = (2, 1, 0)
    self.ABSCO = np.transpose(np.array(arrABSCO), axes=outAx)
    self.layerP = np.array(layPArr)
    self.layerT = np.array(layTArr)
    self.levelT = np.array(levelT)
    self.levelP = np.array(bandABSCO['levelP'])

    return self
  # end arrABSCO()

  def makeNC(self, mol):
    """
    Generate netCDF that conforms to ESDS-RFC-028v1.1.pdf convention
    """

    import netCDF4 as nc

    outNC = '%s/%s_%05d-%05d_v%s_%s.nc' % \
      (self.outDir, mol, self.bands['wn1'][0], \
       self.bands['wn2'][-1], self.version, self.runDesc)
    print('Building %s' % outNC)
    if os.path.exists(outNC): print('WARNING: overwriting %s' % outNC)
    outFP = nc.Dataset(outNC, 'w')
    outFP.set_fill_on()

    pwvStr = 'H2O VMR, '  if mol in self.molH2O else ''
    outFP.description = 'Absorption coefficients for %s ' % mol
    outFP.description += 'as a function of pressure, temperature, '
    outFP.description += '%swavenumber, and band' % pwvStr

    # metadata attribute -- line parameter source
    outFP.source = str(self.dataSource[mol])

    # extract dimensions from data
    inDims = self.ABSCO.shape

    if mol in self.molH2O:
      if mol == 'O2':
        nFreq, nTemp, nLay, nVMR, nVMR = inDims
      else:
        nFreq, nTemp, nLay, nVMR = inDims
      # endif O2
    else:
      nFreq, nTemp, nLay = inDims
    # endif h2o

    nLev = nLay + 1

    # Number of bands actually calculated, whereas self.nBands is the
    # number of bands configured
    numProcBands = self.indBands.shape[0]

    outDimNames = \
      ['nfreq', 'nlev', 'nlay', 'ntemp', 'nranges', 'nranges_lims']
    outDimVals = [nFreq, nLev, nLay, nTemp, numProcBands, 2]

    if mol in self.molH2O:
      outDimNames.append('nvmr')
      outDimVals.append(nVMR)
      abscoDim = ('nfreq', 'ntemp', 'nlay', 'nvmr')
      if mol == 'O2':
        abscoDim = ('nfreq', 'ntemp', 'nlay', 'nvmr', 'nvmr')
    else:
      abscoDim = ('nfreq', 'ntemp', 'nlay')
    # endif WV mol

    for name, val in zip(outDimNames, outDimVals):
      outFP.createDimension(name, val)

    # now onto the variables
    outVar = outFP.createVariable('P_level', float, \
      ('nlev'), zlib=True, complevel=self.compress, fill_value=np.nan)
    outVar[:] = self.levelP
    outVar.units = 'mbar'
    outVar.long_name = 'Pressure Levels'
    outVar.valid_range = (0, 1050)
    outVar.description = 'User-provided layer boundary pressures'

    outVar = outFP.createVariable('P_layer', float, \
      ('nlay', 'ntemp'), zlib=True, complevel=self.compress, \
      fill_value=np.nan)
    outVar[:] = np.ma.array(self.layerP, mask=np.isnan(self.layerP))
    outVar.units = 'mbar'
    outVar.long_name = 'Layer Pressures'
    outVar.valid_range = (0, 1050)
    outVar.description = 'LBLRTM-calculated layer pressures'

    outVar = outFP.createVariable('T_layer', float, \
      ('nlay', 'ntemp'), zlib=True, complevel=self.compress, \
      fill_value=np.nan)
    outVar[:] = np.ma.array(self.layerT, mask=np.isnan(self.layerT))
    outVar.units = 'K'
    outVar.long_name = 'Layer temperatures'
    outVar.valid_range = (180, 320)
    outVar.description = 'LBLRTM-calculated layer temperatures'

    # Use minimum of the frequency dimension size versus the file 
    # configured chunk size
    # All other dimension are small enough to chunk outright
    chunksizes = list(inDims)
    chunksizes[0] = min(inDims[0], self.freq_chunk)

    outVar = outFP.createVariable('Cross_Section', float, abscoDim, \
      zlib=True, complevel=self.compress, chunksizes=chunksizes, \
      fill_value=np.nan)
    outVar[:] = np.ma.array(self.ABSCO, mask=np.isnan(self.ABSCO))
    outVar.units = 'cm2/molecule'
    outVar.long_name = 'Absorption Coefficients'
    outVar.valid_range = (0, 1e-20)
    outVar.description = 'Absorption coefficients (k) calculated ' + \
      'from LBLRTM optical depths and layer amounts'

    outVar = outFP.createVariable('Spectral_Grid', float, \
      ('nfreq'), zlib=True, complevel=self.compress, \
      fill_value=np.nan)
    outVar[:] = self.freq
    outVar.units = 'cm-1'
    outVar.long_name = 'Spectral Points'
    outVar.valid_range = (0, 50000)
    outVar.description = 'Spectral points corresponding to ' + \
      'ABSCOs in a single layer for a single temperature and in ' + \
      'a given spectral range'

    outVar = outFP.createVariable('T_level', float, \
      ('nlev', 'ntemp'), zlib=True, complevel=self.compress, \
      fill_value=np.nan)
    outVar[:] = np.ma.array(self.levelT, mask=np.isnan(self.levelT))
    outVar.units = 'K'
    outVar.long_name = 'Temperature Levels'
    outVar.valid_range = (180, 320)
    outVar.description = 'Applicable temperatures associated ' + \
      'with each layer boundary pressure'

    outVar = outFP.createVariable('Extent_Ranges', float, \
      ('nranges', 'nranges_lims'), zlib=True, \
      complevel=self.compress, fill_value=np.nan)
    # Include only the ranges actually computed
    for bIdx in range(numProcBands):
        iStart = self.indBands[bIdx, 0]
        iEnd = self.indBands[bIdx, 1]
        outVar[bIdx, :] = [self.freq[iStart], self.freq[iEnd]]
    outVar.units = 'cm-1'
    outVar.long_name = 'Spectral Ranges'
    outVar.valid_range = (0, 50000)
    outVar.description = 'Starting and ending spectral points ' + \
      'for each input channel'

    outVar = outFP.createVariable('Extent_Indices', int, \
      ('nranges', 'nranges_lims'), zlib=True, \
      complevel=self.compress, fill_value=sys.maxsize)
    outVar[:] = self.indBands[:]
    outVar.units = 'N/A'
    outVar.long_name = 'Spectral Ranges Reference Indices'
    outVar.valid_range = (0, sys.maxsize)
    outVar.description = 'Pairs of indices defining the start ' + \
      'and end indices of the Cross_Section frequency dimension ' + \
      'for non-continuous calculation regions'

    if mol in self.molH2O:
      outVar = outFP.createVariable('H2O_VMR', float, \
        ('nvmr'), zlib=True, complevel=self.compress, \
        fill_value=np.nan)
      outVar[:] = self.vmrArrWV
      outVar.units = 'ppmv'
      outVar.long_name = 'Water Vapor Mixing Ratio'
      outVar.valid_range = (0, 50000)
      outVar.description = 'Water vapor amount that influences ' + \
        'the continua of [%s] molecules' % (' '.join(self.molH2O))
    # end if WV mol

    if mol == 'O2':
      outVar = outFP.createVariable('O2_VMR', float, \
        ('nvmr'), zlib=True, complevel=self.compress, \
        fill_value=np.nan)
      outVar[:] = self.vmrArrO2
      outVar.units = 'ppmv'
      outVar.long_name = 'Oxygen Mixing Ratio'
      outVar.valid_range = (0, 250000)
      outVar.description = 'Two fixed amounts of O2 [ppmv]. ' + \
        'Used with O2 runs only'
    # end if WV mol

    outFP.close()

    print('Wrote %s' % outNC)
  # end makeNC()
# end makeABSCO()

def combineVMR(inList):
  """
  Combine attributes (ABSCO, vmrWV, and vmrO2) from multiple makeABSCO 
  objects generated with different water vapor or O2 VMRs

  Input
    inList -- list of makeABSCO objects with ABSCO arrays (with same
      dimensions!) and single-value vmr attributes. this should 
      only be a 2-element list, but the function is flexible enough 
      to handle more

  Output
    outObj -- modified makeABSCO object with new ABSCO and vmrWV 
      arrays
  """

  abscoList, vmrList = [], []
  for obj in inList:
    abscoList.append(obj.ABSCO)
  # end object loop

  # the rest of the makeABSCO objects should be the same, so just 
  # replace vmrWV and ABSCO; make the VMR dimension the last one 
  # in the outAxes and convert vmr to ppmv
  outAxes = (1, 2, 3, 0)

  # another dimension if we're working with O2
  if len(np.array(abscoList).shape) == 5: outAxes = (1, 2, 3, 4, 0)

  outObj = inList[0]

  outObj.ABSCO = np.transpose(np.array(abscoList), axes=outAxes)

  return outObj
# end combineVMR()

