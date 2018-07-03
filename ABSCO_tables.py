#!/usr/bin/env python

# standard Python libraries
# for Python 3 compatibility
from __future__ import print_function

import os, sys, glob, argparse
import numpy as np
import subprocess as sub

# path to GIT common submodules (not a Python standard lib)
sys.path.append('common')
import utils
import RC_utils as RC
import lblTools

class configure():
  def __init__(self, inFile):
    """
    Parse the input .ini file (inFile) and return as an object for 
    use in makeABSCO class.  Also do some error checking

    Inputs
      inFile -- string, full path to .ini file that specifies paths 
        and filenames for...
    """

    # all allowed molecule names
    allowed = ['H2O', 'CO2', 'O3', 'N2O', 'CO', 'CH4', 'O2', \
      'NO', 'SO2', 'NO2', 'NH3', 'HNO3', 'OCS', 'CH2O', 'N2', \
      'HCN', 'C2H2', 'HCOOH', 'C2H4', 'CH3OH', 'CCL4', 'CF4', \
      'F11', 'F12', 'F22', 'ISOPRENE', 'PAN', 'HDO', 'BRO', 'O2-O2']

    # standard library, but name depends on Python version
    if sys.version_info.major < 3:
      import ConfigParser
    else:
      import configparser as ConfigParser
    # endif Python version

    errMsg = 'Missing field in %s' % inFile

    cParse = ConfigParser.ConfigParser()
    cParse.read(inFile)
    cpSections = cParse.sections()

    # loop over each field (of all sections) and keep the field and 
    # associated value in returned object (self)
    for iCPS, cps in enumerate(cpSections):
      cItems = cParse.items(cps)
      if cps == 'channels':
        # make channels dictionary with spectral metadata
        # now only allowing 1 channel
        channels = {}
        for cItem in cItems:
          split = cItem[1].split()
          if len(split) > 1: 
            print('More than 1 channel specified, ', end='')
            print('only using the first one.')
          elif len(split) == 0:
            # CONSIDER DEFAULTS
            chanErrMsg = 'No bands specified in %s, returning' % \
              inFile
            sys.exit(chanErrMsg)
          # endif split

          channels[cItem[0]] = float(split[0])
        # end item loop

        # these keys are required
        keys = list(channels.keys())
        for req in ['wn1', 'wn2', 'res']:
          if req not in keys:
            print(errMsg)
            sys.exit('Could not find %s, returning' % req)
          # endif required
        # end required loop

        setattr(self, 'channels', channels)
      elif cps == 'molecules':
        # molecules should be separated by a space
        molecules = []

        molNames = cItems[0][1].upper()
        split = molNames.split()

        # make sure the provided species names are allowed
        # if not, do not include
        for mol in split:
          if mol not in allowed:
            print('Removing %s because it is not allowed' % mol)
            split.remove(mol)
          # endif mol
        # end mol loop

        # what molecules are active to the wn range specified?
        molActive = self.findActiveMol()
        if len(split) == 0:
          print('No molecules specified, ', end='')
          print('finding active molecules in %.2f-%.2f cm-1 range' % \
            (self.channels['wn1'], self.channels['wn2']))
          molecules = list(molActive)
        else:
          # did the user neglect any active molecules?
          molMissed = []
          for active in molActive:
            if active not in split: molMissed.append(active)
          # end active loop

          if len(molMissed):
            prompt = 'The following molecules are active between ' + \
              '%.2f and %.2f cm-1 ' % \
              (self.channels['wn1'], self.channels['wn2']) + \
              'and were not included in the configuration file: ' + \
              '%s, proceed (y/n)? ' % molMissed
            status = input(prompt)
            status = status.upper()

            if status == 'N': sys.exit('Exited without proceeding')
          # endif molMissed

          # did the user include and non-active molecules?
          molExtra = []
          for mol in split:
            if mol not in molActive: molExtra.append(mol)
          # end mol loop

          if len(molExtra):
            prompt = 'The following molecules are not active ' + \
              'between %.2f and %.2f cm-1 ' % \
              (self.channels['wn1'], self.channels['wn2']) + \
              'and were included in the configuration file: ' + \
              '%s, proceed (y/n)? ' % molExtra
            status = input(prompt)
            status = status.upper()

            if status == 'N': sys.exit('Exited without proceeding')
          # endif molExtra

          molecules += split
        # endif no mol

        setattr(self, 'molnames', molecules)
      else:
        for cItem in cItems: setattr(self, cItem[0], cItem[1])
      # endif cps
    # end sections loop

    # in the makeABSCO() class, we expect certain attributes
    # let's make sure they exist in the config file
    reqAtt = ['pfile', 'ptfile', 'channels', 'molnames', 'scale', \
      'lnfl_run_dir', 'lnfl_path', 'tape1_path', 'tape3_dir', \
      'extra_params', 'tape5_dir', 'lbl_path', 'xs_path', 'fscdxs', \
      'lbl_run_dir', 'od_dir', 'absco_dir']

    # loop over all required attributes and do a check
    for req in reqAtt:
      if req not in dir(self):
        print(errMsg)
        sys.exit('Could not find %s attribute, returning' % req)
      # endif req
    # end req loop

    # let's pack all of the files into a single list
    self.paths = [self.pfile, self.ptfile, self.extra_params, \
      self.lnfl_path, self.lbl_path, self.xs_path, self.fscdxs]
    self.outDirs = [self.lnfl_run_dir, self.lbl_run_dir, \
      self.tape3_dir, self.tape5_dir, self.od_dir, self.absco_dir]
  # end constructor

  def findActiveMol(self):
    """
    If the user only specifies a spectral range and no valid 
    molecules, try to determine the molecules to be processed based 
    on the input spectral range
    """

    # "active" spectral regions for each allowed molecule
    # these will only be used if the user does not specify a molecule
    regions = {}
    regions['H2O'] = [[100, 25000]]
    regions['CO2'] = [[500, 12785]]
    regions['O3'] = [[100, 4000], [8500, 24665], [27370, 54000]]
    regions['N2O'] = [[550, 7800]]
    regions['CO'] = [[2000, 2250], [4150, 4350], [6200, 6500], \
      [8200, 8465]]
    regions['CH4'] = [[100, 9100]]
    regions['O2'] = [[0, 150], [1350, 1850], [6200, 6500], \
      [7500, 8500], [9100, 11000], [12990, 13224], [15000, 29870], \
      [36000, 50000]]
    regions['NO'] = [[1750, 2000]]
    regions['SO2'] = [[450, 600], [1050, 1450], [2440, 2450], \
      [3950, 4150], [23995, 43985]]
    regions['NO2'] = [[1550, 1650], [2850, 2950], [15000, 42002]]
    regions['NH3'] = [[750, 1200], [1450, 1800], [3200, 3600]]
    regions['HNO3'] = [[400, 950], [1050, 1450], [1600, 1770]]
    regions['OCS'] = [[500, 2200]]
    regions['CH2O'] = [[25919, 33300]]
    regions['N2'] = [[0, 350], [2000, 2900], [4300, 4950]]
    regions['HCN'] = [[0, 100], [600, 800], [1300, 1500], \
      [3200, 3500]]
    regions['C2H2'] = [[650, 800], [3100, 3400]]
    regions['HCOOH'] = [[950, 1250], [1700, 1890]]
    regions['C2H4'] = [[600, 1175], [1350, 1550], [2900, 3243]]
    regions['CH3OH'] = [[950, 1408], [2600, 3250]]
    regions['CCL4'] = [[740, 820]]
    regions['CF4'] = [[1250, 1300]]
    regions['F11'] = [[830, 860], [1060, 1110]]
    regions['F12'] = [[860, 940], [1080, 1180]]
    regions['F22'] = [[780, 840], [1080, 1150], [1290, 1335]]
    regions['ISOPRENE'] = [[850, 1100], [2800, 3200]]
    regions['PAN'] = [[560, 1400], [1650, 1900]]
    regions['HDO'] = [[1100, 1800], [2500, 3000], [3300, 4300], \
      [4800, 5400], [7000, 7500]]
    regions['BRO'] = [[25927, 34919]]
    regions['O2-O2'] = [[16644, 29785]]

    # find active species inside specified band
    wn1, wn2 = self.channels['wn1'], self.channels['wn2']
    activeMol = []
    for key in regions.keys():
      active = False
      for band in regions[key]:
        if (wn1 >= band[0]) & (wn2 <= band[1]): active = True
      # end band loop
      if active: activeMol.append(key)
    # end key loop

    return activeMol
  # end findActiveMol()
# end configure()

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
      inObj -- configure instance
      molecule -- str, HITRAN molecule name

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

    # gather necessary params from input files
    # [temperature] = K, [pressure] = mbar
    # let's go from surface to TOA
    inP = np.loadtxt(inObj.pfile)

    # is pressure ascending or descending? force descending
    pDiff = np.diff(inP)
    if (pDiff > 0).all():
      inP = inP[::-1]
    elif (pDiff < 0).all():
      pass
    else:
      sys.exit('Please provide monotonic pressures')
    # endif ascend

    allT = np.loadtxt(inObj.ptfile, usecols=(3), unpack=True)
    uniqT = np.unique(allT)

    # set class attributes
    # state, etc. atts
    self.headerOD = inObj.header
    self.cntnmScale = float(inObj.scale)
    self.pLev = np.array(inP)
    self.pressures = np.array(inP)
    self.tLev = np.array(uniqT)
    self.bands = dict(inObj.channels)
    self.nBands = len(inObj.channels['res'])
    self.molNames = list(inObj.molnames)

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

    # all HITRAN molecule names
    """
    # might be useful later...
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
    """

    # some constant LBLRTM TAPE5 records 
    # (see lblrtm_instructions.html for help with each record)
    # record1.2: HI=9: central line contribution omitted
    # CN=6: continuum scale factor for given molecules used
    # OD=1, MG=1: optical depth computation, layer-by-layer
    record12 = ' HI=9 F4=0 CN=6 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 ' + \
      'AM=1 MG=1 LA=0 OD=1 XS=0'

    # multiplicative continuum factors (because CN=6 in record12)
    scales = np.repeat(0.0, 7)

    # records required with IATM=1 (provide some doc on this rec)
    # US Standard atmosphere, path type 2 (slant from H1 to H2), 2 
    # pressure levels, no zero-filling, full printout, 7 molecules,
    record31 = '%5d%5d%5d%5d%5d%5d' % (6, 2, -2, 0, 0, 7)

    for mol in self.molNames:
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

        for itLev, tLev in enumerate(self.tLev):
          for iBand in range(self.nBands):
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
            recs = [record12, record12a, record13, \
              record31, record32, record33b]

            if os.path.exists(outFile):
              print('WARNING: Overwriting %s' % outFile)

            outFP = open(outFile, 'w')
            outFP.write('$ %s\n' % self.headerOD)
            for rec in recs: outFP.write('%s\n' % rec)
            outFP.write('%%%%')
            outFP.close()
          # end band loop
        # end temperature loop
      # end pressure loop
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
  ini = configure(iniFile)

  # ABSCO object instantiation
  absco = makeABSCO(ini)

  if args.run_lnfl: 
    absco.lnflT5()
    absco.runLNFL()
  # end LNFL

  if args.run_lbl:
    absco.lblT5()
    #absco.runLBL()
  # end LBL

  # haven't tested this yet, but no reason it won't work...right?
  if args.end_to_end:
    sys.exit('No e2e yet')
    absco.lblT5()
    absco.runLBL()
  # endif e2e
# end main()

