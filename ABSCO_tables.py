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

class configSetup():
  def __init__(self, inFile):
    """
    Parse the input .ini file (inFile) and return as an object for 
    use in makeABSCO class.  Also do some error checking

    Inputs
      inFile -- string, full path to .ini file that specifies paths 
        and filenames for...
    """

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
        channels = {}
        for cItem in cItems: 
          channels[cItem[0]] = \
            np.array(cItem[1].split()).astype(float)
        # end item loop

        # these keys are required
        keys = list(channels.keys())
        for req in ['wn1', 'wn2', 'res']:
          if req not in keys:
            print(errMsg)
            sys.exit('Could not find %s, returning' % req)
          # endif required
        # end required loop

        # there should be an equal number of starting and ending 
        # wavenumbers and associated spectralresolutions
        for key in keys[1:]:
          if channels[key].size != channels[keys[0]].size:
            chanErrMsg = 'wn1, wn2, and res should have equal ' + \
              'number of elements, returning'
            print('Error in %s' % inFile)
            sys.exit(chanErrMsg)
          # endif channels
        # end key loop

        if channels[keys[0]].size == 0:
          # CONSIDER DEFAULTS
          chanErrMsg = 'No bands specified in %s, returning' % \
            inFile
          sys.exit(chanErrMsg)
        # endif zero

        setattr(self, 'channels', channels)
      elif cps == 'molecules':
        # molecules should be separated by a space
        # CONSIDER DEFAULTS; put list of defaults somewhere in here 
        # and don't do the upper(); stay case sensitive
        molecules = []

        # this is a problem for HOBr and CH3Br
        molNames = cItems[0][1].upper()
        split = molNames.split()
        if len(split) == 0:
          sys.exit('No molecules specified')
        else:
          molecules += split
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
    self.outDirs = [self.tape3_dir, self.tape5_dir, \
      self.lnfl_run_dir, self.lbl_run_dir, \
      self.od_dir, self.absco_dir]
  # end constructor
# end configSetup()

def makeSymLinks(sources, targets):
  """
  Loop over input files and make symbolic links for them
  """

  for source, target in zip(sources, targets):
    os.symlink(source, target)

  return None
# end makeSymLinks()

class makeABSCO():
  """
  - Build TAPE5s (LNFL and LBLRTM) for each molecule of interest
  - Build TAPE3 (binary line file) for specified bands
  - Run LBLRTM to generate ODInt files
  - Use the ODInt files to calculate absorption coefficients and store
    in HDF file

  Generate absorption coefficient tables (ABSCO as a function of 
  wavenumber, pressure, temperature, and band) for specified molecule
  """

  def __init__(self, inObj):

    """
    Inputs
      inObj -- configSetup instance
      molecule -- str, HITRAN molecule name

    Keywords
    """

    # make sure paths exist before proceeding
    #for path in inObj.paths: utils.file_check(path)

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
    self.headerOD = inObj.header
    self.cntnmScale = float(inObj.scale)
    self.pLev = np.array(inP)
    self.pressures = np.array(inP)
    self.tLev = np.array(uniqT)
    self.bands = dict(inObj.channels)
    self.nBands = len(inObj.channels['res'])
    self.molNames = list(inObj.molnames)
    self.pathLNFL = str(inObj.lnfl_path)
    self.runDirLNFL = str(inObj.lnfl_run_dir)
    self.dirExtras = str(inObj.extra_params)
    self.dirT3 = str(inObj.tape3_dir)
    self.runDirLBL = str(inObj.lbl_run_dir)
    self.dirT5 = str(inObj.tape5_dir)
    self.pathLBL = str(inObj.lbl_path)
    self.pathXSDB = str(inObj.xs_path)
    self.pathListXS = str(inObj.fscdxs)
    self.fineOD = str(inObj.od_dir)
    self.coarseOD = str(inObj.absco_dir)

    # all HITRAN molecule names
    # might be useful later...
    """
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
        wvn1 = self.bands['wn1'][iBand]
        wvn2 = self.bands['wn2'][iBand]
        record1 = 'TAPE5 for %s, %.f-%.f' % (mol, wvn1, wvn2)
        record2 = '%10.3f%10.3f' % (wvn1, wvn2)

        molInd = np.array(molIndInit)
        molInd[iMol] = '1'
        record3 = '%47s%4s%40s' % (''.join(list(molInd)), ' ', holInd)

        # write LNFL TAPE5 for this band and molecule
        recs = [record1, record2, record3]

        # making WN1 and WN2 ints just to keep "." out of name
        outFile = '%s/TAPE5_%s_%d-%d' % \
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
    molecule.
    """

    # link to extra broadening and speed dependence parameters
    # equivalent to `ln -s full_path`
    os.chdir(self.runDirLNFL)
    extras = glob.glob(self.dirExtras)
    slExtras = [os.path.basename(extra) for extra in extras]
    makeSymLinks(extras, slExtras)

    for mol in self.molNames:
      inDirT5 = '%s/%s' % (self.dirT5, mol)
      inT5 = sorted(glob.glob('%s/TAPE5_*' % inDirT5))

      for t5 in inT5:
        continue
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
    # print to TAPE7 (not actually needed for final product, probably
    # useful for debugging)
    record31 = '%5d%5d%5d%5d%5d%5d%5d' % (6, 2, -2, 0, 0, 7)

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
            outFile = '%s/TAPE5_%s_P%09.4fmb_T%05.1fK' % \
              (outDirT5, mol, pLev, tLev)

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

    This can be run in parallel for each molecule
    """

    if not os.path.exists(self.workDir): os.mkdir(self.workDir)
    for mol in self.molNames['lblxs']:
      # set up working subdirectory
      workSubDir = '%s/%s' % (self.workDir, mol)
      outDirOD = '%s/%s' % (workSubDir, self.fineOD)
      if not os.path.exists(workSubDir): os.mkdir(workSubDir)
      if not os.path.exists(outDirOD): os.mkdir(outDirOD)
      os.chdir(workSubDir)

      # aliases for symlinks; should correspond to lblFiles
      # these are identical for all LBL runs for this task
      targets = ['lblrtm', 'TAPE3', 'xs', 'FSCDXS']
      lblFiles = \
        [self.pathLBL, self.pathT3, self.pathXS, self.pathFSCDXS]
      for source, target in zip(lblFiles, targets):
        utils.file_check(source)

        # will crash if the link already exists
        if not os.path.islink(target): os.symlink(source, target)
      # end LBL file loop

      # link to the TAPE5s created with lblT5() and run LBL
      globStr = '%s/%s/TAPE5*' % (self.dirT5, mol)
      lblT5 = sorted(glob.glob(globStr))

      if len(lblT5) == 0:
        print('Found no match: %s' % globStr)
        return
      # endif T5

      for t5 in lblT5:
        base = os.path.basename(t5)
        print(base)
        if os.path.islink('TAPE5'): os.remove('TAPE5')
        os.symlink(t5, 'TAPE5')

        # grab extension for use in renaming the ODint LBL output file
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

        os.rename(\
          odStr, '%s/%s' % (self.fineOD, odStr.replace('001', ext)) )
      # end T5 loop
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
  ini = configSetup(iniFile)

  # instantiation
  absco = makeABSCO(ini)

  if args.lnfl_tape5: absco.lnflT5()
  if args.run_lnfl: absco.runLNFL()
  if args.lbl_tape5: absco.lblT5()
  if args.run_lbl: absco.runLBL()

  # haven't tested this yet, but no reason it won't work...right?
  if args.end_to_end:
    absco.lblT5()
    absco.runLBL()
  # endif e2e
# end main()

