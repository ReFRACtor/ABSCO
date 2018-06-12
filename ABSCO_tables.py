#!/usr/bin/env python

# standard Python libraries
# for Python 3 compatibility
from __future__ import print_function

import os, sys, glob, argparse
import numpy as np

# path to GIT common submodules (not a Python standard lib)
sys.path.append('common')
import utils
import RC_utils as RC

class makeABSCO():
  """
  - Build TAPE5s for each molecule of interest
  - Run LBLRTM to generate ODInt files
  - Use the ODInt files to calculate absorption coefficients and store
    in HDF file

  Generate absorption coefficient tables for molecules of interest
  """

  def __init__(self, inObj):

    """
    Inputs
      inObj -- configSetup instance

    Keywords
    """

    for temp in [ptnFile, tBarFile]: utils.file_check(temp)
    for temp in [dirT5, dirDen]:
      if not os.path.exists(temp): os.mkdir(temp)

    # gather necessary params from input files 
    # ([density] = molecules/cm2)
    pBar, tBar, wetAir, h2o, dryAir, o2 = np.loadtxt(\
      ptnFile, usecols=(4, 5, 6, 7, 18, 13), skiprows=5, unpack=True)
    pLev, tLev = np.loadtxt(tBarFile, usecols=(0, 3), unpack=True)

    # set class attributes
    self.ptnFile = ptnFile
    self.TMeanFile = tBarFile
    self.dirT5 = dirT5
    self.dirDensities = dirDen
    self.vmrFile = vmr
    self.pLay = pBar
    self.tLay = tBar
    self.nWetAir = wetAir
    self.nDryAir = dryAir
    self.nH2O = h2o
    self.nO2 = o2
    self.pLev = pLev
    self.tLev = tLev
    self.nLay = pBar.size
    self.nLev = pLev.size
    self.bands = channels
    self.molNames = molecules

    # channel provision: do all keys have the same number of elements
    if not (len(channels['names']) == len(channels['wn1']) == \
      len(channels['wn2']) == len(channels['res'])):
      sys.exit('Channels params do not have equal length')
    self.nBands = len(channels['names'])

    # molecule provision: are we using the same number of ELANOR and 
    # LBL XS names?
    if len(molecules['lblxs']) != len(molecules['elanorxs']):
      print(molecules['lblxs'], molecules['elanorxs'])
      sys.exit('Molecules arrays are not the same size')
    self.nMol = len(molecules['lblxs'])

    self.pathLBL = pathLBL
    self.pathT3 = pathT3
    self.pathXS = pathXS
    self.pathFSCDXS = pathFSCDXS
    self.workDir = dirWork
    self.fineOD = dirFine
    self.coarseOD = dirCoarse

    # i thought this would be used in the eventual binary files, but
    # that's not the case; maybe put it in the TAPE5 header?
    self.headerOD = odHeader

    # for cd'ing back into the cwd
    self.topDir = os.getcwd()
  # end constructor()

  def calcDensity(self):
    """
    Interpolate XS VMRs from LBLATM altitude grid to ABSCO grid
    Calculate densities for each molecule of interest
    Write densities to file for each species
    This will not need to be done everytime
    """
    from scipy import interpolate

    for inFile in [self.vmrFile, self.ptnFile]:
      utils.file_check(inFile)

    # ABSCO altitudes (bottom of layer)
    altABSCO = np.loadtxt(self.ptnFile, unpack=True, \
      usecols=[1], skiprows=5)

    # XS molecule names
    lblHeader = open(self.vmrFile).read().splitlines()[0].split(',')

    # find indices (column numbers) of XS molecules
    lblUseCol = []
    for mol in self.molNames['lblxs']:
      try:
        iXS = lblHeader.index(mol)
        lblUseCol.append(iXS)
      except:
        print('Could not find %s in %s' % (mol, self.vmrFile))
        continue
      # end exception
    # end molecule loop
    
    # loop over each XS and interpolate LBL VMR onto ABSCO altitude
    # grid, then calculate and store associated number densities
    for iXS, lblXS in enumerate(lblUseCol):
      altLBL, vmrLBL = np.loadtxt(self.vmrFile, unpack=True, \
        usecols=(0, lblXS), skiprows=1, delimiter=',')
      vmrInterp = interpolate.interp1d(altLBL, vmrLBL, kind='linear')
      vmrLBL = vmrInterp(altABSCO)
      denXS = vmrLBL * self.nDryAir

      outFile = '%s/%s_densities.txt' % \
        (self.dirDensities, self.molNames['lblxs'][iXS])
      outFP = open(outFile, 'w')
      for p, den in zip(self.pLay, denXS):
        outFP.write('%10.4f%10.3e\n' % (p, den))
      outFP.close()
    # end loop over LBL XS

    return True
  # end calcDensity()

  def makeTAPE5(self):
    """
    Make a TAPE5 for every temperature, band, and molecule
    """

    # some constant LBLRTM TAPE5 records
    # record12 copied straight from make_pan_absco_tape5s.pro
    record12 = ' HI=9 F4=0 CN=0 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 ' + \
      'AM=0 MG=1 LA=0 OD=1 XS=1    0    0'
    record21 = '    1    7   1.00000' # 1 layer, 7 molecules
    record22 = '    1        0' # 1 XS, convolved with P
    record222 = ' 1  1    1  1.000000    XS Values'

    for iMol, mol in enumerate(self.molNames['lblxs']):
      print('Writing %s TAPE5s' % mol)
      xsFile = '%s/%s_densities.txt' % (self.dirDensities, mol)
      utils.file_check(xsFile)
      pMol, denMol = np.loadtxt(xsFile, unpack=True)

      outDirT5 = '%s/%s' % (self.dirT5, mol)
      if not os.path.exists(outDirT5): os.mkdir(outDirT5)

      # TAPE5 records
      record221 = '%10s' % mol

      for iP, pLay in enumerate(self.pLay):
        # find number of temperatures at this pressure
        iTemp = np.where(self.pLev == (iP+1) )

        if iTemp[0].size == 0: continue

        broadener = self.nDryAir[iP] * 0.79 # N2 in molecules/cm2

        # TAPE5 records
        record212 = ''.join(['%10.3e' % 0] * 7)
        record212 += '%10.3e' % broadener
        record224 = '%15.7e' % (denMol[iP])
        for itLev, tLev in enumerate(self.tLev[iTemp]):
          for iBand, band in enumerate(self.bands['names']):
            outFile = '%s/TAPE5_%s_P%09.4fmb_T%05.1fK_%s' % \
              (outDirT5, mol, pLay, tLev, band)

            # TAPE5 records
            # record 1.3 is kinda long...
            record13 = '%10.3e%10.3e' % \
              (self.bands['wn1'][iBand], self.bands['wn2'][iBand])
            # concatenate (NOT append) 6 zeros in scientific notation
            record13 += ''.join(['%10.3e' % 0 for i in range(6)])
            record13 += '%4s%1d%5s%10.3e' % \
              (' ', 0, ' ', self.bands['res'][iBand])

            record211 = '%10.4f%10.4f' % (pLay, tLev)
            record223 = '%15.7e%10.4f' % (pLay, tLev)

            # write the TAPE5 for this set of params
            recs = [record12, record13, record21, record211, \
              record212, record22, record221, record222, \
              record223, record224]

            if os.path.exists(outFile):
              print('WARNING: Overwriting %s' % outFile)

            outFP = open(outFile, 'w')
            outFP.write('$ %s ABSCO %s\n' % \
              (mol, os.path.basename(outFile)))
            for rec in recs: outFP.write('%s\n' % rec)
            outFP.write('%%%%')
            outFP.close()
          # end band loop
        # end temperature loop
      # end pressure loop
    # end molecule loop

    return True
  # end makeTAPE5()

  def runLBL(self):
    """
    Run LBLRTM for each TAPE5 made in makeTAPE5

    This can be run in parallel for each molecule
    """

    # standard library
    import subprocess as sub

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

      # link to the TAPE5s created with makeTAPE5() and run LBL
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
        sub.call(['lblrtm'])
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

# end xsABSCO()

class configSetup():
  def __init__(self, inFile):
    """
    Parse the input .ini file (inFile) and return as a dictionary for 
    use in the rest of this module

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
        # CONSIDER DEFAULTS
        molecules = []
        molNames = cItems[0][1]
        split = molNames.split()
        if len(split) == 0:
          sys.exit('No molecules specified')
        else:
          molecules += split
        setattr(self, 'molnames', molecules)
      else:
        for cItem in cItems:
          if cItem[0] == 'tape5_dir':
            setattr(self, cItem[0], '%s/%s' % (os.getcwd(), cItem[1]))
          else:
            setattr(self, cItem[0], cItem[1])
          # endif t5
        # end item loop
      # endif cps
    # end sections loop

    # in the makeABSCO() class, we expect certain attributes
    # let's make sure they exist in the config file
    reqAtt = ['header', 'channels', 'molnames', \
      'tape5_dir', 'lbl_path', 'tape3_path', 'xs_path', 'fscdxs', \
      'lbl_run_dir', 'od_dir', 'absco_dir']

    # loop over all required attributes and do a check
    for req in reqAtt:
      if req not in dir(self):
        print(errMsg)
        sys.exit('Could not find %s attribute, returning' % req)
      # endif req
    # end req loop

  # end constructor
# end configSetup()

if __name__ == '__main__':

  parser = argparse.ArgumentParser(\
    description='Generate ABSCO tables for user-specified molecules.')
  parser.add_argument('--config_file', type=str, \
    default='ABSCO_config.ini', \
    help='Configuration file that contains the file and ' + \
    'directory names necessary for the makeABSCO class.')

  # argument switches (booleans, no assignments)
  parser.add_argument('-t5', '--make_tape5', action='store_true', \
    help='Make the TAPE5 files for each species, channel, ' + \
    'pressure, and temperature.')
  parser.add_argument('-lbl', '--run_lbl', action='store_true', \
    help='Run LBLRTM for each of the TAPE5s generated and saved ' + \
    'in makeTAPE5().')
  parser.add_argument('-e2e', '--end_to_end', action='store_true', \
    help='Runs the entire process from tape 5 generation to ' + \
    'post-processing (rather than entering all of the keywords ' + \
    'separately).')
  args = parser.parse_args()

  iniFile = args.config_file; utils.file_check(iniFile)
  ini = configSetup(iniFile)
  sys.exit()

  # instantiation; there are lots of keywords here so that we have 
  # some flexibility (i did not wanna force the user to have to use 
  # a configSetup object with the xsABSCO class)
  absco = makeABSCO(ini)

  if args.make_tape5: absco.makeTAPE5()
  if args.run_lbl: absco.runLBL()

  # haven't tested this yet, but no reason it won't work...right?
  if args.end_to_end:
    absco.makeTAPE5()
    absco.runLBL()
  # endif e2e
# end main()

