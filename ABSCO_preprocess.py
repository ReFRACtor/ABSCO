# standard Python libraries
# for Python 3 compatibility
from __future__ import print_function
import os, sys

# path to GIT common submodules
sys.path.append('common')
import utils

# miniconda-installed libs
import numpy as np
import pandas as pd

class configure():
  def __init__(self, inFile, molActiveCheck=True, prompt_user=True):
    """
    Parse the input .ini file (inFile) and return as an object for 
    use in makeABSCO class.  Also do some error checking and other 
    processes that are not directly related to running LNFL and LBLRTM

    Inputs
      inFile -- string, full path to .ini file that specifies paths 
        and filenames for...

    Keywords
      molActiveCheck -- boolean, run a test to determine for which 
        bands LNFL and LBLRTM should be run based on the user-input 
        ranges and whether the given molecules are active in the bands
    """
    # Allow turning of interactive prompts for batch processing
    self.prompt_user = prompt_user

    # extract everything from config (.ini) file (inFile)
    self.configFile = str(inFile)
    self.readConfig()
    self.checkAtts()

    # these intermediate subdirs should be underneath intdir, 
    # which can be an another file system
    self.lnfl_run_dir = os.path.join(self.intdir, self.lnfl_run_dir)
    self.lbl_run_dir = os.path.join(self.intdir, self.lbl_run_dir)
    self.tape3_dir = os.path.join(self.intdir, self.tape3_dir)
    self.outdir = os.path.join(self.intdir, self.outdir)

    # these guys should always be in the "source" directory -- the 
    # directory into which the ABSCO repo is cloned, which is assumed
    # to be the working dir
    gitDir = os.path.dirname(__file__)
    self.pfile = os.path.join(gitDir, self.pfile)
    self.ptfile = os.path.join(gitDir, self.ptfile)
    self.vmrfile = os.path.join(gitDir, self.vmrfile)
    self.xs_lines = os.path.join(gitDir, self.xs_lines)

    # let's pack all of the files into a single list
    self.paths = [self.pfile, self.ptfile, self.vmrfile, \
      self.extra_params, self.lnfl_path, \
      self.lbl_path, self.xs_path, self.fscdxs, self.xs_lines]
    self.outDirs = [self.lnfl_run_dir, self.lbl_run_dir, \
      self.tape3_dir, self.tape5_dir, self.outdir]

    # make sure paths exist before proceeding
    for path in self.paths: utils.file_check(path)

    # cross section molecules will have to be handled differently
    # from HITRAN molecules
    self.xsNames = ['CCL4', 'F11', 'F12', 'F22', 'ISOP', 'PAN', 'BRO']

    # these also have line parameters, but HITRAN recommends to use 
    # their XS coefficients in certain bands (see xsCheck() method).
    self.xsLines = ['CF4', 'SO2', 'NO2', 'HNO3']

    # these guys have continua that are affected by WV amount
    self.molH2O = ['CO2', 'N2', 'O2', 'H2O']

    # and these guys are neither HITRAN or XS molecules
    self.dunno = ['HDO']

    # O2 is a special case where we'll have to fix the VMR
    self.vmrO2 = np.array([0.19, 0.23]) * 1e6

    # read in pressures and do some quality control
    self.processP()

    # XS or line parameter processing?
    self.xsCheck()

    # verification of molecules and associated bands to process
    # this should be called to limit unnecessary runs of LNFL and LBL
    self.molProcess()

    # spectral degradation weight calculation
    self.degradeWeighting()

    # was the line parameter source HITRAN or AER Line Parameter Data?
    self.determineSources()

    self.calcRAM()
  # end constructor

  def readConfig(self):
    """
    Read configuration file and set attributes of configure object
    accordingly
    """

    # all allowed molecule names
    allowed = ['H2O', 'CO2', 'O3', 'N2O', 'CO', 'CH4', 'O2', \
      'NO', 'SO2', 'NO2', 'NH3', 'HNO3', 'OCS', 'H2CO', 'N2', \
      'HCN', 'C2H2', 'HCOOH', 'C2H4', 'CH3OH', 'CCL4', 'CF4', \
      'F11', 'F12', 'F22', 'ISOP', 'PAN', 'HDO', 'BRO']

    # standard library, but name depends on Python version
    if sys.version_info.major < 3:
      import ConfigParser
    else:
      import configparser as ConfigParser
    # endif Python version

    errMsg = 'Missing field in %s' % self.configFile

    cParse = ConfigParser.ConfigParser()
    cParse.read(self.configFile)
    cpSections = cParse.sections()

    # loop over each field (of all sections) and keep the field and 
    # associated value in returned object (self)
    for iCPS, cps in enumerate(cpSections):
      cItems = cParse.items(cps)
      if cps == 'channels':
        # make channels dictionary with spectral metadata
        channels = {}
        for cItem in cItems:
          split = cItem[1].split()

          if len(split) == 0:
            # CONSIDER DEFAULTS
            chanErrMsg = 'No bands specified in %s, returning' % self.configFile 
            sys.exit(chanErrMsg)
          # endif split

          # should only be one unit, but there can be many channels
          if cItem[0] == 'units':
            units = str(split[0]).lower()
          else:
            channels[cItem[0]] = np.array(split).astype(float)
          # endif units
        # end item loop

        # kernel width for degradation
        channels['kernwidth'] = channels['outres']/channels['lblres']

        # these keys are required
        keys = list(channels.keys())
        for req in ['wn1', 'wn2', 'lblres', 'outres', 'kernwidth']:
          if req not in keys:
            sys.exit('Could not find %s, returning' % req)
          # endif required
        # end required loop

        # is "units" defined?
        errMsg = 'Please define "units" field ("cm-1", "um", "nm")'
        if 'units' in locals():
          if units not in ['cm-1', 'um', 'nm']: sys.exit(errMsg)
          setattr(self, 'spectral_units', units)
        else:
          sys.exit(errMsg)
        # endif units

        # there should be an equal number of starting and ending 
        # wavenumbers and associated spectral resolutions
        # need to also check that each wn1 < wn2
        for key in keys[1:]:
          if channels[key].size != channels[keys[0]].size:
            chanErrMsg = 'wn1, wn2, and res should have equal ' + \
              'number of elements, returning'
            print('Error in %s' % self.configFile)
            sys.exit(chanErrMsg)
          # endif channels
        # end key loop

        # make sure degradation is an integer exponent of 2
        if not np.all(np.log2(channels['kernwidth']) % 1 == 0):
          errMsg = 'Please provide an outres/lblres ratio that ' + \
            'is an integer exponent of 2 (2, 4, 8, 16, etc.)'
          sys.exit(errMsg)
        # endif degrade

        # LBLRTM can only handle 2000 cm-1 chunks at a time, so break
        # up any bands that are larger than this
        # also make sure user provides correct WN order and convert 
        # to cm-1 if necessary
        for iBand in range(len(channels['wn1'])):
          wn1, wn2, res, deg, kern = \
            channels['wn1'][iBand], channels['wn2'][iBand], \
            channels['lblres'][iBand], channels['outres'][iBand], \
            channels['kernwidth'][iBand]

          if units in ['um', 'nm']:
            wnConvert = 1e4 if units == 'um' else 1e7
            channels['wn1'][iBand] = wnConvert / wn1
            channels['wn2'][iBand] = wnConvert / wn2
            #channels['lblres'][iBand] = wnConvert / res
            print('Converting band ', end='')
            print('%d (%.3f-%.3f) to cm-1 (%.3f-%.3f)' % \
              (iBand+1, wn1, wn2, channels['wn1'][iBand], \
               channels['wn2'][iBand]))
            wn1 = wnConvert / wn1
            wn2 = wnConvert / wn2
          # endif um

          if wn2 < wn1:
            print('Swapping user-provided WN1 and WN2 for ' + \
              'channel %d' % (iBand+1))
            channels['wn1'][iBand] = wn2
            channels['wn2'][iBand] = wn1
            wn1 = channels['wn1'][iBand]
            wn2 = channels['wn2'][iBand]
          # endif wn2/wn1 switch

          bandwidth = wn2-wn1
          if bandwidth > 2000:
            print('Breaking down channels so they do not ' + \
              'exceed 2000 cm-1 widths')
            subChanWN1, subChanWN2 = [], []
            tempWN1 = float(wn1)

            while (bandwidth/2000 >= 1.0):

              tempWN2 = tempWN1 + 2000
              if tempWN2 >= wn2: tempWN2 = float(wn2) 

              subChanWN1.append(tempWN1)
              subChanWN2.append(tempWN2)

              tempWN1 = float(tempWN2)
              tempWN2 = tempWN1 + 2000
              if tempWN2 >= wn2: tempWN2 = float(wn2)
              bandwidth = tempWN2-tempWN1
            # endif > 1

            subChanWN1.append(tempWN1)
            subChanWN2.append(tempWN2)

            # remove the entire width > 2000 cm-1 band
            for key in ['wn1', 'wn2', 'lblres', 'outres', 'kernwidth']:
              channels[key] = np.delete(channels[key], iBand)

            # reassign channel limits
            channels['wn1'] = \
              np.insert(channels['wn1'], iBand, subChanWN1)
            channels['wn2'] = \
              np.insert(channels['wn2'], iBand, subChanWN2)

            # now apply resolution and degradation from entire 
            # (width > 2000 cm-1) channel to the subchannels
            channels['lblres'] = np.insert(channels['lblres'], \
              iBand, np.repeat(res, len(subChanWN1)))
            channels['outres'] = np.insert(channels['outres'], \
              iBand, np.repeat(deg, len(subChanWN1)))
            channels['kernwidth'] = np.insert(channels['kernwidth'], \
              iBand, np.repeat(kern, len(subChanWN1)))
          # endif > 2000
        # end band loop

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

        setattr(self, 'molnames', split)
      elif cps == 'vmr':
        for cItem in cItems:
          vmrArr = np.array(cItem[1].split()).astype(float)
          if vmrArr.size != 2:
            sys.exit('Please provide 2 values for %s' % cItem[0])
          setattr(self, cItem[0], vmrArr)
        # end cItem loop
      else:
        for cItem in cItems:
          param = cItem[0]
          val = cItem[1]
          if len(val.split()) > 1 and param != 'header':
            errMsg = '%s should only have 1 value assigned to it' % \
              param
            print(errMsg)
            sys.exit(1)
          # endif header

          if param == 'intdir':
            val = os.getcwd() if val == '.' else str(val)
            if not os.path.exists(val): os.makedirs(val)
          # endif 

          setattr(self, param, val)
        # end item loop
      # endif cps
    # end sections loop

    return self

  # end readConfig()

  def checkAtts(self):
    """
    Make sure all of the attributes in the configure object that are 
    necessary for ABSCO_tables.py exist (i.e., that the user did not 
    delete or forget a field in the .ini file)
    """

    # in the makeABSCO() class, we expect certain attributes
    # let's make sure they exist in the config file
    reqAtt = ['pfile', 'ptfile', 'vmrfile', 'channels', 'intdir', \
      'molnames', 'scale', 'lnfl_run_dir', 'lnfl_path', \
      'tape1_path', 'tape3_dir', 'extra_params', 'tape5_dir', \
      'lbl_path', 'xs_path', 'fscdxs', 'lbl_run_dir', 'outdir']

    # loop over all required attributes and do a check
    for req in reqAtt:
      if req not in dir(self):
        errMsg = 'Could not find %s attribute in %s, returning' % \
          (req, self.configFile)
        sys.exit(errMsg)
      # endif req
    # end req loop
  # end checkAtts()

  def processP(self):
    """
    Read in user and standard atmosphere pressures (associated
    with VMRs) and perform quality control -- user-specified pressure 
    grid should be identical to standard atmosphere pressures, and we 
    should make the pressures used in ABSCO_tables.py monotonically 
    decreasing
    """

    inP = np.loadtxt(self.pfile)

    # is pressure ascending or descending? force descending
    # (i.e., go from surface to TOA)
    pDiff = np.diff(inP)
    if (pDiff > 0).all():
      inP = inP[::-1]
    elif (pDiff < 0).all():
      pass
    else:
      sys.exit('Please provide monotonic pressures')
    # endif ascend

    self.pressures = np.array(inP)
    nP = inP.size

    datVMR = pd.read_csv(self.vmrfile)['P'].values
    pErrMsg = 'Inconsistent pressure inputs.  Check %s and %s' % \
      (self.pfile, self.vmrfile) 
    if nP != datVMR.size: sys.exit(pErrMsg)

    return self
  # end processP()

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
    regions['H2CO'] = [[25919, 33300]]
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
    regions['ISOP'] = [[850, 1100], [2800, 3200]]
    regions['PAN'] = [[560, 1400], [1650, 1900]]
    regions['HDO'] = [[1100, 1800], [2500, 3000], [3300, 4300], \
      [4800, 5400], [7000, 7500]]
    regions['BRO'] = [[25927, 34919]]

    # find active species inside specified bands
    # any mol that is active in ANY user-specified band is returned
    wn1, wn2 = self.channels['wn1'], self.channels['wn2']
    activeMol = []
    for iWN, fWN in zip(wn1, wn2):
      print('finding active molecules in %.2f-%.2f cm-1 range' % \
        (iWN, fWN))

      for key in regions.keys():
        active = False
        for band in regions[key]:
          # is there any *overlap*? (user band does not have to be 
          # contained inside active band, just overlap)
          if (fWN >= band[0]) & (iWN <= band[1]): active = True
        # end band loop
        if active: activeMol.append(key)
      # end key loop
    # end WN loop

    self.regActive = dict(regions)
    self.molActive = np.unique(np.array(activeMol))

    return self
  # end findActiveMol()

  def molProcess(self):
    """
    Determine for which bands LNFL and LBLRTM should be run based on 
    the user-input ranges and whether the given molecules are active
    """

    # what molecules are active to the wn range specified?
    self.findActiveMol()
    molActive = np.array(self.molActive)
    molUser = self.molnames

    if len(molUser) == 0:
      print('No active molecules specified.')
      prompt = 'The following molecules are active in the ' + \
        'first input spectral range and can be processed, ' + \
        '%s ' % molActive
      self.show_prompt(prompt)

      molecules = list(molActive)
    else:
      # did the user neglect any active molecules?
      molMissed = []
      for active in molActive:
        if active not in molUser: molMissed.append(active)
      # end active loop

      if len(molMissed):
        prompt = 'The following molecules are active in one ' + \
          'of the input spectral ranges and were not ' + \
          'included in the configuration file: ' + \
          '%s ' % molMissed
        self.show_prompt(prompt)

      # endif molMissed

      # did the user include and non-active molecules?
      molExtra = []
      for mol in molUser:
        if mol not in molActive: molExtra.append(mol)
      # end mol loop

      if len(molExtra):
        prompt = 'The following molecules are not active in ' + \
          'any of the input spectral ranges and were ' + \
          'included in the configuration file: ' + \
          '%s ' % molExtra
        self.show_prompt(prompt)

      # endif molExtra

      molecules = list(molUser)
    # endif no mol

    setattr(self, 'molnames', molecules)

    wn1 = self.channels['wn1']
    wn2 = self.channels['wn2']

    # for each mol to be processed, determine whether it is active in 
    # each of the bands provided by the user
    doBandDict = {}
    for mol in self.molnames:
      doBand = []
      actBands = np.array(self.regActive[mol])

      # user band loop
      for uiWN, ufWN in zip(wn1, wn2):

        # active band loop
        bandOver = []
        for actBand in actBands:
          # any overlap between user band and molecule active band?
          overlap = (ufWN >= actBand[0]) & (uiWN <= actBand[1])
          bandOver.append(overlap)

          # sanity check
          #print(mol, uiWN, ufWN, actBand[0], actBand[1], overlap)
        # end active band loop

        # if there is user overlap with any of the active bands, go 
        # ahead with the band processing
        status = True if any(bandOver) else False
        doBand.append(status)

      # end user band loop
      doBandDict[mol] = list(doBand)
    # end mol loop

    setattr(self, 'doBand', dict(doBandDict))

    return self
  # end molProcess()

  def xsCheck(self):
    """
    Determine if each of the specified molecules should be treated 
    like a XS or line parameter molecule.

    Right now this is not robust enough to handle 2 different statuses
    for a single molecule (e.g. HITRAN recommends XS for 5 of the 6 
    CH3CN bands and line parameters for the remaining band), but that
    feature is also not necessary given the allowed molecules

    Also not robust enough to handle a user-input band spanning line 
    parameter space into XS space
    """

    # user-specified bands
    userWN1 = self.channels['wn1']
    userWN2 = self.channels['wn2']

    # extract the molecules for which line parameters are recommended
    # over the XS coefficients
    csvDat = pd.read_csv(self.xs_lines)
    csvNames = csvDat['species'].values
    xsWN1 = csvDat['wn1'].values
    xsWN2 = csvDat['wn2'].values
    csvFlags = csvDat['flags'].values

    # specify XS processing for each band and each molecule
    doXS = {}
    for mol in self.molnames:
      # determine number of bands in XS database for given molecule
      # by first locating rows corresponding to mol
      iMatch = np.where(mol == csvNames)[0]

      bandXS = []
      for uwn1, uwn2 in zip(userWN1, userWN2):

        if (mol in self.xsNames):
          # no doubter: molecule is only XS
          bandXS.append(True)
        else:
          # this gets complicated...if the molecule has both XS and 
          # line parameters, and there is overlap between the "user"
          # and "XS" bands, and HITRAN recommends XS over lines, doXS
          if (mol in self.xsLines):
            molBandXS = False

            # check if any row in iMatch satisfies criteria
            for iRow in iMatch:
              xwn1, xwn2, flag = \
                xsWN1[iRow], xsWN2[iRow], csvFlags[iRow]

              if (xwn2 >= uwn1) and (xwn1 <= uwn2) and (flag != 1):
                molBandXS = True
                break
              # endif 3-conditional
            # end row loop

            bandXS.append(molBandXS)

          else:
            # another no-doubter: molecule has only line parameters
            bandXS.append(False)
          # endif xsLines
        # endif xsNames
      # end user band loop

      doXS[mol] = list(bandXS)

    # end mol loop

    self.doXS = dict(doXS)

    return self
  # end xsCheck

  def degradeWeighting(self):
    """
    Generate weighting kernels for the spectral degradation
    """

    # the kernels can be of different size, depending on the 
    # degradation factor, so let's use a list instead of arrays
    # 1 kernel per band
    weights = []
    for kWidth in self.channels['kernwidth']:
      # number of kernel elements
      nKernEl = kWidth + 1

      # denominator for weights
      denom = 2**(np.log2(kWidth)-1) + 1

      # not sure what the mathematical term is, but we want to 
      # increase by 1 until a max, then reverse by 1 (e.g. 1 2 3 2 1)
      initKern = np.arange(nKernEl)+1
      revKern = initKern[::-1]
      middle = np.median(initKern)
      initKern[initKern > middle] = revKern[initKern > middle]
      weights.append(np.array(initKern)/denom**2)
    # end kWidth loop

    self.kernel = list(weights)

    return self
  # end degradeWeighting()

  def determineSources(self):
    """
    Determine the source for each band, then combine into a single 
    string for output netCDF
    """

    sources = {}
    for mol in self.molnames:
      molXS = list(self.doXS[mol])
      molSrc = []
      for iBand, doXS in enumerate(molXS):
        src = 'HITRAN2012' if doXS else 'AER v3.6'
        molSrc.append('Band %d: %s' % (iBand+1, src))
      # end band loop
      sources[mol] = '; '.join(molSrc)
    # end mol loop

    self.sources = dict(sources)
  # end determineSources()

  def calcRAM(self):
    """
    Determine the approximate amount of RAM needed for the entirety 
    of the ABSCO computation

    This calculation is simply 8 bytes (np.float64.itemsize) x 
    2 (wavenumber and absco) x number of output wavenumbers (degraded
    spectrum) for all bands combined x number of temperatures x 
    number of layers

    Less HD space is needed because of netCDF compression

    If multiple molecules are provided in the .ini file, they will be
    handled in series, so the by-molecule RAM usage is of greater 
    importance unless the code is run in parallel
    """

    # how many total spectral points?
    nOutWN = (self.channels['wn2']-self.channels['wn1'])
    nOutWN /= self.channels['outres']
    nOutWN = nOutWN.sum()

    estimate = 16 * 2 * nOutWN * 15 * self.pressures.size
    estimateGB = estimate/1e9

    print()
    print('RAM per molecule:')
    print()

    totEst = 0
    for mol in self.molnames:
      # another factor of two to account for 2 H2O VMRs
      scale = 2 if mol in self.molH2O else 1
      est = scale * estimateGB
      totEst += est
      print('\t%s, %.3f GB' % (mol, est))
    # end mol loop

    prompt = 'Full ABSCO table netCDF generation expected to ' + \
      'consume up to %.3f GB of RAM' % totEst
    self.show_prompt(prompt)
  
  def show_prompt(self, message):
    if self.prompt_user:
      message += ', proceed (y/n)? '
      status = input(message).upper()
      if status != 'Y': sys.exit('Exited without proceeding')
    else:
      print(message)

  # end calcRAM()
# end configure()

