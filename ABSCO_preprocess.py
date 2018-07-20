# standard Python libraries
# for Python 3 compatibility
from __future__ import print_function
import sys

# miniconda-installed libs
import numpy as np

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
      'NO', 'SO2', 'NO2', 'NH3', 'HNO3', 'OCS', 'H2CO', 'N2', \
      'HCN', 'C2H2', 'HCOOH', 'C2H4', 'CH3OH', 'CCL4', 'CF4', \
      'F11', 'F12', 'F22', 'ISOP', 'PAN', 'HDO', 'BRO', 'O2-O2']

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
          split = cItem[1].split()
          if len(split) == 0:
            # CONSIDER DEFAULTS
            chanErrMsg = 'No bands specified in %s, returning' % \
              inFile
            sys.exit(chanErrMsg)
          # endif split

          channels[cItem[0]] = np.array(split).astype(float)
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
        # wavenumbers and associated spectral resolutions
        # need to also check that each wn1 < wn2
        for key in keys[1:]:
          if channels[key].size != channels[keys[0]].size:
            chanErrMsg = 'wn1, wn2, and res should have equal ' + \
              'number of elements, returning'
            print('Error in %s' % inFile)
            sys.exit(chanErrMsg)
          # endif channels
        # end key loop

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
      else:
        for cItem in cItems: setattr(self, cItem[0], cItem[1])
      # endif cps
    # end sections loop

    # in the makeABSCO() class, we expect certain attributes
    # let's make sure they exist in the config file
    reqAtt = ['pfile', 'ptfile', 'vmrfile', 'broadfile', 'channels', \
      'molnames', 'scale', 'lnfl_run_dir', 'lnfl_path', \
      'tape1_path', 'tape3_dir', 'extra_params', 'tape5_dir', \
      'lbl_path', 'xs_path', 'fscdxs', 'lbl_run_dir', 'od_dir', \
      'absco_dir']

    # loop over all required attributes and do a check
    for req in reqAtt:
      if req not in dir(self):
        print(errMsg)
        sys.exit('Could not find %s attribute, returning' % req)
      # endif req
    # end req loop

    # let's pack all of the files into a single list
    self.paths = [self.pfile, self.ptfile, self.vmrfile, \
      self.extra_params, self.lnfl_path, self.lbl_path, \
      self.xs_path, self.fscdxs]
    self.outDirs = [self.lnfl_run_dir, self.lbl_run_dir, \
      self.tape3_dir, self.tape5_dir, self.od_dir, self.absco_dir]

    # cross section molecules will have to be handled differently
    # from HITRAN molecules
    self.xsNames = \
      ['CF4', 'CCL4', 'F11', 'F12', 'F22', 'ISOP', 'PAN']

    # and these guys are neither HITRAN or XS molecules
    self.dunno = ['HDO', 'O2-O2', 'BRO']

    # gather necessary params from input files
    # is pressure ascending or descending? force descending
    # (i.e., go from surface to TOA)
    inP = np.loadtxt(self.pfile)

    # is pressure ascending or descending? force descending
    pDiff = np.diff(inP)
    if (pDiff > 0).all():
      inP = inP[::-1]
    elif (pDiff < 0).all():
      pass
    else:
      sys.exit('Please provide monotonic pressures')
    # endif ascend

    self.pressures = np.array(inP)

    # TO DO: ensure user profile has the same amount of pressures
    # could probably do better and return P in profile CSV and then
    # do the 1-1 P check here
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
    regions['O2-O2'] = [[16644, 29785]]

    # find active species inside specified bands
    # any mol that is active in ANY user-specified band is returned
    wn1, wn2 = self.channels['wn1'], self.channels['wn2']
    activeMol = []
    for iWN, fWN in zip(wn1, wn2):
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
    the user-input ranges and whether the given molecules are 
    """

    # what molecules are active to the wn range specified?
    self.findActiveMol()
    molActive = np.array(self.molActive)
    molUser = self.molnames

    if len(molUser) == 0:
      print('No molecules specified, ', end='')
      print('finding active molecules in %.2f-%.2f cm-1 range' % \
        (self.channels['wn1'], self.channels['wn2']))
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
          '%s, proceed (y/n)? ' % molMissed
        status = input(prompt)
        status = status.upper()

        if status == 'N': sys.exit('Exited without proceeding')
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
          '%s, proceed (y/n)? ' % molExtra
        status = input(prompt)
        status = status.upper()

        if status == 'N': sys.exit('Exited without proceeding')
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
# end configure()

