import os, sys
import numpy as np

class constants():
  def __init__(self, cgs=False):
    """
    Just a bunch of physical constants that we'll use at some point
    """

    # all attributes in mks units
    self.kB = 1.38064852e-23 # Boltzman's
    self.g = 9.80665 # acceleration due to gravity
    self.R = 2.87058e2 # dry-air gas constant
    self.nA = 6.022140857e23 # Avogadro's number
    self.c = 2.99792458e8 # speed of light, vaccum
    self.SB = 5.67051e-8 # Stefan-Boltzman (Wm-2K-4)
    self.sPerDay = 60 * 60 * 24 # seconds per day
    self.Rydberg = 8.31 # J mol-1 K-1

    # wet and dry masses are from a KCP script check_wtot2.pro
    # they were used for RFMIP layer density calculations since we
    # wanted to use hydrostatics and not the ideal gas law
    self.mH2O = 1.8016e-2
    self.mDry = 2.8964e-2

    if cgs: self.mks2cgs()

  # end constructor

  def mks2cgs(self):
    self.kB = 1.38064852e-16 # Boltzman's
    self.g = 9.80665e2 # acceleration due to gravity
    self.R = 2.87058e6 # dry-air gas constant
    self.c = 2.99792458e10 # speed of light, vaccum
    self.mH2O = 1.8016e1
    self.mDry = 2.8964e1
    self.SB = 5.67051e-5 # Stefan-Boltzman

    return self
  # end mks2cgs
# end constants

def readTAPE7(inFile, header=True, xs=False):
  """
  Read in a single TAPE7 (LBLRTM profile/layer amounts as calculated
  by LBLATM subroutine) and return parameters in dictionary

  Call
    outDict = readTAPE7(inFile)

  Input
    inFile -- string, path to TAPE7 file

  Output
    outDict -- dictionary with the following keys:
      format: int, pressure format specification
      n_layers: int, number of layers
      n_molecules: int, number of molecules specified in profile
      scale_factor: float array, secant scaling factor (nLayers)
      end_alt: float array, instrument altitude
      obs_alt: float array, observer altitude
      ang: float, viewing angle (nadir = 180)

      p_lay: float array, average pressure of given layer (nLayers)
      T_lay: float array, average temperature of given layer (nLayers)
      type_lay: int array, path type
      path_lay: int array, direction of path (nLayers)

      p_lev: float array, pressure at layer boundaries (nLevels)
      alt_lev: float array, altitudes at layer bounds (nLevels)
      T_lev: float array, temperature at layer bounds (nLevels)

      scale_factor_lay: float array, secant scaling factor for each
        layer (nLayers)
      molDen: float array, molecule density for each gas at each layer
        (nMol-1 x nLayer; the -1 is for the broadening density)
        molecule numbers follow the HITRAN convention --
        https://hitran.org/lbl/ -- but are 0-offset (e.g., molecule
        0 is H2O, 1 is CO2, etc.)

  Keywords
    header -- boolean, does the TAPE7 have the traditional comment
      at the beginning (starts with "$")
    xs -- boolean, extract XS profiles INSTEAD OF line-by-line
      molecule profiles
  """

  def stringSlice(inStr, idxArr):
    """
    Return substring of inStr that spans indices from idxArr
    """

    outStr = inStr[idxArr.min():idxArr.max()+1]

    # if it's an empty string, replace with zero
    outStr = 0 if len(outStr.strip()) == 0 else outStr

    return outStr
  # end stringSlice()

  # skip header of TAPE7, keep everything else
  datT7 = open(inFile).read().splitlines()
  if header: datT7 = datT7[1:]

  record21 = datT7[0]

  iForm = int(stringSlice(record21, np.array([0, 1])))
  nLay = int(stringSlice(record21, np.array([2, 5])))
  nMol = int(stringSlice(record21, np.array([5, 10])))
  secnto = float(stringSlice(record21, np.array([10, 19])))
  h1 = float(stringSlice(record21, np.array([40, 48])))
  h2 = float(stringSlice(record21, np.array([52, 60])))
  ang = float(stringSlice(record21, np.array([65, 73])))

  # the number of "layer lines" is dependent on the number of
  # molecules. the convention is dictated by Record 2.1.2 in the
  # LBLRTM instructions HTML file (8 molecules per line)
  # nLayLines = P/T line + Mixing Ratios lines (records 2.1.1 + 2.1.2)
  if nMol <= 7:
    nLayLines = 1
  else:
    # +1 for the broadener ("molecule" 8)
    nLayLines = np.ceil((nMol+1)/8.0)
  # endif nMol

  # for the P/H/T line
  nLayLines += 1

  # number of molecule/LBL data lines (i.e., non-XS) = layers + r21 + header
  nLinesLBL = int(nLayLines * nLay + 2)

  # either molecules or cross-sections profile extraction
  if xs:
    # have to reassign some variables for XS
    profile = datT7[nLinesLBL-1:]
    nMol = int(profile[0].split()[0])

    # "record 3.7.1" is identical to record21 and redundant; skip XS
    # header, XS names (variable number of lines), and "record 3.7.1"
    # TO DO: have not tested nRec371 past nMol=8
    nRec371 = nMol // 8 + 1
    xsNames = ' '.join(profile[1:1+nRec371]).split()

    # "OTHER" xs is always the broadener (molecules and xs)
    iBroad = xsNames.index('OTHER')
    xsNames[iBroad] = 'broadener'

    profile = profile[2+nRec371:]

    # TO DO: test robustness of this and nRec371 for XS
    nLayLines = 1 if nMol <= 7 else int(np.ceil((nMol+1)/8.0))
    nLayLines += 1
  else:
    profile = datT7[1:nLinesLBL-1]
  # endif xs

  # how the data are read (i.e., array slicing) depends on iForm
  # this is record 2.1.1
  if iForm == 0:
    ipLay = np.array([0, 10])
    itLay = np.array([10, 21])
    iSecant = np.array([21, 30])
    iType = np.array([30, 33])
    iPath = np.array([33, 35])
    iAlt1 = np.array([36, 43])
    ipLev1 = np.array([43, 51])
    itLev1 = np.array([51, 58])
    iAlt2 = np.array([58, 65])
    ipLev2 = np.array([65, 73])
    itLev2 = np.array([73, 80])
  else:
    ipLay = np.array([0, 15])
    itLay = np.array([15, 25])
    iSecant = np.array([25, 35])
    iType = np.array([35, 38])
    iPath = np.array([38, 40])
    iAlt1 = np.array([41, 48])
    ipLev1 = np.array([48, 56])
    itLev1 = np.array([56, 63])
    iAlt2 = np.array([63, 70])
    ipLev2 = np.array([70, 78])
    itLev2 = np.array([78, 85])
  # endif iForm

  # assemble lists for each profile parameter
  pLay, tLay, scaleLay, typeLay, pathLay, altLev, pLev, tLev = \
    ([] for i in range(8))

  # molecule volume mixing ratios will first be a list of lists
  # "sub" den is a subset of all densities for a given layer
  molDen, subMolDen = [], []
  for iLine, line in enumerate(profile):
    if iLine % nLayLines == 0:
      # layer P/T/Z info
      pLay.append(stringSlice(line, ipLay))
      tLay.append(stringSlice(line, itLay))
      scaleLay.append(stringSlice(line, iSecant))
      typeLay.append(stringSlice(line, iType))
      pathLay.append(stringSlice(line, iPath))

      if iLine == 0:
        # the first layer has the info for the first 2 boundaries
        altLev.append(stringSlice(line, iAlt1))
        pLev.append(stringSlice(line, ipLev1))
        tLev.append(stringSlice(line, itLev1))
      # endif iLine

      altLev.append(stringSlice(line, iAlt2))
      pLev.append(stringSlice(line, ipLev2))
      tLev.append(stringSlice(line, itLev2))

      # reset this guy every layer
      subMolDen = []
    else:
      # layer molecule amounts
      subMolDen += line.split()

      # are we on the last line of the layer?
      if iLine % nLayLines == nLayLines-1: molDen.append(subMolDen)
    # end modulo 0
  # end layer loop

  # convert lists to arrays
  pLay, tLay, scaleLay, typeLay, pathLay, altLev, pLev, tLev, molDen = \
    np.array(pLay), np.array(tLay), np.array(scaleLay), \
    np.array(typeLay), np.array(pathLay), np.array(altLev), \
    np.array(pLev), np.array(tLev), np.array(molDen)

  outDict = {'n_layers': nLay, 'n_molecules': nMol, 'format': iForm, \
    'scale_factor': secnto, 'obs_alt': h1, 'end_alt': h2,
    'view_angle': ang}

  # extract broadening density from molecule den and transpose density
  # so it is nMol x nLay
  # molecule 8 is always the broadener (molecules and xs)
  if not xs: iBroad = 7

  broadener = molDen[:, iBroad]
  iDen = np.delete(np.arange(nMol+1), iBroad)
  molDen = molDen[:, iDen].T

  # make a list of all lists, loop through it, and convert all lists
  # to arrays of the proper type and stuff them into outDict
  # TO DO: make outDict['densities'] its own dictionaries with
  # LBL molecule names like we do for XS
  dictKeys = ['p_lay', 'T_lay', 'scale_factor_lay', 'type_lay', \
    'path_lay', 'alt_lev', 'p_lev', 'T_lev', 'densities', 'broadener']
  tempList = [pLay, tLay, scaleLay, typeLay, pathLay, altLev, \
    pLev, tLev, molDen, broadener]
  for iKey, temp in enumerate(tempList):
    if iKey == 'densities' and xs: continue
    outDict[dictKeys[iKey]] = temp.astype(float)

  if xs:
    outDict['densities'] = {}
    xsNames.remove('broadener')

    for iKey, key in enumerate(xsNames):
      outDict['densities'][key] = molDen[iKey].astype(float)
  # endif xs
  return outDict
# end readTAPE7

def readTAPE28(inFile, nSkip=52):
  """
  Read in a single TAPE28 (LBLRTM Brightness Temperature output file)
  and return spectrum in dictionary

  Call
    outDict = readTAPE28(inFile)

  Input
    inFile -- string, path to TAPE28 file

  Output
    outDict -- dictionary with wavenumber and brightness_temperature
      key/value pairs

  Keywords
    nSkip -- int, the number of lines in the header
  """

  # very, very simple function -- maybe just easier to use np.loadtxt
  # directly
  waveNum, bt = np.loadtxt(inFile, unpack=True, skiprows=nSkip)
  outDict = \
    {'wavenumber': waveNum , 'brightness_temperature': bt, 'units': 'K'}

  return outDict
# end readTAPE28

def readTAPE27(inFile, nSkip=52):
  """
  Read in a single TAPE27 (LBLRTM radiance ASCII output file)
  and return spectrum in dictionary

  Call
    outDict = readTAPE27(inFile)

  Input
    inFile -- string, path to TAPE27 file

  Output
    outDict -- dictionary with wavenumber and brightness_temperature
      key/value pairs

  Keywords
    nSkip -- int, the number of lines in the header
  """

  # very, very simple function -- maybe just easier to use np.loadtxt
  # directly
  waveNum, rad = np.loadtxt(inFile, unpack=True, skiprows=nSkip)
  outDict = {'wavenumber': waveNum , 'radiance': rad, \
    'units': '(W cm-2 sr-1)/cm-1'}

  return outDict
# end readTAPE27

def readBinary(inFile, double=True):
  """
  Read LBLRTM binary file (these are special unformatted binary files,
  written in "panel" format)

  Boo...didn't work with ASTI TAPE11...

  Input
    inFile -- string, path to binary TAPE (10, 11, 12, 13)
      output_file from LBLRTM

  Output
    outWN -- float array, wavenumbers spanning spectrum
    param -- float array of ODs, radiances, fluxes, transmittances,
      or whatever other paramter is extracted from inFile

  Keywords
    double -- boolean, is inFile in double precision? defaults to yes
  """

  # probably want something like this, from my ABSCO_diagnostics:
  """
  with open(aFile, "rb") as f:
    # following the procedure i used in xsABSCO.postProcessXS()
    # when i write the binary that i am currently reading

    # wtf variable
    dummy = f.read(4)

    # the first 3 "panel headers" that i wrote
    panelHeader = array.array('d')
    panelHeader.fromfile(f, 178)
    wnDat = array.array('d')
    wnDat.fromfile(f, 3)
    wnDat = np.array(wnDat)
    numFreq = array.array('l')
    numFreq.fromfile(f, 2)
    nFreq = np.array(numFreq)[0]

    # now read in the ABSCO array, which is dependent on nFreq
    abscoArr = array.array('d')
    abscoArr.fromfile(f, nFreq)
    abscoArr = np.array(abscoArr)

    # don't need any of the rest of the "panel header" garbage
    # but we do need to make a wavenumber array associated w/
    # absorption coefficients
    waveNum = np.arange(wnDat[0], wnDat[1]+wnDat[2], wnDat[2])

    # save the spectrum
    abscoList.append(abscoArr)
    wnList.append(waveNum)
  """

  from . import FortranFile
  from .lblTools import readTape12

  outWN, param = readTape12(inFile, double=double)

  return np.array(outWN), np.array(param)
# end readBinary()

def tempIDL(inFile, fType=0, double=True):
  """
  Read in binary TAPE files (inFile), save data to IDL save files,
  and then read them with Python for plotting. this is pretty time-
  consuming

  This is a temporary function until I figure out how to read in
  FORTRAN binary files in Python (looks like it can be done --
  https://stackoverflow.com/questions/37534220/python-read-fortran-binary-file)

  Input
    inFile -- string, path to binary TAPE (10, 11, 12, 13)
      output_file from LBLRTM

  Output

  Keywords
    fType -- int, file type (radiance, transmittance, etc.; see doc in
      /project/rc/rc2/mshep/idl/patbrown/read_lbl_file_dbl.pro)
    double -- boolean, is inFile in double precision? defaults to yes

  """

  from scipy.io.idl import readsav
  from . import utils

  # write_save_file.pro is in this Git repo:
  # https://lex-gitlab.aer.com/rpernak/common_modules
  # and is considered, along with the RC_utils.py and utils.py
  # modules, part of the RC common library
  proFile = 'write_save_file.pro'
  #if not os.path.exists(proFile):
  #  os.symlink('externals/common/%s' % proFile, proFile)

  if double:
    proCall = \
      "write_save_file, \'%s\', file_type=%d, /dbl" % (inFile, fType)
  else:
    proCall = \
      "write_save_file, \'%s\', file_type=%d" % (inFile, fType)
  # endif double

  idlCmd = 'idl -e "%s"' % proCall
  sOut, sErr = utils.spawn(idlCmd)

  # write_save_file.pro always writes a LBLRMT_output.sav file
  # and contains the wavenum and spectrum arrays
  tempSav = 'LBLRTM_output.sav'
  idlDat = readsav(tempSav)
  waveNum, param = idlDat['wavenum'], idlDat['spectrum']
  os.remove(tempSav)
  #os.remove(proFile)

  return {'wavenumber': waveNum, 'spectrum': param}
# end tempIDL()

def radsumRead(inFile):
  """
  Read a single RADSUM output file and return data for a given level
  as a dictionary to be used in radsumPlot()

  Call
    outDict = radsumRead(inFile)

  Input
    inFile -- string, path to RADSUM output file

  Output
    outDict -- dictionary with the following keys
      (with float list values):

      up_flux: upwelling flux (W/m2) as a function of wavenumber
        and level (nLevel x nWavenumber array)
      down_flux: downwelling flux (W/m2) as a function of wavenumber
        and level (nLevel x nWavenumber array)
      net_flux: net flux (W/m2) as a function of wavenumber and level
        (nLevel x nWavenumber array)
      heat_rate: heating rate (K/day) as a function of wavenumber
        and level (nLevel x nWavenumber array)
      wavenumber: spectral points (cm-1) vector (1 x nWavenumber)
      level_pressure: pressure at layer boundaries
        (nLevel-element array)

  Keywords
    None
  """

  inDat = open(inFile).read().splitlines()

  # initialize lists that eventually become
  waveNum1, waveNum2 = [], []
  pLevAll, upFluxAll, dnFluxAll, netFluxAll, heatRateAll = \
    ([] for x in range(5))

  for line in inDat:
    split = line.split()

    # for each band, deduce if we are processing the output or header
    try:
      # if this works, proceed to parsing the rest of the line
      iLev = int(split[0])
    except:
      # header processing -- only wanna extract wavenumber
      if len(split) > 0:
        if split[0] == 'WAVENUMBER':
          # every WAVENUMBER string occurrence implies the start of
          # a new block of RADSUM output, which we break up into
          # 1 x nLev vectors for each parameter, then append to
          # *All lists that eventually become
          # (nWavenumber x nLev) arrays of output
          if 'pLev' in locals():
            # are we past the first block of output (so pLev exists)?
            # otherwise this is unnecessary
            pLevAll.append(pLev)
            upFluxAll.append(upFlux)
            dnFluxAll.append(dnFlux)
            netFluxAll.append(netFlux)
            heatRateAll.append(heatRate)

            pLev, upFlux, dnFlux, netFlux, heatRate = \
              ([] for x in range(5))
          else:
            pLev, upFlux, dnFlux, netFlux, heatRate = \
              ([] for x in range(5))
          # endif pLev len

          # there is no waveNum "array", just an nLev-element vector
          waveNum1.append(float(split[2]))
          waveNum2.append(float(split[4]))
        # endif WAVENUMBER
      # endif split len

      continue
    # end exception

    # RADSUM output processing
    # sometimes radsum has bad pressures because of string formatting
    try:
      pLev.append(float(split[1]))
    except:
      pLev.append(np.nan)

    upFlux.append(float(split[2]))
    dnFlux.append(float(split[3]))
    netFlux.append(float(split[4]))
    heatRate.append(float(split[5]))
  # end loop over lines

  # add last output block to arrays
  pLevAll.append(pLev)
  upFluxAll.append(upFlux)
  dnFluxAll.append(dnFlux)
  netFluxAll.append(netFlux)
  heatRateAll.append(heatRate)

  outDict = {'wavenumber1': np.array(waveNum1), \
    'wavenumber2': np.array(waveNum2), \
    'level_pressure': np.array(pLevAll), \
    'up_flux': np.array(upFluxAll), \
    'down_flux': np.array(dnFluxAll), \
    'net_flux': np.array(netFluxAll), \
    'heat_rate': np.array(heatRateAll)}

  return outDict
# end radsumRead()

def readRRTM(inFile):
  """
  Read RRTM input and return dictionary of parameters (pressure, up
  flux, diffuse down flux, direct down flux, total down flux, net
  flux, and heating rate

  Inputs
    inFile -- str, OUTPUT_RRTM file (output from RRTM run)

  Outputs
    outDict -- dictionary with keys (values are nLev x nBand float
      arrays):

      pressure [mbar]
      up flux [W/m2]
      diffuse down flux [W/m2]
      direct down flux [W/m2]
      total down flux [W/m2]
      net flux [W/m2]

      also contains corresponding broadband arrays (BB is appended)
  """

  dat = open(inFile).read().splitlines()

  # initialize parameter lists
  bandWN1, bandWN2 = [], []
  up, diffuse, direct, down, net, hr = \
    [], [], [], [], [], []
  pBand, upBand, difBand, dirBand, downBand, netBand, hrBand = \
    [], [], [], [], [], [], []

  for line in dat:
    if len(line) == 0 or len(line) == 1: continue
    split = line.split()

    # end of file (don't need anything after this string)
    if split[0] == 'Modules': break

    if split[0] == 'Wavenumbers:':
      # band header
      bandWN1.append(float(split[1]))
      bandWN2.append(float(split[3]))

      # re-initiate this dictionary for every band
      bandDict = {}
      continue
    # endif band header

    if len(split) == 8:
      # extract fluxes for all levels in a given band
      pBand.append(float(split[1]))
      upBand.append(float(split[2]))
      difBand.append(float(split[3]))
      dirBand.append(float(split[4]))
      downBand.append(float(split[5]))
      netBand.append(float(split[6]))
      hrBand.append(float(split[7]))
    # endif flux extract

    if split[0] == '0':
      # surface-level fluxes => end of band
      # save band fluxes
      up.append(upBand)
      diffuse.append(difBand)
      direct.append(dirBand)
      down.append(downBand)
      net.append(netBand)
      hr.append(hrBand)
      pressure = pBand

      # re-initialize flux/hr lists
      pBand, upBand, difBand, dirBand, downBand, netBand, hrBand = \
        [], [], [], [], [], [], []
    # endif surface
  # end dat loop

  up = np.array(up)
  direct = np.array(direct)
  diffuse = np.array(diffuse)
  down = np.array(down)
  net = np.array(net)
  hr = np.array(hr)

  # separate output into by-band and broadband arrays
  # broadband is first in the OUTPUT_RRTM files
  outDict = {'pressure': np.array(pressure), \
    'up': up[1:,:], 'upBB': up[0,:], \
    'direct': direct[1:,:], 'directBB': direct[0,:], \
    'diffuse': diffuse[1:,:], 'diffuseBB': diffuse[0,:], \
    'down': down[1:,:], 'downBB': down[0,:], \
    'net': net[1:,:], 'netBB': net[0,:], \
    'heating_rate': hr[1:,:], 'heating_rate_BB': hr[0,:], \
    'band_lims': np.array([bandWN1[1:], bandWN2[1:]]).T, \
    'band_lims_BB': np.array([bandWN1[0], bandWN2[0]])}

  return outDict
# end readRRTM()

def readXS(inFile, speciesXS):
  """
  Read in the absorption coefficients from HITRAN .xsc files.

  Probably can eventually add some flexibility here for AER LBLRTM
  XS files

  Input
    inFile -- string, HITRAN .xsc file (e.g., CCl4_IR00.xsc)
      eventually: AER xs file (e.g., xs/CCL4) as well
    speciesXS -- string, name of species that is being processed
      (HITRAN "Common Name" convention, e.g., CFC-12, HCFC-22, etc.
       see /nas/project/rc_static/line_files/line_parameters_HITRAN/
       hitran2012/IR-XSect/IRCrossSection_Readme.pdf)

  Output
    eh...working on this. originally thought they'd be float arrays
    but we may need dictionaries because of different sizes of spectra

    outWN -- float array, wavenumbers for spectrum (1-D)
    outK -- float array, absorption coefficient spectrum
      dimensions for both are determined from the number of P/T
      combinations (a proxy for this can be the number of headers in
      the file) and the number of spectral points (which is
  """

  nBlocks = 0
  outWN, outK = {}, {}

  dat = open(inFile).read().splitlines()
  for line in dat:
    split = line.split()

    if speciesXS in line:
      # every header has the species name in it
      # data blocks only have absorption coefficients
      nBlocks += 1

      # when we get to a header, we have to store the previous
      # data block (key should exist by now because it's created
      # from the header)
      if len(outWN.keys()) != 0: outK[key] = np.array(kArr)

      wn1 = float(split[1])
      wn2 = float(split[2])
      nPoints = int(split[3])
      temperature = '%7.2f' % float(split[4])
      pressure = '%6.2f' % float(split[5])
      specRes = (wn2-wn1)/(nPoints-1)
      wnArr = wn1 + np.arange(nPoints) * specRes
      wn1 = '%10.4E' % wn1
      wn2 = '%10.4E' % wn2
      key = '%s/%s/%s/%s' % \
        (wn1.strip(), wn2.strip(), temperature.strip(), pressure.strip())
      key = key.strip()

      outWN[key] = wnArr
      kArr = []
    else:
      # kArr should have been initialized in header processing
      # concatenate (NOT append) to it
      try:
        kList = [float(k) for k in split]
        kArr += kList
      except:
        print('%s may not be the correct species for %s' % \
          (speciesXS, inFile))
        sys.exit(1)
      # end trying
    # endif species check (header)
  # end data loop

  # save the last data block
  outK[key] = np.array(kArr)

  # i suspect that (at least in H16) that sometimes zeroes are used
  # as fill values on the final line such that the number of
  # absorption coefficients is not consistent with the number of
  # spectral points as specified in the header (i ran into this issue
  # with H16 CCl4)
  outK[key] = outK[key][:outWN[key].size]

  if len(outWN.keys()) == 0:
    print('%s not found in %s' % (speciesXS, inFile))
    return {}, {}
  # endif headers

  return outWN, outK
# end readXS()

def rad2BT(inWN, inRad):
  """
  Radiance to Brightness Temperature conversion courtesy of
  https://ncc.nesdis.noaa.gov/data/planck.html

  Input
    inWN -- float array, wavenumbers of spectrum (cm-1)
    inRad -- float array, associated radiance at each wavenumber
      (standard RU used in LBLRTM: W cm-2 sr-1 / cm-1)
  Output
    outBT -- float array, brightness temperatures (K)
  """

  # convert from standard RU to mW m-2 sr-1 / cm-1
  inRad = inRad / 1.0e-7

  num = 1.4387752 * inWN
  denom = np.log( (1.191042e-5 * inWN**3 / inRad) + 1 )
  outBT = num/denom

  return outBT
# end rad2BT()

def colAmt2PWV(amount):
  """
  Convert accumulated molecular amounts for total path (mol/cm2) in
  TAPE6 to precipitable water vapor (PWV, cm)

  The conversion formula was gathered from an email with Vivienne
  Payne ("update to conversion factor for H2O") to the AER RC email
  group on 24-Jul-2006:

  pwv (cm) = \
    (column amount)(1/avogadro)(gram molec wt h2o)(1/sp density h2o)
           = [column amnt  (molec/cm^2) ] x  2.99150e-23 (cm^3/molec)

  Call
    pwv = colAmt2PWV(amount)

  Input
    amount -- float array, column amounts for a given molecule
      (mol/cm2)

  Output
    pwv -- float array, corresponding precipitable water vapor
      (cm)
  """

  return amount * 2.99150e-23
# end colAmt2PWV()

def wvAmtTAPE6(inTAPE6):
  """
  Extract water vapor accumulated (over entire column) amount from
  inTAPE6
  """

  dat = open(inTAPE6).read().splitlines()
  search = '%-55s' % ('0')
  search += 'ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH'
  for iLine, line in enumerate(dat):
    if search in line:
      wvLayAmt = float(dat[iLine+1][57:70])
    else:
      continue
    # endif search

    # break out of the loop as soon after first ACCUMULATED line
    # (otherwise we will grab other molecule densities)
    break
  # end

  return wvLayAmt
# end wvAmtTAPE6

def fluxToHR(flux):
  """
  Flux-to-heating rate calculation using the Stefan Boltzman law

  Input
    flux -- float array, fluxes in Wm-2
    ASSUMED TO BE NLAY X NPROFILE X NBAND! or at least NLAY is first
    dimension

  Output
    hr -- float array, corresponding heating rates in T day-1
  """

  conObj = constants()

  return (np.diff(flux, axis=0) / conObj.SB)**(1/4.) / conObj.sPerDay
# end fluxToHR()
