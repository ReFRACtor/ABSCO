#!/usr/bin/env python

# standard Python libraries
# for Python 3 compatibility
from __future__ import print_function

import os, sys, glob, argparse
import numpy as np

# path to GIT common submodules (not a Python standard lib)
sys.path.append('externals/common')
import utils, lblTools

# GLOBAL VARIABLES 
# defaults for xsABSCO keywords
NAMES = np.array(['2B1', '1B1', '1B2', '2A1', '2A2', '2A3', '2A4', \
  '1A1', '1A2', '1A3', '1A4', '1A5'])
WN1 = np.array([620.0, 780.0, 890.0, 1060.0, 1250.0, 1450.0, \
  1640.0, 1860.0 , 2160.0, 2380.0, 2540.0, 2720.0])
WN2 = np.array([950.0, 1090.0, 1190.0, 1370.0, 1610.0, 1810.0, \
  2000.0, 2280.0, 2500.0, 2720.0, 2880.0, 3080.0])
RES = np.array([0.0001, 0.0001, 0.0002, 0.0002, 0.0002, 0.0002, \
  0.0002, 0.0004, 0.0004, 0.0004, 0.0004, 0.0004])
TESCHANNELS = {'names': NAMES, 'wn1': WN1, 'wn2': WN2, 'res': RES}
LBLXS = ['CCL4', 'F11', 'F12', 'CF4', 'CHCLF2', 'ACET', 'ISOP']
LBLXS = [LBLXS[1]]
ELANORXS = ['CCL4', 'CFC11', 'CFC12', 'CFC14', 'CFC22', \
  'ACET', 'ISOP']
ELANORXS = [ELANORXS[1]]
MOLECULES = {'lblxs': LBLXS, 'elanorxs': ELANORXS}

# paths for LBL runs
RCSTATIC = '/nas/project/rc_static/models'
PATHLBL = RCSTATIC + \
  '/aer_lblrtm/lblrtm_v12.6/lblrtm/lblrtm_v12.6_linux_pgi_dbl'
PATHT3 = '/nas/ASR/LINEFILE_BUILD_TES_V_2.0/lnfl/' + \
  'backup_TAPE3_files/TAPE3_tes_v_2.0'
PATHXS = RCSTATIC + \
  '/aer_line_parameters/AER_line_files/aer_v_3.6/xs_files_v3.6/xs'
PATHFSCDXS = RCSTATIC + \
  '/aer_line_parameters/AER_line_files/aer_v_3.6/xs_files_v3.6/FSCDXS'

class xsABSCO():
  """
  - Build TAPE5s for each XS molecule of interest
  - Run LBLRTM to generate ODInt files
  - etc. etc.

  Generate absorption coefficient files for molecules of interest
  """

  def __init__(self, ptnFile='Pbar_Tbar_Amt.in', \
    tBarFile='Tbar_Grid.in', channels=TESCHANNELS, \
    molecules=MOLECULES, dirT5='TAPE5_dir', vmr='VMR_data.csv', \
    dirDen='XS_densities', pathLBL=PATHLBL, pathT3=PATHT3, \
    pathXS=PATHXS, pathFSCDXS=PATHFSCDXS, dirWork='LBL_runs', \
    dirFine='ABSCO_Fine', dirCoarse='ABSCO_Coarse', odHeader=''):

    """
    Inputs

    Keywords
      ptnFile -- path to ASCII file with layer average pressures and 
        temperatures and associated layer amounts for:
        AIR, H2O, CO2, O3, N2O, CO, CH4, O2, NO, NO2, HNO3, 
        OCS, DRY AIR
      tBarFile -- path to ASCII file with average temperatures as a 
        function of pressure (and corresponding level number)
      dirT5 -- string, path to directory under which subdirectories 
        for each molecule will be made and to which all TAPE5 files 
        will be moved (absolute path is needed)
      vmr -- string, path to file with VMR for many molecules as a 
        function of pressure
      dirDen -- string, name of directory with density text files 
        for each species (output from calcDensity())
      pathLBL -- string, full path to LBLRTM executable to be used
      pathT3 -- string, full path to TAPE3 to use in LBLRTM runs
      pathXS -- string, full path to director with individual XS files
      pathFSCDXS -- string, full path to XS lookup file
      dirWork -- string, working directory name
      dirFine -- string, subdirectory to which LBLRTM ODint files 
        are saved
      dirCoarse -- string, subdirectory to which modified ODint files
        (with coarser spectral resolution) are saved
      odHeader -- string, header that will be placed in 
        each final OD file (on the coarse grid); if it is more than 
        80 characters, only the first 80 will be used

      For channels and molecules, see GLOBAL VARIABLES for examples
        (TESCHANNELS and MOLECULES are the defaults)

      channels -- dictionary with the following keys:
        names -- string array of names for each channel
        wn1 -- float array of starting wavenumber (cm-1) for each band
        wn2 -- float array of ending wavenumber (cm-1) for each band
        res -- float array of spectral resolution of each band
      molecules -- dictionary with the following keys
        lbl -- string list of LBL molecule names
        elanor -- string list of corresponding ELANOR molecule names
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

  def postProcessXS(self):
    """
    Mimic what is done in the process_lblrtm_xs.pro code after optical
    depth files are generated (with runLBL())

      - read in runLBL()-generated ODint files
      - scale ODs with column densities for species
      - sample the spectrum at a coarser resolution

    ASSUMING LITTLE ENDIANNESS!
    """

    # standard library
    import array

    for iMol, mol in enumerate(self.molNames['lblxs']):
      # grab density profile, which will be used for scaling
      denFile = '%s/%s/%s_densities.txt' % \
        (self.topDir, self.dirDensities, mol)
      pProf, denProf = np.loadtxt(denFile, unpack=True)
      pThresh = 0.005

      msgNoOD = 'No OD files for %s, returning' % mol

      # set up working directories
      workSubDir = '%s/%s' % (self.workDir, mol)
      if not os.path.exists(workSubDir):
        print(msgNoOD)
        return
      # endif workSubDir

      outDirOD = '%s/%s' % (workSubDir, self.coarseOD)
      if not os.path.exists(outDirOD): os.mkdir(outDirOD)
      os.chdir(workSubDir)

      # ELANOR molecule name will be used in output files/dirs
      molELANOR = self.molNames['elanorxs'][iMol]

      # grab OD files for given species
      odFiles = sorted(glob.glob(\
        '%s/%s/%s/ODint_*' % (self.topDir, workSubDir, self.fineOD)))
      nFilesOD = len(odFiles)
      if nFilesOD == 0: 
        print(msgNoOD)
        return
      # endif nFilesOD

      for iOD, odFile in enumerate(odFiles):
        print(odFile)
        # we need to extract the pressure associated with this layer
        # and it is located in the naming convention of the OD file
        # see outFile assignment in makeTAPE5(), which propagates to 
        # the OD file naming convention (strip the "P" and "mb" 
        # in the pressure substring)
        base = os.path.basename(odFile)
        split = base.split('_')
        layerP = float(split[2][1:-2])
        layerT = float(split[3][1:-1])
        channel = split[4]
        newName = '%s_P%10.4Emb_T%5.1fK' % (molELANOR, layerP, layerT)

        # grab all "relevant" pressures
        iMatch = np.where(\
          (pProf <= layerP + pThresh) & (pProf >= layerP - pThresh))
        iMatch = iMatch[0]
        if iMatch.size == 0: continue

        # post-processing output directories
        outDirABSCO = '%s/%s/%s/%s' % \
          (self.topDir, self.coarseOD, molELANOR, channel)
        if not os.path.exists(outDirABSCO): os.makedirs(outDirABSCO)

        # convert OD to ABSCO by scaling with density
        scale = denProf[iMatch[0]]
        print('Absorber scaling amout = %.3e' % scale)

        # do a bunch of stuff in create_fine_grid_xs.pro (which is 
        # called by process_lblrtm_xs.pro)
        print('Reading %s, %d of %d' % (odFile, iOD+1, nFilesOD) )

        # scale the ODs
        freq, od = lblTools.readOD(odFile, double=True)
        od /= scale

        # header info
        v1, v2, dv, nFreq = \
          freq.min(), freq.max(), freq[1]-freq[0], freq.size

        # for grid resampling: degrade resolution by a factor of 2 
        # in all channels except 2B1 (which is a factor of 4)
        if ('2B1' in odFile) or ('1B1' in odFile) or ('2A4' in odFile):
          kernel = np.array([1, 2, 3, 2, 1]) / 9.0
        else:
          kernel = np.array([1, 2, 1]) / 4.0
        # endif 2B1

        # convolve optical depth with weighting associated with kernel
        # then resample; most of this taken directly from IDL code
        # resample_grid.pro
        nDegrade = kernel.size - 1
        coarseRes = np.arange(0, nFreq, nDegrade)
        odCoarse = np.convolve(od, kernel)
        odCoarse[0] = od[0]; odCoarse[-1] = od[-1]
        odCoarse = np.array(odCoarse[coarseRes])
        freqCoarse = freq[coarseRes]

        # new header info
        v1, v2, dv, nFreq = freqCoarse.min(), freqCoarse.max(), \
          np.diff(freqCoarse)[0], freqCoarse.size

        # the AER convention of an 80-character header will not work
        #header = b'%80s' % self.headerOD
        header = np.zeros((178))
        header[143] = v1
        header[144] = v2
        header[142] = dv

        outFP = open('%s/%s' % (outDirABSCO, newName), 'wb')

        # ASSUMING LITTLE ENDIAN!
        # for some reason...
        wtf = '1234'
        outFP.write(wtf)
        outFP.write(header)

        # bytearray is Python 2 and 3 compatible
        # sticking with the create_fine_grid_xs.pro convention
        panelHeader = array.array('d', [v1, v2, dv])
        panelHeader.tofile(outFP)
        panelHeader = array.array('l', [nFreq, 0])
        panelHeader.tofile(outFP)
        panelHeader = array.array('d', [0])
        panelHeader.tofile(outFP)

        odOut = array.array('d', odCoarse)
        odOut.tofile(outFP)

        # more create_fine_grid_xs.pro convention
        panelHeader = array.array('d', [0])
        panelHeader.tofile(outFP)
        panelHeader = array.array('d', [0, 0, 0])
        panelHeader.tofile(outFP)
        panelHeader = array.array('l', [-99, -1])
        panelHeader.tofile(outFP)
        panelHeader = array.array('d', [0])
        panelHeader.tofile(outFP)

        outFP.close()
      # end odFile loop
    # end species loop
  # end postProcessXS()
# end xsABSCO()

class configSetup():
  def __init__(self, inFile):
    """
    Parse the input .ini file (inFile) and return as a dictionary for 
    use in the rest of this module

    Inputs
      inFile -- string, full path to .ini file that specifies paths 
        and filenames for ptn_file, tbar_file, tape5_dir, vmr_file, 
        densities_dir, lbl_path, tape3_path, xs_path, fscdxs, and
        lbl_run_dir. 

        there should also be channels and molecules sections, 
        which are handled a bit differently. the two sections are 
        their own dictionaries. channels should have 
        names, wn1, wn2, and res variables. molecules should have 
        lblxs and elanorxs variables. all variable names are case 
        sensitive

    Outputs
      outDict -- dictionary with keys for each of the paths specified 
        in inFile
    """

    # standard library, but name depends on Python version
    if sys.version_info.major < 3:
      import ConfigParser
    else:
      import configparser as ConfigParser
    # endif Python version

    cParse = ConfigParser.ConfigParser()
    cParse.read(inFile)
    cpSections = cParse.sections()

    # loop over each field (of all sections) and keep the field and 
    # associated value in returned object (self)
    for iCPS, cps in enumerate(cpSections):
      cItems = cParse.items(cps)
      if cps == 'channels':
        # make channels dictionary that is like CHANNELS global
        channels = {}
        for cItem in cItems: 
          if cItem[0] in ['wn1', 'wn2', 'res']:
            channels[cItem[0]] = \
              np.array(cItem[1].split()).astype(float)
          else:
            channels[cItem[0]] = cItem[1].split()
          # endif wn
        # end item loop

        setattr(self, 'channels', channels)
      elif cps == 'molecules':
        # make molecules dictionary like MOLECULES global
        molecules = {}
        for cItem in cItems: molecules[cItem[0]] = cItem[1].split()
        setattr(self, 'molecules', molecules)
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

  # end constructor
# end configSetup()

def asciiToCSV(inFile, outFile='VMR_data.csv'):
  """
  Pretty much a one-time script to convert the data blocks in LBLATM
  to a CSV that xsABSCO methods can use

  Inputs
    inFile -- string, full path to text snippet of LBLATM XS data 
      blocks (this should also include the ALTX block)

  Keywords
    outFile -- string, name of CSV file that will be written
  """

  # standard library
  import csv

  dat = open(inFile).read().splitlines()

  datBlock, allBlocks, names = [], [], []
  names.append('ALTX')
  for line in dat:
    # making lots of assumptions here...
    # comments typically represent the start of a new block
    # but there are 
    if '!' in line:
      if line.strip() != '!':
        # keep the XS name and append its data block into the complete
        # list of data blocks
        xsName = line.split()[4]

        # replace F11 and F12 LBLATM aliases with F11 and F12
        # the latter 2 of which are used in FSCDXS and xs/
        if xsName == 'CCL3F': xsName = 'F11'
        if xsName == 'CCL2F2': xsName = 'F12'
        names.append(xsName)

        # combine all lines of data block into a single string, 
        # omitting the metadata in the first and last lines and 
        # excluding the final comma in the data
        allBlocks.append(''.join(datBlock[1:-1])[:-1])
        datBlock = []
        continue
      else:
        continue
      # endif ! only
    # endif new data block

    # now some string manipulation to remove the non-data characters
    line = line.replace('&', '').strip()
    datBlock.append(''.join(line))
  # end line loop

  # process the last data block
  allBlocks.append(''.join(datBlock[1:-1])[:-1])

  # number of pressure levels
  nP = len(allBlocks[0].split(','))

  # now let's take the string list and create float arrays by 
  # utilizing the commas
  # blockArr will be a list of lists, each list being a data block
  blockArr = []
  for iBlock, block in enumerate(allBlocks):
    split = block.split(',')
    if len(split) == 1:
      # empty data blocks just have a "50*-99." string
      # let's replace that with a NaN for each pressure level
      blockArr.append([np.nan] * nP)
    else:
      # real data are separated by commas, so for each data block 
      # (i.e., species), we are creating a list of all data elements
      blockArr.append(block.split(','))
    # endif 
  # end block loop

  # let's make the array [nP x nMol] and true VMRs (not ppmv)
  blockArr = np.array(blockArr).astype(float)
  pressures = blockArr[0, :]
  vmr = blockArr[1:, :] * 1e-6
  blockArr = np.vstack((pressures, vmr)).T

  # write the CSV
  csvFP = open(outFile, 'w')
  csvWrite = csv.writer(csvFP)
  csvWrite.writerow(names)
  for iRow, row in enumerate(blockArr): csvWrite.writerow(row)
  csvFP.close()

  return True
# end asciiToCSV()

if __name__ == '__main__':

  parser = argparse.ArgumentParser(\
    description='Generate ABSCO tables for cross section molecules.')
  parser.add_argument('--config_file', type=str, \
    default='xsABSCO_config.ini', \
    help='Configuration file that contains the file and ' + \
    'directory names necessary for the xsABSCO class.')

  # argument switches (booleans, no assignments)
  parser.add_argument('-csv', '--make_csv', action='store_true', \
    help='Starting from an ASCII text file with LBLATM altitude ' + \
    'and VMR data blocks (assuming it exists, is usable in ' + \
    'asciiToCSV(), and is named "VMR_temp.txt"), condense file ' + \
    'into a more computer-friendly CSV that is used in other ' + \
    'portions of this module. This probably only needs to be ' + \
    'done once at most.')
  parser.add_argument('-den', '--calc_densities', \
    action='store_true', \
    help='Calculate densities from VMR for each XS molecule.')
  parser.add_argument('-t5', '--make_tape5', action='store_true', \
    help='Make the TAPE5 files for each species, channel, ' + \
    'pressure, and temperature.')
  parser.add_argument('-lbl', '--run_lbl', action='store_true', \
    help='Run LBLRTM for each of the TAPE5s generated and saved ' + \
    'in makeTAPE5().')
  parser.add_argument('-re', '--regrid', action='store_true', \
    help='Run the post-processing code that scales the ODs and ' + \
    'writes new, smaller ODint files for each specified species.')
  parser.add_argument('-e2e', '--end_to_end', action='store_true', \
    help='Runs the entire process from tape 5 generation to ' + \
    'post-processing (rather than entering all of the keywords ' + \
    'separately).')
  args = parser.parse_args()

  iniFile = args.config_file; utils.file_check(iniFile)
  ini = configSetup(iniFile)

  csvFile = ini.vmr_file
  if args.make_csv:
    status = asciiToCSV('VMR_temp.txt', outFile=csvFile)
    print('%s was written to working directory' % csvFile)
  # end CSV generation

  # instantiation; there are lots of keywords here so that we have 
  # some flexibility (i did not wanna force the user to have to use 
  # a configSetup object with the xsABSCO class)
  xsObj = xsABSCO(ptnFile=ini.ptn_file, tBarFile=ini.tbar_file, \
    dirT5=ini.tape5_dir, vmr=csvFile, dirDen=ini.densities_dir, \
    channels=ini.channels, molecules=ini.molecules, \
    pathLBL=ini.lbl_path, pathT3=ini.tape3_path, \
    pathXS=ini.xs_path, pathFSCDXS=ini.fscdxs, \
    dirWork=ini.lbl_run_dir, dirFine=ini.od_dir, \
    dirCoarse=ini.absco_dir, odHeader=ini.header)

  if args.calc_densities: xsObj.calcDensity()
  if args.make_tape5: xsObj.makeTAPE5()
  if args.run_lbl: xsObj.runLBL()
  if args.regrid: xsObj.postProcessXS()

  # haven't tested this yet, but no reason it won't work...right?
  if args.end_to_end:
    xsObj.makeTAPE5()
    xsObj.runLBL()
    xsObj.postProcessXS()
  # endif e2e
# end main()

