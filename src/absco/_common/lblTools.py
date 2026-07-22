#!/usr/bin/python

# for Python 3 compatibility
from __future__ import print_function

# -*- coding: utf-8 -*-
import os
import sys
import shutil
import re
import time
import math
import types
import glob
import struct
import numpy as np
import stat
from . import FortranFile
import traceback
import logging

logger = logging.getLogger(__name__)

#from configTools import xmlConfig

OPTICAL_DEPTH = 'opticalDepth'
RADIANCE = 'radiance'
TRANSMISSION = 'transmission'
NUMBER_DENSITY = 'numberDensity'
LAYER_OPTICAL_DEPTHS = 'layerOpticalDepths'
SOLAR = 'solar'
WAVE_NUMBER = 'waveNumber'
SOLAR_OPTICAL_DEPTH = 'solarOD'
SOLAR_RADIANCE = 'solarRadiance'

def interP(p, pin, vin, debug=False, rh=None):
  """
  Interpolate pressure grid

  Input
    p -- float array
    pin -- float array
    vin -- float arr

  Output
    z -- float array
  """

  z = []
  for p1 in p:
    mmin = [-999999, p1]
    mmax = [999999, p1]
    cont = 0
    for k in range(len(pin)):
        if p1 != pin[k]:
            gg = pin[k] - p1
            if gg > 0 and gg < mmax[0]:
                mmax = [gg, pin[k], k]
            elif gg < 0 and gg > mmin[0]:
                mmin = [gg, pin[k], k]
        # endif p1 !=

        if p1 == pin[k]:
            z.append(vin[k])
            cont = 1
            break
        # endif p1 ==

    if cont: continue

    if mmax[1] == p1:
      mmax = [-999999, p1]
      for k in range(len(pin)):
        if p1 != pin[k] and k != mmin[2]:
          gg = pin[k] - p1
          if gg < 0 and gg > mmax[0]: mmax = [gg, pin[k], k]
        # endif p1 and k
    # endif mmax[1]

    if mmin[1] == p1: mmin = [999999, p1]

    for k in range(len(pin)):
      if p1 != pin[k] and k != mmax[2]:
        gg = pin[k] - p1
        if gg > 0 and gg < mmin[0]: mmin = [gg, pin[k], k]
      # end if p1 and k
    # endif mmin[1]

    z1 = math.log(p1) - math.log(mmin[1])
    z2 = math.log(mmax[1]) - math.log(mmin[1])
    z3 = vin[mmax[2]] - vin[mmin[2]]
    z4 = vin[mmin[2]]
    hh = z4 + z1 / z2 * z3

    if (rh is not None) and (p1 < 100) and (hh > 0.003): hh = 0.001

    z.append(hh)
  # end p loop

  return z
# end interP()

def readOD(path, double=False):
  """
  Read in binary LBLRTM ODint files (IMRG=1 and IOD=1, 3, or 4 in
  LBLRTM specifications (Record 1.2 in lblrtm_instructions.html)

  Call
    (ff, od, parms) = readOD(path)

  Input
    path -- string, full path to directory with ODint files

  Output
    ff -- float array, wavenumbers of spectrum
    od -- float array, optical depths of spectrum
    parms -- float array, nLayers x nMolecules, molecular amounts 
      (cm-2) as a funtion of molecule and layer

  Keywords
    double -- boolean, read in double precision OD values
  """
  ff, od = readTape12(path, double=double)

  # this assumes a directory of OD files is provided
  # we will make the function handle one level (i.e., one OD file)
  # as specified by the user in "path"
  """
  # grab list of OD files in path
  ll = glob.glob(os.path.join(path, 'ODint_*'))
  od = None
  ll.sort()

  for i in range(len(ll)):
    # read in binary data like a TAPE12, store into OD array
    (ff, s) = readTape12(ll[i], double=double)
    if od is None:
      od=np.vstack((np.asarray(s),))
    else:
      od=np.vstack((od,np.asarray(s),))
  # end loop over ODint files
  """

  # not all LBLRTM runs have IATM=1 such that LBLATM is called for 
  # density computation
  """
  f = open(os.path.join(path, 'TAPE7'))
  lines = f.readlines()
  f.close()
  z = re.split('\s+', lines[1].strip())
  n = int(z[1])

  # initialize list that will eventually be nMol x nLayers
  # (molecules per cm-2)
  # list of lists, where each sublist contains nMol elements
  parms = []

  # read in ASCII TAPE7 parameters (calculations of molecular amounts)
  # loop over layers (this may be assuming a certain number of 
  # molecules, or at least the number should not exceed a threshold 
  # of 7)
  for i in range(2,len(lines),2):
    p = []
    z = re.split('\s+', lines[i].strip())
    p.append(float(z[0]))

    if len(z)>5:
      a0 = float(z[3])
      p0 = float(z[4])
    else:
      zz = re.split('\.', z[3])
      a0 = float(zz[0] + '.' + zz[1])
      p0 = float('.' + zz[2])
    # endif len z

    if len(z)>7:
      a1 = float(z[6])
      dz = a1 - a0
      a2 = a1
      p1 = float(z[7])
      dp = p0 - p1
      p2 = p1
    elif len(z)>3:
      dz = a0 - a2
      a2 = a0
      dp = p2 - p0
      p2 = p0
    else:
      dz = 0
      dp = 0
    # endif len z

    p.append(dz)
    p.append(dp)

    z = re.split('\s+', lines[i + 1].strip())

    # for i in z:p.append(float(i))
    # concatenatate (NOT append) floating point list of z onto 
    # existing p (mol amounts for all all molecules in a given layer)
    p += map(float, z)
    parms.append(p)
  # end loop over TAPE7
  """

  #return (np.array(ff), np.array(od), np.array(parms))
  return (np.array(ff), np.array(od))
# end readOD()

def readTape7(fName, sList=False):
  """
  Read LBLRTM TAPE7 ASCII file (molecular amounts computed by LBLATM)

  Call
    ll = readTape7(fName)

  Input
    fName -- string, full path to TAPE7 to be read

  Output
    ll -- 

  Keywords
    sList -- 
  """

  def cnvline(l, fmt=','):
    ll = re.split(fmt, l.strip())
    ll = map(float, ll)
    return ll
  # end cnvline()

  f = open(fName)
  l = f.readlines()
  f.close()
  ll = []
  l = l[2:]

  # loop over layers (this may be assuming a certain number of 
  # molecules, or at least the number should not exceed a threshold 
  # of 7)
  for i in range(0, len(l), 2):
    q = re.split('\s+', l[i].strip())
    if not i:
      try:
        if sList:
          q = [float(q[0]), float(q[1]), float(q[6]),
               float(q[6]) - float(q[3])]
        else:
          q = [[float(q[0]), float(q[4]), float(q[7])],
               [float(q[1]), float(q[5]), float(q[8])],
               float(q[6]), float(q[6]) - float(q[3])]
        # endif sList
      except (KeyboardInterrupt, SystemExit):
        raise
      except (IndexError,ValueError):
        q = [float(q[0]), float(q[1]), 0, 0]
    else:
      ind = len(ll) - 1
      if sList:
        q2 = ll[ind][2]
        try:
          q = [float(q[0]), float(q[1]), float(q[3]),
               float(q[3]) - q2]
        except (IndexError,ValueError):
          zz = re.split('\.', q[3])
          q3 = float(zz[0] + zz[1])
          q = [float(q[0]), float(q[1]), q3, q3 - q2]
      else:
        q0 = ll[ind][0][2]
        q1 = ll[ind][1][2]
        q2 = ll[ind][2]
        q = [[float(q[0]), q0, float(q[4])], [float(q[1]), q1,
             float(q[5])], float(q[3]), float(q[3]) - q2]
      # endif sList

    q += cnvline(l[i + 1].strip(), fmt='\s+')
    ll.append(q)
  # end loop over lines

  return ll
# end readTape7

def readTape12(fileName, double=False, fType=0):
  """
  Read in a TAPE12 (or similar TAPE10-13) LBLRTM output binary file
  and return spectrum

  Original: Scott Zaccheo (AER), modified by Rick Pernak (AER, 2017) 

  Call
    waveNumbers, output = readTape12(fileName)

  Input
    fileName -- string, full path to binary file

  Output
    waveNumbers -- float array, wavenumbers of spectrum
    output -- float array, corresponding parameter (radiance, 
      transmittance, optical depth -- see lblrtm_instructions.html 
      file assignments)

  Keywords
    double -- boolean, read in double-precision bits
    fType -- int, file type; from 
      /project/rc/rc2/mshep/idl/patbrown/read_lbl_file.pro:
    
      0: scanned transmittance or radiance, optical depth
      1: radiance from monochromatic radiance calculation
      2: transmittance from monochromatic radiance calculation
      3: Aerosol absorbtance from spectral aerosol transmittance file
      4: Aerosol scattering from spectral aerosol transmittance file
      5: Aerosol asymmetry from spectral aerosol transmittance file
      6: transmittance from monochromatic radiance calculation

      File types 3, 4, and 5 are from TAPE20 files and have not been 
      tested with this function
  """

  def readLBLPanel(ffObj, pHeaderForm):
    """
    Read in a single panel

    Input
      ffObj -- FortranFile object
      pHeaderForm -- string, format of panel header (e.g., 'dddl' for
        3 doubles and a long integer, which would mean wn_start and 
        wn_end are doubles, the spectral resolution is double 
        precision, and the number of points in the panel is a long 
        integer; wn_start and wn_end are always doubles)

    Output
      outWN -- float array, wavenumbers of spectrum for a given panel
      outParam -- float array, corresponding parameter (radiance, 
      transmittance, optical depth -- see lblrtm_instructions.html 
      file assignments) of panel

    Keywords
    """

    # initialization of variables that change with panel
    OK = True;
    outWN = []; outParam = []
    ctr = 1
    while OK:
      buff = fortranFile.getRecord()

      # grab binary data if valid buffer and not a panel header
      #if len(buff) == struct.calcsize(lfmt): continue
      while buff is not None and len(buff) != struct.calcsize(lfmt):
        buff = fortranFile.getRecord()

      # end while buff
      if buff:
        try:
          # read panel header and underlying data
          (v1, v2, dv, nPanel) = struct.unpack(lfmt, buff)
          data = fortranFile.readDoubleVector() if double else \
            fortranFile.readFloatVector()
          #print len(data), ctr
          ctr += 1

          # concatenate (NOT append) onto output
          # for now, this is just radiance -- looks like other
          # parameters like transmittances are in other panels
          outParam += data

          # concatenate wavenumber array (without numpy) based 
          # on spectral resolution and starting wavenumber
          outWN += \
            map(lambda x:v1 + x * dv, range(len(data)))
        except (KeyboardInterrupt, SystemExit):
          raise
        except:
          #print 'Data could not be read, file may be corrupted'
          OK = False
      else:
        # end of panel
        #continue
        OK = False
      # endif buff
    # end while OK
    return outWN, outParam
  # end readLBLPanel()

  # main readTape12()
  iFormat = 'l'
  if struct.calcsize('l') == 8:iFormat = 'i'

  # instantiate FortranFile object
  fortranFile = FortranFile.FortranFile(fileName)

  # read file header
  data = fortranFile.getRecord()

  # format for each panel header
  if double: lfmt = 'dddl'
  else: lfmt = 'ddf%s' % iFormat

  waveNumbers, output = readLBLPanel(fortranFile, lfmt)

  # read file header
  #waveNumbers, output = readLBLPanel(fortranFile, lfmt, header=False)
  #print len(waveNumbers)
  #waveNumbers, output = readLBLPanel(fortranFile, lfmt)

  return np.array(waveNumbers), np.array(output)
# end readTape12()

def rpReadTape12(fileName, double=False, fType=0):
  """
  Read in a TAPE12 (or similar TAPE10-13) LBLRTM output binary file
  and return spectrum

  Original: Scott Zaccheo (AER), modified by Rick Pernak (AER, 2017) 

  Call
    waveNumbers, output = rpReadTape12(fileName)

  Input
    fileName -- string, full path to binary file

  Output
    waveNumbers -- float array, wavenumbers of spectrum
    output -- float array, corresponding parameter (radiance, 
      transmittance, optical depth -- see lblrtm_instructions.html 
      file assignments)

  Keywords
    double -- boolean, read in double-precision bits
    fType -- int, file type; from 
      /project/rc/rc2/mshep/idl/patbrown/read_lbl_file.pro:
    
      0: scanned transmittance or radiance, optical depth
      1: radiance from monochromatic radiance calculation
      2: transmittance from monochromatic radiance calculation
      3: Aerosol absorbtance from spectral aerosol transmittance file
      4: Aerosol scattering from spectral aerosol transmittance file
      5: Aerosol asymmetry from spectral aerosol transmittance file
      6: transmittance from monochromatic radiance calculation

      File types 3, 4, and 5 are from TAPE20 files and have not been 
      tested with this function
  """

  from . import utils

  # main readTape12()
  iFormat = 'l'
  if struct.calcsize('l') == 8:iFormat = 'i'

  # instantiate FortranFile object
  fortranFile = FortranFile.FortranFile(fileName)

  # read file header
  data = fortranFile.getRecord()

  # format for each panel header
  lfmt = 'dddl' if double else 'ddf%s' % iFormat
  headLen = struct.calcsize(lfmt)

  # the number of output parameters differs by file type,
  # how many are expected?
  if fType in [0, 1]:
    paramExp = 1
  elif fType in range(2, 6):
    paramExp = 2
  else:
    paramExp = 3
  # endif fType

  paramStr = ['Scanned', 'Radiance', 'Transmittance']

  while True:
    # skip to first panel header
    buff = fortranFile.getRecord()
    #while buff is None: buff = fortranFile.getRecord()
    #while buff is not None and len(buff) != headLen:
    #  buff = fortranFile.getRecord()

    if buff:
      if buff and len(buff) == headLen:
        (v1, v2, dv, nPtPanel) = struct.unpack(lfmt, buff)
        print(v1, v2, dv)
      else:
        data = fortranFile.readDoubleVector() if double else \
          fortranFile.readFloatVector()
      # endif buf and len
    else:
      data = fortranFile.readDoubleVector() if double else \
        fortranFile.readFloatVector()
      print(len(data))
    # endif buff len

    break
    if data is None: break

    # concatenate (NOT append) onto output
    # for now, this is just radiance -- looks like other
    # parameters like transmittances are in other panels
    outParam += data

    # concatenate wavenumber array (without numpy) based 
    # on spectral resolution and starting wavenumber
    outWN += \
      map(lambda x:v1 + x * dv, range(len(data)))

    #print ctr, v1, v2, len(outParam)
  # endwhile

  sys.exit('GOT HERE')
  if buff:
    try:
      # read panel header and underlying data
      
      data = fortranFile.readDoubleVector() if double else \
        fortranFile.readFloatVector()

      # concatenate (NOT append) onto output
      # for now, this is just radiance -- looks like other
      # parameters like transmittances are in other panels
      outParam += data

      # concatenate wavenumber array (without numpy) based 
      # on spectral resolution and starting wavenumber
      outWN += \
        map(lambda x:v1 + x * dv, range(len(data)))
    except (KeyboardInterrupt, SystemExit):
      raise
    except:
      #print 'Data could not be read, file may be corrupted'
      OK = False
  else:
    # end of panel
    #continue
    OK = False
  # endif buff
# end while OK

  waveNumbers, output = np.array(outWN), np.array(outParam)

  return waveNumbers, output
# end rpReadTape12()

def getOD(tape5, ostream=sys.stdout):
  """
  Read in OD from...TAPE5!?!

  Input
    tape5 -- string, full path to TAPE5 file

  Output
    None
  """

  path = os.path.dirname(tape5)
  return readOD(path, fout=ostream)
# end getOD

def removeFileName(fileName):
  """
  Remove given file

  Input
    fileName -- string, full path to file to be removed
  Output
    None
  """

  try: os.remove(fileName)
  except OSError: pass
# end removeFileName

def generatePressureGrid(pressure, observer, target, nLayers):
  try:
    try:
      mn = int(min(pressure))
      mx = int(max(pressure))
    except (ValueError,IndexError,TypeError):
      mn = int(min([observer, target]))
      mx = int(max([observer, target]))
    # end try

    dx = (mx - mn) / (nLayers - 1)

    layers = map(lambda x: mx - x * dx, range(nLayers))
    return -len(layers), layers
  except (KeyboardInterrupt, SystemExit):
    raise
  except:
    return 0, 0
  # end try
# end generatePressureGrid

def generateHeightGrid(heights, observer, target, nLayers):
    # input heights assumed to be in meters
    try:
        dx = int(max(heights) - min(heights)) / (nLayers - 1)
        layers = map(lambda x: x * dx / 1000. + min(heights) / 1000., range(nLayers))
        return len(layers), layers
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        try:
            hRange = max([observer, target]) - min([observer, target])
            dx = hRange / (nLayers - 1)
            layers = map(lambda x: x * dx + min([observer, target]), range(nLayers))
            return len(layers), layers
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            return 0, 0

def writeTape5(path, parameterDictionary, isFile=False, monoRTM=False,useMeters=False):    
    DEFAULT_NMOL=7
    DEFAULT_CO2=330.
    DEFAULT_CH4=1.7

    if isFile:
        outputFileName = path
    else:
        if monoRTM: baseName = 'MONORTM.IN'
        else: baseName = 'TAPE5'

        outputFileName = os.path.join(path, baseName)
        
        try: 
            os.makedirs(path)
        except OSError:
            pass

    tape5file = open(outputFileName, 'w')

    waveNumber1 = parameterDictionary['v1']
    waveNumber2 = parameterDictionary['v2']
    deltaWaveNumber = parameterDictionary['dv']

    # Define Tangent

    if parameterDictionary.has_key('tangentFlag') and parameterDictionary['tangentFlag']:
        tangent=1
    else:
        tangent=0

    # Define IAERSL

    if parameterDictionary.has_key('aerosols'):
        aerosols = parameterDictionary['aerosols']
    else: 
        aerosols = 0

    # Define IEMIT

    if parameterDictionary.has_key('output'):
        outputType = int(parameterDictionary['output'])
    else: outputType = 0

    if parameterDictionary.has_key('model'):
        model = int(parameterDictionary['model'])
    else: model = None

    if parameterDictionary.has_key('units'):
        units = int(parameterDictionary['units'])
    else: units = None

    if not parameterDictionary.has_key('usePressure'):
        observer = parameterDictionary['h1'] / 1000.
        target = parameterDictionary['h2'] / 1000.
        units = 1
    elif not parameterDictionary['usePressure']:
        observer = parameterDictionary['h1'] / 1000.
        target = parameterDictionary['h2'] / 1000.
        units = 1
    else:
        observer = parameterDictionary['h1']
        target = parameterDictionary['h2']
        units = -1

    if parameterDictionary.has_key('pathL'): pl = parameterDictionary['pathL'] / 1000.
    elif parameterDictionary.has_key('pathLength'): pl = parameterDictionary['pathLength'] / 1000.
    else: pl = 0

    if parameterDictionary.has_key('horz'): 
        if parameterDictionary['horz']:
            horz=1
        else:
            horz=0
    else: horz = 0

    angle = parameterDictionary['angle']

    # define user defined layers/levels
    
    userAltitudes = []
    userPressures = []
    
    if parameterDictionary.has_key('userDefinedLevels'):
        parameterDictionary['udl'] = parameterDictionary['userDefinedLevels']

    if parameterDictionary.has_key('udl'):
        if parameterDictionary['udl']:
            if type(parameterDictionary['udl']) == types.ListType:
                if parameterDictionary.has_key('usePressure') and parameterDictionary['usePressure']:
                    userPressures = sorted(parameterDictionary['udl'],
                            None, None, 1)
                    userDefinedLayers = -len(parameterDictionary['udl'])
                else:
                    userAltitudes = sorted(parameterDictionary['udl'])
                    userDefinedLayers = len(parameterDictionary['udl'])
            else:
                if parameterDictionary['udl'] is True: userDefinedLayers = 100
                elif parameterDictionary['udl'] > 10: userDefinedLayers = parameterDictionary['udl']
                else: userDefinedLayers = 100

                if parameterDictionary.has_key('Height'): height = parameterDictionary['Height']
                else: height = None
                    
                if parameterDictionary.has_key('usePressure'):
                    if parameterDictionary['usePressure']:
                        userDefinedLayers, userPressures = generatePressureGrid(parameterDictionary['Pres'],
                                                                             observer, target,
                                                                             userDefinedLayers)
                    else:                            
                        userDefinedLayers, userAltitudes = generateHeightGrid(height,
                                                                           observer, target, userDefinedLayers)
                else:
                    userDefinedLayers, userAltitudes = generateHeightGrid(height,
                                                                       observer, target, userDefinedLayers)
        else: userDefinedLayers = 0
    else: userDefinedLayers = 0

    if parameterDictionary.has_key('surfaceTerrain'):
        surfaceTerrain = parameterDictionary['surfaceTerrain']
    else: surfaceTerrain = [300, 0.1, 0, 0, 0.9, 0, 0]

    # write record 1.1
    print >> tape5file, '$ TAPE5 by Python, range %f %f %s' % (waveNumber1, waveNumber2, time.asctime())

    if not parameterDictionary.has_key('inFlag'): parameterDictionary['inFlag'] = 0
    if not parameterDictionary.has_key('iotFlag'): parameterDictionary['iotFlag'] = 0

    # solar upwelling
    if outputType == 2 and not monoRTM:
        if parameterDictionary['inFlag'] == 2 and parameterDictionary['iotFlag'] == 2:
            print >> tape5file, ' HI=0 F4=0 CN=0 AE=0 EM=2 SC=0 FI=0 PL=0 TS=0 AM=0 MG=0 LA=0 OD=0 XS=0    0    0'
            print >> tape5file, "%5d%5d  %3d" % (parameterDictionary['inFlag'], parameterDictionary['iotFlag'],
                                              parameterDictionary['solarDay'])
            print >> tape5file, "0.0      0.0"
            print >> tape5file, "-1."
            print >> tape5file, "%"
            tape5file.close()
            return os.path.join(path, 'TAPE5')

    if parameterDictionary.has_key('iodFlag'):
        iodFlag = parameterDictionary['iodFlag']
    else:
        if outputType: iodFlag = 1
        else: iodFlag = 1

    if parameterDictionary.has_key('noContinuum'):
        if parameterDictionary['noContinuum']:
            continuumFlag = 0
        else:
            continuumFlag = 1
    else:
        continuumFlag = 1
            
    # write record 1.2
    if monoRTM:
        iPlot = 1
        iod = 0
        print >> tape5file, "%4s%1i%9s%1i%9s%1i%14s%1i%9s%1i%14s%1i%4s%1i%16s%4i" % ("", 1, "", 1, "", 1, "", iPlot, "", 1, "", iod, "", 0, "", 0)        
    else:
        if outputType < 2:
            print >> tape5file, \
            ' HI=1 F4=1 CN=%0d AE=%0d EM=%0d SC=0 FI=0 PL=0 TS=0 AM=1 MG=0 LA=0 OD=%0d XS=0    0    0' \
              % (continuumFlag,aerosols, outputType, iodFlag)
        else:
            print >> tape5file, \
            ' HI=0 F4=0 CN=0 AE=%0d EM=%0d SC=0 FI=0 PL=0 TS=0 AM=1 MG=0 LA=0 OD=%0d XS=0    0    0' \
              % (aerosols, outputType, iodFlag)

        # write record 1.2a

        if outputType == 2:
            # INFLAG,  IOTFLG,  JULDAT
            #  1-5,    6-10,   13-15
            #   I5,      I5,  2X, I3

            print >> tape5file, "%5d%5d  %3d" % (parameterDictionary['inFlag'], parameterDictionary['iotFlag'],
                                              parameterDictionary['solarDay'])

    # determine molecule scaling
        
    nms = 0
    if parameterDictionary.has_key('co2scale'):
        if parameterDictionary['co2scale']: nms = 6
    if parameterDictionary.has_key('wvScale'):
        if parameterDictionary['wvScale']: nms = 6
    if parameterDictionary.has_key('ch4scale'):
        if parameterDictionary['ch4scale']: nms = 6
            
    if nms>0:
        co2scale = DEFAULT_CO2
        if parameterDictionary.has_key('co2scale'): co2scale = parameterDictionary['co2scale']
        if not co2scale: co2scale = DEFAULT_CO2

        wvScale = 1.0
        if parameterDictionary.has_key('wvScale'): wvScale = parameterDictionary['wvScale']
        if not wvScale: wvScale = 1.0

        ch4scale = DEFAULT_CH4
        if parameterDictionary.has_key('ch4scale'): ch4scale = parameterDictionary['ch4scale']
        if not ch4scale: ch4scale = DEFAULT_CH4

    co2only=False
    o2only=False
    
    if parameterDictionary.has_key('co2only'):
        if parameterDictionary['co2only']:
            nms=DEFAULT_NMOL
            co2only=True
    elif parameterDictionary.has_key('o2only'):
        if parameterDictionary['o2only']:
            nms=DEFAULT_NMOL
            o2only=True
            
    # write record 1.3

    if monoRTM:
        print >> tape5file, '%10.3f%10.3f%10s%10.3e%63s%2i' % (waveNumber1,
                                                       waveNumber2, '', deltaWaveNumber, '', nms)
    else:
        if iodFlag:
            print >> tape5file, '%10.3f%10.3f%70s%10.3e  %2i' % (waveNumber1,
                                                         waveNumber2, '', deltaWaveNumber, nms)
        else:
            print >> tape5file, '%10.3f%10.3f%70s%10.3e  %2i' % (waveNumber1,
                                                         waveNumber2, '', 0, nms)
        
    if nms>0:
        print >> tape5file, '11   1'
        if o2only:
            stringFormat='%15d%15d%15d%15d%15d%15d%15d'
            scaleFactors=(0,0,0,0,0,0,1)
        elif co2only:
            stringFormat='%15d%15d%15d%15d%15d%15d%15d'
            scaleFactors=(0,1,0,0,0,0,0)
        else:
            stringFormat='%15.7e%15.7e%15.7e%15.7e%15.7e%15.7e'
            scaleFactors=(wvScale,co2scale/float(DEFAULT_CO2),1,1,1,ch4scale/float(DEFAULT_CH4))
 
        # write record 1.3a
        print >> tape5file,stringFormat % scaleFactors
        
    # write record 1.4
    if monoRTM:
        print >> tape5file, '%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%5s' % tuple(surfaceTerrain + [''])
    else:
        if outputType == 1:
            if parameterDictionary['angle'] > 90 and parameterDictionary['angle'] <= 180: surfRefl = ['l']
            else: surfRefl = ['s']
            print >> tape5file, '%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%5s' % tuple(surfaceTerrain + surfRefl)

    if observer == target:
        heightType = 1
        userDefinedLayers = 0
    else:
        if tangent:
            heightType = 3
        else:
            heightType = 2

    # write record 3.1

    iFXTYPE = 0
    iMunits = 0
    Re = ''
    hSpace = ''
    vBar = ''

    if parameterDictionary.has_key('refLatitude'): 
        refLat = "%10.3f" % float(parameterDictionary['refLatitude'])
    else: 
        refLat = ''
    
    if horz:
        print >> tape5file, '    %1d    1    0    1    0    7    1%2i %2i%10s%10s%10s%10s%10s' % \
            (model, iFXTYPE, iMunits, Re, hSpace, vBar, '', refLat)
        print >> tape5file, '    0.000                    %10.3f' % pl
        userDefinedLayers = 0
    else:
        print >> tape5file, '    %1d    %1d%5d    1    0    7    1%2i %2i%10s%10s%10s%10s%10s' % (model, heightType, 
                                                                                                  userDefinedLayers, 
                                                                                                  iFXTYPE,
                                                                                                  iMunits, Re, 
                                                                                                  hSpace, vBar, '', 
                                                                                                  refLat)
        print >> tape5file, '%10.3f%10.3f%10.3f%10s%5i' % (observer, target, angle,'',tangent)

    if not model:
        # aptg is altitude, pressure, temperature and gases vector
        
        aptg = parameterDictionary['aptg']

        # build table for the first 6 gases and set default to US standard
        
        for i in range(len(aptg), 9): aptg.append(6)
        if len(aptg) > 9: aptg = aptg[:9]
        
        if parameterDictionary.has_key('Height'): 
            a = parameterDictionary['Height']
        else: 
            a = None

        if parameterDictionary.has_key('Pres'): 
            p = parameterDictionary['Pres']
        else: 
            p = None
            
        if parameterDictionary.has_key('Temp'): 
            t = parameterDictionary['Temp']
        else: 
            t = None

        if parameterDictionary.has_key('WV'): 
            w = parameterDictionary['WV']
        else: 
            w = None

        if aptg[3] > '9': co2 = parameterDictionary['CO2']
        else: co2 = None

        if aptg[4] > '9': o3 = parameterDictionary['O3']
        else: o3 = None
            
        if aptg[5] > '9': n2o = parameterDictionary['N2O']
        else: n2o = None
            
        if aptg[6] > '9': co = parameterDictionary['C0']
        else: co = None
            
        if aptg[7] > '9': ch4 = parameterDictionary['CH4']
        else: ch4 = None

        if aptg[8] > '9': o2 = parameterDictionary['O2']
        else: o2 = None

        if units > 0:
            b = sorted(a)
            z = map(a.index, b)
        else:
            b = sorted(p, reverse=1)
            z = map(p.index, b)

        if a: a = map(lambda x: a[x] / 1000, z)
        if p: p = map(lambda x: p[x], z)

        if t is not None: t = map(lambda x: t[x], z)
        if w is not None: w = map(lambda x: w[x], z)
        if co2 is not None: co2 = map(lambda x: co2[x], z)
        if o3 is not None: o3 = map(lambda x: o3[x], z)
        if n2o is not None: n2o = map(lambda x: n2o[x], z)
        if co is not None: co = map(lambda x: co[x], z)
        if ch4 is not None: ch4 = map(lambda x: ch4[x], z)
        if o2 is not None: o2 = map(lambda x: o2[x], z)

        if userDefinedLayers:
            if userDefinedLayers > 0:
                for i in range(0, len(userAltitudes), 8):
                    print >> tape5file, '',
                    for j in range(8):
                        try: print >> tape5file, '%9.3f' % userAltitudes[i + j],
                        except IndexError: break
                    print >> tape5file
            else:
                for i in range(0, len(userPressures), 8):
                    print >> tape5file, '',
                    for j in range(8):
                        try: 
                            print >> tape5file, '%9.3f' % userPressures[i + j],
                        except IndexError: 
                            break
                    print >> tape5file
        else:
            if not horz:
                print >> tape5file

        if horz:
            print >> tape5file, \
                '    1                 Input from python application max h=%dm' \
                % observer
        else:
            print >> tape5file, \
                '%5d               Input from python application max h=%dm' \
                % (units * len(p), observer)

        fmt = '%10.3f%10.3f%10.3f     %s%s   %s%s%s%s%s%s%s'

        fmt2 = ''

        def addFormatString(variable):
            if variable is not None: return '%10.3e'
            return '          '

        def addVariable(inputVector, variable, index):
            if inputVector is None: inputVector = []
            if variable is not None: inputVector.append(variable[index])
            return inputVector

        fmt2 += addFormatString(w)
        fmt2 += addFormatString(co2)
        fmt2 += addFormatString(o3)
        fmt2 += addFormatString(n2o)
        fmt2 += addFormatString(co)
        fmt2 += addFormatString(ch4)
        fmt2 += addFormatString(o2)
 
        for i in range(len(p)):
            if a: alt = a[i]
            else: alt = 0

            z = [alt, p[i], t[i]]

            # write user defined first row
            for j in aptg: z.append(j)
            print >> tape5file, fmt % tuple(z)

            z = addVariable(None, w, i)
            z = addVariable(z, co2, i)
            z = addVariable(z, o3, i)
            z = addVariable(z, n2o, i)
            z = addVariable(z, co, i)
            z = addVariable(z, ch4, i)
            z = addVariable(z, o2, i)

            # write user defined second row
            print >> tape5file, fmt2 % tuple(z)
    else:
        if userDefinedLayers:
            if userDefinedLayers > 0:
                for i in range(0, len(userAltitudes), 8):
                    print >> tape5file, '',
                    for j in range(8):
                        try:
                            print >> tape5file, '%9.3f' % userAltitudes[i + j],
                        except IndexError:
                            break
                    print >> tape5file
            else:
                for i in range(0, len(userPressures), 8):
                    print >> tape5file, '',
                    for j in range(8):
                        try:
                            print >> tape5file, '%9.3f' % userPressures[i + j],
                        except IndexError:
                            break
                    print >> tape5file
        else:
            if not horz:
                print >> tape5file

    if aerosols and not monoRTM:

        # set record 4.1 IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT
        # Format (6I5,5F10.3)

        iHaze = aerosols
        iSeason = 0
        gndAlt = 0

        aerosolString = '%5d%5d%5d%5d%5d%5d' % (iHaze, iSeason, 0, 0, 0, 0,)
        aerosolString += '%10s%10.3f%10.3f%10.3f%10.3f' % ('', 0.0, 0.0, 0.0, gndAlt)

        # aerosolString+="%10.3f%10.3f%10.3f%10.3f%10.3f"%(vis,0.0,0.0,0.0,gndAlt)

        print >> tape5file, aerosolString

    print >> tape5file, '%'

    tape5file.close()
    return outputFileName

def readARM(fileName):
    nc = pupynere.netcdf_file(fileName)
    h = (nc.variables['alt'])[:]
    h = h - h[0]
    p = (nc.variables['pres'])[:]
    t = (nc.variables['tdry'])[:] + 273
    td = (nc.variables['dp'])[:] + 273
    rh = (nc.variables['rh'])[:]
    uWinds = (nc.variables['u_wind'])[:]
    vWinds = (nc.variables['v_wind'])[:]

    tt = nc.variables['base_time']

    nc.close()
    return (
        h,
        p,
        t,
        td,
        rh,
        uWinds,
        vWinds,
        tt.data,
        )
"""
class LblObject:
    outputList = {OPTICAL_DEPTH: 0, TRANSMISSION: 0, 'Radiance': 1, RADIANCE:1, TRANSMISSION:0, SOLAR:2}
    lblSolarFileName = 'SOLAR.RAD'
    reflFileName = 'SOL.REFLECTANCE'

    modelValues = {
        'User Defined': 0,
        'Tropical': 1,
        'Midlat Summer': 2,
        'Midlat Winter': 3,
        'Arctic Summer': 4,
        'Arctic Winter': 5,
        'U.S. Standard': 6,
         }

    def getModelName(self, number):
        for key in self.modelValues:
            if self.modelValues[key] == int(number): return key
        return None
    
    def getBaseRtParameters(self, dv=0.0001, others=dict()):
        data = {'dv': dv,
              # layers parameters and flags
              'udl': 80,
              'usePressure': False,
              'horz': False,
              # surface parameter
              'surfaceTerrain':[300, 0.8, 0, 0, 0.2, 0, 0],
              # atmospheric parameters
              'aptg': ['A', 'A', 'H', '6', '6', '6', '6', '6', '6'],
              'co2scale': 380.0,
              'Temp': [], 'CO2': [], 'WV': [], 'N2O': [], 'O3': [], 'O2': [], 'CH4': [], 'CO': [],
              'aerosol': 0,
              # solar parameter and flags
              'inFlag':0,
              'iodFlag':0,
              'h2':0,
              'h1':100000                
              }

        for key in others: data[key] = others[key]
        return data
    
    def __init__(self, workPath=None, rtCommand=None, tape3fileName=None, solarFileName=None, debug=True, log=sys.stdout):
        self.workPath = workPath
        self.rtCommand = rtCommand
        self.tape3fileName = tape3fileName
        self.solarFileName = solarFileName
        self.debug = debug
        self.log = log

        self.monoRtmCommand = None
        self.monoRtmSpectraFileName = None

        self.makeWorkPath()
        
    def write(self, buffer):
        print >> self.log, buffer,

    def makeWorkPath(self):
        if self.workPath:
            try:
                os.makedirs(self.workPath)
                print >> self, "creating %s for lblrtm scratch files ..." % (self.workPath)
            except OSError:
                # traceback.print_exc()
                print >> self, "using %s for lblrtm scratch files ..." % (self.workPath)
    
    def loadConfiguration(self, configFileName, section):
        configData = xmlConfig.xmlConfig(configFileName)
        lblrtmConfig = configData.getDictionary(section)['lblrtm'][0]
        self.workPath = configData.parse(lblrtmConfig['lblPath'])
        self.tape3fileName = configData.parse(lblrtmConfig['lblTape3'])
            
        try: 
            self.solarFileName = configData.parse(lblrtmConfig['lblSolarFile'])
        except (KeyboardInterrupt, SystemExit):
            raise
        except: 
            self.solarFileName = None

        try: 
            self.rtCommand = configData.parse(lblrtmConfig['lblCommand'])
        except (KeyboardInterrupt, SystemExit):
            raise
        except: 
            self.rtCommand = configData.parse('$(lblCommand)')

        try: 
            self.monoRtmCommand = configData.parse(lblrtmConfig['monoRtmCommand'])
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            try: self.monoRtmCommand = configData.parse('$(monoRtmCommand)')
            except: self.monoRtmCommand = None

        try: 
            self.monoRtmSpectraFileName = configData.parse(lblrtmConfig['monoRtmSpectraFile'])
        except (KeyboardInterrupt, SystemExit):
            raise
        except: 
            self.monoRtmSpectraFileName = None

        if self.debug:
            print >> self, "using tape3 file %s ..." % (self.tape3fileName)
            print >> self, "using solar file %s ..." % (self.solarFileName)
            print >> self, "using rt command %s ..." % (self.rtCommand)
            print >> self, "using monort command %s ..." % (self.monoRtmCommand)
            print >> self, "using monort spectra file %s ..." % (self.monoRtmSpectraFileName)

        self.makeWorkPath()

    def setLog(self, fileObject):
        self.log = fileObject
        
    def genReflectance(self, wn, coeff):
        return coeff[0] + coeff[1] * wn + coeff[2] * wn ** 2

    def writeReflectance(self, v1, v2, dv, reflCoeff):
    
        maxWN = 24000
        maxBlock = 2400
        nFillHeader = 264
        guardWN = 4

        fortranFile = FortranFile.FortranFile(os.path.join(self.workPath, self.reflFileName), write=True)

        v1 = float(int(v1))
        v2 = float(int(v2))
        
        nWN = int(round((v2 - v1 + 2 * guardWN) / dv) + 1)
        if nWN > maxWN:
            dv = float(v2 - v1 + 2 * guardWN) / (maxWN - 1)
            nWN = int(round((v2 - v1 + 2 * guardWN) / dv) + 1)
 
        startWN = v1 - guardWN
        endWN = startWN + (nWN - 1) * dv
 
        # create formatted header
        
        data = [map(lambda x:0, range(108)), [dv, startWN, endWN, 0, 0, 0], [1, 1], [0, 0, dv, 0, 0], [0, 1, 1, 0]]
        formatString = ['d', 'd', 'i', 'd', 'i', 'd']
        fortranFile.writeFormatVector(data, formatString, 1056)

        # write reflectances 

        for startIndex in range(0, nWN, maxBlock):
            startWN = startWN + startIndex * dv
            endWN = startWN + (maxBlock - 1) * dv

            if endWN > v2 + guardWN: endWN = v2 + guardWN
            nBlock = int((endWN - startWN) / dv + 0.5) + 1
            endWN = startWN + (nBlock - 1) * dv

            fortranFile.writeFormatVector((startWN, endWN, dv, nBlock), 'ddfi')
            refl = map(lambda x:self.genReflectance(startWN + x * dv, reflCoeff), range(nBlock))
            fortranFile.writeFloatVector(refl)
           
        fortranFile.writeFormatVector((endWN + dv, endWN, dv, -99), 'ddfi')
        fortranFile.close()
    
    def run(self, rtParameters,
            opticalDepthFlag=True, radianceFlag=False,
            upwellingFlag=False, downwellingFlag=False,
            layerOpticalDepthFlag=False,
            numberDensityFlag=False,
            thread=True):

        tangentFlag=False
        
        output = dict()

        outputGrid = map(lambda x:x * rtParameters['dv'] + rtParameters['v1'],
                       range(int((rtParameters['v2'] - rtParameters['v1']) / rtParameters['dv'] + 0.5)))
        
        if rtParameters.has_key('opticalDepthFlag'): opticalDepthFlag = rtParameters['opticalDepthFlag']
        if rtParameters.has_key('radianceFlag'): radianceFlag = rtParameters['radianceFlag']
        if rtParameters.has_key('downwellingFlag'): downwellingFlag = rtParameters['downwellingFlag']
        if rtParameters.has_key('upwellingFlag'): upwellingFlag = rtParameters['upwellingFlag']
        if rtParameters.has_key('numberDensityFlag'): numberDensityFlag = rtParameters['numberDensityFlag']
        if rtParameters.has_key('tangentFlag'): tangentFlag = rtParameters['tangentFlag']
        if rtParameters.has_key('layerOpticalDepthFlag'): layerOpticalDepthFlag = rtParameters['layerOpticalDepthFlag']

        if upwellingFlag or downwellingFlag:
            # setup solar file
            try: os.remove(os.path.join(self.workPath, self.lblSolarFileName))
            except OSError: pass

            try: os.symlink(self.solarFileName, os.path.join(self.workPath, self.lblSolarFileName))
            except OSError: shutil.copy(self.solarFileName, os.path.join(self.workPath, self.lblSolarFileName))

        if opticalDepthFlag:
            rtParameters['output'] = self.outputList[OPTICAL_DEPTH]
            tape5 = writeTape5(self.workPath, rtParameters)

            results = run(self.rtCommand, self.tape3fileName, tape5, ostream=self, 
                          readNumberDensity=numberDensityFlag,
                          readLayerOpticalDepths=layerOpticalDepthFlag,
                          thread=thread)

            if results.has_key(OPTICAL_DEPTH):
                if not (downwellingFlag or upwellingFlag): output[WAVE_NUMBER] = outputGrid[:]
                output[OPTICAL_DEPTH] = np.interp(outputGrid, results[OPTICAL_DEPTH][0], results[OPTICAL_DEPTH][1]).tolist()

                if results.has_key(NUMBER_DENSITY):
                    output[NUMBER_DENSITY] = results[NUMBER_DENSITY][:]
                if results.has_key(LAYER_OPTICAL_DEPTHS):
                    output[LAYER_OPTICAL_DEPTHS] = results[LAYER_OPTICAL_DEPTHS][:]
        
        if radianceFlag or downwellingFlag or upwellingFlag:
            rtParameters['output'] = self.outputList[RADIANCE]
            tape5 = writeTape5(self.workPath, rtParameters)
            results = run(self.rtCommand, self.tape3fileName, tape5, ostream=self, thread=thread)
            
            if results.has_key(RADIANCE):
                if radianceFlag: output[WAVE_NUMBER] = outputGrid[:]
                output[RADIANCE] = (np.interp(outputGrid, results[RADIANCE][0], results[RADIANCE][1]) * 1e4).tolist()

        if downwellingFlag or upwellingFlag:
            # lblrtm run number 1
            # Step 1
            # Calculate the downward transmittance and radiance 
            # Make sure TBOUND,SREMISS and SRREFL are 0.0 and
            # H1 is the surface, H2 is TOA and the angle is the solar zenith angle.
            # Note that downward radiance from this calculation is not used.

            localRtParameters = rtParameters.copy()

            if tangentFlag:
                localRtParameters['h2'] = 0
            else:
                localRtParameters['h2'] = 100000
                localRtParameters['h1'] = min([rtParameters['h1'], rtParameters['h2']])
                
            localRtParameters['angle'] = rtParameters['solarZenithAngle']
            localRtParameters['surfaceTerrain'] = [0, 0, 0, 0, 0, 0, 0]
             
            localRtParameters['output'] = self.outputList[OPTICAL_DEPTH]
            tape5 = writeTape5(self.workPath, localRtParameters)
            results = run(self.rtCommand, self.tape3fileName, tape5, ostream=self, thread=thread)
            if results.has_key(OPTICAL_DEPTH):
                output[SOLAR_OPTICAL_DEPTH] = np.interp(outputGrid, results[OPTICAL_DEPTH][0],
                                                      results[OPTICAL_DEPTH][1]).tolist()

            localRtParameters['output'] = self.outputList[RADIANCE]
            tape5 = writeTape5(self.workPath, localRtParameters)
            wn, radiance = run(self.rtCommand, self.tape3fileName, tape5, ostream=self, thread=thread)
            if self.debug: self.copyTapeFile('TAPE5', 'TAPE5.downwelling')

            # run lblrtm a second time (could be combined with first run lbl call)
            
            localRtParameters['output'] = self.outputList[SOLAR]
            localRtParameters['inFlag'] = 0
            localRtParameters['iotFlag'] = 1
            tape5 = writeTape5(self.workPath, localRtParameters)
            results = run(self.rtCommand, self.tape3fileName, tape5, ostream=self, clean=False, readSolar=True,
                        thread=thread)
            
            if results.has_key(SOLAR):
                if not upwellingFlag: output[WAVE_NUMBER] = outputGrid[:]
                # convert from 1/cm2 to 1/m2
                # output[RADIANCE]=(np.interp(outputGrid,wn,radiance)*1e4).tolist()
                
                '''
                The solar code in LBLRTM was originally written to model the direct beam
                reaching the surface and generate output to be compared with radiance
                measuring instruments. Thus the code converts the irradiance to radiance
                by dividing by 6.8-5 (the subangle subtended by the Sun at the Earth). 
                It also does the /m2 to /cm2 conversion.

                In order to get correct reflected radiances you need to consider the
                total incoming radiation and therefore you need to multiply the LBLRTM
                output by 6.8e-5. Furthermore, LBLRTM treats the reflectivity provided
                as correct for the geometry determined by the angles used in the
                downwelling (solar zenith angle) and upwelling (viewer angle) runs. 
                '''

                output[SOLAR_RADIANCE] = (np.interp(outputGrid, results[SOLAR][0], results[SOLAR][1]) * 1e4)
                output[SOLAR_RADIANCE] *= 6.8e-5
                output[SOLAR_RADIANCE] = output[SOLAR_RADIANCE].tolist()
                
        if upwellingFlag:
            localRtParameters = rtParameters.copy()

            if rtParameters.has_key('surfaceTerrain'):   
                surfaceParameters = rtParameters['surfaceTerrain']
            else: 
                surfaceParameters = [300, 0.8, 0, 0, 0.2 / math.pi, 0, 0]

            # write solar reflectance file
                
            self.writeReflectance(rtParameters['v1'], rtParameters['v2'], rtParameters['dv'], surfaceParameters[4:])

            # temporary method
                        
            reflectance = self.genReflectance(np.asarray(outputGrid), surfaceParameters[4:])

            try:
                output[SOLAR_RADIANCE] = np.asarray(output[SOLAR_RADIANCE]) * np.exp(-np.asarray(output[OPTICAL_DEPTH]))
                output[SOLAR_RADIANCE] *= reflectance
                output[SOLAR_RADIANCE] += np.asarray(output[RADIANCE])
                output[SOLAR_RADIANCE] = output[SOLAR_RADIANCE].tolist()
                output[WAVE_NUMBER] = outputGrid[:]
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                pass
                
            ''' original method 

            # lblrtm run number 2
            # Step 2
            # Calculate the upward transmittance and radiance 
            # Make sure H1, H2, ANGLE, TBOUND,SREMISS and SRREFL are set correctly.
            # H1 is the observer height, H2 is the surface and ANGLE is the observer angle (180 is nadir)
            # The TOA radiance obtained from this calculation will be added to the solar transmitted reflected term from step 3.

            localRtParameters['output']= self.outputList[RADIANCE]
            
            tape5=writeTape5(self.workPath,localRtParameters)
            self.copyTapeFile('TAPE12','SOL.PATH.T2')

            wn,radiance=run(self.rtCommand,self.tape3fileName,tape5,clean=False,ostream=self,
                            thread=thread)
            if self.debug: self.copyTapeFile('TAPE5','TAPE5.upwelling')
            
            if wn is None:
                print >> self,'LBLRTM Calculate solar upwelling step #1 failed'
                return output

            self.copyTapeFile('TAPE12','TAPE12.upwelling')

            # lblrtm run number 3
            # Step 3
            # Calculate transmission and reflection of solar radiance
            # Make sure the reflectvity used in Step 2 is consistent with the SOL.REFLECTANCE file needed here. 
            # Final output (radiance at the observer) is in binary file TAPE13.

            localRtParameters['output']= self.outputList[SOLAR]
            localRtParameters['inFlag']= 2
            localRtParameters['iotFlag']= 2
 
            tape5=writeTape5(self.workPath,localRtParameters)

            wn,radiance,solarRadWN,solarRadiance=run(self.rtCommand,self.tape3fileName,tape5,ostream=self,
                                                     clean=False,readSolar=True,
                                                     thread=thread)            
            if self.debug: self.copyTapeFile('TAPE5','TAPE5.solar')
            
            if wn:
                if len(wn) !=len(output['waveNumber']):
                    output['waveNumber']=outputGrid[:]
                    output[SOLAR_RADIANCE]=np.interp(outputGrid,solarRadWN,solarRadiance)*6.8e-5*1e4
                    output[SOLAR_RADIANCE]=output[SOLAR_RADIANCE].tolist()
            else:
                print >> self,'LBLRTM Calculate solar upwelling step #3 failed'
            
            '''

        if not self.debug:
            try: os.remove(os.path.join(self.workPath, self.reflFileName))
            except OSError: pass
            try: os.remove(os.path.join(self.workPath, self.lblSolarFileName))
            except OSError: pass

        return output

    def copyTapeFile(self, inputFile, outputFile):
        shutil.copy(os.path.join(self.workPath, inputFile), os.path.join(self.workPath, outputFile))

    def cleanWorkSpace(self):
        try: shutil.rmtree(self.workPath)
        except OSError: pass
        try: print >> self, "removing %s scratch files and path ..." % (self.workPath)
        except IOError: print "removing %s scratch files and path ..." % (self.workPath)

class MonoRtmObject(LblObject):

    def __init__(self, workPath=None, rtCommand=None, spectraFileName=None, debug=True, log=sys.stdout):
        LblObject.__init__(self, workPath=workPath, rtCommand=rtCommand, tape3fileName=spectraFileName,
                           solarFileName=None, debug=debug, log=log)
 
        self.rtCommand = rtCommand
        self.tape3fileName = spectraFileName
    
    def loadConfiguration(self, configFileName, section):
        configData = xmlConfig.xmlConfig(configFileName)
        rtmConfig = configData.getDictionary(section)['monortm'][0]
        self.workPath = configData.parse(rtmConfig['rtPath'])

        try: 
            self.rtCommand = configData.parse(rtmConfig['monoRtmCommand'])
        except (IndexError,KeyError,RuntimeError): 
            self.rtCommand = configData.parse('$(monoRtmCommand)')

        self.tape3fileName = configData.parse(rtmConfig['monoRtmSpectraFile'])

        self.makeWorkPath()
    
    def run(self, rtParameters, thread=True):
        output = None
        rtParameters['output'] = self.outputList[RADIANCE]
        tape5 = writeTape5(self.workPath, rtParameters, monoRTM=True)
        results = run(self.rtCommand, self.tape3fileName, tape5, \
          ostream=self, thread=thread, monoRTM=True, cwd=self.workPath)

        if results.has_key(RADIANCE):
            output = dict()
            output[WAVE_NUMBER] = results[RADIANCE][0]
            output[RADIANCE] = results[RADIANCE][1]
            output[OPTICAL_DEPTH] = results[OPTICAL_DEPTH][1]
       
        return output
"""
