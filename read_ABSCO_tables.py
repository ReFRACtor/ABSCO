#!/usr/bin/env python

import os, sys, argparse
import numpy as np
import xarray as xa
import netCDF4 as nc
import time

sys.path.append('common')
import utils

class testABSCO():
  def __init__(self, inArgs):
    """
    - Time the difference between netCDF4 and xarray libraries

    - find pressure, temperature, spectral point, and h2o (if 
      necessary) corresponding to user specifications

    - print absorption coefficient at user-specified array coordinate
    """

    utils.file_check(inArgs['ncFile'])
    self.ncFile = str(inArgs['ncFile'])
    self.userP = float(inArgs['in_pressure'])
    self.userT = float(inArgs['in_temp'])
    self.userH2O = float(inArgs['in_h2o'])
    self.userO2 = float(inArgs['in_o2'])
    self.molName = os.path.basename(self.ncFile).split('_')[0]
    self.h2o = True if self.molName in ['O2', 'CO2', 'N2'] else False
    self.o2 = True if self.molName == 'O2' else False

    # ironically, this cannot be frequency
    freq = float(inArgs['in_spectral'][0])
    units = inArgs['in_spectral'][1]
    if units not in ['cm-1', 'um', 'nm']:
      sys.exit('Please provide valid spectral unit [cm-1, um, nm]')

    if units == 'um': wnConvert = 1e4 
    if units == 'nm': wnConvert = 1e7
    self.userWN = wnConvert / freq if units in ['um', 'nm'] else \
      float(freq)

    # tolerance for float equality check in valueLocate() method
    tol = inArgs['tolerance']
    self.tol = 1e-5 if tol is None else float(tol)
  # end constructor

  def valueLocate(self):
    """
    Find array indices that correspond to the user-provided 
    coordinates (P, T, wavenumber, WV VMR), then make sure the values
    from the netCDF are within the user-provided tolerance of the 
    user-provided value
    """

    with xa.open_dataset(self.ncFile) as xaObj:
      ncP = np.array(xaObj.variables['P_level'])
      ncT = np.array(xaObj.variables['Temperature'])
      ncWN = np.array(xaObj.variables['Spectral_Grid'])
      if self.h2o: ncH2O = np.array(xaObj['H2O_VMR'])
      if self.o2: ncO2 = np.array(xaObj['O2_VMR'])
    # endwith

    # find closest values for each coordinate
    idxP = np.nanargmin(np.abs(ncP-self.userP))
    idxT = np.nanargmin(np.abs(ncT[idxP]-self.userT))
    idxWN = np.nanargmin(np.abs(ncWN-self.userWN))

    # pressure input is on levels, but ABSCOs are on layers, so we 
    # need to make sure we don't try to index a layer with a level 
    # index
    if idxP == ncP.size-1:
      sys.exit('Please specify a level P that is not at the TOA')
    
    if self.h2o: idxH2O = np.nanargmin(np.abs(ncH2O-self.userH2O))
    if self.o2: idxO2 = np.nanargmin(np.abs(ncO2-self.userO2))

    # are the closest values close to what the user wants?
    pClose = np.isclose(self.userP, ncP[idxP], rtol=self.tol)
    tClose = np.isclose(self.userT, ncT[idxP, idxT], rtol=self.tol)
    wnClose = np.isclose(self.userWN, ncWN[idxWN], rtol=self.tol)
    closeList = [pClose, tClose, wnClose]
    paramList = ['P', 'T', 'Wavenumber']
    valList = [ncP[idxP], ncT[idxP, idxT], ncWN[idxWN]]

    if self.h2o:
      h2oClose = np.isclose(self.userH2O, ncH2O[idxH2O],rtol=self.tol)
      closeList.append(h2oClose)
      paramList.append('H2O')
      valList.append(ncH2O[idxH2O])
    # endif h2o

    if self.o2:
      o2Close = np.isclose(self.userO2, ncO2[idxO2],rtol=self.tol)
      closeList.append(o2Close)
      paramList.append('O2')
      valList.append(ncO2[idxO2])
    # endif o2

    for iParam, close in enumerate(closeList):
      if not close:
        errMsg = '%s not within %f%% of ' % \
          (paramList[iParam], self.tol*100)
        errMsg += 'closest corresponding value in '
        errMsg += '%s (%f), returning' % \
          (self.ncFile, valList[iParam])
        sys.exit(errMsg)
      # endif close
    # end close loop

    # save for later usage
    self.iP = int(idxP)
    self.iT = int(idxT)
    self.iWN  = int(idxWN)

    if self.h2o: self.iH2O = int(idxH2O)
    if self.o2: self.iO2 = int(idxO2)
  # end valueLocate()

  def readABSCO(self):
    """
    Read in the netCDF (and time the differences in the load by the 
    netCDF4 and xarray libraries)
    """

    # really just get 8% back relative to netCDF4 library, so not a 
    # huge efficiency improvement
    print('Reading %s' % self.ncFile)
    with xa.open_dataset(self.ncFile) as xaObj:
      absco = np.array(xaObj.variables['Cross_Section'])

    if self.o2:
      coord = (self.iWN, self.iT, self.iP, self.iH2O, self.iO2)
      out = absco[self.iWN, self.iT, self.iP, self.iH2O, self.iO2]
    elif self.h2o:
      coord = (self.iWN, self.iT, self.iP, self.iH2O)
      out = absco[self.iWN, self.iT, self.iP, self.iH2O]
    else:
      coord = (self.iWN, self.iT, self.iP)
      out = absco[self.iWN, self.iT, self.iP]
    # endif h2o

    print("Cross_Section indices: %s" % (coord, ))
    print("Cross_Section value: %.6E" % out)

  # end readABSCO
# end readABSCO()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(\
    formatter_class=argparse.ArgumentDefaultsHelpFormatter, \
    description='Read netCDF generated with run_LBLRTM_ABSCO.py ' + \
    'module and print out absorption coefficient (k) at a given ' + \
    'pressure, temperature, spectral point, and water vapor ' + \
    'amount (if molecule continuum is affected by water vapor).')
  parser.add_argument('ncFile', type=str, \
    help='Output netCDF generated by run_LBLRTM_ABSCO.py.')
  parser.add_argument('-p', '--in_pressure', type=float, \
    default=1050.0, \
    help='Reference pressure level [mbar] for which k is ' + \
    'retrieved. There are two pressure boundaries for a given ' + \
    'layer, and this is the lower bound (ie, closest to surface).')
  parser.add_argument('-T', '--in_temp', type=float, \
    default=230.0, \
    help='Reference temperature [K] for which k is retrieved.')
  parser.add_argument('-s', '--in_spectral', nargs=2, \
    default=[500, 'cm-1'], \
    help='Reference spectral point AND units [cm-1, um, or nm] ' + \
    'for which k is retrieved.')
  parser.add_argument('-wv', '--in_h2o', type=float, \
    default=10.0, \
    help='Reference water vapor VMR (ppmv) for which k is ' + \
    'retrieved IF the specified molecule is H2O, CO2, O2, or N2.')
  parser.add_argument('-o2', '--in_o2', type=float, \
    default=190000.0, \
    help='Reference oxygen VMR (ppmv) for which k is ' + \
    'retrieved IF the specified molecule is O2.')
  parser.add_argument('-tol', '--tolerance', type=float, \
    help='Tolerance used when searching for floating point ' + \
    'matches in *each* of the dimensions. This should be a ' + \
    'relative tolerance (e.g. 0.01 would mean P from netCDF is ' + \
    'within 1%% of in_pressure).')
  args = parser.parse_args()

  test = testABSCO(vars(args))
  test.valueLocate()
  test.readABSCO()

# end main()
