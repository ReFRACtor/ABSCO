#!/usr/bin/env python

import os, sys, argparse

# conda modules
import numpy as np

# Git submodules
sys.path.append(os.path.join(os.path.dirname(__file__), 'common'))
import utils

# local modules
import ABSCO_preprocess as preproc
import ABSCO_compute as absco

parser = argparse.ArgumentParser(\
  description='Generate ABSCO tables for user-specified molecules.',
  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--config_file', type=str, \
  default='ABSCO_config.ini', \
  help='Configuration file that contains the file and ' + \
  'directory names necessary for the makeABSCO class.')

# argument switches (booleans, no assignments)
parser.add_argument('-lnfl', '--run_lnfl', action='store_true', \
  help='Run LBLRTM for each of the TAPE5s generated and saved ' + \
  'in lblT5().')
parser.add_argument('-lbl', '--run_lbl', action='store_true', \
  help='Run LBLRTM for each of the TAPE5s generated and saved ' + \
  'in lblT5().')
parser.add_argument('-e2e', '--end_to_end', action='store_true', \
  help='Runs the entire process from tape 5 generation to ' + \
  'post-processing (rather than entering all of the keywords ' + \
  'separately).')
parser.add_argument('-db', '--debug', action='store_true', \
  help='Use for testing (only iterates over a few pressure levels.)')
parser.add_argument('-y', '--no_prompt', action='store_false', dest='prompt_user', default=True,
  help='Do not prompt for continuation at warning messages.')
args = parser.parse_args()

if not (args.run_lbl or args.run_lnfl or args.end_to_end):
  sys.exit('Nothing done. Please set "-e2e", "-lbl", or "-lnfl"')

iniFile = args.config_file; utils.file_check(iniFile)

# configuration object instantiation
ini = preproc.configure(iniFile, prompt_user=args.prompt_user)

if args.run_lnfl or args.end_to_end:
  # don't need to save these objects because runLNFL() will do all 
  # we need (i.e., TAPE3s for LBL run); also VMR is not needed yet
  for mol in ini.molnames: 
    kObj = absco.makeABSCO(ini, mol, vmrWV=np.nan, vmrO2=np.nan)
    kObj.lnflT5(mol)
    kObj.runLNFL()
  # end mol loop
# end LNFL

if args.run_lbl or args.end_to_end:
  for mol in ini.molnames:
    if mol in ini.molH2O:

      # we have to handle water vapor-effected molecules a little
      # differently (they will have an extra dimension in their 
      # output ABSCO arrays)
  
      if mol == 'O2':
        # oxygen is even weirder, because we have to implement its 
        # own vmr dimension
        # we do this by just repeating what we do with the objects
        # for H2O, but running with each O2/H2O VMR pair
        o2ObjList = []
        for o2ppm in ini.vmrO2:
          wvObjList = []
          for wvppm in ini.wv_vmr:
            kObj = absco.makeABSCO(ini, mol, debug=args.debug, \
              vmrWV=wvppm, vmrO2=o2ppm)
            kObj.lblT5(mol)
            kObj.calcABSCO(mol)
            kObj.arrABSCO()
            wvObjList.append(kObj)
          # end H2O VMR loop

          o2ObjList.append(absco.combineVMR(wvObjList))

        # end O2 VMR loop

        kObj = absco.combineVMR(o2ObjList)
      else:
        vmrObjList = []
        for ppm in ini.wv_vmr:
          kObj = absco.makeABSCO(ini,mol,debug=args.debug,vmrWV=ppm)
          kObj.lblT5(mol)
          kObj.calcABSCO(mol)
          kObj.arrABSCO()
          vmrObjList.append(kObj)
        # end VMR loop

        kObj = absco.combineVMR(vmrObjList)
      # endif 

      kObj.makeNC(mol)
    else:
      kObj = absco.makeABSCO(ini, mol, debug=args.debug)
      kObj.lblT5(mol)
      kObj.calcABSCO(mol)
      kObj.arrABSCO()
      kObj.makeNC(mol)
    # endif H2O
  # end mol loop

# end LBL

