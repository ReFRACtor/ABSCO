#!/usr/bin/env python

import os, sys, glob, argparse
import subprocess as sub
from distutils.spawn import find_executable as which

import numpy as np

sys.path.append('common')
import utils

class submodules():
  def __init__(self, inArgs, lnfl=False, lbl=False, lines=False):
    """
    Build models as needed (one at a time) and replace paths in 
    ABSCO_tables.py configuration file if specified
    """

    errOne = 'Please specify one of lnfl, lbl, or lines'
    modBool = np.array([lnfl, lbl, lines])
    if np.all(modBool): sys.exit(errOne)
    if np.where(modBool)[0].size > 1: sys.exit(errOne)

    # preprocessing
    compiler = inArgs['compiler']
    compiler = compiler.lower()
    if compiler not in ['ifort', 'gfortran', 'pgf90']:
      sys.exit('%s is not a supported compiler.' % compiler)

    iniFile = args.config_file
    if iniFile is not None: utils.file_check(iniFile)

    if lnfl:
      lnflDir = args.lnfl_path; utils.file_check(lnflDir)

    if lbl:
      lblDir = args.lblrtm_path; utils.file_check(lblDir)

    if lines:
      linesDir = args.lines_path; utils.file_check(linesDir)

    self.compiler = str(compiler)
    self.iniFile = None if iniFile is None else str(iniFile)
    self.lines = False

    # LNFL: always single precision, LBLRTM: always double
    if lnfl:
      self.modelDir = str(lnflDir)
      self.modelStr = 'LNFL'
      self.doLNFL = True
      self.precision = 'sgl'
      self.makeStr = 'make_lnfl'
      self.pathStr = 'lnfl_path'
    elif lbl:
      self.modelDir = str(lblDir)
      self.modelStr = 'LBLRTM'
      self.doLBL = True
      self.precision = 'dbl'
      self.makeStr = 'make_lblrtm'
      self.pathStr = 'lbl_path'
    elif lines:
      self.modelDir = str(linesDir)
      self.pathStr = ['tape1_path', 'tape2_path', 'extra_params', \
        'xs_path', 'fscdxs']
      self.lines = True
    else:
      sys.exit('No model build chosen')
    # endif model
  # end constructor

  def build(self):
    """
    Build LBLRTM or LNFL
    """

    # OS determination
    # https://docs.python.org/2/library/sys.html#sys.platform
    compPath = which(self.compiler)
    platform = sys.platform
    if platform in ['linux', 'linux2']:
      self.opSys = 'linux'
    elif platform == 'darwin':
      self.opSys = 'osx'
    elif platform == 'win32':
      self.opSys = 'mingw'
    else:
      sys.exit('Could not determine OS, returning')
    # endif OS

    # compiler string generation
    if self.compiler == 'ifort':
      self.compStr = 'INTEL'
    elif self.compiler == 'gfortran':
      self.compStr = 'GNU'
    elif self.compiler == 'pgf90':
      self.compStr = 'PGI'
    # endif compiler

    cmd = '%s%s%s' % (self.opSys, self.compStr, self.precision)

    cwd = os.getcwd()

    os.chdir('%s/build' % self.modelDir)
    print('Building %s' % self.modelStr)
    status = sub.call(['make', '-f', self.makeStr, cmd])
    if status != 0: sys.exit('%s not built' % self.modelStr)
    os.chdir(cwd)

    return self
  # end build()

  def configFile(self):
    """
    Replace paths in ABSCO_tables.py configuration file with the paths
    established in this class
    """

    if self.iniFile is None:
      sys.exit('No configuration file specified, returning')
    else:
      iniDat = open(self.iniFile).read().splitlines()
      outFP = open(self.iniFile, 'w')
      print('Replacing %s in %s' % (self.pathStr, self.iniFile))

      if self.lines:
        # making some assumptions about directory structure here...
        modStr = ['line_file/aer_v_3.6', 'line_file/lncpl_lines', \
          'extra_brd_params', \
          'xs_files_v3.6/xs', 'xs_files_v3.6/FSCDXS']
      else:
        # make_lnfl and make_lblrtm -o arguments
        modStr = '%s/%s_*_%s_%s_%s' % \
          (self.modelDir, self.modelStr.lower(), \
           self.opSys, self.compStr.lower(), self.precision)
        modExe = glob.glob(modStr)[0]
      # endif lines

      for line in iniDat:
        if self.lines:
          for old, new in zip(self.pathStr, modStr):
            if (old in line):
              split = line.split('=')
              line = line.replace(split[1], ' %s' % new)
            # endif LNFL
          # end path loop
        else:
          if (self.pathStr in line):
            split = line.split('=')
            line = line.replace(split[1], ' %s' % modExe)
          # endif LNFL
        # endif lines
        outFP.write('%s\n' % line)
      # end dat loop
      outFP.close()
    # endif ini
  # end configFile()
# end submodules class

if __name__ == '__main__':
  parser = argparse.ArgumentParser(\
    description='Build LBLRTM and LNFL executables for usage in ' + \
    'ABSCO_config.ini and ABSCO_tables.py.')
  parser.add_argument('-c', '--compiler', default='ifort', \
    help='Name of compiler with which user intends to build ' + \
    '([ifort, gfortran, pgf90] are supported). Case-insensitive')
  parser.add_argument('-ini', '--config_file', \
    help='Name of configuration file in which to add exe paths.')
  parser.add_argument('-lnfl', '--lnfl_path', default='LNFL', \
    help='Path of LNFL submodule directory (top level).')
  parser.add_argument('-lbl', '--lblrtm_path', default='LBLRTM', \
    help='Path of LBLRTM submodule directory (top level).')
  parser.add_argument('-lines', '--lines_path', \
    default='AER_Line_File', help='Top-level path of AER line file.')
  parser.add_argument('--only_lines', action='store_true', \
    help='If set, only the line file paths are changed in the ' + \
    'configuration file and no model builds are done.')
  args = parser.parse_args()

  # first replace line file paths
  subObj = submodules(vars(args), lines=True)
  subObj.configFile()
  if args.only_lines:
    sys.exit('Only replaced lines paths in %s' % args.config)

  # now do LNFL -- line file paths will not change
  subObj = submodules(vars(args), lnfl=True)
  subObj.build()
  subObj.configFile()

  # now do LBL -- LNFL build and ini path will not change
  subObj = submodules(vars(args), lbl=True)
  subObj.build()
  subObj.configFile()

# end main()
