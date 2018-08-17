#!/usr/bin/env python

import os, sys, glob, argparse
import subprocess as sub
from distutils.spawn import find_executable as which

sys.path.append('common')
import utils

class submodules():
  def __init__(self, inArgs, lnfl=False, lbl=False):
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

    self.compiler = str(compiler)
    self.iniFile = str(iniFile)

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
    else:
      sys.exit('No model build chosen')
    # endif model
  # end constructor

  def build(self):
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
    if self.iniFile is not None:
      iniDat = open(self.iniFile).read().splitlines()
      outFP = open(self.iniFile, 'w')
      print('Replacing %s in %s' % (self.pathStr, self.iniFile))

      # make_lnfl and make_lblrtm -o arguments
      modStr = '%s/%s_*_%s_%s_%s' % \
        (self.modelDir, self.modelStr.lower(), \
         self.opSys, self.compStr.lower(), self.precision)
      modExe = glob.glob(modStr)[0]

      for line in iniDat:
        if (self.pathStr in line):
          split = line.split('=')
          line = line.replace(split[1], ' %s' % modExe)
        # endif LNFL
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
  args = parser.parse_args()

  subObj = submodules(vars(args), lnfl=True)
  subObj.build()
  subObj.configFile()

  # now do LBL -- LNFL build and ini path will not change
  subObj = submodules(vars(args), lbl=True)
  subObj.build()
  subObj.configFile()

# end main()
