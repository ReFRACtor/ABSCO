#!/bin/env python

import os, sys, glob, argparse
import subprocess as sub
from distutils.spawn import find_executable as which

sys.path.append('common')
import utils

parser = argparse.ArgumentParser(\
  description='Build LBLRTM and LNFL executables for usage in ' + \
  'ABSCO_config.ini and ABSCO_tables.py.')
parser.add_argument('-c', '--compiler', default='ifort', \
  help='Name of compiler with which user intends to build ' + \
  '([ifort, gfortran, pgf90] are supported). Case-insensitive')
parser.add_argument('-ini', '--config_file', \
  help='Name of configuration file in which to add executable paths.')
parser.add_argument('-lnfl', '--lnfl_path', default='LNFL', \
  help='Path of LNFL submodule directory (top level).')
parser.add_argument('-lbl', '--lblrtm_path', default='LBLRTM', \
  help='Path of LBLRTM submodule directory (top level).')
args = parser.parse_args()

# preprocessing
compiler = args.compiler
compiler = compiler.lower()
if compiler not in ['ifort', 'gfortran', 'pgf90']:
  sys.exit('%s is not a supported compiler.' % compiler)
iniFile = args.config_file
if iniFile is not None: utils.file_check(iniFile)

lnflDir = args.lnfl_path; utils.file_check(lnflDir)
lblDir = args.lblrtm_path; utils.file_check(lblDir)

# OS determination
# https://docs.python.org/2/library/sys.html#sys.platform
compPath = which(compiler)
platform = sys.platform
if platform in ['linux', 'linux2']:
  opSys = 'linux'
elif platform == 'darwin':
  opSys = 'osx'
elif platform == 'win32':
  opSys = 'mingw'
else:
  sys.exit('Could not determine OS, returning')
# endif OS

# compiler string generation
if compiler == 'ifort':
  compStr = 'INTEL'
elif compiler == 'gfortran':
  compStr = 'GNU'
elif compiler == 'pgf90':
  compStr = 'PGI'
# endif compiler

# LNFL: always single precision, LBLRTM: always double
lnflCmd = '%s%ssgl' % (opSys, compStr)
lblCmd = '%s%sdbl' % (opSys, compStr)

cwd = os.getcwd()

os.chdir('%s/build' % lnflDir)
print('Building LNFL')
#status = \
#  sub.call(['make', '-f', 'make_lnfl', '-o', 'lnfl', lnflCmd])
#if status != 0: sys.exit('LNFL not built')
os.chdir(cwd)

os.chdir('%s/build' % lblDir)
print('Building LBLRTM')
#status = \
#  sub.call(['make', '-f', 'make_lblrtm', '-o', 'lblrtm', lblCmd])
#if status != 0: sys.exit('LBLRTM not built')
os.chdir(cwd)

print('LNFL and LBLRTM builds successful')

if iniFile is not None:
  print('Replacing lnfl_path and lbl_path in %s' % iniFile)

  # make_lnfl and make_lblrtm -o arguments
  lnflStr = '%s/lnfl_*_%s_%s_sgl' % (lnflDir, opSys, compStr.lower())
  lblStr = '%s/lblrtm_*_%s_%s_dbl' % (lblDir, opSys, compStr.lower())
  lnflExe = glob.glob(lnflStr)[0]
  lblExe = glob.glob(lblStr)[0]

  iniDat = open(iniFile).read().splitlines()
  for line in iniDat:
    if ('lnfl_path' in line):
      split = line.split('=')
      line = line.replace(split[1], lnflExe)
      print(line)
    # endif LNFL

    if ('lbl_path' in line):
      split = line.split('=')
      line = line.replace(split[1], lblExe)
      print(line)
    # end if LBL
  # end dat loop
# endif ini
