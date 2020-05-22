#!/usr/bin/env python

import os, sys
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
      sys.exit('{} is not a supported compiler.'.format(compiler))

    iniFile = args.config_file
    if iniFile is not None: utils.file_check(iniFile)

    if lnfl:
      lnflDir = args.lnfl_path; utils.file_check(lnflDir)

    if lbl:
      lblDir = args.lblrtm_path; utils.file_check(lblDir)

    # we do not expect the line file to exist now that it is no
    # longer a submodule
    if lines: linesDir = args.lines_path

    self.compiler = str(compiler)
    self.iniFile = None if iniFile is None else str(iniFile)
    self.lines = False

    self.topDir = str(inArgs['top_dir'])

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
      self.zRecID = args.record_id
    else:
      sys.exit('No model build chosen')
    # endif model
  # end constructor

  def checkLineFile(self):
    """
    Check if line file dataset from Zenodo needs to be downloaded
    """

    # do we need to dowload the Line File Parameter Database?
    self.lfpdDL = True

    # it is assumed that if AER_Line_File exists, that everything
    # that is needed is underneath the directory
    if os.path.exists(self.modelDir):
      self.lfpdDL = False
      return
    # endif modelDir

    print('Generating list of Zenodo URLs')
    wgetList = 'line_file_list.txt'
    sub.call([self.zGet, str(self.zRecID), '-w', wgetList])
    files = open(wgetList).read().splitlines()

    # we archive a tarball and a license, just need tarball
    for file in files:
      base = os.path.basename(file)
      if '.tar.gz' not in base: continue
      if os.path.exists(base): self.lfpdDL = False
      self.tarBall = base

      # AER convention is that the tarball is just the line file
      # directory name with ".tar.gz" appended
      self.tarDir = base[:-7]
    # end file loop
  # end checkLineFile()

  def getLineFile(self):
    """
    Retrieve line file dataset from Zenodo, extract archive, then
    stage files as expected by LNFL and LBLRTM
    """

    import tarfile
    from zenodo_get.__main__ import zenodo_get as zget

    # what can be downloaded? should just be a tarball and license
    arcList = 'zenodo_archive_list.txt'
    zget([str(self.zRecID), '-w', arcList])
    files = open(arcList).read().splitlines()
    for file in files:
      base = os.path.basename(file)
      if '.tar.gz' not in base: continue
      if os.path.exists(base): self.lfpdDL = False
      self.tarBall = base

      # AER convention is that the tarball is just the line file
      # directory name with ".tar.gz" appended
      self.tarDir = base[:-7]
    # end file loop

    zget([str(self.zRecID)])

    if not os.path.exists(self.modelDir):
      if not os.path.exists(self.tarDir):
        # don't untar if it's not needed
        print('Extracting {}'.format(self.tarBall))
        with tarfile.open(self.tarBall) as tar: tar.extractall()
      # endif tar
      os.rename(self.tarDir, self.modelDir)
    else:
      print('{} already exists, using its Line File contents'.
        format(self.modelDir))
    # endif modelDir
  # end getLineFile()

  def build(self):
    """
    Build LBLRTM or LNFL
    """

    import subprocess as sub

    # OS determination
    # https://docs.python.org/2/library/sys.html#sys.platform
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

    cmd = '{}{}{}'.format(self.opSys, self.compStr, self.precision)

    cwd = os.getcwd()

    os.chdir('{}/build'.format(self.modelDir))
    print('Building {}'.format(self.modelStr))
    status = sub.call(['make', '-f', self.makeStr, cmd])
    if status != 0: sys.exit('{} not built'.format(self.modelStr))
    os.chdir(cwd)

    return self
  # end build()

  def configFile(self):
    """
    Replace paths in ABSCO_tables.py configuration file with the paths
    established in this class
    """

    import glob

    if self.iniFile is None:
      sys.exit('No configuration file specified, returning')
    else:
      iniDat = open(self.iniFile).read().splitlines()
      outFP = open(self.iniFile, 'w')

      if self.lines:
        # making some assumptions about directory structure here...
        modStr = ['line_file/{}'.format(self.tarDir), \
          'line_file/lncpl_lines', 'extra_brd_params', \
          'xs_files/xs', 'xs_files/FSCDXS']
      else:
        # make_lnfl and make_lblrtm -o arguments
        modStr = '{}/{}_*_{}_{}_{}'.format(
          self.modelDir, self.modelStr.lower(), \
          self.opSys, self.compStr.lower(), self.precision)
        modExe = glob.glob(modStr)[0]
      # endif lines

      for line in iniDat:
        if self.lines:
          for old, new in zip(self.pathStr, modStr):
            print('Replacing {} in {}'.format(old, self.iniFile))
            if (old in line):
              split = line.split('=')
              line = line.replace(split[1], \
                ' {}/{}/{}'.format(self.topDir, self.modelDir, new) )
            # endif LNFL
          # end path loop
        else:
          print('Replacing {} in {}'.format(self.pathStr, self.iniFile))
          if (self.pathStr in line):
            split = line.split('=')
            line = line.replace(split[1], ' {}/{}'.format(\
              self.topDir, modExe))
          # endif LNFL
        # endif lines
        outFP.write('{}\n'.format(line))
      # end dat loop
      outFP.close()
    # endif ini
  # end configFile()
# end submodules class

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(\
    formatter_class=argparse.ArgumentDefaultsHelpFormatter, \
    description='Build LBLRTM and LNFL executables for usage in ' + \
    'ABSCO_config.ini and ABSCO_tables.py.')
  parser.add_argument('-c', '--compiler', default='ifort', \
    help='Name of compiler with which user intends to build ' + \
    '([ifort, gfortran, pgf90] are supported). Case-insensitive')
  parser.add_argument('-ini', '--config_file', \
    help='Name of configuration file in which to add exe paths.  ' + \
    'If this is not set, only the models are built and no ' + \
    'configuration file is altered.')
  parser.add_argument('-lnfl', '--lnfl_path', default='LNFL', \
    help='Path of LNFL submodule directory (top level).')
  parser.add_argument('-lbl', '--lblrtm_path', default='LBLRTM', \
    help='Path of LBLRTM submodule directory (top level).')
  parser.add_argument('-lines', '--lines_path', \
    default='AER_Line_File', help='Top-level path of AER line file.')
  parser.add_argument('-record', '--record_id', type=int, \
    default=3837550, help='Zenodo record ID for the line file.')
  parser.add_argument('-no', '--no_build', action='store_true', \
    help='If set, only the line file paths are changed in the ' + \
    'configuration file and no model builds are done.')
  parser.add_argument('-t', '--top_dir', type=str, \
    default=os.getcwd(), \
    help='Full path to top-level directory under which model ' + \
    'directories reside (since .ini file for ABSCO_tables.py ' + \
    'requires absolute paths).')
  args = parser.parse_args()

  # first replace line file paths
  subObj = submodules(vars(args), lines=True)
  subObj.getLineFile()
  subObj.configFile()
  if args.no_build: sys.exit(
    'Only replaced lines paths in {}'.format(args.config_file))

  # now do LNFL -- line file paths will not change
  subObj = submodules(vars(args), lnfl=True)
  subObj.build()
  subObj.configFile()

  # now do LBL -- LNFL build and ini path will not change
  subObj = submodules(vars(args), lbl=True)
  subObj.build()
  subObj.configFile()

# end main()
