# for Python 3 compatibility
from __future__ import print_function

import os, sys
import numpy as np
import subprocess as sub

def log(*message):
  """
  Function from Andre Wehe code in 
  /nas/project/rc/rc2/rpernak/NH3_Retrievals/cris_python.tar

  Modified for better readability

  Used in error/exception handling, specifically file_check() so user
  has a better idea of where error originataed
  """
  import inspect

  frame = inspect.getouterframes(inspect.currentframe())
  (dum, script, topLineNum, dum, dum, dum) = frame[-1]
  (dum, script, funcLineNum, func, dum, dum) = frame[-2]

  outStr = '***\nIn %s (%s()), Line %d\n%s\n***' % \
    (script, func, funcLineNum, ' '.join(message))
  print(outStr)
# end log()

def spawn(cmd, outSplit=True):
  """
  Simplifies the call to a shell command in a Python session

  Newer version of ls(), but I did not want to delete that function 
  because I use it in so many scripts.
  
  Call:
    results = spawn(cmd)

  Input:
    cmd -- a simple string that would be used at the Unix command line 

  Keywords:
    outSplit -- boolean, string-split (space-delimiter) the standard 
      output 

  Returns:
    stOut, stErr -- lists of standard output and standard error
  """

  call = sub.Popen(cmd, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
  callout, callerr = call.communicate()
  if outSplit:
    stOut = callout.split()
  else:
    stOut = callout

  return stOut, callerr
# end spawn()

def file_check(path):
  """
  Quick check if path exists.  Use before reading a file.
  """
  if not os.path.exists(path):
    log('Could not find %s, returning' % path)
    sys.exit()
# end file_check()

def no_overwrite(path):
  """
  Quick check if path exists.  Use before writing to a file.
  """

  if os.path.exists(path):
    sys.exit('%s exists, returning' % path)
# end no_overwrite()

def ls(cmd):
  """
  was designed after lsMap() and lsPPID, but simplifies and 
  generalizes things so we can list more than just MAP and PPID 
  files.
  
  Call:
    lsFiles = ls(cmd)

  Input:
    cmd -- a simple string that would be used at the Unix command 
      line to list a set of files

  Returns:
    lsFiles -- list of files returned from the cmd input
  """

  call = sub.Popen(cmd, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
  callout, callerr = call.communicate()
  lsFiles = callout.split()

  return lsFiles
# end ls()

def pmm(data):
  """
  Prints minimum and maximum of data to screen.

    Call:
      pmm(data)
  
    Inputs:
      data -- list or array of any numeric type

    Returns:
      no value
  """

  data = np.array(data)
  if len(data) == 1:
    return data[0], data[0]
  else:
    return min(data[np.isfinite(data)]), max(data[np.isfinite(data)])
#end pmm()

def value_locate(array, value, lower=False):
  """
  Finds element in array that is closest to value

  Call:
    idx, located_value = value_locate(array, value)

  Inputs:
    array -- array of numbers
    value -- value for which to search
    
  Outputs:
    A 2-element tuple that contains index of the array element 
    that is closest to input value, and the array element

  """

  idx = (np.abs(array-value)).argmin()
  
  return idx, array[idx]

# end value_locate()

def call_CR1():
  """
  Calls a command to use the first AER Corporate Release of Python
  (Python 2.7.1)
  """

  call = sub.Popen('use_py271', shell=True, \
    stdout=sub.PIPE, stderr=sub.PIPE)
  callout, callerr = call.communicate()

  return True  
# end call_CR1

def call_CR2(path='/usr/local/CentOS6/Cr2/misc/set276env'):
  """
  Calls a command to use the second AER Corporate Release of Python
  (Python 2.7.6)
  """

  call = sub.Popen(path, shell=True, \
    stdout=sub.PIPE, stderr=sub.PIPE)
  callout, callerr = call.communicate()
  print(callerr)

  return True
# end call_CR2()

def check_CR2():
  """
  Checks whether the second Corporate Release of Python is activated
  """

  import numpy as np
  if np.__version__ != '1.8.1':
    sys.exit("Please type 'use_py276' to activate the AER " + \
      "Python Corporate Release 2")
  del np
# end check_CR2()

def check_py3():
  """
  Check if use is using Python 3 instead of Python 2 (which is what 
  AER originally was using in the Corporate Release)
  """

  return sys.version_info > (3, 0)
# end check_py3()

