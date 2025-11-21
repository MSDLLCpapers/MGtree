# utility functions for file I/O
# John M. Gaspar
# Version 1.0 2024-10-03

import sys
import gzip

def openRead(filename):
  '''
  Open filename for reading. '-' indicates stdin.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdin
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'rt')
    else:
      f = open(filename, 'r')
  except IOError:
    sys.stderr.write('Error! Cannot open "%s"' % filename
      + ' for reading\n')
    sys.exit(-1)
  return f

def openWrite(filename):
  '''
  Open filename for writing. '-' indicates stdout.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdout
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'wt')
    else:
      f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot open "%s"' % filename
      + ' for writing\n')
    sys.exit(-1)
  return f

def closeFile(f):
  '''
  Close the file
  '''
  if f != sys.stdin and f != sys.stdout:
    f.close()
