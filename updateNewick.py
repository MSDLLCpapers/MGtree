#!/usr/bin/python

# Update a Newick tree:
#   - assign names to all nodes
#   - (optional) add genotypes from given csv file
# Alfredo Gonzalez; John M. Gaspar
# Version 1.0 2024-10-03

import sys
import argparse
import re
from treeMaker import *
from utils import *

def loadNames(filename, verbose):
  '''
  Loads node name and genotype from a headerless csv.
  '''
  f = openRead(filename)
  count = 0
  d = dict()
  for line in f:
    spl = line.rstrip().split(',')
    if len(spl) != 2:
      sys.stderr.write('Warning! Skipping line in CSV file that ' \
        + 'does not have 2 columns: %s' % line)
      continue
    if spl[0] in d:
      sys.stderr.write('Warning! Skipping line in CSV file for ' \
        + 'repeated node: %s\n' % spl[0])
      continue
    # check for forbidden characters
    if re.search(r'[;(),:_]', spl[1]):
      sys.stderr.write('Warning! Skipping line in CSV file that ' \
        + 'has a forbidden character [;(),:_]: %s' % line)
      continue
    d[spl[0]] = spl[1]
    count += 1
  closeFile(f)
  if verbose:
    sys.stderr.write('Node genotypes to be updated: %d\n' % count)
  return d

def main(args):
  '''
  Main.
  '''
  # Read in Tree
  if args.verbose:
    sys.stderr.write('Reading in Newick tree\n')
  f = openRead(args.infile)
  tree = Tree(f.read(), False, args.verbose)
  closeFile(f)

  # load genotype updates
  if args.genfile:
    d = loadNames(args.genfile, args.verbose)
    count = 0
    for name, genotype in d.items():
      if tree.updateGenotype(name, genotype):
        count += 1
    if args.verbose:
      sys.stderr.write('Nodes updated: %d\n' % count)

  # write new Newick tree
  if args.verbose:
    sys.stderr.write('Writing Newick tree\n')
  newick = tree.generateNewick(tree.root)
  f = openWrite(args.outfile)
  f.write(newick)
  closeFile(f)

if __name__ == "__main__":
  '''
  Parse CL args
  '''
  parser = argparse.ArgumentParser(prog=sys.argv[0],
    add_help=False)
  parser._action_groups.pop()

  required = parser.add_argument_group('Required arguments')
  required.add_argument('-i', dest='infile', required=True,
    metavar='<file>', help='Input Newick file')
  required.add_argument('-o', dest='outfile', required=True,
    metavar='<file>', help='Output Newick file')

  other = parser.add_argument_group('Optional arguments')
  other.add_argument('-n', dest='genfile', metavar='<file>',
    help='Input csv file to genotype nodes')
  other.add_argument('-v', '--verbose', dest='verbose',
    action='store_true', help='Run in verbose mode')
  other.add_argument('-h', '--help', dest='help', action='help',
    help='Show help message and exit')

  # parse args, call main()
  args = parser.parse_args()
  main(args)
