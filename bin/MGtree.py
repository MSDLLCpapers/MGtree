#!/usr/bin/env python

# Interpret genotype(s) from a set of read alignments:
#   - Construct a tree from a Newick string
#   - For each read: based on the alignments, add a count to the LCA node
#   - Interpret the counts at the genotype level
# Alfredo Gonzalez; John M. Gaspar
# Version 1.0 2024-10-03

import sys
import argparse
import json
from treeMaker import *
from utils import *

def exportLeafNodeCounts(filename, tree, verbose):
  '''
  Write sorted table of leaf node counts.
  '''
  if verbose:
    sys.stderr.write('Writing table of leaf node counts\n')
  f = openWrite(filename)

  # write rows where readCount > 0
  table = tree.getLeafReadCounts()
  count = 0.0
  for name, readCount in table:
    if readCount > 0:
      f.write('%s\t%.2f\n' % (name, readCount))
      count += readCount

  if verbose:
    sys.stderr.write('  Total leaf node counts: %.1f\n' % count)
  closeFile(f)

def exportGenotypeTable(filename, tree, verbose):
  '''
  Write genotype table.
  '''
  if verbose:
    sys.stderr.write('Writing genotype file\n')

  table, ambigLCA, dupLCA, total = tree.getGenotypeInfo()
  if not total:
    sys.stdout.write('Error! No counts to report\n')
    sys.exit(0)
  f = openWrite(filename)
  f.write('Genotype\tCount\tPercentage\n')
  countCheck = 0
  for genotype, count in table:
    pct = 100.0 * count / total
    f.write('%s\t%d\t%.2f%%\n' % (genotype, count, pct))
    countCheck += count

  if verbose:
    sys.stderr.write('  Total LCA counts: %d\n' % total)
    sys.stderr.write('    Written to table: %d\n' % countCheck)
    if dupLCA:
      sys.stderr.write('    Duplicate LCA counts: %d\n' % dupLCA)
    if ambigLCA:
      sys.stderr.write('    Ambiguous LCA counts: %d\n' % ambigLCA)
  closeFile(f)

def addAlnsToTree(filename, tree, verbose):
  '''
  Process read/fragment alignments from given tsv file.
    Add counts to nodes in tree.
  '''
  if verbose:
    sys.stderr.write('Reading in tsv file of alignments\n')
  f = openRead(filename)
  for line in f:
    spl = line.rstrip().split('\t')
    if len(spl) < 3:
      sys.stderr.write('Warning! Skipping improperly formatted line\n'
        + line)
      continue
    if spl[1] not in ['PE', 'R1', 'R2', 'SE']:
      sys.stderr.write('Warning! Skipping record with unknown aln type:'
        + ' %s\n' % spl[1])
      continue
    if spl[2] in ['unmapped', 'lowMapQ']:
      continue

    # add info to tree
    alns = spl[2].split(',')
    if spl[1] == 'PE':
      tree.increaseLCACount(alns, 2, spl[0])
      for a in alns:
        tree.increaseReadCount(a, 2.0/len(alns))

    else:
      # R1/R2 or SE
      tree.increaseLCACount(alns, 1, spl[0])
      for a in alns:
        tree.increaseReadCount(a, 1.0/len(alns))

  closeFile(f)

def main(args):
  '''
  Main.
  '''
  # construct tree from Newick string
  if args.verbose:
    sys.stderr.write('Reading in Newick tree\n')
  f = openRead(args.nwkfile)
  tree = Tree(f.read(), args.saveGT, args.verbose)
  closeFile(f)

  # load alignment info from tsv file
  addAlnsToTree(args.infile, tree, args.verbose)

  # write genotype output table
  exportGenotypeTable(args.outfile, tree, args.verbose)

  # write table of leaf node info
  if args.tallyfile:
    exportLeafNodeCounts(args.tallyfile, tree, args.verbose)

  # write file of LCA nodes and reads assigned to them
  if args.jsonfile:
    d = dict()
    tree.getRecords(tree.root, d)
    f = openWrite(args.jsonfile)
    json.dump(d, f, sort_keys=True)
    closeFile(f)

if __name__ == "__main__":
  '''
  Parse CL args
  '''
  parser = argparse.ArgumentParser(prog=sys.argv[0],
    add_help=False)
  parser._action_groups.pop()

  required = parser.add_argument_group('Required arguments')
  required.add_argument('-i', dest='infile', required=True, metavar='<file>',
    help='Input tsv file listing read names, aln types, and references')
  required.add_argument('-n', dest='nwkfile', required=True, metavar='<file>',
    help='Input Newick file')
  required.add_argument('-o', dest='outfile', required=True, metavar='<file>',
    help='Output tsv file of genotypes')

  other = parser.add_argument_group('Optional arguments')
  other.add_argument('-t', dest='tallyfile', metavar='<file>',
    help='Output tsv file of leaf node counts')
  other.add_argument('-q', dest='jsonfile', metavar='<file>',
    help='Output JSON file of reads and assigned LCA nodes')
  other.add_argument('-g', dest='saveGT', action='store_false',
    help='Do not interpret leaf node names with \'_\' as having genotypes')
  other.add_argument('-v', '--verbose', dest='verbose',
    action='store_true', help='Run in verbose mode')
  other.add_argument('-h', '--help', dest='help', action='help',
    help='Show help message and exit')

  # parse args, call main()
  args = parser.parse_args()
  main(args)
