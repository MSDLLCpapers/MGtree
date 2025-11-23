#!/usr/bin/env python

# Parse a name-sorted SAM, produce list of ref alns for each read/fragment.
# John M. Gaspar
# Version 1.0 2024-08-26
# Scott S. Norton
# Version 2.0 2025-03-11

import sys
import argparse
from utils import openRead, openWrite, closeFile
try:
  from shlex import join as shlex_join
except ImportError:

  def shlex_join(a):
    '''
    Joins a list of strings, wrapping in quotes elements that have spaces.

    a: List of elements to join

    Return: str, the joined string
    '''
    out = ''
    sep = ''
    for x in a:
      out += sep
      if any(c.isspace() for c in x):
        x = x.replace('"', '\\"')
        out += '"' + x + '"'
      else:
        out += x
      sep = ' '
    return out

def getTag(lis, tag):
  '''
  Get optional tag from a SAM record.

  lis: List of SAM tags, in format XX:Y:ZZZZZ
    XX = tag ID
    Y = type code
    ZZZZZ = value
  tag: str, tag ID to find

  Return:
  optional str: the associated tag value, or None if not found
  '''
  for t in lis:
    spl = t.split(':')
    if spl[0] == tag:
      return spl[-1]
  return None

def loadHeader(f, fSam, unsorted):
  '''
  Parse SAM header. Check sort order.

  f: Input SAM file
  fSam: Output SAM file
  unsorted: True to permit input SAM not in name-sorted order

  Returns: 2-tuple of (first line past the header, total number of SQ)
  '''
  count = 0
  line = f.readline()
  pp = None
  good = False  # queryname sorted?
  while line:
    if line[0] != '@':
      break
    spl = line.rstrip().split('\t')
    if spl[0] == '@HD':
      for s in spl[1:]:
        div = s.split(':')
        if div[0] == 'SO':
          good = div[1] == 'queryname'
    elif spl[0] == '@SQ':
      count += 1
    elif spl[0] == '@PG':
      for x in spl[1:]:
        tag, value = x.split(':')
        if tag == 'ID':
          pp = value
          break
    if fSam:
      fSam.write(line)
    line = f.readline()
  if not good and not unsorted:
    sys.stderr.write('Error! SAM must be name-sorted\n')
    sys.exit(-1)
  if fSam:
    # Add PG line for this program
    pgid = __file__.split('/')[-1]
    if pp is not None:
      fSam.write('@PG\tID:%s\tPN:%s\tVN:2.0\tPP:%s\tCL:%s\n' % (pgid, pgid, pp, shlex_join(sys.argv)))
    else:
      fSam.write('@PG\tID:%s\tPN:%s\tVN:2.0\tCL:%s\n' % (pgid, pgid, shlex_join(sys.argv)))
  return line, count

def processAlns(fOut, fSam, qname, arr, label, asDiff, counts):
  '''
  Process a set of alignments. Print ref names of good alignments
    (within asDiff of best AS).

  fOut: Output TSV
  fSam: Output SAM
  qname: Read query name
  arr: List of parsed SAM records, each entry is:
    unmapped: int, nonzero if unmapped
    mapq_fail: bool, True if MAPQ less than threshold
    proper: bool, True if read is mapped and mate is mapped
    primary: bool, True if read is not secondary, supplementary, duplicate, or QC fail
    rname: str, reference name
    rpos: int, reference start position
    mpos: int, mate reference start position
    as: float, alignment score or -inf if not recorded
    spl: list[str], original SAM record
    spl2: optional list[str], mate SAM record, proper pairs only
  label: Literal['PE', 'SE', 'R1', 'R2']: read label
  asDiff: float, maximum difference between best AS and multimap AS
  counts: list[int], statistics reporting
    total alignments
    primary alignments
    secondary alignments
    unmapped reads
    reads failing MAPQ threshold
  '''
  if not arr:
    return
  fOut.write('%s\t%s\t' % (qname, label))
  if arr[0][0]:
    fOut.write('unmapped')
    counts[3] += 1
  elif arr[0][1]:
    fOut.write('lowMAPQ')
    counts[4] += 1
  else:
    # find max AS
    maxAS = float('-inf')
    for a in arr:
      if a[7] > maxAS:
        maxAS = a[7]
    # print good alns
    first = True
    for a in arr:
      if a[7] >= maxAS - asDiff:
        if not first:
          fOut.write(',')
        fOut.write(a[4])
        first = False
        if fSam:
          fSam.write('\t'.join(a[8]) + '\n')
          if a[9]:
            fSam.write('\t'.join(a[9]) + '\n')
  fOut.write('\n')

def processProper(fOut, fSam, qname, R1, R2, asDiff, counts):
  '''
  Pair up alignments that should be properly paired,
    sum AS's, and pass to processAlns()

  fOut: Output TSV
  fSam: Output SAM
  qname: Read query name
  arr: List of parsed SAM records, each entry is:
    unmapped: int, nonzero if unmapped
    mapq_fail: bool, True if MAPQ less than threshold
    proper: bool, True if read is mapped and mate is mapped
    primary: bool, True if read is not secondary, supplementary, duplicate, or QC fail
    rname: str, reference name
    rpos: int, reference start position
    mpos: int, mate reference start position
    as: float, alignment score or -inf if not recorded
    spl: list[str], original SAM record
    spl2: None for now, filled by this function
  asDiff: float, maximum difference between best AS and multimap AS
  counts: list[int], statistics reporting
    total alignments
    primary alignments
    secondary alignments
    unmapped reads
    reads failing MAPQ threshold
  '''
  # filter R1 alns
  alns = []
  for r in R1:
    # skipped unmapped and not properly paired
    if not r[0] and r[2]:
      alns.append(r)
  if not alns:
    sys.stderr.write('Error! Cannot parse alns for read %s\n' % qname)
    return

  # match R2 alns to R1s
  for r in R2:
    # skipped unmapped and not properly paired
    if not r[0] and r[2]:
      mat = False
      if r[3]:
        # match primary alns
        for s in alns:
          if s[3] and not s[9]:
            s[9] = r[8]
            s[7] += r[7]
            mat = True
      else:
        # match secondary alns based on RNAME and PNEXT
        for s in alns:
          if r[4] == s[4] and r[5] == s[6] and r[6] == s[5] and not s[9]:
            s[9] = r[8]
            s[7] += r[7]
            mat = True

      if not mat:
        sys.stderr.write('Warning! No mate for read %s R2\n' % qname)

  paired = []
  for s in alns:
    if not s[9]:
      sys.stderr.write('Warning! No mate for read %s R1\n' % qname)
      continue
    paired.append(s)
  if not paired:
    sys.stderr.write('Warning! No paired alns for read %s\n' % qname)
    return

  processAlns(fOut, fSam, qname, paired, 'PE', asDiff, counts)

def processRead(fOut, fSam, qname, lines, minMapQ, asDiff, counts):
  '''
  Process a set of alignments for a read/fragment:
    collect R1/single-end and R2 alns, and pass them
    to processPair() or processSingle() as needed.

  fOut: Output TSV
  fSam: Output SAM
  qname: Read query name
  lines: List of SAM records split by tabs
  minMapQ: int, minimum MAPQ score
  asDiff: float, maximum difference between best AS and multimap AS
  counts: list[int], statistics reporting
    total alignments
    primary alignments
    secondary alignments
    unmapped reads
    reads failing MAPQ threshold
  '''
  # separate alignments into R1 and R2 lists
  R1 = []  # list of alns for R1 (or single-end)
  R2 = []  # list of alns for R2
  pe = False
  proper = False
  for spl in lines:
    flag = int(spl[1])
    if flag & 0x1:
      pe = True
    if flag & 0x3 == 0x3:
      proper = True

    counts[1] += 1
    if flag & 0x100:
      counts[2] += 1

    # save alignment score
    score = getTag(spl[11:], 'AS')
    try:
      AS = float(score)
    except (TypeError, ValueError):
      AS = float('-inf')

    if flag & 0x40 or not pe:
      R1.append([flag & 0x4, int(spl[4]) < minMapQ, flag & 0x3 == 0x3,
        not(flag & 0xF00), spl[2], spl[3], spl[7], AS, spl, None])
    else:
      R2.append([flag & 0x4, int(spl[4]) < minMapQ, flag & 0x3 == 0x3,
        not(flag & 0xF00), spl[2], spl[3], spl[7], AS, spl, None])

  if proper:
    processProper(fOut, fSam, qname, R1, R2, asDiff, counts)
  else:
    if not pe:
      processAlns(fOut, fSam, qname, R1, 'SE', asDiff, counts)
    else:
      processAlns(fOut, fSam, qname, R1, 'R1', asDiff, counts)
      processAlns(fOut, fSam, qname, R2, 'R2', asDiff, counts)


def parseSAM(fIn, fOut, fSam, line, minMapQ, asDiff, verbose):
  '''
  Parse SAM file. Collect alns for each read/frag and send them
    to processRead().

  fIn: Input SAM
  fOut: Output TSV
  fSam: Output SAM
  line: First SAM record past the header
  minMapQ: int, minimum MAPQ score
  asDiff: float, maximum difference between best AS and multimap AS
  verbose: bool, whether to print extra logging messages (FIXME: unused)

  Returns: list[int], statistics reporting
    total alignments
    primary alignments
    secondary alignments
    unmapped reads
    reads failing MAPQ threshold
  '''
  lines = []
  counts = [0, 0, 0, 0, 0]
  qname = line.rstrip().split('\t')[0]
  while line:
    if line[0] == '@':
      sys.stderr.write('Warning! Header line found in SAM body\n')
      sys.stderr.write(line)
      continue
    spl = line.rstrip().split('\t')
    if len(spl) < 11:
      sys.stderr.write('Error! Poorly formatted SAM record\n')
      sys.stderr.write(line)
      sys.exit(-1)

    if spl[0] == qname:
      lines.append(spl)
    else:
      # new read: process previous read
      processRead(fOut, fSam, qname, lines, minMapQ, asDiff, counts)
      counts[0] += 1

      # reset for next read
      qname = spl[0]
      lines = [spl]

    line = fIn.readline()

  # process last read
  if lines:
    processRead(fOut, fSam, qname, lines, minMapQ, asDiff, counts)
    counts[0] += 1
  return counts

def main():
  '''
  Main.
  '''
  # set command-line arguments
  parser = argparse.ArgumentParser(prog=sys.argv[0],
    add_help=False)
  parser._action_groups.pop()

  required = parser.add_argument_group('Required arguments')
  required.add_argument('-i', dest='infile', required=True,
    metavar='<file>', help='Input SAM file (name-sorted, with header;'
    + ' use \'-\' for stdin)')
  required.add_argument('-o', dest='outfile', required=True,
    metavar='<file>', help='Output tsv file of reads, aln types, and refs')

  other = parser.add_argument_group('Optional arguments')
  required.add_argument('-O', dest='outsam',
    metavar='<file>', help='Output SAM file of reads filtered by alignment score and proper pairing')
  other.add_argument('-m', dest='minMapQ', type=int, default=0,
    metavar='<int>', help='Minimum MAPQ (def. 0)')
  other.add_argument('-s', dest='asDiff', type=float, default=0.0,
    metavar='<float>', help='Keep sec alns with AS >= bestAS - '
    + '<float> (def. 0)')
  other.add_argument('-S', dest='unsorted', action='store_true',
    help='Skip name-sorting requirement')
  other.add_argument('-v', '--verbose', dest='verbose',
    action='store_true', help='Run in verbose mode')
  other.add_argument('-h', '--help', dest='help', action='help',
    help='Show help message and exit')

  # parse req. arguments, open files
  args = parser.parse_args()
  fIn = openRead(args.infile)
  fOut = openWrite(args.outfile)
  if args.outsam:
    fSam = openWrite(args.outsam)
  else:
    fSam = None

  # check SAM header
  line, count = loadHeader(fIn, fSam, args.unsorted)
  if args.verbose:
    sys.stderr.write('Reference sequences in SAM: %d\n' % count)

  # load info from SAM
  counts = parseSAM(fIn, fOut, fSam, line, args.minMapQ, args.asDiff,
    args.verbose)

  if args.verbose:
    sys.stderr.write('Fragments processed: %d\n' % counts[0])
    sys.stderr.write('  Total alignments analyzed: %d\n' % counts[1])
    sys.stderr.write('    Secondary: %d\n' % counts[2])
    sys.stderr.write('    Unmapped: %d\n' % counts[3])
    if args.minMapQ:
      sys.stderr.write('    MAPQ<%d: %d\n' % (args.minMapQ, counts[4]))

  # close files
  closeFile(fIn)
  closeFile(fOut)
  if fSam:
    closeFile(fSam)

if __name__ == '__main__':
  main()
