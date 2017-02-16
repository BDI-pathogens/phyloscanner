#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''For aligned fasta file arguments, this script prints the
number of columns (i.e. the length) and the number of columns weighted by the
fraction of non-gap charaters (i.e. a column for which 25% of sequences have a
gap has a weight of 0.75).'''

import argparse
import os
import sys
from Bio import AlignIO

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File, nargs='+')
parser.add_argument('-T', '--non-gap-threshold', type=float, \
help='Specify a non-gap fraction threshold, such that a column is given weight'\
+' 1 if its non-gap fraction is at least the threshold, or weight 0 if not.')
args = parser.parse_args()

WeightColumns = args.non_gap_threshold == None

for FastaFile in args.FastaFile:

  # Read in the alignments
  try:
    alignment = AlignIO.read(FastaFile, "fasta")
  except:
    print('Problem reading', FastaFile + ':', file=sys.stderr)
    raise

  # Iterate through columns
  NumSeqs = len(alignment)
  NumCols = alignment.get_alignment_length()
  NumWeightedCols = 0
  for position in range(NumCols):
    column = alignment[:, position]
    NonGapFrac = 1 - float(column.count('-'))/NumSeqs
    if WeightColumns:
      NumWeightedCols += NonGapFrac
    elif NonGapFrac >= args.non_gap_threshold:
      NumWeightedCols += 1

  print(FastaFile, NumCols, NumWeightedCols)
