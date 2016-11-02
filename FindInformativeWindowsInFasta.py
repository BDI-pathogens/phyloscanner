#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script picks out a series of window coordinates
with respect to an alignment, such that each window contains the desired number
of columns with each column weighted by it's non-gap fraction. e.g. a gapless
column counts for 1, a column that's half gaps and half bases counts for 0.5,
etc. The window coordinates are printed as a left coord then the corresponding
right coord, then the next left coord then its right coord, etc. all separated
by commas. For example the output
1,10,5,15
would mean windows 1-10, 5-15.'''

import argparse
import os
import sys
from Bio import AlignIO
import collections

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File)
parser.add_argument('WeightedWindowWidth', type=int, \
help='How wide do you want your windows to be?')
parser.add_argument('Overlap', type=int, help='The overlap between the end of'+\
' one window and the start of the next (measured in numbers of base pairs, or'\
+'columns in the alignment without weighting by non-gap fraction).')
parser.add_argument('-S', '--start-position', type=int, default=1, \
help='The position in the alignment from which to start creating windows (by '+\
'default, 1, i.e. the beginning of the alignment).')
parser.add_argument('-E', '--end-position', type=int, help='''The position in
the alignment to the right of which you don't want any windows. (By default, we
find windows reaching up to the end of the alignment).''')
args = parser.parse_args()

# Read in the alignment. Check there is data.
try:
  alignment = AlignIO.read(args.FastaFile, "fasta")
except:
  print('Problem reading', args.FastaFile + ':', file=sys.stderr)
  raise
NumSeqs = len(alignment)
NumCols = alignment.get_alignment_length()
if NumCols == 0:
  print('All sequences in', args.FastaFile, 'are empty. Quitting.', \
  file=sys.stderr)
  exit(1)

# Argument sanity checks.
if args.WeightedWindowWidth <= 0:
  print('WeightedWindowWidth must be greater than zero. Quitting.', \
  file=sys.stderr)
  exit(1)
StartPos = args.start_position
if StartPos < 1:
  print('The start position argument must be greater than zero. Quitting.', \
  file=sys.stderr)
  exit(1)
elif StartPos >= NumCols:
  print('The start position argument must be less than the alignment length ('+\
  str(NumCols)+'). Quitting.', file=sys.stderr)
  exit(1)

if args.end_position != None:
  if args.end_position <= StartPos:
    print('The end position argument must be greater than the start position.', 
    'Quitting.', file=sys.stderr)
    exit(1)
  EndPos = min(args.end_position, NumCols)
else:
  EndPos = NumCols

def FindWindowEnd(StartPos):
  '''Finds the end of the window, weighting columns by their non-gap fraction.
  Uses zero-based positions.'''
  CurrentWindowWeight = 0
  for position in range(StartPos,EndPos):
    if position == EndPos-1:
      return position
    column = alignment[:, position]
    NonGapFrac = 1 - float(column.count('-'))/NumSeqs
    CurrentWindowWeight += NonGapFrac
    if CurrentWindowWeight >= args.WeightedWindowWidth:
      return position

# Iteratively find windows. We stop when either the right- or left-hand edge of
# the window reaches the end (the latter being possible with negative overlap).
# Windows should always move forward; they won't if there's too much (postive)
# overlap and too narrow a window. We catch this.
WindowLeftEdges = [StartPos]
WindowRightEdges = []
while True:
  LeftEdge = WindowLeftEdges[-1]
  RightEdge = FindWindowEnd(LeftEdge-1) +1  # zero-based indexing
  WindowRightEdges.append(RightEdge)
  if RightEdge == EndPos:
    break
  NextLeftEdge = RightEdge - args.Overlap
  if NextLeftEdge >= EndPos:
    break
  if NextLeftEdge <= LeftEdge:
    print('Error: one of the windows is ', LeftEdge, '-', RightEdge, \
    '; with an overlap of ', args.Overlap, ', the next window would start at ',\
    NextLeftEdge, ', whereas we need each window to be after the previous ',\
    'one. Decrease the overlap and/or increase the WeightedWindowWidth, and ',\
    'try again. Quitting.', sep='', file=sys.stderr)
    exit(1)
  WindowLeftEdges.append(NextLeftEdge)

# Sanity check
if len(WindowLeftEdges) != len(WindowRightEdges):
  print('Malfunction of the code: differing numbers of left and right window',\
  'edges. Quitting.', file=sys.stderr)
  exit(1)

print(','.join(str(WindowLeftEdges[i])+','+str(WindowRightEdges[i]) \
for i in range(len(WindowLeftEdges))))
  

