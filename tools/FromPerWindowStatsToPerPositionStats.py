#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Takes a csv where the first column is the start position
of the window, the second column is the end position of the window, and all
subsequent columns are values associated with that window; it converts these
per-window values to per-position values, by taking the mean of all values
overlapping a particular position. Output is printed to stdout, suitable for
redirection to a new csv file.'''

import os
import copy
import sys
import argparse

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('CsvFile', type=File)
args = parser.parse_args()

StatTotalsByPosition = {}
StatCountsByPosition = {}

with open(args.CsvFile) as f:
  for LineNumMin1, line in enumerate(f):

    # Read in the header
    if LineNumMin1 == 0:
      header = 'position,' + line.split(',',2)[2]
      continue

    # Check for the number of fields
    fields = line.split(',')
    if LineNumMin1 == 1:
      if len(fields) < 3:
        print('Too few columns: need at least three. Quitting.',
        file=sys.stderr)
        exit(1)
      NumStats = len(fields) - 2
    assert len(fields) == 2 + NumStats, 'Line ' + str(LineNumMin1 + 1) + \
    ' has only ' + str(len(fields)) + ' fields. Quitting.'

    # Get the window start and end
    try:
      WindowStart = int(fields[0])
      WindowEnd   = int(fields[1])
      assert WindowEnd > WindowStart
    except ValueError, AssertionError:
      print('Error on line ', LineNumMin1 + 1, ': the first field should be ',
      'the window start, the second the window end. These should be integers, ',
      'with the latter greater than the former. Quitting.', sep='',
      file=sys.stderr)
      exit(1)

    # Get the stats
    try:
      stats = [float(value) for value in fields[2:]]
    except ValueError:
      print('Error on line ', LineNumMin1 + 1, ': unable to understand the ',
      'values after the window coordinates as floats. Quitting.', sep='',
      file=sys.stderr)
      exit(1)

    for pos in range(WindowStart, WindowEnd+1):
      if pos in StatTotalsByPosition:
        for i in range(NumStats):
          StatTotalsByPosition[pos][i] += stats[i]
        StatCountsByPosition[pos] += 1
      else:
        StatTotalsByPosition[pos] = copy.deepcopy(stats)
        StatCountsByPosition[pos] = 1

# Check we have data
if not StatTotalsByPosition:
  print('Found no data in', args.CsvFile + '. Quitting.', file=sys.stderr)
  exit(1)

# Print the output
print(header.rstrip())
for position, StatsTotals in sorted(StatTotalsByPosition.items(),
key=lambda x: x[0]):
  count = StatCountsByPosition[position]
  MeanStats = [float(stat) / count for stat in StatsTotals]
  print(position, ','.join(map(str,MeanStats)), sep=',')

