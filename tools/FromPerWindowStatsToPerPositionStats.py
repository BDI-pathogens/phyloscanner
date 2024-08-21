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
overlapping a particular position. 'nan' values are ignored when taking the mean.
Output is printed to stdout, suitable for redirection to a new csv file.'''

import os
import copy
import sys
import argparse
import math

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('CsvFile', help='The first line should be header', type=File)
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
    except (ValueError, AssertionError):
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

    # For each position in this window and for each stat: increment the total
    # for this stat (obtained when summing over all windows) by its value here,
    # and increment by 1 the count for how many times we've done this.
    # Except for stats with a nan value, then we do neither of these things.
    for pos in range(WindowStart, WindowEnd+1):
      if pos in StatTotalsByPosition:
        for i in range(NumStats):
          if not math.isnan(stats[i]):
            StatTotalsByPosition[pos][i] += stats[i]
            StatCountsByPosition[pos][i] += 1
      else:
        StatTotalsByPosition[pos] = [0 if math.isnan(stat) else stat for stat in stats]
        StatCountsByPosition[pos] = [0 if math.isnan(stat) else 1    for stat in stats]

# Check we have data
if not StatTotalsByPosition:
  print('Found no data in', args.CsvFile + '. Quitting.', file=sys.stderr)
  exit(1)

# Print the output
print(header.rstrip())
for position, StatTotals in sorted(StatTotalsByPosition.items(),
key=lambda x: x[0]):
  StatCounts = StatCountsByPosition[position]
  MeanStats = [float('nan') if count == 0 else float(StatTotals[index]) / count
  for index,count in enumerate(StatCounts)]
  print(position, ','.join(map(str,MeanStats)), sep=',')

