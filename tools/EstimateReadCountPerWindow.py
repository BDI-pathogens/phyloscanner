#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''For each bam file in the list given as input, this
script does the following. The distribution of read lengths, and insert sizes if
reads are found to be paired, is calculated. (Here, length means length of the
mapping reference covered by the read, which will not be the same as the true
read length if there are insertions or deletions.) We then estimate the number
of reads and inserts expected to fully span a window of width W by assuming that
reads are distributed randomly over the genome (i.e. ignoring the actual
location information in the bam). We output this count for each bam file as a
function of W.'''

import os
import sys
import argparse
import pysam
import phyloscanner_funcs as pf
import collections
import numpy as np

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# A class to have new lines in argument help text
class SmartFormatter(argparse.HelpFormatter):
  def _split_lines(self, text, width):
    if text.startswith('R|'):
      return text[2:].splitlines()
    return argparse.HelpFormatter._split_lines(self, text, width)

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage,
formatter_class=SmartFormatter)

# Positional args
parser.add_argument('BamAndRefList', type=File,
help='''R|A csv-format file listing the bam and reference files
(i.e. the fasta-format file containing the sequence to
which the reads were mapped). The first column should
be the bam file, the second column the corresponding
reference file, with a comma separating the two. An
optional third column, if present, will be used to 
rename the bam files in all output. For example:
PatientA.bam,PatientA_ref.fasta,A
PatientB.bam,PatientB_ref.fasta,B''')
parser.add_argument('-N', '--normalise', action='store_true', help='''Normalise
the counts for each bam to the value at a window width of zero, making it easier
to compare the relative decline in number of reads with growing window size
between different bams with different total numbers of reads.''')
parser.add_argument('-O', '--out-filename', help="We'll append '.csv' for the "
"output data file, and '.pdf' for the plot. The default is "
"'EstimatedReadCountsPerWindow'.", default='EstimatedReadCountsPerWindow')
parser.add_argument('-OIS', '--overlapping-insert-sizes', action='store_true',
help='''Just record the insert size distribution for each bam, restricted to
inserts where the mates overlap.''')
parser.add_argument('-DB', '--dont-plot', action='store_true',
help="Don't plot the results.")
parser.add_argument('-MC', '--min-read-count', type=float, help='''Used to
specify a positive number: we'll truncate the x axis when the window width
becomes so large that all bams have a read count per window below this
value. The default is 1.''', default=1)
parser.add_argument('-AS', '--axis-font-size', type=int,
help='For the plot. The default is 15.', default=15)
parser.add_argument('-TS', '--title-font-size', type=int,
help='For the plot. The default is 15.', default=15)
parser.add_argument('-LS', '--legend-font-size', type=int,
help='For the plot. The default is 7.', default=7)
parser.add_argument('-LL', '--legend-location', 
help='''For the plot. The default is 'lower left'. The other options are:
'best', 'upper right', 'upper left', 'lower right', 'right', 'center left',
'center right', 'lower center',' upper center', 'center' ''',
default='lower left')
parser.add_argument('-LY', '--linear-y-axis', 
help='For the plot. The default is logarithmic.', action='store_true')
parser.add_argument('-XM', '--x-min-max', help='The minimum and maximum for '\
'the x axis in the plot, specified together as a comma-separated pair of '\
'numbers.')
parser.add_argument('-YM', '--y-min-max', help='The minimum and maximum for '\
'the y axis in the plot, specified together as a comma-separated pair of '\
'numbers.')
parser.add_argument('--x-samtools', default='samtools', help=\
'Used to specify the command required to run samtools, if it is needed to index'
' the bam files (by default: samtools).')
args = parser.parse_args()

InsertSizesOnly = args.overlapping_insert_sizes

def GetIntPair(arg, ArgName):
  MinMax = arg.split(',')
  if len(MinMax) != 2:
    print(ArgName, 'should be used to specify a comma-separated pair of',
    'numbers. Quitting.', file=sys.stderr)
    exit(1)
  try:
    Min = float(MinMax[0])
    Max = float(MinMax[1])
  except ValueError:
    print(ArgName, 'should be used to specify a comma-separated pair of',
    'numbers. Quitting.', file=sys.stderr)
    exit(1)
  return min(Min, Max), max(Min, Max)

# Get plot limits
if args.x_min_max:
  Xmin, Xmax = GetIntPair(args.x_min_max, '--x-min-max')
if args.y_min_max:
  Ymin, Ymax = GetIntPair(args.y_min_max, '--y-min-max')


# Read in the input bam and ref files
BamFiles, RefFiles, aliases, BamFileBasenames = \
pf.ReadInputCSVfile(args.BamAndRefList)
NumBams = len(BamFiles)

# Make index files for the bam files if needed.
pf.MakeBamIndices(BamFiles, args.x_samtools)

def FindReadCountAsFuncOfWindowWidth(ReadSizeCountDict, RefLength):

  # Return an empty array if there are no reads
  if len(ReadSizeCountDict) == 0:
    return np.zeros(0)

  LargestReadLength = max(ReadSizeCountDict.keys())
  RefLengthPlus1 = RefLength + 1

  # The nth element of this list will eventually contain the number of reads
  # expected to span a window of width n+1 (list is zero-based).
  ReadsCountByWindowWidth = np.zeros(LargestReadLength)

  for ReadLength, count in ReadSizeCountDict.items():

    ReadLengthPlus1 = ReadLength + 1

    # The number of positions at which we could place a window of width W is
    # RefLength - W + 1
    # The number of positions at which we could place a window of width W such
    # that it is wholly inside a read is ReadLength - W + 1
    # Probability of a given read overlapping a window of width W is therefore
    # (ReadLength - W + 1) / (RefLength - W + 1) 
    for W in range(1, ReadLengthPlus1):
      NumSpanningReads = count * \
      float(ReadLengthPlus1 - W) / (RefLengthPlus1 - W)
      ReadsCountByWindowWidth[W-1] += NumSpanningReads

  if args.normalise:
    ReadsCountByWindowWidth = [float(count) / ReadsCountByWindowWidth[0] \
    for count in ReadsCountByWindowWidth]

  return ReadsCountByWindowWidth


ReadLengthCountsByBam = collections.OrderedDict()
InsertSizeCountsByBam = collections.OrderedDict()
InsertSizesOnlyByBam  = collections.OrderedDict()
for i, BamFileName in enumerate(BamFiles):

  alias = aliases[i]
  print('Now counting read and insert sizes for', alias)

  bam = pysam.AlignmentFile(BamFileName, "rb")

  # Find the reference in the bam file; there should only be one.
  AllRefs = bam.references
  if len(AllRefs) != 1:
    print('Expected exactly one reference in', BamFileName + '; found',\
    str(len(AllRefs)) + '.Quitting.', file=sys.stderr)
    exit(1)
  RefName = AllRefs[0]

  # Get the length of the reference.
  AllRefLengths = bam.lengths
  if len(AllRefLengths) != 1:
    print('Pysam error: found one reference but', len(AllRefLengths),
    'reference lengths. Quitting.', file=sys.stderr)
    exit(1)
  RefLength = AllRefLengths[0]

  PairedReadCoords = {}
  ReadLengthCounts = {}
  InsertSizeCounts = {}
  TotalReadCount = 0

  # Iterate through the reads
  for read in bam.fetch(RefName):

    MappedPositions = read.get_reference_positions(full_length=False)

    # Skip unmapped reads
    if not MappedPositions:
      continue

    TotalReadCount += 1

    start = min(MappedPositions[0], MappedPositions[-1])
    end   = max(MappedPositions[0], MappedPositions[-1])
    ReadLength = end - start
    try:
      ReadLengthCounts[ReadLength] += 1
    except KeyError:
      ReadLengthCounts[ReadLength] = 1
      
    # The first time we encounter a mate from a pair, record its start and end.
    # When we encounter its mate, if they overlap, record the insert size; if
    # they don't overlap, record their separate lengths as though they are two
    # different inserts (because phyloscanner won't merge them - they are
    # effectively two separate inserts from the point of view of merging).
    if read.is_paired:
      if read.query_name in PairedReadCoords:
        MateStart, MateEnd, MateFoundBool = PairedReadCoords[read.query_name]
        PairedReadCoords[read.query_name][2] = True
        if start <= MateStart <= end:
          InsertSize = max(end, MateEnd) - start
          try:
            InsertSizeCounts[InsertSize] += 1
          except KeyError:
            InsertSizeCounts[InsertSize] = 1
        elif MateStart <= start <= MateEnd:
          InsertSize = max(end, MateEnd) - MateStart
          try:
            InsertSizeCounts[InsertSize] += 1
          except KeyError:
            InsertSizeCounts[InsertSize] = 1
        else:
          try: 
            InsertSizeCounts[ReadLength] += 1
          except KeyError:
            InsertSizeCounts[ReadLength] = 1
          MateLength = MateEnd - MateStart
          try: 
            InsertSizeCounts[MateLength] += 1
          except KeyError:
            InsertSizeCounts[MateLength] = 1
      else:
        PairedReadCoords[read.query_name] = [start, end, False]

  # For paired reads for which we didn't find a mate, add just the read length
  # to the insert size distribution.
  NumMissingMates = 0
  for start, end, MateFound in PairedReadCoords.values():
    if not MateFound:
      NumMissingMates += 1
      ReadLength = end - start
      try: 
        InsertSizeCounts[ReadLength] += 1
      except KeyError:
        InsertSizeCounts[ReadLength] = 1
  if NumMissingMates > 0:
    print('Info:', NumMissingMates, 'of', TotalReadCount, 'reads in',
    BamFileName, "are flagged as being paired but don't have a mate present.")

  # Skip empty bams
  if TotalReadCount == 0:
    print('Warning: no reads found in', BamFileName + '. Skipping.')
    continue

  if InsertSizesOnly:
    InsertSizesOnlyByBam[alias] = InsertSizeCounts

  ReadLengthCountsByBam[alias] = \
  FindReadCountAsFuncOfWindowWidth(ReadLengthCounts, RefLength)
  InsertSizeCountsByBam[alias] = \
  FindReadCountAsFuncOfWindowWidth(InsertSizeCounts, RefLength)

if InsertSizesOnly:
  with open(args.out_filename + '.csv', 'w') as f:
    f.write('Bam file,Size of overlapping read pair or length of read in ' + \
    'non-overlapping pair,Count\n')
    for alias, InsertSizesOnly in InsertSizesOnlyByBam.items():
      for size, count in sorted(InsertSizesOnly.items(), key=lambda x:x[0]):
        f.write(alias + ',' + str(size) + ',' + str(count) + '\n')
  exit(0)  

# Make a matrix for which the first column is every window size we need to
# consider, in order, and subsequent columns list the number of reads (and
# inserts, if reads are paired) expected to fully span a window of that size,
# for each different bam.
MaxInsertSize = max(len(list_) for list_ in InsertSizeCountsByBam.values())
SomeDataIsPaired = MaxInsertSize > 0
MaxReadOrInsertSize = max(MaxInsertSize,
max(len(list_) for list_ in ReadLengthCountsByBam.values()))
if SomeDataIsPaired:
  matrix = np.zeros((MaxReadOrInsertSize, 2 * NumBams + 1))
else:
  matrix = np.zeros((MaxReadOrInsertSize, NumBams + 1))
matrix[:, 0] = np.arange(1, MaxReadOrInsertSize + 1)
header = 'window width'
if SomeDataIsPaired:
  for alias in aliases:
    header += ',' + 'read count in ' + alias + ',insert size count in ' + alias
  for i, ReadLengthCounts in enumerate(ReadLengthCountsByBam.values()):
    matrix[:len(ReadLengthCounts), 2 * i + 1] = ReadLengthCounts
  for i, InsertSizeCounts in enumerate(InsertSizeCountsByBam.values()):
    matrix[:len(InsertSizeCounts), 2 * i + 2] = InsertSizeCounts
else:
  for alias in aliases:
    header += ',' + 'read count in ' + alias
  for i, ReadLengthCounts in enumerate(ReadLengthCountsByBam.values()):
    matrix[:len(ReadLengthCounts), i + 1] = ReadLengthCounts

# Write the matrix to a csv file.
with open(args.out_filename + '.csv', 'w') as f:
  np.savetxt(f, matrix, delimiter=',', header=header, fmt='%.1f')

if args.dont_plot:
  exit(0)

try:
  import matplotlib.pyplot as plt
except ImportError:
  print("The python library matplotlib does not seem to be installed: you'll "
  "need to plot", args.out_filename + '.csv yourself.' )
  exit(1)

# For plotting: cut off the tail end of the matrix where read counts are too
# small.
LastDesiredRow = 0
for row in range(MaxReadOrInsertSize - 1, -1, -1):
  if max(matrix[row, 1:]) >= args.min_read_count:
    LastDesiredRow = row
    break
if LastDesiredRow == 0:
  print('Warning: no bam has', args.min_read_count, 'reads per window',
  'regardless how small the window is. Ignoring the --min-read-count value.')
  LastDesiredRow = MaxReadOrInsertSize - 1
matrix = matrix[:LastDesiredRow + 1, :]

ax = plt.figure().add_subplot(111)

if args.x_min_max:
  ax.set_xlim(xmin=Xmin, xmax=Xmax)
if args.y_min_max:
  ax.set_ylim(ymin=Ymin, ymax=Ymax)

for i in range(1, matrix.shape[1]):
  if SomeDataIsPaired:
    alias = aliases[(i - 1) / 2]
    if i % 2 == 0:
      label = 'read pairs, ' + alias
      linestyle = '--'
    else:
      label = 'reads, ' + alias
      linestyle = '-'
  else:
    label = aliases[i - 1]
    linestyle = '-'
  plt.plot(matrix[:, 0], matrix[:, i], label=label, linestyle=linestyle)

plt.xlabel('window width', fontsize=args.axis_font_size)
YaxisLabel = 'number of reads'
if args.normalise:
  YaxisLabel += ' relative to\nwhen the window width of zero'
if SomeDataIsPaired:
  title = \
  'Estimating the number of unpaired reads and paired reads (merging\n' + \
  'read in a pair when they overlap) spanning each window, assuming\n' + \
  'reads are randomly distributed over the whole genome'
else:
  title = \
  'Estimating the number of reads spanning each window, assuming\n' + \
  'they are randomly distributed over the whole genome'
plt.ylabel(YaxisLabel, fontsize=args.axis_font_size)
plt.title(title, fontsize=args.title_font_size)
ax.tick_params(axis='both', which='major', labelsize=args.axis_font_size)
if not args.linear_y_axis:
  ax.set_yscale('log')
ax.set_xlim(xmin=0, xmax=LastDesiredRow)
plt.legend(loc=args.legend_location, fontsize=args.legend_font_size)
plt.tight_layout()
plt.savefig(args.out_filename + '.pdf')
