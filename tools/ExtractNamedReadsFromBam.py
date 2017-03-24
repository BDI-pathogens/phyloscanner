#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This extracts named reads from a bam file. The names of
the desired reads can either be passed directly as arguments, or listed in a
file, one per line.'''

import os
import sys
import argparse
import pysam

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('InBamFile', type=File)
parser.add_argument('OutBamFile')
parser.add_argument('-F', '--read-name-file', type=File)
parser.add_argument('-N', '--read-names', nargs='+')
args = parser.parse_args()

ReadNamesAsArgs = args.read_names is not None
ReadNamesAsFile = args.read_name_file is not None
if (ReadNamesAsArgs and ReadNamesAsFile) or ((not ReadNamesAsArgs) and
(not ReadNamesAsFile)):
  print('Exactly one of the --read-name-file and --read-names options should',
  'be used. Quitting.', file=sys.stderr)
  exit(1)

if ReadNamesAsArgs:
  ReadNames = args.read_names
else:
  ReadNames = []
  with open(args.read_name_file, 'r') as f:
    for line in f:
      ReadNames.append(line.strip())


InBam = pysam.AlignmentFile(args.InBamFile, "rb")

# Find the reference in the bam file; there should only be one.
AllReferences = InBam.references
if len(AllReferences) != 1:
  print('Expected exactly one reference in', args.InBamFile+'; found',\
  str(len(AllReferences))+'.Quitting.', file=sys.stderr)
  exit(1)
RefName = AllReferences[0]

OutBam = pysam.AlignmentFile(args.OutBamFile, "wb", template=InBam)

# Hash the desired read names for speed.
ReadNamesDict = {name:False for name in ReadNames}

# Iterate through the reads
for read in InBam.fetch(RefName):
  if read.query_name in ReadNamesDict:
    OutBam.write(read)
    ReadNamesDict[read.query_name] = True

ReadsNotFound = [read for read, found in ReadNamesDict.items() if not found]
if len(ReadsNotFound) != 0:
  print('Error: the following reads were not found in', args.InBamFile + \
  ':\n', ' '.join(ReadsNotFound) + '\nQuitting.', file=sys.stderr)
  exit(1)


