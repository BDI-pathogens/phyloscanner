#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This extracts named reads from a bam file.'''

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
parser.add_argument('ReadName', nargs='+')
args = parser.parse_args()

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
ReadNames = {name:None for name in args.ReadName}

# Iterate through the reads
for read in InBam.fetch(RefName):
  if read.query_name in ReadNames:
   OutBam.write(read)

