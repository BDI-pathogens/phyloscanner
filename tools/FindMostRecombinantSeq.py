#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Considers all triplets of sequences (seqs) in an
alignment and finds the maximum recombination signal.
  
For each possible set of three seqs in the alignment, one seq is considered
the putative recombinant and the other two the parents. For each possible 
'break point' (the point at which recombination occurred), we calculate d_L
as the difference between the Hamming distance from the recombinant to one
parent and the Hamming distance from the recombinant to the other parent,
looking to the left of the break point only; similarly we calculate d_R
looking to the right of the break point only. d_L and d_R are signed integers,
such that their differing in sign indicates that the left and right sides of
the recombinant look like different parents. We maximise the difference
between d_L and d_R (over all possible sets of three sequences and all
possible break points), take the smaller of the two absolute values, and
normalise it by half the alignment length (ignoring sites that are wholly
gaps). This means that the maximum possible score of 1 is obtained if and only
if the two parents disagree at every site, the break point is exactly in the
middle, and either side of the break point the recombinant agrees perfectly
with one of the parents e.g. AAAAAAA, AAAACCC, CCCCCCC.

For speed, Hamming distances are only calculated indirectly - looking only at
informative sites, and considering only changes in distance each time the
break point is slid through the next site. However runtime unavoidably scales
as N^3, where N is the number of sequences.

The script prints to stdout four values separated by spaces: the recombination
metric, the ID of parent 1, the ID of parent 2, and the ID of the recombinant.
If the metric is exactly zero, i.e. no recombination at all, the three sequence
IDs will all be 'None'.'''

import argparse
import os
import sys
from Bio import AlignIO
from phyloscanner_funcs import CalculateRecombinationMetric

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File)
parser.add_argument('-G', '--gap-aware', action='store_true',
help='''By default, when calculating Hamming distances for the recombination
metric, positions with gaps are ignored. This means that e.g. the following
three sequences would have a metric of zero: A-AAAA, A-AAA-A, AAAA-A. With this
option, the gap character counts as a fifth base and so (dis)agreement in gaps
contributes to Hamming distance. This increases sensitivity of the metric to
cases where indels are genuine signals of recombination, but decreases
specificity, since misalignment may falsely suggest recombination.''')
args = parser.parse_args()

alignment = AlignIO.read(args.alignment, "fasta")

result = CalculateRecombinationMetric(alignment, args.gap_aware)

print(' '.join(map(str, result)))