#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''phyloscanner_make_trees.py infers phylogenies containing
both within- and between-host pathogen genetic diversity in windows along the
genome, using mapped reads as input. For each window, for each sample, all reads
mapped to that window are found and identical reads are collected together with
an associated count. Then all reads from all samples in each window are aligned
using MAFFT and a phylogeny is constructed using RAxML. Temporary files and
output files are written to the current working directory; to avoid overwriting
existing files, you should run phyloscanner_make_trees.py from inside an empty
directory.

More information on running phyloscanner can be found in the manual. Options are
explained below; for some options we go into greater detail in the manual,
trying to be as concise as possible here.'''

################################################################################
# The names of some files we'll create.
FileForAlignedRefs = 'RefsAln.fasta'

# Some temporary working files we'll create
TempFileForRefs = 'temp_refs.fasta'
TempFileForPairwiseUnalignedRefs = 'temp_2Refs.fasta'
TempFileForPairwiseAlignedRefs = 'temp_2RefsAln.fasta'
TempFileForReads_basename = 'temp_UnalignedReads'
TempFileForOtherRefs_basename = 'temp_OtherRefs'
TempFileForAllBootstrappedTrees_basename = 'temp_AllBootstrappedTrees'
################################################################################
GapChar = '-'

import os
import collections
import itertools
import subprocess
import sys
import re
import copy
import shutil
import glob
import time
import argparse
import pysam
from Bio import SeqIO
from Bio import Seq
from Bio import AlignIO
from Bio import Align
import tools.phyloscanner_funcs as pf

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Define a function to check directories exist, as a type for the argparse.
def Dir(MyDir):
  if not os.path.isdir(MyDir):
    raise argparse.ArgumentTypeError(MyDir + \
    ' does not exist or is not a directory.')
  return MyDir

# Define a comma-separated integers object, as a type for the argparse.
def CommaSeparatedInts(MyCommaSeparatedInts):
  try:
    ListOfInts = [int(value) for value in MyCommaSeparatedInts.split(',')]
  except:
    raise argparse.ArgumentTypeError('Unable to understand ' +\
    MyCommaSeparatedInts + ' as comma-separated integers.')
  else:
    return ListOfInts

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

parser.add_argument('-Q', '--quiet', action='store_true', help='''Turns off the
small amount of information printed to the terminal (via stdout). We'll still
print warnings and errors (via stderr), and the command you ran, for logging
purposes.''')
parser.add_argument('-V', '--verbose', action='store_true', help='''Print a
little more information than usual. The --time option also provides extra
information on progress.''')

WindowArgs = parser.add_argument_group('Window options - you must choose'
' exactly one of -W, -AW or -E')
WindowArgs.add_argument('-W', '--windows', type=CommaSeparatedInts,
help='''Used to specify a comma-separated series of paired coordinates defining
the boundaries of the windows. e.g. 1,300,301,600,601,900 would define windows
1-300, 301-600, 601-900.''')
WindowArgs.add_argument('-AW', '--auto-window-params',
type=CommaSeparatedInts, help='''Used to specify 2, 3 or 4 comma-separated
integers controlling the automatic creation of regular windows. The first
integer is the width you want windows to be. (If you are not using the
--pairwise-align-to option, columns in the alignment of all references are
weighted by their non-gap fraction; see the manual for clarification.) The
second integer is the overlap between the end of one window and the start of the
next. The third integer, if specified, is the start position for the first
window (if not specified this is 1). The fourth integer, if specified, is the
end position for the last window (if not specified, windows will continue up to
the end of the reference or references).''')
WindowArgs.add_argument('-E', '--explore-window-widths',
type=CommaSeparatedInts, help='''Use this option to explore how the number of
unique reads found in each bam file in each window, all along the genome,
depends on the window width. After this option specify a comma-separated list of
integers. The first integer is the starting position for stepping along the
genome. Subsequent integers are window widths to try. Output is written to the
file specified with the --explore-window-width-file option.''')
WindowArgs.add_argument('-EF', '--explore-window-width-file', help='Used to '
'specify an output file for window width data, when the '
'--explore-window-widths option is used. Output is in in csv format.')

RecommendedArgs = parser.add_argument_group('Options we particularly recommend')
RecommendedArgs.add_argument('-A', '--alignment-of-other-refs', type=File,
help='''Used to specify an alignment of reference sequences for inclusion with
the reads, for comparison. (which need not be
those used to produce the bam files) that will be cut into the same windows as
the bam files and included in the alignment of reads, for comparison. This is
required if phyloscanner is to analyse the trees it produces.''')
RecommendedArgs.add_argument('-2', '--pairwise-align-to', help='''By default,
phyloscanner figures out where corresponding windows are in different bam files
by creating a multiple sequence alignment containing all of the mapping
references used to create the bam files (plus any extra references included with
-A), and window coordinates are intepreted with respect to this alignment.
However using this option, the mapping references used to create the bam files
are each separately pairwise aligned to one of the extra references included
with -A, and window coordinates are interpreted with respect to this
reference. Specify the name of the reference to use after this option. Using
this option is necessary if you want to run the
tools/CalculateTreeSizeInGenomeWindows.py script and feed its output into
phyloscanner_analyse_trees.R via its --normRefFileName option (see the manual
for more details).''')
RAxMLdefaultOptions = "-m GTRCAT -p 1 --no-seq-check"
RaxmlHelp ='''Use this option to specify how RAxML is to be run, including
both the executable (with the path to it if needed), and the options. If you do
not specify anything, we will try to find the fastest RAxML exectuable available
(assuming its path is in your PATH environment variable) and use the
options''' + RAxMLdefaultOptions + '''. -m tells RAxML which evolutionary model
to use, and -p specifies a random number seed for the parsimony inferences; both
are compulsory. You may include any other RAxML options in this command. The set
of things you specify with --x-raxml need to be surrounded with one pair of
quotation marks (so that they're kept together as one option for phyloscanner
and only split up for raxml). If you include a path to your raxml binary, it may
not include whitespace, since whitespace is interpreted as separating raxml
options. Do not include options relating to bootstraps: use phyloscanner's
--num-bootstraps and --bootstrap-seed options instead. Do not include options
relating to the naming of files.'''
RecommendedArgs.add_argument('--x-raxml', help=RaxmlHelp)
RecommendedArgs.add_argument('-P', '--merge-paired-reads', action='store_true',
help='''Relevant only for paired-read data for which the mates in a pair
(sometimes) overlap with each other: merge overlapping mates into a single
(longer) read.''')
RecommendedArgs.add_argument('-CR', '--check-recombination',
action='store_true', help='''Calculates a metric of recombination for each
sample's set of reads in each window. (Recommended only if you're interested, of
course.) Calculation time scales cubically with the number of unique sequences
each sample has per window, and so is turned off by default.''')

QualityArgs = parser.add_argument_group('Options intended to mitigate the '
'impact of poor quality reads')
QualityArgs.add_argument('-I', '--discard-improper-pairs', action='store_true',
help='''For paired-read data, discard all reads that are flagged as improperly
paired.''')
QualityArgs.add_argument('-Q1', '--quality-trim-ends', type=int, help='''Used to
specify a quality threshold for trimming the ends of reads. We trim each read
inwards until a base of this quality is met.''')
QualityArgs.add_argument('-Q2', '--min-internal-quality', type=int, help=\
'''Used to specify an interal quality threshold for reads. Reads are allowed at
most one base below this quality; any read with two or more bases below this
quality are discarded.''')
QualityArgs.add_argument('-MC', '--min-read-count', type=int, default=1, help=\
'''Used to specify a minimum count for each unique read.
Reads with a count (i.e. the number of times that sequence was observed,
after merging if merging is being done) less than this value are discarded.
The default value of 1 means all reads are kept.''')

OtherArgs = parser.add_argument_group('Other assorted options')
OtherArgs.add_argument('-XC', '--excision-coords', type=CommaSeparatedInts,
help='Used to specify a comma-separated set of integer coordinates that will '
'be excised from the aligned reads. Requires the -XR flag.')
OtherArgs.add_argument('-XR', '--excision-ref', help='''Used to specify the name
of a reference (which must be present in the file you specify with -A) with
respect to which the coordinates specified with -XC are interpreted. If you are
also using the --pairwise-align-to option, you must use the same reference there
and here.''')
OtherArgs.add_argument('-MTA', '--merging-threshold-a', type=int, default=0,
help='''Used to specify a similarity threshold for merging similar reads: reads
that differ by a number of bases equal to or less than this are merged, those
with higher counts effectively absorbing those with lower counts. The default
value of 0 means there is no merging. In this version of the merging algorithm,
if B is similar enough to merge into C, and A is similar enough to merge into B
but not into C (A -> B -> C), B will be merged into C and A will be kept
separate. Contrast with --merging-threshold-b.''')
OtherArgs.add_argument('-MTB', '--merging-threshold-b', type=int, default=0,
help='''Similar to --merging-threshold-a, except that in this version of the
merging algorithm, if B is similar enough to merge into C, and A is similar
enough to merge into B but not into C (A -> B -> C), both A and B will be merged
into C.''')
OtherArgs.add_argument('-N', '--num-bootstraps', type=int,
help='Used to specify the number of bootstraps to be calculated for RAxML trees'
' (by default, none, i.e. only the ML tree is calculated).')
OtherArgs.add_argument('-Ns', '--bootstrap-seed', type=int, default=1, help='''
Used to specify the random-number seed for running RAxML with bootstraps. The
default is 1.''')
OtherArgs.add_argument('-OD', '--output-dir', help='''Used to specify the name
of a directory into which output files will be moved.''')
OtherArgs.add_argument('--time', action='store_true',
help='Print the times taken by different steps.')
OtherArgs.add_argument('--x-mafft', default='mafft', help=\
'Used to specify the command required to run mafft (by default: mafft). See '\
'also --x-mafft2.')
OtherArgs.add_argument('--x-mafft2', default='linsi', help='''"If you are using
--pairwise-align-to, by default we will use the command 'linsi' for pairwise
alignment of references regardless of what you specified with --x-mafft. Use
this option to override the linsi default with another command.''')
OtherArgs.add_argument('--x-samtools', default='samtools', help=\
'Used to specify the command required to run samtools, if needed (by default: '
'samtools).')
OtherArgs.add_argument('-KO', '--keep-output-together', action='store_true',
help='''By default, subdirectories are made for different kinds of phyloscanner
output. With this option, all output files will be in the same directory (either
the working directory, or whatever you specify with --output-dir).''')
OtherArgs.add_argument('-KT', '--keep-temp-files', action='store_true', help='''
Keep temporary files we create on the way (these are deleted by default).''')

BioinformaticsArgs = parser.add_argument_group('Options for detailed'
' bioinformatic interrogation of the input bam files (not intended for normal'
' usage)')
BioinformaticsArgs.add_argument('-IO', '--inspect-disagreeing-overlaps',
action='store_true', help='''With --merge-paired-reads, those pairs that 
overlap but disagree are discarded. With this option, these discarded pairs 
are written to a bam file (one per patient, with their reference file copied
to the working directory) for your inspection.''')
BioinformaticsArgs.add_argument('-RN1', '--read-names-1', action='store_true',
help='''Produce a file for each window and each bam, listing the names of the
reads that phyloscanner used.''')
BioinformaticsArgs.add_argument('-RN2', '--read-names-2', action='store_true',
help='''As --read-names-1, except the files will show the correspondence between
read names and which unique sequence they correspond to.''')
BioinformaticsArgs.add_argument('--exact-window-start', action='store_true',
help='''Normally phyloscanner retrieves all reads that fully
overlap a given window, i.e. starting at or anywhere before the window start,
and ending at or anywhere after the window end. If this option is used without
--exact-window-end, the reads that are retrieved are those that start at exactly
the start of the window, and end anywhere (ignoring all the window end
coordinates specified). If this option is used with --exact-window-end, for a
read to be kept it must start at exactly the window start AND end at exactly the
window end.''')
BioinformaticsArgs.add_argument('--exact-window-end', action='store_true',
help='''With this option, the reads that are retrieved are those that end at
exactly the end of the window. Read the --exact-window-start help.''')
BioinformaticsArgs.add_argument('-CE', '--recover-clipped-ends',
action='store_true', help ='''This option recovers (soft-)clipped read ends by
setting the bases they contain to be mapped, assuming no indels inside the
clipped end. WARNING: mapping software clips the ends of reads for a reason.''')

StopEarlyArgs = parser.add_argument_group('Options to only partially run '
'phyloscanner, stopping early or skipping steps')
StopEarlyArgs.add_argument('-AO', '--align-refs-only', action='store_true',
help='''Align the mapping references used to create the bam files, plus any
extra reference sequences specified with -A, then quit without doing anything
else.''')
StopEarlyArgs.add_argument('-RNO', '--read-names-only', action='store_true',
help='''To be combined with --read-names-1 or --read-names-2: quit after writing
the read names to a file (which means the reads are not aligned).''')
StopEarlyArgs.add_argument('-T', '--no-trees', action='store_true',
help='Process and align the reads from each window, then quit without making '
'trees.')
StopEarlyArgs.add_argument('-D', '--dont-check-duplicates', action='store_true',
help="Don't compare reads between samples to find duplicates - a possible "+\
"indication of contamination. (By default this check is done.)")

DeprecatedArgs = parser.add_argument_group('''Deprecated options, left here for
backward compatability or interest.''')
DeprecatedArgs.add_argument('-C', '--contaminant-count-ratio', type=float,
help='''Used to specify a numerical value which is interpreted in the following
'way: if a sequence is found exactly duplicated between any two bam files,
'and is more common in one than the other by a factor at least equal to this
'value, the rarer sequence is deleted and goes instead into a separate
contaminant read fasta file. (This is considered deprecated because including the
exact duplicates in the tree and dealing with them during tree analysis is more
flexible.)''')
DeprecatedArgs.add_argument('-CO', '--flag-contaminants-only',
action='store_true',
help="For each window, just flag contaminant reads then move on (without "
"aligning reads or making a tree). Only makes sense with the -C flag.")
#DeprecatedArgs.add_argument('-CD', '--contaminant-read-dir', type=Dir,
#help='''A directory containing the contaminant read fasta files produced by a
#previous run of phyloscanner using the the -C flag:
#reads flagged as contaminants there will be considered contaminants in this run
#too. The windows must exactly match up between these two runs, most easily
#achieved with the -2 option. The point of this option is to first do a run with
#every bam file you have in which there could conceivably be cross-contamination,
#using the -C flag (and possibly -CO to save time), and then in subsequent runs
#focussing on subsets of bam files you will be able to identify contamination
#from outside that subset.''')
DeprecatedArgs.add_argument('-FR', '--forbid-read-repeats', action='store_true',
help='''Using this option, if a read with the same name is found
to span each of a series of consecutive, overlapping windows, it is only used in
the first window.''')
DeprecatedArgs.add_argument('-O', '--keep-overhangs', action='store_true',
help='''Keep the part of the read that overhangs the edge of the window. (By
default this is trimmed, i.e. only the part of the read inside the window is
kept.)''')
DeprecatedArgs.add_argument('-RC', '--ref-for-coords', help='''By default,
window coordinates are interpreted as alignment coordinates, in the multiple
sequence alignment of all references (which phyloscanner creates). With this
option, window coordinates are interpreted with respect to a named reference in
that alignment. Specify the name of the desired reference after this option.
This is deprecated because the --pairwise-align-to option is expected to perform
better.''')
DeprecatedArgs.add_argument('-RG', '--recombination-gap-aware',
action='store_true',
help='''By default, when calculating Hamming distances for the recombination
metric, positions with gaps are ignored. With this option, the gap character
counts as a fifth base.''')
DeprecatedArgs.add_argument('-RD', '--recombination-norm-diversity',
action='store_true', help='''By default, the normalising constant for the
recombination metric is half the number of sites in the window; with this option
it's half the number of sites in the window that are polymorphic for that bam
file.''')

args = parser.parse_args()

# Shorthand
WindowCoords  = args.windows
UserSpecifiedCoords = args.windows != None
AutoWindows = args.auto_window_params != None
IncludeOtherRefs = args.alignment_of_other_refs != None
QualTrimEnds  = args.quality_trim_ends != None
ImposeMinQual = args.min_internal_quality != None
ExcisePositions = args.excision_coords != None
PairwiseAlign = args.pairwise_align_to != None
FlagContaminants = args.contaminant_count_ratio != None
#RecallContaminants = args.contaminant_read_dir != None
RecallContaminants = False
CheckDuplicates = not args.dont_check_duplicates
ExploreWindowWidths = args.explore_window_widths != None
MergeReadsA = args.merging_threshold_a > 0
MergeReadsB = args.merging_threshold_b > 0
MergeReads = MergeReadsA or MergeReadsB
PrintInfo = not args.quiet
RecombNormToDiv = args.recombination_norm_diversity

# Print how this script was called, for logging purposes.
print('phyloscanner was called thus:\n' + ' '.join(sys.argv))

# Warn if RAxML files exist already.
if glob.glob('RAxML*'):
  print('Warning: RAxML files are present in the working directory. If their',
  'names clash with those that phyloscanner will try to create, RAxML will',
  'fail to run. Continuing.', file=sys.stderr)

TempFiles = set([])

# Try to make the output dir, if desired.
HaveMadeOutputDir = False
if args.output_dir != None:
  if os.path.isdir(args.output_dir):
    HaveMadeOutputDir = True
  else:
    try:
      os.mkdir(args.output_dir)
    except:
      print('Problem creating the directory', args.output_dir+\
      ". Is this a subdirectory inside a directory that doesn't exist? (We can",
      "only create one new directory at a time - not a new directory inside",
      "another new one.) Continuing.", file=sys.stderr)
    else:
      HaveMadeOutputDir = True

# Check only one merging type is specified. Thereafter use its threshold.
if MergeReadsA and MergeReadsB:
  print('You cannot specify both --merging-threshold-a and',
  '--merging-threshold-b. Quitting.', file=sys.stderr)
  exit(1)
if MergeReadsA:
  MergingThreshold = args.merging_threshold_a
elif MergeReadsB:
  MergingThreshold = args.merging_threshold_b

# Check that window coords have been specified either manually or automatically,
# or we're exploring window widths
NumWindowOptions = len([Bool for Bool in [UserSpecifiedCoords, AutoWindows,
ExploreWindowWidths] if Bool == True])
if NumWindowOptions != 1:
  print('Exactly one of the --windows, --auto-window-params,',
  '--explore-window-widths options should specified. Quitting.',
  file=sys.stderr)
  exit(1)

# If using automatic windows (i.e. not specifying any coordinates), the user
# should not specify a reference for their coords to be interpreted with respect
# to, nor use a directory of contaminant reads (since windows must match to make
# use of contaminant reads).
if AutoWindows and args.ref_for_coords != None:
  print('The --ref-for-coords and --auto-window-params',
  'options should not be specified together: the first means your',
  'coordinates should be interpreted with respect to a named reference, and',
  "the second means you're not specfiying any coordinates. Quitting.",
  file=sys.stderr)
  exit(1)
if RecallContaminants and (not UserSpecifiedCoords):
  print('If using the --contaminant-read-dir option you must also specify',
  'windows with the --windows option, because the former requires that the',
  'windows in the current run exactly match up with those from the run that',
  'produces your directory of contaminant reads. Quitting.', file=sys.stderr)
  exit(1)

# If coords were specified with respect to one particular reference, 
# WindowCoords will be reassigned to be the translation of those coords to 
# alignment coordinates. UserCoords are the original coords, which we use for
# labelling things to keep labels intuitive for the user.
UserCoords = WindowCoords

# Find contaminant read files and link them to their windows.
if RecallContaminants:
  ContaminantFilesByWindow = {}
  ContaminantFileRegex = FileForDuplicateSeqs_basename + \
  'InWindow_(\d+)_to_(\d+)\.fasta'
  ContaminantFileRegex2 = re.compile(ContaminantFileRegex)
  LeftWindowEdges  = UserCoords[::2]
  RightWindowEdges = UserCoords[1::2]
  PairedWindowCoords = zip(LeftWindowEdges, RightWindowEdges)
  for AnyFile in os.listdir(args.contaminant_read_dir):
    if ContaminantFileRegex2.match(AnyFile):
      (LeftEdge, RightEdge) = ContaminantFileRegex2.match(AnyFile).groups()
      (LeftEdge, RightEdge) = (int(LeftEdge), int(RightEdge))
      if (LeftEdge, RightEdge) in PairedWindowCoords:
        ContaminantFilesByWindow[(LeftEdge, RightEdge)] = \
        os.path.join(args.contaminant_read_dir, AnyFile)
  if len(ContaminantFilesByWindow) == 0:
    print('Failed to find any files matching the regex', ContaminantFileRegex,
    'in', args.contaminant_read_dir + '. Quitting.', file=sys.stderr)
    exit(1)

# Check the contamination ratio is >= 1
if FlagContaminants and args.contaminant_count_ratio < 1:
  print('The value specified with --contaminant-count-ratio must be greater',
  'than 1. (It is the ratio of the more common duplicate to the less common',
  'one at which we consider the less common one to be contamination; it',
  "should probably be quite a lot larger than 1.) Quitting.", file=sys.stderr)
  exit(1)

# Flagging contaminants requires that we check for duplicates
if args.dont_check_duplicates and FlagContaminants:
  print('The --dont-check-duplicates and --contaminant-count-ratio options',
  'cannot be used together: flagging contaminants requires that we check',
  'duplicates. Quitting.', file=sys.stderr)
  exit(1)

# The -XR and -XC flags should be used together or not all.
if (ExcisePositions and args.excision_ref == None) or \
((not ExcisePositions) and args.excision_ref != None):
  print('The --excision-coords and --excision-ref options require each other:',
  'use both, or neither. Quitting.', file=sys.stderr)
  exit(1)

# --read-names-2 can't be used with read merging or position excising
if args.read_names_2 and (ExcisePositions or MergeReads):
  print('The --read-names-2 option cannot be used with --merging-threshold-a,',
  '--merging-threshold-b or --excision-coords, because they change the',
  'correspondence initially established between unique sequences and reads.',
  'Quitting''', file=sys.stderr)
  exit(1)

# Sanity checks on using the pairwise alignment option.
if PairwiseAlign:
  if args.ref_for_coords != None:
    print('Note that if the --pairwise-align-to option is used, using the',
    '--ref-for-coords as well is redundant.', file=sys.stderr)
    if args.ref_for_coords != args.pairwise_align_to:
      print('Furthermore you have chosen two different values for these flags,'
      , 'indicating some confusion as to their use. Try again.')
      exit(1)
  if ExcisePositions and args.excision_ref != args.pairwise_align_to:
    print('The --pairwise-align-to and --excision-ref options can only be',
    'used at once if the same reference is specified for both. Qutting.',
    file=sys.stderr)
    exit(1)

def CheckMaxCoord(coords, ref):
  '''Check that no coordinate is after the end of the reference with respect to
  which it is supposed to be interpreted.'''
  if max(coords) > len(ref.seq.ungap("-")):
    print('You have specified at least one coordinate (', max(coords),
    ') that is larger than the length of the reference with respect to which',
    ' those coordinates are to be interpreted - ', ref.id, '. Quitting.',
    sep='', file=sys.stderr)
    exit(1)

def SanityCheckWindowCoords(WindowCoords):
  'Check window coordinates come in pairs, all positive, the right > the left.'
  NumCoords = len(WindowCoords)
  if NumCoords % 2 != 0:
    raise ValueError('An even number of WindowCoords must be specified. '+\
    'Quitting.')
  if any(coord < 1 for coord in WindowCoords):
    raise ValueError('All WindowCoords must be greater than zero. Quitting.')
  LeftWindowEdges  = WindowCoords[::2]
  RightWindowEdges = WindowCoords[1::2]
  PairedWindowCoords = zip(LeftWindowEdges, RightWindowEdges)
  for LeftWindowEdge, RightWindowEdge in PairedWindowCoords:
    if LeftWindowEdge >= RightWindowEdge:
      raise ValueError('You specified a window as having left edge ' +\
      str(LeftWindowEdge) +' and right edge ' +str(RightWindowEdge)+\
      '. Left edges should be less than their right edges. Quitting.')
      exit(1)
  return NumCoords

# Sanity checks on user specified WindowCoords
if UserSpecifiedCoords:
  NumCoords = SanityCheckWindowCoords(WindowCoords)

# Sanity checks on auto window parameters
if AutoWindows:
  NumAutoWindowParams = len(args.auto_window_params)
  if not NumAutoWindowParams in [2,3,4]:
    print('The --auto-window-params option requires 2, 3 or 4 integers.',
    'Quitting.', file=sys.stderr)
    exit(1)
  WeightedWindowWidth = args.auto_window_params[0]
  WindowOverlap       = args.auto_window_params[1]
  if NumAutoWindowParams > 2:
    WindowStartPos    = args.auto_window_params[2]
    if WindowStartPos < 1:
      print('The start position for the --auto-window-params option must be',
      'greater than zero. Quitting.', file=sys.stderr)
      exit(1)
    if NumAutoWindowParams == 4:
      WindowEndPos = args.auto_window_params[3]
    else:
      WindowEndPos = float('inf')
  else:
    WindowStartPos = 1
    WindowEndPos = float('inf')
  if WeightedWindowWidth <= 0:
    print('The weighted window width for the --auto-window-params option must',
    'be greater than zero. Quitting.', file=sys.stderr)
    exit(1)

# Sanity checks on window-width exploration parameters
if ExploreWindowWidths:
  if args.explore_window_width_file == None:
    print('The --explore-window-widths option requires the',
    '--explore-window-width-file option. Quitting.', file=sys.stderr)
    exit(1)
  try:
    with open(args.explore_window_width_file, 'w') as f:
      pass
  except:
    print('Unable to open', args.explore_window_width_file, 'for writing. (Is',
    "it a file inside a directory that doesn't exist?). Quitting.",
    file=sys.stderr)
    raise
  if len(args.explore_window_widths) < 2:
    print('The --explore-window-widths option should be used to specify at',
    'least two parameters; use the --help option for more information.',
    'Quitting.', file=sys.stderr)
    exit(1)
  ExploreStart = args.explore_window_widths[0]
  ExploreWidths = args.explore_window_widths[1:]
  if ExploreStart < 1:
    print('The start point for windows when exploring window widths (the '+\
    'first integer specified with --explore-window-widths) cannot be less '+\
    'than 1. Quitting.', file=sys.stderr)
    exit(1)
  ExploreWidths = sorted(ExploreWidths)
  MinExploreWidth = ExploreWidths[0]
  if MinExploreWidth < 2:
    print('The minimum window width specified with --explore-window-widths',
    'should be greater than 1. Quitting.', file=sys.stderr)
    exit(1)
  MaxExploreWidth = ExploreWidths[-1]
  WindowWidthExplorationData = []
  CheckDuplicates = False

def FindExploratoryWindows(EndPoint):
  '''Returns the set of coordinates needed to step across the genome with the
  desired start, end and window width.'''
  # The EndPoint argument should be:
  # * the ref length if --ref-for-coords or --pairwise-align-to is used
  # * the length of the mapping ref if there's only one bam and no extra refs
  # * otherwise, the length of the alignment of all refs
  if EndPoint < ExploreStart + MaxExploreWidth:
    print('With the --explore-window-widths option you specified a start',
    'point of', ExploreStart, 'and your largest window width was',
    str(MaxExploreWidth) + '; one or both of these values should be',
    'decreased since the length of the reference or alignment of references',
    'with respect to which we are interpreting coordinates is only',
    str(EndPoint) + '. We need to be able to fit at least one window in',
    'between the start and end. Quitting.', file=sys.stderr)
    exit(1)
  ExploratoryCoords = []
  for width in ExploreWidths:
    NextStart = ExploreStart
    NextEnd = ExploreStart + width - 1
    while NextEnd <= EndPoint:
      ExploratoryCoords += [NextStart, NextEnd]
      NextStart += width
      NextEnd += width
  return ExploratoryCoords

# Record the names of any external refs being included.
# If we're doing pairwise alignments, we'll also need gappy and gapless copies
# of the ref chosen for pairwise alignment.
# Check that any coordinates that are to be interpreted with respect to a named
# reference do not go past the end of that reference.
ExternalRefNames = []
if IncludeOtherRefs:
  try:
    ExternalRefAlignment = AlignIO.read(args.alignment_of_other_refs, "fasta")
  except:
    print('Problem reading', args.alignment_of_other_refs + ':',
    file=sys.stderr)
    raise
  for ref in ExternalRefAlignment:
    ExternalRefNames.append(ref.id)
    if ref.id == args.pairwise_align_to:
      RefForPairwiseAlnsGappySeq = str(ref.seq)
      RefForPairwiseAlns = copy.deepcopy(ref)
      RefForPairwiseAlns.seq = RefForPairwiseAlns.seq.ungap("-")
      RefForPairwiseAlnsLength = len(RefForPairwiseAlns.seq)
      if UserSpecifiedCoords:
        CheckMaxCoord(WindowCoords, ref)
      elif ExploreWindowWidths:
        MaxCoordForWindowWidthTesting = RefForPairwiseAlnsLength
        WindowCoords = FindExploratoryWindows(MaxCoordForWindowWidthTesting)
        NumCoords = len(WindowCoords)
        UserCoords = WindowCoords
    if ref.id == args.ref_for_coords:
      if UserSpecifiedCoords:
        CheckMaxCoord(WindowCoords, ref)
      elif ExploreWindowWidths:
        MaxCoordForWindowWidthTesting = len(ref.seq.ungap("-"))
        WindowCoords = FindExploratoryWindows(MaxCoordForWindowWidthTesting)
        NumCoords = len(WindowCoords)
        UserCoords = WindowCoords
    if ref.id == args.excision_ref:
      CheckMaxCoord(args.excision_coords, ref)

# Consistency checks on flags that require a ref.
for FlagName, FlagValue in (('--ref-for-coords',  args.ref_for_coords),
('--pairwise-align-to', args.pairwise_align_to),
('--excision-ref', args.excision_ref)):
  #('--ref-for-rooting', args.ref_for_rooting),
  if FlagValue == None:
    continue
  if not IncludeOtherRefs:
    print('The', FlagName, 'flag requires the --alignment-of-other-refs',
    'flag. Quitting.', file=sys.stderr)
    exit(1)
  if not FlagValue in ExternalRefNames:
    print('Reference', FlagValue +', specified with the', FlagName,
    'flag, was not found in', args.alignment_of_other_refs +'. Quitting.',
    file=sys.stderr)
    exit(1)



# Remove duplicated excision coords. Sort from largest to smallest.
if ExcisePositions:
  args.excision_coords = list(set(args.excision_coords))
  args.excision_coords = sorted(args.excision_coords, reverse=True)

TranslateCoordsCode = pf.FindAndCheckCode('TranslateCoords.py')
FindSeqsInFastaCode = pf.FindAndCheckCode('FindSeqsInFasta.py')
FindWindowsCode     = pf.FindAndCheckCode('FindInformativeWindowsInFasta.py')

# Test RAxML works
if not args.no_trees:
  RAxMLargList = pf.TestRAxML(args.x_raxml, RAxMLdefaultOptions, RaxmlHelp)

times = []
if args.time:
  times.append(time.time())

# Read in the input bam and ref files
BamFiles, RefFiles, BamAliases, BamFileBasenames = \
pf.ReadInputCSVfile(args.BamAndRefList)
NumberOfBams = len(BamFiles)

# Don't produce duplication files if there's only one bam.
if NumberOfBams == 1:
  CheckDuplicates = False

# Read in all the reference sequences. Set each seq name to be the corresponding
# alias.
RefSeqs = []
for i,RefFile in enumerate(RefFiles):
  SeqList = list(SeqIO.parse(open(RefFile),'fasta'))
  if len(SeqList) != 1:
    print('There are', len(SeqList), 'sequences in', RefFile+'. There should',
    'be exactly 1.\nQuitting.', file=sys.stderr)
    exit(1)
  SeqList[0].id = BamAliases[i]
  RefSeqs += SeqList

def TranslateCoords(CodeArgs):
  '''Runs TranslateCoordsCode with the supplied args, and returns the results as
  a dict.'''

  try:
    CoordsString = subprocess.check_output([TranslateCoordsCode]+CodeArgs)
  except:
    print('Problem executing', TranslateCoordsCode +'. Quitting.',
    file=sys.stderr)
    raise

  CoordsDict = {}
  for line in CoordsString.splitlines():

    # Trim leading & trailing whitespace and skip blank lines
    line = line.strip()
    if line == '':
      continue

    # Each line in the output of the TranslateCoordsCode should be a sequence 
    # name then the coordinates.
    fields = line.split()
    if len(fields) != NumCoords +1:
      print('Unexpected number of fields in line\n' +line +'\nin the output '+\
      'of ' +TranslateCoordsCode+'\nQuitting.', file=sys.stderr)
      exit(1)
    SeqName = fields[0]
    coords = fields[1:]

    # Convert the coordinates to integers.
    # Where an alignment coordinate is inside a deletion in a particular
    # sequence, TranslateCoords.py returns an integer + 0.5 for the coordinate 
    # with respect to that sequence. Python won't convert such figures directly 
    # from string to int, but we can do so via a float intermediate. This rounds 
    # down, i.e. to the coordinate of the base immediately to the left of the
    # deletion.
    for i in range(len(coords)):
      if coords[i] != 'NaN':
        try:
          coords[i] = int(coords[i])
        except ValueError:
          if '.5' in coords[i]:
            coords[i] = int(float(coords[i]))
          else:
            print('Unable to understand the coordinate', coords[i],
            'as an integer in line\n' +line +'\nin the output of '+\
            TranslateCoordsCode+'\nQuitting.', file=sys.stderr)
            exit(1)
    CoordsDict[SeqName] = coords
  return CoordsDict

def SimpleCoordsFromAutoParams(RefSeqLength):
  WindowEndPosLocal = min(WindowEndPos, RefSeqLength)
  if WindowEndPosLocal < WindowStartPos + WeightedWindowWidth:
    print('With the --auto-window-params option you specified a start', 
    'point of', WindowStartPos, 'and your weighted window width was', 
    str(WeightedWindowWidth) + '; one or both of these values should be', 
    'decreased because the length of your reference or your', 
    'specified end point is only', str(WindowEndPosLocal) + '. We need to be',
    'able to fit at least one window in between the start and end. Quitting.',
    file=sys.stderr)
    exit(1)
  WindowCoords = []
  NextStart = WindowStartPos
  NextEnd = WindowStartPos + WeightedWindowWidth - 1
  while NextEnd <= WindowEndPosLocal:
    WindowCoords += [NextStart, NextEnd]
    NextStart = NextEnd - WindowOverlap + 1
    NextEnd = NextStart + WeightedWindowWidth - 1
  NumCoords = len(WindowCoords)
  UserCoords = WindowCoords

  return WindowCoords, UserCoords, NumCoords, WindowEndPosLocal


# If there is only one bam and no other refs, no coordinate translation
# is necessary - we use the coords as they are, though setting any after the end
# of the reference to be equal to the end of the reference.
if NumberOfBams == 1 and not IncludeOtherRefs:
  if args.align_refs_only:
    print('As you are supplying a single bam file and no external references,',
    "the --align-refs-only option makes no sense - there's nothing to align.",
    "Quitting.", file=sys.stderr)
    exit(1)
  RefSeqLength = len(RefSeqs[0])
  if AutoWindows:
    WindowCoords, UserCoords, NumCoords, WindowEndPos = \
    SimpleCoordsFromAutoParams(RefSeqLength)
  if ExploreWindowWidths:
    MaxCoordForWindowWidthTesting = RefSeqLength
    WindowCoords = FindExploratoryWindows(MaxCoordForWindowWidthTesting)
    NumCoords = len(WindowCoords)
    UserCoords = WindowCoords
  CoordsInRefs = {BamAliases[0] : WindowCoords}

# If there are at least two bam files, or if there is one but we're including
# other refs, we'll be aligning references and translating the user-specified
# coords with respect to each sequence, then storing those coords in a dict
# indexed by the ref's name.
else:

  if args.verbose:
    print('Now determining the correspondence between coordinates in different',
    'bam files.')

  # If we're separately and sequentially pairwise aligning our references to
  # a chosen ref in order to determine window coordinates, do so now.
  if PairwiseAlign:

    # Get the coords if auto coords
    if AutoWindows:
      WindowCoords, UserCoords, NumCoords, WindowEndPos = \
      SimpleCoordsFromAutoParams(RefForPairwiseAlnsLength)

    # Find the coordinates with respect to the chosen ref, in the alignment of
    # just the external refs - we'll need these later.
    ExternalRefWindowCoords = \
    pf.TranslateSeqCoordsToAlnCoords(RefForPairwiseAlnsGappySeq, WindowCoords)

    CoordsInRefs = {}
    for BamRefSeq in RefSeqs:

      # Align
      SeqIO.write([RefForPairwiseAlns,BamRefSeq], TempFileForPairwiseUnalignedRefs,
      "fasta")
      TempFiles.add(TempFileForPairwiseUnalignedRefs)
      with open(TempFileForPairwiseAlignedRefs, 'w') as f:
        try:
          ExitStatus = subprocess.call([args.x_mafft2, '--quiet',
          '--preservecase', TempFileForPairwiseUnalignedRefs], stdout=f)
          assert ExitStatus == 0
        except:
          print('Problem calling mafft. Quitting.', file=sys.stderr)
          raise
      TempFiles.add(TempFileForPairwiseAlignedRefs)

      # Translate.
      # The index names in the PairwiseCoordsDict, labelling the coords found by
      # coord translation, should coincide with the two seqs we're considering.
      PairwiseCoordsDict = TranslateCoords([TempFileForPairwiseAlignedRefs,
      args.pairwise_align_to] + [str(coord) for coord in WindowCoords])
      if set(PairwiseCoordsDict.keys()) != \
      set([BamRefSeq.id,args.pairwise_align_to]):
        print('Malfunction of phylotypes: mismatch between the sequences',
        'found in the output of', TranslateCoordsCode, 'and the two names "' + \
        BamRefSeq.id+'", "'+args.pairwise_align_to +'". Quitting.',
        file=sys.stderr)
        exit(1)
      CoordsInRefs[BamRefSeq.id] = PairwiseCoordsDict[BamRefSeq.id]

  # We're creating a global alignment of all references:
  else:

    # Put all the mapping reference sequences into one file. If an alignment of 
    # other references was supplied, add the mapping references to that 
    # alignment; if not, align the mapping references to each other.
    SeqIO.write(RefSeqs, TempFileForRefs, "fasta")
    TempFiles.add(TempFileForRefs)
    if IncludeOtherRefs:
      FinalMafftOptions = ['--add', TempFileForRefs, args.alignment_of_other_refs]
    else:
      FinalMafftOptions = [TempFileForRefs]
    with open(FileForAlignedRefs, 'w') as f:
      try:
        ExitStatus = subprocess.call([args.x_mafft, '--quiet',
        '--preservecase'] + FinalMafftOptions, stdout=f)
        assert ExitStatus == 0
      except:
        print('Problem calling mafft. Quitting.', file=sys.stderr)
        raise

    if args.align_refs_only:
      print('References aligned in', FileForAlignedRefs+ \
      '. Quitting successfully.')
      exit(0)

    # If we're here and we're exploring window widths, we haven't defined the 
    # coordinates yet (because we haven't known the alignment length), unless
    # --ref-for-coords was specified.
    if ExploreWindowWidths and args.ref_for_coords == None:
      for seq in SeqIO.parse(open(FileForAlignedRefs),'fasta'):
        RefAlignmentLength = len(seq.seq)
        break
      MaxCoordForWindowWidthTesting = RefAlignmentLength
      WindowCoords = FindExploratoryWindows(MaxCoordForWindowWidthTesting)
      NumCoords = len(WindowCoords)
      UserCoords = WindowCoords

    # If window coords were specified with respect to one particular reference, 
    # or if we are excising certain coords, translate to alignment coords.
    if args.ref_for_coords != None or ExcisePositions:
      for seq in SeqIO.parse(open(FileForAlignedRefs),'fasta'):
        if seq.id == args.ref_for_coords:
          WindowCoords = \
          pf.TranslateSeqCoordsToAlnCoords(str(seq.seq), UserCoords)
        if seq.id == args.excision_ref:
          RefForExcisionGappySeq = str(seq.seq)
          AlignmentExcisionCoords = pf.TranslateSeqCoordsToAlnCoords(
          RefForExcisionGappySeq, args.excision_coords)


    # Determine windows automatically if desired
    if AutoWindows:
      command = [FindWindowsCode, FileForAlignedRefs,
      str(WeightedWindowWidth), str(WindowOverlap), '-S', str(WindowStartPos)]
      if NumAutoWindowParams == 4:
        command += ['-E', str(WindowEndPos)]
      try:
        WindowsString = subprocess.check_output(command)
      except:
        print('Problem executing', FindWindowsCode +'. Quitting.',
        file=sys.stderr)
        raise
      try:
        WindowCoords = [int(value) for value in WindowsString.split(',')]
        assert len(WindowCoords) >= 2
      except:
        print('Unable to understand the', FindWindowsCode, 'output -',
        WindowsString, '- as comma-separated integers. Quitting.',
        file=sys.stderr)
        raise
      try:
        NumCoords = SanityCheckWindowCoords(WindowCoords)
      except ValueError:
        print('Problematic output from ' +FindWindowsCode, file=sys.stderr)
        raise
      UserCoords = WindowCoords

    # Translate alignment coordinates to reference coordinates
    CoordsInRefs = TranslateCoords([FileForAlignedRefs, '-A']+\
    [str(coord) for coord in WindowCoords])

    # The index names in the CoordsInSeqs dicts, labelling the coords found by
    # coord translation, should cooincide with all seqs we're considering (i.e.
    # those in FileForAlignedRefs).
    if set(CoordsInRefs.keys()) != set(BamAliases+ExternalRefNames):
      print('Malfunction of phylotypes: mismatch between the sequences found',
      'in the output of', TranslateCoordsCode, 'and those in',
      FileForAlignedRefs +'. Quitting.', file=sys.stderr)
      exit(1)

# Make index files for the bam files if needed.
pf.MakeBamIndices(BamFiles, args.x_samtools)

# Gather some data from each bam file
BamFileRefSeqNames = {}
BamFileRefLengths  = {}
if args.verbose:
  print('Now preparing the bam files for analysis.')
for i,BamFileName in enumerate(BamFiles):

  BamFileBasename = BamFileBasenames[i]
  BamAlias = BamAliases[i]

  # Prep for pysam
  try:
    BamFile = pysam.AlignmentFile(BamFileName, "rb")
  except AttributeError:
    print('Error calling "pysam.AlignmentFile". The AlignmentFile attribute',
    'was introduced in pysam version 0.8.1; are you using an older version',
    'than that? You might be able to update by running\npip install pysam',
    '--upgrade\nfrom the command line. Quitting.', file=sys.stderr)
    exit(1)
  except ValueError:
    print('Error trying to read', BamFileName, 'as a bam file with pysam.',
    'Quitting.', file=sys.stderr)
    raise


  # Find the reference in the bam file; there should only be one.
  AllReferences = BamFile.references
  if len(AllReferences) != 1:
    print('Expected exactly one reference in', BamFileName+'; found',
    str(len(AllReferences))+'.\nQuitting.', file=sys.stderr)
    exit(1)
  BamFileRefSeqNames[BamFileBasename] = AllReferences[0]

  # Get the length of the reference.
  AllReferenceLengths = BamFile.lengths
  if len(AllReferenceLengths) != 1:
    print('Pysam error: found one reference but', len(AllReferenceLengths),
    'reference lengths.\nQuitting.', file=sys.stderr)
    exit(1)
  RefLength = AllReferenceLengths[0]
  BamFileRefLengths[BamFileBasename] = RefLength

  # When translating coordinates, -1 means before the sequence starts; 'NaN'
  # means after it ends. These should be replaced by 1 and the reference length
  # respectively.
  for j,coord in enumerate(CoordsInRefs[BamAlias]):
    if coord == -1:
      CoordsInRefs[BamAlias][j] = 1
    elif coord == 'NaN':
      CoordsInRefs[BamAlias][j] = RefLength

def ProcessReadDict(ReadDict, WhichBam, LeftWindowEdge, RightWindowEdge,
WindowAsStr):
  '''Turns a dict of reads into a list of reads, merging & imposing a minimum
  count.'''

  # For naming things
  BamFileBasename = BamFileBasenames[WhichBam]
  BasenameForReads = BamAliases[WhichBam]

  # Merge similar reads if desired
  if MergeReadsA:
    ReadDict = pf.MergeSimilarStringsA(ReadDict, MergingThreshold)
  if MergeReadsB:
    ReadDict = pf.MergeSimilarStringsB(ReadDict, MergingThreshold)


  # Implement the minimum read count
  if args.min_read_count > 1:
    ReadDict = {read:count for read, count in ReadDict.items() if \
    count >= args.min_read_count}

  # Warn if there are no reads
  if len(ReadDict) == 0 and (not ExploreWindowWidths):
    print('Warning: bam file ', BamFileBasename, ' has no reads (after ',
    'processing) that fully span the window ', WindowAsStr, ".", sep='',
    file=sys.stderr)
    return []

  # Return a list of reads named according to their count.
  reads = []
  for k, (read, count) in \
  enumerate(sorted(ReadDict.items(), key=lambda x: x[1], reverse=True)):
    SeqName = BasenameForReads+'_read_'+str(k+1)+'_count_'+str(count)
    SeqObject = SeqIO.SeqRecord(Seq.Seq(read), id=SeqName, description='')
    reads.append(SeqObject)
  return reads

# This regex matches "_read_" then any integer then "_count_" then any integer,
# constrained to come at the end of the string. We'll need it later.
SampleRegex = re.compile('_read_\d+_count_\d+$')

def ReadAlignedReadsIntoDicts(AlignIOobject, ValuesAreCounts=True):
  '''Collects sample seqs and into dicts, and other seqs into a list.

  The values of the dicts are either the seq count (inferred from the seq name)
  or simply the seq name.'''
  SampleReadCounts = collections.OrderedDict()
  NonSampleSeqs = []
  for seq in AlignIOobject:
    RegexMatch = SampleRegex.search(seq.id)
    if RegexMatch and seq.id[:RegexMatch.start()] in BamAliases:
      SampleName = seq.id[:RegexMatch.start()]
      read = str(seq.seq)
      if ValuesAreCounts:
        value = int(seq.id.rsplit('_',1)[1])
      else:
        value = seq.id
      if SampleName in SampleReadCounts:
        # After excising positions, a sample can have the same read twice:
        if read in SampleReadCounts[SampleName]:
          SampleReadCounts[SampleName][read] += value
        else:
          SampleReadCounts[SampleName][read] = value
      else:
        SampleReadCounts[SampleName] = {read : value}
    else:
      assert seq.id in ExternalRefNames, 'Malfunction of phylotypes: '+\
      'sequence ' + seq.id + ' is not recognised as a read nor as an external'+\
      ' reference.'
      NonSampleSeqs.append(seq)

  return SampleReadCounts, NonSampleSeqs

def RemovePureGapCols(alignment):
  "Removes pure-gap columns from an alignment."
  AlignmentLength = alignment.get_alignment_length()
  for column in reversed(xrange(AlignmentLength)):
    PureGap = True
    for base in alignment[:, column]:
      if base != GapChar:
        PureGap = False
        break
    if PureGap:
      alignment = alignment[:, :column] + alignment[:, column+1:]
  return alignment

def ReMergeAlignedReads(alignment, ForceNoMerging=False):
  '''Splits an alignment object into reads and refs, re-merges the reads,
  renames them, and removes pure-gap columns.'''

  SampleReadCounts, RefSeqsHere = ReadAlignedReadsIntoDicts(alignment)
  NewAlignment = AlignIO.MultipleSeqAlignment([])
  for SampleName in SampleReadCounts:
    if not ForceNoMerging:
      if MergeReadsA:
        SampleReadCounts[SampleName] = \
        pf.MergeSimilarStringsA(SampleReadCounts[SampleName], MergingThreshold)
      if MergeReadsB:
        SampleReadCounts[SampleName] = \
        pf.MergeSimilarStringsB(SampleReadCounts[SampleName], MergingThreshold)
    for k, (read, count) in enumerate(sorted(
    SampleReadCounts[SampleName].items(), key=lambda x: x[1], reverse=True)):
      ID = SampleName+'_read_'+str(k+1)+'_count_'+str(count)
      SeqObject = SeqIO.SeqRecord(Seq.Seq(read), id=ID, description='')
      NewAlignment.append(SeqObject)
  NewAlignment.extend(RefSeqsHere)

  # Merging after alignment means some columns could be pure gap. Remove these.
  if MergeReads:
    NewAlignment = RemovePureGapCols(NewAlignment)
  return NewAlignment

def FindPatientsConsensuses(alignment):
  '''Finds the consensus sequence for each patient appearing in an alignment.'''
  SampleReadCounts, RefSeqsHere = ReadAlignedReadsIntoDicts(alignment)
  ConsensusAlignment = AlignIO.MultipleSeqAlignment([])
  AlignmentLength = alignment.get_alignment_length()
  for SampleName, ReadsAndCounts in SampleReadCounts.items():

    # Count each base seen at each position
    BaseCounterDicts = [{} for pos in range(0,AlignmentLength)]
    TotalCount = 0
    for read, count in ReadsAndCounts.items():
      TotalCount += count
      for pos, base in enumerate(read):
        if base in BaseCounterDicts[pos]:
          BaseCounterDicts[pos][base] += count
        else:
          BaseCounterDicts[pos][base] = count

    # Find the most common 'base' (could be a gap) at each position.
    consensus = ''
    for pos, BaseCounterDict in enumerate(BaseCounterDicts):
      MostCommonBase = None
      HighestCount = 0
      for base, count in BaseCounterDict.items():
        if count > HighestCount:
          MostCommonBase = base
          HighestCount = count
      assert MostCommonBase != None, 'Problem for ' + SampleName + \
      ' at position ' + str(pos)
      consensus += MostCommonBase
    SeqObject = SeqIO.SeqRecord(Seq.Seq(consensus), id=SampleName + \
    '_count_' + str(TotalCount), description='')
    ConsensusAlignment.append(SeqObject)
  for ref in RefSeqsHere:
    ConsensusAlignment.append(ref)
  #return RemovePureGapCols(ConsensusAlignment)
  return ConsensusAlignment

# If we're keeping track list of discarded read pairs for each bam file:
if args.inspect_disagreeing_overlaps:
  DiscardedReadPairsDict = \
  {BamFileBasename:[] for BamFileBasename in BamFileBasenames}

if args.time:
  times.append(time.time())
  LastStepTime = times[-1] - times[-2]
  print('Bam and Reference pre-processing finished. Number of seconds taken:',
  LastStepTime)

AllPatientsReadNamesInThisWindow = {BamFile:set() for BamFile in BamFiles}
ThisWindow = (float('-Inf'), float('-Inf'))

# Subdirectories for output files.
OutputDirs = {}
OutputFilesByDestinationDir = {}
OutputDirs['raxml'] = 'RAxMLfiles'

FileForAlignedReads_basename = 'AlignedReads'
FileForAlignedReads_PositionsExcised_basename = FileForAlignedReads_basename + \
'_PositionsExcised_'
OutputDirs['AlignedReads'] = 'AlignedReads'
OutputFilesByDestinationDir['AlignedReads'] = []

FileForConsensuses_basename = 'Consensuses_'
FileForConsensuses_PositionsExcised_basename = FileForConsensuses_basename + \
'PositionsExcised_'
OutputDirs['Consensuses'] = 'Consensuses'
OutputFilesByDestinationDir['Consensuses'] = []

FileForDiscardedReadPairs_basename = 'DiscardedReads_'
OutputDirs['DiscardedReads'] = 'DiscardedReads'
OutputFilesByDestinationDir['DiscardedReads'] = []

FileForDuplicateReadCountsRaw_basename = 'DuplicateReadCountsRaw_'
FileForDuplicateReadCountsProcessed_basename = 'DuplicateReadCountsProcessed_'
FileForDuplicateSeqs_basename = 'DuplicateReads_contaminants_'
OutputDirs['DupData'] = 'DuplicationData'
OutputFilesByDestinationDir['DupData'] = []

FileForReadNames1_basename = 'ReadNames1_'
FileForReadNames2_basename = 'ReadNames2_'
OutputDirs['ReadNames'] = 'ReadNames'
OutputFilesByDestinationDir['ReadNames'] = []

FileForRecombinantReads_basename = 'RecombinantReads_'
OutputDirs['RecombFiles'] = 'RecombinationData'
OutputFilesByDestinationDir['RecombFiles'] = []

NumMLtreesMade = 0

HaveWarnedNoQualities = False

# Iterate through the windows
for window in range(NumCoords / 2):

  # If coords were specified with respect to one particular reference,
  # WindowCoords is the translation of those coords to alignment coordinates.
  # UserCoords are the original coords, which we use for labelling things to
  # keep labels intuitive for the user.
  UserLeftWindowEdge  = UserCoords[window*2]
  UserRightWindowEdge = UserCoords[window*2 +1]
  ThisWindowSuffix = 'InWindow_'+str(UserLeftWindowEdge)+'_to_'+\
  str(UserRightWindowEdge)
  ThisWindowAsStr = str(UserLeftWindowEdge) + '-' + str(UserRightWindowEdge)

  if PrintInfo:
    print('Now extracting and processing reads in window', ThisWindowAsStr)

  # Prepare some things for checking for reads appearing again the consecutive
  # overlapping windows.
  AllPatientsReadNamesInLastWindow = AllPatientsReadNamesInThisWindow
  AllPatientsReadNamesInThisWindow = {BamFile:set() for BamFile in BamFiles}
  LastWindow = ThisWindow
  ThisWindow = (UserLeftWindowEdge, UserRightWindowEdge)
  OverlapsLastWindow = (LastWindow[0] <= ThisWindow[0] <= LastWindow[1]) or \
                       (ThisWindow[0] <= LastWindow[0] <= ThisWindow[1])

  # Get ready to record reads here from all samples
  AllReadsInThisWindow = []
  if CheckDuplicates:
    AllReadDictsInThisWindow = []

  # Try to find a contamination file for this window. If there is none, that 
  # could be because the user did not put it in the intended directory, or
  # because no contamination was found for this window: warn, but proceed. If
  # there is one, read it in.
  ContaminantReadsInput = {}
  if RecallContaminants:
    try:
      ContaminantFile = \
      ContaminantFilesByWindow[(UserLeftWindowEdge, UserRightWindowEdge)]
    except KeyError:
      print('Warning: no contaminant file found for window', ThisWindowAsStr + \
      '.', file=sys.stderr)
    else:
      for seq in SeqIO.parse(open(ContaminantFile), 'fasta'):
        if seq.id in ContaminantReadsInput:
          ContaminantReadsInput[seq.id].append(str(seq.seq))
        else:
          ContaminantReadsInput[seq.id] = [str(seq.seq)]

  # Iterate through the bam files
  for i,BamFileName in enumerate(BamFiles):

    # Recall some things we've already worked out for this bam file and stored.
    BamFileBasename = BamFileBasenames[i]
    RefSeqName = BamFileRefSeqNames[BamFileBasename]
    RefLength = BamFileRefLengths[BamFileBasename]
    BamAlias = BamAliases[i]
    ThisBamCoords = CoordsInRefs[BamAlias]
    LeftWindowEdge  = ThisBamCoords[window*2]
    RightWindowEdge = ThisBamCoords[window*2 +1]

    if args.verbose:
      print('Now extracting & processing reads from bam', BamAlias + '.')

    # For labelling read name files
    FileForReadNames1_basename_ThisBam = FileForReadNames1_basename + \
    ThisWindowSuffix + '_InBam_'
    FileForReadNames2_basename_ThisBam = FileForReadNames2_basename + \
    ThisWindowSuffix + '_InBam_'

    # Pysam uses zero-based coordinates for positions w.r.t the reference.
    # If we want all reads that start exactly at the window start and end
    # anywhere after, or all reads that end exactly at the window end and start
    # anywhere before, set end=start or start=end respectively, to make sure
    # pysam's fetch function retrieves all the reads we need.
    LeftWindowEdge  = LeftWindowEdge  -1
    RightWindowEdge = RightWindowEdge -1
    LeftWindowEdgeForFetch = LeftWindowEdge
    RightWindowEdgeForFetch = RightWindowEdge
    if args.exact_window_start:
      if not args.exact_window_end:
        RightWindowEdgeForFetch = None
        RightWindowEdge = LeftWindowEdge
    elif args.exact_window_end:
      LeftWindowEdgeForFetch = None
      LeftWindowEdge = RightWindowEdge

    # Find all unique reads in this window and count their occurrences.
    AllReads = {}
    UniqueReads = {}
    ReadNames = []
    ReadNameDict = {}
    BamFile = pysam.AlignmentFile(BamFileName, "rb")
    for read in BamFile.fetch(RefSeqName, LeftWindowEdgeForFetch,
    RightWindowEdgeForFetch):

      # fetch isn't supposed to return unmapped reads, but does.
      if read.is_unmapped:
        continue

      # Skip improperly paired reads if desired
      if args.discard_improper_pairs and read.is_paired and \
      not read.is_proper_pair:
        continue

      if not read.query_qualities:
        read.query_qualities = [106 for base in range(len(read.query_sequence))]
        if not HaveWarnedNoQualities:
          print('Warning: found a read with no information about base',
          'qualities. All bases will be set to have quality 106 for this read',
          'and any others lacking this information found henceforth.',
          file=sys.stderr)
          if QualTrimEnds or ImposeMinQual: 
            print('WARNING: you have specified at least one option relating to',
            'read quality, when read quality information is missing. This is',
            'strongly discouraged. Continuing nonetheless.', file=sys.stderr)
          HaveWarnedNoQualities = True

      if args.merge_paired_reads:

        # If we've seen this read's mate already, merge the pair.
        if read.query_name in AllReads:
          Read1 = AllReads[read.query_name]
          try:
            Read1asPseudoRead = pf.PseudoRead.InitFromRead(Read1)
          except AttributeError:
            # An attribute error will arise if we encounter the same name three
            # times; check if that's the issue and print a more helpful exit
            # message, otherwise just raise.
            ReadCounts = {}
            for NewRead in BamFile.fetch(RefSeqName):
              if NewRead.query_name in ReadCounts:
                if ReadCounts[NewRead.query_name] == 2:
                  print('The name', NewRead.query_name, 'occurs (at least) 3',
                  'times in', BamFileBasename + '; this should never happen -',
                  'the same name should only be found for the two reads in a',
                  'pair. Quitting.', file=sys.stderr)
                  exit(1)
                else:
                  ReadCounts[NewRead.query_name] += 1
              else:
                ReadCounts[NewRead.query_name] = 1
            print('Encountered an error related to converting reads from pysam',
            "format to phyloscanner's format, for this read:", read,
            "\nPlease report to Chris Wymant. Quitting.", file=sys.stderr)
            raise
          Read2 = read
          Read2asPseudoRead = pf.PseudoRead.InitFromRead(read)
          MergedRead = Read1asPseudoRead.MergeReadPairOverWindow(
          Read2asPseudoRead, LeftWindowEdge, RightWindowEdge,
          args.quality_trim_ends, args.min_internal_quality,
          args.recover_clipped_ends)
          if MergedRead == None:
            del AllReads[read.query_name]
            continue
          elif MergedRead == False:
            del AllReads[read.query_name]
            if args.inspect_disagreeing_overlaps:
              DiscardedReadPairsDict[BamFileBasename] += [Read1,Read2]
            continue
          AllReads[read.query_name] = MergedRead

        # We've not come across a read with this name before. Record & move on.
        # Note that we need to save the read, rather than the pseudoread, in
        # case
        else:
          AllReads[read.query_name] = read

      # If we're not merging paired reads, process this read now to save memory.
      # ProcessRead returns None if we don't want to consider this read.
      else:
        ReadAsPseudoRead = pf.PseudoRead.InitFromRead(read)
        seq = ReadAsPseudoRead.ProcessRead(LeftWindowEdge, RightWindowEdge,
          args.quality_trim_ends, args.min_internal_quality,
          args.keep_overhangs, args.recover_clipped_ends,
          args.exact_window_start, args.exact_window_end)
        if seq == None:
          continue

        # We're not merging read pairs here, so we could see the same read name
        # more than once in the same window, in which case skip it. Then, add
        # the read name to this window's list, even if we don't use it
        # because it was in the last window. Otherwise, if a read was in three
        # consecutive windows, we'd skip it in the second and think we were OK
        # to use it again in the third.
        if args.forbid_read_repeats:
          if read.query_name in AllPatientsReadNamesInThisWindow[BamFileName]:
            continue
          AllPatientsReadNamesInThisWindow[BamFileName].add(read.query_name)
          if OverlapsLastWindow and \
          read.query_name in AllPatientsReadNamesInLastWindow[BamFileName]:
            continue

        if seq in UniqueReads:
          UniqueReads[seq] += 1
        else:
          UniqueReads[seq] = 1

        # Record the read name if desired.
        if args.read_names_1:
          ReadNames.append(read.query_name)
        if args.read_names_2:
          if seq in ReadNameDict:
            ReadNameDict[seq].append(read.query_name)
          else:
            ReadNameDict[seq] = [read.query_name]


    # If we did merge paired reads, we now need to process them.
    # AllReads will be a mixture of PseudoRead instances (for merged read pairs)
    # and pysam.AlignedSegment instances (for unmerged single reads). The latter
    # must be converted to PseudoRead instances to be processed.
    if args.merge_paired_reads:
      for read in AllReads.values():
        try:
          seq = read.ProcessRead(LeftWindowEdge, RightWindowEdge,
          args.quality_trim_ends, args.min_internal_quality,
          args.keep_overhangs, args.recover_clipped_ends,
          args.exact_window_start, args.exact_window_end)
        except AttributeError:
          #print(type(read))
          ReadAsPseudoRead = pf.PseudoRead.InitFromRead(read)          
          seq = ReadAsPseudoRead.ProcessRead(LeftWindowEdge, RightWindowEdge,
          args.quality_trim_ends, args.min_internal_quality,
          args.keep_overhangs, args.recover_clipped_ends,
          args.exact_window_start, args.exact_window_end)
          ReadName = read.query_name
        else:
          ReadName = read.name
        if seq == None:
          continue

        # Check if we've seen this merged read pair in the last window.
        if args.forbid_read_repeats:
          AllPatientsReadNamesInThisWindow[BamFileName].add(ReadName)
          if OverlapsLastWindow and \
          ReadName in AllPatientsReadNamesInLastWindow[BamFileName]:
            continue

        if seq in UniqueReads:
          UniqueReads[seq] += 1
        else:
          UniqueReads[seq] = 1

        # Record the read name if desired.
        if args.read_names_1:
          ReadNames.append(ReadName)
        if args.read_names_2:
          if seq in ReadNameDict:
            ReadNameDict[seq].append(ReadName)
          else:
            ReadNameDict[seq] = [ReadName]

    # If we've read in any contaminant reads for this window and this bam,
    # remove them from the read dict. If they're not present in the read dict,
    # warn but proceed. A 1bp slip in alignments between the previous
    # (contaminant-finding) run and the current run could cause such an issue.
    if BamAlias in ContaminantReadsInput:
      HaveWarned = False
      for read in ContaminantReadsInput[BamAlias]:
        try:
          del UniqueReads[read]
        except KeyError:
          if not HaveWarned:
            print('Warning: at least one contaminant read in', ContaminantFile,
            'from', BamAlias, 'was not found in this window in',
            BamFileBasename + '. This could be due to a mismatch in window',
            'coordinates between the run that generated that contamination'
            'file and the present run. Proceeding.', file=sys.stderr)
            HaveWarned = True

    # If we are checking for read duplication between samples, record the file 
    # name and read dict for this sample and move on to the next sample.
    if CheckDuplicates:
      AllReadDictsInThisWindow.append((BamAlias, UniqueReads,
      LeftWindowEdge, RightWindowEdge))

    # If we're not checking for read duplication between samples, process the
    # read dict for this sample now and add it to the list of all reads here.
    else:
      AllReadsInThisWindow += \
      ProcessReadDict(UniqueReads, i, LeftWindowEdge, RightWindowEdge,
      ThisWindowAsStr)

    # Write recorded read names to file if desired.
    if args.read_names_1:
      FileForReadNames1 = FileForReadNames1_basename_ThisBam + BamAlias + '.txt'
      with open(FileForReadNames1, 'w') as f:
        f.write('\n'.join(ReadNames))
      OutputFilesByDestinationDir['ReadNames'].append(FileForReadNames1)
    if args.read_names_2:
      FileForReadNames2 = FileForReadNames2_basename_ThisBam + BamAlias + '.csv'
      with open(FileForReadNames2, 'w') as f:
        for seq, ReadNamesForThatSeq in ReadNameDict.items():
          f.write(seq + ',' + ' '.join(ReadNamesForThatSeq) + '\n')
      OutputFilesByDestinationDir['ReadNames'].append(FileForReadNames2)
  if args.read_names_only:
    continue

  # We've now gathered together reads from all bam files for this window.

  # If we're checking for duplicate reads between samples, do so now.
  # Check every dict against every other dict, and record the ratio of counts
  # for any shared reads.
  if CheckDuplicates:
    if args.verbose:
      print('Now checking for duplication of reads between bam files.')
    DuplicateDetails = []
    ContaminantReadsFound = {}
    for i, (BamFile1Alias, ReadDict1, LeftWindowEdge1, RightWindowEdge1) \
    in enumerate(AllReadDictsInThisWindow):
      for j, (BamFile2Alias, ReadDict2, LeftWindowEdge2, RightWindowEdge2) \
      in enumerate(AllReadDictsInThisWindow[i+1:]):
        DuplicateReadRatios = []
        for read in ReadDict1:
          if read in ReadDict2:
            Bam1Count = ReadDict1[read]
            Bam2Count = ReadDict2[read]
            DuplicateDetails.append(
            (BamFile1Alias, BamFile2Alias, Bam1Count, Bam2Count))

            # Diagnose contaminants
            if FlagContaminants:
              CountRatio = float(Bam1Count) / Bam2Count
              ContaminantAlias = None
              if CountRatio >= args.contaminant_count_ratio:
                ContaminantAlias = BamFile2Alias
              elif CountRatio <= 1. / args.contaminant_count_ratio:
                ContaminantAlias = BamFile1Alias
              if ContaminantAlias != None:
                if ContaminantAlias in ContaminantReadsFound:
                  # It's possible this read for this patient is considered
                  # contamination from more than one source, so check the read
                  # isn't there already before adding it to the list:
                  if not read in ContaminantReadsFound[ContaminantAlias]:
                    ContaminantReadsFound[ContaminantAlias].append(read)
                else:
                  ContaminantReadsFound[ContaminantAlias] = [read]

    if DuplicateDetails != []:
      FileForDuplicateReadCountsRaw = FileForDuplicateReadCountsRaw_basename + \
      ThisWindowSuffix + '.csv'
      with open(FileForDuplicateReadCountsRaw, 'w') as f:
        f.write('"Alias1","Alias2","Count1","Count2"\n')
        f.write('\n'.join(','.join(map(str,data)) for data in DuplicateDetails))
      OutputFilesByDestinationDir['DupData'].append(FileForDuplicateReadCountsRaw)

    # If contaminants are diagnosed, print them and remove them from their
    # ReadDict.
    if ContaminantReadsFound != {}:
      FileForDuplicateSeqs = FileForDuplicateSeqs_basename + \
      ThisWindowSuffix + '.fasta'
      AllContaminants = []
      for alias, reads in ContaminantReadsFound.items():
        for read in reads:
          AllContaminants.append(SeqIO.SeqRecord(Seq.Seq(read), id=alias,
          description=''))
      SeqIO.write(AllContaminants, FileForDuplicateSeqs, "fasta")
      OutputFilesByDestinationDir['DupData'].append(FileForDuplicateSeqs)
      for i, (BamAlias, ReadDict, LeftWindowEdge, RightWindowEdge) \
      in enumerate(AllReadDictsInThisWindow):
        if BamAlias in ContaminantReadsFound:
          for read in ContaminantReadsFound[BamAlias]:
            del AllReadDictsInThisWindow[i][1][read]
    if args.flag_contaminants_only:
      continue

    # Process the read dicts (not yet done if we're checking for duplicates).
    for i, (BamFileBasename, ReadDict, LeftWindowEdge, RightWindowEdge) \
    in enumerate(AllReadDictsInThisWindow):
      AllReadsInThisWindow += \
      ProcessReadDict(ReadDict, i, LeftWindowEdge, RightWindowEdge,
      ThisWindowAsStr)

  # All read dicts have now been processed into the list AllReadsInThisWindow.

  # Skip empty windows.
  if AllReadsInThisWindow == []:
    if ExploreWindowWidths:
      for alias in BamAliases:
        WindowWidthExplorationData.append([UserLeftWindowEdge,
        UserRightWindowEdge, alias, 0])
    else:
      print('WARNING: no bam file had any reads (after a minimum post-merging '+\
      'read count of', args.min_read_count, 'was imposed) in the window',
      ThisWindowAsStr + '. Skipping to the next window.', file=sys.stderr)
    continue

  # Re-define the window edge coords to be with respect to the alignment of refs
  # rather than a bam file.
  LeftWindowEdge  = WindowCoords[window*2]
  RightWindowEdge = WindowCoords[window*2 +1]

  # Create a fasta file with all reads in this window, ready for aligning.
  # If there's only one, we don't need to align (or make trees!).
  TempFileForReadsHere = TempFileForReads_basename + ThisWindowSuffix+\
  '.fasta'
  FileForAlnReadsHere = FileForAlignedReads_basename + \
  ThisWindowSuffix +'.fasta'
  if len(AllReadsInThisWindow) == 1 and not IncludeOtherRefs:
    SeqIO.write(AllReadsInThisWindow, FileForAlnReadsHere, "fasta")
    OutputFilesByDestinationDir['AlignedReads'].append(FileForAlnReadsHere)
    # If we're exploring window widths, record that all bams but one have no
    # reads.
    if ExploreWindowWidths:
      TheReadID = AllReadsInThisWindow[0].id
      RegexMatch = SampleRegex.search(TheReadID)
      if RegexMatch and TheReadID[:RegexMatch.start()] in BamAliases:
        TheBamWithOneRead = TheReadID[:RegexMatch.start()]
      else:
        print('Malfunction of phylotypes: there is only one read in the ',
        'window ', ThisWindowAsStr, ', namely ', TheReadID,
        ", but we can't figure out which bam we got it from. Quitting.", sep='',
        file=sys.stderr)
        exit(1)
      for alias in BamAliases:
        if alias == TheBamWithOneRead:
          count = 1
        else:
          count = 0
        WindowWidthExplorationData.append([UserLeftWindowEdge,
        UserRightWindowEdge, alias, count])
    else:
      if PrintInfo:
        print('There is only one read in this window, written to ' +\
        FileForAlnReadsHere +'. Skipping to the next window.')
    continue
  SeqIO.write(AllReadsInThisWindow, TempFileForReadsHere, "fasta")
  TempFiles.add(TempFileForReadsHere)
  OutputFilesByDestinationDir['AlignedReads'].append(FileForAlnReadsHere)
  FileForTrees = FileForAlnReadsHere

  # If external refs are included, find the part of each one's seq corresponding
  # to this window and put them all in another file.
  # If we did pairwise aligning of refs, we know the coordinates we want in the
  # ExternalRefAlignment object. If we did a global alignment, we slice the
  # desired window out of that alignment.
  if IncludeOtherRefs:
    TempFileForOtherRefsHere = TempFileForOtherRefs_basename + \
    ThisWindowSuffix +'.fasta'
    if PairwiseAlign:
      ExternalRefLeftWindowEdge  = ExternalRefWindowCoords[window*2]
      ExternalRefRightWindowEdge = ExternalRefWindowCoords[window*2 +1]
      RefAlignmentInWindow = ExternalRefAlignment[:,
      ExternalRefLeftWindowEdge-1:ExternalRefRightWindowEdge]
      RefsThatAreNotPureGap = []
      for seq in RefAlignmentInWindow:
        if len(seq.seq.ungap(GapChar)) != 0:
          RefsThatAreNotPureGap.append(seq)
      if len(RefsThatAreNotPureGap) == 0:
        print('Error: all external references are pure gap in this window;',
        'skipping to the next window.', file=sys.stderr)
        continue
      AlignIO.write(Align.MultipleSeqAlignment(RefsThatAreNotPureGap),
      TempFileForOtherRefsHere, 'fasta')
      TempFiles.add(TempFileForOtherRefsHere)
    else:
      with open(TempFileForOtherRefsHere, 'w') as f:
        try:
          ExitStatus = subprocess.call([FindSeqsInFastaCode,
          FileForAlignedRefs, '-B', '-W', str(LeftWindowEdge) + ',' + \
          str(RightWindowEdge), '-v'] + BamAliases, stdout=f)
          assert ExitStatus == 0
        except:
          print('Problem calling', FindSeqsInFastaCode+\
          '. Skipping to the next window.', file=sys.stderr)
          continue

  # Update on time taken if desired
  if args.time:
    times.append(time.time())
    LastStepTime = times[-1] - times[-2]
    print('Read pre-processing in window', ThisWindowAsStr,
    'finished. Number of seconds taken: ', LastStepTime)

  # Align the reads. Prepend 'temp_' to the file name if we'll merge again after
  # aligning.
  if MergeReads:
    FileForReads = 'temp_' + FileForAlnReadsHere
    TempFiles.add(FileForReads)
  else:
    FileForReads = FileForAlnReadsHere
  if IncludeOtherRefs:
    FinalMafftOptions = ['--add', TempFileForReadsHere, TempFileForOtherRefsHere]
  else:
    FinalMafftOptions = [TempFileForReadsHere]
  with open(FileForReads, 'w') as f:
    try:
      ExitStatus = subprocess.call([args.x_mafft, '--quiet', '--preservecase']+\
      FinalMafftOptions, stdout=f)
      assert ExitStatus == 0
    except:
      print('Problem calling mafft. Skipping to the next window.',
      file=sys.stderr)
      continue
    if not os.path.isfile(FileForReads):
      print('Error:', FileForReads +', expected to be produced by mafft, does',
      'not exist. Skipping to the next window.', file=sys.stderr)
      continue


  # Update on time taken if desired
  if args.time:
    times.append(time.time())
    LastStepTime = times[-1] - times[-2]
    print('Read alignment in window', ThisWindowAsStr,
    'finished. Number of seconds taken: ', LastStepTime)

  # Read in the aligned reads.
  try:
    SeqAlignmentHere = AlignIO.read(FileForReads, "fasta")
  except:
    print('Malfunction of phylotypes: problem encountered reading in',
    FileForReads, 'as an alignment. Quitting.', file=sys.stderr)
    raise

  # Do a second round of within-sample read merging now the reads are aligned.
  # Write the output to FileForAlnReadsHere.
  if MergeReads:
    try:
      SeqAlignmentHere = ReMergeAlignedReads(SeqAlignmentHere)
    except:
      print('Problem encountered while analysing', FileForReads +'. Quitting.',
      file=sys.stderr)
      raise
    AlignIO.write(SeqAlignmentHere, FileForAlnReadsHere, 'fasta')

  # Find & write the consensuses.
  ConsensusAlignment = FindPatientsConsensuses(SeqAlignmentHere)
  FileForConsensuses = FileForConsensuses_basename + ThisWindowSuffix +'.fasta'
  AlignIO.write(ConsensusAlignment, FileForConsensuses, 'fasta')
  OutputFilesByDestinationDir['Consensuses'].append(FileForConsensuses)

  # See if there are positions to excise in this window.
  if ExcisePositions:
    FileForAlignedReads_PositionsExcised = \
    FileForAlignedReads_PositionsExcised_basename + ThisWindowSuffix +'.fasta'
    if PairwiseAlign:
      CoordsToExciseInThisWindow = [coord for coord in args.excision_coords \
      if LeftWindowEdge <= coord <= RightWindowEdge]
    else:
      CoordsToExciseInThisWindow = [coord for coord in AlignmentExcisionCoords \
      if LeftWindowEdge <= coord <= RightWindowEdge]
    if CoordsToExciseInThisWindow != []:

      # Define PositionsInUngappedRef to be how far the positions are from the
      # start of the window, in an ungapped version of the ref.
      if PairwiseAlign:
        PositionsInUngappedRef = \
        [coord - LeftWindowEdge + 1 for coord in CoordsToExciseInThisWindow]
        UngappedRefHere = \
        str(RefForPairwiseAlns.seq)[LeftWindowEdge-1:RightWindowEdge]
      else:
        RefInThisWindowGappy = \
        RefForExcisionGappySeq[LeftWindowEdge-1:RightWindowEdge]
        PositionsInUngappedRef = []
        for coord in CoordsToExciseInThisWindow:
          DistanceIntoWindow = coord - LeftWindowEdge
          PositionsInUngappedRef.append(
          len(RefInThisWindowGappy[:DistanceIntoWindow+1].replace(GapChar,'')))
        UngappedRefHere = RefInThisWindowGappy.replace(GapChar,'')

      # Check the ref looks as expected.
      RefInAlignment = None
      for seq in SeqAlignmentHere:
        if seq.id == args.excision_ref:
          RefInAlignment = str(seq.seq)
          break
      if RefInAlignment == None:
        print('Malfunction of phylotypes: unable to find', args.excision_ref,
        'in', FileForAlnReadsHere +'. Quitting.', file=sys.stderr)
        exit(1)
      if RefInAlignment.replace(GapChar,'') != UngappedRefHere:
        print('Malfunction of phylotypes: mismatch between the ref for',
        'excision we expected to find in this window:\n', UngappedRefHere,
        '\nand the ref for excision we actually found in this window:\n',
        RefInAlignment.replace(GapChar,''), '\nQuitting.', file=sys.stderr)
        exit(1)

      # Excise the positions in the aligned set of reads.
      PositionsInAlignment = \
      pf.TranslateSeqCoordsToAlnCoords(RefInAlignment, PositionsInUngappedRef)
      assert PositionsInAlignment == sorted(PositionsInAlignment, reverse=True)
      for pos in PositionsInAlignment:
        SeqAlignmentHere = \
        SeqAlignmentHere[:, :pos-1] + SeqAlignmentHere[:, pos:]

      # Excising positions may have made some sequences identical within a
      # sample, which need to be merged even if the merging parameter is 0.
      try:
        SeqAlignmentHere = ReMergeAlignedReads(SeqAlignmentHere,
        ForceNoMerging=True)
      except:
        print('Problem encountered while analysing', FileForAlnReadsHere + \
        '. Quitting.', file=sys.stderr)
        raise
      AlignIO.write(SeqAlignmentHere, FileForAlignedReads_PositionsExcised,
      'fasta')
      OutputFilesByDestinationDir['AlignedReads'].append(
      FileForAlignedReads_PositionsExcised)
      FileForTrees = FileForAlignedReads_PositionsExcised

      # Find consensuses again after excising positions:
      ConsensusAlignment = FindPatientsConsensuses(SeqAlignmentHere)
      FileForConsensuses_PositionsExcised = \
      FileForConsensuses_PositionsExcised_basename + ThisWindowSuffix +'.fasta'
      AlignIO.write(ConsensusAlignment, FileForConsensuses_PositionsExcised,
      'fasta')
      OutputFilesByDestinationDir['Consensuses'].append(
      FileForConsensuses_PositionsExcised)

  # If we're exploring window widths, we just care how many unique reads
  # were found here. Record & move on.
  if ExploreWindowWidths:
    NumUniqueReadsPerPatient = {alias : 0 for alias in BamAliases}
    for seq in SeqAlignmentHere:
      RegexMatch = SampleRegex.search(seq.id)
      if RegexMatch and seq.id[:RegexMatch.start()] in BamAliases:
        SampleName = seq.id[:RegexMatch.start()]
        NumUniqueReadsPerPatient[SampleName] += 1
    for alias, count in NumUniqueReadsPerPatient.items():
      WindowWidthExplorationData.append([UserLeftWindowEdge,
      UserRightWindowEdge, alias, count])
    continue

  if CheckDuplicates:

    # Find any duplicates
    if len(SeqAlignmentHere) > 1:
      SeqToIDsDict = {}
      DuplicatesDict = {}
      for seq in SeqAlignmentHere:
        SeqAsStr = str(seq.seq) 
        if SeqAsStr in SeqToIDsDict:
          if SeqAsStr in DuplicatesDict:
            DuplicatesDict[SeqAsStr].append(seq.id)
          else:
            DuplicatesDict[SeqAsStr] = [SeqToIDsDict[SeqAsStr], seq.id]
        else:
          SeqToIDsDict[SeqAsStr] = seq.id

      # If we found some duplicated sequences, check that no duplication is 
      # within the same bam file / same alias, then print the duplicate names.
      if DuplicatesDict != {}:
        for SeqNames in DuplicatesDict.values():
          aliases = []
          for SeqName in SeqNames:
            RegexMatch = SampleRegex.search(SeqName)
            if RegexMatch and SeqName[:RegexMatch.start()] in BamAliases:
              alias = SeqName[:RegexMatch.start()]
              aliases.append(alias)
          DuplicatedAliases = [alias for alias, count in \
          collections.Counter(aliases).items() if count > 1]
          if DuplicatedAliases != []:
            print('Malfunction of phylotypes - the each of the following bam',
            'files has more than one copy of the same sequence after '
            'processing:', ' '.join(DuplicatedAliases) + '. Quitting.',
            file=sys.stderr)
            exit(1)
      FileForDuplicateReadCountsProcessed = \
      FileForDuplicateReadCountsProcessed_basename + ThisWindowSuffix + '.csv'
      with open(FileForDuplicateReadCountsProcessed, 'w') as f:
        f.write('\n'.join(','.join(SeqNames) for SeqNames in \
        DuplicatesDict.values()))
      OutputFilesByDestinationDir['DupData'].append(
      FileForDuplicateReadCountsProcessed)

  # Update on time taken if desired
  if args.time:
    times.append(time.time())
    LastStepTime = times[-1] - times[-2]
    if args.check_recombination:
      print('All read processing except the recombination calculation in',
      'window', ThisWindowAsStr, 'finished. Number of seconds taken: ',
      LastStepTime)
    else:
      print('All read processing in window', ThisWindowAsStr,
      'finished. Number of seconds taken: ', LastStepTime)

  # Find the read that looks most like a recombinant for each patient.
  if args.check_recombination:
    SamplesToAlnPosDict = {}
    for i, seq in enumerate(SeqAlignmentHere):
      RegexMatch = SampleRegex.search(seq.id)
      if RegexMatch and seq.id[:RegexMatch.start()] in BamAliases:
        SampleName = seq.id[:RegexMatch.start()]
        if SampleName in SamplesToAlnPosDict:
          SamplesToAlnPosDict[SampleName].append(i)
        else:
          SamplesToAlnPosDict[SampleName] = [i]
    RecombinationResults = []
    for alias, ListOfReadPosInAln in SamplesToAlnPosDict.items():
      ThisAliasAln = Align.MultipleSeqAlignment(SeqAlignmentHere[i] for i in \
      ListOfReadPosInAln)
      #(metric, ParentSeq1, ParentSeq2, RecombinantSeq) = \
      result = (alias, ) + pf.CalculateRecombinationMetric(ThisAliasAln,
      RecombNormToDiv, IncludeGaps=args.recombination_gap_aware)
      RecombinationResults.append(result)
    FileForRecombinantReads = FileForRecombinantReads_basename + \
    ThisWindowSuffix + '.csv'
    with open(FileForRecombinantReads, 'w') as f:
      f.write('Bam file,Recombination metric,Parent seq 1,Parent seq 2,' + \
      'Recombinant seq')
      for result in sorted(RecombinationResults, key=lambda x: x[1],
      reverse=True):
        f.write('\n' + ','.join(map(str, result)) )
      # Add to the recombination data file those bams with no reads here:
      for alias in BamAliases:
        if not alias in SamplesToAlnPosDict:
          f.write('\n' + alias + ',NA,NA,NA,NA')
    OutputFilesByDestinationDir['RecombFiles'].append(FileForRecombinantReads)

    # Update on time taken if desired
    if args.time:
      times.append(time.time())
      LastStepTime = times[-1] - times[-2]
      print('Recombination calculation in window', ThisWindowAsStr,
      'finished. Number of seconds taken: ', LastStepTime)

  if args.no_trees:
    continue

  # Check that there are at least 4 seqs before calling RAxML.
  if len(SeqAlignmentHere) < 4:
    print('Warning: the file of aligned reads in this window,', FileForTrees,
    ', contains only ', len(SeqAlignmentHere), ' sequences; at least 4 are ',
    'needed to make a tree. Skipping to the next window.', sep='',
    file=sys.stderr)
    continue

  # Create the ML tree
  if PrintInfo:
    print('Running RAxML on the processed & aligned reads in window',
    ThisWindowAsStr)
    
  NumMLtreesMade += pf.RunRAxML(FileForTrees, RAxMLargList, ThisWindowSuffix,
  ThisWindowAsStr, UserLeftWindowEdge, UserRightWindowEdge, TempFiles,
  TempFileForAllBootstrappedTrees_basename, args.bootstrap_seed,
  args.num_bootstraps, times)

if ExploreWindowWidths:
  TableHeaders = 'Window start,' + ','.join(sorted(BamAliases))
  # Yes, this is clumsy nesting, but it works:
  # Make a dict indexed by width, of dicts indexed by window, of dicts indexed
  # by bam, with the value being read count.
  ReorganisedData = {}
  for WindowStart, WindowEnd, BamAlias, NumReads in WindowWidthExplorationData:
    width = WindowEnd - WindowStart + 1
    if width in ReorganisedData:
      if WindowStart in ReorganisedData[width]:
        ReorganisedData[width][WindowStart][BamAlias] = NumReads
      else: 
        ReorganisedData[width][WindowStart] = {BamAlias : NumReads}
    else:
      ReorganisedData[width] = {WindowStart : {BamAlias : NumReads}}
  OutputTables = ''
  FirstWidth = True
  for width, DataDictOuter in sorted(ReorganisedData.items(),
  key=lambda x: x[0]):
    if not FirstWidth:
      OutputTables += '\n\n'
    else:
      FirstWidth = False
    OutputTables += 'Number of unique reads per-bam and per-window with ' +\
    'window width = ' + str(width) + ':\n' + TableHeaders
    for WindowStart, DataDictInner in sorted(DataDictOuter.items(),
    key=lambda x: x[0]):
      ReadCountsSortedByBam = [count for bam, count in \
      sorted(DataDictInner.items(), key=lambda x: x[0])]
      OutputTables += '\n' + str(WindowStart) + ',' + \
      ','.join(map(str,ReadCountsSortedByBam))
  with open(args.explore_window_width_file, 'w') as f:  
    f.write(OutputTables)
  exit(0)
    

# Make a bam file of discarded read pairs for each input bam file.
if args.inspect_disagreeing_overlaps:
  for BamFileBasename, DiscardedReadPairs in DiscardedReadPairsDict.items():
    if DiscardedReadPairs != []:
      WhichBamFile = BamFileBasenames.index(BamFileBasename)
      RefFile = RefFiles[WhichBamFile]
      if len(BamFileBasename) >= 4 and (BamFileBasename[-4:] == '.bam' or
        BamFileBasename[-4:] == '.BAM' or BamFileBasename[-4:] == '.Bam'):
        BamNameForRef = BamFileBasename[:-4]
      else:
        BamNameForRef = BamFileBasename
      LocalRefFileName = FileForDiscardedReadPairs_basename + \
      BamNameForRef + '_ref.fasta'
      # Copy the relevant reference file to the working directory, so that it's
      # together with the discarded reads file. This might fail e.g. if the same
      # file exists already - then do nothing.
      try:
        shutil.copy2(RefFile, LocalRefFileName)
      except:
        pass
      else:
        OutputFilesByDestinationDir['DiscardedReads'].append(LocalRefFileName)
      OutFile = FileForDiscardedReadPairs_basename +BamFileBasename
      DiscardedReadPairsOut = pysam.AlignmentFile(OutFile, "wb", template=BamFile)
      for read in DiscardedReadPairs:
        DiscardedReadPairsOut.write(read)
      DiscardedReadPairsOut.close()
      OutputFilesByDestinationDir['DiscardedReads'].append(OutFile)


OutputFilesByDestinationDir['raxml'] = glob.glob('RAxML_*')

# Try to create different directories for each kind of output file we've made.
# Move files if desired.
if not args.keep_output_together:
  for DirKey, files in OutputFilesByDestinationDir.items():
    if len(files) > 0:
      if HaveMadeOutputDir:
        OutputDirs[DirKey] = os.path.join(args.output_dir, OutputDirs[DirKey])
      Dir = OutputDirs[DirKey]
      if not os.path.isdir(Dir):
        try:
          os.mkdir(Dir)
        except:
          print('Problem creating the directory', Dir + \
          ". phyloscanner will create all output files in",
          "the working directory, instead of in file-type-specific",
          "subdirectories.", file=sys.stderr)
          args.keep_output_together = True
          break
  for DirKey, files in OutputFilesByDestinationDir.items():
    Dir = OutputDirs[DirKey]
    for File in files:
      shutil.move(File, os.path.join(Dir, File))
if HaveMadeOutputDir:
  if os.path.isfile(FileForAlignedRefs):
    shutil.move(FileForAlignedRefs, os.path.join(args.output_dir,
    FileForAlignedRefs))
  if args.keep_output_together:
    for File in itertools.chain.from_iterable(
    OutputFilesByDestinationDir.values()):
      shutil.move(File, os.path.join(args.output_dir, File))

# Delete temporary files we've made
if not args.keep_temp_files:
  for TempFile in TempFiles:
    try:
      os.remove(TempFile)
    except:
      print('Failed to delete temporary file', TempFile + '. Leaving it.',
      file=sys.stderr)

# We're only printing info henceforth
if not PrintInfo:
  exit(0)

# Stop if no trees
if NumMLtreesMade == 0:
  print("Info: phyloscanner_make_trees.py has processed all windows but has",
  "not produced any trees, either because you told it not to, or because of a",
  "lack of reads, or because of non-fatal errors. Check earlier",
  "warning/error messages.")
  exit(0)

print("\nphyloscanner_make_trees.py successfully produced some trees! If",
"you're happy with them, step two of phyloscanner is analysing the trees. This",
"is done with phyloscanner_analyse_trees.R; run it with --help to learn more.",
"The analysis requires trees to be rooted. You should either (a) tell",
"phyloscanner_analyse_trees.R (using its --outgroupName option) which sequence",
"in each tree to use as an outgroup for rooting, which should be one of the",
"references you included with the bam files here via the",
"--alignment-of-other-refs option; or (b) manually root the trees before",
"giving them as input to phyloscanner_analyse_trees.R.\n")

def FindFilesForRcode(DescriptionOfFiles, FileBasename, DirKey, ROption,
ExtraText=None):
  '''TODO'''
  if args.keep_output_together:
    if HaveMadeOutputDir:
      Dir = args.output_dir
    else:
      Dir = os.getcwd()
  else:
    Dir = OutputDirs[DirKey]
  Dir = os.path.abspath(Dir)
  FileStart = os.path.join(Dir, FileBasename)
  files = glob.glob(FileStart + '*')
  if not files:
    print("Oops, internally we've lost track of the location of the",
    DescriptionOfFiles, "files we produced. "
    "Expected to find files matching", FileStart + '*\nSorry about that.',
    file=sys.stderr)
  else:
    if ExtraText != None:
      end = ExtraText + '\n'
    else:
      end = '\n'
    print("The", DescriptionOfFiles, "files we've produced can be given to",
    'phyloscanner_analyse_trees.R, via its "' + ROption + '" option, as',
    FileStart, end=end)

  
FindFilesForRcode("tree", 'RAxML_bestTree.', 'raxml', 'tree')
if CheckDuplicates:
  FindFilesForRcode("between-bam duplication data", 
  FileForDuplicateReadCountsProcessed_basename, 'DupData',
  "--duplicateBlacklist", " (be sure to also specify a raw and/or relative "
  "threshold for blacklisting duplicates).")
if args.check_recombination:
  FindFilesForRcode("recombination data", FileForRecombinantReads_basename,
  'RecombFiles', '--recombinationFiles')

