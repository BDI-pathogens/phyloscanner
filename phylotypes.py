#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script creates between-sample alignments of reads.
The reference sequences used to create a number of bam files (i.e. the sequences
to which the reads were mapped) are aligned. User-specified coordinates with
respect to this alignment are translated to coordinates with respect to each
separate reference, dividing each reference into a set of matching windows. For
each sample, for each window: all reads mapped to that window are found,
identical reads are collected together with an associated count, similar reads
are merged together based on the counts, then a minimum count is imposed. Then,
for each window, all reads from all samples are aligned using mafft and a
phylogeny is constructed using RAxML.
Output files are written to the current working directory; to avoid overwriting
existing files, you might to want to call this code from an empty directory.
'''

################################################################################
# USER INPUT

RAxMLseed = 1
RAxMLbootstrapSeed = 1

# Some output files
FileForAlignedReads_basename = 'AlignedReads'
FileForAlignedReads_PositionsExcised_basename = 'AlignedReads_PositionsExcised_'
FileForAlignedRefs = 'RefsAln.fasta'
FileForDiscardedReadPairs_basename = 'DiscardedReads_'
FileForSurvivingDuplicates_basename = 'DuplicateReads_surviving_'
FileForEliminatedDuplicates_basename = 'DuplicateReads_eliminated_'

# Some temporary working files we'll create
FileForRefs = 'temp_refs.fasta'
FileForPairwiseUnalignedRefs = 'temp_2Refs.fasta'
FileForPairwiseAlignedRefs = 'temp_2RefsAln.fasta'
FileForReads_basename = 'temp_UnalignedReads'
FileForOtherRefs_basename = 'temp_OtherRefs'
FileForAllBootstrappedTrees_basename = 'temp_AllBootstrappedTrees'
################################################################################
GapChar = '-'

import os
import collections
import subprocess
import sys
import re
import copy
import glob
from shutil import copy2
import argparse
import pysam
from Bio import SeqIO
from Bio import Seq
from Bio import Phylo
from Bio import AlignIO
from MergeSimilarStrings import MergeSimilarStrings
import phylotypes_funcs as pf
#from Bio import AlignIO
#from matplotlib import pyplot as plt
#from Bio.Phylo.Consenss import bootstrap_trees, majority_consensus, get_support
#from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, \
#DistanceCalculator

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Define a comma-separated integers object, as a type for the argparse.
def CommaSeparatedInts(MyCommaSeparatedInts):
  try:
    ListOfInts = [int(value) for value in MyCommaSeparatedInts.split(',')]
  except:
    raise argparse.ArgumentTypeError('Unable to understand the window '+\
    'coordinates as comma-separated integers.')
  else:
    return ListOfInts

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('MergingThreshold', type=int, help=\
'Reads that differ by a number of bases equal to or less than this are merged'+\
', following the algorithm in MergeSimilarStrings.py. A value equal to or '+\
'less than 0 turns off merging.')
parser.add_argument('MinReadCount', type=int, help=\
'Reads with a count less than this value (after merging) are discarded. A '+\
'value equal to or less than 1 means all reads are kept.')
parser.add_argument('ListOfBamFiles', type=File, help='A file containing the '+\
'names (and paths) of the bam files to be included, one per line. The file '+\
'basenames (i.e. the filename minus the directory) should be unique and free '+\
'of whitespace.')
parser.add_argument('ListOfRefFiles', type=File, help='A file containing the '+\
'names (and paths) of the reference fasta files for the bam files, one per '+\
'line. The file basenames (i.e. the filename minus the directory) should be'+\
' unique and free of whitespace.')
parser.add_argument('-W', '--window-coords', type=CommaSeparatedInts, \
help='A comma-separated series of paired coordinates defining the boundaries '+\
'of the windows. e.g. 1,300,11,310,21,320 would define windows 1-300, 11-310,'+\
' 21-320.')
parser.add_argument('-AW', '--auto-window-params', \
type=CommaSeparatedInts, help='Two or three comma-separated integers '+\
'controlling automatic window finding: the first integer is the length you '+\
'want windows to be, when weighting each column by its non-gap fraction; the '+\
'second is the overlap between the end of one window and the start of the next'\
+'(which can be negative); the optional third integer is the start position '+\
'for the first window (by default, 1).')
parser.add_argument('-A', '--alignment-of-other-refs', type=File,\
help='An alignment of any reference sequences (which need not be those used '+\
'to produce the bam files) to be cut into the same windows as the bam files '+\
'and included in the alignment of reads (e.g. to help root trees).')
parser.add_argument('-C', '--ref-for-coords', help='The coordinates are '\
'specified with respect to the reference named after this flag. (By default '\
+'coordinates are with respect to the alignment of all references.)')
parser.add_argument('-D', '--dont-check-duplicates', action='store_true', \
help="Don't compare reads between samples to find duplicates - a possible "+\
"indication of contamination. (By default this check is done.)")
parser.add_argument('-F', '--renaming-file', type=File, help='Specify a file '\
'with one line per bam file, showing how reads from that bam file should be '\
'named in the output files.')
parser.add_argument('-I', '--discard-improper-pairs', action='store_true', \
help='Any improperly paired reads will be discarded')
parser.add_argument('-M', '--raxml-model', default='GTRCAT',\
help='The evoltionary model used by RAxML (GTRCAT by default).')
parser.add_argument('-N', '--number-of-bootstraps', type=int,\
help='The number of bootstraps to be calculated for RAxML trees (by default, '+\
'none i.e. only the ML tree is calculated).')
parser.add_argument('-O', '--keep-overhangs', action='store_true', \
help='Keep the whole read. (By default, only the part of the read inside the'+\
'window is kept, i.e. overhangs are trimmed.)')
parser.add_argument('-P', '--merge-paired-reads', action='store_true', \
help='Merge overlapping paired reads into a single read.')
parser.add_argument('-Q1', '--quality-trim-ends', type=int, help='Each end of'+\
' the read is trimmed inwards until a base of this quality is met.')
parser.add_argument('-Q2', '--min-internal-quality', type=int, help=\
'Discard reads containing more than one base of a quality below this parameter'\
+'. If used in conjuction with the --quality-trim-ends option, the trimming '+\
'of the ends is done first.')
parser.add_argument('-R', '--ref-for-rooting', help='Used to name a reference'\
'(which must be present in the file you specify with -A) to be an outgroup'+\
' in the tree.')
#parser.add_argument('-S', '--min-support', default=60, type=float, help=\
#'The bootstrap support below which nodes will be collapsed, as a percentage.')
parser.add_argument('--align-refs-only', action='store_true', help='Align the'+\
'references in the bam files (plus any extras specified with -A) then quit '+\
'without parsing the reads.')
parser.add_argument('-T', '--no-trees', action='store_true', help='Generate '+\
'aligned sets of reads for each window then quit without making trees.')
parser.add_argument('-XC', '--excision-coords', type=CommaSeparatedInts, \
help='Used to specify a comma-separated set of integer coordinates that will '+\
'be excised from the aligned reads. Useful for sites of non-neutral '+\
'evolution. Requires the -XR flag.')
parser.add_argument('-XR', '--excision-ref', help='Used to specify the '+\
'name of a reference (which must be present in the file you specify with -A)'+\
'with respect to which the coordinates specified with -XC are interpreted.')
parser.add_argument('-2', '--pairwise-align-to', help='Sequentially, and '+\
'separately, align the bam file references to this reference (which must be '+\
'present in the file you specify with -A) instead of aligning all references '+\
'together. Window coordinates will be interpreted with respect to this '+\
'reference.')
parser.add_argument('--x-raxml', default='raxmlHPC-AVX', help=\
'The command required to invoke RAxML (by default: raxmlHPC-AVX).')
parser.add_argument('--x-mafft', default='mafft', help=\
'The command required to invoke mafft (by default: mafft).')
parser.add_argument('--x-samtools', default='samtools', help=\
'The command required to invoke samtools, if needed (by default: samtools).')

args = parser.parse_args()

# Shorthand
WindowCoords  = args.window_coords
NumBootstraps = args.number_of_bootstraps
UserSpecifiedCoords = args.window_coords != None
AutoWindows = args.auto_window_params != None
IncludeOtherRefs = args.alignment_of_other_refs != None
QualTrimEnds  = args.quality_trim_ends != None
ImposeMinQual = args.min_internal_quality != None
ExcisePositions = args.excision_coords != None
PairwiseAlign = args.pairwise_align_to != None


# Check that window coords have been specified either manually or automatically.
if (UserSpecifiedCoords and AutoWindows) or \
(not UserSpecifiedCoords and not AutoWindows):
  print('Exactly one of the --window-coords or --auto-window-params', \
  'options should specified. (Not neither, not both.) Quitting.', \
  file=sys.stderr)
  exit(1)

# The user shouldn't specify a reference for their coords to be interpreted with
# to, and choose automatic windows (i.e. not specify any coordinates).
if AutoWindows and args.ref_for_coords != None:
  print('The --ref-for-coords and --auto-window-params', \
  'options should not be specified together: the first means your coordinates',\
  'should be interpreted with respect to a named reference, and the second',\
  "means you're not specfiying any coordinates. Quitting.", file=sys.stderr)
  exit(1)

# The -XR and -XC flags should be used together or not all.
if (ExcisePositions and args.excision_ref == None) or \
((not ExcisePositions) and args.excision_ref != None):
  print('The --excision-coords and --excision-ref options require each other:',\
  'use both, or neither. Quitting.', file=sys.stderr)
  exit(1)

# TODO: remove this testing
#read1 = pf.PseudoRead('read1', 'abcdefghij', [1,2,3,4,5,6,7,8,9,10], [30]*10)

# Record the names of any external refs being included.
# If we're doing pairwise reference alignments, we'll need the seqs too. We'll
# also need gappy and gapless copies of the ref chosen for pairwise alignment.
if IncludeOtherRefs:
  try:
    ExternalRefAlignment = AlignIO.read(args.alignment_of_other_refs, "fasta")
  except:
    print('Problem reading', args.alignment_of_other_refs + ':', \
    file=sys.stderr)
    raise
  ExternalRefNames = []
  for ref in ExternalRefAlignment:
    ExternalRefNames.append(ref.id)
    if ref.id == args.pairwise_align_to:
      RefForPairwiseAlnsGappySeq = str(ref.seq)
      RefForPairwiseAlns = copy.deepcopy(ref)
      RefForPairwiseAlns.seq = RefForPairwiseAlns.seq.ungap("-")


# Consistency checks on flags that require a ref.
for FlagName, FlagValue in (('--ref-for-coords',  args.ref_for_coords),\
('--ref-for-rooting', args.ref_for_rooting), \
('--pairwise-align-to', args.pairwise_align_to), \
('--excision-ref', args.excision_ref)):
  if FlagValue == None:
    continue
  if not IncludeOtherRefs:
    print('The', FlagName, 'flag requires the --alignment-of-other-refs',\
    'flag. Quitting.', file=sys.stderr)
    exit(1)
  if not FlagValue in ExternalRefNames:
    print('Reference', FlagValue +', specified with the', FlagName, \
    'flag, was not found in', args.alignment_of_other_refs +'. Quitting.', \
    file=sys.stderr)
    exit(1)

# Sanity checks on using the pairwise alignment option.
if PairwiseAlign:
  if args.ref_for_coords != None:
    print('Note that if the --pairwise-align-to option is used, using the',\
    '--ref-for-coords as well is redundant.', file=sys.stderr)
    if args.ref_for_coords != args.pairwise_align_to:
      print('Furthermore you have chosen two different values for these flags,'\
      , 'indicating some confusion as to their use. Try again.')
      exit(1)
  if AutoWindows:
    print('As you have chosen that references are aligned in a pairwise',\
    'manner, please specify coordinates manually - the automatic option is',\
    "for stepping through a global alignment of all references. Quitting.",
    file=sys.stderr)
    exit(1)
  if ExcisePositions and args.excision_ref != args.pairwise_align_to:
    print('The --pairwise-align-to and --excision-ref options can only be',\
    'used at once if the same reference is specified for both. Qutting.',\
    file=sys.stderr)
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
else:
  # Sanity checks on auto window parameters
  if not len(args.auto_window_params) in [2,3]:
    print('The --auto-window-params option requires either 2 or 3 integers.',\
    'Quitting.', file=sys.stderr)
    exit(1)
  WeightedWindowWidth = args.auto_window_params[0]
  WindowOverlap       = args.auto_window_params[1]
  if len(args.auto_window_params) == 3:
    WindowStartPos    = args.auto_window_params[2]
    if WindowStartPos < 1:
      print('The start position for the --auto-window-params option must be',\
      'greater than zero. Quitting.', file=sys.stderr)
      exit(1)
  else:
    WindowStartPos = 1
  if WeightedWindowWidth <= 0:
    print('The weighted window width for the --auto-window-params option must',\
    'be greater than zero. Quitting.', file=sys.stderr)
    exit(1)

# Sort excision coords, largest to smallest
if ExcisePositions:
  args.excision_coords = sorted(args.excision_coords, reverse=True)

# Check that the bootstrap threshold is between 0 and 100
#if not (0 <= args.min_support <= 100):
#  print('MIN_SUPPORT was given as', str(args.min_support)+'; it should be',
#  'between 0 and 100 inclusive.\nQuitting.', file=sys.stderr)

TranslateCoordsCode = pf.FindAndCheckCode('TranslateCoords.py')
FindSeqsInFastaCode = pf.FindAndCheckCode('FindSeqsInFasta.py')
FindWindowsCode     = pf.FindAndCheckCode('FindInformativeWindowsInFasta.py')

# Test RAxML works, if trees are to be made.
if not args.no_trees:
  FNULL = open(os.devnull, 'w')
  try:
    ExitStatus = subprocess.call([args.x_raxml, '-h'], stdout=FNULL, \
    stderr=subprocess.STDOUT)
    assert ExitStatus == 0
  except:
    print('Problem running', args.x_raxml, '\nQuitting.', file=sys.stderr)
    raise

# Read in lists of bam and reference files
BamFiles, BamFileBasenames = pf.ReadNamesFromFile(args.ListOfBamFiles)
RefFiles, RefFileBasenames = pf.ReadNamesFromFile(args.ListOfRefFiles)
if args.renaming_file != None:
  BamAliases = pf.ReadNamesFromFile(args.renaming_file, False)
else:
  BamAliases = BamFileBasenames

# If the BamFileBasenames are all still unique after removing ".bam" from the
# ends, do so, for aesthetics in output files.
BamlessBasenames = []
for BamFileBasename in BamFileBasenames:
  if len(BamFileBasename) >= 4 and BamFileBasename[-4:] == '.bam':
    BamlessBasenames.append(BamFileBasename[:-4])
  else:
    BamlessBasenames.append(BamFileBasename)
if len(BamlessBasenames) == len(set(BamlessBasenames)):
  BamFileBasenames = BamlessBasenames

# Check that there are the same number of bam and reference files
NumberOfBams = len(BamFiles)
if NumberOfBams != len(RefFiles):
  print('Different numbers of files are listed in', args.ListOfBamFiles, 'and',\
  args.ListOfRefFiles+'.\nQuitting.', file=sys.stderr)
  exit(1)
if args.renaming_file != None and len(BamAliases) != NumberOfBams:
  print('Different numbers of files are listed in', args.ListOfBamFiles, 'and',\
  args.renaming_file+'.\nQuitting.', file=sys.stderr)
  exit(1)

# Read in all the reference sequences. Set each name by the file from which the
# ref comes.
RefSeqs = []
for i,RefFile in enumerate(RefFiles):
  SeqList = list(SeqIO.parse(open(RefFile),'fasta'))
  if len(SeqList) != 1:
    print('There are', len(SeqList), 'sequences in', RefFile+'. There should',\
    'be exactly 1.\nQuitting.', file=sys.stderr)
    exit(1)
  SeqList[0].id = RefFileBasenames[i]
  RefSeqs += SeqList


def TranslateCoords(CodeArgs):
  '''Runs TranslateCoordsCode with the supplied args, and returns the results as
  a dict.'''

  try:
    CoordsString = subprocess.check_output([TranslateCoordsCode]+CodeArgs)
  except:
    print('Problem executing', TranslateCoordsCode +'. Quitting.', \
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
            print('Unable to understand the coordinate', coords[i], \
            'as an integer in line\n' +line +'\nin the output of '+\
            TranslateCoordsCode+'\nQuitting.', file=sys.stderr)
            exit(1)
    CoordsDict[SeqName] = coords
  return CoordsDict


# If there is only one bam and no other refs, no coordinate translation
# is necessary - we use the coords as they are, though setting any after the end
# of the reference to be equal to the end of the reference.
UserCoords = WindowCoords
if NumberOfBams == 1 and not IncludeOtherRefs:
  if AutoWindows:
    print('As you are supplying a single bam file and no external references,',\
    'please specify coordinates manually - the automatic option is designed',\
    "for an alignment of multiple sequences. Quitting.", file=sys.stderr)
    exit(1)
  if args.align_refs_only:
    print('As you are supplying a single bam file and no external references,',\
    "the --align-refs-only option makes no sense - there's nothing to align.",\
    "Quitting.", file=sys.stderr)
    exit(1)
  RefSeqLength = len(RefSeqs[0])
  CoordsToUse = [min(coord, RefSeqLength) for coord in WindowCoords]
  CoordsInRefs = {RefFileBasenames[0] : CoordsToUse}

# If there are at least two bam files, or if there is one but we're including
# other refs, we'll be aligning references and translating the user-specified
# coords with respect to each sequence, then storing those coords in a dict
# indexed by the ref's name.
else:

  # If we're separately and sequentially pairwise aligning our references to
  # a chosen ref in order to determine window coordinates, do so now.
  if PairwiseAlign:

    # Find the coordinates with respect to the chosen ref, in the alignment of
    # just the external refs - we'll need these later.
    ExternalRefWindowCoords = \
    pf.TranslateSeqCoordsToAlnCoords(RefForPairwiseAlnsGappySeq, WindowCoords)

    CoordsInRefs = {}
    for BamRefSeq in RefSeqs:

      # Align
      SeqIO.write([RefForPairwiseAlns,BamRefSeq], FileForPairwiseUnalignedRefs,\
      "fasta")
      with open(FileForPairwiseAlignedRefs, 'w') as f:
        try:
          ExitStatus = subprocess.call([args.x_mafft, '--quiet',  \
          '--preservecase', FileForPairwiseUnalignedRefs], stdout=f)
          assert ExitStatus == 0
        except:
          print('Problem calling mafft. Quitting.', file=sys.stderr)
          raise

      # Translate.
      # The index names in the PairwiseCoordsDict, labelling the coords found by
      # coord translation, should coincide with the two seqs we're considering.
      PairwiseCoordsDict = TranslateCoords([FileForPairwiseAlignedRefs, \
      args.pairwise_align_to] + [str(coord) for coord in WindowCoords])
      if set(PairwiseCoordsDict.keys()) != \
      set([BamRefSeq.id,args.pairwise_align_to]):
        print('Malfunction of phylotypes: mismatch between the sequences',\
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
    SeqIO.write(RefSeqs, FileForRefs, "fasta")
    if IncludeOtherRefs:
      FinalMafftOptions = ['--add', FileForRefs, args.alignment_of_other_refs]
    else:
      FinalMafftOptions = [FileForRefs]
    with open(FileForAlignedRefs, 'w') as f:
      try:
        ExitStatus = subprocess.call([args.x_mafft, '--quiet',  \
        '--preservecase'] + FinalMafftOptions, stdout=f)
        assert ExitStatus == 0
      except:
        print('Problem calling mafft. Quitting.', file=sys.stderr)
        raise

    if args.align_refs_only:
      print('References aligned in', FileForAlignedRefs+ \
      '. Quitting successfully.')
      exit(0)

    # If window coords were specified with respect to one particular reference, 
    # or if we are excising certain coords, translate to alignment coords.
    if args.ref_for_coords != None or ExcisePositions:
      for seq in SeqIO.parse(open(FileForAlignedRefs),'fasta'):
        if seq.id == args.ref_for_coords:
          WindowCoords = \
          pf.TranslateSeqCoordsToAlnCoords(str(seq.seq), WindowCoords)
        if seq.id == args.excision_ref:
          RefForExcisionGappySeq = str(seq.seq)
          AlignmentExcisionCoords = pf.TranslateSeqCoordsToAlnCoords(\
          RefForExcisionGappySeq, args.excision_coords)


    # Determine windows automatically if desired
    if AutoWindows:
      try:
        WindowsString = subprocess.check_output([FindWindowsCode, \
        FileForAlignedRefs, str(WeightedWindowWidth), str(WindowOverlap), \
        '-S', str(WindowStartPos)])
      except:
        print('Problem executing', FindWindowsCode +'. Quitting.', \
        file=sys.stderr)
        raise
      try:
        WindowCoords = [int(value) for value in WindowsString.split(',')]
      except:
        print('Unable to understand the', FindWindowsCode, 'output -', \
        WindowsString, '- as comma-separated integers. Quitting.', \
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
    if set(CoordsInRefs.keys()) != set(RefFileBasenames+ExternalRefNames):
      print('Malfunction of phylotypes: mismatch between the sequences found', \
      'in the output of', TranslateCoordsCode, 'and those in', \
      FileForAlignedRefs +'. Quitting.', file=sys.stderr)
      exit(1)

# Make index files for the bam files if needed.
for BamFileName in BamFiles:
  if not os.path.isfile(BamFileName+'.bai'):
    try:
      ExitStatus = subprocess.call([args.x_samtools, 'index', BamFileName])
      assert ExitStatus == 0
    except:
      print('Problem running samtools index.\nQuitting.', file=sys.stderr)
      raise

# Gather some data from each bam file
BamFileRefSeqNames = {}
BamFileRefLengths  = {}
for i,BamFileName in enumerate(BamFiles):

  BamFileBasename = BamFileBasenames[i]
  RefBasename = RefFileBasenames[i]

  # Prep for pysam
  BamFile = pysam.AlignmentFile(BamFileName, "rb")

  # Find the reference in the bam file; there should only be one.
  AllReferences = BamFile.references
  if len(AllReferences) != 1:
    print('Expected exactly one reference in', BamFileName+'; found',\
    str(len(AllReferences))+'.\nQuitting.', file=sys.stderr)
    exit(1)
  BamFileRefSeqNames[BamFileBasename] = AllReferences[0]

  # Get the length of the reference.
  AllReferenceLengths = BamFile.lengths
  if len(AllReferenceLengths) != 1:
    print('Pysam error: found one reference but', len(AllReferenceLengths), \
    'reference lengths.\nQuitting.', file=sys.stderr)
    exit(1)
  RefLength = AllReferenceLengths[0]
  BamFileRefLengths[BamFileBasename] = RefLength

  # When translating coordinates, -1 means before the sequence starts; 'NaN'
  # means after it ends. These should be replaced by 1 and the reference length
  # respectively.
  for j,coord in enumerate(CoordsInRefs[RefBasename]):
    if coord == -1:
      CoordsInRefs[RefBasename][j] = 1
    elif coord == 'NaN':
      CoordsInRefs[RefBasename][j] = RefLength

def ProcessReadDict(ReadDict, WhichBam, LeftWindowEdge, RightWindowEdge):
  '''Turns a dict of reads into a list of reads, merging & imposing a minimum
  count.'''

  # For naming things
  BamFileBasename = BamFileBasenames[WhichBam]
  BasenameForReads = BamAliases[WhichBam]

  # Merge similar reads if desired
  if args.MergingThreshold > 0:
    ReadDict = MergeSimilarStrings(ReadDict, args.MergingThreshold)

  # Implement the minimum read count
  if args.MinReadCount > 1:
    ReadDict = {read:count for read, count in ReadDict.items() if \
    count >= args.MinReadCount}

  # Warn if there are no reads
  if len(ReadDict) == 0:
    print('Warning: bam file', BamFileBasename, 'has no reads in window', \
    str(LeftWindowEdge+1)+'-'+   str(RightWindowEdge+1), file=sys.stderr)
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

def ReadFastaOfReadsIntoDicts(FastaFile, ValuesAreCounts=True):
  '''Collects sample seqs and into dicts, and other seqs into a list.

  The values of the dicts are either the seq count (inferred from the seq name)
  or simply the seq name.'''
  SampleReadCounts = collections.OrderedDict()
  NonSampleSeqs = []
  for seq in SeqIO.parse(open(FastaFile),'fasta'):
    RegexMatch = SampleRegex.search(seq.id)
    if RegexMatch and seq.id[:RegexMatch.start()] in BamAliases:
      SampleName = seq.id[:RegexMatch.start()]
      read = str(seq.seq)
      if ValuesAreCounts:
        value = int(seq.id.rsplit('_',1)[1])
      else:
        value = seq.id
      if SampleName in SampleReadCounts:
        if read in SampleReadCounts[SampleName]:
          print('Malfunction of phylotypes:', FastaFile, \
          'contains two identical sequences for sample', SampleName+\
          '. This should not happen. Quitting.', file=sys.stderr)
          exit(1)
        SampleReadCounts[SampleName][read] = value
      else:
        SampleReadCounts[SampleName] = {read : value}
    else:
      if not seq.id in ExternalRefNames:
        print('Malfunction of phylotypes: sequence', seq.id, 'in', \
        FileForReads, 'not recognised as a read nor as an external reference.',\
        'Quitting.', file=sys.stderr)
        exit(1)
      NonSampleSeqs.append(seq)
  return SampleReadCounts, NonSampleSeqs

# We'll keep a list of discarded read pairs for each bam file:
DiscardedReadPairsDict = \
{BamFileBasename:[] for BamFileBasename in BamFileBasenames}

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

  print('Now processing window ', UserLeftWindowEdge, '-', UserRightWindowEdge,\
  sep='')

  # Get ready to record reads here from all samples
  AllReadsInThisWindow = []
  if not args.dont_check_duplicates:
    AllReadDictsInThisWindow = []

  # Iterate through the bam files
  for i,BamFileName in enumerate(BamFiles):

    # Recall some things we've already worked out for this bam file and stored.
    BamFileBasename = BamFileBasenames[i]
    RefSeqName = BamFileRefSeqNames[BamFileBasename]
    RefLength = BamFileRefLengths[BamFileBasename]
    RefBasename = RefFileBasenames[i]
    ThisBamCoords = CoordsInRefs[RefBasename]
    LeftWindowEdge  = ThisBamCoords[window*2]
    RightWindowEdge = ThisBamCoords[window*2 +1]

    # Find all unique reads in this window and count their occurrences.
    # NB pysam uses zero-based coordinates for positions w.r.t the reference
    AllReads = {}
    UniqueReads = {}
    LeftWindowEdge  = LeftWindowEdge  -1
    RightWindowEdge = RightWindowEdge -1
    BamFile = pysam.AlignmentFile(BamFileName, "rb")
    for read in BamFile.fetch(RefSeqName, LeftWindowEdge, RightWindowEdge):

      # Skip improperly paired reads if desired
      if args.discard_improper_pairs and read.is_paired and \
      not read.is_proper_pair:
        continue

      ReadAsPseudoRead = pf.PseudoRead.InitFromRead(read)

      if args.merge_paired_reads:

        # We've seen this read's mate already. Merge the pair.
        if read.query_name in AllReads:
          Read1 = AllReads[read.query_name]
          Read1asPseudoRead = pf.PseudoRead.InitFromRead(Read1)
          Read2 = read
          Read2asPseudoRead = ReadAsPseudoRead
          MergedRead = Read1asPseudoRead.MergeReadPairOverWindow( \
          Read2asPseudoRead, LeftWindowEdge, RightWindowEdge, \
          args.quality_trim_ends, args.min_internal_quality)
          if MergedRead == None:
            del AllReads[read.query_name]
            continue
          elif MergedRead == False:
            DiscardedReadPairsDict[BamFileBasename] += [Read1,Read2]
            del AllReads[read.query_name]
            continue
          AllReads[read.query_name] = MergedRead

        # We've not come across a read with this name before. Record & move on.
        else:
          AllReads[read.query_name] = read

      # If we're not merging reads, process this read now to save memory.
      # ProcessRead returns None if we don't want to consider this read.
      else:
        seq = ReadAsPseudoRead.ProcessRead(LeftWindowEdge, RightWindowEdge, \
          args.quality_trim_ends, args.min_internal_quality, \
          args.keep_overhangs)
        if seq == None:
          continue
        if seq in UniqueReads:
          UniqueReads[seq] += 1
        else:
          UniqueReads[seq] = 1

    # If we did merge paired reads, we now need to process them.
    # AllReads will be a mixture of PseudoRead instances (for merged read pairs)
    # and pysam.AlignedSegment instances (for unmerged single reads). The latter
    # must be converted to PseudoRead instances to be processed.
    if args.merge_paired_reads:
      for read in AllReads.values():
        try:
          seq = read.ProcessRead(LeftWindowEdge, RightWindowEdge, \
          args.quality_trim_ends, args.min_internal_quality, \
          args.keep_overhangs)
        except AttributeError:
          #print(type(read))
          ReadAsPseudoRead = pf.PseudoRead.InitFromRead(read)          
          seq = ReadAsPseudoRead.ProcessRead(LeftWindowEdge, RightWindowEdge, \
          args.quality_trim_ends, args.min_internal_quality, \
          args.keep_overhangs)
        if seq == None:
          continue
        if seq in UniqueReads:
          UniqueReads[seq] += 1
        else:
          UniqueReads[seq] = 1

    # If we are checking for read duplication between samples, record the file 
    # name and read dict for this sample and move on to the next sample.
    if not args.dont_check_duplicates:
      AllReadDictsInThisWindow.append((BamFileBasename, UniqueReads, \
      LeftWindowEdge, RightWindowEdge))

    # If we're not checking for read duplication between samples, process the
    # read dict for this sample now and add it to the list of all reads here.
    else:
      AllReadsInThisWindow += \
      ProcessReadDict(UniqueReads, i, LeftWindowEdge, RightWindowEdge)

  # We've now gathered together reads from all bam files for this window.

  # If we're checking for duplicate reads between samples, do so now.
  # Check every dict against every other dict, and record the ratio of counts
  # for any shared reads...
  if not args.dont_check_duplicates:
    DuplicateDetails = []
    for i, (BamFile1Basename, ReadDict1, LeftWindowEdge1, RightWindowEdge1) \
    in enumerate(AllReadDictsInThisWindow):
      for j, (BamFile2Basename, ReadDict2, LeftWindowEdge2, RightWindowEdge2) \
      in enumerate(AllReadDictsInThisWindow[i+1:]):
        DuplicateReadRatios = []
        for read in ReadDict1:
          if read in ReadDict2:
            BamFile1Alias = BamAliases[BamFileBasenames.index(BamFile1Basename)]
            BamFile2Alias = BamAliases[BamFileBasenames.index(BamFile2Basename)]
            Bam1Count = ReadDict1[read]
            Bam2Count = ReadDict2[read]
            DuplicateDetails.append(\
            (BamFile1Alias, BamFile2Alias, read, Bam1Count, Bam2Count))
    # ... and process the read dicts.
    for i, (BamFileBasename, ReadDict, LeftWindowEdge, RightWindowEdge) \
    in enumerate(AllReadDictsInThisWindow):
      AllReadsInThisWindow += \
      ProcessReadDict(ReadDict, i, LeftWindowEdge, RightWindowEdge)

  # All read dicts have now been processed into the list AllReadsInThisWindow.

  # Skip empty windows.
  if AllReadsInThisWindow == []:
    print('WARNING: no bam file had any reads (after a minimum post-merging '+\
    'read count of', args.MinReadCount, 'was imposed) in the window', \
    str(UserLeftWindowEdge)+'-'+str(UserRightWindowEdge)+'. Skipping to the', \
    'next window.', file=sys.stderr)
    continue

  # Re-define the window edge coords to be with respect to the alignment of refs
  # rather than a bam file.
  LeftWindowEdge  = WindowCoords[window*2]
  RightWindowEdge = WindowCoords[window*2 +1]

  # Create a fasta file with all reads in this window, ready for aligning.
  # If there's only one, we don't need to align (or make trees!).
  FileForReadsHere = FileForReads_basename + ThisWindowSuffix+\
  '.fasta'
  FileForAlnReadsHere = FileForAlignedReads_basename + \
  ThisWindowSuffix +'.fasta'
  if len(AllReadsInThisWindow) == 1 and not IncludeOtherRefs:
    SeqIO.write(AllReadsInThisWindow, FileForAlnReadsHere, "fasta")
    print('There is only one read in this window, written to ' +\
    FileForAlnReadsHere +'. Skipping to the next window.')
    continue
  SeqIO.write(AllReadsInThisWindow, FileForReadsHere, "fasta")

  # If external refs are included, find the part of each one's seq corresponding
  # to this window and put them all in another file.
  # If we did pairwise aligning of refs, we know the coordinates we want in the
  # ExternalRefAlignment object. If we did a global alignment, we slice the
  # desired window out of that alignment.
  if IncludeOtherRefs:
    FileForOtherRefsHere = FileForOtherRefs_basename + \
    ThisWindowSuffix +'.fasta'
    if PairwiseAlign:
      ExternalRefLeftWindowEdge  = ExternalRefWindowCoords[window*2]
      ExternalRefRightWindowEdge = ExternalRefWindowCoords[window*2 +1]
      RefAlignmentInWindow = ExternalRefAlignment[:, \
      ExternalRefLeftWindowEdge-1:ExternalRefRightWindowEdge]
      AlignIO.write(RefAlignmentInWindow, FileForOtherRefsHere, 'fasta')
    else:
      with open(FileForOtherRefsHere, 'w') as f:
        try:
          ExitStatus = subprocess.call([FindSeqsInFastaCode, \
          FileForAlignedRefs, '-B', '-W', str(LeftWindowEdge) + ',' + \
          str(RightWindowEdge), '-v'] + RefFileBasenames, stdout=f)
          assert ExitStatus == 0
        except:
          print('Problem calling', FindSeqsInFastaCode+\
          '. Skipping to the next window.', file=sys.stderr)
          continue

  # Align the reads. Prepend 'temp_' to the file name if we'll merge again after
  # aligning.
  if args.MergingThreshold > 0:
    FileForReads = 'temp_'+FileForAlnReadsHere
  else:
    FileForReads = FileForAlnReadsHere
  if IncludeOtherRefs:
    FinalMafftOptions = ['--add', FileForReadsHere, FileForOtherRefsHere]
  else:
    FinalMafftOptions = [FileForReadsHere]
  with open(FileForReads, 'w') as f:
    try:
      ExitStatus = subprocess.call([args.x_mafft, '--quiet', '--preservecase']+\
      FinalMafftOptions, stdout=f)
      assert ExitStatus == 0
    except:
      print('Problem calling mafft. Skipping to the next window.', \
      file=sys.stderr)
      continue

  # Do a second round of within-sample read merging now the reads are aligned. 
  # Make a dict (indexed by sample name) of dicts (indexed by the sequences
  # themselves) of read counts. Those sequences that are from a sample are 
  # found by matching the RegexMatch '_read_\d+_count_\d+$'; other sequences
  # should be external references the user included, and are not processed.
  if args.MergingThreshold > 0:
    SampleReadCounts, RefSeqs = ReadFastaOfReadsIntoDicts(FileForReads)
    SampleSeqsToPrint = []
    for SampleName in SampleReadCounts:
      SampleReadCounts[SampleName] = \
      MergeSimilarStrings(SampleReadCounts[SampleName], args.MergingThreshold)
      for k, (read, count) in enumerate(sorted(\
      SampleReadCounts[SampleName].items(), key=lambda x: x[1], reverse=True)):
        ID = SampleName+'_read_'+str(k+1)+'_count_'+str(count)
        SeqObject = SeqIO.SeqRecord(Seq.Seq(read), id=ID, description='')
        SampleSeqsToPrint.append(SeqObject)
    AllSeqsToPrint = SampleSeqsToPrint + RefSeqs
    if AllSeqsToPrint == []:
      print('Malfunction of phylotypes: no sequences were found in', \
      FileForReads +'. Quitting.', file=sys.stderr)
      exit(1)

    # Merging after alignment means some columns could be pure gap.
    # Remove these.
    PureGapColumns = []
    FirstSeq = str(AllSeqsToPrint[0].seq)
    for position,base in enumerate(FirstSeq):
      if base == '-':
        PureGapColumns.append(position)
    if PureGapColumns != []:
      for seq in AllSeqsToPrint[1:]:
        SeqAsString = str(seq.seq)
        GapColsToRemove = []
        for i,col in enumerate(PureGapColumns):
          if SeqAsString[col] != '-':
            GapColsToRemove.append(col)
        PureGapColumns = \
        [col for col in PureGapColumns if not col in GapColsToRemove]
        if PureGapColumns == []:
          break
      if PureGapColumns != []:
        for i,seq in enumerate(AllSeqsToPrint):
          SeqAsString = str(seq.seq)
          for position in PureGapColumns[::-1]:
            SeqAsString = SeqAsString[:position]+SeqAsString[position+1:]
          AllSeqsToPrint[i].seq = Seq.Seq(SeqAsString)
    SeqIO.write(AllSeqsToPrint, FileForAlnReadsHere, "fasta")
    FileForTrees = FileForAlnReadsHere

  # If we're checking for duplicates, read in the file of aligned reads (which
  # may have gone through a second round of merging), remove gaps, and see which
  # of the duplicates we found before merging (and imposing a minimum count) are
  # are still there.
  if not args.dont_check_duplicates and DuplicateDetails != []:
    SampleReadNames, RefSeqs = \
    ReadFastaOfReadsIntoDicts(FileForAlnReadsHere, False)
    UngappedReadNames = {}
    for SampleName,ReadDict in SampleReadNames.items():
      UngappedReadNames[SampleName] = \
      {read.replace('-',''):SeqName for read,SeqName in ReadDict.items()}
    SurvivingDuplicates  = []
    EliminatedDuplicates = []
    for BamFile1Alias, BamFile2Alias, read, Bam1Count, Bam2Count in \
    DuplicateDetails:
      for WhichBam,alias in enumerate([BamFile1Alias, BamFile2Alias]):
        if WhichBam == 0:
          ThisCount, TheOtherCount = Bam1Count, Bam2Count
          ThisAlias, TheOtherAlias = BamFile1Alias, BamFile2Alias
        else:
          ThisCount, TheOtherCount = Bam2Count, Bam1Count
          ThisAlias, TheOtherAlias = BamFile2Alias, BamFile1Alias
        CountRatio = str(float(ThisCount)/TheOtherCount)
        # The read may no longer exist for this alias; also this alias may have
        # no reads present in the file at all (if all merged reads failed the
        # minimum count). Either will cause a KeyError, meaning this is an
        # eliminated duplicate.
        try:
          DuplicateSeqName = UngappedReadNames[alias][read]
        except KeyError:
          EliminatedDuplicates.append((ThisAlias, TheOtherAlias, CountRatio))
        else:
          SurvivingDuplicates.append((DuplicateSeqName, TheOtherAlias, \
          CountRatio))
    FileForSurvivingDuplicates = FileForSurvivingDuplicates_basename + \
    ThisWindowSuffix + '.csv'
    with open(FileForSurvivingDuplicates, 'w') as f:
      f.write('"SeqName","SharedWith","SeqCount/SharedCount"\n')
      f.write('\n'.join(','.join(data) for data in SurvivingDuplicates))
    FileForEliminatedDuplicates = FileForEliminatedDuplicates_basename + \
    ThisWindowSuffix + '.csv'
    with open(FileForEliminatedDuplicates, 'w') as f:
      f.write('"BamFile","SharedWith","SeqCount/SharedCount"\n')
      f.write('\n'.join(','.join(data) for data in EliminatedDuplicates))

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
    print('CoordsToExciseInThisWindow', CoordsToExciseInThisWindow)
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
          PositionsInUngappedRef.append(\
          len(RefInThisWindowGappy[:DistanceIntoWindow+1].replace('-','')))
        UngappedRefHere = RefInThisWindowGappy.replace('-','')

      # Read in the aligned reads, and check the ref looks as expected.
      RefInAlignment = None
      SeqAlignmentHere = AlignIO.read(FileForAlnReadsHere, "fasta")
      for seq in SeqAlignmentHere:
        if seq.id == args.excision_ref:
          RefInAlignment = str(seq.seq)
          break
      if RefInAlignment == None:
        print('Malfunction of phylotypes: unable to find', args.excision_ref, \
        'in', FileForAlnReadsHere +'. Quitting.', file=sys.stderr)
        exit(1)
      if RefInAlignment.replace('-','') != UngappedRefHere:
        print('Malfunction of phylotypes: mismatch between the ref for',\
        'excision we expected to find in this window:\n', UngappedRefHere,\
        '\nand the ref for excision we actually found in this window:\n',\
        RefInAlignment.replace('-',''), '\nQuitting.', file=sys.stderr)
        exit(1)

      # Excise the positions in the aligned set of reads.
      PositionsInAlignment = \
      pf.TranslateSeqCoordsToAlnCoords(RefInAlignment, PositionsInUngappedRef)
      assert PositionsInAlignment == sorted(PositionsInAlignment, reverse=True)
      for pos in PositionsInAlignment:
        SeqAlignmentHere = \
        SeqAlignmentHere[:, :pos-1] + SeqAlignmentHere[:, pos:]
      AlignIO.write(SeqAlignmentHere, FileForAlignedReads_PositionsExcised, \
      'fasta')
      FileForTrees = FileForAlignedReads_PositionsExcised

  if args.no_trees:
    continue

  # Create the ML tree
  MLtreeFile = 'RAxML_bestTree.' +ThisWindowSuffix +'.tree'
  RAxMLcall = [args.x_raxml, '-m', args.raxml_model, '-p', str(RAxMLseed),\
  '-s', FileForTrees, '-n', ThisWindowSuffix+'.tree']
  if args.ref_for_rooting != None:
    RAxMLcall += ['-o', args.ref_for_rooting]
  try:
    ExitStatus = subprocess.call(RAxMLcall)
    assert ExitStatus == 0
  except:
    print('Problem making the ML tree with RAxML.\nSkipping to the next', \
    'window.', file=sys.stderr)
    continue
  if not os.path.isfile(MLtreeFile):
    print(MLtreeFile +', expected to be produced by RAxML, does not exist.'+\
    '\nSkipping to the next window.', file=sys.stderr)
    continue

  # If desired, make bootstrapped alignments
  if NumBootstraps != None:
    try:
      ExitStatus = subprocess.call([args.x_raxml, '-m', args.raxml_model, '-p',\
      str(RAxMLseed), '-b', str(RAxMLbootstrapSeed), '-f', 'j', '-#', \
      str(NumBootstraps), '-s', FileForTrees, '-n', ThisWindowSuffix+\
      '_bootstraps'])
      assert ExitStatus == 0
    except:
      print('Problem generating bootstrapped alignments with RAxML', \
      '\nSkipping to the next window.', file=sys.stderr)
      continue
    BootstrappedAlignments = [FileForTrees+'.BS'+str(bootstrap) for \
    bootstrap in range(NumBootstraps)]
    if not all(os.path.isfile(BootstrappedAlignment) \
    for BootstrappedAlignment in BootstrappedAlignments):
      print('At least one of the following files, expected to be produced by'+\
      ' RAxML, is missing:\n', ' '.join(BootstrappedAlignments)+\
      '\nSkipping to the next window.', file=sys.stderr)
      continue

    # Make a tree for each bootstrap
    for bootstrap,BootstrappedAlignment in enumerate(BootstrappedAlignments):
      try:
        ExitStatus = subprocess.call([args.x_raxml, '-m', args.raxml_model, \
        '-p', str(RAxMLseed), '-s', BootstrappedAlignment, '-n', \
        ThisWindowSuffix+'_bootstrap_'+str(bootstrap)+'.tree'])
        assert ExitStatus == 0
      except:
        print('Problem generating a tree with RAxML for bootstrap', \
        str(bootstrap), '\Breaking.', file=sys.stderr)
        break
    BootstrappedTrees = ['RAxML_bestTree.' +ThisWindowSuffix +'_bootstrap_' +\
    str(bootstrap) +'.tree' for bootstrap in range(NumBootstraps)]
    if not all(os.path.isfile(BootstrappedTree) \
    for BootstrappedTree in BootstrappedTrees):
      print('At least one of the following files, expected to be produced by'+\
      ' RAxML, is missing:\n', ' '.join(BootstrappedTrees)+\
      '\nSkipping to the next window.', file=sys.stderr)
      continue

    # Collect the trees from all bootstraps into one file
    AllBootstrappedTreesFile = FileForAllBootstrappedTrees_basename +\
    ThisWindowSuffix+'.tree'
    with open(AllBootstrappedTreesFile, 'w') as outfile:
      for BootstrappedTree in BootstrappedTrees:
        with open(BootstrappedTree, 'r') as infile:
          outfile.write(infile.read())

    # Collect the trees from all bootstraps onto the ML tree
    MainTreeFile = 'MLtreeWbootstraps' +ThisWindowSuffix +'.tree'
    try:
      ExitStatus = subprocess.call([args.x_raxml, '-m', args.raxml_model, '-p',\
      str(RAxMLseed), '-f', 'b', '-t', MLtreeFile, '-z', \
      AllBootstrappedTreesFile, '-n', MainTreeFile])
      assert ExitStatus == 0
    except:
      print('Problem collecting all the bootstrapped trees onto the ML tree', \
      'with RAxML.\nSkipping to the next window.', file=sys.stderr)
      continue
    MainTreeFile = 'RAxML_bipartitions.' +MainTreeFile
    if not os.path.isfile(MainTreeFile):
      print(MainTreeFile +', expected to be produced by RAxML, does not '+\
      'exist.\nSkipping to the next window.', file=sys.stderr)
      continue

  # With no bootstraps, just use the ML tree:
  else:
    MainTreeFile = MLtreeFile

  #MainTree = Phylo.read(MainTreeFile, 'newick')
  #for TipOrMonoSampleClade in ResolveTree(MainTree):
  #  print(TipOrMonoSampleClade)

  #MainTree.collapse_all(lambda c: c.confidence is not None and \
  #c.confidence < args.min_support)
  #for clade in MainTree.find_clades(order='level'):
  #  node_path = MainTree.get_path(clade)
  #  if len(node_path) == 0:     
  #    print('whole tree?')
  #    parent = 'N/A'
  #  elif len(node_path) == 1: 
  #    parent = MainTree.root 
  #  else:
  #    parent = node_path[-2]
  #  if  len(node_path) == 1: 
  #    print(clade.is_terminal(), MainTree.get_path(clade))
  #    print('parent:', parent)
  #    print(' '.join([tip.name for tip in clade.get_terminals()]))
  #  continue
  #  if not clade.is_terminal():
  #    print('Subclade:')
  #    for clade2 in clade.find_clades(order='level'):
  #      print(clade2.is_terminal())
  #      print('MainTree.get_path(clade):', MainTree.get_path(clade2))
  #      print('clade.get_path(clade):', clade.get_path(clade2))
  #      print(' '.join([tip.name for tip in clade2.get_terminals()]))
  #  print()

  #  #if clade.name == None:
  #  #  for clade2 in clade.find_clades():
  #  #print(clade2.name, clade2.confidence, clade2.count_terminals(), \
  #  #clade2.is_preterminal(), '\n', clade2, '\n\n')
  #  if clade2.is_preterminal()

  #MainTree.ladderize()   # Flip branches so deeper clades are displayed at top
  #with open(MainTreeFile+'_image.txt', 'w') as f:
  #  Phylo.draw_ascii(MainTree, file=f, column_width=1000)

  #plt.ion()
  #Phylo.draw(MainTree)
  #plt.savefig('foo.pdf')

# Make a bam file of discarded read pairs for each input bam file.
DiscardedReadPairsFiles = []
for BamFileBasename, DiscardedReadPairs in DiscardedReadPairsDict.items():
  if DiscardedReadPairs != []:
    WhichBamFile = BamFileBasenames.index(BamFileBasename)
    RefFile = RefFiles[WhichBamFile]
    LocalRefFileName = BamFileBasename+'_ref.fasta'
    # Copy the relevant reference file to the working directory, so that it's
    # together with the discarded reads file. This might fail e.g. if the same
    # file exists already - then do nothing.
    try:
      copy2(RefFile, LocalRefFileName)
    except:
      pass
    if len(BamFileBasename) >= 4 and BamFileBasename[-4:] == '.bam':
      OutFile = FileForDiscardedReadPairs_basename +BamFileBasename
    else:
      OutFile = FileForDiscardedReadPairs_basename +BamFileBasename +'.bam'
    DiscardedReadPairsOut = pysam.AlignmentFile(OutFile, "wb", template=BamFile)
    for read in DiscardedReadPairs:
      DiscardedReadPairsOut.write(read)
    DiscardedReadPairsOut.close()
    DiscardedReadPairsFiles.append(OutFile)
if DiscardedReadPairsFiles != []:
  print('Info: read pairs that overlapped but disagreed on the overlap were',\
  'found. These have been written to', ' '.join(DiscardedReadPairsFiles) +'.')


# Some code not being used at the moment:
'''DuplicateReadRatios.append(float(ReadDict1[read])/ReadDict2[read])
        if DuplicateReadRatios != []:
          DuplicateDetails.append([BamFile1Basename, BamFile2Basename] + \
          DuplicateReadRatios)
    if DuplicateDetails != []:
      DuplicateDetails.sort(key=lambda entry: len(entry), reverse=True)
      FileForDuplicates = FileForDuplicates_basename + ThisWindowSuffix + '.csv'
      with open(FileForDuplicates, 'w') as f:
        f.write('"BamFile1","BamFile2","BamFile1Count/BamFile2Count'+\
        ' for each duplicated read"\n')
        f.write('\n'.join(','.join(map(str,data)) for data in DuplicateDetails))'''
