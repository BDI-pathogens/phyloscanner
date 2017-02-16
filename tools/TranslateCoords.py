#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251 
##
## Overview: this script translates coordinates with respect to one sequence to
## coordinates with respect to all other sequences in an alignment.
## Usage: call it from the command line thus, for reference-based coordinates:
## ./TranslateCoords.py MyAlignmentFile MyChosenReference coord1 [coord2...]
## or thus, for alignment-based coordinates:
## ./TranslateCoords.py MyAlignmentFile -A coord1 [coord2...]
## The translated coordinates are reported in the order in which they were
## specified.
##
################################################################################
## USER INPUT
# Characters that indicates a gap (missing base)
GapChars = '-.?'
################################################################################

# Import what's needed
import sys, os.path, collections
from optparse import OptionParser
from Bio import SeqIO

# Define the arguments and options
parser = OptionParser()
parser.add_option("-A", action="store_true", dest="AlignmentCoords",
default=False,
help="specify that the coordinates are with respect to the alignment")
(options, args) = parser.parse_args()

# Check this file is called from the command line with the correct number of
# arguments, and that the specified file(s) exist.
if options.AlignmentCoords:
  if len(args) < 2:
    print('At least two arguments are required with the -A option: the',\
    'alignment file and at least one (integer) coordinate.\nQuitting.', \
    file=sys.stderr)
    exit(1)
  AlignmentFile = args[0]
  coords        = args[1:]
else:
  if len(args) < 3:
    print('At least three arguments are required: firstly the alignment file,',\
    'secondly the chosen reference therein, then at least one (integer)',\
    'coordinate for that reference.\nQuitting.', file=sys.stderr)
    exit(1)
  AlignmentFile = args[0]
  ChosenRef     = args[1]
  coords        = args[2:]
if not os.path.isfile(AlignmentFile):
  print(AlignmentFile, 'does not exist or is not a file. Quitting.', \
  file=sys.stderr)
  exit(1)

# Try to understand the coordinates as integers. Check they're positive. Sort.
for i in range(0,len(coords)):
  try:
    coords[i] = int(coords[i])
  except ValueError:
    print('Unable to understand coordinate', coords[i], 'as an integer.'+\
    '\nQuitting.', file=sys.stderr)
    exit(1)
if any(coord < 1 for coord in coords):
  print('All coordinates must be greater than zero. Quitting.', file=sys.stderr)
  exit(1)

# Read in the sequences from the alignment file (into an ordered dictionary)
SeqDict = collections.OrderedDict()
for seq in SeqIO.parse(open(AlignmentFile),'fasta'):
  if seq.id in SeqDict:
    print('Two (or more) sequences in', AlignmentFile, 'are called', seq.id+\
    '. Sequence names should be unique. Quitting.', file=sys.stderr)
    exit(1)
  SeqDict[seq.id] = str(seq.seq)

if len(SeqDict) == 0:
  print("There are no sequences in", AlignmentFile+". Quitting.", \
  file=sys.stderr)
  exit(1)

# Check all sequences have the same length
AlignmentLength = len(SeqDict.values()[0])
if any(len(seq) != AlignmentLength for seq in SeqDict.values()):
  print("The sequences in", AlignmentFile, "are not all of the same length - ",\
  "it's supposed to be an alignment file. Quitting.", file=sys.stderr)

# Check all sequences have at least one base
for SeqName, seq in SeqDict.items():
  NoBases = True
  for base in seq:
    if not base in GapChars:
      NoBases = False
      break
  if NoBases:
    print(SeqName, "has no bases, it's just one big gap.\nQuitting.", \
    file=sys.stderr)
    exit(1)

# If coordinates were specified with respect to the alignment:
if options.AlignmentCoords:

  # Check for coordinates after the end of the alignment
  TooLargeCoords = [coord for coord in coords if coord > AlignmentLength]
  if TooLargeCoords != []:
    print('Coordinates', ', '.join(map(str,TooLargeCoords)), 'occur after the',\
    'end of the alignment ('+str(AlignmentLength), 'bases long).\nQuitting.', \
    file=sys.stderr)
    exit(1)
  CoordsInAlignment_ZeroBased = [coord-1 for coord in coords]

# Coordinates were specified with respect to a chosen reference:
else:

  # Check that the reference is in the alignment
  if not ChosenRef in SeqDict:
    print('Could not find', ChosenRef, 'in', AlignmentFile+'.\nQuitting.', \
    file=sys.stderr)
    exit(1)
  ChosenRefSeq = SeqDict[ChosenRef]

  # Find the coordinates.
  PositionInRef=0
  CoordsInAlignment_ZeroBased = [-1 for coord in coords]
  for PositionMin1,base in enumerate(ChosenRefSeq):
    if not base in GapChars:
      PositionInRef += 1
      for i,coord in enumerate(coords):
        if coord == PositionInRef:
          CoordsInAlignment_ZeroBased[i] = PositionMin1
      if not -1 in CoordsInAlignment_ZeroBased:
        break

  # Check that all coordinates were found
  MissingCoords = \
  [coords[i] for i,coord in enumerate(CoordsInAlignment_ZeroBased) if coord == -1]
  if len(MissingCoords) != 0:
    print('Coordinates', ', '.join(map(str,MissingCoords)), 'occur after the',\
    'end of', ChosenRef, '('+str(PositionInRef), 'bases long).\nQuitting.', \
    file=sys.stderr)
    exit(1)


# Translate those coordinates to the other sequences. Put -1 for a coordinate
# occuring before the start and NaN for a coordinate occuring after the end.
# Put a half-integer (an integer + 0.5) if a coordinate occurs inside a gap.
CoordsDict = collections.OrderedDict()
for SeqName, seq in SeqDict.items():

  # Find the start and the end of this sequence
  StartOfSeq = 0
  while seq[StartOfSeq] in GapChars:
    StartOfSeq += 1
  EndOfSeq = AlignmentLength-1
  while seq[EndOfSeq] in GapChars:
    EndOfSeq -= 1

  # Find the coordinates in this sequence.
  CoordsInThisSeq = [0 for coord in coords]
  PositionInThisSeq = 0
  for PositionMin1,base in enumerate(seq):
    BaseHere = not base in GapChars
    if BaseHere:
      PositionInThisSeq += 1
    for i,coord in enumerate(CoordsInAlignment_ZeroBased):
      if coord == PositionMin1:
        if coord < StartOfSeq:
          CoordsInThisSeq[i] = -1
        elif coord > EndOfSeq:
          CoordsInThisSeq[i] = 'NaN'
        elif BaseHere:
          CoordsInThisSeq[i] = PositionInThisSeq
        else:
          CoordsInThisSeq[i] = PositionInThisSeq + 0.5
    if not 0 in CoordsInThisSeq:
      break
  if 0 in CoordsInThisSeq:
    MissingCoords = [coords[i] for i in range(len(coords))\
    if CoordsInThisSeq[i] == 0]
    print('Internal malfunction of the code: the user-specified coordinates',\
    ' '.join(map(str, MissingCoords)), "were not found in sequence", SeqName+\
    '\nQuitting.', file=sys.stderr)
    exit(1)

  # If coordinates were specified with respect to a chosen reference, check that
  # this process for the chosen ref should recover the user's input coordinates.
  if not options.AlignmentCoords and SeqName == ChosenRef and \
  CoordsInThisSeq != coords:
    print('Internal malfunction of the code: converting the chosen',\
    "reference's coordinates to the alignment coordinates and back again",\
    'gives a different result.\nQuitting.', file=sys.stderr)
    exit(1)

  CoordsDict[SeqName] = CoordsInThisSeq

# Print the output
for SeqName,CoordsInThisSeq in CoordsDict.items():
  print(SeqName, ' '.join(map(str,CoordsInThisSeq)))
