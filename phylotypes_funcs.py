from __future__ import print_function
import os
import sys
import subprocess

# Test that we can run code we'll need
DirectoryOfThisScript = os.path.dirname(os.path.realpath(__file__))
def FindAndCheckCode(CodeBasename):
  '''Checks that code exists in the same directory as this script, that it's
  executable with a -h flag, and returns its path.'''
  CodeFullPath = os.path.join(DirectoryOfThisScript, CodeBasename)
  if not os.path.isfile(CodeFullPath):
    print(CodeBasename, 'is not in the same directory as', sys.argv[0] +\
    '\nQuitting', file=sys.stderr)
    exit(1)
  FNULL = open(os.devnull, 'w')
  try:
    ExitStatus = subprocess.call([CodeFullPath, '-h'], stdout=FNULL, \
    stderr=subprocess.STDOUT)
    assert ExitStatus == 0
  except:
    print('Problem running', CodeFullPath+'.\nTry running\nchmod u+x ', \
    CodeFullPath+'\nIt might help...', file=sys.stderr)
    exit(1)
  return CodeFullPath


def ReadFilenamesFromFile(TheFile):
  '''Reads the filenames in a file, and also their basenames, into two lists.

  The files are checked to exist and to have unique basenames without whitespace
  in them. No filenames being present causes a non-zero exit.'''

  files = []
  basenames = []
  with open(TheFile, 'r') as f:
    for line in f:
      filename = line.strip()
      if filename == '':
        continue
      if not os.path.isfile(filename):
        print(filename, 'does not exist or is not a file.', file=sys.stderr)
        exit(1)
      basename = os.path.basename(filename)
      if len(basename.split(None,1)) > 1:
        print('File', filename, 'in', TheFile, 'contains whitespace in the',\
        'basename. Rename to avoid this and try again. Quitting.',\
        file=sys.stderr)
        exit(1)
      if basename in basenames:
        print('Multiple files in', TheFile, 'have basename', basename+\
        '. Basenames should be unique, as they are used as labels. Quitting.',\
        file=sys.stderr)
        exit(1)
      files.append(filename)
      basenames.append(basename)
  if files == []:
    print(TheFile, 'contains no filenames. This is assumed to be an error.\n'+\
    'Quitting.', file=sys.stderr)
    exit(1)
  return files, basenames



class PseudoRead:
  "A class similar to pysam.AlignedSegment. Writable, and with extra features."

  def __init__(self, name, sequence, positions, qualities):
    "A manual constructor. At least one position must not be None."
    assert len(sequence) == len(positions) == len(qualities) > 0
    assert any(value != None for value in positions)
    self.name = name
    self.sequence = sequence
    self.positions = positions
    self.qualities = qualities

  @classmethod
  def InitFromRead(cls, read):
    '''A constructor for pysam.AlignedSegments.
    Not a true constructor, but a decorated class method: the pythonic 
    work-around for multiple constructors.'''
    positions = read.get_reference_positions(full_length=True)
    if not len(read.query_sequence) == len(positions) == \
    len(read.query_qualities) > 0:
      print('Unexpected attribute properties for pysam.AlignedSegment\n', read,\
      '\nQuitting',  file=sys.stderr)
      exit(1)
    return cls(read.query_name, read.query_sequence, positions, \
    read.query_qualities)

  def __repr__(self):
    'Defining how a PseudoRead can be printed'
    return 'name: %s\nseq: %s\npositions: %s\nqualities: %s' % (self.name, \
    self.sequence, ' '.join(map(str,self.positions)), \
    ' '.join(map(str,self.qualities)))

  def SpansWindow(self, LeftWindowEdge, RightWindowEdge):
    "Returns True if the read fully spans the specified window."
    assert LeftWindowEdge <= RightWindowEdge
    LeftMostMappedBase = 0
    try:
      while self.positions[LeftMostMappedBase] == None:
        LeftMostMappedBase += 1
    except IndexError:
      return False
    if self.positions[LeftMostMappedBase] > LeftWindowEdge:
      return False
    RightMostMappedBase = len(self.positions)-1
    while self.positions[RightMostMappedBase] == None:
      RightMostMappedBase -= 1
    return self.positions[RightMostMappedBase] >= RightWindowEdge

  def QualityTrimEnds(self, MinQual):
    '''Trims inwards until a base of sufficient quality is met.
    Returns a blank read if no bases are of sufficient quality.'''
    FirstHighQBase = 0
    try:
      while self.qualities[FirstHighQBase] < MinQual:
        FirstHighQBase += 1
    except IndexError:
      self.sequence = ''
      self.positions = []
      self.qualities = []
    else:
      LastHighQBase = len(self.positions)-1
      while self.qualities[LastHighQBase] < MinQual:
        LastHighQBase -= 1
      self.sequence  = self.sequence[FirstHighQBase:LastHighQBase+1]
      self.positions = self.positions[FirstHighQBase:LastHighQBase+1]
      self.qualities = self.qualities[FirstHighQBase:LastHighQBase+1]

  def IsLowQual(self, MinQual):
    '''Returns True if two or more bases are below the quality threshold, False
    otherwise.'''
    NumLowQbases = 0
    for qual in self.qualities:
      if qual < MinQual:
        NumLowQbases += 1
        if NumLowQbases == 2:
          return True
    return False

  def ProcessRead(self, LeftWindowEdge, RightWindowEdge, MinQualForEnds, \
  MinInternalQual, KeepOverhangs):
    '''Returns reads that span a given window.
    Overhangs & low-Q bases are trimmed if desired. None is returned for reads
    that do not span the window. The coordinates of the window edges should be
    zero-based.'''

    # Skip reads that only partially overlap the window
    if not self.SpansWindow(LeftWindowEdge, RightWindowEdge):
      return None

    # Trim low-Q ends if desired. Skip if that leaves only a partial overlap.
    if MinQualForEnds != None:
      self.QualityTrimEnds(MinQualForEnds)
      if not self.SpansWindow(LeftWindowEdge, RightWindowEdge):
        return None

    # Skip reads containing more than one low-quality base
    if MinInternalQual != None and \
    (self.IsLowQual(MinInternalQual):
      return None

    SeqToReturn = self.sequence

    # Trim the part of the read overhanging the window if desired.
    if not KeepOverhangs:
      try:
        LeftEdgePositionInRead = 0
        while self.positions[LeftEdgePositionInRead] == None or \
        self.positions[LeftEdgePositionInRead] < LeftWindowEdge:
          LeftEdgePositionInRead += 1
        RightEdgePositionInRead = len(self.positions)-1
        while self.positions[RightEdgePositionInRead] == None or \
        self.positions[RightEdgePositionInRead] > RightWindowEdge:
          RightEdgePositionInRead -= 1
        assert LeftEdgePositionInRead <= RightEdgePositionInRead
      except (IndexError, AssertionError):
        print('Unexpected behaviour for read', read.name+', which',\
        'maps to the following positions in the reference:\n'+ \
        ' '.join(map(str,self.positions)) +'\nUnable to determine ',\
        'where the window edges ('+str(LeftWindowEdge+1), 'and', \
        str(RightWindowEdge+1)+') are in this read. Quitting.', file=sys.stderr)
        exit(1)
      SeqToReturn = \
      SeqToReturn[LeftEdgePositionInRead:RightEdgePositionInRead+1]

    return SeqToReturn


  def MergeReadPairOverWindow(self, other, LeftWindowEdge, RightWindowEdge, \
  MinQualForEnds, MinInternalQual):
    '''TODO:
    Returns the value None if the pair do not overlap each other and span the
    window. Returns the value False if the pair overlap but disagree on the
    overlap.'''

    assert self.name == other.name

    # Trim low-Q ends if desired. If either read has no sequence left after 
    # trimming, return None.
    if MinQualForEnds != None:
      self.QualityTrimEnds(MinQualForEnds)
      other.QualityTrimEnds(MinQualForEnds)
      if self.sequence == '' or other.sequence == '':
        return None

    # Check that the pair overlap and span the window. If not, return None.
    SelfLeftEdge  = next(pos for pos in self.positions if pos != None)
    SelfRightEdge = next(pos for pos in reversed(self.positions) if pos != None)
    OtherLeftEdge  = next(pos for pos in other.positions if pos != None)
    OtherRightEdge = next(pos for pos in reversed(other.positions) \
    if pos != None)
    if OtherRightEdge < SelfLeftEdge or OtherLeftEdge > SelfRightEdge or \
    min(SelfLeftEdge, OtherLeftEdge) > LeftWindowEdge or \
    max(SelfRightEdge, OtherRightEdge) < RightWindowEdge:
      return None

    if SelfLeftEdge < OtherLeftEdge:
      LeftRead = self
      RightRead = other
    else:
      LeftRead = other
      RightRead = self
    Length_LeftRead  = len(LeftRead.positions)
    Length_RightRead = len(RightRead.positions)

    # Slide the reads along each other until we find a position such that 
    # they agree perfectly on the overlap - both on the bases it contains,
    # and on the positions mapped to in the reference. 
    # If no such position is found, they disagree: return False.
    OverlapStartInLeftRead = None
    for j in range(Length_LeftRead):
      this_j_works = True
      for k in range(min(Length_RightRead, Length_LeftRead -j)):
        if LeftRead.positions[j+k] != RightRead.positions[k] or \
        LeftRead.sequence[j+k] != RightRead.sequence[k]:
          this_j_works = False
          break
      if this_j_works:
        OverlapStartInLeftRead = j
        break
    if OverlapStartInLeftRead == None:
      return False

    # Merge the two reads, conservatively taking the quality of each base in the
    # overlap to be the larger from the two reads. (As opposed to a quality
    # corresponding to the probability that both reads independently made the
    # same miscall, for example.)
    merged_sequence = LeftRead.sequence[:OverlapStartInLeftRead] +\
    RightRead.sequence +\
    LeftRead.sequence[OverlapStartInLeftRead+Length_RightRead:]
    merged_positions = LeftRead.positions[:OverlapStartInLeftRead] +\
    RightRead.positions +\
    LeftRead.positions[OverlapStartInLeftRead+Length_RightRead:]
    merged_qualities = LeftRead.sequence[:OverlapStartInLeftRead]
    for j in range(\
    max(Length_LeftRead - OverlapStartInLeftRead, Length_RightRead)):
      try:
        BaseQLeftRead = LeftRead.sequence[OverlapStartInLeftRead+j]
      except IndexError:
        BaseQLeftRead  = float('-inf')
      try:
        BaseQRightRead = RightRead.sequence[j]
      except IndexError:
        BaseQRightRead = float('-inf')
      BestBaseQ = max(BaseQLeftRead, BaseQRightRead)
      merged_qualities += BestBaseQ
    assert len(merged_qualities) == len(merged_sequence)
    MergedRead = PseudoRead(self.name, merged_sequence, \
    merged_positions, merged_qualities)
    return MergedRead


def IsMonoSampleClade(clade, SampleRegex):
  '''Checks whether all tips inside this clade come from the same sample.

  We check that all tip names match the regex, and have the same string
  before the regex, and that string is one of the BamFileBasenames. If so, a 
  list of the tip names is returned; if not, the value False is returned.'''

  TipNames = []
  for tip in clade.get_terminals():
    TipName = tip.name
    RegexMatch = SampleRegex.search(TipName)
    if RegexMatch and TipName[:RegexMatch.start()] in BamFileBasenames:
      SampleName = TipName[:RegexMatch.start()]
    else:
      return False
    if TipNames == []:
      FirstSample = SampleName
    elif SampleName != FirstSample:
      return False
    TipNames.append(TipName)
  return TipNames


# TODO: how to make this work for unrooted trees? find_clades assumes a root, I 
# think.
def ResolveTree(tree):
  '''Resolves a tree into a list, each element of which is either a) a list of
  tip names where these tips form a mono-sample clade, or b) a tip name.'''
  ResolvedCladesAtThisLevel = []
  for clade in tree.find_clades(order='level'):
    CladeLevel = len(tree.get_path(clade))
    if CladeLevel == 0:
      continue
    elif CladeLevel == 1:
      if clade.is_terminal():
        ResolvedCladesAtThisLevel.append(clade.name)
      else:
        TipNamesIfMonoSample = IsMonoSampleClade(clade)
        if TipNamesIfMonoSample:
          ResolvedCladesAtThisLevel.append(TipNamesIfMonoSample)
        else:
          ResolvedCladesAtThisLevel += ResolveTree(clade)
    else:
      break
  return ResolvedCladesAtThisLevel
