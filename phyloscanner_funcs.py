from __future__ import print_function
import os
import sys
import subprocess
import itertools

GapChar = '-'

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


def ReadNamesFromFile(TheFile, IsFile=True):
  '''Reads names from each line of a file into a list.

  Names are checked to be unique, as they will be used as labels, and may not
  contain whitespace as they may be used as sequence IDs in fasta files.
  If IsFile=True, each name is checked to be a file that exists, and a second
  list - containing the files' basenames - is also returned.
  No names being present causes a non-zero exit.'''

  assert IsFile == True or IsFile == False

  NamesChecked = []
  files = []
  with open(TheFile, 'r') as f:
    for line in f:
      name = line.strip()
      if name == '':
        continue
      if IsFile:
        if not os.path.isfile(name):
          print(name, 'does not exist or is not a file.', file=sys.stderr)
          exit(1)
        NameToCheck = os.path.basename(name)
      else:
        NameToCheck = name
      if len(NameToCheck.split(None,1)) > 1:
        print('Name', NameToCheck, 'in', TheFile, 'contains whitespace.',\
        'Rename to avoid this and try again. Quitting.',\
        file=sys.stderr)
        exit(1)
      if NameToCheck in NamesChecked:
        print('Encountered name', NameToCheck, 'multiple times in', TheFile +\
        '. names should be unique, as they are used as labels. Quitting.',\
        file=sys.stderr)
        exit(1)
      NamesChecked.append(NameToCheck)
      if IsFile:
        files.append(name)
  if NamesChecked == []:
    print(TheFile, 'contains no names. This is assumed to be an error.\n'+\
    'Quitting.', file=sys.stderr)
    exit(1)
  if IsFile:
    return files, NamesChecked
  else:
    return NamesChecked


def TranslateSeqCoordsToAlnCoords(seq, coords):
  '''Takes a sequence that contains gaps (in general), and a set of coordinates
  specified with a respect to that sequence without gaps. The coordinates are
  translated to their positions in the gappy version of the sequence.
  e.g. called with the arguments "-a--cg-t-" and [1,2,3], we return [2,5,6].
  '''
  assert type(seq) == type('abc')
  assert all(type(coord) == int for coord in coords)
  assert all(coord > 0 for coord in coords)
  TranslatedCoords = [-1 for coord in coords]
  PositionInSeq = 0
  for GappyPostitionMin1,base in enumerate(seq):
    if base != GapChar:
      PositionInSeq += 1
      for i,coord in enumerate(coords):
        if coord == PositionInSeq:
          TranslatedCoords[i] = GappyPostitionMin1+1
      if not -1 in TranslatedCoords:
        break
  assert not -1 in TranslatedCoords
  assert len(TranslatedCoords) == len(coords)
  return TranslatedCoords


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

  def RecoverClippedEnds(self):
    '''Replaces any 'None' mapped positions at read ends to continuous ints.

    Recovers clipped ends by considering any bases at the ends of the read
    that are mapped to 'None' to be mapped instead to 1 more than the base to
    the left (at the right end) or 1 less than the base to the right (at the
    end left). e.g. a 9bp read mapped to positions
    None,None,10,11,13,14,None,None,None
    (i.e. clipped on the left by 2bp, and on the right by 3bp, with a 1bp
    deletion in the middle), on processing with this function, is taken to be
    mapped instead to positions
    8,9,10,11,13,14,15,16,17.
    Note that bases overhanging the left end of the references will be
    considered mapped to negative positions, and bases overhanging the right end
    of the reference will be considered mapped to positions greater than the
    length of the reference.'''

    LeftMostMappedBase = 0
    try:
      while self.positions[LeftMostMappedBase] == None:
        LeftMostMappedBase += 1
    except IndexError:
      # Every base mapped to None.
      pass
    else:
      ReadLength = len(self.positions)
      RightMostMappedBase = ReadLength - 1
      while self.positions[RightMostMappedBase] == None:
        RightMostMappedBase -= 1
      if LeftMostMappedBase > 0:
        #print(self.positions)
        RefPosOfLeftEdge = self.positions[LeftMostMappedBase]
        for i in range(LeftMostMappedBase):
          self.positions[i] = RefPosOfLeftEdge - LeftMostMappedBase + i
      if RightMostMappedBase < ReadLength - 1:
        #print(self.positions)
        RefPosOfRightEdge = self.positions[RightMostMappedBase]
        for i in range(RightMostMappedBase + 1, ReadLength):
          self.positions[i] = RefPosOfRightEdge + i - RightMostMappedBase


  def ProcessRead(self, LeftWindowEdge, RightWindowEdge, MinQualForEnds, \
  MinInternalQual, KeepOverhangs, RecoverClippedEnds):
    '''Returns reads that span a given window.
    Overhangs & low-Q bases are trimmed if desired. None is returned for reads
    that do not span the window. The coordinates of the window edges should be
    zero-based.'''

    if RecoverClippedEnds:
      self.RecoverClippedEnds()

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
    (self.IsLowQual(MinInternalQual)):
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
        print('Unexpected behaviour for read', self.name+', which',\
        'maps to the following positions in the reference:\n'+ \
        ' '.join(map(str,self.positions)) +'\nUnable to determine ',\
        'where the window edges ('+str(LeftWindowEdge+1), 'and', \
        str(RightWindowEdge+1)+') are in this read. Skipping it.', \
        file=sys.stderr)
        return None
      SeqToReturn = \
      SeqToReturn[LeftEdgePositionInRead:RightEdgePositionInRead+1]

    return SeqToReturn


  def MergeReadPairOverWindow(self, other, LeftWindowEdge, RightWindowEdge, \
  MinQualForEnds, MinInternalQual, RecoverClippedEnds):
    '''TODO:
    Returns the value None if the pair do not overlap each other and span the
    window. Returns the value False if the pair overlap but disagree on the
    overlap.'''

    assert self.name == other.name

    if RecoverClippedEnds:
      self.RecoverClippedEnds()
      other.RecoverClippedEnds()

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

    # Find where, relative to the reference, the left edge of each read is,
    # *including unmapped reads*. This is the position of the left-most mapped
    # base, minus the length of any unmapped reads to its left. Find which of
    # the two reads has this position more to the left. 
    SelfStartIncClipping  = SelfLeftEdge  - self.positions.index(SelfLeftEdge)
    OtherStartIncClipping = OtherLeftEdge - other.positions.index(OtherLeftEdge)
    if SelfStartIncClipping < OtherStartIncClipping:
      LeftRead = self
      RightRead = other
    else:
      LeftRead = other
      RightRead = self
    Length_LeftRead  = len(LeftRead.positions)
    Length_RightRead = len(RightRead.positions)

    # Slide the reads along each other until we find a position such that 
    # they agree perfectly on the overlap - both on the bases it contains,
    # and on the positions mapped to in the reference. At least one position in
    # the overlap must be mapped, i.e. they can't all be mapped to 'None'.
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
        AtLeastOnePosMapped = False
        for k in range(min(Length_RightRead, Length_LeftRead -j)):
          if LeftRead.positions[j+k] != None:
            AtLeastOnePosMapped = True
            break
        if AtLeastOnePosMapped:
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

def MergeSimilarStrings(DictOfStringCounts, SimilarityThreshold=1):
  '''Absorbs those strings with lower counts into those with higher counts.

  Takes a dict for which the keys are the strings and the values are the counts
  associated with those strings. Sorts it. We iterate forwards - from the rarest
  string (with the lowest count) to the most common string (with the highest
  count) - with index i. For each i, we iterate backwards from the most common
  string with index j, until j reaches i. While iterating backwards, if a
  sufficiently similar string is encountered, we add the count of the rarer
  string to that of the more common string, remove the rarer string, and skip to
  the next i.
  In the resulting dict, every string differs from every other string by an
  amount greater than the SimilarityThreshold.
  '''

  # Check that the keys of the dict are strings.
  for String in DictOfStringCounts:
    if type(String) != type('foo'):
      print('The function MergeSimilarStrings was called with a dict',\
      "containing a key that's not a string.\nQuitting.", file=sys.stderr)
      exit(1)

  # Support duck-typing: the dict's values need not be numbers, we just need to 
  # be able add them together.
  try:
    TotalStringCount = sum(DictOfStringCounts.values())
  except TypeError:
    print('The function MergeSimilarStrings was called with a dict',\
    "containing values of types that cannot be added together.\nQuitting.",\
    file=sys.stderr)
    exit(1)
      
  # Nothing needs to be done to dicts with fewer than two entries.
  NumberOfUniqueStrings = len(DictOfStringCounts)
  if NumberOfUniqueStrings < 2:
    return DictOfStringCounts

  # Sort the strings by their counts
  SortedDict = sorted(DictOfStringCounts.items(), key=lambda x: x[1])

  # Iterate i forwards through the strings, rarest first.
  MergedDict = {}
  for i,(RareString,RareCount) in enumerate(SortedDict):

    # Iterate j backwards through the strings, from the most common one to the
    # current one.
    MatchesAnotherString = False
    MatchingString = None
    for j in reversed(xrange(i+1,NumberOfUniqueStrings)):

      CommonString, CommonCount = SortedDict[j]

      # Compare the two strings. Initialise the number of differences as the
      # difference in string length, then do a pairwise comparison of characters
      # over the length of the shorter string. As soon as the number of 
      # differences exceeds the threshold, we know we shouldn't be merging them. 
      NumDifferingBases = abs(len(RareString) - len(CommonString))
      if NumDifferingBases > SimilarityThreshold:
        break
      for base1, base2 in itertools.izip(RareString, CommonString):
        if base1 != base2:
          NumDifferingBases += 1
          if NumDifferingBases > SimilarityThreshold:
            break
      # If we're merging these two strings, we don't have to check smaller j,
      # so we break.
      if NumDifferingBases <= SimilarityThreshold:
        MatchingString = CommonString
        MatchingCount  = CommonCount
        break

    # NB either string (rare or matching) will already be in MergedDict if and
    # only if we have already merged an even rarer string into it, meaning we
    # should use its count from MergedDict and not from the original dict.

    # If, after looping through j, we found a string to match this i:
    if MatchingString != None:
      if RareString in MergedDict:
        RareCount = MergedDict[RareString]
        del MergedDict[RareString]
      if MatchingString in MergedDict:
        MergedDict[MatchingString] += RareCount
      else:
        MergedDict[MatchingString]  = RareCount + MatchingCount

    # If we did not find a string to match this i: put it in the updated dict,
    # if it's not there already.
    elif not RareString in MergedDict:
      MergedDict[RareString] = RareCount

  # Check the total count is the same before and after.
  TotalStringCountAfterMerging = sum(MergedDict.values())
  if TotalStringCountAfterMerging != TotalStringCount:
    print('FATAL ERROR: the number of strings after merging, '+\
    str(TotalStringCountAfterMerging)+', != the number before merging, '+\
    str(TotalStringCount)+'. Please report this to Chris Wymant',
    'c.wymant@imperial.ac.uk.\nQuitting.', file=sys.stderr)
    exit(1)

  return MergedDict


# for testing:
import random
def GenerateRandomSequence(length, bases='ACGT'):
  '''Generates a random string of the specified size using the specified
  characters (ACGT by default).'''
  return ''.join(random.choice(bases) for _ in range(length))

