from __future__ import print_function
import os
import sys
import subprocess
import itertools
import csv
import time

GapChar = '-'

# Test that we can run code we'll need
DirectoryOfThisScript = os.path.dirname(os.path.realpath(__file__))
def FindAndCheckCode(CodeBasename):
  '''Checks that code exists in the same directory as this script, that it's
  executable with a -h flag, and returns its path.'''
  CodeFullPath = os.path.join(DirectoryOfThisScript, CodeBasename)
  if not os.path.isfile(CodeFullPath):
    print(CodeBasename, 'is not in the same directory as', __file__ +\
    '\nQuitting', file=sys.stderr)
    exit(1)
  FNULL = open(os.devnull, 'w')
  try:
    ExitStatus = subprocess.call([CodeFullPath, '-h'], stdout=FNULL,
    stderr=subprocess.STDOUT)
    assert ExitStatus == 0
  except:
    print('Problem running', CodeFullPath+'.\nTry running\nchmod u+x ',
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
        print('Name', NameToCheck, 'in', TheFile, 'contains whitespace.',
        'Rename to avoid this and try again. Quitting.',
        file=sys.stderr)
        exit(1)
      if NameToCheck in NamesChecked:
        print('Encountered name', NameToCheck, 'multiple times in', TheFile +\
        '. names should be unique, as they are used as labels. Quitting.',
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


def ReadInputCSVfile(TheFile):
  '''Reads in a csv file listing the bams, refs and optionally aliases.
  
  Bam and ref files are checked to exist; bam file base names (i.e. the file
  name after stripping the path), and aliases if present, are required to be
  unique.'''

  assert os.path.isfile(TheFile), TheFile + \
  ' does not exist or is not a file. Quitting.'

  with open(TheFile, 'r') as f:
    reader = csv.reader(f, delimiter=',', quotechar='"')
    BamFiles = []
    RefFiles = []
    aliases = []
    BamBaseNames = []
    for LineNumberMin1, fields in enumerate(reader):

      # Check for the correct number of fields
      if LineNumberMin1 == 0:
        NumFields = len(fields)
        if NumFields != 2 and NumFields != 3:
          print(TheFile, 'should contain either 2 or 3 comma-separated fields;',
          'found', NumFields, 'on the first line. Quitting.', file=sys.stderr)
          exit(1)
        RenamingColPresent = NumFields == 3
      elif len(fields) != NumFields:
        print('Line', LineNumberMin1 + 1, 'of', TheFile, 'contains',
        len(fields), 'fields, but the first line contains', str(NumFields) +
        '. All lines should have the same number of fields. Quitting.',
        file=sys.stderr)
        exit(1)

      # Check the bam and ref files exist, and don't contain whitespace.
      BamFile = fields[0].strip()
      RefFile = fields[1].strip()
      for FileToCheck in (BamFile, RefFile):
        if not os.path.isfile(FileToCheck):
          print(FileToCheck + ', specified in ' + TheFile + \
          ', does not exist or is not a file. Quitting.', file=sys.stderr)
          exit(1)

      # Check the bam basename, and alias if present, is unique and
      # whitespace-free.
      BamBaseName = os.path.basename(BamFile)
      if len(BamBaseName.split(None, 1)) > 1:
        print(BamBaseName, 'named in', TheFile, 'contains whitespace.',
        'Rename to avoid this and try again. Quitting.',
        file=sys.stderr)
        exit(1)
      if BamBaseName in BamBaseNames:
        print('A bam file with base name', BamBaseName, 'was multiply specified in', TheFile + \
        '. Bam file base names should be unique. Quitting.')
        exit(1)
      BamBaseNames.append(BamBaseName)
      if RenamingColPresent:
        alias = fields[2].strip()
        if len(alias.split(None, 1)) > 1:
          print(alias, 'named in', TheFile, 'contains whitespace.',
          'Rename to avoid this and try again. Quitting.',
          file=sys.stderr)
          exit(1)
        if alias in aliases:
          print('The alias', alias, 'was found a second time on line',
          LineNumberMin1 + 1, 'in', TheFile +
          '. Aliases should be unique. Quitting.')
          exit(1)
        aliases.append(alias)

      BamFiles.append(BamFile)
      RefFiles.append(RefFile)

  if not RenamingColPresent:
    aliases = BamBaseNames

  return BamFiles, RefFiles, aliases, BamBaseNames


def TestRAxML(ArgString, DefaultFlags, HelpMessage):
  '''Runs RAxML with the desired options and --flag-check.'''

  # The user has specified how to call RAxML. Try it.
  if ArgString != None:
    ArgList = ArgString.split()
    out = None
    err = None
    try:
      proc = subprocess.Popen(ArgList + ['--flag-check'],
      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
      success = False
      raise
    else:
      out, err = proc.communicate()
      success = proc.returncode == 0
    finally:
      if not success:
        print('Error: could not run the command "', ArgString,
        ' --flag-check".', sep='', file=sys.stderr)
        if out != None and err != None:
          print('RAxML produced the error messages: ', out + '\n' + err, sep='',
          file=sys.stderr)
        else:
          print('If RAxML is not installed, please install it first. If it is',
          'installed, it seems you need to specify a different executable',
          'and/or set of options. Quitting.', file=sys.stderr)
        print('Quitting.', file=sys.stderr)
        exit(1)

  # The user has not specified how to call RAxML. Try different executables.
  else:
    FlagList = DefaultFlags.split()
    success = False
    out = None
    err = None
    ExesToTry = ['raxmlHPC-AVX', 'raxmlHPC-SSE3', 'raxmlHPC']
    for exe in ExesToTry:
      ArgList = [exe] + FlagList
      try:
        proc = subprocess.Popen(ArgList + ['--flag-check'],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      except:
        continue
      else:
        out, err = proc.communicate()
        if proc.returncode != 0:
          continue
      success = True
      break
    if not success:
      print("Error: could not successfully call RAxML using any of the ",
      "commands ", ', '.join(ExesToTry), ', with the options "', DefaultFlags,
      '".', sep='', file=sys.stderr)
      if out != None and err != None:
        print('RAxML produced the error messages: ', out + '\n' + err, sep='',
        file=sys.stderr)
      else:
        print('If RAxML is not installed, please install it first. If it is ',
        'installed, try adding the path containing its executable files to the',
        ' PATH environment variable of your terminal, or rerunning',
        ' using the --x-raxml option:\n', HelpMessage, sep='',
        file=sys.stderr)
      print('Quitting.', file=sys.stderr)
      exit(1)

  return ArgList

def RunRAxML(alignment, RAxMLargList, WindowSuffix, WindowAsStr, LeftEdge,
RightEdge, TempFilesSet, TempFileForAllBootstrappedTrees_basename,
BootstrapSeed=None, NumBootstraps=None, TimesList=[]):
  '''Runs RAxML on aligned sequences in a window, with bootstraps if desired.

  Returns 1 if an ML tree was produced (regardless of whether any subsequent
  bootstrapping worked), 0 if not.'''

  # Update on times if we weren't given an empty list
  UpdateTimes = TimesList != []

  MLtreeFile = 'RAxML_bestTree.' + WindowSuffix + '.tree'
  RAxMLcall = RAxMLargList + ['-s', alignment, '-n',
  WindowSuffix+'.tree']
  proc = subprocess.Popen(RAxMLcall, stdout=subprocess.PIPE,
  stderr=subprocess.PIPE)
  out, err = proc.communicate()
  ExitStatus = proc.returncode
  if ExitStatus != 0:
    print('Problem making the ML tree with RAxML in window ', WindowAsStr,
    '. It returned an exit code of ', ExitStatus, ', printed this to stdout:\n',
    out, '\nand printed this to stderr:\n', err,
    '\nSkipping to the next window.', sep='', file=sys.stderr)
    return 0
  if not os.path.isfile(MLtreeFile):
    print(MLtreeFile +', expected to be produced by RAxML, does not exist.'+\
    '\nSkipping to the next window.', file=sys.stderr)
    return 0

  # Update on time taken if desired
  if UpdateTimes:
    TimesList.append(time.time())
    LastStepTime = TimesList[-1] - TimesList[-2]
    print('ML tree in window', WindowAsStr,
    'finished. Number of seconds taken: ', LastStepTime)

  # If desired, make bootstrapped alignments
  if NumBootstraps != None:
    try:
      ExitStatus = subprocess.call(RAxMLargList + ['-b',
      str(BootstrapSeed), '-f', 'j', '-#', str(NumBootstraps), '-s',
      alignment, '-n', WindowSuffix + '_bootstraps'])
      assert ExitStatus == 0
    except:
      print('Problem generating bootstrapped alignments with RAxML in window ',
      WindowAsStr, '. Skipping to the next window.', sep='', file=sys.stderr)
      return 1
    BootstrappedAlignments = [alignment+'.BS'+str(bootstrap) for \
    bootstrap in range(NumBootstraps)]
    if not all(os.path.isfile(BootstrappedAlignment) \
    for BootstrappedAlignment in BootstrappedAlignments):
      print('At least one of the following files, expected to be produced by'+\
      ' RAxML, is missing:\n', ' '.join(BootstrappedAlignments)+\
      '\nSkipping to the next window.', file=sys.stderr)
      return 1

    # Make a tree for each bootstrap
    for bootstrap,BootstrappedAlignment in enumerate(BootstrappedAlignments):
      try:
        ExitStatus = subprocess.call(RAxMLargList + ['-s',
        BootstrappedAlignment, '-n', WindowSuffix + '_bootstrap_' + \
        str(bootstrap)+'.tree'])
        assert ExitStatus == 0
      except:
        print('Problem generating a tree with RAxML for bootstrap',
        str(bootstrap), '. Skipping subsequent bootstraps.', file=sys.stderr)
        break
    BootstrappedTrees = ['RAxML_bestTree.' +WindowSuffix +'_bootstrap_' +\
    str(bootstrap) +'.tree' for bootstrap in range(NumBootstraps)]
    if not all(os.path.isfile(BootstrappedTree) \
    for BootstrappedTree in BootstrappedTrees):
      print('At least one of the following files, expected to be produced by'+\
      ' RAxML, is missing:\n', ' '.join(BootstrappedTrees)+\
      '\nSkipping to the next window.', file=sys.stderr)
      return 1

    # Collect the trees from all bootstraps into one file
    TempAllBootstrappedTreesFile = TempFileForAllBootstrappedTrees_basename +\
    WindowSuffix+'.tree'
    with open(TempAllBootstrappedTreesFile, 'w') as outfile:
      for BootstrappedTree in BootstrappedTrees:
        with open(BootstrappedTree, 'r') as infile:
          outfile.write(infile.read())
    TempFilesSet.add(TempAllBootstrappedTreesFile)

    # Collect the trees from all bootstraps onto the ML tree
    MainTreeFile = 'MLtreeWbootstraps' +WindowSuffix +'.tree'
    try:
      ExitStatus = subprocess.call(RAxMLargList + ['-f', 'b', '-t', MLtreeFile,
       '-z', TempAllBootstrappedTreesFile, '-n', MainTreeFile])
      assert ExitStatus == 0
    except:
      print('Problem in window', WindowAsStr, 'trying to collect all the',
      'bootstrapped trees onto the ML tree with RAxML. Skipping to the next',
      'window.', file=sys.stderr)
      return 1
    MainTreeFile = 'RAxML_bipartitions.' +MainTreeFile
    if not os.path.isfile(MainTreeFile):
      print(MainTreeFile +', expected to be produced by RAxML, does not '+\
      'exist.\nSkipping to the next window.', file=sys.stderr)
      return 1

    # Update on time taken if desired
    if UpdateTimes:
      TimesList.append(time.time())
      LastStepTime = TimesList[-1] - TimesList[-2]
      print('Bootstrapped trees in window', WindowAsStr,
      'finished. Number of seconds taken: ', LastStepTime)
  return 1




def TranslateSeqCoordsToAlnCoords(seq, coords):
  '''Takes a sequence that contains gaps (in general), and a set of coordinates
  specified with a respect to that sequence without gaps. The coordinates are
  translated to their positions in the gappy version of the sequence.
  e.g. called with the arguments "-a--cg-t-" and [1,2,3], we return [2,5,6].
  '''
  assert type(seq) == type('abc'), \
  'TranslateSeqCoordsToAlnCoords called with a sequence not of string type.'
  assert all(type(coord) == int for coord in coords), \
  'TranslateSeqCoordsToAlnCoords called with coords not of int type.'
  assert all(coord > 0 for coord in coords), \
  'TranslateSeqCoordsToAlnCoords called with at least one non-positive coord.'
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
  assert not -1 in TranslatedCoords, \
  'TranslateSeqCoordsToAlnCoords failed to find at least one coord.'
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
      print('Unexpected attribute properties for pysam.AlignedSegment\n', read,
      '\nSpecifically, expected equal numbers of bases, mapped positions and ',
      'base qualities, but found ', len(read.query_sequence), ', ', 
      len(positions), ' and ', len(read.query_qualities),
      ' respectively. Quitting.', sep='', file=sys.stderr)
      exit(1)
    return cls(read.query_name, read.query_sequence, positions,
    read.query_qualities)

  def __repr__(self):
    'Defining how a PseudoRead can be printed'
    return 'name: %s\nseq: %s\npositions: %s\nqualities: %s' % (self.name,
    self.sequence, ' '.join(map(str,self.positions)),
    ' '.join(map(str,self.qualities)))

  def SpansWindow(self, LeftWindowEdge, RightWindowEdge, ExactWindowStart,
    ExactWindowEnd):
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
    if ExactWindowStart and self.positions[LeftMostMappedBase] != \
    LeftWindowEdge:
      return False
    RightMostMappedBase = len(self.positions)-1
    while self.positions[RightMostMappedBase] == None:
      RightMostMappedBase -= 1
    if ExactWindowEnd and self.positions[RightMostMappedBase] != \
    RightWindowEdge:
      return False
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


  def ProcessRead(self, LeftWindowEdge, RightWindowEdge, MinQualForEnds,
  MinInternalQual, KeepOverhangs, RecoverClippedEnds, ExactWindowStart,
  ExactWindowEnd):
    '''Returns reads that span a given window.
    Overhangs & low-Q bases are trimmed if desired. None is returned for reads
    that do not span the window. The coordinates of the window edges should be
    zero-based.'''

    if RecoverClippedEnds:
      self.RecoverClippedEnds()

    # Skip reads that only partially overlap the window
    if not self.SpansWindow(LeftWindowEdge, RightWindowEdge, ExactWindowStart,
    ExactWindowEnd):
      return None

    # Trim low-Q ends if desired. Skip if that leaves only a partial overlap.
    if MinQualForEnds != None:
      self.QualityTrimEnds(MinQualForEnds)
      if not self.SpansWindow(LeftWindowEdge, RightWindowEdge, ExactWindowStart,
    ExactWindowEnd):
        return None

    # Skip reads containing more than one low-quality base
    if MinInternalQual != None and \
    (self.IsLowQual(MinInternalQual)):
      return None

    SeqToReturn = self.sequence

    # Now we trim the part of the read overhanging the window if desired.
    # If ExactWindowStart == True, we have already selected only those reads
    # whose first mapped position is exactly the left window edge. Ditto
    # ExactWindowEnd and the right window edge. If either of them are true, all
    # we want to do is trim unmapped bases off the read ends, no more. That's
    # because if ExactWindowStart == True && ExactWindowEnd == False, the
    # left-most mapped base is the start of what we want to keep and we don't
    # care where the read end is, ditto for the inverted statement, and finally
    # if both are true we do care about the start and end but we've already
    # preselected the reads to start and end at exactly the points of interest.
    if not KeepOverhangs:
      if ExactWindowStart or ExactWindowEnd:
        # Do this so that the only bases to get trimmed are unmapped ones at the
        # edges.
        LeftWindowEdge = float('-Inf')
        RightWindowEdge = float('Inf')
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
        print('Unexpected behaviour for read', self.name+', which',
        'maps to the following positions in the reference:\n'+ \
        ' '.join(map(str,self.positions)) +'\nUnable to determine ',
        'where the window edges ('+str(LeftWindowEdge+1), 'and',
        str(RightWindowEdge+1)+') are in this read. Skipping it.',
        file=sys.stderr)
        return None
      SeqToReturn = \
      SeqToReturn[LeftEdgePositionInRead:RightEdgePositionInRead+1]

    return SeqToReturn


  def MergeReadPairOverWindow(self, other, LeftWindowEdge, RightWindowEdge,
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
    # *including unmapped bases*. This is the position of the left-most mapped
    # base, minus the number of unmapped bases to its left. Find which of
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
    for j in range(
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
    MergedRead = PseudoRead(self.name, merged_sequence,
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

def MergeSimilarStringsB(DictOfStringCounts, SimilarityThreshold=1):
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
      print('The function MergeSimilarStringsB was called with a dict',
      "containing a key that's not a string.\nQuitting.", file=sys.stderr)
      exit(1)

  # Support duck-typing: the dict's values need not be numbers, we just need to 
  # be able add them together.
  try:
    TotalStringCount = sum(DictOfStringCounts.values())
  except TypeError:
    print('The function MergeSimilarStringsB was called with a dict',
    "containing values of types that cannot be added together.\nQuitting.",
    file=sys.stderr)
    exit(1)

  # Nothing needs to be done if the SimilarityThreshold is zero
  if SimilarityThreshold == 0:
    return DictOfStringCounts
      
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


def CalculateRecombinationMetric(SeqAlignment, NormaliseToDiversity, IncludeGaps=False):
  '''Considers all triplets of seqs and finds the maximum recombination signal.
  
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
  possible break points) and take the smaller of the two absolute values.
  If NormaliseToDiversity=True, we normalise by dividing by
  half the number of informative sites (i.e. ignoring sites where all sequences
  have the same base). This means that the maximum possible score of 1 is
  obtained if and only the two parents disagree at every informative site, the
  break point is exactly in the middle of all informative sites, and either side
  of the break point the recombinant agrees perfectly with one of the parents
  e.g.
  TATATATATA
  TATATATCTC
  TCTCTCTCTC
  If NormaliseToDiversity=False, we normalise by half the alignment length
  (ignoring sites that are wholly gaps). This means that the maximum possible
  score of 1 is obtained if and only if the two parents disagree at every site,
  the break point is exactly in the middle, and either side of the break point
  the recombinant agrees perfectly with one of the parents e.g.
  AAAAAAA
  AAAACCC
  CCCCCCC
  For either normaliser, after dividing the appropriate number of sites (either
  just informative ones, or all sites) by 2, we take the floor of the result
  i.e. round half-integers down. This means that the perfect score of 1 can be
  obtained with an odd number of sites.
  
  With the default value of False for the IncludeGaps argument, any position
  where any of the three sequences has the gap character '-' will be ignored.
  This means that e.g. the following three sequences would have a metric of
  zero: A-AAAA, A-AAA-A, AAAA-A. Setting IncludeGaps=True, the gap character
  will be treated the same as any other character, so that (dis)agreements in
  gaps count towards Hamming distance in exactly the same way as point
  mutations. This increases sensitivity of the metric to cases where indels are
  genuine signals of recombination, but decreases specificity, since
  misalignment may falsely suggest recombination.
  
  For speed, Hamming distances are calculated indirectly - looking only at
  informative sites, and considering only changes in distance each time the
  break point is slid through the next such site. However, runtime necessarily
  scales as N^3, where N is the number of sequences.

  The function returns a tuple of length four: (metric, ID of parent 1, ID of
  parent 2, ID of recombinant). If the metric is exactly zero, i.e. no
  recombination at all, the three sequence IDs will all be the 'None' value.
  '''

  MaxScoreAndSeqs = (0, None, None, None)

  # Recombination requires at least three sequences
  NumSeqs = len(SeqAlignment)
  if NumSeqs < 3:
    return MaxScoreAndSeqs

  # Strip alignment positions where every sequence has the same thing.
  # Also calculate the number of positions at which at least one sequence does
  # not have a gap.
  AlignmentLength = SeqAlignment.get_alignment_length()
  NumPureGapCols = 0
  for column in reversed(xrange(AlignmentLength)):
    RemoveThisCol = True
    FirstBaseSeen = None
    for base in SeqAlignment[:, column]:
      if FirstBaseSeen == None:
        FirstBaseSeen = base
      elif base != FirstBaseSeen:
        RemoveThisCol = False
        break
    if RemoveThisCol:
      SeqAlignment = SeqAlignment[:, :column] + SeqAlignment[:, column + 1:]
      if FirstBaseSeen == '-':
        NumPureGapCols += 1
  NumColsWithABase = AlignmentLength - NumPureGapCols
  NumInformativeSites = SeqAlignment.get_alignment_length()

  # Recombination requires sequences at least 2bp long
  ReducedAlignmentLength = SeqAlignment.get_alignment_length()
  if ReducedAlignmentLength < 2:
    return MaxScoreAndSeqs

  # Convert all the seqs to strings
  SeqsAsStrings = [str(seq.seq) for seq in SeqAlignment]

  # Calculate the score for sequences i and j being the two original seqs, and
  # k the recombinant between them.
  for i in range(NumSeqs):
    for j in range(i+1, NumSeqs):
      iSeq = SeqsAsStrings[i]
      jSeq = SeqsAsStrings[j]
      DisagreeingPositions = []
      for pos in range(ReducedAlignmentLength):
        if iSeq[pos] != jSeq[pos] and (IncludeGaps or
        (iSeq[pos] != '-' and jSeq[pos] != '-')):
          DisagreeingPositions.append(pos)

      NumDisagreeingPositions = len(DisagreeingPositions)
      if NumDisagreeingPositions == 0:
        continue
      for k in range(NumSeqs):
        if k == i or k == j:
          continue
        kSeq = SeqsAsStrings[k]

        # At each position where i and j disagree, record the 'loyalty' value 1
        # if k agrees with i, -1 if k agrees with j, and 0 otherwise.
        loyalties = []
        for pos in DisagreeingPositions:
          if kSeq[pos] == iSeq[pos]:
            loyalties.append(1)
          elif kSeq[pos] == jSeq[pos]:
            loyalties.append(-1)
          else:
            loyalties.append(0)

        # We consider all possible 'break points' for the list of loyalties,
        # i.e. points at which we split it into a left part and a right part.
        # The sum of the loyalties to the left or to the right equals (the
        # Hamming distance of i to k) minus (the Hamming distance of i to j) for
        # that part of the sequence (i.e. for the left or for the right).
        # We choose the break point that maximises the difference (between left
        # and right) of these two Hamming distance differences, i.e. the
        # difference between the sum of the loyalties to the left and the sum of
        # the loyalties to the right. For this break point we record the smaller
        # absolute value of the two loyalty sums.
        MaxLoyalty = 0
        MaxLoyaltyDiff = 0
        LoyaltyLeftOfBreakPoint = 0
        LoyaltyRightOfBreakPoint = sum(loyalties)
        for BreakPointPos in range(1, NumDisagreeingPositions):
          LoyaltyLeftOfBreakPoint  += loyalties[BreakPointPos - 1]
          LoyaltyRightOfBreakPoint -= loyalties[BreakPointPos - 1]
          #MeanLoyaltyLeft  = LoyaltyLeftOfBreakPoint / BreakPointPos
          #MeanLoyaltyRight = LoyaltyRightOfBreakPoint / \
          #(NumDisagreeingPositions - BreakPointPos)
          #LoyaltyDiff = abs(MeanLoyaltyLeft - MeanLoyaltyRight)
          LoyaltyDiff = abs(LoyaltyLeftOfBreakPoint - LoyaltyRightOfBreakPoint)
          if LoyaltyDiff > MaxLoyaltyDiff:
            MaxLoyaltyDiff = LoyaltyDiff
            MaxLoyalty = min(abs(LoyaltyLeftOfBreakPoint),
            abs(LoyaltyRightOfBreakPoint))
        if MaxLoyalty > MaxScoreAndSeqs[0]:
          MaxScoreAndSeqs = (MaxLoyalty, i, j, k)

  MaxScore = MaxScoreAndSeqs[0]
  if MaxScore == 0:
    return MaxScoreAndSeqs
  if NormaliseToDiversity:
    denominator = NumInformativeSites / 2
  else:
    denominator = NumColsWithABase / 2
  assert MaxScore <= denominator, \
  'Error: recombination metric > 1 found. Please report to Chris Wymant.'
  return (float(MaxScore) / denominator, ) + \
  tuple([SeqAlignment[i].id for i in MaxScoreAndSeqs[1:]])

def MakeBamIndices(BamFiles, SamtoolsCommand):
  '''Tries to run samtools index on bam files that don't have a .bai file.'''
  for BamFileName in BamFiles:
    if not os.path.isfile(BamFileName+'.bai'):
      IndexCommandPieces = [SamtoolsCommand, 'index', BamFileName]
      try:
        ExitStatus = subprocess.call(IndexCommandPieces)
        assert ExitStatus == 0
      except:
        print('Warning: encountered a problem running the command "' + \
        ' '.join(IndexCommandPieces) + '" (which we tried to run because',
        'there is no file', BamFileName + '.bai, i.e. it seems the bam file is',
        'not indexed). This may prevent the bam file from being readable later',
        'in the code. Continuing...', file=sys.stderr)


def MergeSimilarStringsA(DictOfStringCounts, SimilarityThreshold=1,
RecordCorrespondence=False):
  '''Absorbs those strings with lower counts into those with higher counts.
  TODO
  '''

  # Check that the keys of the dict are strings.
  for String in DictOfStringCounts:
    if type(String) != type('foo'):
      print('The function MergeSimilarStringsA was called with a dict',
      "containing a key that's not a string.\nQuitting.", file=sys.stderr)
      exit(1)

  # Sum the dict's values
  try:
    TotalStringCount = sum(DictOfStringCounts.values())
  except TypeError:
    print('The function MergeSimilarStringsA was called with a dict',
    "containing values of types that cannot be added together.\nQuitting.",
    file=sys.stderr)
    exit(1)


  if RecordCorrespondence:
    AfterToBeforeDict = \
    {string:[string] for string in DictOfStringCounts.keys()}

  # Nothing needs to be done to dicts with fewer than two entries, or if the
  # SimilarityThreshold is zero
  NumberOfUniqueStrings = len(DictOfStringCounts)
  if NumberOfUniqueStrings < 2 or SimilarityThreshold == 0:
    if RecordCorrespondence:
      return DictOfStringCounts, AfterToBeforeDict
    else:
      return DictOfStringCounts

  # Sort the strings by their counts
  SortedDict = sorted(DictOfStringCounts.items(), key=lambda x: x[1])

  PositionsOfStringsThatGetAbsorbed = set([])
  MergedDict = {}

  for i in range(NumberOfUniqueStrings - 1, -1, -1):

    if i in PositionsOfStringsThatGetAbsorbed:
      continue

    CommonString, CommonCount = SortedDict[i]
    CountForJsToMergeToThisI = 0

    for j in range(0, i):

      if j in PositionsOfStringsThatGetAbsorbed:
        continue

      RareString, RareCount = SortedDict[j]
      if RareCount == CommonCount:
        break

      # Compare the two strings. Initialise the number of differences as the
      # difference in string length, then do a pairwise comparison of characters
      # over the length of the shorter string. As soon as the number of 
      # differences exceeds the threshold, we know we shouldn't be merging them. 
      NumDifferingBases = abs(len(RareString) - len(CommonString))
      if NumDifferingBases > SimilarityThreshold:
        continue
      for base1, base2 in itertools.izip(RareString, CommonString):
        if base1 != base2:
          NumDifferingBases += 1
          if NumDifferingBases > SimilarityThreshold:
            break
      if NumDifferingBases <= SimilarityThreshold:
        CountForJsToMergeToThisI += RareCount
        PositionsOfStringsThatGetAbsorbed.add(j)
        if RecordCorrespondence:
          AfterToBeforeDict[CommonString].append(RareString)
          del AfterToBeforeDict[RareString]

    MergedDict[CommonString] = CommonCount + CountForJsToMergeToThisI

  # Check the total count is the same before and after.
  TotalStringCountAfterMerging = sum(MergedDict.values())
  if TotalStringCountAfterMerging != TotalStringCount:
    print('FATAL ERROR: the number of strings after merging, '+\
    str(TotalStringCountAfterMerging)+', != the number before merging, '+\
    str(TotalStringCount)+'. Please report this to Chris Wymant. Quitting.',
    file=sys.stderr)
    exit(1)

  if RecordCorrespondence:
    return MergedDict, AfterToBeforeDict
  else:
    return MergedDict
