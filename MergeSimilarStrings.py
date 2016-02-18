from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251

import sys
import itertools

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
  Note that by construction, calling the function twice - feeding the output
  back in as input - is exactly the same as calling it once (except slower).
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

