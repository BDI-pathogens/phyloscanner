#!/usr/bin/env bash

set -u
set -o pipefail

ExpectedNumArgs=2
if [ "$#" -ne "$ExpectedNumArgs" ]; then
    echo "This script should be run from the command line with"\
    "$ExpectedNumArgs arguments: the config file, and the directory containing"\
    "the RAxML_bestTree.* files (one per window) and the BamIDs.txt file"\
    "produced by phyloscanner. Quitting."
    exit 1
fi
ConfFile=$1
TreeDir=$2

if [[ ! -f "$ConfFile" ]]; then
  echo "$ConfFile does not exist or is not a regular file. Quitting."
  exit 1
fi
if [[ ! -d "$TreeDir" ]]; then
  echo "$TreeDir does not exist or is not a directory. Quitting."
  exit 1
fi

ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ToolsDir="$ThisDir"/tools
if [[ ! -d "$ToolsDir" ]]; then
  echo "Expected to find a subdirectory named tools in the same directory as"\
  "this script; found no such subdirectory. Quitting."
  exit 1
fi


source "$ConfFile"

TreeDir=$(cd "$TreeDir"; pwd)
ToolsDir=$(cd "$ToolsDir"; pwd)
PatientIDfile="$TreeDir"/BamIDs.txt
if [[ ! -f "$PatientIDfile" ]]; then
  echo "$PatientIDfile does not exist or is not a regular file. Quitting."
  exit 1
fi

################################################################################
# Files we'll produce
RoguesPrefix='RogueBlacklist_'
DualsPrefix='MultipleInfections_'
FinalBlacklistPrefix='FinalBlacklist_'
SubgraphsPrefix='subgraphs_'
SummaryPrefix='summary'
ClassPrefix='Classification_'
TransmissionSummary='TransmissionSummary.csv'
NormalisationReference='groupM_reference_trees_grubbs1_stats.csv'
RawNormalisationLookup='raw_normalisations.csv'
ProcessedNormalisationLookup='processed_normalisations.csv'
################################################################################

# Install any missing packages

Rscript "$ToolsDir"/PackageInstall.R

Rscript "$ToolsDir"/NormalisationLookupWriter.R "$TreeDir"/'RAxML_bestTree.' "$NormalisationReference" "$RawNormalisationLookup" "MEDIAN_PWD" -D "$ToolsDir" --standardize

# Find rogue reads and, if desired, reads that look like they're part of a dual
# infection.
if [[ "$ExcludeDuals" == "true" ]]; then
  Rscript "$ToolsDir"/ParsimonyBasedBlacklister.R "$SubgraphMinCount" \
  "$SubgraphMinRatio" "$SankhoffK" "$TreeDir"/'RAxML_bestTree.' "$RoguesPrefix" -x "$regex" -D \
  "$ToolsDir" -r "$root" -d "$DualsPrefix" -n "$RawNormalisationLookup" || { echo \
  'Problem running ParsimonyBasedBlacklister.R. Quitting.' ; exit 1 ; }
else
  Rscript "$ToolsDir"/ParsimonyBasedBlacklister.R "$SubgraphMinCount" \
  "$SubgraphMinRatio" "$SankhoffK" "$TreeDir"/'RAxML_bestTree.' "$FinalBlacklistPrefix" -x "$regex" -D \
  "$ToolsDir" -r "$root" -n "$RawNormalisationLookup" || { echo \
  'Problem running ParsimonyBasedBlacklister.R. Quitting.' ; exit 1 ; }
fi

# Find patients who look dual in enough windows, and add all of their reads from
# all windows to the blacklists, IF we're removing duals.
if [[ "$ExcludeDuals" == "true" ]]; then
  Rscript "$ToolsDir"/DualPatientBlacklister.R $FractionOfWindowsToCallDual \
  "$TreeDir"/'RAxML_bestTree.' "$DualsPrefix" "$FinalBlacklistPrefix" -b "$RoguesPrefix" -D "$ToolsDir" || { echo \
  'Problem running DualPatientBlacklister.R. Quitting.' ; exit 1 ; }
fi

# Split patients into their subgraphs
Rscript "$ToolsDir"/SplitPatientsToSubgraphs.R "$TreeDir"/'RAxML_bestTree.' "$RunLabel" -R -r "$root" -b "$FinalBlacklistPrefix" -x "$regex" -s "$SplitsRule" -t "$TiesRule" -k "$SankhoffK" -m "$MultifurcationThreshold" -D "$ToolsDir" -n "$RawNormalisationLookup" -pw 20 -ph 0.5 || { echo \
  'Problem running SplitPatientsToSubgraphs.R. Quitting.' ; exit 1 ; }

Rscript "$ToolsDir"/NormalisationLookupWriter.R "$TreeDir"/'ProcessedTree_' "$NormalisationReference" "$ProcessedNormalisationLookup" "MEDIAN_PWD" -D "$ToolsDir" --standardize

# Generate summary stats over all windows
Rscript "$ToolsDir"/SummaryStatistics.R "$PatientIDfile" 'ProcessedTree_'"$SplitsRule"'_'"$RunLabel" "$SubgraphsPrefix$SplitsRule"'_'"$RunLabel" \
"$SummaryPrefix" -b "$FinalBlacklistPrefix" -x "$regex" -D "$ToolsDir" || { echo \
  'Problem running SummaryStatistics.R. Quitting.' ; exit 1 ; }

# Classify relationships between patients in each window
Rscript "$ToolsDir"/ClassifyRelationships.R 'ProcessedTree_'"$SplitsRule"'_'"$RunLabel" "$SubgraphsPrefix$SplitsRule"'_'"$RunLabel" "$ClassPrefix$SplitsRule" -c -D "$ToolsDir" -n "$ProcessedNormalisationLookup" || { echo \
  'Problem running ClassifyRelationships.R. Quitting.' ; exit 1 ; }

# Summarise relationships across all windows
Rscript "$ToolsDir"/TransmissionSummary.R "$PatientIDfile" "$ClassPrefix$SplitsRule"'_classification_' "$TransmissionSummary" -D "$ToolsDir" -s "$SummaryPrefix"'_patStatsFull.csv' -m "$MinWindowsForTransmissionLink" -c "$MaxDistanceForTransmissionLink" || { echo \
  'Problem running TransmissionSummary.R. Quitting.' ; exit 1 ; }


