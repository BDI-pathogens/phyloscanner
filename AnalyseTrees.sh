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
################################################################################


# Find rogue reads and, if desired, reads that look like they're part of a dual
# infection.
if [[ "$ExcludeDuals" == "true" ]]; then
  Rscript "$ToolsDir"/ParsimonyBasedBlacklister.R "$SubgraphMinCount" \
  "$SubgraphMinRatio" "$SankoffK" "$TreeDir"/'RAxML_bestTree.' "$RoguesPrefix" -x "$regex" -D \
  "$ToolsDir" -r "$root" -d "$DualsPrefix"
else
  Rscript "$ToolsDir"/ParsimonyBasedBlacklister.R "$SubgraphMinCount" \
  "$SubgraphMinRatio" "$SankoffK" "$TreeDir"/'RAxML_bestTree.' "$FinalBlacklistPrefix" -x "$regex" -D \
  "$ToolsDir" -r "$root"
fi

# Find patients who look dual in enough windows, and add all of their reads from
# all windows to the blacklists, IF we're removing duals.
if [[ "$ExcludeDuals" == "true" ]]; then
  Rscript "$ToolsDir"/DualPatientBlacklister.R $FractionOfWindowsToCallDual \
  "$TreeDir"/'RAxML_bestTree.' "$DualsPrefix" "$FinalBlacklistPrefix" -b "$RoguesPrefix" -D \
  "$ToolsDir"
fi

# Split patients into their subgraphs
Rscript "$ToolsDir"/SplitPatientsToSubgraphs.R "$TreeDir"/'RAxML_bestTree.' "$RunLabel" -r "$root" -b "$FinalBlacklistPrefix" -x "$regex" -s "$SplitsRule" -k "$SankoffK" -D "$ToolsDir" -pw 20 -ph 0.5

# Generate summary stats over all windows
Rscript "$ToolsDir"/SummaryStatistics.R "$PatientIDfile" 'ProcessedTree_'"$SplitsRule"'_'"$RunLabel" "$SubgraphsPrefix$SplitsRule"'_'"$RunLabel" \
"$SummaryPrefix" -b "$FinalBlacklistPrefix" -x "$regex" -D "$ToolsDir"

# Classify relationships between patients in each window
Rscript "$ToolsDir"/NewClassifyRelationships.R 'ProcessedTree_'"$SplitsRule"'_'"$RunLabel" "$SubgraphsPrefix$SplitsRule"'_'"$RunLabel" "$ClassPrefix$SplitsRule" -c -D "$ToolsDir"

# Summarise relationships across all windows
Rscript "$ToolsDir"/NewTransmissionSummary.R "$PatientIDfile" "$ClassPrefix$SplitsRule"'_classification_' "$TransmissionSummary" -D "$ToolsDir" -s "$SummaryPrefix"'_patStatsFull.csv' -m "$MinWindowsForTransmissionLink" -c "$MaxDistanceForTransmissionLink"


