#!/usr/bin/env bash

#PBS -l walltime=0:29:59
#PBS -l select=1:ncpus=1:mem=3GB
#PBS -J 1-155

# This script requires input variables chosen by the user, and stored in the
# following file - make sure it exists with your chosen values! The PBJ -J
# directive above should run from 1 to the number of lines in the "$WindowFile"
# variable you set there.
source "$HOME/phylotypes_HPC_config.bash"


################################################################################
# INITIALISATION

# Exit this script if an undefined variable is encountered:
set -u

# Load modules
module load anaconda/2.3.0 &&
module load samtools &&
module load mpi &&
module load raxml/8.2.4 &&
module load mafft/7 || \
{ echo 'Failed to load required modules. Quitting.' >&2 ; exit 1;  }

# Check required files exist
for i in "$PhylotypesCode" "$ListOfBamFiles" "$ListOfRefFiles" "$WindowFile"; do
  if [ ! -f "$i" ]; then
    echo "$i" 'does not exist. Quitting.' >&2
    exit 1
  fi
done

# Make the output directory if it doesn't exist
if [ ! -d "$OutputDir" ]; then
  mkdir "$OutputDir"
fi || { echo 'Unable to create the specified output directory. Quitting.' >&2 ;\
exit 1; }

# Remove trailing slashes from user-specified directories, if present.
OutputDir=$(cd "$OutputDir"; pwd)

# Each line in the list of all windows is a set of windows, so count the number
# of lines.
NumSetsOfWindows=$(wc -l "$WindowFile" | awk '{print $1}')

# Check that the current array index is not larger than the number of input
# files.
if [ "$PBS_ARRAY_INDEX" -gt "$NumSetsOfWindows" ]; then
  echo 'Error: job array index' "$PBS_ARRAY_INDEX" 'is greater than the number'\
  'of lines in' "$WindowFile"'. Quitting.' >&2
  exit 1
fi

# Find the window(s) corresponding to this array index.
WindowsForThisJob=$(sed -n "$PBS_ARRAY_INDEX"'p' "$WindowFile")

# Check all bam files and ref files exist, copy them to the current working dir,
# and create new lists of bam and ref files that point to the local copies.
function MakeRemoteFilesLocal {
  RemoteFileList="$1"
  LocalFileList="$2"
  echo -n > "$LocalFileList"
  while read File; do
    if [ ! -f "$File" ]; then
      echo "$File" 'does not exist. Quitting.' >&2
      exit 1
    fi
    cp "$File" .
    basename "$File" >> "$LocalFileList"
  done < "$RemoteFileList"
}
MakeRemoteFilesLocal "$ListOfBamFiles" "$ListOfLocalBams"
MakeRemoteFilesLocal "$ListOfRefFiles" "$ListOfLocalRefs"

################################################################################


# THE MAIN CODE:

"$PhylotypesCode" $MergingParam $MinReadCount "$ListOfLocalBams" \
"$ListOfLocalRefs" $WindowsForThisJob $ExtraArgs

cp \
AlignedReadsInWindow*.fasta \
*_ref.fasta \
DiscardedReads*.bam \
RAxML_* \
-t "$OutputDir"
