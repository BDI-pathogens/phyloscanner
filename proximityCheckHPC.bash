#!/usr/bin/env bash

#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=1:mem=6GB
#PBS -q pqeelab
#PBS -k oe

# The above commands say this script should not take longer than 2 seconds,
# requests 1 node with 1 CPU per node and 3mb memory per node.
# The three main directories for use are $HOME, $WORK (for your output) and
# $TMPDIR (for intermediate files that you don't want to keep).

# Put this near the top of all your bash scripts, always. It means that an
# undefined variable raises an error, instead of being assumed to be null.
set -u

module load intel-suite 
module load R

################################################################################
# INPUT

# An input file you want for processing. NB you should ensure this exists before
# submitting this script!

let "STARTWINDOW = 600 + PBS_ARRAY_INDEX*200"
let "ENDWINDOW = STARTWINDOW + 300"

BlacklistInputFile="$HOME/phylotypes/run20160211/DuplicateReads_surviving_InWindow_${STARTWINDOW}_to_${ENDWINDOW}.csv"
BlacklistOutputFile="$WORK/Blacklist_InWindow_${STARTWINDOW}_to_${ENDWINDOW}.csv"

# This will be created
TreeInputFile="$HOME/phylotypes/run20160211/RAxML_bestTree.InWindow_${STARTWINDOW}_to_${ENDWINDOW}.tree"
ProxOutputFile="$WORK/ProximityCheck_${STARTWINDOW}_to_${ENDWINDOW}.csv"

OutputDirectory="$WORK"

################################################################################



################################################################################
# INITIALISATION

# Try to make the output directory if it doesn't exist. Quit if we are unable to
# (e.g. if we don't have permission, or the specified output dir is a subdir of
# a dir that does not exist).

if [ ! -d "$OutputDirectory" ]; then
  mkdir "$OutputDirectory"
fi || { echo 'Unable to create the specified output directory. Quitting.' >&2 ;\
exit 1; }

################################################################################



################################################################################
# THE MAIN PROGRAM

echo "Making blacklist"

cd $TMPDIR
Rscript $HOME/phylotypes/MakeBlacklist.R "$BlacklistInputFile" "$BlacklistOutputFile"

echo "Running LikelyTransmissions"

if [ -e "$BlacklistOutputFile" ]
then
echo "Found blacklist"
Rscript $HOME/phylotypes/ProximityCheck.R -f "$TreeInputFile" -o "$ProxOutputFile" -r "B.FR.83.HXB2_LAI_IIIB_BRU.K03455" -p 1 -s "_" -b "$BlacklistOutputFile" -t 0.08 
fi
echo "No blacklist found"
Rscript $HOME/phylotypes/ProximityCheck.R -f "$TreeInputFile" -o "$ProxOutputFile" -r "B.FR.83.HXB2_LAI_IIIB_BRU.K03455" -p 1 -s "_" -t 0.08 

rm "$BlacklistOutputFile"
cp "$ProxOutputFile" $WORK

################################################################################



################################################################################
# SOME EXTRA GUIDANCE (NOT PART OF THE SCRIPT)

## If your Imperial username is ab123, copy this file into your home directory
## thus (supplying your password when prompted):
# scp MyQsubTest.bash ab123@login.cx1.hpc.ic.ac.uk:/home/ab123/MyQsubTest.bash
## then login (supplying your password when prompted):
# ssh ab123@login.cx1.hpc.ic.ac.uk
## See that this script has been copied here:
# ls
## Submit this job:
# qsub MyQsubTest.bash
## Then wait for $OutputFile to be created wherever you specified, and stdout
## and stderr files to be created in the directory you were in when you 
## submitted your job (I always do this from $HOME so I don't forget where 
## these files will be).
## While a job is queued or running, check its progress with
# qstat
## Abort a queued or running job (whose jobID will be displayed by qstat) with
# qdel jobID
##
## To copy everything from $WORK to a local directory /path/to/HPCoutput,
# rsync -az ab123@login.cx1.hpc.ic.ac.uk:'$WORK/' /path/to/HPCoutput
## where the -z flag uses compression to speed up the copy. rsync has many
## possible options; try
# man rsync
## (another example is --ignore-existing to skip files you've copied already).
##

