

for a in `seq 43` 

do
let "STARTWINDOW = 600 + a*200"
let "ENDWINDOW = STARTWINDOW + 300"

BlacklistInputFile="/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160211/DuplicateReads_surviving_InWindow_${STARTWINDOW}_to_${ENDWINDOW}.csv"
BlacklistOutputFile="/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160211/Blacklist_InWindow_${STARTWINDOW}_to_${ENDWINDOW}.csv"

# This will be created
TreeInputFile="/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160211/RAxML_bestTree.InWindow_${STARTWINDOW}_to_${ENDWINDOW}.tree"
TransOutputFile="/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160211/LikelyTransmissions_${STARTWINDOW}_to_${ENDWINDOW}.csv"

echo "Making blacklist"

Rscript /Users/twoseventwo/Documents/phylotypes/MakeBlacklist.R "$BlacklistInputFile" "$BlacklistOutputFile"

echo "Running LikelyTransmissions"

if [ -e "$BlacklistOutputFile" ]
then
echo "Found blacklist"
Rscript /Users/twoseventwo/Documents/phylotypes/LikelyTransmissions.R -f "$TreeInputFile" -o "$TransOutputFile" -r "C.BW.00.00BW07621.AF443088" -p 1 -s "_" -b "$BlacklistOutputFile" -t 0.08 
fi
echo "No blacklist found"
Rscript /Users/twoseventwo/Documents/phylotypes/LikelyTransmissions.R -f "$TreeInputFile" -o "$TransOutputFile" -r "C.BW.00.00BW07621.AF443088" -p 1 -s "_" -t 0.08


done
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

