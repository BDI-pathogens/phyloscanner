#!/usr/bin/env bash

# This file gets sourced (i.e. the values of the variables specified here are
# used) by phylotypes_HPC.bash, if you put this file in the place where
# phylotypes_HPC.bash is expecting it to be (see the top of that file).

# Where does your copy of phylotypes.py live?
PhylotypesCode="$HOME/PhylotypesCode/phylotypes.py"

MergingParam=1
MinReadCount=2

# Two files listing the bam and reference files, respectively, to be used.
# In both cases, no two files should have the same basename (i.e. the file name
# minus the path).
#ListOfBamFiles="/home/cw109/JobInputs/test2_BamFiles.txt"
#ListOfRefFiles="/home/cw109/JobInputs/test2_RefFiles.txt"
ListOfBamFiles="/home/cw109/JobInputs/LongReadBamFiles_BEEHIVE.txt"
ListOfRefFiles="/home/cw109/JobInputs/LongReadBamFiles_BEEHIVE_refs.txt"


ExtraArgs='--x-raxml raxml -Q1 23 -Q2 23 -P -A /home/cw109/JobInputs/B.FR.83.HXB2_LAI_IIIB_BRU.K03455.fasta --renaming-file /home/cw109/JobInputs/LongReadBamFiles_BEEHIVE_PatientIDs.txt'
#ExtraArgs='--x-raxml raxml -Q1 23 -Q2 23 -P -A /home/cw109/JobInputs/B.FR.83.HXB2_LAI_IIIB_BRU.K03455.fasta --renaming-file /home/cw109/JobInputs/test2_PatientIDs.txt'

# NB if you include --auto-window-params in the ExtraArgs, the following file is 
# not needed and the following bool should be set to false:
# A file in which each line is an even number of integers specifying
# (alternately) left- and right-hand edges of a window, in the usual phylotypes
# argument format; each line will run as its own job.
# To create such a file easily from a list of space-separated numbers all on the
# same line, see this example:
# IntsPerLine=2; count=0; for i in 1 100 101 200 201 300; do echo -n "$i "; 
# count=$((count+1)); if (( $count % $IntsPerLine == 0 )); then echo; fi;
# done > MyWindowFile.txt
WindowFile="$HOME/JobInputs/AllLongMiseqPatientsAlignedToHXB2_SmartWindows.txt"
UseWindowFile=true

OutputDir="$WORK/PhylotypesOutput_BEEHIVE"

################################################################################
# Some temporary intermediate files we'll create:
ListOfLocalBams='temp_ListOfLocalBams.txt'
ListOfLocalRefs='temp_ListOfLocalRefs.txt'
