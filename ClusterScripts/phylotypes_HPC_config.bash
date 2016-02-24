#!/usr/bin/env bash

# This file gets sourced (i.e. the values of the variables specified here are
# used) by phylotypes_HPC.bash, if you put this file in the place where
# phylotypes_HPC.bash is expecting it to be (see the top of that file).

# Where does your copy of phylotypes.py live?
PhylotypesCode="$HOME/PhylotypesCode/phylotypes.py"

MergingParam=3
MinReadCount=2

# Two files listing the bam and reference files, respectively, to be used.
# In both cases, no two files should have the same basename (i.e. the file name
# minus the path).
ListOfBamFiles="/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/ptyruns/ptyr2_bam.txt"
ListOfRefFiles="/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/ptyruns/ptyr2_ref.txt"

#ExtraArgs='-Q1 18 -Q2 2 -P -T --keep-overhangs --x-raxml raxml -D'
ExtraArgs='-T --x-raxml raxml -D'

# NB if you include --auto-window-params in the ExtraArgs, the following file is 
# not needed and the following bool should be set to false:
# A file in which each line is an even number of comma-separated integers 
# specifying (alternately) left- and right-hand edges of a window, in the usual
# phylotypes argument format; each line will run as its own job.
# To create such a file easily from a list of space-separated numbers all on the
# same line, see this example:
# IntsPerLine=2; count=0; for i in 1 100 101 200 201 300; do echo -n "$i "; 
# count=$((count+1)); if (( $count % $IntsPerLine == 0 )); then echo; fi;
# done > MyWindowFile.txt
WindowFile="$HOME/JobInputs/PhylotypesWindows.txt"
UseWindowFile=true

OutputDir="$WORK/PhylotypesOutput"

################################################################################
# Some temporary intermediate files we'll create:
ListOfLocalBams='temp_ListOfLocalBams.txt'
ListOfLocalRefs='temp_ListOfLocalRefs.txt'
