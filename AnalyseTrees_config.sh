#!/bin/bash

root='C.BW.00.00BW07621.AF443088' # The name of the reference seq used for rooting
RunLabel='MyRun' # a label for this particular phyloscanner run
ExcludeDuals=false # Exclude patients who seem to be dually infected?


# Not worth changing these for the moment:

MinWindowsForTransmissionLink=5
MaxDistanceForTransmissionLink=0.05
SubgraphMinCount=3
SubgraphMinRatio=0.005
SankhoffK=20
SankhoffP=0
FractionOfWindowsToCallDual=0.15
MultifurcationThreshold=1E-5
#s=Sankhoff (slow, rigorous), r=Romero-Severson (quick, less rigorous with >2
# patients):
SplitsRule='s'   
TiesRule='b'
regex="^(.*)_read_([0-9]+)_count_([0-9]+)$" # use this unless you need to merge different bams into one
