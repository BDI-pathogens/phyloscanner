#!/bin/bash

# This script currently expects folders named ReadNames and BAMs

# saner programming env: these switches turn some bugs into errors
set -o errexit -o pipefail -o noclobber -o nounset

! getopt --test > /dev/null 
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo "I'm sorry, `getopt --test` failed in this environment."
    exit 1
fi

# Where does phyloscanner live normally?

phyDir="~/phyloscanner"

# Most of these are to be passed, but a couple are unique and this script also needs to know about a few of them

fullOpts="$@"

OPTIONS=m:b:v:x:y:p:i:
LONGOPTS=phyloscannerPath:,baminput:,og:,outgroupName:,multifurcationThreshold:,blacklist:,od:,outputDir:,verbose:,tipRegex:,tfe:,treeFileExtension:,cfe:,csvFileExtension:,sd:,seed:,ow,overwrite,rda,outputRDA,nr:,normRefFileName:,ns,normStandardiseGagPol,nc:,normalisationConstants:,db:,duplicateBlacklist:,pbk:,parsimonyBlacklistK:,rwt:,rawBlacklistThreshold:,rtt:,ratioBlacklistThreshold:,ub,dualBlacklist,rcm,readCountsMatterOnZeroLengthBranches,blr,blacklistReport,dsl:,maxReadsPerHost:,dsb,blacklistUnderrepresented


! PARSED=$(getopt -a --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    exit 2
fi
# read getopt's output this way to handle the quoting right:

eval set -- "$PARSED"

outDir="."
overwrite=false
# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
        -p|--phyloscannerPath)
            phyDir="$2"
            shift 2
            ;;
        -i|--baminput)
            bamInputFile="$2"
            shift 2
            ;;
        --od|--outputDir)
            outDir="$2"
            shift 2
            ;;
        --ow|--overwrite)
            overwrite=true
            shift
            ;;
        --ow|--overwrite|--rda|--outputRDA|--ns|--normStandardiseGagPol|--ub|--dualBlacklist|--rcm|--readCountsMatterOnZeroLengthBranches|--blr|--blacklistReport|--dsb|--blacklistUnderrepresented)
        	shift
        	;;
	    --)
            shift
            break
            ;;
        *)
            shift 2
            ;;
    esac
done



# This script takes the same arguments as phyloscanner_clean_alignments.R except -i and -p. Delete these

pcaopts=$(echo $fullOpts sed s/-[ip]\ [^\ ]*\ //g)
pcaopts=$(echo $pcaopts sed s/--phyloscannerPath\ [^\ ]*\ //g)
pcaopts=$(echo $pcaopts sed s/--baminput\ [^\ ]*\ //g)

runName=$3

Rscript ${phyDir}/tools/phyloscanner_clean_alignments.R $pcaopts

mkdir -p ${outDir}/cleaning_temp

for wrf in ReadNames/ReadNames2_*
do
    fn=$(basename $wrf)
    windowstring=${fn:20}
    windowstring=${windowstring::-4}
    if [ "$overwrite" = true ]
    then
        python ${phyDir}/tools/FindAllNonBlacklistedReads.py ${outDir}/${runName}_blacklistReport.csv $wrf ${outDir}/cleaning_temp/${runName}_keptReads_${windowstring} --overwrite
    else 
        python ${phyDir}/tools/FindAllNonBlacklistedReads.py ${outDir}/${runName}_blacklistReport.csv $wrf ${outDir}/cleaning_temp/${runName}_keptReads_${windowstring}
    fi
done

Rscript ${phyDir}/tools/collect_kept_reads_by_BAM.R ${outDir}/cleaning_temp/

myregex="${runName}_keptReads_allWindows_(.*)\.txt"

mkdir -p ${outDir}/cleaned_BAMs

for mf in ${outDir}/cleaning_temp/*allWindows*
do
    f=$(basename $mf)

    if [[ $f =~ $myregex ]]
    then
        name="${BASH_REMATCH[1]}"
        echo "Cleaning ${name}.bam to ${name}_${runName}_cleaned.bam"
        python ${phyDir}/tools/ExtractNamedReadsFromBam.py BAMs/${name}.bam ${outDir}/cleaned_BAMs/${name}_${runName}_cleaned.bam -F $mf

    fi
done

rm -r ${outDir}/cleaning_temp
