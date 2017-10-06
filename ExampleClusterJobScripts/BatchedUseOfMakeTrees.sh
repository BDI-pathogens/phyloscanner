#!/usr/bin/env bash

#PBS -l walltime=9:59:59
#PBS -l select=1:ncpus=1:mem=6000MB
#PBS -J 1-27

# Above are the job specs which you'll obviously need to change; -J should run 
# from 1 to the number of lines in the FileListingAllBatchPairs variable below.

# Splitting up a large set of bams into N batches and then doing one
# phyloscanner_make_trees.py run for each possible pair of batches means that 
# every possible pair of patients has been considered in at least one run (in 
# many different runs if they're in the same batch), and each run contains a 
# smaller total number of patients, which may improve MAFFT and RAxML accuracy. 
# The total computational burden will increase, but slower than linearly with N: 
# on the one hand each patient is now being processed again and again in N-1 
# different runs, but on the other hand because RAxML runtime increases faster
# than linearly with the number of sequences, there's a gain from splitting
# the total set of sequences into smaller groups that are processed separately.

# TODO: work out a modification of phyloscanner_analyse_trees.R that averages
# results for individual patients (such as their root-to-tip distance) over all
# runs in which that patient appears, and averages results for pairs of patients
# (namely their topological & distance relationships in each window) over all
# runs in which that pair appears.

## Before using this job script, you need to have defined your batches of bams.
## Say MyUsualBamRefIDlist.csv is your usual phyloscanner_make_trees.py input
## file, i.e. csv format with columns bam file name, reference file name, ID.
## Splitting this file into a number of separate files each with the same number
## of lines is trivial, e.g. splitting to size 100 with
#split MyUsualBamRefIDlist.csv --lines 100 batch_
## However in my example, IDs are either the patient ID, or of the format
## PatientID_SomethingElse where SomethingElse distinguishes one of a number of
## bam files associated with the same patient. I want batches with equal numbers
## of patients (not necessarily equal numbers of bams). First of all, sort the
## input file by ID:
#sort -t, -k3,3 MyUsualBamRefIDlist.csv > MyUsualBamRefIDlist_sorted.csv
## Then split up the sorted list into batches of 100 patients thus:
#BatchSize=100;
#LastID='';
#CurrentBatch=1;
#CurrentPosInBatch=0;
#while read line; do
#  ID=$(echo $line | awk -F, '{print $3}' | awk -F_ '{print $1}');
#  if [[ "$ID" == "$LastID" ]]; then
#    echo $line >> batch_${CurrentBatch}.csv;
#    continue;
#  fi;
#  CurrentPosInBatch=$((CurrentPosInBatch+1));
#  if [[ $CurrentPosInBatch -gt $BatchSize ]]; then
#    CurrentBatch=$((CurrentBatch+1));
#    CurrentPosInBatch=1;
#  fi;
#  if [[ $CurrentPosInBatch -eq 1 ]]; then
#    echo $line > batch_${CurrentBatch}.csv;
#  else
#    echo $line >> batch_${CurrentBatch}.csv;
#  fi
#  LastID="$ID";
#done < MyUsualBamRefIDlist_sorted.csv
## This creates a bunch of files named batch_1.csv, batch_2.csv etc. each of
## which is a phyloscanner-style input file for 100 patients (except the last 
## file which will have less than 100 of the total number of patients is not a
## multiple of 100). I then combine all possible pairs of these files into a new
## set of files, named batches_1_and_2.csv, batches_1_and_3.csv etc. thus:
#NumBatches=$(ls batch_*.csv | wc -l)
#for i in $(seq 1 $NumBatches); do
#  for j in $(seq $((i+1)) $NumBatches); do
#    cat batch_$i.csv batch_$j.csv > batches_${i}_and_${j}.csv
#  done
#done
## I then make one master file each line of which is the name of one of those
## batch-pair files:
#ls batches_*_and_*.csv > FileListingAllBatchPairs.txt
## The variable FileListingAllBatchPairs just below should be set to this file.
## This array-job script will run one array element for each line in that file,
## i.e. one phyloscanner run for each pair of batches.

FileListingAllBatchPairs=~/JobInputs/AllBEEHIVEPhyloscannerBatches.txt

# A file that contains, all in one line, the comma-separated list of window
# coordinates to be used as the argument to --windows:
WindowFileAllInOne="$HOME/JobInputs/HXB2_NewWindows_320w_160i_AllInOne.txt"

PhyloscannerCode="$HOME/phyloscanner/phyloscanner_make_trees.py"

# TODO: adjust the raxml binary and -T to match the ncpus specified
raxmlargs='raxmlHPC-SSE3 -m GTRCAT -p 1 --no-seq-check'

ExtraArgs="-Q1 25 -Q2 25 -P -A /home/cw109/JobInputs/2refs_HXB2_C.BW.fasta -2 B.FR.83.HXB2_LAI_IIIB_BRU.K03455 -XR B.FR.83.HXB2_LAI_IIIB_BRU.K03455 -XC "\
'823,824,825,892,893,894,907,908,909,1012,1013,1014,1156,1157,1158,1384,1385,1386,1444,1445,1446,1930,1931,1932,1957,1958,1959,2014,2015,2016,2023,2024,2025,2080,2081,2082,2134,2135,2136,2191,2192,2193,2280,2281,2282,2283,2284,2285,2298,2299,2300,2310,2311,2312,2316,2317,2318,2319,2320,2321,2322,2323,2324,2340,2341,2342,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2373,2374,2375,2379,2380,2381,2385,2386,2387,2388,2389,2390,2391,2392,2393,2394,2395,2396,2400,2401,2402,2409,2410,2411,2412,2413,2414,2415,2416,2417,2424,2425,2426,2430,2431,2432,2436,2437,2438,2439,2440,2441,2442,2443,2444,2457,2458,2459,2460,2461,2462,2463,2464,2465,2469,2470,2471,2472,2473,2474,2478,2479,2480,2481,2482,2483,2496,2497,2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,2514,2515,2516,2517,2518,2519,2520,2521,2522,2526,2527,2528,2529,2530,2531,2535,2536,2537,2670,2671,2672,2679,2680,2681,2703,2704,2705,2709,2710,2711,2733,2734,2735,2742,2743,2744,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2769,2770,2771,2772,2773,2774,2778,2779,2780,2811,2812,2813,2814,2815,2816,2817,2818,2819,2823,2824,2825,2841,2842,2843,2847,2848,2849,2850,2851,2852,2856,2857,2858,2865,2866,2867,2871,2872,2873,2892,2893,2894,2895,2896,2897,2901,2902,2903,2904,2905,2906,2952,2953,2954,2961,2962,2963,3000,3001,3002,3015,3016,3017,3018,3019,3020,3030,3031,3032,3042,3043,3044,3084,3085,3086,3090,3091,3092,3099,3100,3101,3111,3112,3113,3117,3118,3119,3135,3136,3137,3171,3172,3173,3177,3178,3179,3180,3181,3182,3189,3190,3191,3192,3193,3194,3204,3205,3206,3210,3211,3212,3222,3223,3224,3228,3229,3230,3237,3238,3239,3246,3247,3248,3249,3250,3251,3255,3256,3257,3261,3262,3263,3396,3397,3398,3501,3502,3503,3546,3547,3548,3705,3706,3707,4425,4426,4427,4449,4450,4451,4503,4504,4505,4518,4519,4520,4590,4591,4592,4641,4642,4643,4647,4648,4649,4656,4657,4658,4668,4669,4670,4671,4672,4673,4692,4693,4694,4722,4723,4724,4782,4783,4784,4974,4975,4976,5016,5017,5018,5067,5068,5069,7863,7864,7865,7866,7867,7868,7869,7870,7871,7872,7873,7874,7875,7876,7877,7881,7882,7883,7884,7885,7886'\
" --merging-threshold-a 1 --min-read-count 2 -T"

# We'll make a subdirectory in here for each batch pair:
OutputBaseDir=$WORK/PhyloscannerOutput_BEEHIVE


################################################################################
# INITIALISATION

# Exit this script if an undefined variable is encountered:
set -u

# Load modules
module load anaconda/2.3.0 &&
module load samtools &&
module load mpi &&
module load raxml/8.2.9 &&
module load mafft/7 || \
{ echo 'Failed to load required modules. Quitting.' >&2 ; exit 1;  }


# Check required files exist
for i in "$PhyloscannerCode" "$FileListingAllBatchPairs" "$WindowFileAllInOne"; do
  if [ ! -f "$i" ]; then
    echo "$i" 'does not exist. Quitting.' >&2
    exit 1
  fi
done

# Each line in the input file is a batch pair to analyse, so count the number
# of lines.
NumBatchPairs=$(wc -l "$FileListingAllBatchPairs" | awk '{print $1}')

# Check that the current array index is not larger than the number of input
# files.
if [ "$PBS_ARRAY_INDEX" -gt "$NumBatchPairs" ]; then
  echo 'Error: job array index' "$PBS_ARRAY_INDEX" 'is greater than the number'\
  'of lines in' "$FileListingAllBatchPairs"'. Quitting.' >&2
  exit 1
fi

# Find the window(s) corresponding to this array index.
BatchPairFile=$(sed -n "$PBS_ARRAY_INDEX"'p' "$FileListingAllBatchPairs")
BatchPairID=$(basename "${BatchPairFile%.csv}")

# Check required files exist
if [ ! -f "$BatchPairFile" ]; then
  echo "$BatchPairFile" 'does not exist. Quitting.' >&2
  exit 1
fi

# Make the output directory if it doesn't exist
OutputDir="$OutputBaseDir/$BatchPairID"
if [ ! -d "$OutputDir" ]; then
  mkdir -p "$OutputDir"
fi || { echo 'Unable to create the specified output directory. Quitting.' >&2 ;\
exit 1; }

################################################################################


"$PhyloscannerCode" "$BatchPairFile" \
-W $(cat "$WindowFileAllInOne") \
$ExtraArgs \
--x-raxml "$raxmlargs" \
--output-dir "$OutputDir" \
> "$BatchPairID"Log.out 2> "$BatchPairID"Log.err

cp "$BatchPairID"Log.out "$BatchPairID"Log.err -t "$OutputDir"

