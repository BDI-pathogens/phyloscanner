#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
explanatory_message = '''This script analyses the files produced by
phyloscanner_make_trees.py's --read-names-2 option together with the file
phyloscanner_analyse_trees.R's --blacklistReport option. The former contain the
correspondence between a tip name (from one of the trees) and the reads that
went into that tip; the latter says which of these tips were blacklisted and
which were not. In general, for each read we have multiple results for whether
it is blacklisted or not, due to the read appearing in multiple windows. This
script collects together all of the blacklisting results for each read, and
decides whether the read should be kept or not. We output a list of kept reads,
one per bam; you can then extract those reads from their bam by giving that list
as the --read-name-file argument of
~/phyloscanner/tools/ExtractNamedReadsFromBam.py.'''

import argparse
import os
import sys
import re
from collections import defaultdict, Counter

blacklist_window_regex_string = "(\d+_to_\d+)$"
blacklist_window_regex = re.compile(blacklist_window_regex_string)
read_names_window_regex_string = "ReadNames2_InWindow_(\d+_to_\d+).csv$"
read_names_window_regex = re.compile(read_names_window_regex_string)
tip_regex_string = "_read_\d+_count_\d+$"
tip_regex = re.compile(tip_regex_string)

def read_blacklist_report(blacklist_report_file):

  blacklists_by_window = defaultdict(dict)

  with open(blacklist_report_file, 'r') as f:
    for line_num_min_1, line in enumerate(f):

      # Remove leading & trailing whitespace
      line = line.strip()

      # Check for the expected header line
      if line_num_min_1 == 0:
        expected_header = "tree.id,tip,kept,status"
        if line == expected_header:
          continue
        else:
          print("Error: unexpected header line '" + line + "' in",
          blacklist_report_file + "; expected '" + expected_header + \
          "'. Quitting.", file=sys.stderr)
          exit(1)

      # Split into fields
      fields = line.split(",")
      if len(fields) != 4:
        print("Error: expected four fields, found", len(fields), "on line",
        line_num_min_1 + 1, "in", blacklist_report_file + ". Quitting.",
        file=sys.stderr)
        exit(1)
      tree_id, tip, kept, status = fields

      # Find which window this line is for.
      tree_match = blacklist_window_regex.search(tree_id)
      if not tree_match:
        print("Error: the tree_id", tree_id, "on line", line_num_min_1 + 1,
        "in", blacklist_report_file + " does not match the expected pattern '"\
         + blacklist_window_regex_string + "'. Quitting.", file=sys.stderr)
        exit(1)
      window = tree_match.groups()[0]

      # kept should be a bool
      if kept == "TRUE":
        kept = True
      elif kept == "FALSE":
        kept = False
      else:
        print("Error: unexpected value", kept, "for 'kept' field on line",
        line_num_min_1 + 1, "in", blacklist_report_file + ". Quitting.",
        file=sys.stderr)
        exit(1)

      # The blacklist report in general will contain external references as well
      # reads; skip these.
      tip_match = tip_regex.search(tip)
      if not tip_match:
        continue

      # Check we haven't seen this tip before, then record it.
      if tip in blacklists_by_window[window]:
        print("Error: encountered tip", tip, "twice for window", window,
        "on line", line_num_min_1 + 1, "in", blacklist_report_file + \
        ". Quitting.", file=sys.stderr)
        exit(1)
      blacklists_by_window[window][tip] = kept

  return blacklists_by_window

def update_blacklists_by_bam_by_read(blacklists_by_bam_by_read,
tips_to_read_names_file, blacklists_by_window):

  # Extract the window from the file name
  tree_match = read_names_window_regex.search(tips_to_read_names_file)
  if not tree_match:
    print("Error: the file name", tips_to_read_names_file, "does not match the",
    "expected pattern '" + read_names_window_regex_string + "'. Quitting.",
    file=sys.stderr)
    exit(1)
  window = tree_match.groups()[0]

  # Retrieve the blacklist for this window
  if window in blacklists_by_window:
    blacklist = blacklists_by_window[window]
  else:
    print("Error: the file name", tips_to_read_names_file, "corresponds to a",
    "window", window, "that we did not find in the blacklist report. Quitting.",
    file=sys.stderr)
    exit(1)

  tip_to_reads_dict = {}

  with open(tips_to_read_names_file, 'r') as f:
    for line_num_min_1, line in enumerate(f):

      # Remove leading & trailing whitespace
      line = line.strip()

      # Split into fields
      fields = line.split(",")
      if len(fields) < 2:
        print("Error: expected at least 2 fields, found", len(fields),
        "on line", line_num_min_1 + 1, "in", tips_to_read_names_file + \
        ". Quitting.", file=sys.stderr)
        exit(1)
      tip = fields[0]
      read_names = fields[1:]

      if tip in tip_to_reads_dict:
        print("Error: encountered tip", tip, "twice in", \
        tips_to_read_names_file + ". Quitting.", file=sys.stderr)
        exit(1)
      tip_to_reads_dict[tip] = read_names

  all_tips_here = set(tip_to_reads_dict.keys())
  all_tips_before = set(blacklist.keys())
  if all_tips_here != all_tips_before:
    extra_tips = all_tips_here - all_tips_before
    if len(extra_tips) > 0:
      print("Error: these tips, found in", tips_to_read_names_file + \
      ", were not found in the blacklist report:", " ".join(extra_tips),
      file=sys.stderr)
    missing_tips = all_tips_before - all_tips_here
    if len(missing_tips) > 0:
      print("Error: these tips, found in the blacklist report, were not",
      "found in", tips_to_read_names_file + ":", " ".join(missing_tips),
      file=sys.stderr)
    print("Quitting.", file=sys.stderr)
    exit(1)

  # For each tip, record whether that tip was blacklisted or not for all reads
  # associated with the tip. Record results by which bam file that tip came
  # from.
  for tip, read_names in tip_to_reads_dict.items():
    kept = blacklist[tip]
    bam = tip[:tip_regex.search(tip).start()]
    for read_name in read_names:
      blacklists_by_bam_by_read[bam][read_name][kept] += 1




if __name__ == '__main__':

  # Define a function to check files exist, as a type for the argparse.
  def file_type(file_):
    if not os.path.isfile(file_):
      raise argparse.ArgumentTypeError(file_ + \
      ' does not exist or is not a file.')
    return file_

  # Set up the arguments for this script
  parser = argparse.ArgumentParser(description=explanatory_message)
  parser.add_argument('blacklist_report', type=file_type)
  parser.add_argument('tips_to_read_names_csv', type=file_type, nargs="+")
  parser.add_argument('output_file_stem', help='''We will create per-bam output
  files by appending the bam file IDs to this stem.''')
  parser.add_argument('--discarded_reads', action="store_true", help='''With
  this option we will also output the names of the reads that were not kept.''')
  parser.add_argument('--keep_criterion', help='''In general, for each read we
  will have multiple results for whether it is blacklisted or not, due to the
  read appearing in multiple windows. By default we go with the consensus of
  the results: we keep the read if greater than half of the results say so. Use
  this option to specify one of two alternatives: "strict" (keep the read only
  if all results say so) or "permissive" (keep the read if any of the results
  say so).''')
  parser.add_argument('--overwrite', action="store_true", help='''By default, if
  an output file exists already we will exit without overwriting it. With this
  option we will overwrite it.''')
  args = parser.parse_args()

  # Sanity check on the --keep_criterion arg.
  strict = args.keep_criterion == "strict"
  permissive = args.keep_criterion == "permissive"
  if args.keep_criterion != None and (not strict) and (not permissive):
    print("Error: --keep_criterion should be used to specify either 'strict'",
    "or 'permissive'; you specified", args.keep_criterion + ". Quitting.",
    file=sys.stderr)
    exit(1)

  blacklist_report = read_blacklist_report(args.blacklist_report)
  
  # Get the set of blacklist results for each read.
  blacklists_by_bam_by_read = defaultdict(lambda: defaultdict(Counter))
  for tips_to_read_names_file in args.tips_to_read_names_csv:
    update_blacklists_by_bam_by_read(blacklists_by_bam_by_read,
    tips_to_read_names_file, blacklist_report)

  for bam, per_read_blacklists in blacklists_by_bam_by_read.items():

    # Set up the output files for this bam.
    out_file = args.output_file_stem + "_" + bam + ".txt"
    if not args.overwrite and os.path.isfile(out_file):
      print(out_file, "exists already. Quitting to prevent overwriting. (Be",
      "aware of the --overwrite option.)", file=sys.stderr)
      exit(1)
    if args.discarded_reads:
      out_file_discarded = args.output_file_stem + "_" + bam + "_discarded.txt"
      if not args.overwrite and os.path.isfile(out_file_discarded):
        print(out_file_discarded, "exists already. Quitting to prevent",
        "overwriting. (Be aware of the --overwrite option.)", file=sys.stderr)
        exit(1)

    # Decide whether to keep each read based on its blacklist results.
    keep_reads = set([])
    discard_reads = set([])
    for read, keep_and_discard_counts in per_read_blacklists.items():
      keep_count = keep_and_discard_counts[True]
      discard_count = keep_and_discard_counts[False]
      at_least_one_keep = keep_count > 0
      at_least_one_discard = discard_count > 0
      assert at_least_one_keep or at_least_one_discard
      if strict:
        if at_least_one_discard:
          discard_reads.add(read)
        else:
          keep_reads.add(read)
      elif permissive:
        if at_least_one_keep > 0:
          keep_reads.add(read)
        else:
          discard_reads.add(read)
      else:
        total_count = keep_count + discard_count
        if 2 * keep_count > total_count:
          keep_reads.add(read)
        else:
          discard_reads.add(read)

    with open(out_file, "w") as f:
      f.write("\n".join(keep_reads) + "\n")

    if args.discarded_reads:
      with open(out_file_discarded, "w") as f:
        f.write("\n".join(discard_reads) + "\n")
