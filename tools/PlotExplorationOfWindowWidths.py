#!/usr/bin/python2
# TODO: change shebang
from __future__ import print_function

## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''TODO'''

import os
import sys
import argparse
import pysam
import collections
import re
import numpy as np
import matplotlib.pyplot as plt

def read_window_width_file(file_):
  '''Reads in a file of the format produced by phyloscanner_make_trees.py's
  --explore-window-widths option'''

  title_regex_string = "Number of unique reads per-bam and per-window with " + \
  "window width = (\d+):"
  title_regex = re.compile(title_regex_string)

  counts_by_window_width = collections.defaultdict(list)
  current_window_width = None

  with open(file_) as f:
    for line_num_min_1, line in enumerate(f):
      line = line.strip()

      # Blank lines mark the end of data for one window width
      if line == "":
        current_window_width = None
        continue
        
      # Skip the header of each csv block
      if line.startswith("Window start,"):
        continue

      # For error messages
      error_suffix = "on line " + str(line_num_min_1 + 1) + " in " + file_ + \
      ". Quitting."

      # Extract the window width from relevant lines
      if title_regex.match(line):
        try:
          current_window_width = int(title_regex.match(line).groups()[0])
        except:
          print("Internal error extracting the window width from '" + line + \
          "' by comparison with regex '" + title_regex_string + "',",
          error_suffix, file=sys.stderr)
          exit(1)
        continue

      # If we're at a data line without having met a window width line yet,
      # something is wrong.
      if current_window_width == None:
        print("Error: encountered line which is neither blank nor the start",
        "of a new block of data for a particular window width",
        error_suffix, file=sys.stderr)
        exit(1)

      # Data lines = comma-separated ints
      fields = line.split(",")
      try:
        window_and_counts = [int(field) for field in fields]
      except ValueError:
        print("Could not understand the fields as comma-separated integers",
        error_suffix, file=sys.stderr)
        exit(1)

      # There should be at least one count. Record them.
      counts = window_and_counts[1:]
      if len(counts) == 0:
        print("No count values found", error_suffix, file=sys.stderr)
        exit(1)
      counts_by_window_width[current_window_width] += counts
      
  # Check we found at least one window width
  if len(counts_by_window_width) == 0:
    print("No count values found in", file_ + ". Quitting.", file=sys.stderr)
    exit(1)
       
  return counts_by_window_width
    
def calculate_percentiles(counts_by_window_width, percentiles):
  '''For each window width we calculate the desired percentiles for that set of
  counts. Returns a matrix where each row is one percentile and each column is
  one window width.'''

  matrix = np.empty((len(percentiles), len(counts_by_window_width)))
  windows = []
  for i, (window, counts) in enumerate(sorted(counts_by_window_width.items(),
  key=lambda x:x[0])):
    windows.append(window)
    matrix[:, i] = np.percentile(counts, percentiles)
  return windows, matrix
  
def plot_percentiles_by_window_width(windows, matrix, args):
  '''TODO'''

  ax = plt.figure().add_subplot(111)

  for row_number in range(matrix.shape[0]):
    percentile = args.percentiles[row_number]
    this_percentile_for_all_window_widths = matrix[row_number, :]
    plt.plot(windows, this_percentile_for_all_window_widths, label=percentile)

  if args.x_min_max:
    ax.set_xlim(xmin=Xmin, xmax=Xmax)
  if args.y_min_max:
    ax.set_ylim(ymin=Ymin, ymax=Ymax)

  handles, labels = ax.get_legend_handles_labels()
  labels, handles = \
  zip(*sorted(zip(labels, handles), key=lambda t: t[0], reverse=True))

  plt.xlabel('window width', fontsize=args.axis_font_size)
  plt.ylabel('number of unique reads', fontsize=args.axis_font_size)
  title = \
  'Characterising the distribution, over all bam files and all genomic\n' + \
  'windows, of the number of unique reads as a function of window width.\n'
  plt.title(title, fontsize=args.title_font_size)
  ax.tick_params(axis='both', which='major', labelsize=args.axis_font_size)
  plt.legend(handles, labels, title="percentile", loc=args.legend_location,
  fontsize=args.legend_font_size)
  plt.tight_layout()
  plt.savefig(args.OutputPDF)


if __name__ == '__main__':

  # Define a function to check files exist, as a type for the argparse.
  def File(MyFile):
    if not os.path.isfile(MyFile):
      raise argparse.ArgumentTypeError(MyFile + \
      ' does not exist or is not a file.')
    return MyFile

  # Define a comma-separated integers object, as a type for the argparse.
  def CommaSeparatedPercentages(MyCommaSeparatedPercentages):
    try:
      ListOfPercentages = [float(value) for value in
      MyCommaSeparatedPercentages.split(',')]
    except:
      raise argparse.ArgumentTypeError('Unable to understand ' +\
      MyCommaSeparatedPercentages + ' as comma-separated floats.')
    else:
      try:
        assert all(0 <= value <= 100 for value in MyCommaSeparatedPercentages)
      except AssertionError:
        raise argparse.ArgumentTypeError('Percentiles must not be less than ' +\
        'zero or greater than 100.')
      return list(set(ListOfPercentages))

  # Set up the arguments for this script
  parser = argparse.ArgumentParser(description=ExplanatoryMessage)
  parser.add_argument('DataFile', type=File, help='''TODO''')
  parser.add_argument('OutputPDF', help="")
  parser.add_argument('-P', '--percentiles', type=CommaSeparatedPercentages,
  help='''The percentiles to plot, as a comma-separated list. Default:
  10,20,30,40,50,60,70,80,90''', default=[10,20,30,40,50,60,70,80,90])
  parser.add_argument('-AS', '--axis-font-size', type=int,
  help='For the plot. The default is 15.', default=15)
  parser.add_argument('-TS', '--title-font-size', type=int,
  help='For the plot. The default is 15.', default=15)
  parser.add_argument('-LS', '--legend-font-size', type=int,
  help='For the plot. The default is 8.', default=8)
  parser.add_argument('-LL', '--legend-location', 
  help='''For the plot. The default is 'upper right'. The other options are:
  'best', 'upper right', 'upper left', 'lower right', 'right', 'center left',
  'center right', 'lower center',' upper center', 'center' ''',
  default='upper right')
  parser.add_argument('-XM', '--x-min-max', help='The minimum and maximum for '\
  'the x axis in the plot, specified together as a comma-separated pair of '\
  'numbers.')
  parser.add_argument('-YM', '--y-min-max', help='The minimum and maximum for '\
  'the y axis in the plot, specified together as a comma-separated pair of '\
  'numbers.')
  args = parser.parse_args()

  counts_by_window_width = read_window_width_file(args.DataFile)
  
  windows, matrix = calculate_percentiles(counts_by_window_width,
  args.percentiles)
  
  plot_percentiles_by_window_width(windows, matrix, args)
