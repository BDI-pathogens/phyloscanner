#!/usr/bin/env Rscript

options("warn"=1)

list.of.packages <- c("argparse", "kimisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Please run package_install.R to continue\n")
  quit(save="no", status=1)
}

suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(kimisc, quietly=TRUE, warn.conflicts=FALSE))

arg_parser = ArgumentParser(description="Identify phylogeny tips for blacklisting as representing suspected contaminants, based on being exact duplicates of more numerous reads from other hosts")

arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                        help="Regular expression identifying tips from the dataset. Three capture groups: host ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the host ID will be the BAM file name.")
arg_parser$add_argument("-b", "--blacklist", action="store", help="A blacklist to be applied before this script is run.")
arg_parser$add_argument("rawThreshold", action="store", type="double", help="Raw threshold; tips with read counts less than this that are identical to a tip from another host will be blacklisted, regardless of the count of the other read.")
arg_parser$add_argument("ratioThreshold", action="store", type="double", help="Ratio threshold; tips will be blacklisted if the ratio of their tip count to that of another, identical tip from another host is less than this value.")
arg_parser$add_argument("inputFileName", action="store", help="A file (comma-separated) outlining groups of tips that have identical sequences, each forming a single line.")
arg_parser$add_argument("outputFileName", action="store", help="The file to write the output to, a list of tips to be blacklisted.")
arg_parser$add_argument("-D", "--scriptDir", action="store", help="Full path of the /tools directory.")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what the script doing.")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")

# Parse arguments

args                 <- arg_parser$parse_args()

if(!is.null(args$scriptDir)){
  script.dir          <- args$scriptDir
} else {
  script.dir          <- dirname(thisfile())
  if(!dir.exists(script.dir)){
    stop("Cannot detect the location of the /phyloscanner/tools directory. Please specify it at the command line with -D.")
  }
}

verbose              <- args$verbose
tree.fe              <- args$treeFileExtension
csv.fe               <- args$csvFileExtension
raw.threshold        <- args$rawThreshold
ratio.threshold      <- args$ratioThreshold
tip.regex            <- args$tipRegex
input.string         <- args$inputFileName
output.string        <- args$outputFileName
blacklist.input      <- args$blacklist

# Load necessary functions

source(file.path(script.dir, "tree_utility_functions.R"))
source(file.path(script.dir, "general_functions.R"))
source(file.path(script.dir, "blacklist_functions.R"))

# Decide if we're in single or batch mode

all.tree.info <- list()

if(file.exists(input.string)){
  
  tree.info <- list()
  
  tree.info$input.file.name <- input.string
  tree.info$output.file.name <- output.string

  if(!is.null(blacklist.input)){
    tree.info$blacklist.file.name <- blacklist.input
  }
  
  all.tree.info[[input.string]] <- tree.info
  
} else {

  file.regex <- paste0(basename(input.string), '.*\\.',csv.fe,'$')
  input.file.names	<- list.files.mod(dirname(input.string), pattern=file.regex, full.names=TRUE)
  if (length(input.file.names) == 0) {
    stop(paste0("Error: failed to find read duplication data. Searched for file names matching the regex ", file.regex, " in the directory ", dirname(input.string), ". Quitting.\n"))
  }
  suffixes <- substr(input.file.names, nchar(input.string) + 1, nchar(input.file.names)-nchar(csv.fe)-1)
  
  for(suffix.no in 1:length(suffixes)){
    tree.info <- list()
    
    tree.info$suffix <- suffixes[suffix.no]
    tree.info$input.file.name <- input.file.names[suffix.no]
    
    tree.info$output.file.name <- paste0(output.string, "_", suffixes[suffix.no], ".", csv.fe)
    
    if(!is.null(blacklist.input)){
      tree.info$blacklist.file.name <- paste0(blacklist.input, suffixes[suffix.no], ".", csv.fe)
    }
    
    all.tree.info[[suffixes[suffix.no]]] <- tree.info
  }
  
  
}

for(tree.info in all.tree.info){
  input.file.name <- tree.info$input.file.name

  if(!is.null(tree.info$blacklist.file.name)){
    if(!file.exists( tree.info$blacklist.file.name)){
      warning(paste0("Specified blacklist file ", tree.info$blacklist.file.name," does not exist; will be ignored.\n"))
      already.blacklisted <- vector()
    } else {
      already.blacklisted <- readLines( tree.info$blacklist.file.name, warn=F)
    } 
  } else {
    already.blacklisted <- vector()
  }
  
  if(verbose) cat("DuplicateBlacklister.R run on: ", input.file.name, "\n", sep="")
  
  # Read the lines of each input file and process them as a list
  
  tree.info$duplicate.tips <- strsplit(readLines(input.file.name, warn=F),",")
  
  blacklisted <- blacklist.exact.duplicates(tree.info, raw.threshold, ratio.threshold, tip.regex, verbose)
  
  if (verbose) cat(length(blacklisted)," newly blacklisted tips.\n",sep="")

  blacklisted <- unique(c(blacklisted, already.blacklisted))
  
  write.table(blacklisted, tree.info$output.file.name, sep=",", row.names=FALSE, col.names=FALSE, quote=F)
}
