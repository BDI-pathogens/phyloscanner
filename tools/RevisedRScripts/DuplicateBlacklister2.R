list.of.packages <- "argparse"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))

arg_parser = ArgumentParser(description="Identify phylogeny tips for blacklisting as representing suspected contaminants, based on being exact duplicates of more numerous reads from other patients")

arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                        help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
arg_parser$add_argument("-b", "--blacklist", action="store", help="A blacklist to be applied before this script is run.")
arg_parser$add_argument("rawThreshold", action="store", type="double", help="Raw threshold; tips with read counts less than this that are identical to a tip from another patient will be blacklisted, regardless of the count of the other read.")
arg_parser$add_argument("ratioThreshold", action="store", type="double", help="Ratio threshold; tips will be blacklisted if the ratio of their tip count to that of another, identical tip from another patient is less than this value.")
arg_parser$add_argument("inputFileName", action="store", help="A file (comma-separated) outlining groups of tips that have identical sequences, each forming a single line.")
arg_parser$add_argument("outputFileName", action="store", help="The file to write the output to, a list of tips to be blacklisted.")
arg_parser$add_argument("-D", "--scriptDir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")

# Parse arguments

args                 <- arg_parser$parse_args()

script.dir           <- args$scriptDir
verbose              <- args$verbose
tree.fe              <- args$treeFileExtension
csv.fe               <- args$csvFileExtension
raw.threshold        <- args$rawThreshold
ratio.threshold      <- args$ratioThreshold
tip.regex            <- args$tipRegex
input.file.name      <- args$inputFileName
output.file.name     <- args$outputFileName
blacklist.file.name  <- args$blacklist

# Load necessary functions

source(file.path(script.dir, "TreeUtilityFunctions.R"))
source(file.path(script.dir, "GeneralFunctions.R"))
source(file.path(script.dir, "RevisedRScripts/BlacklistFunctions.R"))

# Decide if we're in single or batch mode

if(file.exists(input.file.name)){
  
  input.file.names <- input.file.name
  output.file.names <- output.file.name
  if(!is.null(blacklist.file.name)){
    blacklist.file.names <- blacklist.file.name
  }
  
} else {

  file.regex <- paste0(basename(input.file.name), '.*\\.',csv.fe,'$')
  input.file.names	<- list.files.mod(dirname(input.file.name), pattern=file.regex, full.names=TRUE)
  if (length(input.file.names) == 0) {
    stop(paste0("Error: failed to find read duplication data. Searched for file names matching the regex ", file.regex, " in the directory ", dirname(input.file.name), ". Quitting.\n"))
  }
  suffixes <- substr(input.file.names, nchar(input.file.name) + 1, nchar(input.file.names)-nchar(csv.fe))
  output.file.names <- paste0(output.file.name, suffixes, csv.fe)
  if(!is.null(blacklist.file.name)){
    blacklist.file.names <- paste0(blacklist.file.name, suffixes, csv.fe)
  }
  
}



for(file.no in 1:length(input.file.names)){
  input.file.name <- input.file.names[file.no]
  blacklist.file.name <- blacklist.file.names[file.no]
  if(!file.exists(blacklist.file.name)){
    warning(paste0("Specified blacklist file ",blacklist.file.name," does not exist; will be ignored.\n"))
    already.blacklisted <- vector()
  } else {
    already.blacklisted <- readLines(blacklist.file.name, warn=F)
  }
  
  if(verbose) cat("DuplicateBlacklister.R run on: ", input.file.name, "\n", sep="")
  
  # Read the lines of each input file and process them as a list
  
  entries <- strsplit(readLines(input.file.name, warn=F),",")
  
  blacklisted <- blacklist.exact.duplicates(entries, raw.threshold, ratio.threshold, tip.regex, verbose)
  
  if (verbose) cat(length(blacklisted)," newly blacklisted tips.\n",sep="")

  blacklisted <- unique(c(blacklisted, already.blacklisted))
  
  write.table(blacklisted, output.file.names[file.no], sep=",", row.names=FALSE, col.names=FALSE, quote=F)
}
