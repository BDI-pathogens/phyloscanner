#!/usr/bin/env Rscript

list.of.packages <- c("argparse", "ape", "kimisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Please run package_install.R to continue\n")
  quit(save="no", status=1)
}

suppressMessages(require(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ape, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(kimisc, quietly=TRUE, warn.conflicts=FALSE))

# 	Define arguments

arg_parser		     <- ArgumentParser(description="Write a file of host IDs from all tree tips in the input")

arg_parser$add_argument("tree.file.root", action="store", type="character", help="Start of tree file names, or single tree file name")
arg_parser$add_argument("output.file.name", action="store", type="character", help="Output file name")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what the script is doing.")
arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. Three capture groups: host ID, read ID, and read count; if the latter two groups are missing then read information will not be used. If absent, input will be assumed to be from the phyloscanner pipeline, and the host ID will be the BAM file name.")
arg_parser$add_argument("-tfe", "--treeFileExtension", action="store", default="tree", help="The file extension for tree files (default .tree).")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")
arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the /tools directory.")

# Parse arguments

args                <- arg_parser$parse_args()

if(!is.null(args$scriptDir)){
  script.dir          <- args$scriptDir
} else {
  script.dir          <- dirname(thisfile())
  if(!dir.exists(script.dir)){
    stop("Cannot detect the location of the /phyloscanner/tools directory. Please specify it at the command line with -D.")
  }
}

tree.file.root 	    <- args$tree.file.root
output.file.name    <- args$output.file.name
tip.regex           <- args$tipRegex
verbose             <- args$verbose
tree.fe             <- args$treeFileExtension
csv.fe              <- args$csvFileExtension

source(file.path(script.dir, "tree_utility_functions.R"))
source(file.path(script.dir, "general_functions.R"))

if(file.exists(tree.file.root)){
  tree.file.names <- tree.file.root
} else {
  tree.file.names <- list.files.mod(dirname(tree.file.root), pattern=paste(basename(tree.file.root),'.*',tree.fe,'$',sep=''), full.names=TRUE)
}

host.names <- vector()

for(tree.file.name in tree.file.names){
  tree <- read.tree(tree.file.name)
  new.names <- sapply(tree$tip.label, function(x) host.from.label(x, tip.regex))
  new.names <- na.omit(new.names)
  host.names <- unique(c(host.names, new.names))
}

host.names <- host.names[order(host.names)]

if(verbose) cat("Identified ",length(host.names), " hosts from these trees.")

write.table(host.names, output.file.name, sep=",", row.names=FALSE, col.names=FALSE, quote=F)