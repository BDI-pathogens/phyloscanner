#!/usr/bin/env Rscript


suppressMessages(library(prodlim, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(phyloscannerR, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(reshape2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(gdata, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ggplot2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(data.table, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(kimisc, quietly=TRUE, warn.conflicts=FALSE))

suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
tmp	<- "Summarise topological relationships suggesting direction of transmission across windows. Outputs a .csv file of relationships between patient IDs."
arg_parser = ArgumentParser(description=tmp)
arg_parser$add_argument("-m", "--minThreshold", action="store", default=0, type="double", help="Relationships between two patients will only appear in output if they are within the distance threshold and ajacent to each other in more than this proportion of trees many tree (default 0). High numbers are useful for drawing figures in e.g. Cytoscape with few enough arrows to be comprehensible.")
arg_parser$add_argument("-c", "--distanceThreshold", action="store", default=-1, help="Maximum distance threshold on a window for a relationship to be reconstructed between two patients on that window. If absent, no such threshold will be applied.")
arg_parser$add_argument("-p", "--allowMultiTrans", action="store_true", default=FALSE, help="If absent, directionality is only inferred between pairs of patients where a single clade from one patient is nested in one from the other; this is more conservative")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")
arg_parser$add_argument("inputFiles", action="store", help="Either a list of all input files (output from classify_relationships.R), separated by colons, or a single string that begins every input file name.")
arg_parser$add_argument("outputFile", action="store", help="A .csv file to write the output to.")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")

args <- arg_parser$parse_args()

output.file              <- args$outputFile
verbose                  <- args$verbose
csv.fe                   <- args$csvFileExtension
min.threshold            <- as.numeric(args$minThreshold)
dist.threshold           <- as.numeric(args$distanceThreshold)
if(dist.threshold==-1){
  dist.threshold         <- Inf
}
if(is.null(min.threshold)){
  split.threshold        <- 1L
}
allow.mt                 <- args$allowMultiTrans
input.file.name          <- args$inputFiles

input.files <- list.files.mod(dirname(input.file.name), pattern=paste(basename(input.file.name)), full.names=TRUE)

if(length(input.files)==0){
  stop("No input files found.")
}

all.tree.info <- list()

for(file in input.files){
  tree.info <- list()
  tree.info$classification.file.name <- file
  
  suffix <- get.suffix(file, input.file.name, csv.fe)
  tree.info$suffix <- suffix
  
  all.tree.info[[file]] <- tree.info
}

results <- summarise.classifications(all.tree.info, min.threshold*length(all.tree.info), dist.threshold, allow.mt, F, verbose)

if (verbose) cat('Writing summary to file',output.file,'\n')
write.csv(results, file=output.file, row.names=FALSE, quote=FALSE)

