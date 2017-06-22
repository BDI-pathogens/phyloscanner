#!/usr/bin/env Rscript

list.of.packages <- c("prodlim","reshape2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos="http://cran.ma.imperial.ac.uk/")

suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
tmp	<- "Summarise topological relationships suggesting direction of transmission across windows. Outputs a .csv file of relationships between patient IDs. Relationship types: 'anc' indicates a topology that suggests the patient in column 1 infected the patient in column 2, 'cher' where there is only one subtree per patient and the MRCAs of both have the same parent in the phylogeny, 'unint' where cher is not present reads from the patients form sibling clades from an unsampled common ancestor patient (and that neither has any sampled ancestors occuring later in the transmission chain than that common ancestor), 'int' one where reads from both patients are intermingled and the direction of transmission cannot be established. (More documentation to follow.)"
arg_parser = ArgumentParser(description=tmp)
arg_parser$add_argument("-m", "--minThreshold", action="store", default=1, type="integer", help="Relationships between two patients will only appear in output if they are within the distance threshold and ajacent to each other least these many windows (default 1). High numbers are useful for drawing figures in e.g. Cytoscape with few enough arrows to be comprehensible. The script is extremely slow if this is set to 0.")
arg_parser$add_argument("-c", "--distanceThreshold", action="store", default=-1, help="Maximum distance threshold on a window for a relationship to be reconstructed between two patients on that window.")
arg_parser$add_argument("-p", "--allowMultiTrans", action="store_true", default=FALSE, help="If absent, directionality is only inferred between pairs of patients where a single clade from one patient is nested in one from the other; this is more conservative")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")
arg_parser$add_argument("idFile", action="store", help="A file containing a list of the IDs of all the patients to calculate and display statistics for.")
arg_parser$add_argument("inputFiles", action="store", help="Either (if -l is present) a list of all input files (output from ClassifyRelationships.R), separated by colons, or (if not) a single string that begins every input file name.")
arg_parser$add_argument("outputFile", action="store", help="A .csv file to write the output to.")
arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")


args <- arg_parser$parse_args()


summary.file             <- args$summaryFile
id.file                  <- args$idFile
output.file              <- args$outputFile
script.dir               <- args$scriptdir  
verbose                  <- args$verbose
csv.fe                   <- args$csvFileExtension
min.threshold            <- as.numeric(args$minThreshold)
dist.threshold           <- as.numeric(args$distanceThreshold)
if(dist.threshold==-1){
  dist.threshold         <- Inf
}
detailed.output          <- args$detailed
if(is.null(min.threshold)){
  split.threshold        <- 1L
}
allow.mt                 <- args$allowMultiTrans
input.file.name          <- args$inputFiles

source(file.path(script.dir, "TreeUtilityFunctions.R"))
source(file.path(script.dir, "GeneralFunctions.R"))
source(file.path(script.dir, "CollapsedTreeMethods2.R"))

suppressMessages(library(prodlim, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(reshape2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(gdata, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ggplot2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(data.table, quietly=TRUE, warn.conflicts=FALSE))


input.files <- list.files.mod(dirname(input.file.name), pattern=paste(basename(input.file.name)), full.names=TRUE)

if(length(input.files)==0){
  stop("No input files found.")
}


# Read in the IDs. Remove duplicates. Shuffle their order if desired.

if(verbose) cat("Reading patient IDs...\n")

patient.ids	<- unique(scan(id.file, what="", sep="\n", quiet=TRUE))
if(!length(patient.ids)){
  stop(paste("No IDs found in ", id.file, ". Quitting.\n", sep=""))  
}

all.tree.info <- list()

for(file in input.files){
  tree.info <- list()
  tree.info$classification.file.name <- file
  
  suffix <- get.suffix(file, input.file.name, csv.fe)
  tree.info$suffix <- suffix
  
  all.tree.info[[file]] <- tree.info
}


results <- summarise.classifications(all.tree.info, patient.ids, min.threshold, dist.threshold, allow.mt, csv.fe, verbose)

if (verbose) cat('Writing summary to file',output.file,'\n')
write.csv(results, file=output.file, row.names=FALSE, quote=FALSE)

