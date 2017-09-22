#!/usr/bin/env Rscript

list.of.packages <- c("argparse", "data.table", "kimisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Missing dependencies; replacing the path below appropriately, run\n[path to your phyloscanner code]/tools/package_install.R\nthen try again.\n")
  quit(save="no", status=1)
}

suppressMessages(require(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(data.table, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(kimisc, quietly=TRUE, warn.conflicts=FALSE))

# 	Define arguments

arg_parser		     <- ArgumentParser(description="Calculate normalising constants for each tree file from reference table.")

arg_parser$add_argument("tree.file.root", action="store", type="character", help="Start of tree file names.")
arg_parser$add_argument("norm.file.name", action="store", type="character", help="File name of reference table.")
arg_parser$add_argument("output.file.name", action="store", type="character", help="Output file name of csv file that has two columns, the base tree file name and the corresponding normalising constant.")
arg_parser$add_argument("-s", "--standardiseGagPol", action="store_true", default=FALSE, help="HIV-specific; if true, the normalising constants are standardized so that the average on gag+pol equals 1, rather than that the average on the whole genome equals 1. The normalised branch lengths are interpretable as typical distances on gag+pol")
arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the /tools directory.")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what the script is doing.")
arg_parser$add_argument("-tfe", "--treeFileExtension", action="store", default="tree", help="The file extension for tree files (default .tree).")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")

# Parse arguments

args              <- arg_parser$parse_args()

if(!is.null(args$scriptDir)){
  script.dir          <- args$scriptDir
} else {
  script.dir          <- dirname(thisfile())
  if(!dir.exists(script.dir)){
    stop("Cannot detect the location of the /phyloscanner/tools directory. Please specify it at the command line with -D.")
  }
}

tree.file.root 	    <- args$tree.file.root
norm.ref.file.name 	<- args$norm.file.name
output.file.name    <- args$output.file.name
norm.var 		        <- args$norm.var
norm.standardise    <- args$standardiseGagPol
verbose             <- args$verbose
tree.fe             <- args$treeFileExtension
csv.fe              <- args$csvFileExtension

# Load necessary functions

source(file.path(script.dir, "general_functions.R"))
source(file.path(script.dir, "normalisation_functions.R"))

#	Load reference table, define normalising constant

if (verbose) cat('Loading normalising constants reference file ', norm.ref.file.name, "\n", sep="")

norm.table	<- as.data.table(read.csv(norm.ref.file.name, stringsAsFactors=FALSE))

if(ncol(norm.table)!=2){
  stop(paste0(norm.ref.file.name," is not formatted as expected for a normalisation lookup file; expecting two columns.\n"))
} else {
  setnames(norm.table, 1, 'POSITION')
  setnames(norm.table, 2, 'NORM_CONST')
  
  if(norm.standardise){
    #	Standardize to mean of 1 on gag+pol ( prot + first part of RT in total 1300bp )
    
    #790 - 3385
    if (verbose) cat('Standardising normalising constants to 1 on the gag+pol (prot + first part of RT in total 1300bp pol) region\n')
    
    tmp		<- subset(norm.table, POSITION>=790L & POSITION<=3385L)
    
    if(nrow(tmp)<=0){
      stop(paste0("No positions from gag+pol present in file ",norm.ref.file.name,"; unable to standardise"))
    }
    
    tmp		<- tmp[, mean(NORM_CONST)]
    
    if(!is.finite(tmp)){
      stop(paste0("Standardising constant is not finite"))
    }
    
    set(norm.table, NULL, 'NORM_CONST', norm.table[, NORM_CONST/tmp])
  } else {
    if (verbose) cat('Standardising normalising constants to 1 on the whole genome\n')
    
    tmp		<- norm.table[, mean(NORM_CONST)]
    if(!is.finite(tmp)){
      stop(paste0("Standardising constant is not finite"))
    }
    
    set(norm.table, NULL, 'NORM_CONST', norm.table[, NORM_CONST/tmp])
  }
}

#	guess window coordinates of tree files from file name
if (verbose) cat('Reading tree files starting with ', tree.file.root, "\n")

wdf		<- list.files.mod(dirname(tree.file.root), pattern=paste('^',basename(tree.file.root),'.*\\.tree$',sep=''), full.names=FALSE)

if(length(wdf)==0){
  stop("ERROR: No tree files found.\n")
}

tree.info <- lapply(wdf, function(x){
  suffix <- gsub(paste0('\\.',tree.fe),'',substring(x, nchar(basename(tree.file.root)) + 1L))
  coords <- get.window.coords(suffix)
  
  list(name =  gsub(paste0('\\.',tree.fe),'',substring(x, nchar(basename(tree.file.root)) + 1L)),
       tree.file.name = x,
       window.coords = coords 
       )
})

nc.lookup <- function(name) {
  sapply(name, function(x) {
    lookup.normalisation.for.tree(tree.info[[x]], norm.table, lookup.column = "NORM_CONST")})
}

wdf		<- data.table(FILE_NAME=sort(wdf))

wdf[, NORM.CONST := nc.lookup(FILE_NAME) ]


#	write to file

if(verbose) cat('Writing normalising constants to file', output.file.name,'\n')
write.table(wdf, output.file.name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
