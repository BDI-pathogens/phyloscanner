suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(data.table, quietly=TRUE, warn.conflicts=FALSE))


# 	Define arguments

arg_parser		<- ArgumentParser(description="Calculate normalising constants for each tree file from reference table.")
arg_parser$add_argument("tree.file.root", action="store", type="character", help="Start of tree file names.")
arg_parser$add_argument("norm.file.name", action="store", type="character", help="File name of reference table.")
arg_parser$add_argument("output.file.name", action="store", type="character", help="Output file name of csv file that has two columns, the base tree file name and the corresponding normalising constant.")
arg_parser$add_argument("norm.var", action="store", type="character", help="Column name in reference table that is to be used as normalising constant.")
arg_parser$add_argument("-s", "--standardise", action="store_true", default=FALSE, help="If true, the normalising constants are standardized so that the average on gag+pol equals 1. This way the normalised branch lengths are interpretable as typical distances on gag+pol")
arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what the script is doing.")
arg_parser$add_argument("-tfe", "--treeFileExtension", action="store", default="tree", help="The file extension for tree files (default .tree).")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")

# Parse arguments

args              <- arg_parser$parse_args()
script.dir        <- args$scriptdir
tree.file.root 	  <- args$tree.file.root
norm.file.name 	  <- args$norm.file.name
output.file.name  <- args$output.file.name
norm.var 		      <- args$norm.var
norm.standardise  <- args$standardise
verbose           <- args$verbose
tree.fe           <- args$treeFileExtension
csv.fe            <- args$csvFileExtension

# Load necessary functions

source(file.path(script.dir, "GeneralFunctions.R"))
source(file.path(script.dir, "RevisedRScripts/NormalisationFunctions.R"))

#	Load reference table, define normalising constant

if (verbose) cat('Loading normalising constants reference file ', norm.file.name, "\n", sep="")

if(grepl(paste0(csv.fe, "$"), norm.file.name)){
	norm.table	<- as.data.table(read.csv(norm.file.name, stringsAsFactors=FALSE))
} else if(grepl('rda$', norm.file.name)){
	tmp			<- load(norm.file.name)
	if(length(tmp)!=1)	stop("Expected one R data.table in file ",norm.file.name)
	eval(parse(text=paste("norm.table<- ",tmp,sep='')))
} else {
  stop(paste0("Unknown input file format.\n"))
}

if(any(!(c('W_FROM','W_TO',norm.var) %in% colnames(norm.table)))){
  stop(paste0("One or more of W_FROM, W_TO and ",norm.var," columns not present in file ",norm.file.name))
}

setnames(norm.table, norm.var, 'NORM_CONST')

#	Standardize to mean of 1 on gag+pol ( prot + first part of RT in total 1300bp )

if(norm.standardise){
  norm.table[, W_MID:= (W_FROM+W_TO)/2]
	#790 - 3385
	if (verbose) cat('Standardising normalising constants to 1 on the gag+pol (prot + first part of RT in total 1300bp pol) region\n')
  
	tmp		<- subset(norm.table, W_MID>=790L & W_MID<=3385L)
	
	if(nrow(tmp)<=0){
	  stop(paste0("No positions from gag+pol present in file ",norm.file.name,"; unable to standardise"))
	}
	
	tmp		<- tmp[, mean(NORM_CONST)]
	if(!is.finite(tmp)){
	  stop(paste0("Standardising constant is not finite"))
	}

	set(norm.table, NULL, 'NORM_CONST', norm.table[, NORM_CONST/tmp])
}

norm.table	         <- subset(norm.table, select=c("W_FROM", "W_TO", "NORM_CONST"))

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
