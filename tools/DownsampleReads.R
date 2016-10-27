list.of.packages <- c("argparse", "phangorn")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

library(argparse)
library(phangorn)

arg_parser = ArgumentParser(description="Examine all patients from the input tree that match a given regexp and blacklist tips with low read counts that are distant from the main clades of tips (and hence likely contaminants).")

arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                        help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
arg_parser$add_argument("-b", "--blacklist", action="store", help="A blacklist to be applied before this script is run.")
arg_parser$add_argument("maxReadsPerPatient", type="double", action="store", help="The upper limit for the number of reads to be included from each patient")
arg_parser$add_argument("inputFile", metavar="inputTreeFileName", help="Tree file name. Alternatively, a base name that identifies a group of tree file names can be specified. Tree files are assumed to end in .tree.")  
arg_parser$add_argument("outputFile", metavar="outputFileName", help="The file to write the output to, a list of tips to be blacklisted.")  
arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")


args <- arg_parser$parse_args()
script.dir <- args$scriptdir
tip.regex <- args$tipRegex
input.file.name <- args$inputFile
output.file.name <- args$outputFile
blacklist.file.name <- args$blacklist
root.name <- args$outgroupName
max.reads <- args$maxReadsPerPatient

source(file.path(script.dir, "TransmissionUtilityFunctions.R"))
source(file.path(script.dir, "SubtreeMethods.R"))

#
#	additional function definitions
#

# This returns the ones to get rid of (i.e. blacklist) per patient for a given tree
downsample.patient <- function(patient, tree, number, tip.regex, patient.ids){
	tips.to.keep <- which(patient.ids==patient)
	labels.to.keep <- tree$tip.label[tips.to.keep]
	
	tip.counts <- as.numeric(sapply(labels.to.keep, function(tip) read.count.from.label(tip, tip.regex))) 
	
	total.reads <- sum(tip.counts)
	
	if(total.reads <= number){
		return(vector())
	}
	
	tip.props <- tip.counts/total.reads
	
	sample <- rmultinom(1, size = number, prob = tip.props)
	
	sampled.names <- labels.to.keep[which(sample > 0)]
	
	return(setdiff(labels.to.keep, sampled.names))
}
# This returns the ones to get rid of (i.e. blacklist) for a given tree
downsample.tree<- function(input.file.name, blacklist.file.name, output.file.name, root.name, tip.regex, max.reads)
{
	tree <- read.tree(input.file.name)
	tree <- unroot(tree)
	tree <- di2multi(tree, tol = 1E-5)
	if(!is.null(root.name)){
		tree <- root(tree, outgroup = which(tree$tip.label==root.name), resolve.root = T)
	}
	
	blacklist <- vector()
	
	if(!is.null(blacklist.file.name)){
		if(file.exists(blacklist.file.name)){
			cat("Reading blacklist file",blacklist.file.name,'\n')
			blacklisted.tips <- read.table(blacklist.file.name, sep=",", header=F, stringsAsFactors = F, col.names="read")
			if(nrow(blacklisted.tips)>0){
				blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
			}
		} else {
			warning(paste("File ",blacklist.file.name," does not exist; skipping.",paste=""))
		}	
	}
	
	tree.1 <- drop.tip(tree, blacklist)
	
	patients.present <- sapply(tree.1$tip.label, function(tip) patient.from.label(tip, tip.regex))
	patients.present <- patients.present[!is.na(patients.present)]
	patients.present <- unique(patients.present)
	
	patient.ids <- sapply(tree.1$tip.label, function(x) patient.from.label(x, tip.regex))
	
	
	
	excluded <- unlist(lapply(patients.present, downsample.patient, tree=tree.1, number = max.reads, tip.regex=tip.regex, patient.ids=patient.ids))
	
	new.blacklist <- c(tree$tip.label[blacklist], excluded)
	write.table(new.blacklist, output.file.name, sep=",", row.names=FALSE, col.names=FALSE, quote=F)
}
#
#	run script 
#

#	option 1: input.file.name specifies a single input file
if(file.exists(input.file.name))
{
	downsample.tree(input.file.name, blacklist.file.name, output.file.name, root.name, tip.regex, max.reads)
}
#	option 2: input.file.name specifies a regular expression of several input files
if(!file.exists(input.file.name))
{
	input.file.names		<- sort(list.files(dirname(input.file.name), pattern=paste(basename(input.file.name),'.*\\.tree$',sep=''), full.names=TRUE))
	blacklist.file.names	<- NULL
	if(!is.null(blacklist.file.name))
		blacklist.file.names<- sort(list.files(dirname(blacklist.file.name), pattern=paste(basename(blacklist.file.name),'.*\\.csv$',sep=''), full.names=TRUE))	
	stopifnot(length(input.file.names)==length(blacklist.file.names))
	for(tree.i in seq_along(input.file.names))
	{
		input.file.name		<- input.file.names[tree.i]
		blacklist.file.name	<- blacklist.file.names[tree.i]		
		tmp					<- paste(output.file.name, gsub('tree$','csv',regmatches(input.file.name,regexpr('InWindow.*',input.file.name))),sep='')		
		downsample.tree(input.file.name, blacklist.file.name, tmp, root.name, tip.regex, max.reads)		
	}
}