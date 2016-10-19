list.of.packages <- c("argparse", phangorn)
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
max.reads <- args$maxReadsPerPatient

setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161007_couples_w270_rerun/ptyr1_trees_newick/")
output.dir <- getwd()
script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools/"
tree.file.name <- "ptyr1_InWindow_850_to_1099.tree"
blacklist.file.name <- NULL
root.name <- "REF_CPX_AF460972"
tip.regex <- "^(.*)_read_([0-9]+)_count_([0-9]+)$"
max.reads <- 100

source(file.path(script.dir, "TransmissionUtilityFunctions.R"))
source(file.path(script.dir, "SubtreeMethods.R"))

tree <- read.tree(tree.file.name)
tree <- unroot(tree)
tree <- di2multi(tree, tol = 1E-5)
if(!is.null(root.name)){
  tree <- root(tree, outgroup = which(tree$tip.label==root.name), resolve.root = T)
}

blacklist <- vector()

if(!is.null(blacklist.file.name)){
  if(file.exists(blacklist.file.name)){
    cat("Reading BlackList file",blacklist.file.name,'\n')
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

# This returns the ones to get rid of (i.e. blacklist)

downsample.me <- function(patient, tree, number){
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

excluded <- unlist(lapply(patients.present, downsample.me, tree=tree.1, number = max.reads))

new.blacklist <- c(tree$tip.label[blacklist], new.blacklist)

write.table(new.blacklist, output.file.name, sep=",", row.names=FALSE, col.names=FALSE, quote=F)