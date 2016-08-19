# Quick script to gather all patient IDs matching a given regexp from one or more tree files

require(ape)
require(optparse)

source("/Users/twoseventwo/Documents/phylotypes/TransmissionUtilityFunctions.R")

option_list = list(
  make_option("--treeFileRoot", type="character", default=NULL, 
              help="A string which every tree file begins with."),
  make_option("--treeFileList", type="character", default=NULL, 
              help="A list of tree files (single string, with file names separated by with colons). Ignored if
              --treeFileRoot is present"),
  make_option(c("-x", "--tipRegex"), type="character", default=NULL, 
              help="Regular expression identifying tips from the dataset. Three groups: patient ID,
              read ID, and read count."),
  make_option(c("-o", "--outputFileName"), type="character", default=NULL, help="File name for output")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

tree.file.root <- opt$treeFileRoot

if(!is.null(tree.file.root)){
  tree.files <- sort(list.files(getwd(), pattern=paste(tree.file.root,".*\\.tree",sep="")))
} else {
  tree.files <- unlist(strsplit(opt$treeFileList, ":"))
}

tip.regex <- opt$tipRegex
out.file.name <- opt$outputFileName

names <- vector()

for(tree.file in tree.files){
  tree <- read.tree(tree.file)
  tip.names <- tree$tip.label
  names <- c(names, sapply(tip.names, function(x) patient.from.label(x, tip.regex)))
}

names <- names[which(!is.na(names))]
names <- unique(names)
names <- sort(names)

write.table(names, out.file.name, quote=F, row.names = F, col.names = F)
