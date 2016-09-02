command.line <- T
list.of.packages <- c("phangorn", "argparse", "phytools", "ggplot2", "gdata", "mvtnorm", "expm")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

if(!("ggtree" %in% installed.packages()[,"Package"])){
  source("https://bioconductor.org/biocLite.R")
  biocLite("ggtree")
}

if(command.line){
  require(argparse)
  
  arg_parser = ArgumentParser(description="Split the tips from each patient up, according to one of several schemes, such that the tree can subsequently be partitioned into connected subtrees, one per split.")
  
  arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                          help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
  arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
  arg_parser$add_argument("-z", "--zeroLengthTipsCount", action="store_true", default=FALSE, help="If present, a zero length terminal branch associates the patient at the tip with its parent node, interrupting any inferred transmission pathway between another pair of hosts that goes through the node.")
  arg_parser$add_argument("-s", "--splitsRule", action="store", default="r", help="The rules by which the sets of patients are split into groups in order to ensure that all groups can be members of connected subtrees without causing conflicts. Currently available: c=conservative, r=Romero-Severson (default).")
  arg_parser$add_argument("-b", "--blacklist", action="store", help="A .csv file listing tips to ignore.")
  arg_parser$add_argument("inputFile", help="Input tree file name", metavar="inputTreeFileName")
  arg_parser$add_argument("outputFileIdentifier", help="A string identifying output files.")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
  arg_parser$add_argument("-O", "--outputdir", action="store", help="Full path of the directory for output; if absent, current working directory")
  arg_parser$add_argument("-pw", "--pdfwidth", action="store", default=100, help="Width of tree pdf in inches.")
  arg_parser$add_argument("-ph", "--pdfrelheight", action="store", default=0.15, help="Relative height of tree pdf.")
  arg_parser$add_argument("-db", "--debug", action="store_true", default=FALSE, help="Debugging mode.")
  
  args <- arg_parser$parse_args()
  script.dir <- args$scriptdir
  zero.length.tips.count <- args$zeroLengthTipsCount
  file.name <- args$inputFile
  output.dir <- args$outputDir
  if(is.null(output.dir)){
    output.dir <- getwd()
  }
  out.identifier <- args$outputFileIdentifier
  blacklist.file <- args$blacklist
  root.name <- args$outgroupName
  tip.regex <- args$tipRegex
  mode <- args$splitsRule
  pdf.hm <- as.numeric(args$pdfrelheight)
  pdf.w <- as.numeric(args$pdfwidth)
  debug <- args$debug
  
  if(!(mode %in% c("c", "r"))){
    stop(paste("Unknown split classifier: ", mode, "\n", sep=""))
  }
#  if(debug)
#	  cat(paste(script.dir, zero.length.tips.count, file.name, out.identifier, blacklist.file, root.name, tip.regex, mode, pdf.hm, pdf.w, sep='\n'))
  
  
} else {
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517_clean/")
  output.dir <- getwd()
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools/"
  file.name <- "RAxML_bestTree.InWindow_800_to_1150.tree"
  blacklist.file <- "FullBlacklist_InWindow_800_to_1150.csv"
  out.identifier <- "test"
  root.name <- "C.BW.00.00BW07621.AF443088"
  tip.regex <- "^(.*)_read_([0-9]+)_count_([0-9]+)$"
  mode <- "c"
  zero.length.tips.count <- F
  if(0)
  {
	  script.dir		<- '/Users/Oliver/git/phylotypes/tools'
	  file.name 		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptoutput/ptyr115_InWindow_800_to_1049.tree'
	  blacklist.file	<- NULL
	  #blacklist.file <- "FullBlacklist_InWindow_6800_to_7150.csv"
	  out.identifier 			<- "/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptoutput/ptyr115_InWindow_800_to_1049"
	  root.name 		<- "REF_CPX_AF460972_read_1_count_0"
	  tip.regex 		<- "^(.*)_read_([0-9]+)_count_([0-9]+)$"	  
	  mode 				<- "r"
	  pdf.hm 			<- 0.15
	  pdf.w 			<- 30	  
	  zero.length.tips.count	<- FALSE
  }
}

require(phangorn, quietly=T)
require(argparse, quietly=T)
require(phytools, quietly=T)
require(ggplot2, quietly=T)
require(ggtree, quietly=T)

source(file.path(script.dir, "TransmissionUtilityFunctions.R"))
source(file.path(script.dir, "SubtreeMethods.R"))

cat("SplitPatientsToSubtrees.R run on: ", file.name,", rules = ",mode,"\n", sep="")

tree <- read.tree(file.name)
tree <- unroot(tree)
tree <- di2multi(tree, tol = 1E-5)
if(!is.null(root.name)){
  tree <- root(tree, outgroup = which(tree$tip.label==root.name), resolve.root = T)
}

tip.labels <- tree$tip.label

blacklist <- vector()

if(!is.null(blacklist.file)){
  if(file.exists(blacklist.file)){
    blacklisted.tips <- read.table(blacklist.file, sep=",", header=F, stringsAsFactors = F, col.names="read")
    if(nrow(blacklisted.tips)>0){
      blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
    }
  } else {
    warning(paste("File ",blacklist.file.name," does not exist; skipping.",paste=""))
  }
} 



patient.ids <- sapply(tip.labels, function(x) patient.from.label(x, tip.regex))

if(length(blacklist)>0){
  patients <- unique(patient.ids[-blacklist])
} else {
  patients <- unique(patient.ids)
}

patients <- patients[!is.na(patients)]

patient.tips <-
  sapply(patients, function(x)  setdiff(which(patient.ids==x), blacklist))

patient.mrcas <- lapply(patient.tips, function(node) mrca.phylo.or.unique.tip(tree, node, zero.length.tips.count))

results <- split.and.annotate(tree, patients, patient.tips, patient.mrcas, blacklist, tip.regex, mode)

node.shapes <- rep(FALSE, length(tree$tip.label) + tree$Nnode)
for(mrca in results$first.nodes){
  node.shapes[mrca] <- TRUE
}

cat("Drawing tree...\n")

temp.ca <- rep(NA, length(tree$tip.label) + tree$Nnode)

for(item in seq(1, length(results$assocs))){
  if(!is.null(results$assocs[[item]])){
    if(results$assocs[[item]] != "*"){
      temp.ca[item] <- results$assocs[[item]]
    }
  }
}

temp.ca.pat <- sapply(temp.ca, function(x) unlist(strsplit(x, "-S"))[1] )
names(temp.ca.pat) <- NULL
temp.ca.pat <- factor(temp.ca.pat, levels = sample(levels(as.factor(temp.ca.pat))))
attr(tree, 'INDIVIDUAL') <- temp.ca.pat
attr(tree, 'NODE_SHAPES') <- node.shapes

tree.display <- ggtree(tree, aes(color=INDIVIDUAL)) +
  geom_point2(shape = 16, size=3, aes(subset=NODE_SHAPES)) +
  scale_fill_hue(na.value = "black") +
  scale_color_hue(na.value = "black") +
  theme(legend.position="none") +
  geom_tiplab(aes(col=INDIVIDUAL))

tree.display

ggsave(file.path(output.dir,paste('Tree_',mode,'_',out.identifier,'.pdf',sep='')), device="pdf", 
       height = pdf.hm*length(tree$tip.label), width = pdf.w, limitsize = F)

cat("Writing output...\n")

orig.patients <- vector()
patient.splits <- vector()
tip.names <- vector()

for(patient.plus in results$split.patients){
  tips <- results$split.tips[[patient.plus]]
  orig.patients <- c(orig.patients, rep(unlist(strsplit(patient.plus, "-S"))[1], length(tips)))
  patient.splits <- c(patient.splits, rep(patient.plus, length(tips)))
  tip.names <- c(tip.names, tree$tip.label[tips])
}

rs.subtrees <- data.frame(orig.patients, patient.splits, tip.names)


rs.file.name <- file.path(output.dir, paste('Subtrees_',mode,'_',out.identifier,'.csv',sep=''))
write.csv(rs.subtrees, rs.file.name, row.names = F, quote=F)
rs.file.name <- file.path(output.dir,paste('Subtrees_',mode,'_',out.identifier,'.rda',sep=''))
save(rs.subtrees, tree, file=rs.file.name)
  
