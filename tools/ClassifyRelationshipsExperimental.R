#!/usr/bin/env Rscript

command.line <- T

list.of.packages <- c("phangorn", "argparse", "phytools", "OutbreakTools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

if(!("ggtree" %in% installed.packages()[,"Package"])){
  source("https://bioconductor.org/biocLite.R")
  biocLite("ggtree")
}

if(command.line){
  suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
  
  arg_parser = ArgumentParser(description="Identify the likely direction of transmission between all patients present in a phylogeny.")
  arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                          help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name. Necessary only if -s is specified.")
  arg_parser$add_argument("-t", "--dualInfectionThreshold", type="double", help="Length threshold at which a branch of the subtree constructed just from reads from a single patient is considered long enough to indicate a dual infection (such a patient will be ignored, at present). If absent, keep all patients.")
  arg_parser$add_argument("-c", "--collapsedTree", action="store", help="If present, the collapsed tree (in which all adjacent nodes with the same assignment are collapsed to one) is output as a .csv file to the path specified.")
  arg_parser$add_argument("treeFileName", action="store", help="Tree file name. Alternatively, a base name that identifies a group of tree file names can be specified. Tree files are assumed to end in '.tree'. This must be the 'processed' tree produced by SplitTreesToSubtrees.R.")
  arg_parser$add_argument("splitsFileName", action="store", help="Splits file name. Alternatively, a base name that identifies a group of split file names can be specified. Split files are assumed to end in '_subtree_[a-z].csv'.")
  arg_parser$add_argument("outputFileName", action="store", help="Output file name (.csv format)")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
  
  args <- arg_parser$parse_args()
  
  tree.file.names <- args$treeFileName
  splits.file.names <- args$splitsFileName
  collapsed.file.names <- args$collapsedTree
  script.dir <- args$scriptdir
  output.name <- args$outputFileName
  tip.regex <- args$tipRegex
  split.threshold <- args$dualInfectionThreshold
  if(is.null(split.threshold)){
    split.threshold <- Inf
  }
  zero.length.tips.count <- args$zeroLengthTipsCount
  
} else {
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools"
  
  # BEEHIVE example
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517_clean/")
  tree.file.names <- "ProcessedTree_r_run20160517_inWindow_800_to_1150.tree"
  splits.file.names <- "Subtrees_r_run20160517_inWindow_800_to_1150.csv"
  
  output.name <- "hi.csv"
  collapsed.file.names <- "collapsed.csv"
  split.threshold <- NA
  
  # Rakai example
  
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161007_couples_w270_rerun/")
  tree.file.names <- "ProcessedTree_r_test_pytr5.tree"
  splits.file.names <- "Subtrees_r_test_pytr5.csv"
  
  output.name <- "hi.csv"
  collapsed.file.names <- "collapsed.csv"
  split.threshold <- NA
  
  #other example
  
  # setwd("/Users/twoseventwo/Documents/Work/pty_16-11-23-15-05-35/")
  # tree.file.names <- "ProcessedTree_r_ptyr22_InWindow_"
  # splits.file.names <- "Subtrees_r_ptyr22_InWindow_"
  # 
  # output.name <- "ptyr22_"
  # collapsed.file.names <- NULL
  # split.threshold <- NA
  
  # MRSA example
  
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/Thai MRSA 6/Matthew/")
  tree.file.names <- "ProcessedTree_s_mrsa_k50.tree"
  splits.file.names <- "Subtrees_s_mrsa_k50.csv"
  output.name <- "LT_s_mrsa_k.csv"
  collapsed.file.names <- "collapsed_s_mrsa_k50.csv"
  split.threshold <- NA
  tip.regex <- "^([ST][0-9][0-9][0-9])_[A-Z0-9]*_[A-Z][0-9][0-9]$"
  
  
  if(0)
  {
    script.dir				<- "/Users/Oliver/git/phylotypes/tools"
    tree.file.names			<- '/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptoutput/ptyr115_'
    splits.file.names		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptoutput/ptyr115_'
    output.name 			<- '/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptoutput/ptyr115_'
    root.name				<- "REF_CPX_AF460972_read_1_count_0"
    tip.regex 				<- "^(.*)_read_([0-9]+)_count_([0-9]+)$"	  		
    split.threshold			<- Inf 
    romero.severson 		<- TRUE
    zero.length.tips.count 	<- FALSE
  }
}

suppressMessages(library(phytools, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(phangorn, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ggtree, quietly=TRUE, warn.conflicts=FALSE))
#
#	load external functions
#
source(file.path(script.dir, "TransmissionUtilityFunctions.R"))
source(file.path(script.dir, "SubtreeMethods.R"))
#
#	define internal functions
#
likely.transmissions<- function(tree.file.name, splits.file.name, tip.regex, split.threshold, romero.severson, zero.length.tips.count, tt.file.name = NULL)
{	
  cat("Opening file: ", tree.file.name, "...\n", sep = "")
  
  pseudo.beast.import <- read.beast(tree.file.name)
  
  tree <- attr(pseudo.beast.import, "phylo")
  
  cat("Reading splits file",splits.file.name,"...\n")
  
  splits <- read.csv(splits.file.name, stringsAsFactors = F)
  
  cat("Collecting tips for each patient...\n")
  
  patients <- unique(splits$orig.patients)
  patient.tips <- lapply(patients, function(x) which(tree$tip.label %in% splits[which(splits$orig.patients==x), "tip.names"]))
  names(patient.tips) <- patients
  all.splits <- unique(splits$patient.splits)
  
  cat("Reading annotations...\n")
  
  annotations <- attr(pseudo.beast.import, "stats")
  annotations$INDIVIDUAL <- sapply(as.character(annotations$INDIVIDUAL), function(x) substr(x, 2, nchar(x)-1))
  annotations$SPLIT <- sapply(as.character(annotations$SPLIT), function(x) substr(x, 2, nchar(x)-1))
  
  was.split <- unique(annotations$INDIVIDUAL[which(annotations$SPLIT!=annotations$INDIVIDUAL)])
  
  in.order <- match(seq(1, length(tree$tip.label) + tree$Nnode), annotations$node)
  
  assocs <- annotations$SPLIT[in.order]
  assocs.ind <- annotations$INDIVIDUAL[in.order]
  
  assocs <- lapply(assocs, function(x) replace(x,is.na(x),"none"))
  assocs.ind <- lapply(assocs.ind, function(x) replace(x,is.na(x),"none"))
  
  splits.for.patients <- lapply(patients, function(x) unique(splits$patient.splits[which(splits$orig.patients==x)] ))
  names(splits.for.patients) <- patients
  
  patients.for.splits <- lapply(all.splits, function(x) unique(splits$orig.patients[which(splits$patient.splits==x)] ))
  names(patients.for.splits) <- all.splits
  
  # get rid of what look like dual infections
  
  ignore.because.split <- vector()
  
  # TODO currently this checks only split patients and it should probably check all patients?
  if(!is.na(split.threshold)){
    for(split.patient in was.split){
      tips <- patient.tips[[split.patient]]
      subtree <- drop.tip(tree, tree$tip.label[-patient.tips[[split.patient]]])
      max.length <- max(subtree$edge.length)
      if(max.length > split.threshold){
        ignore.because.split <- c(ignore.because.split, split.patient)
      }
    }
  }
  
  patients.included <- setdiff(patients, ignore.because.split)
  
  total.pairs <- (length(patients.included) ^ 2 - length(patients.included))/2
  
  cat("Collapsing subtrees...\n")
  
  tt <- output.trans.tree(tree, assocs, tt.file.name)
  
  cat("Identifying pairs of unblocked splits...\n")
  
  adjacent <- subtrees.unblocked(tt, all.splits)
  
  cat("Calculating pairwise distances between splits...\n")
  
  split.distances <- all.subtree.distances(tree, tt, all.splits, assocs)
  
  cat("Testing pairs...\n")
  
  count <- 0
  contiguous.matrix <- matrix(NA, length(patients.included), length(patients.included))
  path.matrix <- matrix(NA, length(patients.included), length(patients.included))
  nodes.12.matrix <- matrix(NA, length(patients.included), length(patients.included))
  nodes.21.matrix <- matrix(NA, length(patients.included), length(patients.included))
  dir.12.matrix <- matrix(NA, length(patients.included), length(patients.included))
  dir.21.matrix <- matrix(NA, length(patients.included), length(patients.included))
  mean.distance.matrix <- matrix(NA, length(patients.included), length(patients.included))
  
  tree.distances <- dist.nodes(tree)
  
  for(pat.1 in seq(1, length(patients.included))){
    for(pat.2 in  seq(1, length(patients.included))){
      if (pat.1 < pat.2) {
        count <- count + 1
        
        pat.1.id <- patients.included[pat.1]
        pat.2.id <- patients.included[pat.2]
        
        #	      cat(pat.1.id, pat.2.id, "\n")
        
        nodes.1 <- splits.for.patients[[pat.1.id]]
        nodes.2 <- splits.for.patients[[pat.2.id]]
        
        all.nodes <-  c(nodes.1, nodes.2)
        
        # we want that they form a perfect blob with no intervening nodes for any other patients
        OK <- check.contiguous(tt, c(pat.1.id, pat.2.id), splits.for.patients, patients.for.splits)		
        
        contiguous.matrix[pat.1, pat.2] <- OK
      
        count.12 <- 0
        count.21 <- 0
        
        for(node.2 in nodes.2){
          ancestors <- get.tt.ancestors(tt, node.2)
          if(length(intersect(ancestors, nodes.1)) > 0){
            count.12 <- count.12 + 1
          }
        }
        
        for(node.1 in nodes.1){
          ancestors <- get.tt.ancestors(tt, node.1)
          if(length(intersect(ancestors, nodes.2)) > 0){
            count.21 <- count.21 + 1
          }
        }
        
        prop.12 <- count.12/length(nodes.2)
        prop.21 <- count.21/length(nodes.1)
        
        nodes.12.matrix[pat.1, pat.2] <- length(nodes.1)
        nodes.21.matrix[pat.1, pat.2] <- length(nodes.2)
        
        dir.12.matrix[pat.1, pat.2] <- count.12
        dir.21.matrix[pat.1, pat.2] <- count.21
        
        if(count.12 == 0 & count.21 == 0){
          path.matrix[pat.1, pat.2] <- "none"
        } else if(count.12 != 0 & count.21 == 0 & prop.12 == 1) {
          if(count.12 == 1){
            path.matrix[pat.1, pat.2] <- "anc"
          } else {
            path.matrix[pat.1, pat.2] <- "multiAnc"
          }
        } else if(count.21 != 0 & count.12 == 0 & prop.21 == 1) {
          if(count.21 == 1){
            path.matrix[pat.1, pat.2] <- "desc"
          } else {
            path.matrix[pat.1, pat.2] <- "multiDesc"
          }
        } else {
          path.matrix[pat.1, pat.2] <- "conflict"
        }
        
        pairwise.distances <- 0
        
        for(node.1 in nodes.1){
          for(node.2 in nodes.2){
            if(adjacent[node.1, node.2]){
              pairwise.distances <- c(pairwise.distances, split.distances[node.1, node.2])
            }
          }
        }
        
        mean.distance.matrix[pat.1, pat.2] <- mean(pairwise.distances)
        
        
        if (count %% 100 == 0) {
          cat(paste("Done ",format(count, scientific = F)," of ",format(total.pairs, scientific = F)," pairwise calculations\n", sep = ""))
        }
      }
    }
  }
  
  contiguous.table <- as.table(contiguous.matrix)
  dir.12.table <- as.table(dir.12.matrix)
  dir.21.table <- as.table(dir.21.matrix)
  nodes.12.table <- as.table(nodes.12.matrix)
  nodes.21.table <- as.table(nodes.21.matrix)
  path.table <- as.table(path.matrix)
  mean.distance.table <- as.table(mean.distance.matrix)
  
  colnames(contiguous.table) <- patients.included
  rownames(contiguous.table) <- patients.included
  
  colnames(path.table) <- patients.included
  rownames(path.table) <- patients.included
  
  colnames(dir.12.table) <- patients.included
  rownames(dir.12.table) <- patients.included
  
  colnames(dir.21.table) <- patients.included
  rownames(dir.21.table) <- patients.included
  
  colnames(nodes.12.table) <- patients.included
  rownames(nodes.12.table) <- patients.included
  
  colnames(nodes.21.table) <- patients.included
  rownames(nodes.21.table) <- patients.included
  
  colnames(mean.distance.table) <- patients.included
  rownames(mean.distance.table) <- patients.included
  
  cdf <- as.data.frame(contiguous.table)
  
  keep <- complete.cases(cdf)
  
  cdf <- cdf[keep,]

  pdf <- as.data.frame(path.table)
  pdf <- pdf[keep,]
  
  d12df <- as.data.frame(dir.12.table)
  d21df <- as.data.frame(dir.21.table)
  d12df <- d12df[keep,]
  d21df <- d21df[keep,]
  
  n12df <- as.data.frame(nodes.12.table)
  n21df <- as.data.frame(nodes.21.table)
  n12df <- n12df[keep,]
  n21df <- n21df[keep,]
  
  mddf <- as.data.frame(mean.distance.table)
  mddf <- mddf[keep,]
  
  cdf <- cbind(cdf, d12df[,3], d21df[,3], n12df[,3], n21df[,3], pdf[,3], mddf[,3])
  
  colnames(cdf) <- c("Patient_1", "Patient_2", "contiguous", "paths.12", "paths21", "nodes12", "nodes21", "path.classification", "mean.distance.between.subtrees")
  cdf
}

#
#	start script
#

#	check if 'tree.file.names' is tree

single.file	<- file.exists(tree.file.names)

if(single.file)
{
  #	if 'tree.file.names' is tree, process just one tree
  tree.file.name		<- tree.file.names[1]
  splits.file.name	<- splits.file.names[1]
  dddf				<- likely.transmissions(tree.file.name, splits.file.name, tip.regex, split.threshold, romero.severson, zero.length.tips.count, collapsed.file.names)
  cat("Write to file",output.name,"...\n")
  write.table(dddf, file = output.name, sep = ",", row.names = FALSE, col.names = TRUE, quote=F)	
} else {
  #	if 'tree.file.names' is not tree, process multiple trees
  tree.file.names		<- sort(list.files(dirname(tree.file.names), pattern=paste(basename(tree.file.names),'.*\\.tree$',sep=''), full.names=TRUE))
  splits.file.names	<- sort(list.files(dirname(splits.file.names), pattern=paste(basename(splits.file.names),'.*\\.csv$',sep=''), full.names=TRUE))	
  stopifnot(length(tree.file.names)==length(splits.file.names))
  for(tree.i in seq_along(tree.file.names))
  {
    tree.file.name		<- tree.file.names[tree.i]
    splits.file.name	<- splits.file.names[tree.i]
    output.name			<- gsub('\\.tree','_LikelyTransmissions.csv', tree.file.name)		
    collapsed.file.name	<- gsub('\\.tree','_collapsed.csv', tree.file.name)
    dddf				<- likely.transmissions(tree.file.name, splits.file.name, tip.regex, split.threshold, romero.severson, zero.length.tips.count, collapsed.file.name)
    cat("Write to file",output.name,"...\n")
    write.table(dddf, file = output.name, sep = ",", row.names = FALSE, col.names = TRUE, quote=F)
  }
}