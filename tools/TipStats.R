#!/usr/bin/env Rscript

command.line <- T

if(command.line){
  suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
  
  arg_parser = ArgumentParser(description="Calculates per-tip statistics, currently just RTT and SRTT (subgraph root-to-tip)")
  arg_parser$add_argument("treeFileName", action="store", help="Tree file name. Alternatively, a base name that identifies a group of tree file names can be specified. Tree files are assumed to end in '.tree'. This must be the 'processed' tree produced by SplitTreesToSubtrees.R.")
  arg_parser$add_argument("splitsFileName", action="store", help="Splits file name. Alternatively, a base name that identifies a group of split file names can be specified. Split files are assumed to end in '_subtree_[a-z].csv'.")
  arg_parser$add_argument("outputFileName", action="store", help="Output file name (.csv format) or a base name identifying multiple output files.")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
  
  args <- arg_parser$parse_args()
  
  tree.file.names <- args$treeFileName
  splits.file.names <- args$splitsFileName
  script.dir <- args$scriptdir
  output.name <- args$outputFileName

} else {
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools"
  
  # BEEHIVE example
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20161013/")
  tree.file.names <- "ProcessedTrees_r/ProcessedTree_r_run20161013_RS_inWindow_4550_to_4900.tree"
  splits.file.names <- "SubtreeFiles_r/Subtrees_r_run20161013_RS_inWindow_4550_to_4900.csv"

  output.name <- "outTest.csv"

}

suppressMessages(library(phytools, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(phyloscannerR, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(phangorn, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ggtree, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(data.table, quietly=TRUE, warn.conflicts=FALSE))
#
#	load external functions
#
#
#	define internal functions
#
tip.stats <- function(tree.file.name, splits.file.name)
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
  
  stats	<- data.table(TIP=rep(NA_character_,length(tree$tip.label)), SPLIT = rep(NA_character_,length(tree$tip.label)), PATIENT = rep(NA_character_,length(tree$tip.label)), RTT=rep(NA_real_,length(tree$tip.label)), SRTT=rep(NA_real_,length(tree$tip.label)))
  
  dd	<- lapply(seq(1, length(tree$tip.label)), function(index){
    rtt <- 0
    current.node <- index
    
    while(!is.root(tree, current.node)){
      rtt <- rtt + get.edge.length(tree, current.node)
      current.node <- Ancestors(tree, current.node, type="parent")
    }
    
    if(assocs[[index]]=="none"){
      stats[index, 1:5 := list(tree$tip.label[index], "blacklisted", "blacklisted", rtt, NA)]
    } else {
      srtt <- 0
      current.node <- index
      split <- assocs[[index]]
      individual <- assocs.ind[[index]]
      
      while(!is.root(tree, current.node) & assocs[current.node]==split){
        srtt <- srtt + get.edge.length(tree, current.node)
        current.node <- Ancestors(tree, current.node, type="parent")
      } 
      stats[index, 1:5 := list(tree$tip.label[index], split, individual, rtt, srtt)]
      
      

    }
  })
  
  stats
}

#
#	start script
#

#	check if 'tree.file.names' is tree

single.file	<- file.exists(tree.file.names)

if(single.file){
  #	if 'tree.file.names' is tree, process just one tree
  tree.file.name		<- tree.file.names[1]
  splits.file.name	<- splits.file.names[1]
  dddf	<- tip.stats(tree.file.name, splits.file.name)
  cat("Write to file",output.name,"...\n")
  write.csv(dddf, file = output.name, row.names = FALSE, quote=F)	
} else {
  #	if 'tree.file.names' is not tree, process multiple trees
  tree.file.names		<- sort(list.files(dirname(tree.file.names), pattern=paste(basename(tree.file.names),'.*\\.tree$',sep=''), full.names=TRUE))
  splits.file.names	<- sort(list.files(dirname(splits.file.names), pattern=paste(basename(splits.file.names),'.*\\.csv$',sep=''), full.names=TRUE))	
  stopifnot(length(tree.file.names)==length(splits.file.names))
  for(tree.i in seq_along(tree.file.names))
  {
    tree.file.name		<- tree.file.names[tree.i]
    splits.file.name	<- splits.file.names[tree.i]
    output.file.name			<- gsub('\\.tree',paste('_',output.name,'.csv',sep=''), tree.file.name)
	tryCatch({
				dddf	<- tip.stats(tree.file.name, splits.file.name)
				cat("Write to file",output.name,"...\n")
				write.csv(dddf, file = output.file.name, row.names = FALSE, quote=F)			
			},
			error=function(e){ print(e$message); cat("No output written to file.\n") })    
  }
}