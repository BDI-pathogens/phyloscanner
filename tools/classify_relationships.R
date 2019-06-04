#!/usr/bin/env Rscript

suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(phyloscannerR, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(phytools, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(phangorn, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ggtree, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(kimisc, quietly=TRUE, warn.conflicts=FALSE))

arg_parser = ArgumentParser(description="Classify relationships between study subjects based on their relative positions in the phylogeny")
arg_parser$add_argument("-c", "--collapsedTree", action="store_true", help="If present, the collapsed tree (in which all adjacent nodes with the same assignment are collapsed to one) is output as a CSV file or files.")
arg_parser$add_argument("-n", "--branchLengthNormalisation", action="store", help="If present and a number, a normalising constant for all branch lengths in the tree or trees. If present and a file, the path to a .csv file with two columns: tree file name and normalising constant")
arg_parser$add_argument("treeFileName", action="store", help="Tree file name. Alternatively, a base name that identifies a group of tree file names. Tree files are assumed to end in '.tree'. This must be the 'processed' tree produced by SplitTreesToSubtrees.R.")
arg_parser$add_argument("splitsFileName", action="store", help="Splits file name. Alternatively, a base name that identifies a group of split file.")
arg_parser$add_argument("outputFileName", action="store", help="A path and root character string identifiying all output files.")
arg_parser$add_argument("-tfe", "--treeFileExtension", action="store", default="tree", help="The file extension for tree files (default .tree).")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what the script is doing.")
args <- arg_parser$parse_args()

tree.file.name           <- args$treeFileName
splits.file.name         <- args$splitsFileName
do.collapsed             <- args$collapsedTree
output.root              <- args$outputFileName
has.normalisation        <- !is.null(args$branchLengthNormalisation)
normalisation.argument   <- args$branchLengthNormalisation
tree.fe                  <- args$treeFileExtension
csv.fe                   <- args$csvFileExtension
verbose                  <- args$verbose

#	load external functions

all.tree.info                        <- list()

if(file.exists(tree.file.name)){

  tree.info                          <- list()
  tree.info$tree.file.name           <- tree.file.name
  tree.info$splits.file.names        <- splits.file.name
  tree.info$classification.file.name <- paste(output.root, "_classification", ".", csv.fe, sep="")
  if(do.collapsed){
    tree.info$collapsed.file.name    <- paste(output.root, "_collapsedTree", ".", csv.fe, sep="")
  }
  
  if(!has.normalisation){
    tree.info$normalisation.constant <- 1
  } else {
    nc <- suppressWarnings(as.numeric(normalisation.argument))
    
    if(is.na(nc)){
      norm.table <- read.csv(normalisation.argument, stringsAsFactors = F, header = F)
      if(nrow(norm.table[which(norm.table[,1]==basename(tree.file.name)),])==1){
        tree.info$normalisation.constant <- as.numeric(norm.table[which(norm.table[,1]==basename(tree.file.name)),2])
      } else if (nrow(norm.table[which(norm.table[,1]==basename(tree.file.name)),])==0){
        warning(paste("No normalisation constant given for tree file name ", basename(tree.file.name), " in file ",normalisation.argument,"; normalised distances will not be given", sep=""))
        tree.info$normalisation.constant <- 1
      } else {
        warning(paste("Two or more entries for normalisation constant for tree file name ",basename(tree.file.name)," in file ",normalisation.argument,"; taking the first", sep=""))
        tree.info$normalisation.constant <- as.numeric(norm.table[which(norm.table[,1]==basename(tree.file.name))[1],2])
      }
    } else {
      tree.info$normalisation.constant <- nc
    }
  }
  
  all.tree.info[[1]] <- tree.info
} else {
  # Assume we are dealing with a group of files
  
  tree.file.names		<- list.files.mod(dirname(tree.file.name), pattern=paste(basename(tree.file.name),'.*',tree.fe,'$',sep=''), full.names=TRUE)
  splits.file.names		<- list.files.mod(dirname(splits.file.name), pattern=paste(basename(splits.file.name),'.*',csv.fe,'$',sep=''), full.names=TRUE)
  
  if(length(tree.file.names)!=length(splits.file.names)){
    stop("Numbers of tree files and subgraph files differ\n")
  }
  
  if(length(tree.file.names)==0){
    stop("No input trees found.\n Quitting.\n")
  }
  
  suffixes <- substr(tree.file.names, nchar(tree.file.name) + 1, nchar(tree.file.names) - nchar(tree.fe) - 1)
  suffixes <- suffixes[order(suffixes)]
  
  tree.file.names <- tree.file.names[order(suffixes)]
  splits.file.names <- splits.file.names[order(suffixes)]
  
  classification.file.names <- paste(output.root, "_classification_", suffixes, ".", csv.fe, sep="")
  if(do.collapsed){
    collapsed.file.names <- paste(output.root, "_collapsedTree_",suffixes, ".", csv.fe, sep="")
  }

  if(has.normalisation){
    normalisation.constant <- suppressWarnings(as.numeric(normalisation.argument))    
    if(is.na(normalisation.constant)){
      norm.table <- read.csv(normalisation.argument, stringsAsFactors = F, header = F)	  
	  normalisation.constants <- sapply(basename(tree.file.names), function(x){				  
		tmp <- 	which(norm.table[,1]==x)
		if(length(tmp)==0)
			tmp	<- which(norm.table[,1]==gsub('ProcessedTree_[a-z]_','',x))	  
        if(length(tmp)>=1){
          if(length(tmp)>1){
            warning(paste("Two or more entries for normalisation constant for tree file name ",x, " in file ",normalisation.argument,"; taking the first", sep=""))
          }
          return(norm.table[tmp[1],2])
        } else {
          warning(paste("No normalisation constant for tree file ",x," found in file, ",normalisation.argument,"; this tree will not have branch lengths normalised.", sep=""))
          return(1)
        }
        
      })        
    } else {
      # it's a number to be used in every tree
      normalisation.constants <- rep(normalisation.constant, length(tree.file.names))
    }
  } else {
    normalisation.constants <- rep(1, length(tree.file.names))
  }
  
  fn.df <- data.frame(row.names = suffixes, 
                      tree.file.name = tree.file.names, 
                      splits.file.name = splits.file.names, 
                      normalisation.constant = normalisation.constants, 
                      classification.file.name = classification.file.names,
                      stringsAsFactors = F)
  if(do.collapsed){
    fn.df$collapsed.file.name <- collapsed.file.names
  }
  all.tree.info <- split(fn.df, rownames(fn.df))
}
  
for(tree.info in all.tree.info){	

  tree.info <- as.list(tree.info)
  
  results <- classify(tree.info, verbose)
  if (verbose) cat("Writing to file",tree.info$classification.file.name,"...\n")
  write.table(results$classification, file = tree.info$classification.file.name, sep = ",", row.names = FALSE, col.names = TRUE, quote=F)	
  
  if(!is.null(tree.info$collapsed.file.name)){
    # clean it
    nice.tt.output <- prune.unsampled.tips(results$collapsed)
    write.table(nice.tt.output, file = tree.info$collapsed.file.name, sep = ",", row.names = FALSE, col.names = TRUE, quote=F)	
  }
} 
