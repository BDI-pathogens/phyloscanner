#!/usr/bin/env Rscript

command.line <- T

list.of.packages <- c("phangorn", "argparse", "phytools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

if(!("ggtree" %in% installed.packages()[,"Package"])){
  source("https://bioconductor.org/biocLite.R")
  biocLite("ggtree")
}

tree.fe <- ".tree"
csv.fe <- ".csv"

if(command.line){
  suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
  
  arg_parser = ArgumentParser(description="Classify relationships between study subjects based on their relative positions in the phylogeny")
  arg_parser$add_argument("-c", "--collapsedTree", action="store_true", default=T, help="If present, the collapsed tree (in which all adjacent nodes with the same assignment are collapsed to one) is output as a .csv file or files.")
  arg_parser$add_argument("-n", "--branchLengthNormalisation", action="store", help="If present and a number, a normalising constant for all branch lengths in the tree or trees. If present and a file, the path to a .csv file with two columns: tree file name and normalising constant")
  arg_parser$add_argument("treeFileName", action="store", help="Tree file name. Alternatively, a base name that identifies a group of tree file names. Tree files are assumed to end in '.tree'. This must be the 'processed' tree produced by SplitTreesToSubtrees.R.")
  arg_parser$add_argument("splitsFileName", action="store", help="Splits file name. Alternatively, a base name that identifies a group of split file.")
  arg_parser$add_argument("outputFileName", action="store", help="A path and root character string identifiying all output files.")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
  arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")
  args <- arg_parser$parse_args()
  
  tree.file.name<- args$treeFileName
  splits.file.name <- args$splitsFileName
  do.collapsed <- args$collapsedTree
  script.dir <- args$scriptdir
  output.root <- args$outputFileName
  has.normalisation <- !is.null(args$branchLengthNormalisation)
  normalisation.argument <- args$branchLengthNormalisation
  verbose <- args$verbose
  
} else {
  script.dir <- "/Users/mdhall/phylotypes/tools"
  
  # BEEHIVE example
  setwd("/Users/mdhall/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20161013/")
  tree.file.names <- "ProcessedTrees_s/ProcessedTree_s_run20161013D_inWindow_800_to_1150.tree"
  splits.file.names <- "SubtreeFiles_s/Subtrees_s_run20161013D_inWindow_800_to_1150.csv"
  
  output.name <- "outTest.csv"
  collapsed.file.names <- "outCollapsed.csv"
  
  # Rakai example
  
  # setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161007_couples_w270_rerun/")
  # tree.file.names <- "ProcessedTree_r_test_pytr5.tree"
  # splits.file.names <- "Subtrees_r_test_pytr5.csv"
  # 
  # output.name <- "hi.csv"
  # collapsed.file.names <- "collapsed.csv"
  
  #other example
  
  # setwd("/Users/twoseventwo/Documents/Work/pty_16-11-23-15-05-35/")
  # tree.file.names <- "ProcessedTree_r_ptyr22_InWindow_"
  # splits.file.names <- "Subtrees_r_ptyr22_InWindow_"
  # 
  # output.name <- "ptyr22_"
  # collapsed.file.names <- NULL
  
  # MRSA example
  
  setwd("/Users/mdhall/Dropbox (Infectious Disease)/Thai MRSA 6/Matthew/refinement")
  tree.file.names <- "ProcessedTree_s_mrsa_k10_bp_yetanother.tree"
  splits.file.names <- "Subtrees_s_mrsa_k10_bp_yetanother.csv"
  output.file.ID <- "LT_s_mrsa_k10_yetanother.csv"
  collapsed.file.ID <- "collapsed_s_mrsa_k10_yetanother.csv"
  
  
  if(0)
  {
    setwd("/Users/twoseventwo/Library/Containers/com.apple.mail/Data/Library/Mail Downloads/8D40D980-06B2-46EE-9BA2-789D541F19E2")
    script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools"
    tree.file.names			<- '/Users/twoseventwo/Library/Containers/com.apple.mail/Data/Library/Mail Downloads/8D40D980-06B2-46EE-9BA2-789D541F19E2/pty_17-03-22-16-12-38/ProcessedTree_s_ptyr22_InWindow_'
    splits.file.names		<- '/Users/twoseventwo/Library/Containers/com.apple.mail/Data/Library/Mail Downloads/8D40D980-06B2-46EE-9BA2-789D541F19E2/pty_17-03-22-16-12-38/subgraphs_s_ptyr22_InWindow_'
    output.name 			<- 'testpit'
    has.normalisation <- F
  }
}

suppressMessages(library(phytools, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(phangorn, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ggtree, quietly=TRUE, warn.conflicts=FALSE))
#
#	load external functions
#
source(file.path(script.dir, "TreeUtilityFunctions.R"))
source(file.path(script.dir, "ParsimonyReconstructionMethods.R"))
source(file.path(script.dir, "CollapsedTreeMethods.R"))
source(file.path(script.dir, "GeneralFunctions.R"))

#	define internal functions
#
classify <- function(tree.file.name, splits.file.name, normalisation.constant = 1, tt.file.name = NULL) {	
  if (verbose) cat("Reading tree file ", tree.file.name, "...\n", sep = "")
  
  pseudo.beast.import <- read.beast(tree.file.name)
  tree <- attr(pseudo.beast.import, "phylo")
  
  if (verbose) cat("Reading splits file",splits.file.name,"...\n")
  
  splits <- read.csv(splits.file.name, stringsAsFactors = F)

  # this should go eventually
  
  colnames(splits) <- c("patient", "subgraph", "split")
  
  if (verbose) cat("Collecting tips for each patient...\n")
  
  patients <- unique(splits$patient)

  all.splits <- unique(splits$subgraph)
  
  if (verbose) cat("Reading annotations...\n")
  
  annotations <- attr(pseudo.beast.import, "stats")

  annotations$INDIVIDUAL <- as.character(annotations$INDIVIDUAL)
  annotations$SPLIT <- as.character(annotations$SPLIT)
  
  in.order <- match(seq(1, length(tree$tip.label) + tree$Nnode), annotations$node)
  
  assocs <- annotations$SPLIT[in.order]
  assocs <- lapply(assocs, function(x) replace(x,is.na(x),"none"))

  splits.for.patients <- lapply(patients, function(x) unique(splits$subgraph[which(splits$patient==x)] ))
  names(splits.for.patients) <- patients
  
  patients.for.splits <- lapply(all.splits, function(x) unique(splits$patient[which(splits$subgraph==x)] ))
  names(patients.for.splits) <- all.splits
  
  patients.included <- patients
  
  total.pairs <- (length(patients.included) ^ 2 - length(patients.included))/2
  
  if (verbose) cat("Collapsing subtrees...\n")
  
  tt <- output.trans.tree(tree, assocs, tt.file.name)
  
  if (verbose) cat("Identifying pairs of unblocked splits...\n")
  
  collapsed.adjacent <- subtrees.unblocked(tt, all.splits)
  
  if (verbose) cat("Calculating pairwise distances between splits...\n")
  
  split.distances <- tryCatch(
    all.subtree.distances(tree, tt, all.splits, assocs), warning=function(w){return(NULL)}, error=function(e){return(NULL)})
  
  if(is.null(split.distances)){
    split.distances <- all.subtree.distances(tree, tt, all.splits, assocs, TRUE)
  }

  if (verbose) cat("Testing pairs...\n")
  
  count <- 0
  adjacency.matrix <- matrix(NA, length(patients.included), length(patients.included))
  contiguity.matrix <- matrix(NA, length(patients.included), length(patients.included))
  path.matrix <- matrix(NA, length(patients.included), length(patients.included))
  nodes.1.matrix <- matrix(NA, length(patients.included), length(patients.included))
  nodes.2.matrix <- matrix(NA, length(patients.included), length(patients.included))
  dir.12.matrix <- matrix(NA, length(patients.included), length(patients.included))
  dir.21.matrix <- matrix(NA, length(patients.included), length(patients.included))
  min.distance.matrix <- matrix(NA, length(patients.included), length(patients.included))
  
  for(pat.1 in seq(1, length(patients.included))){
    for(pat.2 in  seq(1, length(patients.included))){
      if (pat.1 < pat.2) {

        count <- count + 1
        
        pat.1.id <- patients.included[pat.1]
        pat.2.id <- patients.included[pat.2]
        
        nodes.1 <- splits.for.patients[[pat.1.id]]
        nodes.2 <- splits.for.patients[[pat.2.id]]
        
        all.nodes <-  c(nodes.1, nodes.2)
        
        adjacency.matrix[pat.1, pat.2] <- check.adjacency(tt, c(pat.1.id, pat.2.id), splits.for.patients)		
        contiguity.matrix[pat.1, pat.2] <- check.contiguous(tt, c(pat.1.id, pat.2.id), splits.for.patients, patients.for.splits)		

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
        
        nodes.1.matrix[pat.1, pat.2] <- length(nodes.1)
        nodes.2.matrix[pat.1, pat.2] <- length(nodes.2)
        
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
        
        pairwise.distances <- vector()
        
        for(node.1 in nodes.1){
          for(node.2 in nodes.2){
            if(collapsed.adjacent[node.1, node.2]){
              pairwise.distances <- c(pairwise.distances, split.distances[node.1, node.2])
            }
          }
        }
        
        min.distance.matrix[pat.1, pat.2] <- min(pairwise.distances)
        
        
        if (count %% 100 == 0) {
          if (verbose) cat(paste("Done ",format(count, scientific = F)," of ",format(total.pairs, scientific = F)," pairwise calculations\n", sep = ""))
        }
      }
    }
  }
  

  normalised.distance.matrix<- min.distance.matrix/normalisation.constant
  
  adjacency.table <- as.table(adjacency.matrix)
  contiguity.table <- as.table(contiguity.matrix)
  dir.12.table <- as.table(dir.12.matrix)
  dir.21.table <- as.table(dir.21.matrix)
  nodes.1.table <- as.table(nodes.1.matrix)
  nodes.2.table <- as.table(nodes.2.matrix)
  path.table <- as.table(path.matrix)
  min.distance.table <- as.table(min.distance.matrix)
  

  normalised.distance.table <- as.table(normalised.distance.matrix)

  
  colnames(adjacency.table) <- patients.included
  rownames(adjacency.table) <- patients.included
  
  colnames(contiguity.table) <- patients.included
  rownames(contiguity.table) <- patients.included
  
  colnames(path.table) <- patients.included
  rownames(path.table) <- patients.included
  
  colnames(dir.12.table) <- patients.included
  rownames(dir.12.table) <- patients.included
  
  colnames(dir.21.table) <- patients.included
  rownames(dir.21.table) <- patients.included
  
  colnames(nodes.1.table) <- patients.included
  rownames(nodes.1.table) <- patients.included
  
  colnames(nodes.2.table) <- patients.included
  rownames(nodes.2.table) <- patients.included
  
  colnames(min.distance.table) <- patients.included
  rownames(min.distance.table) <- patients.included
  
  colnames(normalised.distance.table) <- patients.included
  rownames(normalised.distance.table) <- patients.included
  
  adf <- as.data.frame(adjacency.table)
  keep <- complete.cases(adf)
  adf <- adf[keep,]
  
  cdf <- as.data.frame(contiguity.table)
  keep <- complete.cases(cdf)
  cdf <- cdf[keep,]
  
  pdf <- as.data.frame(path.table)
  pdf <- pdf[keep,]
  
  d12df <- as.data.frame(dir.12.table)
  d21df <- as.data.frame(dir.21.table)
  d12df <- d12df[keep,]
  d21df <- d21df[keep,]
  
  n1df <- as.data.frame(nodes.1.table)
  n2df <- as.data.frame(nodes.2.table)
  n1df <- n1df[keep,]
  n2df <- n2df[keep,]
  
  mddf <- as.data.frame(min.distance.table)
  mddf <- mddf[keep,]
  
  if(normalisation.constant!=1){
    nddf <- as.data.frame(normalised.distance.table)
    nddf <- nddf[keep,]
  }
  
  adf <- cbind(adf, cdf[,3], d12df[,3], d21df[,3], n1df[,3], n2df[,3], pdf[,3], mddf[,3])
  
  column.names <- c("Patient_1", "Patient_2", "adjacent", "contiguous", "paths12", "paths21", "nodes1", "nodes2", "path.classification", "min.distance.between.subtrees")
  
  if(normalisation.constant!=1){
    adf <- cbind(adf, nddf[,3])
    column.names <- c(column.names, "normalised.min.distance.between.subtrees")
  }
  
  colnames(adf) <- column.names
  adf
}

#
#	start script
#

#	check if 'tree.file.names' is tree


if(file.exists(tree.file.name)){
  tree.file.names <- tree.file.name
  splits.file.names <- splits.file.name
  output.file.names <- paste(output.root, "_classification", csv.fe, sep="")
  if(do.collapsed){
    collapsed.file.names <- paste(output.root, "_collapsedTree", csv.fe, sep="")
  }
  
  if(!has.normalisation){
    normalisation.constants <- 1
  } else {
    normalisation.constants <- suppressWarnings(as.numeric(normalisation.argument))
    
    if(is.na(normalisation.constants)){
      norm.table <- read.csv(normalisation.argument, stringsAsFactors = F, header = F)
      if(nrow(norm.table[which(norm.table[,1]==basename(tree.file.name)),])==1){
        normalisation.constants <- as.numeric(norm.table[which(norm.table[,1]==basename(tree.file.name)),2])
      } else if (nrow(norm.table[which(norm.table[,1]==basename(tree.file.name)),])==0){
        warning(paste("No normalisation constant given for tree file name ", basename(tree.file.name), " in file ",normalisation.argument,"; normalised distances will not be given", sep=""))
        normalisation.constants <- 1
      } else {
        warning(paste("Two or more entries for normalisation constant for tree file name ",basename(tree.file.name)," in file ",normalisation.argument,"; taking the first", sep=""))
        normalisation.constants <- as.numeric(norm.table[which(norm.table[,1]==basename(tree.file.name))[1],2])
      }
    }
  }
} else {
  # Assume we are dealing with a group of files
  
  tree.file.names		<- list.files.mod(dirname(tree.file.name), pattern=paste(basename(tree.file.name),'.*\\',tree.fe,'$',sep=''), full.names=TRUE)
  splits.file.names		<- list.files.mod(dirname(splits.file.name), pattern=paste(basename(splits.file.name),'.*\\',csv.fe,'$',sep=''), full.names=TRUE)
  
  if(length(tree.file.names)!=length(splits.file.names)){
    cat("Numbers of tree files and subgraph files differ\n")
    quit(save="no", status=1)
  }
  
  if(length(tree.file.names)==0){
    cat("No input trees found.\n Quitting.\n")
    quit(save="no", status=1)
  }
  
  suffixes <- substr(tree.file.names, nchar(tree.file.name) + 1, nchar(tree.file.names) - nchar(tree.fe))
  suffixes <- suffixes[order(suffixes)]
  
  tree.file.names <- tree.file.names[order(suffixes)]
  splits.file.names <- splits.file.names[order(suffixes)]
  
  output.file.names <- paste(output.root, "_classification_", suffixes, csv.fe, sep="")
  if(do.collapsed){
    collapsed.file.names <- paste(output.root, "_collapsedTree_",suffixes, csv.fe, sep="")
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
}


for(i in 1:length(tree.file.names)){	
  if(do.collapsed){
    dddf <- classify(tree.file.names[i], splits.file.names[i], normalisation.constants[i], collapsed.file.names[i])
  } else {
    dddf <- classify(tree.file.names[i], splits.file.names[i], normalisation.constants[i])
  }
  if (verbose) cat("Writing to file",output.file.names[i],"...\n")
  write.table(dddf, file = output.file.names[i], sep = ",", row.names = FALSE, col.names = TRUE, quote=F)	
} 
