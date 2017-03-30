#!/usr/bin/env Rscript

command.line <- T

list.of.packages <- c("phangorn", "argparse", "phytools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

if(!("ggtree" %in% installed.packages()[,"Package"])){
  source("https://bioconductor.org/biocLite.R")
  biocLite("ggtree")
}

if(command.line){
  suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
  
  arg_parser = ArgumentParser(description="Identify the likely direction of transmission between all patients present in a phylogeny.")
  arg_parser$add_argument("-c", "--collapsedTree", action="store", help="If present, the collapsed tree (in which all adjacent nodes with the same assignment are collapsed to one) is output as a .csv file to the path specified.")
  arg_parser$add_argument("-n", "--branchLengthNormalisation", action="store", help="If present and a number, a normalising constant for all branch lengths in the tree or trees. If present and a file, the path to a .csv file with two columns: tree file name and normalising constant")
  arg_parser$add_argument("treeFileName", action="store", help="Tree file name. Alternatively, a base name that identifies a group of tree file names. Tree files are assumed to end in '.tree'. This must be the 'processed' tree produced by SplitTreesToSubtrees.R.")
  arg_parser$add_argument("splitsFileName", action="store", help="Splits file name. Alternatively, a base name that identifies a group of split file.")
  arg_parser$add_argument("outputFileName", action="store", help="Output file name (.csv format)")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
  args <- arg_parser$parse_args()
  
  tree.file.names <- args$treeFileName
  splits.file.names <- args$splitsFileName
  collapsed.file.names <- args$collapsedTree
  script.dir <- args$scriptdir
  output.name <- args$outputFileName
  has.normalisation <- !is.null(args$branchLengthNormalisation)
  normalisation.argument <- args$branchLengthNormalisation
  
} else {
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools"
  
  # BEEHIVE example
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20161013/")
  tree.file.names <- "ProcessedTrees_r/ProcessedTree_r_run20161013_RS_inWindow_6050_to_6400.tree"
  splits.file.names <- "SubtreeFiles_r/Subtrees_r_run20161013_RS_inWindow_6050_to_6400.csv"
  
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
  
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/Thai MRSA 6/Matthew/refinement")
  tree.file.names <- "ProcessedTree_s_mrsa_k10_bp_yetanother.tree"
  splits.file.names <- "Subtrees_s_mrsa_k10_bp_yetanother.csv"
  output.name <- "LT_s_mrsa_k10_yetanother.csv"
  collapsed.file.names <- "collapsed_s_mrsa_k10_yetanother.csv"
  
  
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
#
#	define internal functions
#
classify <- function(tree.file.name, splits.file.name, normalisation.constant = 1, tt.file.name = NULL) {	
  cat("Reading tree file ", tree.file.name, "...\n", sep = "")
  
  pseudo.beast.import <- read.beast(tree.file.name)
  tree <- attr(pseudo.beast.import, "phylo")
  
  cat("Reading splits file",splits.file.name,"...\n")
  
  splits <- read.csv(splits.file.name, stringsAsFactors = F)
  
  colnames(splits) <- c("patient", "subgraph", "split")
  
  cat("Collecting tips for each patient...\n")
  
  patients <- unique(splits$patient)
  patient.tips <- lapply(patients, function(x) which(tree$tip.label %in% splits[which(splits$patient==x), "tip.names"]))
  names(patient.tips) <- patients
  all.splits <- unique(splits$subgraph)
  
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
  
  splits.for.patients <- lapply(patients, function(x) unique(splits$subgraph[which(splits$patient==x)] ))
  names(splits.for.patients) <- patients
  
  patients.for.splits <- lapply(all.splits, function(x) unique(splits$patient[which(splits$subgraph==x)] ))
  names(patients.for.splits) <- all.splits
  
  patients.included <- patients
  
  total.pairs <- (length(patients.included) ^ 2 - length(patients.included))/2
  
  cat("Collapsing subtrees...\n")
  
  tt <- output.trans.tree(tree, assocs, tt.file.name)
  
  cat("Identifying pairs of unblocked splits...\n")
  
  adjacent <- subtrees.unblocked(tt, all.splits)
  
  cat("Calculating pairwise distances between splits...\n")
  
  split.distances <- tryCatch(
    all.subtree.distances(tree, tt, all.splits, assocs), warning=function(w){return(NULL)}, error=function(e){return(NULL)})
  
  if(is.null(split.distances)){
    split.distances <- all.subtree.distances(tree, tt, all.splits, assocs, TRUE)
  }

  cat("Testing pairs...\n")
  
  count <- 0
  adjacency.matrix <- matrix(NA, length(patients.included), length(patients.included))
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
        
        # we want that they form a perfect blob with no intervening nodes for any other patients
        OK <- check.adjacency(tt, c(pat.1.id, pat.2.id), splits.for.patients)		
        
        adjacency.matrix[pat.1, pat.2] <- OK

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
            if(adjacent[node.1, node.2]){
              pairwise.distances <- c(pairwise.distances, split.distances[node.1, node.2])
            }
          }
        }
        
        min.distance.matrix[pat.1, pat.2] <- min(pairwise.distances)
        
        
        if (count %% 100 == 0) {
          cat(paste("Done ",format(count, scientific = F)," of ",format(total.pairs, scientific = F)," pairwise calculations\n", sep = ""))
        }
      }
    }
  }
  

  normalised.distance.matrix<- min.distance.matrix/normalisation.constant
  
  adjacency.table <- as.table(adjacency.matrix)
  dir.12.table <- as.table(dir.12.matrix)
  dir.21.table <- as.table(dir.21.matrix)
  nodes.1.table <- as.table(nodes.1.matrix)
  nodes.2.table <- as.table(nodes.2.matrix)
  path.table <- as.table(path.matrix)
  min.distance.table <- as.table(min.distance.matrix)
  
  if(normalisation.constant!=1){
    normalised.distance.table <- as.table(normalised.distance.matrix)
  }
  
  colnames(adjacency.table) <- patients.included
  rownames(adjacency.table) <- patients.included
  
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
  
  if(!is.null(normalisation.constant)){
    colnames(normalised.distance.table) <- patients.included
    rownames(normalised.distance.table) <- patients.included
  }
  
  
  adf <- as.data.frame(adjacency.table)
  
  keep <- complete.cases(adf)
  
  adf <- adf[keep,]
  
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
    nddf <- mddf[keep,]
  }
  
  adf <- cbind(adf, d12df[,3], d21df[,3], n1df[,3], n2df[,3], pdf[,3], mddf[,3])
  
  column.names <- c("Patient_1", "Patient_2", "adjacent", "paths12", "paths21", "nodes1", "nodes2", "path.classification", "min.distance.between.subtrees")
  
  if(!is.null(normalisation.constant)){
    adf <- cbind(adf, nddf[,3])
    column.names <- c(column.names, "normalised.mean.distance.between.subtrees")
  }
  
  colnames(adf) <- column.names
  adf
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
  
  normalisation.constant <- NULL
  
  if(has.normalisation){
    normalisation.constant <- suppressWarnings(as.numeric(normalisation.argument))
    
    if(is.na(normalisation.constant)){
      norm.table <- read.csv(normalisation.argument, stringsAsFactors = F, header = F)
      if(nrow(norm.table[which(norm.table[,1]==tree.file.name),])==1){
        normalisation.constant <- as.numeric(norm.table[which(norm.table[,1]==tree.file.name),2])
      } else if (nrow(norm.table[which(norm.table[,1]==tree.file.name),])==0){
        warning(paste("No normalisation constant given for tree file name ",tree.file.name, " in file ",normalisation.argument,"; normalised distances will not be given", sep=""))
        normalisation.constant <- NULL
      } else {
        warning(paste("Two or more entries for normalisation constant for tree file name ",tree.file.name, " in file ",normalisation.argument,"; taking the first", sep=""))
        normalisation.constant <- as.numeric(norm.table[which(norm.table[,1]==tree.file.name)[1],2])
      }
    }
  }
  
  dddf				<- classify(tree.file.name, splits.file.name, normalisation.constant = normalisation.constant, tt.file.name = collapsed.file.names)
  cat("Write to file",output.name,"...\n")
  write.table(dddf, file = output.name, sep = ",", row.names = FALSE, col.names = TRUE, quote=F)	
} else {
  #	if 'tree.file.names' is not tree, process multiple trees
  tree.file.names		<- sort(list.files(dirname(tree.file.names), pattern=paste(basename(tree.file.names),'.*\\.tree$',sep=''), full.names=TRUE))
  splits.file.names	<- sort(list.files(dirname(splits.file.names), pattern=paste(basename(splits.file.names),'.*\\.csv$',sep=''), full.names=TRUE))	
  stopifnot(length(tree.file.names)==length(splits.file.names))
  
  normalisation.constant <- 1
  
  if(has.normalisation){
    normalisation.constant <- suppressWarnings(as.numeric(normalisation.argument))
    
    if(is.na(normalisation.constant)){
      norm.table <- read.csv(normalisation.argument, stringsAsFactors = F, header = F)
    }
  }
  
  
  for(tree.i in seq_along(tree.file.names)){
    tree.file.name		<- tree.file.names[tree.i]
    splits.file.name	<- splits.file.names[tree.i]
    output.name			<- gsub('\\.tree','_LikelyTransmissions.csv', tree.file.name)		
    collapsed.file.name	<- gsub('\\.tree','_collapsed.csv', tree.file.name)

    if(has.normalisation){
      if(has.normalisation & is.na(normalisation.constant)) {
        if(nrow(norm.table[which(norm.table[,1]==basename(tree.file.name)),])==1){
          normalisation.constant <- as.numeric(norm.table[which(norm.table[,1]==basename(tree.file.name)),2])
        } else if (nrow(norm.table[which(norm.table[,1]==basename(tree.file.name)),])==0){
          warning(paste("No normalisation constant given for tree file name ",basename(tree.file.name), " in file ",basename(normalisation.argument),"; normalised distances will not be given", sep=""))
          normalisation.constant <- 1
        } else {
          warning(paste("Two or more entries for normalisation constant for tree file name ",basename(tree.file.name), " in file ",basename(normalisation.argument),"; taking the first", sep=""))
          normalisation.constant <- as.numeric(norm.table[which(norm.table[,1]==basename(tree.file.name))[1],2])
        }
      }
    }
    
    tryCatch({
      dddf	<- classify(tree.file.name, splits.file.name, normalisation.constant = normalisation.constant, collapsed.file.name)
      cat("Write to file",output.name,"...\n")
      write.table(dddf, file = output.name, sep = ",", row.names = FALSE, col.names = TRUE, quote=F)			
    },
    error=function(e){ print(e$message); cat("No output written to file.\n") })    
  }
}