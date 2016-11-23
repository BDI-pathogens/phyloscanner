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
  tree.file.names <- "ProcessedTree_r_test_r.tree"
  splits.file.names <- "Subtrees_r_test_r.csv"
  
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
  romero.severson <- T
  
  
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
#	OR changelog:
#	Put all calculations per tree into function, so that this script can  
#	also loop over several windows. This allows to combine the SummaryStats.R and LikelyTransmissions.R into
#	a single UNIX script

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
  
  assocs <- lapply(assocs, function(x) replace(x,is.na(x),"none"))
  
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
	
	cat("Testing pairs\n")
	
	count <- 0
	direct.descendant.matrix <- matrix(NA, length(patients.included), length(patients.included))
	distance.matrix <- matrix(NA, length(patients.included), length(patients.included))
	branches.longer.than.root.matrix <- matrix(NA, length(patients.included), length(patients.included))
	details.matrix <- matrix(NA, length(patients.included), length(patients.included))
	
	for(pat.1 in seq(1, length(patients.included))){
	  for(pat.2 in  seq(1, length(patients.included))){
	    if (pat.1 < pat.2) {
	      
	      count <- count + 1
	      
	      pat.1.id <- patients.included[pat.1]
	      pat.2.id <- patients.included[pat.2]
	      
	      all.nodes <-  c(splits.for.patients[[pat.1.id]], splits.for.patients[[pat.2.id]])
	      
	      #we want that they form a perfect blob with no intervening nodes for any other patients
	      OK <- TRUE
	      
	      for(node.1 in seq(1, length(all.nodes))){
	        for(node.2 in seq(1, length(all.nodes))){
	          if(node.1 < node.2){
	            node.1.id <- all.nodes[node.1]
	            node.2.id <- all.nodes[node.2]
	            path <- get.tt.path(tt, node.1.id, node.2.id)
	            for(node in path){
	              if(!startsWith(node, "none")){
	                if(!(patients.for.splits[[node]] %in% c(pat.1.id, pat.2.id))){
	                  OK <- FALSE
	                  break
	                }
	              }
	            }
	          }
	          if(!OK){
	            break
	          }
	        }
	        if(!OK){
	          break
	        }
	      }
	      
	      # then we want to see if any node from one patient interrupts a path between two nodes from the other
	      if(OK){
	        entangled <- F
	        
	        if(length(splits.for.patients[[pat.1.id]])>1){
	          combns.1 <- t(combn(splits.for.patients[[pat.1.id]], 2))
	          for(comb in seq(1, nrow(combns.1))){
	            path <- get.tt.path(tt, combns.1[comb, 1], combns.1[comb, 2])
	            for(node in path){
	              if(!startsWith(node, "none")){
	                if(patients.for.splits[[node]]==pat.2.id){
	                  entangled <- T
	                  break
	                }
	              }
	            }
	            if(entangled){
	              break
	            }
	          }
	        }
	        if(!entangled & length(splits.for.patients[[pat.2.id]])>1){
	          combns.2 <- t(combn(splits.for.patients[[pat.2.id]], 2))
	          for(comb in seq(1, nrow(combns.2))){
	            path <- get.tt.path(tt, combns.2[comb, 1], combns.2[comb, 2])
	            for(node in path){
	              if(!startsWith(node, "none")){
	                if(patients.for.splits[[node]]==pat.1.id){
	                  entangled <- T
	                  break
	                }
	              }
	            }
	            if(entangled){
	              break
	            }
	          }
	        }
	        rel.determined <- F
	        
	        if(!entangled){
	          for(pat.1.splt in splits.for.patients[[pat.1.id]]){
	            ancestors <- get.tt.ancestors(tt, pat.1.splt)
	            if(length(intersect(ancestors, splits.for.patients[[pat.2.id]])) > 0){
	              direct.descendant.matrix[pat.1, pat.2] <- "desc"
	              if(length(all.nodes)==2){
	                length <- 0
	                current.node <- tt$root.nos[which(tt$unique.splits==pat.1.id)]
	                while(assocs[[current.node]]!=pat.2.id){
	                  length <- length + get.edge.length(tree, current.node)
	                  current.node <- Ancestors(tree, current.node, type="parent")
	                }
	                distance.matrix[pat.1, pat.2] <- length
	              } else {
	                distance.matrix[pat.1, pat.2] <- NA
	              }
	              rel.determined <- T
	              break
	            }
	          }
	          if(!rel.determined){
	            for(pat.2.splt in splits.for.patients[[pat.2.id]]){
	              ancestors <- get.tt.ancestors(tt, pat.2.splt)
	              if(length(intersect(ancestors, splits.for.patients[[pat.1.id]]))>0){
	                direct.descendant.matrix[pat.1, pat.2] <- "anc"
	                if(length(all.nodes)==2){
	                  length <- 0
	                  current.node <- tt$root.nos[which(tt$unique.splits==pat.2.id)]
	                  while(assocs[[current.node]]!=pat.1.id){
	                    length <- length + get.edge.length(tree, current.node)
	                    current.node <- Ancestors(tree, current.node, type="parent")
	                  }
	                  distance.matrix[pat.1, pat.2] <- length
	                } else {
	                  distance.matrix[pat.1, pat.2] <- NA
	                }
	                rel.determined <- T
	                break
	              }
	            }
	          }
	          if(!rel.determined){
	            subtree <- extract.subtrees.for.patients(tree, c(pat.1.id, pat.2.id), splits)
	            max.length.prop <- 0
	            details.string <- ""
	            first <- T
	            for(node in all.nodes){
	              length.prop <- prop.internal.longer.than.root(subtree, node, splits)
	              if(length.prop > max.length.prop){
	                max.length.prop <- length.prop
	              }
	              if(first){
	                first <- F
	              } else {
	                details.string <- paste(details.string, ";", sep="")
	              }
	              details.string <- paste(details.string,node,":",length.prop, sep="")
	            }
	            branches.longer.than.root.matrix[pat.1, pat.2] <- max.length.prop
	            details.matrix[pat.1, pat.2] <- details.string
	            
	            if(length(all.nodes)==2){
	              root.1 <- tt$root.nos[which(tt$unique.splits==pat.1.id)]
	              root.2 <- tt$root.nos[which(tt$unique.splits==pat.2.id)]
	              ancestor.node.1 <- Ancestors(tree, root.1, type="parent")
	              if(length(ancestor.node.1) == 0){
	                ancestor.node.1 <- "root"
	              }
	              ancestor.node.2 <- Ancestors(tree, root.2, type="parent")
	              if(length(ancestor.node.2) == 0){
	                ancestor.node.2 <- "root"
	              }
	              if(ancestor.node.1 == ancestor.node.2){
	                direct.descendant.matrix[pat.1, pat.2] <- "cher"
	              } else {
	                direct.descendant.matrix[pat.1, pat.2] <- "unint"
	              }
	              distance.matrix[pat.1, pat.2] <- dist.nodes(tree)[root.1, root.2]
	            } else {
	              direct.descendant.matrix[pat.1, pat.2] <- "unint"
	              distance.matrix[pat.1, pat.2] <- NA
	            }
	          }
	        } else {
	          current.node <- all.nodes[1]
	          for(node in all.nodes[2:length(all.nodes)]){
	            current.node <- get.tt.mrca(tt, current.node, node)
	          }
	          if(startsWith(current.node, "none")){
	            direct.descendant.matrix[pat.1, pat.2] <- "trueInt"
	          } else if(patients.for.splits[[current.node]]==pat.1.id){
	            direct.descendant.matrix[pat.1, pat.2] <- "intAnc"
	          } else if(patients.for.splits[[current.node]]==pat.2.id){
	            direct.descendant.matrix[pat.1, pat.2] <- "intDesc"
	          } else {
	            stop("Classification failure (email the authors)")
	          }
	          distance.matrix[pat.1, pat.2] <- NA
	        }
	      }
	      
	      if (count %% 100 == 0) {
	        cat(
	          paste(
	            "Done ",format(count, scientific = F)," of ",format(total.pairs, scientific = F)," pairwise calculations\n", sep = ""
	          )
	        )
	      }
	    }
	  }
	}
	
	direct.descendant.table <- as.table(direct.descendant.matrix)
	branches.longer.than.root.table <- as.table(branches.longer.than.root.matrix)
	distance.table <- as.table(distance.matrix)
	details.table <- as.table(details.matrix)
	
	colnames(direct.descendant.table) <- patients.included
	rownames(direct.descendant.table) <- patients.included
	
	colnames(distance.table) <- patients.included
	rownames(distance.table) <- patients.included
	
	colnames(branches.longer.than.root.table) <- patients.included
	rownames(branches.longer.than.root.table) <- patients.included
	
	colnames(details.table) <- patients.included
	rownames(details.table) <- patients.included
	
	dddf <- as.data.frame(direct.descendant.table)
	
	keep <- complete.cases(dddf)
	
	dddf <- dddf[keep,]
	
	dmf <- as.data.frame(distance.table)
	dmf <- dmf[keep,]
	
	blmf <- as.data.frame(branches.longer.than.root.table)
	blmf <- blmf[keep,]
	
	dtmf <- as.data.frame(details.table)
  dtmf <- dtmf[keep,]
		
	dddf <- cbind(dddf, dmf[,3], blmf[,3], dtmf[,3])
	
	colnames(dddf) <- c("Patient_1", "Patient_2", "Relationship", "Distance", "Prop.internal.branches.longer.than.root.branch", "details")
	dddf
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