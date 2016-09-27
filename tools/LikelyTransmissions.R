#!/usr/bin/env Rscript

command.line <- T

list.of.packages <- c("phangorn", "argparse", "phytools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

if(command.line){
  library(argparse, quietly=TRUE, warn.conflicts=FALSE)
  
  arg_parser = ArgumentParser(description="Identify the likely direction of transmission between all patients present in a phylogeny.")
  arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                          help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name. Necessary only if -s is specified.")
  arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
  arg_parser$add_argument("-z", "--zeroLengthTipsCount", action="store_true", default=FALSE, help="If present, a zero length terminal branch associates the patient at the tip with its parent node, interrupting any inferred transmission pathway between another pair of hosts that goes through the node.")
  arg_parser$add_argument("-t", "--dualInfectionThreshold", type="double", help="Length threshold at which a branch of the subtree constructed just from reads from a single patient is considered long enough to indicate a dual infection (such a patient will be ignored, at present). If absent, keep all patients.")
  arg_parser$add_argument("-s", "--romeroSeverson", default=FALSE, action="store_true", help="If present, the Romero-Severson classification will be repeated on tree to annotate internal nodes beyond the MRCAs of each split. Recommended if R-S was used to determine splits.")
  arg_parser$add_argument("-c", "--collapsedTree", action="store", help="If present, the collapsed tree (in which all adjacent nodes with the same assignment are collapsed to one) is output as a .csv file to the path specified.")
  arg_parser$add_argument("treeFileName", action="store", help="Tree file name. Alternatively, a base name that identifies a group of tree file names can be specified. Tree files are assumed to end in '.tree'.")
  arg_parser$add_argument("splitsFileName", action="store", help="Splits file name. Alternatively, a base name that identifies a group of split file names can be specified. Split files are assumed to end in '_subtree_[a-z].csv'.")
  arg_parser$add_argument("outputFileName", action="store", help="Output file name (.csv format)")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
  
  args <- arg_parser$parse_args()
  
  tree.file.names <- args$treeFileName
  splits.file.names <- args$splitsFileName
  collapsed.file.names <- args$collapsedTree
  script.dir <- args$scriptdir
  output.name <- args$outputFileName
  root.name <- args$outgroupName
  tip.regex <- args$tipRegex
  split.threshold <- args$dualInfectionThreshold
  if(is.null(split.threshold)){
    split.threshold <- Inf
  }
  zero.length.tips.count <- args$zeroLengthTipsCount
  romero.severson <- args$romeroSeverson
  
} else {
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/")
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools"
  tree.file.names <- "ptyr78_trees_newick/ptyr78_InWindow_850_to_1099.tree"
  splits.file.names <- "ptyr78_subtrees_r_csv/Subtrees_r_ptyr78_InWindow_850_to_1099.csv"
  output.name <- "hi.csv"
  root.name <- "REF_CPX_AF460972"
  split.threshold <- 0.08
  zero.length.tips.count <- F
  romero.severson <- T
  tip.regex <- "^(.*)_read_([0-9]+)_count_([0-9]+)$"
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

library(phytools, quietly=TRUE, warn.conflicts=FALSE)
library(phangorn, quietly=TRUE, warn.conflicts=FALSE)
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
	
	tree <-read.tree(tree.file.name)
	tree <- root(tree, root.name)
	
	cat("Making tree multiphyletic...\n")
	
	tree <- di2multi(tree, tol = 1E-5)
	
	cat("Reading splits file",splits.file.name,"...\n")
	
	splits <- read.csv(splits.file.name, stringsAsFactors = F)
	
	cat("Collecting tips for each patient...\n")
	
	patients <- unique(splits$orig.patients)
	patient.tips <- lapply(patients, function(x) which(tree$tip.label %in% splits[which(splits$orig.patients==x), "tip.names"]))
	names(patient.tips) <- patients
	all.splits <- unique(splits$patient.splits)
	if(romero.severson){
	  
	  # The reason the procedure is repeated is that the RS classification for internal nodes cannot simply be determined 
	  # from the set of tips in each subtree (unlike the conservative classification). While it would be possible to export 
	  # annotations in SplitPatientsToSubtrees.R, the concern is over anything that might renumber the tree nodes.
	  
	  # Really the only reason to load the splits in is to avoid needing the blacklist yet again. Maybe rethink this.
	  
	  patient.mrcas <- lapply(patient.tips, function(node) mrca.phylo.or.unique.tip(tree, node, zero.length.tips.count))
	  
	  # custom blacklist
	  
	  blacklist <- which(!(tree$tip.label %in% splits$tip.names))
	  
	  classification <- split.and.annotate(tree, patients, patient.tips, patient.mrcas, blacklist, tip.regex, "r")

	  split.ids <- classification$split.patients
	  
	  was.split <- unique(patients[!(patients %in% split.ids)])
	  
	  split.tips <- classification$split.tips
	  assocs <- classification$assocs
	  
	} else {
	  cat("Collecting tips for each split...\n")
	  
	  split.tips <- lapply(all.splits, function(x) which(tree$tip.label %in% splits[which(splits$patient.splits==x), "tip.names"]))
	  names(split.tips) <- all.splits
	  
	  tip.assocs.splits <- lapply(tree$tip.label, function(x) if(x %in% splits$tip.names) return(splits$patient.splits[which(splits$tip.names==x)]) else return("*"))
	  
	  split.mrcas <- lapply(split.tips, function(node) mrca.phylo.or.unique.tip(tree, node, zero.length.tips.count))
	  
	  node.assocs <- annotate.internal(tree, all.splits, split.tips, split.mrcas)
	  
	  was.split <- unique(splits$orig.patients[which(splits$orig.patients!=splits$patient.splits)])
	  
	  assocs <- node.assocs$details
	  
	}
	
	splits.for.patients <- lapply(patients, function(x) unique(splits$patient.splits[which(splits$orig.patients==x)] ))
	names(splits.for.patients) <- patients
	
	patients.for.splits <- lapply(all.splits, function(x) unique(splits$orig.patients[which(splits$patient.splits==x)] ))
	names(patients.for.splits) <- all.splits
	
	# get rid of what look like dual infections
	
	ignore.because.split <- vector()
	
	for(split.patient in was.split){
	  tips <- patient.tips[[split.patient]]
	  subtree <- drop.tip(phy = tree, tip = tree$tip.label[-patient.tips[[split.patient]]])
	  max.length <- max(subtree$edge.length)
	  if(max.length > split.threshold){
	    ignore.because.split <- c(ignore.because.split, split.patient)
	  }
	}
	
	patients.included <- setdiff(patients, ignore.because.split)
	
	total.pairs <- (length(patients.included) ^ 2 - length(patients.included))/2
	
	cat("Collapsing subtrees...\n")
	
	tt <- output.trans.tree(tree, assocs, tt.file.name)
	
	cat("Testing pairs\n")
	
	count <- 0
	direct.descendant.matrix <- matrix(NA, length(patients.included), length(patients.included))
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
	            if(length(intersect(ancestors, splits.for.patients[[pat.2.id]]))>0){
	              direct.descendant.matrix[pat.1, pat.2] <- "desc"
	              rel.determined <- T
	              break
	            }
	          }
	          if(!rel.determined){
	            for(pat.2.splt in splits.for.patients[[pat.2.id]]){
	              ancestors <- get.tt.ancestors(tt, pat.2.splt)
	              if(length(intersect(ancestors, splits.for.patients[[pat.1.id]]))>0){
	                direct.descendant.matrix[pat.1, pat.2] <- "anc"
	                rel.determined <- T
	                break
	              }
	            }
	          }
	          if(!rel.determined){
	            direct.descendant.matrix[pat.1, pat.2] <- "sib"
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
	
	colnames(direct.descendant.table) <- patients.included
	rownames(direct.descendant.table) <- patients.included
	
	
	dddf <- as.data.frame(direct.descendant.table)
	dddf <- dddf[complete.cases(dddf),]
	
	colnames(dddf) <- c("Patient_1", "Patient_2", "Relationship")
	dddf
}

#
#	start script
#

#	check if 'tree.file.names' is tree
options(show.error.messages = FALSE)		
can.read.tree		<- try(suppressWarnings(read.tree(tree.file.names)))
options(show.error.messages = TRUE)	

if(!inherits(can.read.tree, "try-error"))
{
	#	if 'tree.file.names' is tree, process just one tree
	tree.file.name		<- tree.file.names[1]
	splits.file.name	<- splits.file.names[1]
	dddf				<- likely.transmissions(tree.file.name, splits.file.name, tip.regex, split.threshold, romero.severson, zero.length.tips.count, collapsed.file.names)
	cat("Write to file",output.name,"...\n")
	write.table(dddf, file = output.name, sep = ",", row.names = FALSE, col.names = TRUE, quote=F)	
}

if(inherits(can.read.tree, "try-error"))
{
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