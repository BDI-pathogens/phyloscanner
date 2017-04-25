list.of.packages <- c("argparse", "ggplot2", "ff", "ggtree", "phangorn", "ape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  cat("Please run PackageInstall.R to continue\n")
  quit(save='no')
}

command.line <- T

suppressMessages(library(phangorn, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ggplot2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ggtree, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(data.table, quietly=TRUE, warn.conflicts=FALSE))

tree.fe <- ".tree"
csv.fe <- ".csv"

if(command.line){
  suppressMessages(require(argparse, quietly=TRUE, warn.conflicts=FALSE))

  arg_parser = ArgumentParser(description="Annotate the internal nodes of the phylogeny, by an ancestral state reconstruction scheme, such that the tree can subsequently be partitioned into connected subgraphs, one per split. Input and output file name arguments are either a single file, or the root name for a group of files, in which case all files matching that root will be processed, and matching output files generated.")  
  arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
  arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
  arg_parser$add_argument("-s", "--splitsRule", action="store", default="r", help="The rules by which the sets of patients are split into groups in order to ensure that all groups can be members of connected subgraphs without causing conflicts. Currently available: s=Sankhoff (slow, rigorous), r=Romero-Severson (quick, less rigorous with >2 patients).")
  arg_parser$add_argument("-b", "--blacklist", action="store", help="A .csv file listing tips to ignore. Alternatively, a base name that identifies blacklist files. Blacklist files are assumed to end in .csv.")
  arg_parser$add_argument("-k", "--kParam", action="store", default=0, help="The k parameter in the cost matrix for Sankhoff reconstruction (see documentation)")
  arg_parser$add_argument("-P", "--pruneBlacklist", action="store_true", help="If present, all blacklisted and references tips (except the outgroup) are pruned away before starting.")
  arg_parser$add_argument("-p", "--proximityThreshold", action="store", default=0)
  arg_parser$add_argument("-R", "--readCountsMatterOnZeroBranches", default = FALSE, action="store_true")
  arg_parser$add_argument("-m", "--multifurcationThreshold", help="If specified, short branches in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. If a number, this number will be used as the threshold. If 'g', it will be guessed from the branch lengths (use this only if you've checked by eye that the tree does indeed have multifurcations).")
  arg_parser$add_argument("-n", "--branchLengthNormalisation", default = 1, action="store", help="If present and a number, a normalising constant for all branch lengths in the tree or trees. If present and a file, the path to a .csv file with two columns: tree file name and normalising constant")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
  arg_parser$add_argument("-OD", "--outputdir", action="store", help="Full path of the directory for output; if absent, current working directory")
  arg_parser$add_argument("-pw", "--pdfwidth", action="store", default=100, help="Width of tree pdf in inches.")
  arg_parser$add_argument("-ph", "--pdfrelheight", action="store", default=0.15, help="Relative height of tree pdf.")
  arg_parser$add_argument("-ff", "--useff", action="store_true", default=FALSE, help="Use ff to store parsimony reconstruction matrices. Use if you run out of memory.")
  arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")
  arg_parser$add_argument("inputFileName", action="store", help="The file or file root for the input tree in Newick format")
  arg_parser$add_argument("outputFileID", action="store", help="A string shared by all output file names.")

  args <- arg_parser$parse_args()
  
  script.dir <- args$scriptdir
  
  source(file.path(script.dir, "TreeUtilityFunctions.R"))
  source(file.path(script.dir, "ParsimonyReconstructionMethods.R"))
  source(file.path(script.dir, "CollapsedTreeMethods.R"))
  source(file.path(script.dir, "WriteAnnotatedTrees.R"))
  source(file.path(script.dir, "GeneralFunctions.R"))
  
  tree.file.name <- args$inputFile

  sankhoff.k <- as.numeric(args$kParam)
  sankhoff.p <- as.numeric(args$proximityThreshold)
  output.file.ID <- args$outputFileID
  blacklist.file.name <- args$blacklist
  root.name <- args$outgroupName
  read.counts.matter <- args$readCountsMatterOnZeroBranches
  verbose <- args$verbose
  if(!is.null(args$outputDir)){
    output.dir <- args$outputDir
  } else {
    output.dir <- getwd()
  }
  
  prune.blacklist <- args$pruneBlacklist
  
  use.m.thresh <- !is.null(args$multifurcationThreshold)
  if(use.m.thresh){
    if(args$multifurcationThreshold=="g"){
      m.thresh <- NA
    } else if(!is.na(as.numeric(args$multifurcationThreshold))){
      m.thresh <- as.numeric(args$multifurcationThreshold)
    } else {
      cat("Unknown argument for -m specified\n")
      quit(save="no")
    }
  }
  
  if(!use.m.thresh & read.counts.matter){
    cat("Please specify a multifurcation collapse threshold with -m to count reads on zero length branches in the reconstruction")
    quit(save="no")
  }
  
  if(is.null(root.name)){
    cat("No outgroup name given; will assume the tree is already rooted\n")
  }
  
  tip.regex <- args$tipRegex
  mode <- args$splitsRule
  pdf.hm <- as.numeric(args$pdfrelheight)
  pdf.w <- as.numeric(args$pdfwidth)
  useff = args$useff
  
  has.normalisation <- !is.null(args$branchLengthNormalisation)
  normalisation.argument <- args$branchLengthNormalisation
  
  if(useff){
    suppressMessages(library(ff, quietly=TRUE, warn.conflicts=FALSE))
  }
  

  if(!(mode %in% c("r", "s", "f"))){
    stop(paste("Unknown split classifier: ", mode, "\n", sep=""))
  }


} else {
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools/"
  pdf.w = 10
  pdf.hm = 0.2
  
  # # BEEHIVE example
  # 
  # setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517_clean/")
  # output.dir <- getwd()
  # tree.file.names <- "RAxML_bestTree.InWindow_1550_to_1900.tree"
  # blacklist.files <- "PatientBlacklist_InWindow_1550_to_1900.csv"
  # out.identifier <- "test_s"
  # mode <- "s"
  # root.name <- "C.BW.00.00BW07621.AF443088"
  # tip.regex <- "^(.*)-[0-9].*_read_([0-9]+)_count_([0-9]+)$"
  # sankhoff.k <- 25
  # break.ties.unsampled <- TRUE
  # 
  # Rakai example
  # 
  #   setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr3_trees_collapsed.zip")
  #   output.dir <- getwd()
  #   tree.file.names <- "/Users/Oliver/duke/tmp/pty_17-03-23-20-42-16/ptyr122_InWindow_3200_to_3449.tree"
  #   blacklist.files <- "/Users/Oliver/duke/tmp/pty_17-03-23-20-42-16/ptyr122_blacklistdwns_InWindow_3200_to_3449.csv"
  #   out.identifier <- "ptyr1_"
  #   root.name <- "REF_CPX_AF460972"
  #   tip.regex <- "^(.*)_read_([0-9]+)_count_([0-9]+)$"
  #   mode <- "s"
  #   zero.length.tips.count <- F
  #   sankhoff.k <- 20
  #   break.ties.unsampled <- F
  #  
  # MRSA example
  
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/Thai MRSA 6/Matthew/")
  output.dir <- "/Users/twoseventwo/Dropbox (Infectious Disease)/Thai MRSA 6/Matthew/refinement"
  tree.file.names <- "RAxML_bipartitions.ST239_no_bootstraps_T056corr.tree"
  blacklist.file.name <- NULL
  output.file.IDs <- "mrsa_k10_bp_test_finaliteration"
  tip.regex <- "^([ST][0-9][0-9][0-9])_[A-Z0-9]*_[A-Z][0-9][0-9]$"
  root.name <- "TW20"
  mode <- "f"
  sankhoff.k <- 10
  sankhoff.p <- 0.1
  useff  <- F
  
  setwd("/Users/twoseventwo/Documents/Croucher alignments/")
  output.dir <- "/Users/twoseventwo/Dropbox (Infectious Disease)/Thai MRSA 6/Matthew/refinement"
  tree.file.names <- "MBconsensus.tre"
  blacklist.file.name <- NULL
  output.file.IDs <- "test"
  tip.regex <- "^(ARI-[0-9][0-9][0-9][0-9]_[A-Z]*)_[0-9][0-9][0-9][0-9]_[0-9]_[0-9][0-9]?$"
  root.name <- NULL
  mode <- "s"
  sankhoff.k <- 10
  useff  <- F
  
  setwd("/Users/twoseventwo/Documents/Croucher alignments/")
  output.dir <- "/Users/twoseventwo/Dropbox (Infectious Disease)/Thai MRSA 6/Matthew/refinement"
  tree.file.names <- "MBconsensus.tre"
  blacklist.file.name <- NULL
  output.file.IDs <- "test"
  tip.regex <- "^(ARI-[0-9][0-9][0-9][0-9]_[A-Z]*)_[0-9][0-9][0-9][0-9]_[0-9]_[0-9][0-9]?$"
  root.name <- NULL
  mode <- "s"
  sankhoff.k <- 10
  useff  <- F
  
}


#
#	define internal functions
#
split.patients.to.subgraphs<- function(tree.file.name, normalisation.constant = 1, mode, 
                                       blacklist.file, root.name, tip.regex, sankhoff.k, 
                                       sankhoff.p, useff, count.reads){
  
  if (verbose) cat("SplitPatientsToSubgraphs.R run on: ", tree.file.name,", rules = ",mode,"\n", sep="")
  
  # Read, root and multifurcate the tree

  tree <- read.tree(tree.file.name)
  tree <- unroot(tree)
  
  if(use.m.thresh){
    if(is.na(m.thresh)){
      min.bl <- min(tree$edge.length)
      m.thresh <- min.bl*1.001
    } 
    tree <- di2multi(tree, tol = m.thresh)
  }
  
  if(!is.null(root.name)){
    tree <- root(tree, outgroup = which(tree$tip.label==root.name), resolve.root = T)
  }
  
  # Load the blacklist if it exists
  
  blacklist <- vector()
  
  if(!is.null(blacklist.file)){
    if(file.exists(blacklist.file)){
      if (verbose) cat("Reading blacklist file",blacklist.file,'\n')
      blacklisted.tips <- read.table(blacklist.file, sep=",", header=F, stringsAsFactors = F, col.names="read")
      if(nrow(blacklisted.tips)>0){
        blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
      }
    } else {
      cat("WARNING: File ",blacklist.file," does not exist; skipping.\n",sep="")
    }
  } 
  
  if(prune.blacklist){
    non.patient.tips <- which(is.na(sapply(tree$tip.label, function(name) patient.from.label(name, tip.regex))))
    outgroup <- which(tree$tip.label==root.name)
    tips.to.go <- c(non.patient.tips, blacklist)
    tips.to.go <- setdiff(tips.to.go, outgroup)
    
    tree <- drop.tip(tree, tips.to.go)
    
    # clear the blacklist
    
    blacklist <- vector()
  }
  
  orig.tree <- tree
  
  tree$edge.length <- tree$edge.length/normalisation.constant
  
  tip.labels <- tree$tip.label
  
  if (verbose) cat("Identifying tips with patients...\n")
  
  # Find patient IDs from each tip
  
  patient.ids <- sapply(tip.labels, function(x) patient.from.label(x, tip.regex))
  
  # Find the complete list of patients present in this tree minus blacklisting
  
  if(length(blacklist)>0){
    patients <- unique(na.omit(patient.ids[-blacklist]))
  } else {
    patients <- unique(na.omit(patient.ids))
  }
  
  patients <- patients[order(patients)]
  
  if(length(patients)==0){
    cat("No patient IDs detected, nothing to do.\n")
    quit(save="no")
  } 

  
  # patient.tips is a list of tips that belong to each patient
  
  patient.tips <- lapply(patients, function(x) setdiff(which(patient.ids==x), blacklist))
  names(patient.tips) <- patients
  
  # patient.mrcas is a list of MRCA nodes for each patient
  
  patient.mrcas <- lapply(patient.tips, function(node) mrca.phylo.or.unique.tip(tree, node))
  
  # Do the main function
  
  results <- split.and.annotate(tree, patients, patient.tips, patient.mrcas, blacklist, tip.regex, mode, sankhoff.k, sankhoff.p, root.name, m.thresh, useff = useff, count.reads, verbose)
  
  # Where to put the node shapes that display subgraph MRCAs
  
  node.shapes <- rep(FALSE, length(tree$tip.label) + tree$Nnode)
  node.shapes[unlist(results$first.nodes)] <- TRUE

  # For display purposes. "unsampled" nodes should have NA as annotations
  
  split.annotation <- rep(NA, length(tree$tip.label) + tree$Nnode)
  
  for(item in seq(1, length(results$assocs))){
    if(!is.null(results$assocs[[item]])){
      if(!(results$assocs[[item]] %in% c("*", "unsampled", "unsampled.rt", "unsampled.int"))) {
        split.annotation[item] <- results$assocs[[item]]
      }
    }
  }
  
  # This is the annotation for each node by patient
  
  patient.annotation <- sapply(split.annotation, function(x) unlist(strsplit(x, "-S"))[1] )
  names(patient.annotation) <- NULL
  patient.annotation <- factor(patient.annotation, levels = sample(levels(as.factor(patient.annotation))))
  
  # For annotation
  # We will return the original tree with the original branch lengths
  
  attr(orig.tree, 'SPLIT') <- factor(split.annotation)
  attr(orig.tree, 'INDIVIDUAL') <- patient.annotation
  attr(orig.tree, 'NODE_SHAPES') <- node.shapes
  
  rs.subgraphs <- data.table(subgraph=results$split.patients)
  rs.subgraphs <- rs.subgraphs[, list(patient= unlist(strsplit(subgraph, "-S"))[1], tip= tree$tip.label[ results$split.tips[[subgraph]] ]	), by='subgraph']
  rs.subgraphs <- as.data.frame(subset(rs.subgraphs, select=c(patient, subgraph, tip)))
  
  list(tree=orig.tree, rs.subgraphs=rs.subgraphs)
}
#
#	start script
#

if(file.exists(tree.file.name)){

  file.details <- list()
  
  file.name.list <- list()
  file.name.list$tree.input <- tree.file.name

  tree.file.names <- tree.file.name
  if(!is.null(blacklist.file.name)){
    file.name.list$blacklist.input <- blacklist.file.name
  }
  
  file.name.list$output.ID <- output.file.ID
  
  normalisation.constant <- suppressWarnings(as.numeric(normalisation.argument))
  
  if(is.na(normalisation.constant)){
    norm.table <- read.csv(normalisation.argument, stringsAsFactors = F, header = F)
    if(nrow(norm.table[which(norm.table[,1]==basename(tree.file.name)),])==1){
      normalisation.constant <- as.numeric(norm.table[which(norm.table[,1]==basename(tree.file.name)),2])
    } else if (nrow(norm.table[which(norm.table[,1]==basename(tree.file.name)),])==0){
      warning(paste("No normalisation constant given for tree file name ", basename(tree.file.name), " in file ",normalisation.argument,"; normalised distances will not be given", sep=""))
      normalisation.constant <- 1
    } else {
      warning(paste("Two or more entries for normalisation constant for tree file name ",basename(tree.file.name)," in file ",normalisation.argument,"; taking the first", sep=""))
      normalisation.constant <- as.numeric(norm.table[which(norm.table[,1]==basename(tree.file.name))[1],2])
    }
  }

  file.name.list$normalisation.constant <- normalisation.constant
  
  file.details[[tree.file.name]] <- file.name.list
} else {
  # Assume we are dealing with a group of files
  tree.file.names	<- list.files.mod(dirname(tree.file.name), pattern=paste(basename(tree.file.name),'.*\\',tree.fe,'$',sep=''), full.names=TRUE)

  if(length(tree.file.names)==0){
    cat("No input trees found.\n")
    quit(save="no")
  }
  
  suffixes <- substr(tree.file.names, nchar(tree.file.name) + 1, nchar(tree.file.names)-nchar(tree.fe))
  if(!is.null(blacklist.file.name)){  
    blacklist.file.names <- paste(blacklist.file.name, suffixes, csv.fe, sep="")
  }
  output.file.IDs <- paste(output.file.ID, suffixes, sep="")
  
  
  normalisation.constant <- suppressWarnings(as.numeric(normalisation.argument))
  
  if(is.na(normalisation.constant)){
    norm.table <- read.csv(normalisation.argument, stringsAsFactors = F, header = F)
    normalisation.constants <- sapply(basename(tree.file.names), function(x){
      if(x %in% norm.table[,1]){
        if(length(which(norm.table[,1]==x))>1){
          warning(paste("Two or more entries for normalisation constant for tree file name ",x, " in file ",normalisation.argument,"; taking the first", sep=""))
        }
        return(norm.table[which(norm.table[,1]==x)[1],2])
      } else {
        warning(paste("No normalisation constant for tree file ",x," found in file, ",normalisation.argument,"; this tree will not have branch lengths normalised.", sep=""))
        return(1)
      }
    })  
  } else {
    # it's a number to be used in every tree
    normalisation.constants <- rep(normalisation.constant, length(tree.file.names))
  }
  fn.df <- data.frame(row.names = suffixes, tree.input = tree.file.names, output.ID = output.file.IDs, normalisation.constant = normalisation.constants, stringsAsFactors = F)
  if(!is.null(blacklist.file.name)){
    fn.df$blacklist.input <- blacklist.file.names
  }
  file.details <- split(fn.df, rownames(fn.df))
}

for(i in file.details){
  #	if 'tree.file.names' is tree, process just one tree
  tree.file.name		<- i$tree.input	
  if(!is.null(blacklist.file.name)){
    blacklist.file	<- i$blacklist.input
  } else {
    blacklist.file <- NULL
  }
  output.string <- i$output.ID 
  
  tmp					<- split.patients.to.subgraphs(tree.file.name, i$normalisation.constant, mode, blacklist.file, root.name, tip.regex, sankhoff.k, sankhoff.p, useff, read.counts.matter)
  tree				<- tmp[['tree']]	
  rs.subgraphs			<- tmp[['rs.subgraphs']]
  #
  #	write rda file (potentially before plotting fails so we can recover)
  #
  tmp 				<- file.path(output.dir,paste('subgraphs_',mode,'_',output.string,'.rda',sep='')) 
  if (verbose) cat("Writing output to file",tmp,"...\n")
  save(rs.subgraphs, tree, file=tmp)		
  #
  #	plot tree
  #
  if (verbose) cat("Drawing tree...\n")
  tree.display 		<- ggtree(tree, aes(color=INDIVIDUAL)) +
    geom_point2(shape = 16, size=3, aes(subset=NODE_SHAPES)) +
    scale_fill_hue(na.value = "black") +
    scale_color_hue(na.value = "black") +
    theme(legend.position="none") +
    geom_tiplab(aes(col=INDIVIDUAL)) + 
    geom_treescale(width=0.01, y=-5, offset=1.5)
  
  x.max <- ggplot_build(tree.display)$layout$panel_ranges[[1]]$x.range[2]
  
  tree.display <- tree.display + ggplot2::xlim(0, 1.1*x.max)
  tree.display	

  tmp	<- file.path(output.dir,paste('Tree_',output.string,'.pdf',sep=''))
  if (verbose) cat("Plot to file",tmp,"...\n")
  ggsave(tmp, device="pdf", height = pdf.hm*length(tree$tip.label), width = pdf.w, limitsize = F)
  
  tmp <- file.path(output.dir, paste('ProcessedTree_',mode,'_',output.string, tree.fe,sep=''))
  if (verbose) cat("Writing rerooted, multifurcating, annotated tree to file",tmp,"...\n")	
  write.ann.nexus(tree, file=tmp, annotations = c("INDIVIDUAL", "SPLIT"))
  
  #
  #	write csv file
  #
  
  tmp 				<- file.path(output.dir, paste('subgraphs_',mode,'_',output.string,csv.fe,sep=''))
  if (verbose) cat("Writing output to file",tmp,"...\n")	
  write.csv(rs.subgraphs, file=tmp, row.names = F, quote=F)	
  
}
