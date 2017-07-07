#!/usr/bin/env Rscript

list.of.packages <- c("argparse", "ggplot2", "ff", "ggtree", "phangorn", "ape", "kimisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  cat("Please run package_install.R to continue\n")
  quit(save="no", status=1)
}

suppressMessages(library(phangorn, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ggplot2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ggtree, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(data.table, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(kimisc, quietly=TRUE, warn.conflicts=FALSE))

arg_parser = ArgumentParser(description="Annotate the internal nodes of the phylogeny, by an ancestral state reconstruction scheme, such that the tree can subsequently be partitioned into connected subgraphs, one per split. Input and output file name arguments are either a single file, or the root name for a group of files, in which case all files matching that root will be processed, and matching output files generated.")  
arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. Three groups: host ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the host ID will be the BAM file name.")
arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
arg_parser$add_argument("-s", "--splitsRule", action="store", default="r", help="The rules by which the sets of hosts are split into groups in order to ensure that all groups can be members of connected subgraphs without causing conflicts. Currently available: s=Sankoff (slow, rigorous), r=Romero-Severson (quick, less rigorous with >2 hosts).")
arg_parser$add_argument("-b", "--blacklist", action="store", help="A .csv file listing tips to ignore. Alternatively, a base name that identifies blacklist files. Blacklist files are assumed to end in .csv.")
arg_parser$add_argument("-k", "--kParam", action="store", default=0, help="The k parameter in the cost matrix for Sankoff reconstruction (see documentation)")
arg_parser$add_argument("-P", "--pruneBlacklist", action="store_true", help="If present, all blacklisted and references tips (except the outgroup) are pruned away before starting.")
arg_parser$add_argument("-p", "--proximityThreshold", action="store", default=0, help="The parameter determining when nodes return to the unsampled state; p for splitsRule = s, q for splitsRule = f; see documentation.")
arg_parser$add_argument("-i", "--idFile", action="store", help="Optionally, a full list of hosts; this will be used to standardise colours in PDF output across multiple runs for the same data")
arg_parser$add_argument("-R", "--readCountsMatterOnZeroBranches", help="If present, read counts will be taken into account in parsimony reconstructions at the parents of zero-length branches. Not applicable for the Romero-Severson-like reconstruction method.", default = FALSE, action="store_true")
arg_parser$add_argument("-m", "--multifurcationThreshold", help="If specified, short branches in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. If a number, this number will be used as the threshold. If 'g', it will be guessed from the branch lengths (use this only if you've checked by eye that the tree does indeed have multifurcations).")
arg_parser$add_argument("-n", "--branchLengthNormalisation", default = 1, action="store", help="If present and a number, a normalising constant for all branch lengths in the tree or trees. If present and a file, the path to a .csv file with two columns: tree file name and normalising constant")
arg_parser$add_argument("-D", "--scriptDir", action="store", help="Full path of the /tools directory.")
arg_parser$add_argument("-OD", "--outputdir", action="store", help="Full path of the directory for output; if absent, current working directory")
arg_parser$add_argument("-pw", "--pdfwidth", action="store", default=100, help="Width of tree pdf in inches.")
arg_parser$add_argument("-ph", "--pdfrelheight", action="store", default=0.15, help="Relative height of tree pdf.")
arg_parser$add_argument("-ff", "--useff", action="store_true", help="Use ff to store parsimony reconstruction matrices. Use if you run out of memory.")
arg_parser$add_argument("-v", "--verbose", action="store_true", help="Talk about what I'm doing.")
arg_parser$add_argument("-tfe", "--treeFileExtension", action="store", default="tree", help="The file extension for tree files (default .tree).")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")
arg_parser$add_argument("-ORDA", "--outputAsRDA", action="store_true", help="Store output in rda format.")
arg_parser$add_argument("inputFileName", action="store", help="The file or file root for the input tree in Newick format")
arg_parser$add_argument("outputFileID", action="store", help="A string shared by all output file names.")

args <- arg_parser$parse_args()

if(!is.null(args$scriptDir)){
  script.dir          <- args$scriptDir
} else {
  script.dir          <- dirname(thisfile())
  if(!dir.exists(script.dir)){
    stop("Cannot detect the location of the /phyloscanner/tools directory. Please specify it at the command line with -D.")
  }
}

tree.file.name         <- args$inputFile
sankoff.k              <- as.numeric(args$kParam)
sankoff.p              <- as.numeric(args$proximityThreshold)
output.file.ID         <- args$outputFileID
blacklist.file.name    <- args$blacklist
root.name              <- args$outgroupName
read.counts.matter     <- args$readCountsMatterOnZeroBranches
verbose                <- args$verbose
outputAsRDA	           <- args$outputAsRDA
prune.blacklist        <- args$pruneBlacklist
tip.regex              <- args$tipRegex
mode                   <- args$splitsRule
pdf.hm                 <- as.numeric(args$pdfrelheight)
pdf.w                  <- as.numeric(args$pdfwidth)
useff                  <- args$useff
tree.fe                <- args$treeFileExtension
csv.fe                 <- args$csvFileExtension
has.normalisation      <- !is.null(args$branchLengthNormalisation)
normalisation.argument <- args$branchLengthNormalisation

if(!is.null(args$outputdir)){
  output.dir         <- args$outputdir
} else {
  output.dir         <- getwd()
}

use.m.thresh <- !is.null(args$multifurcationThreshold)
if(use.m.thresh){
  if(args$multifurcationThreshold=="g"){
    m.thresh          <- NA
  } else if(!is.na(as.numeric(args$multifurcationThreshold))){
    m.thresh          <- as.numeric(args$multifurcationThreshold)
  } else {
    stop("Unknown argument for --multifurcationThreshold specified\n")
  }
}

if(!use.m.thresh & read.counts.matter){
  stop("Please specify a multifurcation collapse threshold with -m to count reads on zero length branches in the reconstruction")
}

if(is.null(root.name)){
  cat("No outgroup name given; will assume the tree is already rooted\n")
}

source(file.path(script.dir, "tree_utility_functions.R"))
source(file.path(script.dir, "parsimony_reconstruction_methods.R"))
source(file.path(script.dir, "collapsed_tree_methods.R"))
source(file.path(script.dir, "write_annotated_trees.R"))
source(file.path(script.dir, "general_functions.R"))

if(useff){
  suppressMessages(library(ff, quietly=TRUE, warn.conflicts=FALSE))
}

if(!(mode %in% c("r", "s", "f"))){
  stop(paste("Unknown split classifier: ", mode, "\n", sep=""))
}

if(!is.null(args$idFile)){
  host.master.list <- read.csv(args$idFile, header = F)[,1]
} else {
  host.master.list <- NULL
}

#
#	start script
#

if(file.exists(tree.file.name)){

  all.tree.info <- list()
  
  tree.info <- list()
  tree.info$tree.input <- tree.file.name

  tree.file.names <- tree.file.name
  if(!is.null(blacklist.file.name)){
    tree.info$blacklist.input <- blacklist.file.name
  }
  
  tree.info$output.ID <- output.file.ID
  
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

  tree.info$normalisation.constant <- normalisation.constant
  
  all.tree.info[[tree.file.name]] <- tree.info
  
} else {
  
  # Assume we are dealing with a group of files
  
  tree.file.names	<- list.files.mod(dirname(tree.file.name), pattern=paste('^',basename(tree.file.name),".*",tree.fe,'$',sep=''), full.names=TRUE)

  if(length(tree.file.names)==0){
    cat("No input trees found.\n")
    quit(save="no", status=1)
  }
  
  suffixes <- substr(tree.file.names, nchar(tree.file.name) + 1, nchar(tree.file.names)-nchar(tree.fe)-1)
  if(!is.null(blacklist.file.name)){  
    blacklist.file.names <- paste(blacklist.file.name, suffixes, ".", csv.fe, sep="")
  }
  output.file.IDs <- paste(output.file.ID, "_", suffixes, sep="")
  
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
  all.tree.info <- split(fn.df, rownames(fn.df))
}



for(i in all.tree.info){
  
  # we now do all the tree processing here
  
  tree.file.name		<- i$tree.input
  if (verbose) cat("Processing file ",tree.file.name,"...\n",sep="")
  
  if(!is.null(blacklist.file.name)){
    blacklist.file	<- i$blacklist.input
  } else {
    blacklist.file <- NULL
  }
  
  output.string <- i$output.ID 
  
  tree <- read.tree(tree.file.name)
  
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
    tips.to.go <- blacklist
    blacklist <- vector()
  } else {
    tips.to.go <- vector()
  }
  
  if(use.m.thresh){
    if(is.na(m.thresh)){
      min.bl <- min(tree$edge.length)
      m.thresh <- min.bl*1.001
    } 
  } else {
    m.thresh <- -1
  }

  tree <- process.tree(tree, root.name, m.thresh, tips.to.go, i$normalisation.constant)
  
  tmp					    <- split.hosts.to.subgraphs(tree, blacklist, mode, tip.regex, sankoff.k, sankoff.p, useff, read.counts.matter, host.master.list, verbose)
  tree					  <- tmp[['tree']]	
  rs.subgraphs		<- tmp[['rs.subgraphs']]
  
  # reconvert the tree to starting branch lengths
  
  tree$edge.length <- tree$edge.length * i$normalisation.constant
  
  #
  #	write rda file (potentially before plotting fails so we can recover)
  #
  if(outputAsRDA){
	  tmp 				<- file.path(output.dir,paste('subgraphs_',mode,'_',output.string,'.rda',sep='')) 
  	if (verbose) cat("Writing output to file",tmp,"...\n")
  	save(rs.subgraphs, tree, file=tmp)  
  }	
  #
  #	plot tree
  #
  if(!is.null(rs.subgraphs)){		
	  if (verbose) cat("Drawing tree...\n")
	  tree.display <- ggtree(tree, aes(color=BRANCH_COLOURS)) +
							    geom_point2(shape = 16, size=1, aes(color=INDIVIDUAL)) +
							    scale_fill_hue(na.value = "black", drop=F) +
							    scale_color_hue(na.value = "black", drop=F) +
							    theme(legend.position="none") +
							    geom_tiplab(aes(col=INDIVIDUAL)) + 
							    geom_treescale(width=0.01, y=-5, offset=1.5)	  
	  x.max <- ggplot_build(tree.display)$layout$panel_ranges[[1]]$x.range[2]	  
	  tree.display <- tree.display + ggplot2::xlim(0, 1.1*x.max)
	  tree.display		
	  tmp	<- file.path(output.dir,paste('Tree_',output.string,'.pdf',sep=''))
	  if (verbose) cat("Plot to file",tmp,"...\n")
	  ggsave(tmp, device="pdf", height = pdf.hm*length(tree$tip.label), width = pdf.w, limitsize = F)
  }
  #
  # write nexus file
  #
  if(!is.null(rs.subgraphs)){
	  tmp <- file.path(output.dir, paste('ProcessedTree_',mode,'_',output.string,".", tree.fe,sep=''))
	  if (verbose) cat("Writing rerooted, multifurcating, annotated tree to file",tmp,"...\n")	
	  write.ann.nexus(tree, file=tmp, annotations = c("INDIVIDUAL", "SPLIT"))
  }
  #
  #	write csv file
  #
  if(!is.null(rs.subgraphs)){
	  tmp 				<- file.path(output.dir, paste('subgraphs_',mode,'_',output.string,".",csv.fe,sep=''))
	  if (verbose) cat("Writing output to file",tmp,"...\n")	
	  write.csv(rs.subgraphs, file=tmp, row.names = F, quote=F)		  
  }
}
