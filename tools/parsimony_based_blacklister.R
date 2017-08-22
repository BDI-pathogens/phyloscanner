#!/usr/bin/env Rscript

list.of.packages <- c("argparse", "ape", "phangorn", "kimisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Please run package_install.R to continue\n")
  quit(save="no", status=1)
}

suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
# Define arguments

arg_parser = ArgumentParser(description="Identify phylogeny tips for blacklisting as representing suspected contaminants, based on a Sankoff parsimony reconstruction. Input and output file name arguments are either a single file, or the root name for a group of files, in which case all files matching that root will be processed, and matching output files generated.")
arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. Three capture groups, in order: host ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the host ID will be the BAM file name.")
arg_parser$add_argument("-D", "--scriptDir", action="store", help="Full path of the /tools directory.")
arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
arg_parser$add_argument("-b", "--blacklist", action="store", help="A blacklist to be applied before this script is run.")
arg_parser$add_argument("-c", "--noReadCounts", action="store_true", help="If present, read counts are not taken from tip labels and each tip is assumed to represent one read")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what the script is doing.")
arg_parser$add_argument("-n", "--branchLengthNormalisation", action="store", help="If present and a number, a normalising constant for all branch lengths in the tree or trees. If present and a file, the path to a .csv file with two columns: tree file name and normalising constant")
arg_parser$add_argument("-m", "--multifurcationThreshold", help="If specified, short branches in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. If a number, this number will be used as the threshold. If 'g', it will be guessed from the branch lengths (use this only if you've checked by eye that the tree does indeed have multifurcations).")
arg_parser$add_argument("-d", "--dualsOutputFile", action="store", help="A file or file root to write the set of hosts which seem to be dually infected according to the parameters of this run; if unspecified, do not output this.")
arg_parser$add_argument("-tfe", "--treeFileExtension", action="store", default="tree", help="The file extension for tree files (default .tree).")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")
arg_parser$add_argument("-dpr", "--developmentParsimonyReconstruction", default=F, action="store_true")
arg_parser$add_argument("rawThreshold", action="store", type="double", help="Raw threshold; subgraphs with read counts less than this will be blacklisted, regardless of the count of any other subgraphs from the same host")
arg_parser$add_argument("ratioThreshold", action="store", type="double", help="Ratio threshold; subgraphs will be blacklisted if the ratio of their tip count to that of another subgraph from the same host is less than this.")
arg_parser$add_argument("sankoffK", action="store", type="double", help="The k parameter in the cost matrix for Sankoff reconstruction (see documentation)")
arg_parser$add_argument("-p", "--sankoffP", action="store", type="double", default=0, help="The p parameter in the cost matrix for Sankoff reconstruction (see documentation)")
arg_parser$add_argument("inputFileName", action="store", help="The file or file root for the input tree in Newick or Nexus format")
arg_parser$add_argument("blacklistOutputFileName", action="store", help="The file or file root where the list of tips to be blacklisted will be written")


# Parse arguments
args                   <- arg_parser$parse_args()

if(!is.null(args$scriptDir)){
  script.dir          <- args$scriptDir
} else {
  suppressMessages(library(kimisc, quietly=TRUE, warn.conflicts=FALSE))
  script.dir          <- dirname(thisfile())
  if(!dir.exists(script.dir)){
    stop("Cannot detect the location of the /phyloscanner/tools directory. Please specify it at the command line with -D.")
  }
}

raw.threshold          <- args$rawThreshold
ratio.threshold        <- args$ratioThreshold
sankoff.k              <- args$sankoffK
sankoff.p              <- args$sankoffP
sankoff.method         <- if(args$developmentParsimonyReconstruction==T) "f" else "s"
tip.regex              <- args$tipRegex
input.name             <- args$inputFileName
b.output.name          <- args$blacklistOutputFileName
d.output.name          <- args$dualsOutputFile
root.name              <- args$outgroupName
blacklist.file.name    <- args$blacklist
no.read.counts         <- args$noReadCounts
verbose                <- args$verbose
normalisation.argument <- args$branchLengthNormalisation
tree.fe                <- args$treeFileExtension
csv.fe                 <- args$csvFileExtension

has.normalisation      <- !is.null(normalisation.argument)

use.m.thresh           <- !is.null(args$multifurcationThreshold)

if(use.m.thresh){
  if(args$multifurcationThreshold=="g"){
    m.thresh          <- NA
  } else if(!is.na(as.numeric(args$multifurcationThreshold))){
    m.thresh          <- as.numeric(args$multifurcationThreshold)
  } else {
    stop("Unknown argument for --multifurcationThreshold specified\n")
  }
} else {
  m.thresh            <- -1       
}

if(verbose) cat("Loading libraries...\n")

suppressMessages(library(ape, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(phangorn, quietly=TRUE, warn.conflicts=FALSE))

if(verbose) cat("Reading functions...\n")

source(file.path(script.dir, "tree_utility_functions.R"))
source(file.path(script.dir, "parsimony_reconstruction_methods.R"))
source(file.path(script.dir, "general_functions.R"))
source(file.path(script.dir, "blacklist_functions.R"))

if(is.null(root.name)){
  cat("No outgroup name given; will assume the tree is rooted and use a random other tip as an outgroup for each host\n")
}

all.tree.info <- list()

if(file.exists(input.name)){
  
  # single mode
  
  tree.info <- list()
  tree.info$tree.input <- input.name
  tree.info$blacklist.input <- blacklist.file.name
  tree.info$duals.output <- d.output.name
  tree.info$blacklist.output <- b.output.name

  if(!has.normalisation){
    tree.info$normalisation.constant <- 1
  } else {
    normalisation.constant <- suppressWarnings(as.numeric(normalisation.argument))

    if(is.na(normalisation.constant)){
      norm.table <- read.csv(normalisation.argument, stringsAsFactors = F, header = F)
      if(nrow(norm.table[which(norm.table[,1]==basename(input.name)),])==1){
        tree.info$normalisation.constant <- as.numeric(norm.table[which(norm.table[,1]==basename(input.name)),2])
      } else if (nrow(norm.table[which(norm.table[,1]==basename(input.name)),])==0){
        cat(paste("WARNING: No normalisation constant given for tree file name ", basename(input.name), " in file ",normalisation.argument,"; normalised distances will not be given", sep=""))
        tree.info$normalisation.constant <- 1
      } else {
        cat(paste("Two or more entries for normalisation constant for tree file name ",basename(input.name)," in file ",normalisation.argument,"; taking the first", sep=""))
        tree.info$normalisation.constant <- as.numeric(norm.table[which(norm.table[,1]==basename(input.name))[1],2])
      }
    } else {
      tree.info$normalisation.constant <- normalisation.constant
    }
  }
  
  all.tree.info[[input.name]] <- tree.info
  
} else {
  
  # batch mode
  
  tree.input.names <- list.files.mod(dirname(input.name), pattern=paste0("^",basename(input.name),'.*',tree.fe,'$'), full.names=TRUE)
  
  if(length(tree.input.names)==0){
    stop("No tree files found.")
  }
  
  suffixes <- substr(tree.input.names, nchar(input.name) + 1, nchar(tree.input.names) - nchar(tree.fe) - 1)
  
  b.output.names <- paste(b.output.name, "_", suffixes, ".", csv.fe, sep="")
  if(!is.null(d.output.name)){
    d.output.names <- paste0(d.output.name, "_", suffixes, ".", csv.fe)
  }
  if(!is.null(blacklist.file.name)){
    b.input.names <- paste0(blacklist.file.name, suffixes, ".", csv.fe)
  }
  
  if(has.normalisation){
    normalisation.constant <- suppressWarnings(as.numeric(normalisation.argument))
    if(is.na(normalisation.constant)){
      norm.table <- read.csv(normalisation.argument, stringsAsFactors = F, header = F)
      normalisation.constants <- sapply(basename(tree.input.names), function(x){
        if(x %in% norm.table[,1]){
          if(length(which(norm.table[,1]==x))>1){
            cat(paste("WARNING: Two or more entries for normalisation constant for tree file name ",x, " in file ",normalisation.argument,"; taking the first\n", sep=""))
          }
          return(norm.table[which(norm.table[,1]==x)[1],2])
        } else {
          cat(paste("WARNING: No normalisation constant for tree file ",x," found in file, ",normalisation.argument,"; this tree will not have branch lengths normalised.\n", sep=""))
          return(1)
        }
      })  
      
    } else {
      # it's a number to be used in every tree
      normalisation.constants <- rep(normalisation.constant, length(tree.input.names))
    }
  } else {
    normalisation.constants <- rep(1, length(tree.input.names))
  }
  
  for(suffix.no in 1:length(suffixes)){
    suffix <- suffixes[suffix.no]
    
    tree.info <- list()
    tree.info$tree.input <- tree.input.names[suffix.no]
    tree.info$blacklist.output <- b.output.names[suffix.no]
    tree.info$normalisation.constant <- normalisation.constants[suffix.no]
    if(!is.null(d.output.name)){
      tree.info$duals.output <- d.output.names[suffix.no]
    }
    if(!is.null(blacklist.file.name)){
      tree.info$blacklist.input <- b.input.names[suffix.no]
    }
    all.tree.info[[suffix]] <- tree.info
  }

}

for(tree.info in all.tree.info){
  
  if (verbose) cat(paste("Reading tree (",tree.info$tree.input,")...\n",sep=""))
  
  tree <- read.tree(tree.info$tree.input)
  
  if(tree.info$normalisation.constant!=1){
    if (verbose) cat("Normalising tree branch length (nc = ",tree.info$normalisation.constant,")\n", sep="")
  }
  
  tree$edge.length <- tree$edge.length/tree.info$normalisation.constant
  
  # Re-root the tree
  
  tree <- unroot(tree)
  
  if(!is.null(root.name)){
    outgroup.no <- which(tree$tip.label==root.name)
    tree <- root(tree, outgroup = outgroup.no, resolve.root = T)
  }

  if(is.na(m.thresh)){
    m.thresh                          <- min(tree$edge.length*1.0001) 
  }
  
  tip.labels <- tree$tip.label
  
  tree.info$tree <- tree
  
  # Import existing blacklist

  blacklist <- vector()
  
  if(!is.null(blacklist.file.name)){
    if(file.exists(tree.info$blacklist.input)){
      if (verbose) cat("Reading blacklist file",tree.info$blacklist.input,'\n')
      blacklisted.tips <- read.table(tree.info$blacklist.input, sep=",", header=F, stringsAsFactors = F, col.names="read")
      if(nrow(blacklisted.tips)>0){
        blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
      }
    } else {
      cat(paste("WARNING: File ",tree.info$blacklist.input," does not exist; skipping.\n",sep=""))
    }
  } 
  
  if (verbose) cat("Collecting tips for each host...\n")
  
  tip.hosts <- sapply(tip.labels, function(x) host.from.label(x, tip.regex))
  tip.hosts[blacklist] <- NA

  
  hosts <- unique(na.omit(tip.hosts))
  
  new.blacklist <- tree$tip.label[blacklist]
  
  # Re-order the hosts (mostly for the sake of getting a sense of progress in screen output)
  
  hosts <- hosts[order(hosts)]
  
  tree <- process.tree(tree, root.name, m.thresh)
  
  results <- lapply(hosts, function(x) get.splits.for.host(x, tip.hosts, tree, root.name, raw.threshold, ratio.threshold, sankoff.method, sankoff.k, sankoff.p, !is.null(tree.info$duals.output), no.read.counts, verbose))

  if (verbose) cat("Finished\n")
  
  new.blacklist <- c(new.blacklist, unlist(lapply(results, "[[", 2)))
  
  # Write the output
  
  write.table(new.blacklist, tree.info$blacklist.output, sep=",", row.names=FALSE, col.names=FALSE, quote=F)
  
  if(!is.null(tree.info$duals.output)){
    
    which.are.duals <- which(unlist(lapply(results, "[[", 6)))
    
    if(length(which.are.duals) > 0) {
      mi.df <- data.frame(host = unlist(sapply(results[which.are.duals], function (x) rep(x$id, length(x$tip.names)) )), 
                          tip.name = unlist(lapply(results[which.are.duals], "[[", 3)),
                          reads.in.subtree = unlist(lapply(results[which.are.duals], "[[", 4)),
                          tips.in.subtree = unlist(lapply(results[which.are.duals], "[[", 5))
      )
      
      write.csv(mi.df, tree.info$duals.output, row.names = F)
    } else {
      # Other scripts behave better when these files exist even if they are empty
      
      cat("No multiple infections detected; writing empty file.\n")
      file.create(tree.info$duals.output)
    }
  }
}
