list.of.packages <- c("argparse", "ape", "phangorn")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Please run PackageInstall.R to continue\n")
  quit(save='no')
}

suppressMessages(library(ape, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(phangorn, quietly=TRUE, warn.conflicts=FALSE))

cat("Reading functions...\n")

tree.fe <- ".tree"
csv.fe <- ".csv"

command.line <- T
if(command.line){
  suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
  # Define arguments
  
  arg_parser = ArgumentParser(description="Identify phylogeny tips for blacklisting as representing suspected contaminants, based on a Sankhoff parsimony reconstruction. Input and output file name arguments are either a single file, or the root name for a group of files, in which case all files matching that root will be processed, and matching output files generated.")
  arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. Three groups, in order: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.")
  arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
  arg_parser$add_argument("-b", "--blacklist", action="store", help="A blacklist to be applied before this script is run.")
  arg_parser$add_argument("-c", "--noReadCounts", action="store_true", help="If present, read counts are not taken from tip labels and each tip is assumed to represent one read")
  arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what the script is doing.")
  arg_parser$add_argument("-n", "--branchLengthNormalisation", action="store", help="If present and a number, a normalising constant for all branch lengths in the tree or trees. If present and a file, the path to a .csv file with two columns: tree file name and normalising constant")
  arg_parser$add_argument("-m", "--multifurcationThreshold", help="If specified, short branches in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. If a number, this number will be used as the threshold. If 'g', it will be guessed from the branch lengths (use this only if you've checked by eye that the tree does indeed have multifurcations).")
  arg_parser$add_argument("-d", "--dualsOutputFile", action="store", help="A file or file root to write the set of patients which seem to be dually infected according to the parameters of this run; if unspecified, do not output this.")
  arg_parser$add_argument("rawThreshold", action="store", type="double", help="Raw threshold; subgraphs with read counts less than this will be blacklisted, regardless of the count of any other subgraphs from the same patient")
  arg_parser$add_argument("ratioThreshold", action="store", type="double", help="Ratio threshold; subgraphs will be blacklisted if the ratio of their tip count to that of another subgraph from the same patient is less than this.")
  arg_parser$add_argument("sankhoffK", action="store", type="double", help="The k parameter in the cost matrix for Sankhoff reconstruction (see documentation)")
  arg_parser$add_argument("inputFileName", action="store", help="The file or file root for the input tree in Newick format")
  arg_parser$add_argument("blacklistOutputFileName", action="store", help="The file or file root to write a list of tips to be blacklisted to.")
  
  # Parse arguments
  args <- arg_parser$parse_args()
  raw.threshold <- args$rawThreshold
  ratio.threshold <- args$ratioThreshold
  sankhoff.k <- args$sankhoffK
  tip.regex <- args$tipRegex
  input.name <- args$inputFileName
  b.output.name <- args$blacklistOutputFileName
  d.output.name <- args$dualsOutputFile
  script.dir <- args$scriptdir
  root.name <- args$outgroupName
  blacklist.file.name <- args$blacklist
  no.read.counts <- args$noReadCounts
  verbose <- args$verbose
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
  has.normalisation <- !is.null(args$branchLengthNormalisation)
  normalisation.argument <- args$branchLengthNormalisation
  
  source(file.path(script.dir, "TreeUtilityFunctions.R"))
  source(file.path(script.dir, "ParsimonyReconstructionMethods.R"))
  source(file.path(script.dir, "GeneralFunctions.R"))
  
  if(is.null(root.name)){
    cat("No outgroup name given; will assume the tree is rooted and use a random other tip as an outgroup for each patient\n")
  }
  
  file.details <- list()
  
  if(file.exists(input.name)){
    
    file.name.list <- list()
    file.name.list$tree.input <- input.name
    file.name.list$blacklist.input <- blacklist.file.name
    file.name.list$duals.output <- d.output.name
    file.name.list$blacklist.output <- b.output.name
    file.name.list$blacklist.input <- blacklist.file.name
    
    if(!has.normalisation){
      file.name.list$normalisation.constant <- 1
    } else {
      normalisation.constant <- suppressWarnings(as.numeric(normalisation.argument))
      
      if(is.na(normalisation.constant)){
        norm.table <- read.csv(normalisation.argument, stringsAsFactors = F, header = F)
        if(nrow(norm.table[which(norm.table[,1]==basename(input.name)),])==1){
          file.name.list$normalisation.constant <- as.numeric(norm.table[which(norm.table[,1]==basename(input.name)),2])
        } else if (nrow(norm.table[which(norm.table[,1]==basename(input.name)),])==0){
          cat(paste("WARNING: No normalisation constant given for tree file name ", basename(input.name), " in file ",normalisation.argument,"; normalised distances will not be given", sep=""))
          file.name.list$normalisation.constant <- 1
        } else {
          cat(paste("Two or more entries for normalisation constant for tree file name ",basename(input.name)," in file ",normalisation.argument,"; taking the first", sep=""))
          file.name.list$normalisation.constant <- as.numeric(norm.table[which(norm.table[,1]==basename(input.name))[1],2])
        }
      } else {
        file.name.list$normalisation.constant <- normalisation.constant
      }
    }

    file.details[[input.name]] <- file.name.list

  } else {
    # Assume we are dealing with a group of files
    
    tree.input.names <- list.files.mod(dirname(input.name), pattern=paste(basename(input.name),'.*\\',tree.fe,'$',sep=''), full.names=TRUE)

    if(length(tree.input.names)==0){
      cat("No tree files found.\nQuitting.\n",)
      quit(save="no")
    }
    
    suffixes <- substr(tree.input.names, nchar(input.name) + 1, nchar(tree.input.names))
    suffixes <- gsub(paste('\\', tree.fe, sep="") , csv.fe ,suffixes)
    
    b.output.names <- paste(b.output.name, suffixes, sep="")
    if(!is.null(d.output.name)){
      d.output.names <- paste(d.output.name, suffixes, sep="")
    }
    if(!is.null(blacklist.file.name)){
      b.input.names <- paste(blacklist.file.name, suffixes, sep="")
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
    
    fn.df <- data.frame(row.names = suffixes, tree.input = tree.input.names, blacklist.output = b.output.names, normalisation.constant = normalisation.constants, stringsAsFactors = F)
    if(!is.null(d.output.name)){
      fn.df$duals.output <- d.output.names
    }
    if(!is.null(blacklist.file.name)){
      fn.df$blacklist.input <- b.input.names
    }
    
    file.details <- split(fn.df, rownames(fn.df))
    
  }
  
} else {
  setwd("/Users/mdhall/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/Final6MiseqForPaper/trees/")
  script.dir <- "/Users/mdhall/phylotypes/tools/"
  
  input.name <- "RAxML_bestTree.InWindow_4200_to_4519.tree"
  
  #anything already blacklisted
  blacklist.file.name <- NULL
  
  root.name <- "C.BW.00.00BW07621.AF443088"
  tip.regex <- "^(.*)-[0-9].*_read_([0-9]+)_count_([0-9]+)$"
  
  sankhoff.k <- 20
  
  # Threshold for propotion of reads that come from a patient to be in a split for it to be considered a dual infection, not a
  # rogue
  
  raw.threshold <- 3
  ratio.threshold <- 0.005
  
  b.output.names <- "test1"
  d.output.names <- "test2"
  
  no.read.counts <- F
  verbose <- T
}


# The main function for getting splits for a given patient ID

get.splits.for.patient <- function(patient, tip.patients, tree, root.name, raw.threshold, ratio.threshold, d.output.name, verbose){
  cat("Identifying splits for patient ", patient, "\n", sep="")
  
  blacklist.items <- vector()
  
  dual <- F
  
  if(length(which(tip.patients==patient))>1){
    
    # Keep only the tips from this patient and the outgroup
    
    if(!is.null(root.name)){
      outgroup.no <- which(tree$tip.label==root.name)
    } else {
      patient.mrca <- mrca.phylo(tree, which(tip.patients==patient))
      if(patient.mrca == getRoot(tree)){
        message("No suitable outgroup found for patient ",patient,"; try adding one and specifying its tip name with --outgroupName", sep="")
        quit(save = F)
      }
      
      outgroup.no <- sample(1:length(tree$tip.label), 1)
      while(patient.mrca %in% Ancestors(tree, a.random.tip)){
        outgroup.no <- sample(1:length(tree$tip.label), 1)
      }
    }
    
    tip.nos <- c(outgroup.no, which(tip.patients==patient))
    
    subtree <- drop.tip(tree, setdiff(1:length(tree$tip.label), tip.nos))
    
    # Collapse multifurcations
    
    if(use.m.thresh){
      if(is.na(m.thresh)){
        min.bl <- min(tree$edge.length)
        subtree <- di2multi(subtree, tol = min.bl*1.001)
      } else {
        subtree <- di2multi(subtree, tol = m.thresh)
      }
    }
    
    # Re-find the root
    
    st.outgroup.no <- which(subtree$tip.label==root.name)
    
    # di2multi actually unroots the tree!
    
    subtree <- root(subtree, st.outgroup.no, resolve.root = T)
    
    # Set up split.and.annotate (which takes a list)
    
    patient.tips <- list()
    patient.tips[[patient]] <- setdiff(1:length(subtree$tip.label), st.outgroup.no)
    
    # Perform split.and.annotate; get a list of splits
    
    split.results <- split.and.annotate(subtree, patient, patient.tips, NULL, NULL, tip.regex, "s", sankhoff.k, 0, root.name, "u", useff = F)
    # vector of of split IDs
    
    patient.split.ids <- split.results$split.patients
    
    # list of tips that belong to each split ID
    
    tips.for.splits <- split.results$split.tips
    
    # The total number of reads is either got from the tips or assumed to be one per tip
    
    if(!no.read.counts){
      reads.per.tip <- sapply(subtree$tip.label, function(x) read.count.from.label(x, tip.regex))
    } else {
      reads.per.tip <- rep(1, length(subtree$tip.label))
      reads.per.tip[root.no] <- NA
    }
    
    total.reads <- sum(reads.per.tip[!is.na(reads.per.tip)])
    
    if(length(patient.split.ids)==1){
      
      # if the algorithm didn't split the tips
      
      if(verbose){
        cat("No splits for patient ", patient, "\n", sep="")
      }
      if(total.reads < raw.threshold){
        # blacklist the whole patient if the read count from its only subgraph is below the raw threshold
        
        if(verbose){
          cat("Blacklisting ",patient,"; not enough reads in total.\n", sep="")
        }
        blacklist.items <- subtree$tip.label[which(subtree$tip.label!=root.name)]
      }
      
    } else {
      
      # at least two subgraphs
      
      cat(length(patient.split.ids), " splits for patient ", patient, "\n", sep="")
      
      props <- vector()
      too.small <- vector()
      
      tips.vector <- vector()
      read.count.vector <- vector()
      tip.count.vector <- vector()
      
      too.small <- sapply(patient.split.ids, function(x) check.read.count.for.split(x, tips.for.splits, raw.threshold, ratio.threshold, reads.per.tip, total.reads))
      
      if(length(which(!too.small))>1 & !is.null(d.output.name)){
        # there are at least two subgraphs which are too big to be contaminants
        
        if(verbose){
          cat(patient, " looks like a dual infection\n", sep="")
        }
        
        dual <- T
      }
      
      # blacklist the tips from small subgraphs
      
      if(length(which(too.small))>0){
        blacklist.items <- subtree$tip.label[unlist(tips.for.splits[patient.split.ids[which(too.small)]])]
      }
      
    }
    
    tips.vector <- subtree$tip.label[unlist(tips.for.splits[patient.split.ids])]
    read.count.vector <- unlist(lapply(patient.split.ids, function(x){
      rep(sum(read.count.from.label(subtree$tip.label[tips.for.splits[[x]]], tip.regex)), length(tips.for.splits[[x]]))
    }))
    tip.count.vector <- unlist(lapply(patient.split.ids, function(x){
      rep(length(tips.for.splits[[x]]), length(tips.for.splits[[x]]))
    }))
  } else {
    # this covers the situation where there is only one patient tip (because ape has trouble with two-tip trees)
    
    tip.no <- which(tip.patients==patient)
    
    reads <- read.count.from.label(tree$tip.label[tip.no], tip.regex)
    
    if(reads < raw.threshold){
      if(verbose){
        cat("Blacklisting ",patient,"; not enough reads in total.\n", sep="")
      }
      blacklist.items <- tree$tip.label[tip.no]
    }
    
    tips.vector <- tree$tip.label[tip.no]
    read.count.vector <- reads
    tip.count.vector <- 1
  }
  
  
  list(id = patient, blacklist.items = blacklist.items, tip.names = tips.vector, read.counts = read.count.vector, tip.counts = tip.count.vector, dual=dual)
}

# Function that returns the third group from a regex as a numeric (the read count from a tip label)

get.count <- function(string){
  if(length(grep(tip.regex, string)>0)) {
    return(as.numeric(sub(tip.regex, "\\3", string)))
  } else {
    return(NA)
  }
}

# Return whether the number of reads in this split is below one of the thresholds

check.read.count.for.split <- function(split, tips.for.splits, raw.threshold, ratio.threshold, reads.per.tip, total.reads){
  # get the read counts
  
  count.in.split <- sum(reads.per.tip[tips.for.splits[[split]]])
  
  # find what proportion of all reads from this patient are in this subgraph and check against the thresholds
  
  prop.in.split <- count.in.split/total.reads
  return (count.in.split < raw.threshold | prop.in.split < ratio.threshold)
}


for(i in file.details){
  
  cat(paste("Reading tree (",i$tree.input,")...\n",sep=""))

  tree <- read.tree(i$tree.input)
  
  if(i$normalisation.constant!=1){
    cat("Normalising tree branch length (nc = ",i$normalisation.constant,")\n", sep="")
  }
  
  tree$edge.length <- tree$edge.length/i$normalisation.constant
  
  # Re-root the tree
  
  tree <- unroot(tree)
  
  if(!is.null(root.name)){
    outgroup.no <- which(tree$tip.label==root.name)
    tree <- root(tree, outgroup = outgroup.no, resolve.root = T)
  }

  tip.labels <- tree$tip.label
  
  # Import existing blacklist
  
  blacklist <- vector()
  
  if(!is.null(blacklist.file.name)){
    if(file.exists(i$blacklist.input)){
      cat("Reading blacklist file",i$blacklist.input,'\n')
      blacklisted.tips <- read.table(i$blacklist.input, sep=",", header=F, stringsAsFactors = F, col.names="read")
      if(nrow(blacklisted.tips)>0){
        blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
      }
    } else {
      cat(paste("WARNING: File ",i$blacklist.input," does not exist; skipping.\n",sep=""))
    }
  } 
  
  cat("Collecting tips for each patient...\n")
  
  tip.patients <- sapply(tip.labels, function(x) patient.from.label(x, tip.regex))
  tip.patients[blacklist] <- NA
  
  patients <- unique(na.omit(tip.patients))
  
  new.blacklist <- tree$tip.label[blacklist]
  
  # Re-order the patients (mostly for the sake of getting a sense of progress in screen output)
  
  patients <- patients[order(patients)]
  
  results <- lapply(patients, function(x) get.splits.for.patient(x, tip.patients, tree, root.name, raw.threshold, ratio.threshold, i$duals.output, verbose))
  
  cat("Finished\n")
  
  new.blacklist <- c(new.blacklist, unlist(lapply(results, "[[", 2)))
  
  # Write the output
  
  write.table(new.blacklist, i$blacklist.output, sep=",", row.names=FALSE, col.names=FALSE, quote=F)
  
  if(!is.null(i$duals.output)){
    
    which.are.duals <- which(unlist(lapply(results, "[[", 6)))
    
    if(length(which.are.duals) > 0) {
      mi.df <- data.frame(patient = unlist(sapply(results[which.are.duals], function (x) rep(x$id, length(x$tip.names)) )), 
                          tip.name = unlist(lapply(results[which.are.duals], "[[", 3)),
                          reads.in.subtree = unlist(lapply(results[which.are.duals], "[[", 4)),
                          tips.in.subtree = unlist(lapply(results[which.are.duals], "[[", 5))
      )
      
      write.csv(mi.df, i$duals.output, row.names = F)
    } else {
      # Other scripts behave better when these files exist even if they are empty
      
      cat("No multiple infections detected; writing empty file.\n")
      file.create(i$duals.output)
    }
  }
}