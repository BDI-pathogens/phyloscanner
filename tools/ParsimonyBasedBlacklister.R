list.of.packages <- c("argparse", "phangorn")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))

command.line <- T
if(command.line){
  # Define arguments
  
  arg_parser = ArgumentParser(description="Identify phylogeny tips for blacklisting as representing suspected contaminants, based on a Sankhoff parsimony reconstruction. Input and output file name arguments are either a single file, or the root name for a group of files, in which case all files matching that root will be processed, and matching output files generated.")
  arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. Three groups, in order: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.")
  arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
  arg_parser$add_argument("-b", "--blacklist", action="store", help="A blacklist to be applied before this script is run.")
  arg_parser$add_argument("-c", "--noReadCounts", action="store_true", help="If present, read counts are not taken from tip labels and each tip is assumed to represent one read")
  arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what the script is doing.")
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
  
  if(is.null(root.name)){
    cat("No outgroup name given; will assume the tree is rooted and use a random other tip as an outgroup for each patient\n")
  }
  
  if(file.exists(input.name)){
    input.names <- input.name
    d.output.names <- d.output.name
    b.output.names <- b.output.name
  } else {
    # Assume we are dealing with a group of files
    
    input.names		<- sort(list.files(dirname(input.name), pattern=paste(basename(input.name),'.*\\.tree$',sep=''), full.names=TRUE))
    
    
    if(length(input.names)==0){
      cat("No input trees found,\n")
      quit()
    }
    
    suffixes <- substr(input.names, nchar(input.name) + 1, nchar(input.names))
    suffixes <- gsub('\\.tree','.csv',suffixes)
    
    b.output.names <- paste(b.output.name, suffixes, sep="")
    if(!is.null(d.output.name)){
      d.output.names <- paste(d.output.name, suffixes, sep="")
    }
  }
  
} else {
  
}

# The main function for getting splits for a given patient ID

get.splits.for.patient <- function(patient, tip.patients, tree, root.name, raw.threshold, ratio.threshold, d.output.name, verbose){
  cat("Identifying splits for patient ", patient, "\n")
  
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
    
    subtree <- di2multi(subtree, tol = 1E-5)
    
    # Re-find the root
    
    st.outgroup.no <- which(subtree$tip.label==root.name)
    
    # Set up split.and.annotate (which takes a list)
    
    patient.tips <- list()
    patient.tips[[patient]] <- setdiff(1:length(subtree$tip.label), st.outgroup.no)
    
    # Perform split.and.annotate; get a list of splits
    
    split.results <- split.and.annotate(subtree, patient, patient.tips, NULL, NULL, tip.regex, "s", sankhoff.k, "u", sankhoff.mode = "total.lengths", useff = F)
    
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

cat("Reading functions...\n")

source(file.path(script.dir, "TreeUtilityFunctions.R"))
source(file.path(script.dir, "ParsimonyReconstructionMethods.R"))

for(i in 1:length(input.names)){
  cat(paste("Reading tree (",input.names[i],")...\n",sep=""))
  
  tree <- read.tree(input.names[i])
  
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
    if(file.exists(blacklist.file.name)){
      cat("Reading blacklist file",blacklist.file.name,'\n')
      blacklisted.tips <- read.table(blacklist.file.name, sep=",", header=F, stringsAsFactors = F, col.names="read")
      if(nrow(blacklisted.tips)>0){
        blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
      }
    } else {
      warning(paste("File ",blacklist.file.name," does not exist; skipping.",paste=""))
    }
  } 
  
  cat("Collecting tips for each patient...\n")
  
  tip.patients <- sapply(tip.labels, function(x) patient.from.label(x, tip.regex))
  tip.patients[blacklist] <- NA
  
  patients <- unique(na.omit(tip.patients))
  
  new.blacklist <- tree$tip.label[blacklist]
  
  # Re-order the patients (mostly for the sake of getting a sense of progress in screen output)
  
  patients <- patients[order(patients)]
  
  results <- lapply(patients, function(x) get.splits.for.patient(x, tip.patients, tree, root.name, raw.threshold, ratio.threshold, d.output.names[i], verbose))
  
  cat("Finished\n")
  
  new.blacklist <- c(new.blacklist, unlist(lapply(results, "[[", 2)))
  
  # Write the output
  
  write.table(new.blacklist, b.output.names[i], sep=",", row.names=FALSE, col.names=FALSE, quote=F)
  
  if(!is.null(d.output.names[i])){
    
    which.are.duals <- which(unlist(lapply(results, "[[", 6)))
    
    if(length(which.are.duals) > 0) {
      mi.df <- data.frame(patient = unlist(sapply(results[which.are.duals], function (x) rep(x$id, length(x$tip.names)) )), 
                          tip.name = unlist(lapply(results[which.are.duals], "[[", 3)),
                          reads.in.subtree = unlist(lapply(results[which.are.duals], "[[", 4)),
                          tips.in.subtree = unlist(lapply(results[which.are.duals], "[[", 5))
      )
      
      write.csv(mi.df, d.output.names[i], row.names = F)
    } else {
      # Other scripts behave better when these files exist even if they are empty
      
      cat("No multiple infections detected; writing empty file.\n")
      file.create(d.output.name)
    }
  }
}