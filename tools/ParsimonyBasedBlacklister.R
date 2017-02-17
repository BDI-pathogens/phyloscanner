list.of.packages <- c("argparse", "phangorn")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))

command.line <- T

if(command.line){
  arg_parser = ArgumentParser(description="Identify phylogeny tips for blacklisting as representing suspected contaminants, based on a Sankhoff parsimony reconstruction")
  
  arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                          help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
  arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
  arg_parser$add_argument("-b", "--blacklist", action="store", help="A blacklist to be applied before this script is run.")
  arg_parser$add_argument("-c", "--noReadCounts", action="store_true", help="If present, each tip is assumed to represent one read")
  arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")
  arg_parser$add_argument("rawThreshold", action="store", type="double", help="Raw threshold; tips with read counts less than this that are identical to a tip from another patient will be blacklisted, regardless of the count of the other read.")
  arg_parser$add_argument("ratioThreshold", action="store", type="double", help="Ratio threshold; tips will be blacklisted if the ratio of their tip count to that of of another, identical tip from another patient is less than this value.")
  arg_parser$add_argument("sankhoffK", action="store", type="double", help="The k parameter in the cost matrix for Sankhoff reconstruction (see documentation)")
  arg_parser$add_argument("inputFileName", action="store", help="A CSV file outlining groups of tips that have identical sequences, each forming a single line.")
  arg_parser$add_argument("blacklistOutputFileName", action="store", help="The file to write a list of tips to be blacklisted to.")
  arg_parser$add_argument("dualCandidatesOutputFileName", action="store", help="The file to write the set of patients who seem to be dually infected according to the parameters of this run.")
  
  
  # Parse arguments
  
  args <- arg_parser$parse_args()
  
  raw.threshold <- args$rawThreshold
  ratio.threshold <- args$ratioThreshold
  sankhoff.k <- args$sankhoffK
  tip.regex <- args$tipRegex
  input.name <- args$inputFileName
  b.output.name <- args$blacklistOutputFileName
  d.output.name <- args$dualCandidatesOutputFileName
  script.dir <- args$scriptdir
  root.name <- args$outgroupName
  blacklist.file.name <- args$blacklist
  no.read.counts <- args$noReadCounts
  verbose <- args$verbose
  
} else {
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20161013/AllTrees/")
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools/"
  
  input.name <- "RAxML_bestTree.InWindow_1550_to_1900.tree"
  
  #anything already blacklisted
  blacklist.file.name <- "AmpliconBlacklist_InWindow_1550_to_1900.csv"
  
  root.name <- "C.BW.00.00BW07621.AF443088"
  tip.regex <- "^(.*)-[0-9].*_read_([0-9]+)_count_([0-9]+)$"
  
  sankhoff.k <- 20
  
  # Threshold for propotion of reads that come from a patient to be in a split for it to be considered a dual infection, not a
  # rogue
  
  raw.threshold <- 3
  ratio.threshold <- 0.005
}

get.count <- function(string){
  if(length(grep(tip.regex, string)>0)) {
    return(as.numeric(sub(tip.regex, "\\3", string)))
  } else {
    return(NA)
  }
}

cat("Reading functions...\n")
source(file.path(script.dir, "TransmissionUtilityFunctions.R"))
source(file.path(script.dir, "SubtreeMethods.R"))

cat(paste("Reading tree (",input.name,")...\n",sep=""))

tree <- read.tree(input.name)
tree <- unroot(tree)
outgroup.no <- which(tree$tip.label==root.name)

if(!is.null(root.name)){
  tree <- root(tree, outgroup = outgroup.no, resolve.root = T)
}

tip.labels <- tree$tip.label

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

patient.ids <- sapply(tip.labels, function(x) patient.from.label(x, tip.regex))
patient.ids[blacklist] <- NA

patients <- unique(na.omit(patient.ids))

new.blacklist <- tree$tip.label[blacklist]

mi.patient.column <- vector()
mi.tip.names.column <- vector()
mi.read.count.column <- vector()
mi.tip.count.column <- vector()

for(patient in patients[order(patients)]){
  
  if(length(which(patient.ids==patient))>1){
    tip.nos <- c(outgroup.no, which(patient.ids==patient))
    
    subtree <- drop.tip(tree, setdiff(1:length(tree$tip.label), tip.nos))
    subtree <- di2multi(subtree, tol = 1E-5)
    
    root.no <- which(subtree$tip.label==root.name)
    
    patient.tips <- list()
    patient.tips[[patient]] <- setdiff(1:length(subtree$tip.label), root.no)
    
    split.results <- split.and.annotate(subtree, patient, patient.tips, NULL, NULL, tip.regex, "s", sankhoff.k, "u", sankhoff.mode = "total.lengths", useff = F)
    
    if(!no.read.counts){
      total.reads <- sum(sapply(subtree$tip.label[setdiff(1:length(subtree$tip.label), root.no)], function(x) read.count.from.label(x, tip.regex)))
    } else {
      total.reads <- length(subtree$tip.label)
    }
    
    if(length(split.results$split.patients)==1){
      if(verbose){
        cat("No splits for patient ", patient, "\n", sep="")
      }
      if(total.reads < raw.threshold){
        if(verbose){
          cat("Blacklisting ",patient,"; not enough reads in total.\n", sep="")
        }
        new.blacklist <- c(new.blacklist, subtree$tip.label[which(subtree$tip.label!=root.name)])
      }
      
    } else {
      
      cat(length(split.results$split.patients), " splits for patient ", patient, "\n", sep="")
      
      
      
      props <- vector()
      too.small <- vector()
      
      tips.vector <- vector()
      read.count.vector <- vector()
      tip.count.vector <- vector()
      
      for(split in split.results$split.patients){
        if(!no.read.counts){
          count.in.split <- sum(sapply(subtree$tip.label[split.results$split.tips[[split]]], function(x) read.count.from.label(x, tip.regex)))
        } else {
          count.in.split <- length(subtree$tip.label[split.results$split.tips[[split]]])
        }
        prop.in.split <- count.in.split/total.reads
        props <- c(props, prop.in.split)
        too.small <- c(too.small, count.in.split < raw.threshold | prop.in.split < ratio.threshold)
        
        tips.vector <- c(tips.vector, subtree$tip.label[split.results$split.tips[[split]]])
        read.count.vector <- c(read.count.vector, rep(count.in.split, length(split.results$split.tips[[split]])))
        tip.count.vector <- c(tip.count.vector, rep(length(split.results$split.tips[[split]]), length(split.results$split.tips[[split]])))
        
      }
      if(length(which(!too.small))>1){
        if(verbose){
          cat(patient, " looks like a dual infection\n", sep="")
        }
        mi.patient.column <- c(mi.patient.column, rep(patient, length(subtree$tip.label)-1))
        mi.tip.names.column <- c(mi.tip.names.column, tips.vector)
        mi.read.count.column <- c(mi.read.count.column, read.count.vector)
        mi.tip.count.column <- c(mi.tip.count.column, tip.count.vector)
        
      }
      for(small.group in which(too.small)){
        cat("Tips from ", split.results$split.patients[small.group], " blacklisted\n", sep="")
        new.blacklist <- c(new.blacklist, subtree$tip.label[split.results$split.tips[[split.results$split.patients[small.group]]]])
      }
    }
  } 
}


write.table(new.blacklist, b.output.name, sep=",", row.names=FALSE, col.names=FALSE, quote=F)

if(length(mi.patient.column>0)){
  mi.df <- data.frame(patient = mi.patient.column, tip.name = mi.tip.names.column, reads.in.subtree = mi.read.count.column,
                      tips.in.subtree = mi.tip.count.column)
  
  write.csv(mi.df, d.output.name, row.names = F)
} else {
  cat("No multiple infections detected.\n")
}
