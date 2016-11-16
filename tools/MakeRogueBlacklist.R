list.of.packages <- c("argparse", "phangorn")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

suppressMessages(library(argparse))
suppressMessages(library(phangorn))

arg_parser = ArgumentParser(description="Examine all patients from the input tree that match a given regexp and blacklist tips with low read counts that are distant from the main clades of tips (and hence likely contaminants).")

arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                        help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
arg_parser$add_argument("-b", "--blacklist", action="store", help="A blacklist to be applied before this script is run.")
arg_parser$add_argument("dropProportion", type="double", action="store", help="Distant tips are to be blacklisted if the ratio of reads on that tip to the total number of reads from this patient is less than this threshold.")
arg_parser$add_argument("longestBranchLength", type="double", action="store", help="Tips are classified as distant if, in the subtree obtained by dropping every tip not from the patient in question, the longest branch is longer than this.")
arg_parser$add_argument("inputFile", metavar="inputTreeFileName", help="Tree file name. Alternatively, a base name that identifies a group of tree file names can be specified. Tree files are assumed to end in .tree.")  
arg_parser$add_argument("outputFile", metavar="outputFileName", help="The file to write the output to, a list of tips to be blacklisted.")  
arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")

args <- arg_parser$parse_args()
script.dir <- args$scriptdir
tip.regex <- args$tipRegex
input.file.name <- args$inputFile
output.file.name <- args$outputFile
blacklist.file.name <- args$blacklist
root.name <- args$outgroupName
drop.prop <- args$dropProportion
if(drop.prop>=0.5){
  stop("Drop proportion threshold must be less than 0.5")
}
branch.limit <- args$longestBranchLength

if(F){
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/")
  drop.prop <- 0.01
  branch.limit <- 0.1
  tip.regex <- "^(.*)_read_([0-9]+)_count_([0-9]+)$"
  outgroup <- "REF_CPX_AF460972"
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools/"
  input.file.name <- "ptyr42_trees_newick/ptyr42_InWindow_1225_to_1474.tree"
  blacklist.file.name <- "ptyr42_otherstuff/blacklist_test.txt"
  output.file.name <- "ptyr42_otherstuff/blacklist_test_2.txt"
}

source(file.path(script.dir, "TransmissionUtilityFunctions.R"))
source(file.path(script.dir, "SubtreeMethods.R"))

tree <- read.tree(input.file.name)
tree <- unroot(tree)
tree <- di2multi(tree, tol = 1E-5)
if(!is.null(root.name)){
  tree <- root(tree, outgroup = which(tree$tip.label==root.name), resolve.root = T)
}

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

tree.1 <- drop.tip(tree, blacklist)

patients.present <- sapply(tree.1$tip.label, function(tip) patient.from.label(tip, tip.regex))
patients.present <- patients.present[!is.na(patients.present)]
patients.present <- unique(patients.present)

patient.ids <- sapply(tree.1$tip.label, function(x) patient.from.label(x, tip.regex))

rogue.hunt <- function(tree, patient, patient.ids, length.threshold, read.prop.threshold){
  tips.to.keep <- which(patient.ids==patient)
  new.blacklist <- vector()
  if(length(tips.to.keep) > 1){
    tree.2 <- drop.tip(tree, tip=tree$tip.label[setdiff(seq(1, length(tree$tip.label)), tips.to.keep)])
    if(length(tree.2$tip.label)>2){
      tree.2 <- unroot(tree.2)
      total.reads <- sum(as.numeric(sapply(tree.2$tip.label, function(tip) read.count.from.label(tip, tip.regex))))
      longest.branch <- max(tree.2$edge.length)
      if(longest.branch >= length.threshold){
        long.edge.ends <- tree.2$edge[which(tree.2$edge.length>=length.threshold),]
        which.end <- vector()
        groups.identified <- rep(FALSE, length(tree.2$tip.label))
        counter <- 1
        for(tip.no in seq(1, length(tree.2$tip.label))){
          if(!groups.identified[tip.no]){
            tips.in.group <- tips.reachable(tree.2, tip.no, long.edge.ends)
            groups.identified[tips.in.group] <- TRUE
            which.end[tips.in.group] <- counter
            counter <- counter + 1
          }
        }
        tips.going <- vector()
        for(group in unique(which.end)){
          group.tips <- tree.2$tip.label[which(which.end==group)]
          if(length(group.tips)>0){
            reads.in.group <- sum(as.numeric(sapply(group.tips, function(tip) read.count.from.label(tip, tip.regex))))
            if(reads.in.group/total.reads < read.prop.threshold){
              tips.going <- c(tips.going, group.tips)
            }
          }
        }
        new.blacklist <- c(new.blacklist, tips.going)
      }
    } else {
      total.length = sum(tree.2$edge.length)
      if(total.length > length.threshold){
        reads.1 <- as.numeric(read.count.from.label(tree.2$tip.label[1], tip.regex))
        reads.2 <- as.numeric(read.count.from.label(tree.2$tip.label[2], tip.regex))
        if(reads.1/reads.2 < read.prop.threshold){
          new.blacklist <- c(new.blacklist,1)
        }
        if(reads.2/reads.1 < read.prop.threshold){
          new.blacklist <- c(new.blacklist,2)
        }
      }
    }
  }
  return(new.blacklist)
}

new.blacklists <- lapply(patients.present, function(pat) rogue.hunt(tree.1, pat, patient.ids, branch.limit, drop.prop))
new.blacklist <- unlist(new.blacklists)
new.blacklist <- c(tree$tip.label[blacklist], new.blacklist)
write.table(new.blacklist, output.file.name, sep=",", row.names=FALSE, col.names=FALSE, quote=F)
