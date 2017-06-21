list.of.packages <- c("argparse", "phangorn")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

suppressMessages(library(argparse))
suppressMessages(library(phangorn))

tree.fe <- ".tree"
csv.fe <- ".csv"

arg_parser = ArgumentParser(description="Examine all patients from the input tree that match a given regexp and blacklist tips with low read counts that are distant from the main clades of tips (and hence likely contaminants).")

arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                        help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
arg_parser$add_argument("-b", "--blacklist", action="store", help="A blacklist to be applied before this script is run.")
arg_parser$add_argument("-p", "--patients", action="store", help="A file listing which patients to downsample (on separate lines, with no header). If absent, then downsample every patient present.")
arg_parser$add_argument("-s", "--seed", action="store", help="Random number seed (integer value). Can be ommitted.")
arg_parser$add_argument("-m", "--multifurcationThreshold", help="If specified, short branches in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. If a number, this number will be used as the threshold. If 'g', it will be guessed from the branch lengths (use this only if you've checked by eye that the tree does indeed have multifurcations).")
arg_parser$add_argument("-R", "--rename", action="store", help="If present, tip labels will be rewritten to include the number of reads that were sampled, not the total in the starting data. A new tree file will be output to this argument or, if there are multiple tree files, to this base name. This is recommended if this script is used and the -R option is used in SplitPatientsToSubtrees.R, which should then be given the relabelled tree as input.")
arg_parser$add_argument("-e", "--excludeUnderrepresented", action="store", help="If present, hosts without enough reads will be blacklisted entirely from this tree. If absent, all reads from that host will be included.")
arg_parser$add_argument("maxReadsPerPatient", type="double", action="store", help="The upper limit for the number of reads to be included from each patient")
arg_parser$add_argument("inputFile", metavar="inputTreeFileName", help="Tree file name. Alternatively, a base name that identifies a group of tree file names can be specified. Tree files are assumed to end in .tree.")  
arg_parser$add_argument("outputFile", metavar="outputFileName", help="The file to write the output to, a list of tips to be blacklisted.")  
arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what the script is doing.")

args <- arg_parser$parse_args()
script.dir <- args$scriptdir
tip.regex <- args$tipRegex
input.file.name <- args$inputFile
output.file.name <- args$outputFile
blacklist.file.name <- args$blacklist
root.name <- args$outgroupName
max.reads <- args$maxReadsPerPatient
seed <- ifelse(is.null(args$seed), NA_real_,args$seed)  
patients.file.name <- args$patients
rename <- !is.null(args$rename)
if(rename){
  renamed.file.name <- args$rename
}
exclude.underrepresented <- args$excludeUnderrepresented

verbose <- args$verbose
use.m.thresh <- !is.null(args$multifurcationThreshold)
if(use.m.thresh){
  if(args$multifurcationThreshold=="g"){
    m.thresh <- NA
  } else if(!is.na(as.numeric(args$multifurcationThreshold))){
    m.thresh <- as.numeric(args$multifurcationThreshold)
  } else {
    cat("Unknown argument for -m specified\n")
    quit(save="no", status=1)
  }
} else {
  m.thresh <- -1
}

source(file.path(script.dir, "TreeUtilityFunctions.R"))
source(file.path(script.dir, "ParsimonyReconstructionMethods.R"))
source(file.path(script.dir, "DownsamplingFunctions.R"))

all.tree.info <- list()

if(file.exists(input.file.name)){
  
  #	option 1: input.file.name specifies a single input file
  
  tree.info <- list()
  tree.info$tree.input <- input.file.name
  tree.info$blacklist.input <- blacklist.file.name
  
  blacklist <- vector()
  
  if(!is.null(blacklist.file.name)){
    if(file.exists(blacklist.file.name)){
      if(verbose) cat("Reading existing blacklist file",blacklist.file.name,'\n')
      blacklisted.tips <- read.table(blacklist.file.name, sep=",", header=F, stringsAsFactors = F, col.names="read")
      if(nrow(blacklisted.tips)>0){
        blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
      }
    } else {
      warning(paste("File ",blacklist.file.name," does not exist; skipping.",paste=""))
    }	
  }
  
  tree.info$blacklist        <- blacklist
  
  tree.info$blacklist.output <- output.file.name
  if(rename){
    tree.info$renamed.tree.output <- renamed.file.name
  }
  
  
  tree <- read.tree(input.file.name)
  
  tree <- process.tree(tree, root.name, m.thresh, tree.info$blacklist)
  
  tree.info$tree <- tree
  
  all.tree.info[[input.file.name]] <- tree.info
} else {
  #	option 2: input.file.name specifies a root name for several input files
  input.file.names <- sort(list.files(dirname(input.file.name), pattern=paste(basename(input.file.name),'.*\\.tree$',sep=''), full.names=TRUE))
  
  if(length(input.file.names)==0){
    cat("No tree files found.\nQuitting.\n")
    quit(save="no", status=1)
  }
  
  suffixes <- substr(input.file.names, nchar(input.file.name) + 1, nchar(input.file.names)-nchar(tree.fe) )
  
  b.output.names <- paste(output.file.name, suffixes, csv.fe, sep="")
  if(rename){
    rt.output.names <- paste(renamed.file.name, suffixes, tree.fe, sep="")
  }
  if(!is.null(blacklist.file.name)){
    b.input.names <- paste(blacklist.file.name, suffixes, csv.fe, sep="")
  }
  
  fn.df <- data.frame(row.names = suffixes, tree.input = input.file.names, blacklist.output = b.output.names, stringsAsFactors = F)
  if(rename){
    fn.df$renamed.tree.output <- rt.output.names
  }
  if(!is.null(blacklist.file.name)){
    fn.df$blacklist.input <- b.input.names
  }
  
  all.tree.info <- split(fn.df, rownames(fn.df))
  
  all.tree.info <- sapply(all.tree.info, function(tree.info){
    
    blacklist <- vector()
    
    if(!is.null(tree.info$blacklist.file.name)){
      if(file.exists(tree.info$blacklist.file.name)){
        if(verbose) cat("Reading existing blacklist file",tree.info$blacklist.file.name,'\n')
        blacklisted.tips <- read.table(tree.info$blacklist.file.name, sep=",", header=F, stringsAsFactors = F, col.names="read")
        if(nrow(blacklisted.tips)>0){
          blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
        }
      } else {
        warning(paste("File ",tree.info$blacklist.file.name," does not exist; skipping.",paste=""))
      }	
    }
    
    tree.info$blacklist        <- blacklist
    
    tree <- read.tree(tree.info$tree.file.name)
    
    tree <- process.tree(tree, root.name, m.thresh, blacklist)
    
    tree.info$tree        <- tree
    
    tree.info
  }, simplify = F, USE.NAMES = T)
}

hosts <- read.table(patients.file.name, header = F, stringsAsFactors = F)

for(tree.info in all.tree.info){
  tree.info <- downsample.tree(tree.info, hosts, max.reads, rename, exclude.underrepresented, seed, verbose)
  
  if(rename){
    if(verbose){
      cat("Writing new tree to ",renamed.tree.file.name,"...\n",sep = "")
    }
    write.tree(tree.info$tree, renamed.tree.file.name)
  }
  
  if(verbose){
    cat("Writing new blacklist to ",output.file.name,"...\n",sep = "")
  }
  write.table(tree.info$tip.label[tree.info$blacklist],  tree.info$blacklist.output, sep=",", row.names=FALSE, col.names=FALSE, quote=F)
  
}

