#!/usr/bin/env Rscript

list.of.packages <- c("argparse", "phangorn", "kimisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Missing dependencies; replacing the path below appropriately, run\n[path to your phyloscanner code]/tools/package_install.R\nthen try again.\n")
  quit(save="no", status=1)
}
options("warn"=1)

suppressMessages(library(argparse))
suppressMessages(library(phangorn))
suppressMessages(library(kimisc))

arg_parser = ArgumentParser(description="Examine all hosts from the input tree that match a given regexp and blacklist tips with low read counts that are distant from the main clades of tips (and hence likely contaminants).")

arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. Three groups: host ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the host ID will be the BAM file name.")
arg_parser$add_argument("-b", "--blacklist", action="store", help="A blacklist to be applied before this script is run.")
arg_parser$add_argument("-H", "--hosts", action="store", help="A file listing which hosts to downsample (on separate lines, with no header). If absent, then downsample every host present.")
arg_parser$add_argument("-s", "--seed", action="store", help="Random number seed (integer value).")
arg_parser$add_argument("-nrc", "--noReadCounts", action="store_true", help="If present, treat each tip as a single read")
arg_parser$add_argument("-R", "--rename", action="store", help="If present, tip labels will be rewritten to include the number of reads that were sampled, not the total in the starting data. A new tree file will be output to this argument or, if there are multiple tree files, to this base name. This is recommended if this script is used and the -R option is used in split_hosts_to_subtrees.R, which should then be given the relabelled tree as input.")
arg_parser$add_argument("-e", "--excludeUnderrepresented", action="store_true", help="If present, hosts without enough reads will be blacklisted entirely from this tree. If absent, all reads from that host will be included.")
arg_parser$add_argument("maxReadsPerHost", type="double", action="store", help="The upper limit for the number of reads to be included from each host")
arg_parser$add_argument("inputFile", metavar="inputTreeFileName", help="Tree file name. Alternatively, a base name that identifies a group of tree file names can be specified. Tree files are assumed to end in .tree.")  
arg_parser$add_argument("outputFile", metavar="outputFileName", help="The file to write the output to, a list of tips to be blacklisted.")  
arg_parser$add_argument("-D", "--scriptDir", action="store", help="Full path of the script directory.")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what the script is doing.")
arg_parser$add_argument("-tfe", "--treeFileExtension", action="store", default="tree", help="The file extension for tree files (default .tree).")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")

args <- arg_parser$parse_args()

if(!is.null(args$scriptDir)){
  script.dir          <- args$scriptDir
} else {
  script.dir          <- dirname(thisfile())
  if(!dir.exists(script.dir)){
    stop("Cannot detect the location of the /phyloscanner/tools directory. Please specify it at the command line with -D.")
  }
}

tip.regex <- args$tipRegex
input.file.name <- args$inputFile
output.file.name <- args$outputFile
blacklist.file.name <- args$blacklist
max.reads <- args$maxReadsPerHost
seed <- ifelse(is.null(args$seed), NA_real_,args$seed)  
hosts.file.name <- args$hosts
no.read.counts <- args$noReadCounts
rename <- !is.null(args$rename)
if(rename){
  renamed.file.name <- args$rename
}
exclude.underrepresented <- args$excludeUnderrepresented
tree.fe <- args$treeFileExtension
csv.fe <- args$csvFileExtension

verbose <- args$verbose
source(file.path(script.dir, "tree_utility_functions.R"))
source(file.path(script.dir, "parsimony_reconstruction_methods.R"))
source(file.path(script.dir, "downsampling_functions.R"))
source(file.path(script.dir, "general_functions.R"))

all.tree.info <- list()


if(file.exists(input.file.name)){
  
  #	option 1: input.file.name specifies a single input file
  
  tree.info <- list()
  tree.info$tree.input <- input.file.name
  
  if(verbose) cat("Reading tree ", input.file.name, "\n", sep="")
  tree <- read.tree(input.file.name)
  
  tree.info$tree <- tree
  
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
      warning("File ",blacklist.file.name," does not exist; skipping.",paste="")
    }	
  }
  
  tree.info$blacklist        <- blacklist
  
  tree.info$blacklist.output <- output.file.name
  if(rename){
    tree.info$renamed.tree.file.name <- renamed.file.name
  }

  
  all.tree.info[[input.file.name]] <- tree.info
} else {
  #	option 2: input.file.name specifies a root name for several input files

  input.file.names <- sort(list.files.mod(dirname(input.file.name), pattern=paste(basename(input.file.name),'.*\\.', tree.fe,'$',sep=''), full.names=TRUE))

  if(length(input.file.names)==0){
    cat("No tree files found.\nQuitting.\n")
    quit(save="no", status=1)
  }  
  suffixes <- substr(basename(input.file.names), nchar(basename(input.file.name)) + 1, nchar(basename(input.file.names))-nchar(tree.fe)-1)
  
  b.output.names <- paste(output.file.name, suffixes, ".", csv.fe, sep="")  
  if(rename){
    rt.output.names <- paste(renamed.file.name, suffixes, ".", tree.fe, sep="")
  }
  if(!is.null(blacklist.file.name)){
    b.input.names <- paste(blacklist.file.name, suffixes, ".", csv.fe, sep="")
  }
  
  fn.df <- data.frame(row.names = suffixes, suffix = suffixes, tree.input = input.file.names, blacklist.output = b.output.names, stringsAsFactors = F)
  if(rename){
    fn.df$renamed.tree.file.name <- rt.output.names
  }
  if(!is.null(blacklist.file.name)){
    fn.df$blacklist.input <- b.input.names
  }
  
  all.tree.info <- split(fn.df, rownames(fn.df))
  
  all.tree.info <- sapply(all.tree.info, function(tree.info){
    tree.info <- as.list(tree.info)
    
    if(verbose) cat("Reading tree ", tree.info$tree.input, "\n", sep="")
    tree <- read.tree(tree.info$tree.input)
    
    blacklist <- vector()
    
    if(!is.null(tree.info$blacklist.input)){
      
      if(file.exists(tree.info$blacklist.input)){
        if(verbose) cat("Reading existing blacklist file",tree.info$blacklist.input,'\n')
        blacklisted.tips <- read.table(tree.info$blacklist.input, sep=",", header=F, stringsAsFactors = F, col.names="read")
        if(nrow(blacklisted.tips)>0){
          blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
        }
      } else {
        warning(paste("File ",tree.info$blacklist.input," does not exist; skipping.",paste=""))
      }	
    }
  
    tree.info$tree             <- tree
    tree.info$blacklist        <- blacklist
    
    tree.info
  }, simplify = F, USE.NAMES = T)
  

}


if(!is.null(hosts.file.name)){
  hosts <- read.table(hosts.file.name, header = F, stringsAsFactors = F)[,1]
} else {
  hosts <- NULL
}

for(tree.info in all.tree.info){
  tree.info <- downsample.tree(tree.info, hosts, max.reads, rename, exclude.underrepresented, no.read.counts, seed, verbose)
  
  if(rename){
    if(verbose){
      cat("Writing new tree to ",tree.info$renamed.tree.file.name,"...\n",sep = "")
    }
    write.tree(tree.info$tree, tree.info$renamed.tree.file.name)
  }
  
  if(verbose){
    cat("Writing new blacklist to ",tree.info$blacklist.output,"...\n",sep = "")
  }
  
  write.table(tree.info$tree$tip.label[tree.info$blacklist],  tree.info$blacklist.output, sep=",", row.names=FALSE, col.names=FALSE, quote=F)
  
}

