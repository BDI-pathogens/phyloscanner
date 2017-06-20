list.of.packages <- c("argparse", "ape", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Please run PackageInstall.R to continue\n")
  quit(save="no", status=1)
}

suppressMessages(library(data.table, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ape, quietly=TRUE, warn.conflicts=FALSE))

suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))

arg_parser = ArgumentParser(description="Look at identified multiple infections across all windows, and identify which are to be ignored entirely or reduced to largest subtree only. Requires the output of ParsimonyBasedBlacklister.R.")

arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. Three groups, in order: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")
arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.")
arg_parser$add_argument("-b", "--existingBlacklistsPrefix", action="store", help="A file path and initial string identifying existing (.csv) blacklists whose output is to be appended to. Only such files whose suffix (the string after the prefix, without the file extension) matches a suffix from the dual reports will be appended to.")
arg_parser$add_argument("-s", "--summaryFile", action="store", help="A file to write a summary of the results to, detailing window counts for all patients")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")
arg_parser$add_argument("threshold", action="store", help="The proportion of windows that a dual infection needs to appear in to be considered genuine. Those patients whose dual infections are considered genuine are blacklisted entirely. Those who are not are reduced to their largest subtrees only.")
arg_parser$add_argument("treePrefix", action="store", help="A file path and initial string identifying read trees of all analyzed windows.")
arg_parser$add_argument("dualReportsPrefix", action="store", help="A file path and initial string identifying all dual infection files output from ParsimonyBasedBlacklister.R.")
arg_parser$add_argument("newBlacklistsPrefix", action="store", help="A file path and initial string for all output blacklist files.")

# Read in the arguments

args                   <- arg_parser$parse_args()

script.dir             <- args$scriptdir
tree.prefix            <- args$treePrefix
existing.bl.prefix     <- args$existingBlacklistsPrefix
duals.prefix           <- args$dualReportsPrefix
output.prefix          <- args$newBlacklistsPrefix 
threshold              <- args$threshold
summary.file           <- args$summaryFile
verbose                <- args$verbose
tip.regex              <- args$tipRegex
csv.fe                 <- args$csvFileExtension

source(file.path(script.dir, "TreeUtilityFunctions.R"))
source(file.path(script.dir, "GeneralFunctions.R"))
source(file.path(script.dir, "BlacklistFunctions.R"))

dual.files  <- list.files.mod(dirname(duals.prefix), pattern=paste('^',basename(duals.prefix),sep=""), full.names=TRUE)
tree.files	<- list.files.mod(dirname(tree.prefix), pattern=paste0('^',basename(tree.prefix)), full.names=TRUE)

suffixes	  <- substr(dual.files, nchar(duals.prefix)+1, nchar(dual.files)-nchar(csv.fe)-1)

suffixes    <- suffixes[order(suffixes)]
tree.files  <- tree.files[order(suffixes)]
dual.files  <- dual.files[order(suffixes)]

output.bls  <- sapply(suffixes, function(x) paste0(output.prefix, x, ".", csv.fe, sep=""))

if(!is.null(existing.bl.prefix)){
  expected.blacklists <- paste(existing.bl.prefix, suffixes, ".", csv.fe, sep="")
  observed.bl.files <- list.files.mod(dirname(existing.bl.prefix), pattern=paste('^',basename(existing.bl.prefix),sep=""), full.names=TRUE)
  expected.but.not.seen <- setdiff(expected.blacklists, observed.bl.files)
  seen.but.not.expected <- setdiff(observed.bl.files, expected.blacklists)
  if(length(expected.but.not.seen)>0){
    for(ebns in expected.but.not.seen){
      cat("Blacklist file ",ebns," not found. Will make a new blacklist file with this extension.\n",sep="")
    }
  }
}

all.tree.info <- list()

hosts <- vector()

if(verbose) cat("Collecting window data...\n")

for(suffix.index in 1:length(suffixes)){
  suffix               <- suffixes[suffix.index]
  out                  <- list()
  out$suffix           <- suffix
  out$tree             <- read.tree(tree.files[suffix.index])
  out$hosts.for.tips   <- sapply(out$tree$tip.label, function(x) host.from.label(x, tip.regex))
  duals.table          <- fread(dual.files[suffix.index], stringsAsFactors = F)
  out$duals.info       <- duals.table
  out$bl.output.name   <- output.bls[suffix.index]
  hosts                <- unique(c(hosts, duals.table$host))
  if(!is.null(existing.bl.prefix)){
    if(file.exists(expected.blacklists[suffix.index])){
      out$existing.blacklist <- read.csv(expected.blacklists[suffix.index], stringsAsFactors = F, header = F)
    } else {
      out$existing.blacklist <- vector()
    }
  }
  all.tree.info[[suffix]]  <- out
}

hosts <- hosts[order(hosts)]

for(suffix.index in 1:length(suffixes)){
  all.tree.info[[suffix]]$dual.infection  <- sapply(hosts, function(x) x%in% duals.table$host, simplify = F, USE.NAMES = T)
}


#	Count number of potential dual windows by patient

if(verbose) cat("Making new blacklists...\n")

results <- blacklist.duals(all.tree.info, hosts, threshold, summary.file, verbose)

for(tree.info in all.tree.info){
  if(verbose) cat("Writing new blacklist file ", tree.info$bl.output.name, "...\n")
  tree.info$blacklist <- results[tree.info$name]
  write.table(tree.info$blacklist, tree.info$bl.output.name, sep=",", row.names=FALSE, col.names=FALSE, quote=F)
}

