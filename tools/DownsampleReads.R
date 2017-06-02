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
} 

source(file.path(script.dir, "TreeUtilityFunctions.R"))
source(file.path(script.dir, "ParsimonyReconstructionMethods.R"))

#
#	additional function definitions
#

# This returns the ones to get rid of (i.e. blacklist) per patient for a given tree
downsample.patient <- function(patient, tree, number, tip.regex, patient.ids, rename=F, verbose=F){
  if(verbose) cat("Downsampling tips for patient ", patient, "...\n", sep="")
	tips.to.keep <- which(patient.ids==patient)
	labels.to.keep <- tree$tip.label[tips.to.keep]
	
	read.counts <- as.numeric(sapply(labels.to.keep, function(tip) read.count.from.label(tip, tip.regex))) 
	
	total.reads <- sum(read.counts)
	
	# if you have less reads than you are trying to sample, return them all
	
	if(total.reads <= number){
	  if(verbose) cat("Insufficient reads for downsampling at a count of ",number,", returning all.\n", sep="")
	  return(list (map=NULL, blacklist=vector()))
	}
	
	tip.props <- read.counts/total.reads
	
	sample <- rmultinom(1, size = number, prob = tip.props)
	
	if(rename){
	  if(verbose) cat("Renaming tree tips with new read counts...\n", sep="")
	  
	  # need a new regex to do the find and replace. This is probably preferable to asking the user for one.
	  
	  # remove the first opening bracket
	  new.tip.regex <- sub("\\(", "", tip.regex)
	  # remove the second opening bracket
	  new.tip.regex <- sub("\\(", "", new.tip.regex)
	  # remove the first closing bracket
	  new.tip.regex <- sub("\\)", "", new.tip.regex)
	  # remove the second closing bracket
	  new.tip.regex <- sub("\\)", "", new.tip.regex)
	  # Start a new capture group start after the third one
	  new.tip.regex <- sub("\\)", "\\)\\(", new.tip.regex)
	  # End a new capture group start before the third one
	  new.tip.regex <- sub("\\(", "\\)\\(", new.tip.regex)
	  # Deal with the beginning
	  if(substr(new.tip.regex,1,1)=="^"){
	    new.tip.regex <- gsub("\\^", "\\^\\(", new.tip.regex)
	  } else {
	    new.tip.regex <- paste0("^(", new.tip.regex)
	  }
	  # Deal with the end
	  if(substr(new.tip.regex,nchar(new.tip.regex),nchar(new.tip.regex))=="$"){
	    new.tip.regex <- gsub("\\$", "\\)\\$", new.tip.regex)
	  } else {
	    new.tip.regex <- paste0(new.tip.regex,')$')
	  }
	  
	  label.map <- lapply(1:length(labels.to.keep), function(x){
	    gsub(new.tip.regex, paste0("\\1", sample[x], "\\3"), labels.to.keep[x])
	  })
	  
	  names(label.map) <- labels.to.keep
	  
	  labels.to.keep <- unlist(label.map[labels.to.keep])
	} else {
	  label.map <- NULL
	}
	
	sampled.names <- labels.to.keep[which(sample > 0)]
  
	return(list(blacklist = setdiff(labels.to.keep, sampled.names),map=label.map))
}

# This returns the ones to get rid of (i.e. blacklist) for a given tree
downsample.tree<- function(input.file.name, blacklist.file.name, output.file.name, renamed.tree.file.name=NULL, patients.file.name, root.name, tip.regex, max.reads, use.m.thresh, seed=NA, multifurcation.threshold = NA, verbose=F) {
  
  if(verbose){
    cat("Downsampling reads on tree ",i$tree.input,"...\n",sep = "")
  }
  
  rename = !is.null(renamed.tree.file.name)
  
	if(!is.na(seed))
		set.seed(seed)
  
  if(verbose){
    cat("Loading tree...\n",sep = "")
  }
  
	tree <- read.tree(input.file.name)
	tree <- unroot(tree)
	
	
	if(use.m.thresh){
	  if(verbose){
	    cat("Collapsing short branches...\n",sep = "")
	  }
	  if(is.na(multifurcation.threshold)){
	    min.bl <- min(tree$edge.length)
	    multifurcation.threshold <- min.bl*1.001
	  } 
	  tree <- di2multi(tree, tol = multifurcation.threshold)
	}
	
	if(!is.null(root.name)){
		tree <- root(tree, outgroup = which(tree$tip.label==root.name), resolve.root = T)
	}
	
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
	
	tree.1 <- drop.tip(tree, blacklist)
	
	patients.present <- sapply(tree.1$tip.label, function(tip) patient.from.label(tip, tip.regex))
	patients.present <- patients.present[!is.na(patients.present)]
	patients.present <- unique(patients.present)
	
	if(!is.null(patients.file.name)){
	  patients.to.include <- read.table(patients.file.name, header = F, stringsAsFactors = F)
	  if(interaction(setdiff(patients.to.include, patients.present)) == 0){
	    stop(paste("No entries in ", patients.file.name, " are actually present in the tree", sep=""))
	  }
	  if(length(setdiff(patients.to.include, patients.present)) > 0){
	    warning(paste("Not all patients in file ", patients.file.name, " are present in the tree", sep=""))
	  }
	  patients.to.include <- intersection(patients.to.include, patients.present)
	  
	} else {
	  patients.to.include <- patients.present
	}
	
	patient.ids <- sapply(tree.1$tip.label, function(x) patient.from.label(x, tip.regex))
	
	new.tip.labels <- tree$tip.label
	excluded <- unlist(lapply(patients.to.include, function(x){
	                          result <- downsample.patient(x, tree=tree.1, number = max.reads, tip.regex=tip.regex, patient.ids=patient.ids, rename, verbose)
	                          if(rename & !is.null(result$map)){
	                            new.tip.labels <<- sapply(new.tip.labels, function(y){
	                              if(y %in% names(result$map)){
	                                result$map[[y]]
	                              } else {
	                                y
	                              }
	                            })
	                          }
	                          return(result$blacklist)
	                   }))
	
	if(rename){
	  tree$tip.label <- new.tip.labels
	  if(verbose){
	    cat("Writing new tree to ",renamed.tree.file.name,"...\n",sep = "")
	  }
	  write.tree(tree, renamed.tree.file.name)
	}
	
	new.blacklist <- c(tree$tip.label[blacklist], excluded)
	if(verbose){
	  cat("Writing new blacklist to ",output.file.name,"...\n",sep = "")
	}
	
	write.table(new.blacklist, output.file.name, sep=",", row.names=FALSE, col.names=FALSE, quote=F)
}
#
#	run script 
#

file.details <- list()

if(file.exists(input.file.name)){
  
  #	option 1: input.file.name specifies a single input file
  
  file.name.list <- list()
  file.name.list$tree.input <- input.file.name
  file.name.list$blacklist.input <- blacklist.file.name
  file.name.list$blacklist.output <- output.file.name
  if(rename){
    file.name.list$renamed.tree.output <- renamed.file.name
  }
  file.details[[input.file.name]] <- file.name.list

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
  
  file.details <- split(fn.df, rownames(fn.df))
}

for(i in file.details){
  downsample.tree(i$tree.input, i$blacklist.input, i$blacklist.output, i$renamed.tree.output, patients.file.name, root.name, tip.regex, max.reads, use.m.thresh, seed, m.thresh, verbose)
}

