#!/usr/bin/env Rscript

options("warn"=1)

suppressMessages(require(ape, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phyloscannerR, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phytools, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ggtree, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(reshape, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gtable, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(grid, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gridExtra, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(RColorBrewer, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(scales, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(dtplyr, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(tidyverse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(kimisc, quietly=TRUE, warn.conflicts=FALSE))

arg_parser = ArgumentParser(description="Summarise patient statistics from each phylogeny across all windows.")

arg_parser$add_argument("-c", "--noReadCounts", action="store_true", default=FALSE, help="If present, tip labels have no read counts and each tip is assumed to represent a single read") 
arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                        help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
arg_parser$add_argument("-b", "--blacklists", help="An optional file path and initial string identifying blacklist files (file extension must be .csv).")
arg_parser$add_argument("-w", "--windowCoords", action="store", help="The path for a .csv file describing the position in the genome which each tree file represents (e.g. the midpoint of a window); used for the x-axis in output. The first column in the file should the tree file name, the second the coordinate. If not given, the script will attempt to obtain these from the file names (which should work if the input is from the phyloscanner pipeline).")
arg_parser$add_argument("idFile", action="store", help="A file containing a list of the IDs of all the patients to calculate and display statistics for.")
arg_parser$add_argument("treeFiles", action="store", help="A file path and initial string identifying all tree files (processed tree files output by SplitPatientsToSubgraphs.R file extension must be .tree).")
arg_parser$add_argument("splitsFiles", action="store",help="A file path and initial string identifying all splits files (file extension must be .csv).")
arg_parser$add_argument("outputBaseName", action="store", help="A path and string to begin the names of all output files.")
arg_parser$add_argument("-R", "--recombinationFiles", action="store", help="A file path and initial string identifying all recombination data files.")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about the script is doing.")
arg_parser$add_argument("-tfe", "--treeFileExtension", action="store", default="tree", help="The file extension for tree files (default .tree).")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")

# Read in the arguments

args                     <- arg_parser$parse_args()

id.file                  <- args$idFile
verbose                  <- args$verbose
tip.regex                <- args$tipRegex
tree.file.root           <- args$treeFiles
splits.file.root         <- args$splitsFiles
blacklist.file.root      <- args$blacklists  
output.root              <- args$outputBaseName
window.coords.file       <- args$windowCoords
recomb.file.root         <- args$recombinationFiles
tree.fe                  <- args$treeFileExtension
csv.fe                   <- args$csvFileExtension
recomb.files.exist       <- !is.null(recomb.file.root)
no.read.counts           <- args$noReadCounts

# Find the input files

tree.files <- list.files.mod(dirname(tree.file.root), pattern=paste('^',basename(tree.file.root),".*",tree.fe,'$',sep=""), full.names=TRUE)
splits.files <- list.files.mod(dirname(splits.file.root), pattern=paste('^',basename(splits.file.root),".*",csv.fe,'$',sep=""), full.names=TRUE)	  

blacklist.files <- NULL
if(!is.null(blacklist.file.root)){
  blacklist.files <- list.files.mod(dirname(blacklist.file.root), pattern=paste('^',basename(blacklist.file.root),".*",csv.fe,"$",sep=""), full.names=TRUE)
}

if(recomb.files.exist) { 
  recomb.files <- list.files.mod(dirname(recomb.file.root), pattern=paste('^',basename(recomb.file.root),".*",csv.fe,'$',sep=""), full.names=TRUE)
}

if(length(tree.files)==0){
  cat("No tree files found. \n Quitting. \n")
  quit(save="no", status=1)
}
if(length(splits.files)==0){
  cat("No subgraph files found. \n Quitting. \n")
  quit(save="no", status=1)
}



# Get the suffixes

tree.suffixes	  <- sapply(tree.files, function(x) get.suffix(x, tree.file.root, tree.fe))
splits.suffixes	<- sapply(splits.files, function(x) get.suffix(x, splits.file.root, csv.fe))

# Only take suffixes that have both a tree file and a subgraph file
ts.both.present <- intersect(tree.suffixes, splits.suffixes)
ts.both.present <- ts.both.present[order(ts.both.present)]
if(length(ts.both.present) < length(tree.files)){
  missing.suffixes <- setdiff(tree.suffixes, ts.both.present)
  no.partner.files <- paste(tree.file.root, missing.suffixes, tree.fe, sep="")
  
  for(nps in no.partner.files){
    warning(paste("No subgraph file found for tree file ", nps, "; will ignore this file", sep=""))
  }
} 
if(length(ts.both.present) < length(splits.files)){
  missing.suffixes <- setdiff(tree.suffixes, ts.both.present)
  no.partner.files <- paste(splits.file.root, missing.suffixes, csv.fe, sep="")
  
  for(nps in no.partner.files){
    warning(paste("No tree file found for subgraph file ", nps, "; will ignore this file", sep=""))
  }
} 

all.tree.info <- list()

for(suffix in ts.both.present){
  tree.info <- list()
  tree.info$tree.id <- suffix
  tree.info$tree.file.name <- paste0(tree.file.root, suffix, ".", tree.fe)
  tree.info$splits.file.name <- paste0(splits.file.root, suffix, ".", csv.fe)
  tree.info$normalisation.constant <- 1
  all.tree.info[[suffix]] <- tree.info
}

# Check which suffixes from the current list have a blacklist file and vice versa

if(!is.null(blacklist.file.root)){
  blacklist.suffixes <- sapply(blacklist.files, function(x) get.suffix(x, blacklist.file.root, csv.fe))
  
  # the only remaining tree files will have a splits file, no need to check the splits files separately
  bt.both.present <- intersect(ts.both.present, blacklist.suffixes)	
  
  for(suffix in bt.both.present){
    all.tree.info[[suffix]]$blacklist.file.name <- paste(blacklist.file.root, suffix, ".", csv.fe, sep="")
  }
  
  if(length(bt.both.present) < length(ts.both.present)){
    missing.suffixes <- setdiff(ts.both.present, bt.both.present)
    no.partner.files <- paste(tree.file.root, missing.suffixes, ".", tree.fe, sep="")
    
    for(nps in no.partner.files){
      warning(paste("No blacklist file found for tree file ",nps,"; assuming no tips were blacklisted", sep=""))
    }
  } 
  if(length(bt.both.present) < length(blacklist.suffixes)){
    missing.suffixes <- setdiff(blacklist.suffixes, bt.both.present)
    no.partner.files <- paste(blacklist.file.root, missing.suffixes, ".", csv.fe, sep="")
    
    for(nps in no.partner.files){
      warning(paste("Tree or subgraph file missing for blacklist file ",nps,"; ignoring this file", sep=""))
    }
  } 
}

# Check there's a one-to-one correspondence between tree files and
# recomb files.
if (recomb.files.exist) { 
  recomb.suffixes	<- sapply(recomb.files, function(x) get.suffix(x, recomb.file.root, csv.fe))
  recomb.but.no.tree.suffixes <- setdiff(recomb.suffixes, ts.both.present)
  tree.but.no.recomb.suffixes <- setdiff(ts.both.present, recomb.suffixes)
  
  if (length(recomb.but.no.tree.suffixes) > 0) {
    missing.files <- paste(recomb.file.root, recomb.but.no.tree.suffixes, ".", csv.fe, sep="")
    stop("Error: no tree file found for the following recombination files:", paste(missing.files, collapse=' '), "\nQuitting.\n")
  }
  if (length(tree.but.no.recomb.suffixes) > 0) {
    missing.files <- paste(tree.file.root, tree.but.no.recomb.suffixes, ".", csv.fe, sep="")
    stop("Error: no recombination file found for the following tree files:", paste(missing.files, collapse=' '), "\nQuitting.\n")
  }
  for(suffix in recomb.suffixes){
    all.tree.info[[suffix]]$recombination.file.name <-  paste(recomb.file.root, suffix, ".", csv.fe, sep="")
  }
}

# From now on we work with suffixes only

suffixes <- ts.both.present
suffixes <- suffixes[order(suffixes)]

if(!is.null(window.coords.file)){
  if (verbose) cat("Reading genome coordinates from file ", window.coords.file, "\n", sep="")
  trees.to.be.worked.with <- paste(basename(tree.file.root), ts.both.present, tree.fe, sep="")
  
  window.coords <- read.csv(window.coords.file, header = F, stringsAsFactors = F)
  window.coords[,2] <- as.numeric(window.coords[,2])
  files.expected <- window.coords[,1]
  # do all trees actually have a window?
  if(length(setdiff(trees.to.be.worked.with, files.expected))){
    for(missing in setdiff(trees.to.be.worked.with, files.expected)){
      cat("No genome position information found in file ", args$windowCoords, " for tree file ", missing, "\n", sep="")
    }
    quit(save="no", status=1)
  }
  
  # Try to find sensible start and end x-axis values.
  
  range <- max(window.coords[,2]) - min(window.coords[,2])
  increment <- range/nrow(window.coords)
  
  ews <- min(window.coords[,2]) - 0.45*increment
  lwe <- max(window.coords[,2]) + 0.45*increment
  
  # translate the first column to suffixes
  
  window.coords[,1] <- sapply(window.coords[,1], function(x) get.suffix(x, tree.file.root, tree.fe ))
  
  for(row.no in 1:nrow(window.coords)){
    if(window.coords[row.no,1] %in% suffixes){
      all.tree.info[[window.coords[row.no,1]]]$xcoord <- window.coords[row.no,2]
    }
  }
  
} else {
  if (verbose) cat("Attempting to read genome coordinates from file names...", sep="")
  
  regex <- "^\\D*([0-9]+)_to_([0-9]+).*$"
  
  starts <- sapply(suffixes, function(x) if(length(grep(regex, x))>0) as.numeric(sub(regex, "\\1", x)) else NA)
  ends <- sapply(suffixes, function(x) if(length(grep(regex, x))>0) as.numeric(sub(regex, "\\2", x)) else NA)
  
  if(any(is.na(starts)) | any(is.na(ends)) | length(starts)==0 | length(ends)==0){
    cat("Cannot obtain window coordinates from some file names\n")
    quit(save="no", status=1)
  }
  
  if (verbose) cat(" OK.\n")
  
  ews <- min(starts)
  lwe <- max(ends)
  
  middles <- (starts+ends)/2
  
  window.coords <- data.frame(suffix = suffixes, coordinate=middles)
  
  for(row.no in 1:nrow(window.coords)){
    all.tree.info[[window.coords[row.no,1]]]$xcoord <- window.coords[row.no,2]
  }
}

x.limits <- c(ews, lwe)

# Read in the IDs. Remove duplicates. Alphabeticise.

hosts <- scan(id.file, what="", sep="\n", quiet=TRUE)
hosts <- unique(hosts)
hosts <- hosts[order(hosts)]

num.hosts <- length(hosts)

if (num.hosts == 0) {
  stop(paste("No host IDs found in ", id.file, ". Quitting.\n", sep=""))
}

# Load the splits first, since these tables get reused

for(suffix in suffixes){
  
  splits.file.name <- all.tree.info[[suffix]]$splits.file.name
  
  if(!file.exists(splits.file.name)){
    cat("Subgraph file ",splits.file.name," does not exist.\n", sep="")
    quit(save="no", status=1)
  }
  
  splits.table <- read_csv(splits.file.name)
  
  colnames(splits.table) <- c("host", "subgraph", "tip")
  
  if (verbose) cat('Read subgraph file ',splits.file.name,'\n', sep="")
  
  if(!no.read.counts){
    splits.table$reads <- sapply(splits.table$tip, function(x) as.numeric(read.count.from.label(x, tip.regex)))
  } else {
    splits.table$reads <- 1
  }
  all.tree.info[[suffix]]$splits.table <- splits.table
  
  tree.file.name <- all.tree.info[[suffix]]$tree.file.name
  
  if(!file.exists(tree.file.name)){
    cat("Tree ",tree.file.name," does not exist.\n", sep="")
    quit(save="no", status=1)
  }
  
  pseudo.beast.import <- read.beast(tree.file.name)
  tree <- attr(pseudo.beast.import, "phylo")
  
  if (verbose) cat('Read tree file ',tree.file.name,'\n', sep="")
  
  all.tree.info[[suffix]]$tree <- tree
  
  blacklist.file.name <- all.tree.info[[suffix]]$blacklist.file.name
  blacklist <- vector()
  
  if(!is.null(blacklist.file.name)){
    if(file.exists(blacklist.file.name)){
      # There should already have been a warning if this file does not exist. It isn't fatal.
      
      blacklisted.tips <- read.table(blacklist.file.name, sep=",", stringsAsFactors = F, col.names="read")
      
      
      if (verbose) cat('Read blacklist file ',blacklist.file.name,'\n', sep="")
      
      if(nrow(blacklisted.tips)>0){
        blacklist <- sapply(blacklisted.tips, get.tip.no, tree=tree)
        all.tree.info[[suffix]]$blacklist <- blacklist
      }
    }
  }
  
  clade.results <- resolve.tree.into.host.clades(tree, hosts, tip.regex, blacklist, no.read.counts)
  
  all.tree.info[[suffix]]$clades.by.host      <- clade.results$clades.by.patients
  all.tree.info[[suffix]]$clade.mrcas.by.host <- clade.results$clade.mrcas.by.patient
  
}


# Calculate all the statistics apart from the subgraph proportions

pat.stats <- lapply(all.tree.info, function(x) calc.all.stats.in.window(x, hosts, tip.regex, verbose))
pat.stats <- bind_rows(pat.stats)

# If you're looking at absolutely nothing, I'm not drawing you any graphs

if(length(which(pat.stats$tips > 0))==0){
  stop("No patients from ID file ",id.file," present in any tree; stopping. Is your regex correct?\n", sep="")
  quit(save="no", status=1)
}

# The proportions of reads in each patient subgraph in each window

read.proportions <- lapply(setNames(suffixes, suffixes), function(y) lapply(setNames(hosts, hosts), function(x) get.read.proportions(x, y, all.tree.info[[y]]$splits.table)))

# Get the max split count over every window and patient (the exact number of columns depends on this)

max.splits <- max(sapply(read.proportions, function(x) max(sapply(x, function(y) length(y)  ) )))

read.prop.columns <- lapply(setNames(suffixes, suffixes), function(x){
  out <- lapply(setNames(hosts, hosts), function(y){
    window.props <- read.proportions[[x]]
    result <- window.props[[y]]
    
    if(is.na(result[1])){
      result <- rep(NA, max.splits)
      
    } else {
      if(length(result)<max.splits){
        result[(length(result)+1):max.splits] <- 0
      }
    }
    
    result
  })
  out <- as_tibble(do.call(rbind, out))
  colnames(out) <- paste("prop.gp.",seq(1,max.splits),sep="")
  out
})


read.prop.columns <- bind_rows(read.prop.columns)

pat.stats         <- cbind(pat.stats, read.prop.columns)

pat.stats$reads   <- as.numeric(pat.stats$reads)
pat.stats$tips    <- as.numeric(pat.stats$tips)

# Output the pat.stats table to file

tmp	<- file.path(paste(output.root,"patStatsFull.csv",sep=""))
if (verbose) cat("Writing output to file ",tmp,"...\n",sep="")
write.csv(pat.stats, tmp, quote = F, row.names = F)

pat.stats$prop.reads.largest.subtree <- pat.stats$prop.gp.1

# Output a summary of the pat.stats table to file
if (recomb.files.exist) { 
  tmp <- as_tibble(pat.stats) %>%
    filter(reads>0) %>%
    select(host.id, tips, reads, subgraphs, clades, overall.rtt, largest.rtt, max.pat.distance, prop.reads.largest.subtree, max.branch.length, subgraph.mean.pat.distance, global.mean.pat.distance, recombination.metric)
} else {
  tmp <- as_tibble(pat.stats) %>%
    filter(reads>0) %>%
    select(host.id, tips, reads, subgraphs, clades, overall.rtt, largest.rtt, max.pat.distance, prop.reads.largest.subtree, max.branch.length, subgraph.mean.pat.distance, global.mean.pat.distance)
}

pat.stats.summary <- tmp %>%
  group_by(host.id) %>%
  summarise_all(list(mean, sd), na.rm = TRUE)
   
tmp <- file.path(paste(output.root,"patStatsSummary.csv",sep=""))
if (verbose) cat("Writing output to file ",tmp,"...\n",sep="")
write.csv(pat.stats.summary, tmp, quote = F, row.names = F)

# Draw the graphs

tmp <- paste(output.root,"patStats.pdf",sep="")
if (verbose) cat("Plotting to file ",tmp,"...\n",sep="")

# Set up the boundaries of each window's region on the x-axis

xcoords <- unique(pat.stats$xcoord)

range <- max(xcoords) - min(xcoords)
increment <- range/length(xcoords)

ews <- min(xcoords) - 0.45*increment
lwe <- max(xcoords) + 0.45*increment

x.limits <- c(ews, lwe)

if(length(xcoords)>1){
  xcoords <- xcoords[order(xcoords)]
  
  missing.window.data <- find.gaps(xcoords, x.limits)
  
  xcoords <- missing.window.data$x.coordinates
  regular.gaps <- missing.window.data$regular.gaps
  rectangles.for.missing.windows <- missing.window.data$rectangles.for.missing.windows
  bar.width <- missing.window.data$width
  
  produce.pdf.graphs(tmp, pat.stats, hosts, xcoords, x.limits, rectangles.for.missing.windows, bar.width, regular.gaps, readable.coords = !is.null(window.coords.file), verbose = verbose)
} else {
  cat("Only one window present; not plotting summary statistics")
}
