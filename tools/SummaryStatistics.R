list.of.packages <- c("argparse","phytools", "dplyr", "ggplot2", "ggtree", "reshape", "dtplyr", "gtable", "grid", "gridExtra", "RColorBrewer", "scales", "pegas")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Please run PackageInstall.R to continue\n")
  quit(save='no')
}

suppressMessages(require(ape, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phytools, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ggplot2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ggtree, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(reshape, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gtable, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(grid, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gridExtra, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(RColorBrewer, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(scales, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(dplyr, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(dtplyr, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(data.table, quietly=TRUE, warn.conflicts=FALSE))

#	constants

tree.fe <- ".tree" 
csv.fe <- ".csv"

# This function grabs the suffix which is used to uniquely identify each window through all input files

get.suffix <- function(file.name, prefix, extension){

  file.name <- basename(file.name)
  prefix <- basename(prefix)

  substr(file.name, nchar(prefix)+1, nchar(file.name)-nchar(extension))
}

command.line <- T

if (command.line) {
  require(argparse, quietly=TRUE, warn.conflicts=FALSE)
  
  arg_parser = ArgumentParser(description="Summarise patient statistics from each phylogeny across all windows.")
  
  arg_parser$add_argument("-c", "--noReadCounts", action="store_true", default=FALSE, help="If present, tip labels have no read counts and each tip is assumed to represent a single read") 
  arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                          help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
  arg_parser$add_argument("-b", "--blacklists", help="An optional file path and initial string identifying blacklist files (file extension must be .csv).")
  arg_parser$add_argument("-w", "--windowCoords", action="store", help="The path for a .csv file describing the position in the genome which each tree file represents (e.g. the midpoint of a window); used for the x-axis in output. The first column in the file should the tree file name, the second the coordinate. If not given, the script will attempt to obtain these from the file names (which should work if the input is from the phyloscanner pipeline).")
  arg_parser$add_argument("idFile", action="store", help="A file containing a list of the IDs of all the patients to calcualte and display statistics for.")
  arg_parser$add_argument("treeFiles", action="store", help="A file path and initial string identifying all tree files (processed tree files output by SplitPatientsToSubgraphs.R file extension must be .tree).")
  arg_parser$add_argument("splitsFiles", action="store",help="A file path and initial string identifying all splits files (file extension must be .csv).")
  arg_parser$add_argument("outputBaseName", action="store", help="A path and string to begin the names of all output files.")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.")
  
  # Read in the arguments
  
  args <- arg_parser$parse_args()
  
  id.file <- args$idFile
  script.dir <- args$scriptdir

  # The scripts are necessary here
  
  source(file.path(script.dir, "TreeUtilityFunctions.R"))
  source(file.path(script.dir, "ParsimonyReconstructionMethods.R"))
  source(file.path(script.dir, "SummariseTrees_funcs.R"))
  
  tip.regex <- args$tipRegex
  tree.file.root <- args$treeFiles
  splits.file.root <- args$splitsFiles
  blacklist.file.root <- args$blacklists  
  output.root <- args$outputBaseName
  window.coords.file <- args$windowCoords

  no.read.counts <- args$noReadCounts
  
  # Find the input files
  
  tree.files <- sort(list.files.mod(dirname(tree.file.root), pattern=paste('^',basename(tree.file.root),".*",tree.fe,sep=""), full.names=TRUE))
  
  splits.files <- sort(list.files.mod(dirname(splits.file.root), pattern=paste('^',basename(splits.file.root),".*",csv.fe,sep=""), full.names=TRUE))	  
  blacklist.files <- NULL
  if(!is.null(blacklist.file.root)){
    blacklist.files <- sort(list.files.mod(dirname(blacklist.file.root), pattern=paste('^',basename(blacklist.file.root),".*\\.csv",sep=""), full.names=TRUE))
  }
  
  # Get the suffixes
  
  tree.suffixes	<- sapply(tree.files, function(x) get.suffix(x, tree.file.root, tree.fe))
  
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
  
  # Check which suffixes from the current list have a blacklist file and vice versa
  
  if(!is.null(blacklist.file.root)){
    blacklist.suffixes <- sapply(blacklist.files, function(x) get.suffix(x, blacklist.file.root, csv.fe))


    # the only remaining tree files will have a splits file, no need to check the splits files separately
    bt.both.present <- intersect(ts.both.present, blacklist.suffixes)	
    if(length(bt.both.present) < length(ts.both.present)){
      missing.suffixes <- setdiff(ts.both.present, bt.both.present)
      no.partner.files <- paste(tree.file.root, missing.suffixes, tree.fe, sep="")
      
      for(nps in no.partner.files){
        warning(paste("No blacklist file found for tree file ",nps,"; assuming no tips were blacklisted", sep=""))
      }
    } 
    if(length(bt.both.present) < length(blacklist.suffixes)){
      missing.suffixes <- setdiff(blacklist.suffixes, bt.both.present)
      no.partner.files <- paste(blacklist.file.root, missing.suffixes, csv.fe, sep="")
      
      for(nps in no.partner.files){
        warning(paste("Tree or subgraph file missing for blacklist file ",nps,"; ignoring this file", sep=""))
      }
    } 
  }
  
  # From now on we work with suffixes only
  
  suffixes <- ts.both.present
  
  if(!is.null(window.coords.file)){
    cat("Reading genome coordinates from file ", window.coords.file, "\n", sep="")
    trees.to.be.worked.with <- paste(basename(tree.file.root), ts.both.present, tree.fe, sep="")
    
    
    window.coords <- read.csv(window.coords.file, header = F, stringsAsFactors = F)
    window.coords[,2] <- as.numeric(window.coords[,2])
    files.expected <- window.coords[,1]
    # do all trees actually have a window?
    if(length(setdiff(trees.to.be.worked.with, files.expected))){
      for(missing in setdiff(trees.to.be.worked.with, files.expected)){
        cat("No genome position information found in file ", args$windowCoords, " for tree file ", missing, "\n", sep="")
      }
      quit(save="no")
    }
    
    # Try to find sensible start and end x-axis values. Untested.
    
    range <- max(window.coords[,2]) - min(window.coords[,2])
    increment <- range/nrow(window.coords)
    
    ews <- min(window.coords[,2]) - 0.45*increment
    lwe <- max(window.coords[,2]) + 0.45*increment
    
    # translate the first column to suffixes
    
    window.coords[,1] <- sapply(window.coords[,1], function(x) get.suffix(x, tree.file.root, tree.fe ))
    colnames(window.coords) <- c("suffix", "coordinate")
    
    
  } else {
    cat("Attempting to read genome coordinates from file names...", sep="")
    
    regex <- "^\\D*([0-9]+)_to_([0-9]+).*$"
    
    starts <- sapply(suffixes, function(x) if(length(grep(regex, x))>0) as.numeric(sub(regex, "\\1", x)) else NA)
    ends <- sapply(suffixes, function(x) if(length(grep(regex, x))>0) as.numeric(sub(regex, "\\2", x)) else NA)
  
    if(any(is.na(starts)) | any(is.na(ends)) | length(starts)==0 | length(ends)==0){
      cat("Cannot obtain window coordinates from some file names\n")
      quit(save="no")
    }
    
    cat(" OK.\n")
    
    ews <- min(starts)
    lwe <- max(ends)
    
    middles <- (starts+ends)/2

    window.coords <- data.frame(suffix = suffixes, coordinate=middles)
  }
  
  window.coords <- as.list(setNames(window.coords$coordinate, window.coords$suffix ))

} else {

  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20161013")
  tree.file.root <- "AllTrees/RAxML_bestTree.InWindow_"
  blacklist.file.root <- "Blacklists/Duals blacklists/DualBlacklist.InWindow_"
  splits.file.root <- "SubtreeFiles_s/Subtrees_s_run20161013_inWindow_"
  tip.regex <- "^(.*)-[0-9].*_read_([0-9]+)_count_([0-9]+)$"
  id.file <- "patientIDListShortWDubious.txt"
  
  setwd("/Users/twoseventwo/Documents/Croucher alignments/")
  id.file <- "patients.txt"
  tree.file.root <- "ProcessedTree_s_pneumo_"
  splits.file.root <- "subgraphs_s_pneumo_"
  blacklist.file.root <- NULL
  window.coords.file <- "samplecoords.csv"
  tip.regex <- "^(ARI-[0-9][0-9][0-9][0-9]_[A-Z]*)_[0-9][0-9][0-9][0-9]_[0-9]_[0-9][0-9]?$"
  no.read.counts <- T
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools"
  
  setwd("/Users/twoseventwo/Downloads/bams/")
  id.file <- "PatientIDfile.txt"
  tree.file.root <- "ProcessedTree_s_SimulatedDataRefs_InWindow_"
  splits.file.root <- "subgraphs_s_SimulatedDataRefs_InWindow_"
  blacklist.file.root <- "FinalBlacklist.InWindow_"
  no.read.counts <- F
  tip.regex <- "^(.*)_read_([0-9]+)_count_([0-9]+)$"
}

# Align the graph output

AlignPlots <- function(...) {
  LegendWidth <- function(x) x$grobs[[15]]$grobs[[1]]$widths[[4]]
  
  plots.grobs <- lapply(list(...), ggplotGrob)
  
  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]], unit(abs(diff(c(LegendWidth(x), max.legends.width))), "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}

# Add coloured rectangles to the graphs to represent areas with no coverage

add.no.data.rectangles <- function(graph, rectangle.coords, log = F){
  y.limits <- ggplot_build(graph)$layout$panel_ranges[[1]]$y.range
  
  if(nrow(rectangle.coords)>0){
    for(rect.no in seq(1, nrow(rectangle.coords))){
      if(log){
        graph <- graph + annotate("rect", xmin=rectangle.coords$start[rect.no], xmax = rectangle.coords$end[rect.no], 
                                  ymin=10^(y.limits[1]), ymax=10^(y.limits[2]), fill=rectangle.coords$colour[rect.no], alpha=0.5)
      } else {
        graph <- graph + annotate("rect", xmin=rectangle.coords$start[rect.no], xmax = rectangle.coords$end[rect.no], 
                                  ymin=y.limits[1], ymax=y.limits[2], fill=rectangle.coords$colour[rect.no], alpha=0.5)
      }
    }
  }
  
  return(graph)
}

# Calculate the coordinates to build the rectangles for. missing.coords are x-coordinates of windows with no coverage. all.coords is the
# full list of window coordinates.

form.rectangles <- function(missing.coords, all.coords, colour = "grey"){
  if(length(missing.coords)>0){
    gap <-  unique(all.coords[2:length(all.coords)] - all.coords[1:length(all.coords)-1])
    if(length(gap)>1){
      cat("Fatal error drawing rectangles for missing data")
      quit(save="no")
    }
    temp.starts <- missing.coords - gap/2
    temp.ends <- missing.coords + gap/2
    
    final.starts <- setdiff(temp.starts, temp.ends)
    final.ends <- setdiff(temp.ends, temp.starts)
    
    return(data.frame(starts = final.starts, ends = final.ends, colour=colour))
  } 
  stop("No missing coordinates given")
}

unfactorDataFrame <- function(x) {
  x <- data.frame(lapply(x, as.character), stringsAsFactors = F)
  x <- data.frame(lapply(x, function(y) type.convert(y, as.is=T)),
                  stringsAsFactors = F)
}

# Collects a variety of statistics about a single patient in a single tree

calc.subtree.stats <- function(id, suffix, tree, tips.for.patients, splits.table, no.read.counts, verbose = F){
  if(verbose) cat("Calculating detailed statistics for patient ",id,".\n", sep="")
  
  subgraphs <- length(unique(splits.table$subgraph[which(splits.table$patient==id)]))
  
  all.tips <- tips.for.patients[[id]]
  
  subtree.all <- NULL
  
  if(length(all.tips)>1){
    
    subtree.all <- drop.tip(tree, tip=tree$tip.label[!(tree$tip.label %in% all.tips)])
    
    # just in case pruning changes the tip order
    
    if(!no.read.counts){
      reads.per.tip <- sapply(subtree.all$tip.label, function(x) as.numeric(read.count.from.label(x, tip.regex)))
    } else {
      reads.per.tip <- rep(1, length(subtree.all$tip.label))
    }
    
    names(reads.per.tip) <- subtree.all$tip.label
    
    overall.rtt <- calcMeanRootToTip(subtree.all, reads.per.tip)
    
    if(length(subtree.all$tip.label)>2){
      unrooted.subtree <- unroot(subtree.all)
      max.branch.length <- max(unrooted.subtree$edge.length)
    } else {
      
      # ape won't unroot two-tip trees
      max.branch.length <- sum(subtree.all$edge.length)
    }
    
    pat.distances <- cophenetic(subtree.all)
    max.pat.distance <- max(pat.distances)
    
    branch.to.pat.ratio <- max.branch.length/max.pat.distance
    
    if(subgraphs==1){
      mean.pat.distance <- mean(pat.distances[upper.tri(pat.distances)])
      largest.rtt <- overall.rtt
    } else {
      relevant.reads <- splits.table[which(splits.table$patient==id),]
      splits <- unique(relevant.reads$subgraph)
      reads.per.split <- sapply(splits, function(x) sum(splits.table$reads[which(splits.table$subgraph==x)] ) )
      
      winner <- splits[which(reads.per.split==max(reads.per.split))]

      winner.tips <- relevant.reads$tip[which(relevant.reads$subgraph==winner)]
      if(length(winner.tips)==1){
        largest.rtt <- 0
        mean.pat.distance <- NA
      } else {
        subtree <- drop.tip(tree, tip=tree$tip.label[!(tree$tip.label %in% winner.tips)])
        
        # just in case pruning changes the tip order
        
        if(!no.read.counts){
          reads.per.tip <- sapply(subtree$tip.label, function(x) as.numeric(read.count.from.label(x, tip.regex)))
        } else {
          reads.per.tip <- rep(1, length(subtree$tip.label))
        }
        
        names(reads.per.tip) <- subtree$tip.label
        
        largest.rtt <- calcMeanRootToTip(subtree, reads.per.tip)
        pat.distances <- cophenetic(subtree)
        mean.pat.distance <- mean(pat.distances[upper.tri(pat.distances)])
      }
    }
    
  } else {
    
    # A patient with zero or one tips has no branches and patristic distances are meaningless.

    max.branch.length <- NA
    max.pat.distance <- NA
    branch.to.pat.ratio <- NA
    mean.pat.distance <- NA
    
    if(length(all.tips)==1){
      # A patient with one tip does have a root
      overall.rtt <- 0
      largest.rtt <- 0
    } else {
      # An absent patient has no root or tips
      overall.rtt <- NA
      largest.rtt <- NA
    }
  }
  return(list(overall.rtt = overall.rtt, largest.rtt = largest.rtt, max.branch.length = max.branch.length, max.pat.distance = max.pat.distance,
              branch.to.pat.ratio  = branch.to.pat.ratio, mean.pat.distance = mean.pat.distance))
}

# Calculates all statistics (apart from read proportions) for all patients in a given window

calc.all.stats.in.window <- function(suffix, verbose = F){
  
  # Make the file names anew from the suffixes (probably the cleanest way to do this)
  
  tree.file.name <- paste(tree.file.root, suffix, tree.fe, sep="")
  
  if(!is.null(blacklist.file.root)){
    blacklist.file.name <- paste(blacklist.file.root, suffix, csv.fe, sep="")
  } else {
    blacklist.file.name = NULL
  }
  
  # The following should not happen in practice; we've already checked they exist
  
  if(!file.exists(tree.file.name)){
    cat("Tree file ",tree.file.name," does not exist.\n", sep="")
    quit(save="no")
  }
  
  # Read the tree
  
  pseudo.beast.import <- read.beast(tree.file.name)
  tree <- attr(pseudo.beast.import, "phylo")
  
  cat('Read tree file ',tree.file.name,'\n', sep="")
  
  # Load the blacklist if it exists
  
  blacklist <- vector()
  
  if(!is.null(blacklist.file.name)){
    if(file.exists(blacklist.file.name)){
      # There should already have been a warning if this file does not exist. It isn't fatal.
      
      blacklisted.tips <- read.table(blacklist.file.name, sep=",", stringsAsFactors = F, col.names="read")
      
      cat('Read blacklist file ',blacklist.file.name,'\n', sep="")
      
      if(nrow(blacklisted.tips)>0){
        blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
      }
    }
  }
  
  # Get the splits
  
  splits.table <- all.splits.table[[suffix]]
  
  # Find the clades
  clade.results <- resolveTreeIntoPatientClades(tree, ids, tip.regex, blacklist, no.read.counts)
  clade.mrcas.by.patient <- clade.results$clade.mrcas.by.patient
  all.clades.by.patient <- clade.results$clades.by.patient
  
  # A vector of patients for each tip
  
  patients.for.tips <- sapply(tree$tip.label, function(x) patient.from.label(x, tip.regex))
  patients.for.tips[blacklist] <- NA
  
  # If no patients from the input file are actually here, give a warning
  
  patients.present <- intersect(ids, unique(patients.for.tips))
  if(length(patients.present)==0){
    warning(paste("No patients from file ",id.file," appear in tree ",tree.file.name,"\n",sep=""))
  }
  
  # A list of tips for each patient 
  
  tips.for.patients <- lapply(setNames(ids, ids), function(x) tree$tip.label[which(patients.for.tips==x)])
  
  # Make the data.table
  
  window.table <- data.table(id=ids)
  window.table <- window.table[, file.suffix := suffix]
  window.table <- window.table[, xcoord := window.coords[[suffix]]]
  window.table <- window.table[, tips :=  sapply(ids, function(x) length(tips.for.patients[[x]]))]
  window.table <- window.table[, reads :=  sapply(ids, function(x){
    if(length(tips.for.patients[[x]])==0){
      return(0)
    } else {
      if(!no.read.counts){
        return(sum(sapply(tips.for.patients[[x]], function(y) as.numeric(read.count.from.label(y, tip.regex)))))
      } else {
        return(length(tips.for.patients[[x]]))
      }
    }
  } 
  )]
  window.table <- window.table[, subgraphs :=  sapply(ids, function(x) length(unique(splits.table[which(splits.table$patient==x),]$subgraph)))]
  window.table <- window.table[, clades := sapply(ids, function(x) length(all.clades.by.patient[[x]]))  ]
  
  new.cols <- sapply(ids, function(x) calc.subtree.stats(x, suffix, tree, tips.for.patients, splits.table, no.read.counts, verbose))
  new.cols <- as.data.table(t(new.cols))
  new.cols <- sapply(new.cols, as.numeric)
  window.table <- cbind(window.table, new.cols) 
  
  window.table
}

# Get the proportion of reads from a given patient in each of its subgraphs in a single tree

get.read.proportions <- function(id, suffix, splits.table){
  this.pat.splits <- splits.table[which(splits.table$patient==id),]
  
  if(nrow(this.pat.splits)>0){
    this.pat.reads.by.split <- aggregate(this.pat.splits$reads, by=list(Category=this.pat.splits$subgraph), sum)
    
    this.pat.reads.by.split <- this.pat.reads.by.split[order(this.pat.reads.by.split$x, decreasing = T),]
    return(this.pat.reads.by.split$x/sum(this.pat.reads.by.split$x))

  } else {
    return(NA)
  }
}


# Read in the IDs. Remove duplicates. Alphabeticise.

ids <- scan(id.file, what="", sep="\n", quiet=TRUE)
ids <- unique(ids)
ids <- ids[order(ids)]

num.ids <- length(ids)

if (num.ids == 0) {
  cat(paste("No IDs found in ", id.file, ". Quitting.\n", sep=""))
  quit("no", 1)
}

# Load the splits first, since these tables get reused

all.splits.table <- lapply(setNames(suffixes, suffixes), function(x){
  splits.file.name <- paste(splits.file.root, x, csv.fe, sep="")
  
  if(!file.exists(splits.file.name)){
    cat("Subgraph file ",splits.file.name," does not exist.\n", sep="")
    quit(save="no")
  }
  
  splits.table <- fread(splits.file.name)
  
  # todo get rid of this
  colnames(splits.table) <- c("patient", "subgraph", "tip")
  
  cat('Read subgraph file ',splits.file.name,'\n', sep="")
  
  if(!no.read.counts){
    splits.table$reads <- sapply(splits.table$tip, function(x) as.numeric(read.count.from.label(x, tip.regex)))
  } else {
    splits.table$reads <- 1
  }
  splits.table$suffix <- x
  splits.table
})

# Calculate all the statistics apart from the subgraph proportions

pat.stats <- lapply(suffixes, function(x) calc.all.stats.in.window(x, F))
pat.stats <- rbindlist(pat.stats)

# If you're looking at absolutely nothing, I'm not drawing you any graphs

if(length(which(pat.stats$tips > 0))==0){
  cat("No patients from ID file ",id.file," present in any tree; stopping. Is your regex correct?\n", sep="")
  quit(save="no")
}

# The proportions of reads in each patient subgraph in each window

read.proportions <- lapply(setNames(suffixes, suffixes), function(y) lapply(setNames(ids, ids), function(x) get.read.proportions(x, y, all.splits.table[[y]])))

# Get the max split count over every window and patient (the exact number of columns depends on this)

max.splits <- max(sapply(read.proportions, function(x) max(sapply(x, function(y) length(y)  ) )))

read.prop.columns <- lapply(setNames(suffixes, suffixes), function(x){
  out <- lapply(setNames(ids, ids), function(y){
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
  out <- as.data.table(do.call(rbind, out))
  colnames(out) <- paste("prop.gp.",seq(1,max.splits),sep="")
  out
})
  

read.prop.columns <- rbindlist(read.prop.columns)
  
pat.stats <- cbind(pat.stats, read.prop.columns)

# Output the pat.stats table to file

tmp	<- file.path(paste(output.root,"_patStatsFull.csv",sep=""))
cat("Writing output to file ",tmp,"...\n",sep="")
write.csv(pat.stats, tmp, quote = F, row.names = F)

pat.stats$prop.reads.largest.subtree <- pat.stats$prop.gp.1

mean.na.rm <- function(x) mean(x, na.rm = T)

# Output a summary of the pat.stats table to file

pat.stats.temp <- pat.stats[which(pat.stats$reads>0), c("id", "tips", "reads", "subgraphs", "clades", "overall.rtt", "largest.rtt", "max.pat.distance",
                                                        "prop.reads.largest.subtree", "max.branch.length", "mean.pat.distance", "branch.to.pat.ratio")] 


by.patient <- pat.stats.temp %>% group_by(id)
pat.stats.summary <-
  as.data.frame(by.patient %>% summarise_each(funs(mean.na.rm)))

tmp <- file.path(paste(output.root,"_patStatsSummary.csv",sep=""))
cat("Writing output to file ",tmp,"...\n",sep="")
write.csv(pat.stats.summary, tmp, quote = F, row.names = F)

# Draw the graphs

tmp <- paste(output.root,"_patStats.pdf",sep="")
cat("Plotting to file ",tmp,"...\n",sep="")

# Set up the boundaries of each window's region on the x-axis

xcoords <- unique(pat.stats$xcoord)
xcoords <- xcoords[order(xcoords)]

# Try to detect if the windows are evenly spaced. If they aren't, skip the rectangle-drawing, it's too complex.

gaps <- xcoords[2:length(xcoords)]-xcoords[1:length(xcoords)-1]

if(length(unique(gaps))==1){
  # perfect
  regular.gaps <- T
  xcoords.reg <- xcoords
  missing.window.rects <- NULL
  
} else {
  # if the length of every gap is a multiple of the smallest such gap length, it suggests missing windows
  smallest.gap <- min(gaps)
  if(length(which(gaps/smallest.gap != floor(gaps/smallest.gap)))==0){
    regular.gaps <- T
    xcoords.reg <- seq(min(xcoords), max(xcoords), smallest.gap)
    missing.window.rects <- form.rectangles(setdiff(xcoords.reg, xcoords), xcoords.reg, "grey20")
  } else {
    regular.gaps <- F
    missing.window.rects <- NULL
  }
}

range <- max(xcoords.reg) - min(xcoords.reg)
bar.width.5 <- range/(1.5*(length(xcoords.reg)+1))

pdf(file=tmp, width=8.26772, height=11.6929)

for (i in seq(1, num.ids)) {
  
  patient <- ids[i]
  
  this.pat.stats <- pat.stats[which(pat.stats$id==patient),]
  
  this.pat.stats <- this.pat.stats[order(this.pat.stats$xcoord),]
  this.pat.stats$reads <- as.numeric(this.pat.stats$reads)
  this.pat.stats$tips <- as.numeric(this.pat.stats$tips)
  
  this.pat.stats.temp <- this.pat.stats[which(this.pat.stats$reads>0),]
  if(nrow(this.pat.stats.temp)>0){
    cat("Drawing graphs for patient ",patient,"\n", sep="")
    
    #Get the zero runs for the grey rectangles
    
    missing.read.rects <- NULL

    if(regular.gaps & length(which(this.pat.stats$reads==0))){
      missing.read.rects <- form.rectangles(this.pat.stats$xcoord[which(this.pat.stats$reads==0)], xcoords.reg, "grey")
    }
    
    missing.rects <- rbind(missing.window.rects, missing.read.rects)

    
    #graph 1: read count and tip count
    
    this.pat.stats.1col <- melt(this.pat.stats.temp[,c("xcoord","id","tips","reads")], id.vars=c("id", "xcoord"))
    
    graph.1 <- ggplot(this.pat.stats.1col, aes(x=xcoord, y=value, col=variable))
    
    graph.1 <- graph.1 + geom_point(na.rm=TRUE) +
      theme_bw() + 
      scale_y_log10(breaks=c(1,10,100,1000,10000)) +
      ylab("Count") +
      xlab("Window centre") +
      scale_x_continuous(limits=c(ews, lwe)) +
      scale_color_discrete(name="Variable", labels=c("Tips", "Reads")) + 
      theme(text = element_text(size=7))
    
    if(regular.gaps & !is.null(missing.rects)){
      graph.1 <- add.no.data.rectangles(graph.1, missing.rects,  TRUE)
    }
    
    
    #graph 2: subgraph and clade counts
    
    this.pat.stats.1col <- melt(this.pat.stats.temp[,c("xcoord", "id", "subgraphs", "clades")], id.vars=c("id", "xcoord"))
    
    graph.2 <- ggplot(this.pat.stats.1col, aes(x=xcoord, y=value))
    
    graph.2 <- graph.2 +
      geom_point(aes(shape=variable, size=variable), na.rm=TRUE) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Count") +
      xlab("Window centre") +
      scale_x_continuous(limits=c(ews, lwe)) +
      scale_shape_manual(values=c(1,19), name="Variable", labels=c("Subgraphs", "Clades")) +  
      scale_size_manual(values=c(2,1), name="Variable", labels=c("Subgraphs", "Clades")) +		
      scale_color_discrete(name="Variable", labels=c("Subgraphs", "Clades")) + 
      theme(text = element_text(size=7))
    
    if(max(this.pat.stats.1col$value)==1){
      graph.2 <- graph.2 + scale_y_continuous(breaks = c(0,1))
    } else if(max(this.pat.stats.1col$value)==2){
      graph.2 <- graph.2 + expand_limits(y=0)+ scale_y_continuous(breaks = c(0,1,2)) 
    } else {
      graph.2 <- graph.2 + expand_limits(y=0) + scale_y_continuous(breaks = pretty_breaks()) 
    }
    
    if(regular.gaps & !is.null(missing.rects)){
      graph.2 <- add.no.data.rectangles(graph.2, missing.rects)
    }
    
    this.pat.stats.1col <- melt(this.pat.stats.temp[,c("xcoord","id","overall.rtt","largest.rtt")], id.vars=c("id", "xcoord"))
   
    # graph 3: root-to-tip distances for all reads and largest subgraph
    
    graph.3 <- ggplot(this.pat.stats.1col, aes(x=xcoord, y=value))
    
    graph.3 <- graph.3 +
      geom_point(aes(shape=variable, size=variable), na.rm=TRUE) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Mean root-to-tip-distance\n(read-weighted)") +
      xlab("Window centre") +
      scale_x_continuous(limits=c(ews, lwe)) +
      expand_limits(y=0) + 
      scale_color_discrete(name="Tip set", labels=c("All", "Largest subgraph")) + 
      scale_shape_manual(values=c(1,19), name="Tip set", labels=c("All", "Largest subgraph")) +
      scale_size_manual(values=c(2,1), name="Tip set", labels=c("All", "Largest subgraph")) +
      theme(text = element_text(size=7))
    
    if(regular.gaps & !is.null(missing.rects)){
      graph.3 <- add.no.data.rectangles(graph.3, missing.rects)
    }
    
    # graph 4: largest patristic distance in largest subgraph
    
    graph.4 <- ggplot(this.pat.stats.temp, aes(x=xcoord, y=max.pat.distance))
    
    graph.4 <- graph.4 +
      geom_point(alpha=0.5, na.rm=TRUE) +
      theme_bw() + 
      ylab("Largest patristic distance in\nlargest connected subgraph") +
      xlab("Window centre") +
      scale_x_continuous(limits=c(ews, lwe)) +
      expand_limits(y=0) +
      theme(text = element_text(size=7))
    
    if(regular.gaps & !is.null(missing.rects)){
      graph.4 <- add.no.data.rectangles(graph.4, missing.rects)
    }
    

    proportion.column.names <- colnames(pat.stats)[which(substr(colnames(pat.stats), 1, 8)=="prop.gp.")]
    
  
    # only columns with nonzero entries should be kept
    proportion.columns <- pat.stats[which(pat.stats$id==patient), proportion.column.names, with=F]
    
    sums <- colSums(proportion.columns, na.rm=T)
    proportion.column.names <- proportion.column.names[which(sums>0)]
    
    splits.props <- pat.stats[which(pat.stats$id==patient),c("id", "xcoord", proportion.column.names), with=F]
    
    splits.props.1col <- melt(splits.props, id=c("id", "xcoord"))
    splits.props.1col <- splits.props.1col[!is.na(splits.props.1col$value),]
    
    splits.props.1col$ngroup <- sapply(splits.props.1col$variable, function(x) which(levels(splits.props.1col$variable)==x))
    
    splits.props.1col$fgroup <- as.factor(splits.props.1col$ngroup)
    
    colourCount = length(unique(splits.props.1col$ngroup))
    getPalette = colorRampPalette(brewer.pal(9, "Greens"))
    
    # graph 5: read proportions in each subgraph
    
    graph.5 <- ggplot(splits.props.1col, aes(x=xcoord, weight=value, fill=reorder(fgroup, rev(order(splits.props.1col$ngroup)))))
  
    graph.5 <- graph.5 +
      geom_bar(width=bar.width.5, colour="black", size=0.25) +
      theme_bw() + 
      ylab("Proportion of reads\nin discrete subraphs") +
      xlab("Window centre") +
      scale_x_continuous(limits=c(ews, lwe)) +
      scale_fill_manual(values = getPalette(colourCount)) +
      theme(text = element_text(size=7)) + 
      guides(fill = guide_legend(title = "Subgraph rank\n(by tip count)", keywidth = 1, keyheight = 0.4))
    
    if(regular.gaps & !is.null(missing.rects)){
      graph.5 <- add.no.data.rectangles(graph.5, missing.rects)
    }
    
    # graph 6: longest branch to largest patristic distance ratios
    
    graph.6 <- ggplot(this.pat.stats.temp, aes(x=xcoord, y=branch.to.pat.ratio))
    
    graph.6 <- graph.6 +
      geom_point(alpha = 0.5, na.rm=TRUE) +
      theme_bw() + 
      ylab("Ratio of longest branch to greatest\n patristic distance between tips") +
      xlab("Window centre") +
      scale_x_continuous(limits=c(ews, lwe)) +
      expand_limits(y=0) +
      #      scale_color_discrete(name="Tip set", labels=c("Longest branch", "Greatest patristic distance")) + 
      theme(text = element_text(size=7))
    
    if(regular.gaps & !is.null(missing.rects)){
      graph.6 <- add.no.data.rectangles(graph.6, missing.rects)
    }
  
    all.plots <- AlignPlots(graph.1, graph.2, graph.3, graph.4, graph.5, graph.6)

    if(i!=1){
      grid.newpage()
    }
    
    pushViewport(viewport(layout = grid.layout(7, 1, heights = unit(c(0.25, rep(1,6)), "null") )))
    grid.text(patient, gp=gpar(fontsize=20), vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
    for(plot.no in 1:6){
      plot <- all.plots[[plot.no]]
      plot$vp = viewport(layout.pos.col = 1, layout.pos.row = plot.no + 1)
      grid.draw(plot)
    }
  
  } else {
    cat("Skipping graphs for patient ",patient," as no reads are present and not blacklisted\n", sep="")
  }
}

dev.off()