list.of.packages <- c("argparse", "data.table", "ape", "ff", "phangorn", "ggtree", "phytools", "scales", "RColorBrewer", "gtable", "grid", "gridExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Please run PackageInstall.R to continue\n")
  quit(save="no", status=1)
}

suppressMessages(require(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(data.table, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ape, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phangorn, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ggtree, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phytools, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(scales, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(RColorBrewer, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gtable, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(grid, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gridExtra, quietly=TRUE, warn.conflicts=FALSE))

command.line <- T

if(command.line){

  arg_parser		     <- ArgumentParser()

  arg_parser$add_argument("tree", action="store", help="A path and string that begins all the tree file names.")
  arg_parser$add_argument("outputString", action="store", help="This string identifies all output files.")

  arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
  arg_parser$add_argument("-m", "--multifurcationThreshold", help="If specified, short branches in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. If a number, this number will be used as the threshold. If 'g', it will be guessed from the branch lengths (use this only if you've checked by eye that the tree does indeed have multifurcations).")
  
  arg_parser$add_argument("-b", "--blacklist", action="store", help="A path and string that begins all the file names for pre-existing blacklist files.")
  
  # General, bland options
  
  arg_parser$add_argument("-od", "--outputDir", action="store", help="All output will be written to this directory. If absent, current working directory.")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.")
  arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what the script is doing.")
  arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. Three capture groups: patient ID, read ID, and read count; if the last is missing then read information will not be used. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
  arg_parser$add_argument("-y", "--fileNameRegex", action="store", default="^\\D*([0-9]+)_to_([0-9]+).*$", help="Regular expression identifying window coordinates. Two capture groups: start and end. If absent, input will be assumed to be from the phyloscanner pipeline.")
  arg_parser$add_argument("-tfe", "--treeFileExtension", action="store", default="tree", help="The file extension for tree files (default tree).")
  arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default csv).")
  arg_parser$add_argument("-pw", "--pdfwidth", action="store", default=50, help="Width of tree pdf in inches.")
  arg_parser$add_argument("-ph", "--pdfrelheight", action="store", default=0.15, help="Relative height of tree pdf.")
  
  # Normalisation options
  
  arg_parser$add_argument("-nr", "--normRefFileName", action="store", help="File name of reference for window normalising constants. If absent, no normalisation will be performed.")
  arg_parser$add_argument("-nv", "--normVar", action="store", default="MEDIAN_PWD", help="Column name in reference table that is to be used as normalising constant.")
  arg_parser$add_argument("-ns", "--normStandardise", action="store_true", default=FALSE, help="If true, the normalising constants are standardized so that the average on gag+pol equals 1. This way the normalised branch lengths are interpretable as typical distances on gag+pol")
  arg_parser$add_argument("-nc", "--normalisationConstants", action="store", help="Either a CSV file listing the file name for each tree (column 1) and the normalisation constant (column 2) or a single numerical normalisation constant to be applied to each window.")
  
  # Blacklisting
  
  arg_parser$add_argument("-db", "--duplicateBlacklist", action="store", help="Perform duplicate blacklisting for likely contaminants; argument here is input file produced by Python output.")
  arg_parser$add_argument("-pb", "--parsimonyBlacklist", action="store_true", help="Perform parsimony-based blacklisting for likely contaminants.")
  arg_parser$add_argument("-rwt", "--rawBlacklistThreshold", action="store", default=0, help="Raw threshold for blacklisting; subgraphs or exact duplicate tips with read counts less than this will be blacklisted, regardless of the count of any other subgraphs from the same host or identical tips from another host. Default 0; one or both of this and -rtt must be specified and >0 for -db or -pb to do anything.")
  arg_parser$add_argument("-rtt", "--ratioBlacklistThreshold", action="store", default=0, help="Ratio threshold for blacklisting; subgraphs or exact duplicate tips will be blacklisted if the ratio of their tip count to that of another subgraph from the same host or tip from another host is less than this. Default 0; one or both of this and -rwt must be specified and >0 for -db or -pb to do anything.")
  arg_parser$add_argument("-ub", "--dualBlacklist", action="store_true", default=F, help="Blacklist all reads from the minor subgraphs for all patients established as dual by parsimony blacklisting.")
  
  # Parsimony reconstruction
  
  arg_parser$add_argument("-ff", "--useff", action="store_true", default=FALSE, help="Use ff to store parsimony reconstruction matrices. Use if you run out of memory.")
  arg_parser$add_argument("-s", "--splitsRule", action="store", default="r", help="The rules by which the sets of patients are split into groups in order to ensure that all groups can be members of connected subgraphs without causing conflicts. Currently available: s=Sankoff with optional within-host diversity penalty (slow, rigorous), r=Romero-Severson (quick, less rigorous with >2 patients), f=Sankoff with continuation costs (experimental).")
  arg_parser$add_argument("-k", "--sankoffK", action="store", default = 0, help="The k parameter in the cost matrix for Sankoff reconstruction (see documentation)")
  arg_parser$add_argument("-p", "--sankoffP", action="store", default = 0, help="For the Sankoff recontruction with within-host diversity penalty, this is the (normalised) branch length threshold at which lineages transfer back to the unsampled state rather than stay in their existing host. For the reconstruction with continuation costs, this is the (normalised) branch length threshold at which an internal node is sufficiently far from all its neighbours to be reconstructed as unsampled. See documentation for both.")
  
  arg_parser$add_argument("-P", "--pruneBlacklist", action="store_true", help="If present, all blacklisted and references tips (except the outgroup) are pruned away before starting parsimony-based reconstruction")
  arg_parser$add_argument("-rcm", "--readCountsMatterOnZeroBranches", default = FALSE, action="store_true", help="If present, read counts will be taken into account in parsimony reconstructions at the parents of zero-length branches. Not applicable for the Romero-Severson-like reconstruction method.")
  
  arg_parser$add_argument("-tp", "--outputPDFTree", action="store_true", help="Output annotated trees in PDF format, via ggtree.")
  arg_parser$add_argument("-tn", "--outputNexusTree", action="store_true", help="Output annotated trees in Nexus format.")
  
  # Summary statistics
  
  arg_parser$add_argument("-ssg", "--makeSummaryStatisticsGraph", action="store_true", help="Produce summary statistic graphs.")
  arg_parser$add_argument("-ssf", "--makeSummaryStatisticsFile", action="store_true", help="Write summary statistic output to file.")
  arg_parser$add_argument("-R", "--recombinationFiles", action="store", help="An optional file path and initial string identifying all recombination data files.")
  
  # Classification
  
  arg_parser$add_argument("-cd", "--allClassifications", action="store_true", default=T, help="If present, the per-window host relationships will be writted to a separate CSV file for each window.")
  arg_parser$add_argument("-ct", "--collapsedTree", action="store_true", default=T, help="If present, the collapsed tree (in which all adjacent nodes with the same assignment are collapsed to one) is output as a CSV file or files.")
  
  # Classification summary
  
  arg_parser$add_argument("-swt", "--windowThreshold", action="store", default=0, type="integer", help="Relationships between two patients will only appear in output if they are within the distance threshold and ajacent to each other least this proportion of windows (default 0).")
  arg_parser$add_argument("-sdt", "--distanceThreshold", action="store", default=-1, help="Maximum distance threshold on a window for a relationship to be reconstructed between two patients on that window.")
  arg_parser$add_argument("-amt", "--allowMultiTrans", action="store_true", default=FALSE, help="If absent, directionality is only inferred between pairs of patients where a single clade from one patient is nested in one from the other; this is more conservative")
  
  
  args                  <- arg_parser$parse_args()
  
  verbose               <- args$verbose
  
  tree.input            <- args$tree
  blacklist.input       <- args$blacklist
  
  output.dir            <- args$outputDir
  if(is.null(output.dir)){
    output.dir          <- getwd()
  }
  output.string         <- args$outputString
  
  outgroup.name         <- args$outgroupName
  use.m.thresh          <- !is.null(args$multifurcationThreshold)
  tree.fe               <- args$treeFileExtension
  csv.fe                <- args$csvFileExtension
  
  pdf.hm                <- as.numeric(args$pdfrelheight)
  pdf.w                 <- as.numeric(args$pdfwidth)
  
  if(use.m.thresh){
    if(args$multifurcationThreshold=="g"){
      m.thresh          <- NA
    } else if(!is.na(as.numeric(args$multifurcationThreshold))){
      m.thresh          <- as.numeric(args$multifurcationThreshold)
    } else {
      stop("Unknown argument for --multifurcationThreshold specified\n")
    }
  } else {
    m.thresh            <- -1       
  }
  
  norm.ref.file.name    <- args$normRefFileName
  norm.column.name      <- args$normVar
  norm.standardise      <- args$normStandardise
  norm.constants.input  <- args$normalisationConstants

  
  tip.regex             <- args$tipRegex
  file.name.regex       <- args$fileNameRegex
  
  if(!is.null(norm.ref.file.name) & !is.null(norm.constants.input)){
    warning("Normalisation reference file name and predetermined normalisation constants both specified. Only the latter will be used.")
  }
  
  do.dup.blacklisting   <- !is.null(args$duplicateBlacklist)
  dup.input.file.name   <- args$duplicateBlacklist
  do.par.blacklisting   <- args$parsimonyBlacklist
  do.dual.blacklisting  <- args$dualBlacklist
  
  if(do.dual.blacklisting & !do.par.blacklisting){
    warning("Dual blacklisting requires parsimony blacklisting. Turning dual blacklisting off.")
    do.dual.blacklisting <- F
  }
  
  bl.raw.threshold      <- as.numeric(args$rawBlacklistThreshold)
  bl.ratio.threshold    <- as.numeric(args$ratioBlacklistThreshold)
  
  useff                 <- args$useff
  if(useff){
    suppressMessages(require(ff, quietly=TRUE, warn.conflicts=FALSE))
  }
  
  reconstruction.mode   <- args$splitsRule
  
  if(!(reconstruction.mode %in% c("r", "s", "f"))){
    stop(paste("Unknown split classifier: ", mode, "\n", sep=""))
  }
  
  sankoff.k             <- as.numeric(args$sankoffK)
  
  if(sankoff.k == 0 & reconstruction.mode == "f"){
    stop("k must be >0 for continuation cost reconstruction")
  }
  
  sankoff.p             <- as.numeric(args$sankoffP)
  read.counts.matter    <- args$readCountsMatterOnZeroBranches
  prune.blacklist       <- args$pruneBlacklist
  
  output.pdf            <- args$outputPDFTree
  output.nexus          <- args$outputNexusTree
  
  output.ssg            <- args$makeSummaryStatisticsGraph
  output.ssf            <- args$makeSummaryStatisticsFile
  
  do.summary.statistics <- output.ssg | output.ssf
  
  recomb.input          <- args$recobinationFiles
  do.recomb             <- !is.null(recomb.input)
  
  do.collapsed          <- args$collapsedTree
  do.class.detail       <- args$allClassifications
  
  win.threshold         <- args$windowThreshold 
  dist.threshold        <- args$distanceThreshold
  allow.mt              <- args$allowMultiTrans
  
  script.dir            <- args$scriptdir
  
  source(file.path(script.dir, "GeneralFunctions.R"))
  source(file.path(script.dir, "NormalisationFunctions.R"))
  source(file.path(script.dir, "TreeUtilityFunctions.R"))
  source(file.path(script.dir, "TreeUtilityFunctions2.R"))
  source(file.path(script.dir, "BlacklistFunctions.R"))
  source(file.path(script.dir, "ParsimonyReconstructionMethods2.R"))
  source(file.path(script.dir, "CollapsedTreeMethods.R"))
  source(file.path(script.dir, "WriteAnnotatedTrees.R"))
  source(file.path(script.dir, "SummariseTrees_funcs.R"))
  source(file.path(script.dir, "SummaryStatsFunctions.R"))
  source(file.path(script.dir, "PlottingFunctions.R"))
  
} else {
  
  setwd("/Users/mdhall/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20161013/Phyloscanner2Test/")
  script.dir <- ("/Users/mdhall/phylotypes/tools")
  
  tree.input          <- "RAxML_bestTree.InWindow_"
  blacklist.input     <- "AmpliconBlacklist_InWindow_"
  reconstruction.mode <- "r"
  
  root.name <- "C.BW.00.00BW07621.AF443088"
  
  do.par.blacklisting <- T
  bl.raw.threshold <- 3
  bl.ratio.threshold <- 0.005
  
  dup.input.file.name <- "DuplicateReadCountsProcessed_InWindow_"
  
  output.string <- "p2test"
  sankoff.k <- 20
  
  output.pdf <- T
  output.nexus <- T
  output.dir <- getwd()
  
  setwd("/Users/mdhall/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/PossibleTransmissionClusters/")
  script.dir <- ("/Users/mdhall/phylotypes/tools")
  
  tree.input          <- "RAxML_bestTree.InWindow_"
  blacklist.input     <- NULL
  reconstruction.mode <- "r"
  
  root.name <- "C.BW.00.00BW07621.AF443088"
  
  do.par.blacklisting <- F
  bl.raw.threshold <- 3
  bl.ratio.threshold <- 0.005
  
  dup.input.file.name <- "DuplicateReadCountsProcessed_InWindow_"
  
  output.string <- "p2test"
  sankoff.k <- 20
  
  output.pdf <- T
  output.nexus <- T
  output.dir <- getwd()
  
}

# todo If there's only one tree, skip summary statistics and transmission summary
# todo All functions that get passed tree info should check that the lists have what they need. If they have file names but not the contents, they should load the contents in.

# 1. Make the big list

all.tree.info <- list()

if(file.exists(tree.input)){
  
  tree.info                   <- list()
  tree.info$tree.file.name    <- tree.input
  tree.info$output.string     <- output.string
  tree.info$recomb.string     <- recomb.input
  #todo check blacklist exists
  
  tree.info$blacklist.input   <- blacklist.input
  
} else {
  
  tree.file.names <- list.files.mod(dirname(tree.input), pattern=paste(basename(tree.input),'.*',tree.fe,'$',sep=''), full.names=TRUE)
  
  if(length(tree.file.names)==0){
    stop("No tree files found.")
  }
  
  if(!is.null(blacklist.input)){
    blacklist.file.names <- list.files.mod(dirname(blacklist.input), pattern=paste(basename(blacklist.input),'.*\\.',csv.fe,'$',sep=''), full.names=TRUE)
  }
  
  if(do.recomb){
    recomb.file.names <- list.files.mod(dirname(recomb.input), pattern=paste(basename(recomb.input),'.*\\.',csv.fe,'$',sep=''), full.names=TRUE)
  }
  
  suffixes <- sapply(tree.file.names, get.suffix, tree.input, tree.fe)
  
  for(suffix.no in 1:length(suffixes)){
    suffix                        <- suffixes[suffix.no] 
    
    tree.info                     <- list()
    tree.info$suffix              <- suffix
    tree.info$tree.file.name      <- tree.file.names[suffix.no]
    tree.info$output.string       <- paste0(output.string, "_", suffix)
    
    expected.blacklist.file.name  <- paste0(blacklist.input, suffix , ".", csv.fe)
    
    if(file.exists(expected.blacklist.file.name)){
      tree.info$prexisting.blacklist.file.name <- expected.blacklist.file.name
    }
    
    if(do.recomb){
      expected.recomb.file.name  <- paste0(recomb.input, suffix , ".", csv.fe)
      
      if(file.exists(expected.recomb.file.name)){
        tree.info$recomb.file.name <- expected.blacklist.file.name
      } else {
        warning("Not all expected recombination files exist; will not report recombination metric")
        do.recomb <- F
      }
    }
    
    all.tree.info[[suffix]] <- tree.info
  }
}


# 2. Read the trees

all.tree.info <- sapply(all.tree.info, function(tree.info) {
  if(verbose){
    cat("Reading tree file",tree.info$tree.file.name,'\n')
  }
  
  tree                           <- read.tree(tree.info$tree.file.name)
  tree.info$tree                 <- tree
  tree.info
}, simplify = F, USE.NAMES = T)


# 3. Find some things out - are there read counts at the tips? Are there window coordinates in suffixes?

read.counts.check <- sapply(all.tree.info, function(tree.info){
  tip.labels   <- tree.info$tree$tip.label
  read.counts  <- sapply(tip.labels, read.count.from.label, tip.regex)
  return(all(is.na(read.counts)))
})

if(any(read.counts.check) & !all(read.counts.check)){
  warning("Read counts are not present in some trees; ignoring them throughout.\n")
}
no.read.counts <- all(read.counts.check)

all.tree.info <- sapply(all.tree.info, function(tree.info) {
  tryCatch({
    coords                  <- get.window.coords(tree.info$suffix, file.name.regex)
    tree.info$window.coords <- coords
    tree.info$xcoord        <- (coords$end + coords$start)/2
    tree.info
  }, error = function(e){
    warning("Cannot obtain genome window coordinates from file names. Summary stats, if requested, will be plotted in file order.")
    tree.info$xcoord        <- which(names(all.tree.info) == suffix)
    tree.info
  })
}, simplify = F, USE.NAMES = T)


# 4. Root tree, collapse multifurcations, get host for each tip

all.tree.info <- sapply(all.tree.info, function(tree.info){
  if(verbose) cat("Processing tree for suffix ",tree.info$suffix,".\n", sep="")
  
  tree <- tree.info$tree
  if(is.na(m.thresh)){
    m.thresh                          <- min(tree$edge.length*1.0001) 
  }
  
  new.tree <- process.tree(tree, outgroup.name, m.thresh)
  
  tree.info$tree <- new.tree
  tree.info$original.tip.labels       <- new.tree$tip.label
  
  hosts.for.tips                      <- sapply(tree$tip.label, function(x) host.from.label(x, tip.regex))
  tree.info$hosts.for.tips            <- hosts.for.tips
  
  tree.info
  
}, simplify = F, USE.NAMES = T)




# 5. Read the blacklists

all.tree.info <- sapply(all.tree.info, function(tree.info) {
  if(!is.null(tree.info$prexisting.blacklist.file.name)){
    if(file.exists(tree.info$prexisting.blacklist.file.name)){
      if (verbose) cat("Reading blacklist file ",tree.info$prexisting.blacklist.file.name,'\n',sep="")
      blacklisted.tips                    <- read.table(tree.info$prexisting.blacklist.file.name, sep=",", header=F, stringsAsFactors = F, col.names="read")
      blacklist                           <- vector()
      if(nrow(blacklisted.tips)>0){
        blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree.info$tree))
      }
      tree.info$hosts.for.tips[blacklist] <- NA 
      
      tree.info$blacklist                 <- blacklist
    } else {
      cat(paste("WARNING: File ",i$blacklist.input," does not exist; skipping.\n",sep=""))
    }
  } 
  
  tree.info
}, simplify = F, USE.NAMES = T)


# 6. Rename tips from the prexisting blacklist

all.tree.info <- sapply(all.tree.info, function(tree.info) {
  tree <- tree.info$tree
  
  if(is.null(tree.info$tree)){
    stop("No tree for suffix ",tree.info$suffix,"\n")
  }
  
  if(!is.null(tree.info$blacklist)){
    old.tip.labels                         <- tree$tip.label
    if(length(tree.info$blacklist)>0){
      new.tip.labels                       <- old.tip.labels
      new.tip.labels[tree.info$blacklist]  <- paste0(new.tip.labels[tree.info$blacklist], "_X_USER")
      tree$tip.label                       <- new.tip.labels
      tree.info$tree                       <- tree
      
    }
  }
  tree.info
}, simplify = F, USE.NAMES = T)


# 7. Get the normalisation constants

if(!is.null(norm.constants.input)){
  if(!is.na(suppressWarnings(as.numeric(norm.constants.input)))){
    nc <- as.numeric(norm.constants.input)
    all.tree.info <- sapply(all.tree.info, function(tree.info) {
      tree.info$normalisation.constant  <- nc
      tree.info
    }, simplify = F, USE.NAMES = T)
  } else if(file.exists(norm.constants.input)){
    nc.df   <- read.csv(norm.constants.input, stringsAsFactors = F, header = F)
    all.tree.info <- sapply(all.tree.info, function(tree.info) {
      rows  <- which(nc.df[,1]==basename(tree.info$tree.file.name))
      
      if(length(rows)>0){
        # Take the first appearance of the file name
        
        row <- rows[1]
        nc  <- as.numeric(nc.df[row, 2])
        if(is.na(nc)){
          warning("Tree with suffix ",tree.info$suffix," given a normalisation constant in file ",norm.constants.input," which cannot be processed as a number; this tree will not have normalised branch lengths\n")
          nc <- 1
        }
      } else {
        warning("Tree with suffix ",tree.info$suffix," has no normalisation constant in this lookup file and will not have normalised branch lengths\n")
        nc <- 1
      }
      
      tree.info$normalisation.constant  <- nc
      tree.info
    }, simplify = F, USE.NAMES = T)
    
  } else {
    warning("Cannot parse -nc argument; tree branch lengths will not be normalised.\n")
  }
} else if(!is.null(norm.ref.file.name)){
  if (verbose) cat('Loading normalising constants reference file ', norm.ref.file.name, "\n", sep="")
  
  if(grepl(paste0(csv.fe, "$"), norm.ref.file.name)){
    norm.table	<- as.data.table(read.csv(norm.ref.file.name, stringsAsFactors=FALSE))
    
    if(any(!(c('W_FROM','W_TO',norm.column.var) %in% colnames(norm.table)))){
      warning(paste0("One or more of W_FROM, W_TO and ",norm.var," columns not present in file ",norm.ref.file.name," skipping normalisation.\n"))
      all.tree.info <- sapply(all.tree.info, function(tree.info) {
        tree.info$normalisation.constant  <- 1
        tree.info
      }, simplify = F, USE.NAMES = T)
    } else {
      setnames(norm.table, norm.column.var, 'NORM_CONST')
      
      if(norm.standardise){
        #	Standardize to mean of 1 on gag+pol ( prot + first part of RT in total 1300bp )
        
        norm.table[, W_MID := (W_FROM+W_TO)/2]
        #790 - 3385
        if (verbose) cat('Standardising normalising constants to 1 on the gag+pol (prot + first part of RT in total 1300bp pol) region\n')
        
        tmp		<- subset(norm.table, W_MID>=790L & W_MID<=3385L)
        
        if(nrow(tmp)<=0){
          stop(paste0("No positions from gag+pol present in file ",norm.ref.file.name,"; unable to standardise"))
        }
        
        tmp		<- tmp[, mean(NORM_CONST)]
        if(!is.finite(tmp)){
          stop(paste0("Standardising constant is not finite"))
        }
        
        set(norm.table, NULL, 'NORM_CONST', norm.table[, NORM_CONST/tmp])
      }
      all.tree.info <- sapply(all.tree.info, function(tree.info){
        tree.info$normalisation.constant <- lookup.normalisation.for.tree(tree.info, norm.table, lookup.column = "NORM_CONST")
        tree.info
      }, simplify = F, USE.NAMES = T)
    }
  } else {
    warning(paste0("Unknown input file format for normalisation file; tree branch lengths will not be normalised.\n"))
  }
} else {
  all.tree.info <- sapply(all.tree.info, function(tree.info) {
    tree.info$normalisation.constant  <- 1
    tree.info
  }, simplify = F, USE.NAMES = T)
}

# 8. Apply the normalisation constants

all.tree.info <- sapply(all.tree.info, function(tree.info) {
  tree              <- tree.info$tree
  tree$edge.length  <- tree$edge.length/tree.info$normalisation.constant
  tree.info$tree    <- tree
  tree.info
}, simplify = F, USE.NAMES = T)


# 9. Duplicate blacklisting

if(do.dup.blacklisting){

  dup.file.names <- list.files.mod(dirname(dup.input.file.name), pattern=paste(basename(dup.input.file.name),'.*',csv.fe,'$',sep=''), full.names=TRUE)
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) {
    if(file.exists(paste0(dup.input.file.name, tree.info$suffix, ".", csv.fe))){
      tree.info$duplicate.tips <- strsplit(readLines(paste0(dup.input.file.name, tree.info$suffix, ".", csv.fe), warn=F),",")
    } else {
      warning("No duplicates file found for tree suffix ",tree.info$suffix, "; skipping duplicate blacklisting.")
    }
    tree.info
  }, simplify = F, USE.NAMES = T)
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) {
    tree <- tree.info$tree
    

    if(!is.null(tree.info$duplicate.tips)){

      duplicated                                   <- blacklist.exact.duplicates(tree.info, bl.raw.threshold, bl.ratio.threshold, tip.regex, verbose)
    
      duplicate.nos                                <- which(tree.info$original.tip.labels %in% duplicated)
      
      newly.blacklisted                            <- setdiff(duplicate.nos, tree.info$blacklist) 
      
      tree.info$hosts.for.tips[newly.blacklisted]  <- NA
      
      old.tip.labels    <- tree.info$tree$tip.label
      
      if(length(newly.blacklisted)>0){
        new.tip.labels                             <- old.tip.labels
        new.tip.labels[newly.blacklisted]          <- paste0(new.tip.labels[newly.blacklisted], "_X_DUPLICATE")
        tree$tip.label                             <- new.tip.labels
        tree.info$tree                             <- tree
      }
      
      
      tree.info$blacklist <- unique(c(tree.info$blacklist, duplicate.nos))
      tree.info$blacklist <- tree.info$blacklist[order(tree.info$blacklist)]
    }
    tree.info
  }, simplify = F, USE.NAMES = T)
}


# 10. Parsimony blacklisting

if(do.par.blacklisting){
  all.tree.info <- sapply(all.tree.info, function(tree.info){
    
    tree <- tree.info$tree
    tip.hosts <- sapply(tree$tip.label, function(x) host.from.label(x, tip.regex))
    tip.hosts[tree.info$blacklist] <- NA
    
    hosts <- unique(na.omit(tip.hosts))
    
    hosts <- hosts[order(hosts)]
    
    results <- sapply(hosts, function(x) get.splits.for.host(x, tip.hosts, tree, outgroup.name, bl.raw.threshold, bl.ratio.threshold, sankoff.k, T, no.read.counts, verbose), simplify = F, USE.NAMES = T)
    
    contaminant                                 <- unlist(lapply(results, "[[", 2))
    contaminant.nos                             <- which(tree.info$tree$tip.label %in% contaminant)
    
    newly.blacklisted                           <- setdiff(contaminant.nos, tree.info$blacklist) 
    
    tree.info$hosts.for.tips[newly.blacklisted] <- NA
    
    old.tip.labels                              <- tree.info$tree$tip.label
    
    if(length(newly.blacklisted)>0){
      new.tip.labels                            <- old.tip.labels
      new.tip.labels[newly.blacklisted]         <- paste0(new.tip.labels[newly.blacklisted], "_X_CONTAMINANT")
      tree$tip.label                            <- new.tip.labels
      tree.info$tree                            <- tree
    }
    
    tree.info$blacklist                         <- unique(c(tree.info$blacklist, contaminant.nos))
    tree.info$blacklist                         <- tree.info$blacklist[order(tree.info$blacklist)]
    
    which.are.duals                             <- which(unlist(lapply(results, "[[", 6)))
    
    if(length(which.are.duals) > 0) {
      mi.df                                     <- data.frame(host = unlist(sapply(results[which.are.duals], function (x) rep(x$id, length(x$tip.names)) )), 
                                                              tip.name = unlist(lapply(results[which.are.duals], "[[", 3)),
                                                              reads.in.subtree = unlist(lapply(results[which.are.duals], "[[", 4)),
                                                              tips.in.subtree = unlist(lapply(results[which.are.duals], "[[", 5)), 
                                                              stringsAsFactors = F
      )
      
      tree.info$duals.info <- mi.df
    } 
    
    tree.info
  }, simplify = F, USE.NAMES = T)
}

# 11. All IDs can be safely gathered and mapped now; dual blacklisting as performed by this script cannot eliminate patients entirely so we know who's in and who's out

if(verbose) cat("Gathering patient IDs...\n")

hosts <- lapply(all.tree.info, "[[" , "hosts.for.tips")
hosts <- unique(unlist(hosts))
hosts <- hosts[!is.na(hosts)]
hosts <- hosts[order(hosts)]

all.tree.info <- sapply(all.tree.info, function(tree.info){
  tips.for.hosts <- sapply(hosts, function(x){
    
    which(tree.info$hosts.for.tips == x)
    
  }, simplify = F, USE.NAMES = T)
  tree.info$tips.for.hosts <- tips.for.hosts
  tree.info
}, simplify = F, USE.NAMES = T)

# 12. Dual blacklisting

if(do.dual.blacklisting){
  hosts.that.are.duals <- lapply(all.tree.info, function(tree.info){
    tree.info$duals.info$host
  })
  hosts.that.are.duals <- unique(unlist(hosts.that.are.duals))
  hosts.that.are.duals <- hosts.that.are.duals[order(hosts.that.are.duals)]
  
  dual.results <- blacklist.duals(all.tree.info, hosts.that.are.duals, summary.file = NULL, verbose)
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) {
    tree <- tree.info$tree
    
    if(!is.null(dual.results[[tree.info$suffix]])){
      dual                                        <- dual.results[[tree.info$suffix]]
      dual.nos                                    <- which(tree.info$original.tip.labels %in% dual)
      
      newly.blacklisted                           <- setdiff(dual.nos, tree.info$blacklist) 
      
      tree.info$hosts.for.tips[newly.blacklisted] <- NA
      
      old.tip.labels                              <- tree.info$tree$tip.label
      
      if(length(newly.blacklisted)>0){
        new.tip.labels                            <- old.tip.labels
        new.tip.labels[newly.blacklisted]         <- paste0(new.tip.labels[newly.blacklisted], "_X_DUAL")
        tree$tip.label                            <- new.tip.labels
        tree.info$tree                            <- tree
      }
      
      tree.info$blacklist                         <- unique(c(tree.info$blacklist, dual.nos))
      tree.info$blacklist                         <- tree.info$blacklist[order(tree.info$blacklist)]
    }
    tree.info
  }, simplify = F, USE.NAMES = T)
}


# 13. Prune away the blacklist if so requested

if(prune.blacklist){
  all.tree.info <- sapply(all.tree.info, function(tree.info) {
    tree                <- tree.info$tree
    
    new.tree            <- process.tree(tree, outgroup.name, m.thresh = -1, tree.info$blacklist)
    
    tree.info$tree      <- new.tree
    tree.info$blacklist <- vector()
    
    tree.info
    
  }, simplify = F, USE.NAMES = T)
}


# 14. Parsimony reconstruction

all.tree.info <- sapply(all.tree.info, function(tree.info) {
  # Do the reconstruction
  
  tmp					     <- split.patients.to.subgraphs(tree.info$tree, tree.info$blacklist, reconstruction.mode, tip.regex, sankoff.k, sankoff.p, useff, read.counts.matter, hosts, verbose)
  tree					   <- tmp[['tree']]	
  
  # trees are annotated from now on
  
  tree.info$ tree  <- tree
  rs.subgraphs		 <- tmp[['rs.subgraphs']]
  
  # Convert back to unnormalised branch lengths for output
  
  tree$edge.length <- tree$edge.length * tree.info$normalisation.constant
  
  if(!is.null(rs.subgraphs)){		
    
    # Write PDF output
    
    if(output.pdf){
      if (verbose) cat("Drawing PDF tree...\n")
      tree.display <- ggtree(tree, aes(color=BRANCH_COLOURS)) +
        geom_point2(shape = 16, size=1, aes(color=INDIVIDUAL)) +
        scale_fill_hue(na.value = "black", drop=F) +
        scale_color_hue(na.value = "black", drop=F) +
        theme(legend.position="none") +
        geom_tiplab(aes(col=INDIVIDUAL)) + 
        geom_treescale(width=0.01, y=-5, offset=1.5)	  
      x.max <- ggplot_build(tree.display)$layout$panel_ranges[[1]]$x.range[2]	  
      tree.display <- tree.display + ggplot2::xlim(0, 1.1*x.max)
      tree.display		
      tmp	<- file.path(output.dir, paste0('ProcessedTree_', tree.info$output.string, '.pdf'))
      if (verbose) cat("Plot to file",tmp,"...\n")
      ggsave(tmp, device="pdf", height = pdf.hm*length(tree$tip.label), width = pdf.w, limitsize = F)
    }
    
    # Write nexus output

    if(output.nexus){
      tmp <- file.path(output.dir, paste0('ProcessedTree_', tree.info$output.string, ".", tree.fe))
      if (verbose) cat("Writing rerooted, multifurcating, annotated tree to file",tmp,"...\n")	
      write.ann.nexus(tree, file=tmp, annotations = c("INDIVIDUAL", "SPLIT"))
    }
    
    #	Save the table
    
    if(!no.read.counts){
      rs.subgraphs$reads <- sapply(rs.subgraphs$tip, function(x) as.numeric(read.count.from.label(x, tip.regex)))
    } else {
      rs.subgraphs$reads <- 1
    }
    
    tree.info$splits.table  <- rs.subgraphs
  }
  
  tree.info
}, simplify = F, USE.NAMES = T)


# 14. Summary statistics

if(do.summary.statistics){
  
  coordinates <- lapply(all.tree.info, "[[" , "window.coords")
  
  if(all(!is.null(coordinates))){
    starts <- sapply(coordinates, "[[", "start")
    ends <- sapply(coordinates, "[[", "end")
    
    ews <- min(starts)
    lwe <- max(ends)
  } else {
    coordinates <- sapply(all.tree.info, "[[" , "xcoord")
    range <- max(coordinates) - min(coordinates)
    increment <- range/length(coordinates)
    
    ews <- min(coordinates) - 0.45*increment
    lwe <- max(coordinates) + 0.45*increment
  }
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) {
    clade.results                 <- resolveTreeIntoPatientClades(tree.info$tree, hosts, tip.regex, tree.info$blacklist, no.read.counts)
    
    tree.info$clades.by.host      <- clade.results$clades.by.patients
    tree.info$clade.mrcas.by.host <- clade.results$clade.mrcas.by.patient
    
    tree.info
  }, simplify = F, USE.NAMES = T)
  
  pat.stats <- lapply(all.tree.info, function(x) calc.all.stats.in.window(x, hosts, tip.regex, verbose))
  pat.stats <- rbindlist(pat.stats)
  
  read.proportions <- lapply(all.tree.info, function(y) sapply(hosts, function(x) get.read.proportions(x, y$suffix, y$splits.table), simplify = F, USE.NAMES = T))
  
  # Get the max split count over every window and patient (the exact number of columns depends on this)
  
  max.splits <- max(sapply(read.proportions, function(x) max(sapply(x, function(y) length(y)))))
  
  read.prop.columns <- lapply(all.tree.info, function(x){
    out <- lapply(hosts, function(y){
      window.props <- read.proportions[[x$suffix]]
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
  
  pat.stats         <- cbind(pat.stats, read.prop.columns)
  
  pat.stats$tips    <- as.numeric(pat.stats$tips)
  
  if(output.ssf){
    tmp	<- file.path(paste0(output.dir, "/", output.string,"_patStatsFull.csv"))
    if (verbose) cat("Writing output to file ",tmp,"...\n",sep="")
    write.csv(pat.stats, tmp, quote = F, row.names = F)
  }
  if(output.ssg){
    tmp <- file.path(paste0(output.dir, "/", output.string,"_patStats.pdf"))
    if (verbose) cat("Plotting to file ",tmp,"...\n",sep="")
    
    # Set up the boundaries of each window's region on the x-axis
    
    xcoords <- unique(pat.stats$xcoord)
    xcoords <- xcoords[order(xcoords)]
    
    missing.window.data <- find.gaps(xcoords)
    
    xcoords <- missing.window.data$x.coordinates
    regular.gaps <- missing.window.data$regular.gaps
    rectangles.for.missing.windows <- missing.window.data$rectangles.for.missing.windows
    bar.width <- missing.window.data$width
    
    produce.pdf.graphs(tmp, pat.stats, hosts, xcoords, rectangles.for.missing.windows, bar.width, regular.gaps, verbose = verbose)
  }
}

# 15. Individual window classifications

all.tree.info <- sapply(all.tree.info, function(tree.info) {

  tree.info$classification.results <- classify(tree.info, verbose)
  
  if(do.collapsed){
    tree.info$collapsed.file.name <- file.path(output.dir, paste0("CollapsedTree_",tree.info$output.string,".","csv.fe"))
    write.csv(tree.info$classification.results$collapsed, tree.info$collapsed.file.name, quote=F, row.names = F)
  }
  if(do.class.detail){
    tree.info$classification.file.name <- file.path(output.dir, paste0("Classification_",tree.info$output.string,".","csv.fe"))
    write.csv(tree.info$classification.results$classification, tree.info$classification.file.name, quote=F, row.names = F)
  }
  
  tree.info
}, simplify = F, USE.NAMES = T)


# 16. Transmission summary

results <- summarise.classifications(all.tree.info, hosts, win.threshold*length(all.tree.info), dist.threshold, allow.mt, csv.fe, verbose)

if (verbose) cat('Writing summary to file',output.file,'\n')
write.csv(results, file=file.path(output.dir, paste0("TransmissionSummary_",output.string,".",csv.fe)), row.names=FALSE, quote=FALSE)