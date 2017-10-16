#!/usr/bin/env Rscript

list.of.packages <- c("argparse", "data.table", "ape", "ff", "phangorn", "phytools", "scales", "RColorBrewer", "gtable", "grid", "gridExtra", "kimisc", "GGally", "network", "sna")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Missing dependencies; replacing the path below appropriately, run\n[path to your phyloscanner code]/tools/package_install.R\nthen try again.\n")
  quit(save="no", status=1)
}

options("warn"=1)

suppressMessages(require(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(data.table, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ape, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phangorn, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ggtree, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phytools, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(network, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(scales, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(RColorBrewer, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gtable, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(grid, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gridExtra, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(kimisc, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(GGally, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(sna, quietly=TRUE, warn.conflicts=FALSE))

arg_parser		     <- ArgumentParser()

arg_parser$add_argument("tree", action="store", help="A string that begins the file names (including the path) of all input trees. e.g. path/to/RAxML_bestTree.InWindow_")
arg_parser$add_argument("outputString", action="store", help="A string that will be used to label all output files.")

arg_parser$add_argument("-og", "--outgroupName", action="store", help="The name of the tip in the phylogeny/phylogenies to be used as outgroup (if unspecified, trees will be assumed to be already rooted). This should be sufficiently distant to any sequence obtained from a host that it can be assumed that the MRCA of the entire tree was not a lineage present in any sampled individual.")
arg_parser$add_argument("-m", "--multifurcationThreshold", help="If specified, short branches in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. If a number, this number will be used as the threshold, with all branches strictly smaller collapsed. If 'g', it will be guessed from the branch lengths and the width of the genomic window (if appropriate). It is recommended that trees are examined by eye to check that they do appear to have multifurcations if 'g' is used.")
arg_parser$add_argument("-b", "--blacklist", action="store", help="A path and string that begins all the file names for pre-existing blacklist files.")

# General, bland options

arg_parser$add_argument("-od", "--outputDir", action="store", help="All output will be written to this directory. If absent, current working directory.")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what the script is doing.")
arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. This expects up to three capture groups, for host ID, read ID, and read count (in that order). If the latter two groups are missing then read information will not be used. The default is '^(.*)_read_([0-9]+)_count_([0-9]+)$', which matches input from the phyloscanner pipeline where the host ID is the BAM file name.")
arg_parser$add_argument("-y", "--fileNameRegex", action="store", default="^\\D*([0-9]+)_to_([0-9]+)\\D*$", help="Regular expression identifying window coordinates. Two capture groups: start and end; if the latter is missing then the first group is a single numerical identifier for the window. The default is '^(.*)_read_([0-9]+)_count_([0-9]+)$', which matches input from the phyloscanner pipeline.")
arg_parser$add_argument("-tfe", "--treeFileExtension", action="store", default="tree", help="The file extension for tree files (default tree).")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default csv).")
arg_parser$add_argument("-pw", "--pdfWidth", action="store", default=50, help="Width of tree PDF in inches.")
arg_parser$add_argument("-ph", "--pdfRelHeight", action="store", default=0.15, help="Relative height of tree PDF")
arg_parser$add_argument("-psb", "--pdfScaleBarWidth", action="store", default=0.01, help="Width of the scale bar in the PDF output (in branch length units)")
arg_parser$add_argument("-rda", "--outputRDA", action="store_true", help="Write the final R workspace image to file.")
arg_parser$add_argument("-sd", "--seed", action="store_true", help="Random number seed; used by the downsampling process, and also ties in some parsimony reconstructions can be broken randomly.")
arg_parser$add_argument("-D", "--toolsDir", action="store", help="Full path of the /tools/ directory.")
arg_parser$add_argument("-ow", "--overwrite", action="store_true", help="Overwrite existing output files with the same names.")

# Normalisation options

arg_parser$add_argument("-nr", "--normRefFileName", action="store", help="Name of a file giving a normalisation constant for every genome position. Cannot be used simultaneously with -nc. If both -nr and -nc are absent then no normalisation will be performed.")
arg_parser$add_argument("-ns", "--normStandardiseGagPol", action="store_true", default=FALSE, help="Use in conjunction with -nr only. An HIV-specific option: if true, the normalising constants are standardised so that the average on gag+pol equals 1. Otherwise they are standardised so the average on the whole genome equals 1.")
arg_parser$add_argument("-nc", "--normalisationConstants", action="store", help="Either a CSV file listing the file name for each tree (column 1) and the normalisation constant (column 2) or a single numerical normalisation constant to be applied to each window. If both -nr and -nc are absent then no normalisation will be performed.")

# Blacklisting

arg_parser$add_argument("-db", "--duplicateBlacklist", action="store", help="Perform blacklisting for likely contamination between samples in this data set, suggested by exact duplicate reads found in different samples. Use this option to specify a string that begins the file names (including their path) of the duplication data produced by phyloscanner_make_trees.py. For example, 'path/to/DuplicateReadCountsProcessed_InWindow_'. This option must be used in conjunction with the --rawBlacklistThreshold and/or --ratioBlacklistThreshold option.")
arg_parser$add_argument("-pbk", "--parsimonyBlacklistK", action="store", type="double", help="Perform parsimony-based blacklisting for likely contaminant sequences (including those originating from a sample outside the current data set). Use this option to specify the value of the within-host diversity penalty used (corresponding to the k parameter in the Sankoff option of the splitsRule argument). This option must be used in conjunction with the --rawBlacklistThreshold and/or --ratioBlacklistThreshold option.")
arg_parser$add_argument("-rwt", "--rawBlacklistThreshold", action="store", default=0, help="Used to specify a read count to be used as a raw threshold for blacklisting. --parsimonyBlacklistK and/or --duplicateBlacklist must also be used. If --parsimonyBlacklistK is used, subgraph with a read count strictly less than this threshold will be blacklisted. If --duplicateBlacklist is used, duplicate reads with a count strictly less than this threshold will be blacklisted. The default value of 0 means nothing is blacklisted.")
arg_parser$add_argument("-rtt", "--ratioBlacklistThreshold", action="store", default=0, help="Used to specify a read count ratio (between 0 and 1) to be used as a threshold for blacklisting. --parsimonyBlacklistK and/or --duplicateBlacklist must also be used. If --parsimonyBlacklistK is used, subgraphs will be blacklisted if the ratio of its read count to the total read count from the same host is strictly less than this threshold. If --duplicateBlacklist is used, duplicate reads will be blacklisted if the ratio of their count to the count of the duplicate (from another host) is strictly less than this threshold.")
arg_parser$add_argument("-ub", "--dualBlacklist", action="store_true", default=F, help="Blacklist all reads from the minor subgraphs for all hosts established as dual by parsimony blacklisting.")

# Downsampling

arg_parser$add_argument("-dsl", "--maxReadsPerHost", action="store", type="integer", help="If given, blacklist to downsample read counts (or tip counts if no read counts are identified) from each host to this number.")
arg_parser$add_argument("-dsb", "--blacklistUnderrepresented", action="store_true", help="If present and -dsl is given, blacklist hosts from trees where their total tip count does not reach the maximum.")

# Parsimony reconstruction

arg_parser$add_argument("-ff", "--useff", action="store_true", default=FALSE, help="Use ff to store parsimony reconstruction matrices. Use if you run out of memory.")
arg_parser$add_argument("splitsRule", action="store", help="The rules by which the sets of hosts are split into groups in order to ensure that all groups can be members of connected subgraphs without causing conflicts. This takes one string as an argument which is itself a variable number of arguments, comma-separated. The first dictates the algorithm: s=Sankoff with optional within-host diversity penalty (slow, rigorous, recommended), r=Romero-Severson (quick, less rigorous with >2 hosts), f=Sankoff with continuation costs (experimental). For 'r' no further arguments are expected. For 's' and 'f' the k parameter in the Sankoff reconstruction, which penalises within-host diversity, must also be given. There is an optional third argument for 's' and 'f'. For 's' this is the branch length threshold at which a lineage reconstructed as infecting a host will transition to the unsampled state. For 'f' this is the branch length at which an node is reconstructed as unsampled if all its neighbouring nodes are a greater distance away. Both defaults are 0.")
arg_parser$add_argument("-P", "--pruneBlacklist", action="store_true", help="If present, all blacklisted and references tips (except the outgroup) are pruned away before starting parsimony-based reconstruction.")
arg_parser$add_argument("-rcm", "--readCountsMatterOnZeroLengthBranches", default = FALSE, action="store_true", help="If present, read counts will be taken into account in parsimony reconstructions at the parents of zero-length branches. Not applicable for the Romero-Severson-like reconstruction method.")
arg_parser$add_argument("-tn", "--outputNexusTree", action="store_true", help="Standard output of annotated trees are in PDF format. If this option is present, output them as NEXUS instead.")

# Summary statistics

arg_parser$add_argument("-R", "--recombinationFiles", action="store", help="Include in the summary plots the recombination metric calculated by phyloscanner_make_trees.py. Use this option to specify a string that begins the file names (including their path) of the recombination data files. For example, 'path/to/RecombinantReads_InWindow_'.")

# Classification

arg_parser$add_argument("-cd", "--allClassifications", action="store_true", help="If present, the per-window host relationships will be writted to a separate CSV file for each window.")
arg_parser$add_argument("-ct", "--collapsedTrees", action="store_true", help="If present, the collapsed tree (in which all adjacent nodes with the same assignment are collapsed to one) is output as a CSV file or files.")

# Classification summary

arg_parser$add_argument("-swt", "--windowThreshold", action="store", default=0.5, type="double", help="Relationships between two hosts will only appear in output if they are within the distance threshold and ajacent to each other in more than this proportion of windows (default 0.5).")
arg_parser$add_argument("-sdt", "--distanceThreshold", action="store", default=-1, type="double", help="Maximum distance threshold on a window for a relationship to be reconstructed between two hosts on that window. If tree branchs lengths were normalised this will be applied to those normalised lengths. If absent then no such threshold will be applied.")
arg_parser$add_argument("-amt", "--allowMultiTrans", action="store_true", help="If absent, directionality is only inferred between pairs of hosts where a single clade from one host is nested in one from the other; this is more conservative")

# Classification simplification

arg_parser$add_argument("-sat", "--directionThreshold", action="store", default=0.33, type="double", help="In the simplified graph diagram, links will be shown as arrows if direction of transmission was inferred in at least this proportion of windows (default 0.33). Must be less than or equal to --windowThreshold.")
arg_parser$add_argument("-spd", "--summaryPlotDimensions", action="store", default=25, type="double", help="Width and height of the simplified graph PDF file in inches. Default is 25. If this output is too crowded, try increaing this.")
arg_parser$add_argument("-sks", "--skipSummaryGraph", action="store_true", help="If present, do not output a simplified relationship graph")

args                  <- arg_parser$parse_args()

verbose               <- args$verbose

overwrite             <- args$overwrite

tree.input            <- args$tree
blacklist.input       <- args$blacklist

output.dir            <- args$outputDir
if(is.null(output.dir)){
  output.dir          <- getwd()
}
output.string         <- args$outputString

if(!overwrite & file.exists(paste0(output.dir, "/", output.string,"_patStats.csv"))){
  stop("Previous output with this output string (",output.string,") detected. Please re-run with --overwrite if you wish to overwrite this.")
}

outgroup.name         <- args$outgroupName

if(is.null(output.string)){
  warning("No outgroup name provided. Trees are assumed to be correctly rooted.")
}

use.m.thresh          <- !is.null(args$multifurcationThreshold)
tree.fe               <- args$treeFileExtension
csv.fe                <- args$csvFileExtension

pdf.hm                <- as.numeric(args$pdfRelHeight)
pdf.w                 <- as.numeric(args$pdfWidth)
pdf.scale.bar.width   <- as.numeric(args$pdfScaleBarWidth)

seed                  <- args$seed

output.rda            <- args$outputRDA

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
norm.standardise.gp   <- args$normStandardiseGagPol
norm.constants.input  <- args$normalisationConstants

tip.regex             <- args$tipRegex
file.name.regex       <- args$fileNameRegex

if(!is.null(norm.ref.file.name) & !is.null(norm.constants.input)){
  warning("Normalisation reference file name and predetermined normalisation constants both specified. Only the latter will be used.")
}

do.dup.blacklisting   <- !is.null(args$duplicateBlacklist)
dup.input.file.name   <- args$duplicateBlacklist
do.par.blacklisting   <- !is.null(args$parsimonyBlacklistK)
par.blacklisting.k    <- args$parsimonyBlacklistK
do.dual.blacklisting  <- args$dualBlacklist

if(do.dual.blacklisting & !do.par.blacklisting){
  warning("Dual blacklisting requires parsimony blacklisting. Turning dual blacklisting off.")
  do.dual.blacklisting <- F
}

bl.raw.threshold      <- as.numeric(args$rawBlacklistThreshold)
bl.ratio.threshold    <- as.numeric(args$ratioBlacklistThreshold)

if(do.par.blacklisting & bl.raw.threshold == 0 & bl.ratio.threshold == 0){
  stop("Parsimony blacklisting requested but no thresholds specified with -rwt or -rtt")
}


if(do.dup.blacklisting & bl.raw.threshold == 0 & bl.ratio.threshold == 0){
  stop("Duplicate blacklisting requested but no thresholds specified with -rwt or -rtt")
}

useff                 <- args$useff
if(useff){
  suppressMessages(require(ff, quietly=TRUE, warn.conflicts=FALSE))
}

reconst.mode.arg      <- args$splitsRule

reconst.mode.arg      <- unlist(strsplit(reconst.mode.arg, ","))

if(!(reconst.mode.arg[1] %in% c("r", "s", "f"))){
  stop(paste("Unknown split classifier: ", reconst.mode.arg[1], "\n", sep=""))
}

reconstruction.mode   <- reconst.mode.arg[1]

if(reconstruction.mode == "r" & length(reconst.mode.arg)>1){
  warning("Romero-Severson reconstuction takes no additional arguments; ignoring everything after 'r'")
}

if(reconstruction.mode!="r" & length(reconst.mode.arg)==1){
  stop("At least one additional argument is required for Sankoff reconstruction")
}

if(reconstruction.mode!="r"){
  sankoff.k             <- as.numeric(reconst.mode.arg[2])
  if(sankoff.k == 0 & reconstruction.mode=="f"){
    stop("k=0 for continuation costs parsimony is not supported (for simple parsimony use 's,0')")
  }
  
  if(length(reconst.mode.arg) > 2){
    sankoff.p           <- as.numeric(reconst.mode.arg[3])
  } else {
    sankoff.p           <- 0
  }
  
  if(is.na(sankoff.k) | is.na(sankoff.p)){
    stop("Expected numerical arguments for reconstruction algorithm parameters.")
  }
} else {
  sankoff.k <- NA
  sankoff.p <- NA
}

downsample            <- !is.null(args$maxReadsPerHost)
downsampling.limit    <- args$maxReadsPerHost
blacklist.ur          <- args$blacklistUnderrepresented

read.counts.matter    <- args$readCountsMatterOnZeroLengthBranches
prune.blacklist       <- args$pruneBlacklist

output.nexus          <- args$outputNexusTree
output.pdf            <- !output.nexus

if(!("package:ggtree" %in% search())){
  warning("ggtree is not installed; annotated trees will be output in NEXUS format")
  output.nexus        <- T
  output.pdf          <- F
}

recomb.input          <- args$recombinationFiles
do.recomb             <- !is.null(recomb.input)

do.collapsed          <- args$collapsedTrees
do.class.detail       <- args$allClassifications

win.threshold         <- args$windowThreshold 
dist.threshold        <- args$distanceThreshold
arrow.threshold       <- args$directionThreshold

if(arrow.threshold >= win.threshold){
  stop("Direction threshold cannot be larger than window threshold")
}

do.simplified.graph   <- !args$skipSummaryGraph
simp.plot.dim         <- args$summaryPlotDimensions


if(dist.threshold == -1){
  dist.threshold      <- Inf
}
allow.mt              <- args$allowMultiTrans

if(!is.null(args$toolsDir)){
  tools.dir           <- args$toolsDir
} else {
  tools.dir           <- paste0(dirname(thisfile()),"/tools")
  if(!dir.exists(tools.dir)){
    stop("Cannot detect the location of the /phyloscanner/tools directory. Please specify it at the command line with -D.")
  }
}

source(file.path(tools.dir, "general_functions.R"))
source(file.path(tools.dir, "normalisation_functions.R"))
source(file.path(tools.dir, "tree_utility_functions.R"))
source(file.path(tools.dir, "blacklist_functions.R"))
source(file.path(tools.dir, "parsimony_reconstruction_methods.R"))
source(file.path(tools.dir, "collapsed_tree_methods.R"))
source(file.path(tools.dir, "write_annotated_trees.R"))
source(file.path(tools.dir, "clade_functions.R"))
source(file.path(tools.dir, "summary_stats_functions.R"))
source(file.path(tools.dir, "plotting_functions.R"))
source(file.path(tools.dir, "downsampling_functions.R"))

# todo All functions that get passed tree info should check that the lists have what they need. If they have file names but not the contents, they should load the contents in.

# 1. Make the big list

all.tree.info <- list()

if(file.exists(tree.input)){
  
  tree.info                   <- list()
  tree.info$tree.file.name    <- tree.input
  tree.info$output.string     <- output.string
  tree.info$suffix            <- "only.tree"
  
  if(do.recomb){
    warning("Only one window; recombination metric files will be ignored")
  }
  
  tree.info$prexisting.blacklist.file.name   <- blacklist.input
  
  if(!do.class.detail){
    do.class.detail       <- T
  }
  
  all.tree.info[["only.tree"]] <- tree.info
  
  single.input <- T
  
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
  
  if(any(suffixes=="")){
    stop("Some trees have identifying suffixes of length zero. The tree file argument should be a path and a string that begins every tree file, but there must be something after that string (excluding the file extension).")
  }
  
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
        tree.info$recombination.file.name <- expected.recomb.file.name
      } else {
        warning("Not all expected recombination files exist; will not report recombination metric")
        do.recomb <- F
      }
    }
    
    all.tree.info[[suffix]] <- tree.info
  }
  
  single.input <- F
}

if(length(all.tree.info)==1 & !single.input){
  warning("Only a single input tree file detected, summary statistics will not be plotted and transmission summary will be skipped.")
}

# single.file is TRUE if there is just one input file, whereas single.input is if the user specified just one file (they might have 
# specified a tree prefix which matches just one file)

single.file <- single.input | length(all.tree.info)==1

# 2. Read the trees

all.tree.info <- sapply(all.tree.info, function(tree.info) {
  if(verbose){
    cat("Reading tree file",tree.info$tree.file.name,'\n')
  }
  
  first.line                     <- readLines(tree.info$tree.file.name, n=1)
  
  if(first.line == "#NEXUS"){
    tree                           <- read.nexus(tree.info$tree.file.name)
  } else {
    tree                           <- read.tree(tree.info$tree.file.name)
  }
  
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

readable.coords <- T

all.tree.info <- sapply(all.tree.info, function(tree.info) {
  tryCatch({
    if(single.file){
      coords                <- get.window.coords(tree.info$tree.file.name, file.name.regex)
    } else {
      coords                <- get.window.coords(tree.info$suffix, file.name.regex)
    }
    tree.info$window.coords <- coords
    tree.info$xcoord        <- (coords$end + coords$start)/2
    tree.info
  }, error = function(e){
    readable.coords <<- F
    tree.info$xcoord        <- which(names(all.tree.info) == tree.info$suffix)
    tree.info
  })
}, simplify = F, USE.NAMES = T)

if(!readable.coords & !single.file){
  warning("Cannot obtain genome window coordinates from file names. Summary statistics will be plotted in alphabetical order of file name.")
}

# 4. Root tree, collapse multifurcations, get host for each tip

if(use.m.thresh & is.na(m.thresh)){
  warning("Attempting to guess a branch length threshold for multifurcations from the tree. Please visually examine the tree or trees for multifurcations before using the results of this analysis.")
}

all.tree.info <- sapply(all.tree.info, function(tree.info){
  if(verbose){
    if(!single.file){
      cat("Processing tree for suffix ",tree.info$suffix,"...\n", sep="")
    } else {
      cat("Processing tree...\n", sep="")
    }
  }
  
  tree <- tree.info$tree
  
  if(is.na(m.thresh)){
    minimum.bl                        <- min(tree$edge.length)
    if(readable.coords){
      window.width = tree.info$window.coords$end - tree.info$window.coords$start + 1
      one.snp <- 1/window.width
      
      if(minimum.bl > 0.25*one.snp){
        if(verbose){
          if(!single.file){
            cat("In window suffix ",tree.info$suffix," the minimum branch length is ",minimum.bl,", which is equivalent to ",minimum.bl/one.snp," SNPs. Assuming this tree has no multifurcations.\n", sep="")
          } else {
            cat("The minimum branch length is ",minimum.bl,", which is equivalent to ",minimum.bl/one.snp," SNPs. Assuming this tree has no multifurcations.\n", sep="")
          }
        }
        m.thresh                      <- -1
      } else {
        if(verbose) {
          if(!single.file){
            cat("In window suffix ",tree.info$suffix," the minimum branch length is ",minimum.bl,", which is equivalent to ",minimum.bl/one.snp," SNPs. Using this branch length as a multifurcation threshold.\n", sep="")
          } else {
            cat("The minimum branch length is ",minimum.bl,", which is equivalent to ",minimum.bl/one.snp," SNPs. Using this branch length as a multifurcation threshold.\n", sep="")
          }
        }
        if(minimum.bl==0){
          m.thresh                    <- 1E-9
        } else {
          m.thresh                    <- minimum.bl*1.0001
        }
      }

    } else {
      warning("Attempting to guess a branch length threshold for multifurcations from the tree. Please ensure that the tree has multifurcations before using the results of this analysis.")
      if(verbose){
        if(!single.file){
          cat("In window suffix ",tree.info$suffix," the minimum branch length is ",minimum.bl,". Using this as a multifurcation threshold.\n", sep="")
        } else {
          cat("The minimum branch length is ",minimum.bl,". Using this as a multifurcation threshold.\n", sep="")
        }
      } 
      if(minimum.bl==0){
        m.thresh                    <- 1E-9
      } else {
        m.thresh                    <- minimum.bl*1.0001
      }
    }
  }
  
  new.tree <- process.tree(tree, outgroup.name, m.thresh)
  
  tree.info$tree <- new.tree
  tree.info$original.tip.labels       <- new.tree$tip.label
  
  hosts.for.tips                      <- sapply(tree$tip.label, function(x) host.from.label(x, tip.regex))
  tree.info$hosts.for.tips            <- hosts.for.tips
  
  tree.info
  
}, simplify = F, USE.NAMES = T)

# sanity check

all.tree.info <- sapply(all.tree.info, function(tree.info){
  if(all(is.na(tree.info$hosts.for.tips))){
    warning("For tree suffix ",tree.info$suffix," no non-blacklisted tips remain; this window will be removed from the analysis.")
    NULL 
  } else {
    tree.info
  }
}, simplify = F, USE.NAMES = T)

all.tree.info[sapply(all.tree.info, is.null)] <- NULL

if(length(all.tree.info)==0){
  stop("Cannot find any hosts on any tree that match this regex. Please check that it is correct.")
}


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
      if(any(is.na(blacklist))){
        warning("Some tips listed in blacklist file ",tree.info$prexisting.blacklist.file.name," are not tips of tree ",tree.info$tree.file.name, sep="")
      }
      blacklist <- blacklist[!is.na(blacklist)]
      
      tree.info$hosts.for.tips[blacklist] <- NA 
      
      if(verbose & length(blacklist)>0) {
        if(!single.file){
          cat(length(blacklist), " tips pre-blacklisted for tree suffix ",tree.info$suffix, ".\n", sep="")
        } else {
          cat(length(blacklist), " tips pre-blacklisted.\n", sep="")
        }
      }
      
      tree.info$blacklist                 <- blacklist
    } else {
      cat(paste("WARNING: File ",tree.info$blacklist.input," does not exist; skipping.\n",sep=""))
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
  if(readable.coords){
    if (verbose) cat('Loading normalising constants reference file ', norm.ref.file.name, "\n", sep="")
    
    if(grepl(paste0(csv.fe, "$"), norm.ref.file.name)){
      norm.table	<- as.data.table(read.csv(norm.ref.file.name, stringsAsFactors=FALSE))
      
      if(ncol(norm.table)!=2){
        stop(paste0(norm.ref.file.name," is not formatted as expected for a normalisation lookup file; expecting two columns.\n"))
      } else {
        setnames(norm.table, 1, 'POSITION')
        setnames(norm.table, 2, 'NORM_CONST')
        
        if(norm.standardise.gp){
          #	Standardize to mean of 1 on gag+pol ( prot + first part of RT in total 1300bp )
          
          #790 - 3385
          if (verbose) cat('Standardising normalising constants to 1 on the gag+pol (prot + first part of RT in total 1300bp pol) region\n')
          
          tmp		<- subset(norm.table, POSITION>=790L & POSITION<=3385L)
          
          if(nrow(tmp)<=0){
            stop(paste0("No positions from gag+pol present in file ",norm.ref.file.name,"; unable to standardise"))
          }
          
          tmp		<- tmp[, mean(NORM_CONST)]
          
          if(!is.finite(tmp)){
            stop(paste0("Standardising constant is not finite"))
          }
          
          set(norm.table, NULL, 'NORM_CONST', norm.table[, NORM_CONST/tmp])
        } else {
          if (verbose) cat('Standardising normalising constants to 1 on the whole genome\n')
          
          tmp		<- norm.table[, mean(NORM_CONST)]
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
      all.tree.info <- sapply(all.tree.info, function(tree.info) {
        tree.info$normalisation.constant  <- 1
        tree.info
      }, simplify = F, USE.NAMES = T)
    }
  } else {
    warning(paste0("Cannot normalise branch lengths from file without window cooardinates in file suffixes; tree branch lengths will not be normalised.\n"))
    all.tree.info <- sapply(all.tree.info, function(tree.info) {
      tree.info$normalisation.constant  <- 1
      tree.info
    }, simplify = F, USE.NAMES = T)
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
  
  if(!single.input){
    all.tree.info <- sapply(all.tree.info, function(tree.info) {
      if(file.exists(paste0(dup.input.file.name, tree.info$suffix, ".", csv.fe))){
        tree.info$duplicate.tips <- strsplit(readLines(paste0(dup.input.file.name, tree.info$suffix, ".", csv.fe), warn=F),",")
      } else {
        warning("No duplicates file found for tree suffix ",tree.info$suffix, "; skipping duplicate blacklisting.")
      }
      tree.info
    }, simplify = F, USE.NAMES = T)
  } else {
    all.tree.info[[1]]$duplicate.tips <- strsplit(readLines(dup.input.file.name),",")
  }
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) {
    tree <- tree.info$tree
    
    if(!is.null(tree.info$duplicate.tips)){
      
      duplicated                                   <- blacklist.exact.duplicates(tree.info, bl.raw.threshold, bl.ratio.threshold, tip.regex, verbose)
      
      duplicate.nos                                <- which(tree.info$original.tip.labels %in% duplicated)
      
      newly.blacklisted                            <- setdiff(duplicate.nos, tree.info$blacklist) 
      
      if(verbose & length(newly.blacklisted > 0)) cat(length(newly.blacklisted), " tips blacklisted as duplicates for tree suffix ",tree.info$suffix, "\n", sep="")
      
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
    
    results <- sapply(hosts, function(x) get.splits.for.host(x, tip.hosts, tree, outgroup.name, bl.raw.threshold, bl.ratio.threshold, "s", par.blacklisting.k, 0, T, no.read.counts, verbose), simplify = F, USE.NAMES = T)

    contaminant                                 <- unlist(lapply(results, "[[", 2))
    contaminant.nos                             <- which(tree.info$tree$tip.label %in% contaminant)
    
    newly.blacklisted                           <- setdiff(contaminant.nos, tree.info$blacklist) 
    
    if(verbose & length(newly.blacklisted)>0) cat(length(newly.blacklisted), " tips blacklisted as probable contaminants by parsimony reconstruction for tree suffix ",tree.info$suffix, "\n", sep="")
    
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
    mi.count                                    <- unlist(lapply(results, "[[", 7))
    
    multiplicity.table                          <- data.frame(host = hosts, count = mi.count, stringsAsFactors = F)
    tree.info$dual.detection.splits             <- multiplicity.table
    
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

# 11. Dual blacklisting

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
      
      if(verbose & length(newly.blacklisted)>0) cat(length(newly.blacklisted), " tips blacklisted for belonging to minor subgraphs in tree suffix ",tree.info$suffix, "\n", sep="")
      
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


# 12. Downsampling

if(downsample){
  all.tree.info <- sapply(all.tree.info, function(tree.info){
    tree.info <- downsample.tree(tree.info, NULL, downsampling.limit, T, blacklist.ur, no.read.counts, seed, verbose)
    tree.info$hosts.for.tips[tree.info$blacklist] <- NA 
    
    tree.info
  }, simplify = F, USE.NAMES = T)
}

# 13. All IDs can be safely gathered and mapped now. Windows with no patients should be removed.

if(verbose) cat("Gathering host IDs...\n")

hosts <- lapply(all.tree.info, "[[" , "hosts.for.tips")
hosts <- unique(unlist(hosts))
hosts <- hosts[!is.na(hosts)]
hosts <- hosts[order(hosts)]

if(length(hosts)==1){
  warning("Only one host detected in any tree, will skip classification of topological relationships between hosts.")
}

all.tree.info <- sapply(all.tree.info, function(tree.info){
  if(all(is.na(tree.info$hosts.for.tips))){
    warning("For tree suffix ",tree.info$suffix," no non-blacklisted tips remain; this window will be removed from the analysis.")
    NULL 
  } else {
    tips.for.hosts <- sapply(hosts, function(x){
      
      which(tree.info$hosts.for.tips == x)
      
    }, simplify = F, USE.NAMES = T)
    tree.info$tips.for.hosts <- tips.for.hosts
    tree.info
  }
}, simplify = F, USE.NAMES = T)

all.tree.info[sapply(all.tree.info, is.null)] <- NULL

if(length(all.tree.info)==0){
  stop("All hosts have been blacklisted from all trees; nothing to do.")
}

if(length(all.tree.info)==1 & !single.file){
  warning("Only one window with any hosts is present after blacklisting, summary statistics will not be plotted and transmission summary will be skipped.")
  single.file <- T
}

# 14. Prune away the blacklist if so requested

if(prune.blacklist){
  all.tree.info <- sapply(all.tree.info, function(tree.info) {
    tree                <- tree.info$tree
    
    new.tree            <- process.tree(tree, outgroup.name, m.thresh = -1, tree.info$blacklist)
    
    tree.info$tree      <- new.tree
    tree.info$blacklist <- vector()
    
    tree.info
    
  }, simplify = F, USE.NAMES = T)
}


# 15. Parsimony reconstruction

all.tree.info <- sapply(all.tree.info, function(tree.info) {
  # Do the reconstruction
  
  if(!single.file){
    if(verbose) cat("Reconstructing internal node hosts on tree suffix ",tree.info$suffix, "\n", sep="")
  } else {
    if(verbose) cat("Reconstructing internal node hosts\n", sep="")
  }
  
  tmp					     <- split.hosts.to.subgraphs(tree.info$tree, tree.info$blacklist, reconstruction.mode, tip.regex, sankoff.k, sankoff.p, useff, read.counts.matter, hosts, verbose)
  tree					   <- tmp[['tree']]	
  
  # trees are annotated from now on
  
  rs.subgraphs		 <- tmp[['rs.subgraphs']]
  
  if(is.null(rs.subgraphs)){
    warning("No non-blacklisted tips for suffix ",tree.info$suffix,"; this window will not be analysed further.")
  }
  
  # Convert back to unnormalised branch lengths for output
  
  tree$edge.length <- tree$edge.length * tree.info$normalisation.constant
  
  if(!is.null(rs.subgraphs)){		
    
    # Write PDF output
    
    if(output.pdf){
      if (verbose) cat("Drawing PDF tree...\n")
      tree.display <- ggtree(tree, aes(color=BRANCH_COLOURS)) +
        geom_point2(aes(subset=SUBGRAPH_MRCA, color=INDIVIDUAL), shape = 23, size = 3, fill="white") +
        geom_point2(aes(color=INDIVIDUAL), shape=16, size=1) +
        scale_fill_hue(na.value = "black", drop=F) +
        scale_color_hue(na.value = "black", drop=F) +
        theme(legend.position="none") +
        geom_tiplab(aes(col=INDIVIDUAL)) + 
        geom_treescale(width=pdf.scale.bar.width, y=-5, offset=1.5)	  
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
  
  tree.info$tree <- tree
  
  tree.info
  
}, simplify = F, USE.NAMES = T)


# 16. Summary statistics

all.tree.info <- sapply(all.tree.info, function(tree.info) {
  clade.results                 <- resolveTreeIntoPatientClades(tree.info$tree, hosts, tip.regex, tree.info$blacklist, no.read.counts)
  
  tree.info$clades.by.host      <- clade.results$clades.by.patient
  tree.info$clade.mrcas.by.host <- clade.results$clade.mrcas.by.patient
  
  tree.info
}, simplify = F, USE.NAMES = T)

pat.stats <- lapply(all.tree.info, function(x) calc.all.stats.in.window(x, hosts, tip.regex, verbose))
pat.stats <- rbindlist(pat.stats)

read.proportions <- lapply(all.tree.info, function(y) sapply(hosts, function(x) get.read.proportions(x, y$suffix, y$splits.table), simplify = F, USE.NAMES = T))

# Get the max split count over every window and host (the exact number of columns depends on this)

max.splits <- max(sapply(read.proportions, function(x) max(sapply(x, function(y) length(y)))))

read.prop.columns <- lapply(all.tree.info, function(x){
  window.props <- read.proportions[[x$suffix]]
  out <- lapply(hosts, function(y){
    
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

tmp	<- file.path(paste0(output.dir, "/", output.string,"_patStats.csv"))
if (verbose) cat("Writing output to file ",tmp,"...\n",sep="")
write.csv(pat.stats, tmp, quote = F, row.names = F)

if(!single.file){
  coordinates <- lapply(all.tree.info, "[[" , "window.coords")
  
  if(readable.coords){
    coordinates <- lapply(all.tree.info, "[[" , "window.coords")
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
  
  x.limits <- c(ews, lwe)
  
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
  
  produce.pdf.graphs(tmp, pat.stats, hosts, xcoords, x.limits, rectangles.for.missing.windows, bar.width, regular.gaps, readable.coords = readable.coords, verbose = verbose)
}

# for(tree.info in all.tree.info){
#   splits.file.name <- file.path(output.dir, paste0("Splits_",tree.info$output.string,".",csv.fe))
#   write.csv(tree.info$splits.table, splits.file.name, quote=F, row.names = F)
# }


# 17. Individual window classifications

if(length(hosts)>1){
  all.tree.info <- sapply(all.tree.info, function(tree.info) {
    
    if(!single.file){
      if(verbose) cat("Classifying pairwise host relationships for tree suffix ",tree.info$suffix, ".\n", sep="")
    } else {
      if(verbose) cat("Classifying pairwise host relationships.\n", sep="")
    }
    
    tree.info$classification.results <- classify(tree.info, verbose)
    
    if(do.collapsed){
      tree.info$collapsed.file.name <- file.path(output.dir, paste0("CollapsedTree_",tree.info$output.string,".",csv.fe))
      write.csv(tree.info$classification.results$collapsed[,1:4], tree.info$collapsed.file.name, quote=F, row.names = F)
    }
    if(do.class.detail){
      tree.info$classification.file.name <- file.path(output.dir, paste0("Classification_",tree.info$output.string,".",csv.fe))
      write.csv(tree.info$classification.results$classification, tree.info$classification.file.name, quote=F, row.names = F)
    }
    
    tree.info
  }, simplify = F, USE.NAMES = T)
}


# 18. Transmission summary

if(!single.file & length(hosts)>1){
  results <- summarise.classifications(all.tree.info, win.threshold*length(all.tree.info), dist.threshold, allow.mt, verbose)
  
  if (verbose) cat('Writing summary to file', paste0(output.string,"_hostRelationshipSummary.",csv.fe),'\n')
  write.csv(results, file=file.path(output.dir, paste0(output.string,"_hostRelationshipSummary.",csv.fe)), row.names=FALSE, quote=FALSE)
  
  results$ancestry <- as.character(results$ancestry)
  
  if(do.simplified.graph){
    simplified.graph <- simplify.summary(results, arrow.threshold, length(all.tree.info), plot = T)
    simplified.graph$simp.diagram
    ggsave(file = file.path(output.dir, paste0(output.string,"_simplifiedRelationshipGraph.pdf")), width=simp.plot.dim, height=simp.plot.dim)
  }
}

# 19. Workspace image

if(output.rda){
  if (verbose) cat('Saving R workspace image to file', paste0("workspace_",output.string,".rda"),'\n')
  save(all.tree.info, file=file.path(output.dir, paste0("workspace_",output.string,".rda")))
}
