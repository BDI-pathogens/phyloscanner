#!/usr/bin/env Rscript

options("warn"=1)

suppressMessages(require(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phyloscannerR, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(network, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ggplot2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(sna, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(scales, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(readr, quietly=TRUE, warn.conflicts=FALSE))

arg_parser		     <- ArgumentParser()


arg_parser$add_argument("tree", action="store", help="A string that begins the file names (including the path) of all input trees. e.g. path/to/RAxML_bestTree.InWindow_, a directory containing all trees, or a single input tree file.")
arg_parser$add_argument("outputString", action="store", help="A string that will be used to label all output files.")

arg_parser$add_argument("-og", "--outgroupName", action="store", help="The name of the tip in the phylogeny/phylogenies to be used as outgroup (if unspecified, trees will be assumed to be already rooted). This should be sufficiently distant to any sequence obtained from a host that it can be assumed that the MRCA of the entire tree was not a lineage present in any sampled individual.")
arg_parser$add_argument("-m", "--multifurcationThreshold", help="If specified, short branches in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. If a number, this number will be used as the threshold, with all branches strictly smaller collapsed. If 'g', it will be guessed from the branch lengths and the width of the genomic window (if appropriate). It is recommended that trees are examined by eye to check that they do appear to have multifurcations if 'g' is used.")
arg_parser$add_argument("-b", "--userBlacklist", action="store", help="A path and string that begins all the file names for pre-existing blacklist files, or the path of a directory containing all those files.")
arg_parser$add_argument("-aln", "--alignment", action="store", help="A path and string that begins all the file names for alignments, or a directory containing all alignment files. Needed if ancestral state reconstruction is desired at a later point.")

# General, bland options

arg_parser$add_argument("-od", "--outputDir", action="store", help="All output will be written to this directory. If absent, current working directory.")
arg_parser$add_argument("-v", "--verbose", action="store", type="integer", default=1, help="The level of verbosity in command line output. 0=none, 1=minimal, 2=detailed")
arg_parser$add_argument("-npb", "--noProgressBars", action="store_true", default=FALSE, help="If --verbose, do not display progress bars")
arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. This expects up to three capture groups, for host ID, read ID, and read count (in that order). If the latter two groups are missing then read information will not be used. If absent, the default is '^(.*)_read_([0-9]+)_count_([0-9]+)$', which matches input from the phyloscanner pipeline where the host ID is the BAM file name.")
arg_parser$add_argument("-y", "--fileNameRegex", action="store", default="^(?:.*\\D)?([0-9]+)_to_([0-9]+).*$", help="Regular expression identifying window coordinates in tree file names. Two capture groups: start and end; if the latter is missing then the first group is a single numerical identifier for the window. If absent, input will be assumed to be from the phyloscanner pipeline.")
arg_parser$add_argument("-tfe", "--treeFileExtension", action="store", default=".tree", help="The file extension for tree files (default .tree). This includes the dot; use '' if there is no file extension.")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default=".csv", help="The file extension for table files (default .csv). This includes the dot; use '' if there is no file extension..")
arg_parser$add_argument("-afe", "--alignmentFileExtension", action="store", default=".fasta", help="The file extension for nucleotide alignment files (default .fasta). This includes the dot; use '' if there is no file extension.")
arg_parser$add_argument("-aft", "--alignmentFormat", action="store", default="fasta", help="File format for alignment files. The options are the same as ape::read.dna.")
arg_parser$add_argument("-pw", "--pdfWidth", action="store", default=50, help="Width of tree PDF in inches.")
arg_parser$add_argument("-ph", "--pdfRelHeight", action="store", default=0.15, help="Relative height of tree PDF")
arg_parser$add_argument("-psb", "--pdfScaleBarWidth", action="store", default=0.01, help="Width of the scale bar in the PDF output (in branch length units)")
arg_parser$add_argument("-rda", "--outputRDA", action="store_true", help="Write the final R workspace image to file.")
arg_parser$add_argument("-rdo", "--RDAonly", action="store_true", help="Suppress all non-RDA output. Overrides --blacklistReport, --allClassifications and --collapsedTrees.")
arg_parser$add_argument("-sd", "--seed", action="store", help="Random number seed; used by the downsampling process, and also ties in some parsimony reconstructions can be broken randomly.")
arg_parser$add_argument("-ow", "--overwrite", action="store_true", help="Overwrite existing output files with the same names.")

# Normalisation options

arg_parser$add_argument("-nr", "--normRefFileName", action="store", help="Name of a file giving a normalisation constant for every genome position. Cannot be used simultaneously with -nc. If both -nr and -nc are absent then no normalisation will be performed. If you are running this script on trees produced by phyloscanner_make_trees.py, and you want to use this option with a file produced by tools/CalculateTreeSizeInGenomeWindows.py, then phyloscanner_make_trees.py must have been run with the --pairwise-align-to option and the reference specified there must match the reference whose coordinates were used for the output of tools/CalculateTreeSizeInGenomeWindows.py. This ensures the coordinates match for correct normalisation of the trees.")
arg_parser$add_argument("-ns", "--normStandardiseGagPol", action="store_true", default=FALSE, help="Use in conjunction with -nr only. An HIV-specific option: if true, the normalising constants are standardised so that the average on gag+pol equals 1. Otherwise they are standardised so the average on the whole genome equals 1.")
arg_parser$add_argument("-nc", "--normalisationConstants", action="store", help="Either a CSV file listing the file name for each tree (column 1) and the normalisation constant (column 2) or a single numerical normalisation constant to be applied to each window. If both -nr and -nc are absent then no normalisation will be performed.")

# Blacklisting

arg_parser$add_argument("-db", "--duplicateBlacklist", action="store", help="Perform blacklisting for likely contamination between samples in this data set, suggested by exact duplicate reads found in different samples. Use this option to specify a string that begins the file names (including their path) of the duplication data produced by phyloscanner_make_trees.py (for example, 'path/to/DuplicateReadCountsProcessed_InWindow_') or a directory containing those files. This option must be used in conjunction with the --rawBlacklistThreshold and/or --ratioBlacklistThreshold option.")
arg_parser$add_argument("-pbk", "--parsimonyBlacklistK", action="store", type="double", help="Perform parsimony-based blacklisting for likely contaminant sequences (including those originating from a sample outside the current data set). Use this option to specify the value of the within-host diversity penalty used (corresponding to the k parameter in the Sankoff option of the splitsRule argument). This option must be used in conjunction with the --rawBlacklistThreshold and/or --ratioBlacklistThreshold option.")
arg_parser$add_argument("-rwt", "--rawBlacklistThreshold", action="store", default=0, help="Used to specify a read count to be used as a raw threshold for blacklisting. --parsimonyBlacklistK and/or --duplicateBlacklist must also be used. If --parsimonyBlacklistK is used, subgraphs with a read count strictly less than this threshold will be blacklisted. If --duplicateBlacklist is used, duplicate reads with a count strictly less than this threshold will be blacklisted. The default value of 0 means nothing is blacklisted.")
arg_parser$add_argument("-rtt", "--ratioBlacklistThreshold", action="store", default=0, help="Used to specify a read count ratio (between 0 and 1) to be used as a threshold for blacklisting. --parsimonyBlacklistK and/or --duplicateBlacklist must also be used. If --parsimonyBlacklistK is used, a subgraph will be blacklisted if the ratio of its read count to the total read count from the same host (in that tree) is strictly less than this threshold. If --duplicateBlacklist is used, duplicate reads will be blacklisted if the ratio of their count to the count of the duplicate (from another host) is strictly less than this threshold.")
arg_parser$add_argument("-ub", "--dualBlacklist", action="store_true", help="Blacklist all reads from the minor subgraphs for all hosts established as dual by parsimony blacklisting.")
arg_parser$add_argument("-mr", "--minReadsPerHost", action="store", type="integer", default = 1, help="Blacklist hosts from trees where they have less than this number of reads.")
arg_parser$add_argument("-mt", "--minTipsPerHost", action="store",  type="integer", default = 1, help="Blacklist hosts from trees where they have less than this number of tips.")
arg_parser$add_argument("-phcb", "--postHocCountBlacklisting", action="store_true", help="Perform minimum read and tip based blacklisting as a separate step at the end of the analysis. (A legacy option).")

# Downsampling

arg_parser$add_argument("-dsl", "--maxReadsPerHost", action="store", type="integer", help="If given, blacklist to downsample read counts (or tip counts if no read counts are identified) from each host to this number.")
arg_parser$add_argument("-dsb", "--blacklistUnderrepresented", action="store_true", help="If present and -dsl is given, blacklist hosts from trees where their total tip count does not reach the maximum.")

# Blacklisting report

arg_parser$add_argument("-blr", "--blacklistReport", action="store_true", help="If present, output a CSV file of blacklisted tips from each tree.")

# Parsimony reconstruction

arg_parser$add_argument("-ff", "--useff", action="store_true", default=FALSE, help="Use ff to store parsimony reconstruction matrices. Use if you run out of memory.")
arg_parser$add_argument("splitsRule", action="store", help="The rules by which the sets of hosts are split into groups in order to ensure that all groups can be members of connected subgraphs without causing conflicts. This takes one string as an argument which is itself a variable number of arguments, comma-separated. The first dictates the algorithm: s=Sankoff with optional within-host diversity penalty (slow, rigorous, recommended), r=Romero-Severson (quick, less rigorous with >2 hosts). For 'r' no further arguments are expected. For 's' the k parameter in the Sankoff reconstruction, which penalises within-host diversity, must also be given.")
arg_parser$add_argument("-P", "--pruneBlacklist", action="store_true", help="If present, all blacklisted and references tips (except the outgroup) are pruned away before starting parsimony-based reconstruction")
arg_parser$add_argument("-rcm", "--readCountsMatterOnZeroLengthBranches", default = FALSE, action="store_true", help="If present, read counts will be taken into account in parsimony reconstructions at the parents of zero-length branches. Not applicable for the Romero-Severson-like reconstruction method.")
arg_parser$add_argument("-tn", "--outputNexusTree", action="store_true", help="Standard output of annotated trees are in PDF format. If this option is present, output them as NEXUS instead.")

# Summary statistics

arg_parser$add_argument("-R", "--recombinationFiles", action="store", help="Include in the summary plots the recombination metric calculated by phyloscanner_make_trees.py. Use this option to specify a string that begins the file names (including their path) of the recombination data files (for example, 'path/to/RecombinantReads_InWindow_') or a directory containing those files.")
arg_parser$add_argument("-ast", "--alignmentStatistics", action = "store_true", help = "Adds columns for nucleotide diversity and cumulative minor allele frequency to the summary statistics file. Requires --alignment." )

# Classification

arg_parser$add_argument("-cd", "--allClassifications", action="store_true", help="If present, the per-window host relationships will be writted to a separate CSV file for each window.")
arg_parser$add_argument("-ct", "--collapsedTrees", action="store_true", help="If present, the collapsed tree (in which all adjacent nodes with the same assignment are collapsed to one) is output as a CSV file or files.")

# Classification summary

arg_parser$add_argument("-swt", "--windowThreshold", action="store", default=0.5, type="double", help="Relationships between two hosts will only appear in output if they are within the distance threshold and ajacent to each other in more than this proportion of windows (default 0.5).")
arg_parser$add_argument("-sdt", "--distanceThreshold", action="store", nargs="+", default=-1, type="double", help="Maximum distance threshold on a window for a relationship to be reconstructed between two hosts on that window. If tree branchs lengths were normalised this will be applied to those normalised lengths. If absent then no such threshold will be applied. If --multinomial is also specified and a second value is given, the two values will be used as the 'close threshold' and 'distant threshold' in that model. Any arguments beyond the first will otherwise be ignored.")
arg_parser$add_argument("-amt", "--allowMultiTrans", action="store_true", help="If absent, directionality is only inferred between pairs of hosts where a single clade from one host is nested in one from the other; this is more conservative")
arg_parser$add_argument("-rla", "--relaxedAncestry", action="store_true", help="If absent, directionality can be inferred so long as at least one subraph from one host is descended from one from the other, and no pair of subgraphs exist in the opposite direction. Otherwise it is required that every subgraph from one host is descended from one from the other.")
arg_parser$add_argument("-mlt", "--multinomial", action="store_true", help="Use the adjustment for missing and overlapping windows as described in Ratmann et al., Nature Communications, 2019.")

# Classification simplification

arg_parser$add_argument("-sat", "--directionThreshold", action="store", default=0.33, type="double", help="In the simplified graph diagram, links will be shown as arrows if direction of transmission was inferred in at least this proportion of windows (default 0.33). Must be less than or equal to --windowThreshold.")
arg_parser$add_argument("-spd", "--summaryPlotDimensions", action="store", default=25, type="double", help="Width and height of the simplified graph PDF file in inches. Default is 25. If this output is too crowded, try increaing this.")
arg_parser$add_argument("-sks", "--skipSummaryGraph", action="store_true", help="If present, do not output a simplified relationship graph")

args                            <- arg_parser$parse_args()

# make packages less chatty

options(readr.num_columns = 0)

# basics
verbosity                       <- args$verbose

if(!(verbosity %in% 0:2)){
  stop("Invalid verbosity option: ",verbosity)
}

no.progress.bars                <- args$noProgressBars
overwrite                       <- args$overwrite
tree.fe                         <- args$treeFileExtension
re.tree.fe                      <- gsub("\\.", "\\\\.", tree.fe)
csv.fe                          <- args$csvFileExtension
re.csv.fe                       <- gsub("\\.", "\\\\.", csv.fe)
alignment.fe                    <- args$alignmentFileExtension
re.alignment.fe                 <- gsub("\\.", "\\\\.", alignment.fe)
alignment.format                <- args$alignmentFormat

# tree input
tree.input                      <- args$tree

if(!file.exists(tree.input)){
  tree.directory                <- dirname(tree.input)

  tree.file.regex               <- paste0("^", basename(tree.input), "(.*)", re.tree.fe, "$")
  single.tree                   <- F

} else {
  if(file.info(tree.input)[['isdir']]){
    tree.directory              <- tree.input
    tree.file.regex             <- paste0("^(.*)",re.tree.fe,"$")
    single.tree                 <- F
  } else {
    single.tree                 <- T
  }
}

# user blacklist
blacklist.input                 <- args$userBlacklist

if(!is.null(blacklist.input)){
  if(!file.exists(blacklist.input)){
    user.blacklist.directory    <- dirname(blacklist.input)
    user.blacklist.file.regex   <- paste0("^", basename(blacklist.input), "(.*)", re.csv.fe,"$")
  } else {
    if(file.info(blacklist.input)[['isdir']]){
      user.blacklist.directory  <- blacklist.input
      user.blacklist.file.regex <- paste0("^(.*)",re.csv.fe,"$")
    } 
  }
} else {
  user.blacklist.directory      <- NULL
  user.blacklist.file.regex     <- NULL
}

# alignment files
alignment.input                 <- args$alignment

if(!is.null(alignment.input)){
  if(!file.exists(alignment.input)){
    alignment.directory         <- dirname(alignment.input)
    alignment.file.regex        <- paste0("^", basename(alignment.input), "(.*)", re.alignment.fe,"$")
  } else {
    if(file.info(alignment.input)[['isdir']]){
      alignment.directory       <- alignment.input
      
      # this is a heinous hack. PMT has two types of alignment and if one exists then it is the right one. If it does not, then the other one is.
      if(length(list.files(alignment.directory, pattern = "^AlignedReadsInWindow_[0-9].*")) > 0){
        # this is PMT output, we assume
        alignment.file.regex      <- paste0("^AlignedReadsInWindow_([0-9]+_to_[0-9]+)",re.alignment.fe,"$")
      } else {
        # it isn't, just take the whole folder
        alignment.file.regex      <- paste0("^(.*)",re.alignment.fe,"$")
      }
      
      
    } else {
      alignment.directory         <- dirname(alignment.input)
      alignment.file.regex        <- alignment.input
    }
  }
} else {
  alignment.directory           <- NULL
  alignment.file.regex          <- NULL
}

# output files
output.dir                      <- args$outputDir
if(is.null(output.dir)){
  output.dir                    <- getwd()
}

if(!dir.exists(output.dir)){
  stop("Specified output directory (",output.dir,") does not exist.")
}

output.string                   <- args$outputString


if(!overwrite & length(list.files(path=output.dir, pattern=paste0("^", output.string, ".*")))>0 ){
  stop("Previous output with this output string (",output.string,") detected. Please re-run with --overwrite if you wish to overwrite this.")
}

# outgroup name
outgroup.name                   <- args$outgroupName

if(is.null(outgroup.name)){
  warning("No outgroup name provided. Trees are assumed to be correctly rooted.")
}

# multifurcation threshold
use.m.thresh                    <- !is.null(args$multifurcationThreshold)

if(use.m.thresh){
  if(args$multifurcationThreshold=="g"){
    m.thresh                    <- NA
  } else if(!is.na(as.numeric(args$multifurcationThreshold))){
    m.thresh                    <- as.numeric(args$multifurcationThreshold)
  } else {
    stop("Unknown argument for --multifurcationThreshold specified\n")
  }
} else {
  m.thresh                      <- -1
}

# PDF tree dimensions
pdf.hm                        <- as.numeric(args$pdfRelHeight)
pdf.w                         <- as.numeric(args$pdfWidth)
pdf.scale.bar.width           <- as.numeric(args$pdfScaleBarWidth)

# Random number seed
seed                          <- args$seed
if(is.null(seed)){
  seed                        <- sample.int(1000000000, 1)
} else {
  seed                        <- as.numeric(seed)
}
if(verbosity!=0) cat("Random number seed is",seed,"\n")
set.seed(seed)

# Output RDA?
output.rda                    <- args$outputRDA | args$RDAonly
output.files                  <- !args$RDAonly

# Normalisation options
norm.ref.file.name            <- args$normRefFileName
norm.standardise.gp           <- args$normStandardiseGagPol
norm.constants.input          <- args$normalisationConstants
if(!is.null(norm.ref.file.name) & !is.null(norm.constants.input)){
  warning("Normalisation reference file name and predetermined normalisation constants both specified. Only the latter will be used.")
}

# Regexes for suffixes and window coordinates
tip.regex                     <- args$tipRegex
file.name.regex               <- args$fileNameRegex

# What blacklisting to do and how
do.dup.blacklisting           <- !is.null(args$duplicateBlacklist)
dup.input.file.name           <- args$duplicateBlacklist
if(!is.null(dup.input.file.name)){
  if(!file.exists(dup.input.file.name)){
    duplicate.file.directory    <- dirname(dup.input.file.name)
    duplicate.file.regex        <- paste0("^", basename(dup.input.file.name), "(.*)", re.csv.fe,"$")
  } else {
    if(file.info(dup.input.file.name)[['isdir']]){
      duplicate.file.directory  <- dup.input.file.name
      duplicate.file.regex      <- paste0("^(.*)",re.csv.fe,"$")
    } 
  }
} else {
  duplicate.file.directory    <- NULL
  duplicate.file.regex        <- NULL
}

do.par.blacklisting           <- !is.null(args$parsimonyBlacklistK)
par.blacklisting.k            <- args$parsimonyBlacklistK
if(is.null(par.blacklisting.k)){
  par.blacklisting.k          <- 0
}
do.dual.blacklisting          <- args$dualBlacklist

output.blacklisting.report    <- args$blacklistReport
if(output.blacklisting.report & !output.files){
  warning("You asked for the blacklist report but also for no output files. No file will be produced.")
}

bl.raw.threshold              <- as.numeric(args$rawBlacklistThreshold)
bl.ratio.threshold            <- as.numeric(args$ratioBlacklistThreshold)

if(do.par.blacklisting & bl.raw.threshold == 0 & bl.ratio.threshold == 0){
  stop("Parsimony blacklisting requested but no thresholds specified with -rwt or -rtt")
}


if(do.dup.blacklisting & bl.raw.threshold == 0 & bl.ratio.threshold == 0){
  stop("Duplicate blacklisting requested but no thresholds specified with -rwt or -rtt")
}

# Slow runs using less memory
useff                         <- args$useff
# if(useff){
#   suppressMessages(require(ff, quietly=TRUE, warn.conflicts=FALSE))
# }

# read and tip count blacklisting
min.reads.per.host            <- args$minReadsPerHost
min.tips.per.host             <- args$minTipsPerHost

post.hoc.min.counts           <- args$postHocCountBlacklisting
# Downsampling
downsampling.limit            <- args$maxReadsPerHost
if(is.null(downsampling.limit)){
  downsampling.limit          <- Inf
}
blacklist.ur                  <- args$blacklistUnderrepresented

# Miscellaneous options
read.counts.matter            <- args$readCountsMatterOnZeroLengthBranches
prune.blacklist               <- args$pruneBlacklist

# What reconstruction mode to use
reconst.mode.arg              <- args$splitsRule
reconst.mode.arg              <- unlist(strsplit(reconst.mode.arg, ","))

if(!(reconst.mode.arg[1] %in% c("r", "s"))){
  stop(paste("Unknown split classifier: ", reconst.mode.arg[1], "\n", sep=""))
}

reconstruction.mode           <- reconst.mode.arg[1]

if(reconstruction.mode == "r" & length(reconst.mode.arg)>1){
  warning("Romero-Severson reconstuction takes no additional arguments; ignoring everything after 'r'")
}

if(reconstruction.mode=="s" & length(reconst.mode.arg)==1){
  stop("One additional argument is required for Sankoff reconstruction")
}

if(reconstruction.mode=="s"){
  sankoff.k                  <- as.numeric(reconst.mode.arg[2])
  sankoff.p                  <- 0
  
  if(is.na(sankoff.k)){
    stop("Expected numerical arguments for reconstruction algorithm parameters.")
  }
} else {
  sankoff.k                  <- NA
  sankoff.p                  <- NA
}

# Whether to output trees
if(args$outputNexusTree){
  tree.output.format         <- "nex"
  if(!output.files){
    warning("You asked for Nexus trees but also for no output files. No trees will be produced.")
  }
} else {
  tree.output.format         <- "pdf"
}

recomb.input                   <- args$recombinationFiles
do.recomb                      <- !is.null(recomb.input)
do.alignment.statistics        <- args$alignmentStatistics

if(do.alignment.statistics & is.null(alignment.input)){
  warning("You must provide alignment files to generate alignment statistics; these will not be produced.")
}

if(do.recomb){
  if(!file.exists(recomb.input)){
    recomb.file.directory      <- dirname(recomb.input)
    recomb.file.regex          <- paste0("^", basename(recomb.input), "(.*)",re.csv.fe,"$")
  } else {
    if(file.info(recomb.input)[['isdir']]){
      recomb.file.directory    <- recomb.input
      recomb.file.regex        <- paste0("^(.*)",re.csv.fe,"$")
    } 
  }
} else {
  recomb.file.directory      <- NULL
  recomb.file.regex          <- NULL
}


# Output classification and collapsed tree files?
do.collapsed                   <- args$collapsedTrees
if(do.collapsed & !output.files){
  warning("You asked for collapsed tree files but also for no output files. No files will be produced.")
}
do.class.detail                <- args$allClassifications
if(do.class.detail & !output.files){
  warning("You asked for classification files but also for no output files. No files will be produced.")
}

multinomial                    <- args$multinomial

# Thresholds for summaries
win.threshold                  <- args$windowThreshold
dist.thresholds                <- as.list(args$distanceThreshold)
if(length(dist.thresholds) == 1){
  if(!multinomial){
    if(dist.thresholds[[1]] == -1){
      dist.threshold <- Inf
    } else {
      dist.threshold <- dist.thresholds[[1]]
    }
  } else {
    if(dist.thresholds[[1]] == -1){
      stop("A finite distance threshold must be specified for the multinomial model")
    } else {
      close.threshold <- dist.thresholds[[1]]
      distant.threshold <- dist.thresholds[[1]]
    }
  } 
} else {
  if(!multinomial){
    warning("More arguments than expected to --distanceThreshold. Ignoring all but the first.")
    if(dist.thresholds[[1]] == -1){
      dist.threshold <- Inf
    } else {
      dist.threshold <- dist.thresholds[[1]]
    }
  } else {
    if(length(dist.thresholds) > 2){
      warning("More arguments than expected to --distanceThreshold. Ignoring all but the first two")
    } 
    close.threshold <- dist.thresholds[[1]]
    distant.threshold <- dist.thresholds[[2]]
    if(close.threshold == -1 | distant.threshold == -1){
      stop("Finite distance thresholds must be specified for the multinomial model")
    }
    if(close.threshold > distant.threshold){
      stop("Threshold for close pairs must be smaller than or equal to threshold for distant pairs.")
    }
  }
}

arrow.threshold                <- args$directionThreshold

if(arrow.threshold >= win.threshold){
  stop("Direction threshold cannot be larger than window threshold")
}

allow.mt                       <- args$allowMultiTrans
relaxed.ancestry               <- args$relaxedAncestry

# Do the simplified plot?
do.simplified.graph            <- !args$skipSummaryGraph
simp.plot.dim                  <- args$summaryPlotDimensions

if(single.tree){
  phyloscanner.trees <- phyloscanner.analyse.tree(
    tree.input,
    reconstruction.mode,
    sankoff.k,
    sankoff.unassigned.switch.threshold = 0,
    continuation.unassigned.proximity.cost = 1000,
    outgroup.name,
    m.thresh,
    is.na(m.thresh),
    blacklist.input,
    dup.input.file.name,
    recomb.input,
    alignment.input,
    alignment.format,
    tip.regex,
    file.name.regex,
    seed,
    norm.ref.file.name,
    norm.standardise.gp,
    norm.constants.input,
    allow.mt,
    relaxed.ancestry,
    par.blacklisting.k,
    bl.raw.threshold,
    bl.ratio.threshold,
    do.dual.blacklisting,
    downsampling.limit,
    blacklist.ur,
    ifelse(post.hoc.min.counts, 1, min.reads.per.host),
    ifelse(post.hoc.min.counts, 1, min.tips.per.host),
    useff,
    prune.blacklist,
    read.counts.matter,
    verbosity,
    no.progress.bars)
} else {	
  phyloscanner.trees <- phyloscanner.analyse.trees(
    tree.directory,
    tree.file.regex,
    reconstruction.mode,
    sankoff.k,
    sankoff.unassigned.switch.threshold = 0,
    continuation.unassigned.proximity.cost = 1000,
    outgroup.name,
    m.thresh,
    is.na(m.thresh),
    user.blacklist.directory,
    user.blacklist.file.regex,
    duplicate.file.directory,
    duplicate.file.regex,
    recomb.file.directory,
    recomb.file.regex,
    alignment.directory,
    alignment.file.regex,
    alignment.format,
    tip.regex,
    file.name.regex,
    seed,
    norm.ref.file.name,
    norm.standardise.gp,
    norm.constants.input,
    allow.mt,
    relaxed.ancestry,
    par.blacklisting.k,
    bl.raw.threshold,
    bl.ratio.threshold,
    do.dual.blacklisting,
    downsampling.limit,
    blacklist.ur,
    ifelse(post.hoc.min.counts, 1, min.reads.per.host),
    ifelse(post.hoc.min.counts, 1, min.tips.per.host),
    useff,
    prune.blacklist,
    read.counts.matter,
    verbosity,
    no.progress.bars)
}

hosts <- all.hosts.from.trees(phyloscanner.trees)

if(length(hosts)<=1 & (do.collapsed | do.class.detail)){
  do.collapsed <- F
  do.class.detail <- F
}

if(verbosity!=0 & output.files){
  cat("Writing annotated trees in .",tree.output.format," format...\n", sep="")
}

if(output.files){
  silent <- sapply(phyloscanner.trees, function(ptree){
    if(single.tree){
      file.name <- paste0(output.string, "_processedTree.", tree.output.format)
    } else {
      file.name <- paste0(output.string, "_processedTree_", ptree$id, ".", tree.output.format)
    }
    write.annotated.tree(ptree, file=file.path(output.dir, file.name), tree.output.format, pdf.scale.bar.width, pdf.w, pdf.hm, verbosity == 2)
  }, simplify = F, USE.NAMES = T)
}

if(do.collapsed & output.files){
  if(verbosity!=0){
    cat("Writing collapsed trees to ",csv.fe," files...\n", sep="")
  }
  
  silent <- sapply(phyloscanner.trees, function(ptree){
    if(single.tree){
      file.name <- paste0(output.string, "_collapsedTree", csv.fe)
    } else {
      file.name <- paste0(output.string, "_collapsedTree_", ptree$id, csv.fe)
    }
    if(verbosity==2) cat("Writing collapsed tree for tree ID ",ptree$id," to file ",file.name, "...\n", sep="")
    
    write_csv(ptree$classification.results$collapsed[,1:4], file.path(output.dir, file.name))
    
  }, simplify = F, USE.NAMES = T)
}


if(do.class.detail & output.files){
  if(verbosity!=0){
    cat("Writing host pairwise classifications to ",csv.fe," files...\n", sep="")
  }
  
  silent <- sapply(phyloscanner.trees, function(ptree){
    if(single.tree){
      file.name <- paste0(output.string, "_classification", csv.fe)
    } else {
      file.name <- paste0(output.string, "_classification_", ptree$id, csv.fe)
    }
    if(verbosity==2) cat("Writing relationship classifications for tree ID ",ptree$id," to file ",file.name, "...\n", sep="")
    write_csv(ptree$classification.results$classification, file.path(output.dir, file.name))
    
  }, simplify = F, USE.NAMES = T)
}

if(length(phyloscanner.trees)>1){
  
  if(verbosity!=0){
    cat("Calculating per-host summary statistics...\n", sep="")
  }
  
  summary.stats <- gather.summary.statistics(phyloscanner.trees, tip.regex = tip.regex, do.alignment.stats = do.alignment.statistics, verbose = verbosity==2)
  
  ss.csv.fn <- paste0(output.string,"_patStats", csv.fe)
  
  if(verbosity!=0 & output.files){
    cat("Writing summary statistics to file ",ss.csv.fn,"...\n", sep="")
  }
  
  if(output.files){
    write_csv(summary.stats, file.path(output.dir, ss.csv.fn))
  }
  
  ss.graphs.fn <- paste0(output.string,"_patStats.pdf")
  
  if(output.files){
    if(verbosity!=0){
      cat("Graphing summary statistics to file ",ss.graphs.fn,"...\n", sep="")
    }
    silent <- multipage.summary.statistics(phyloscanner.trees, summary.stats, file.name = file.path(output.dir, ss.graphs.fn), verbose = verbosity==2)
  }
  
  hosts <- all.hosts.from.trees(phyloscanner.trees)
  
  if(length(hosts)>1){
    
    if(!multinomial){
      if(post.hoc.min.counts){
        ts <- transmission.summary(phyloscanner.trees, win.threshold, dist.threshold, tip.regex, min.tips.per.host, min.reads.per.host, close.sib.only = F, verbosity==2)
      } else {
        ts <- transmission.summary(phyloscanner.trees, win.threshold, dist.threshold, tip.regex, 1, 1, close.sib.only = F, verbosity==2)
      }
      
      
      if (output.files){
        if (verbosity!=0) cat('Writing transmission summary to file ', paste0(output.string,"_hostRelationshipSummary", csv.fe),'...\n', sep="")
        write_csv(ts, file.path(output.dir, paste0(output.string,"_hostRelationshipSummary", csv.fe)))
      }
      
    } else {
      relationship.types	<- c('proximity.3.way',
                              'any.ancestry',
                              'close.x.contiguous',			
                              'close.and.adjacent',					
                              'close.and.adjacent.and.directed',					
                              'close.and.adjacent.and.ancestry',			
                              'close.and.contiguous',					
                              'close.and.contiguous.and.directed',					
                              'close.and.contiguous.and.ancestry')	
      
      dwin <- classify.pairwise.relationships(phyloscanner.trees, 
                                              close.threshold = close.threshold, 
                                              distant.threshold = distant.threshold,
                                              relationship.types = relationship.types, 
                                              verbosity == 2)			
      
      if(post.hoc.min.counts){
        tmp <- select.windows.by.read.and.tip.count(phyloscanner.trees, dwin, tip.regex, min.reads.per.host, min.tips.per.host)
      } else {
        tmp <- select.windows.by.read.and.tip.count(phyloscanner.trees, dwin, tip.regex, 1, 1)
      }
      dc  <- count.pairwise.relationships(tmp) 
      
      if (verbosity!=0 & output.files) cat('Writing multinomial model output to file ', paste0(output.string,"_multinomialOutput", csv.fe),'...\n', sep="")
      write_csv(dc, file.path(output.dir, paste0(output.string,"_multinomialOutput", csv.fe)))
    }
    
    if(do.simplified.graph){

      if(!multinomial){
        if(nrow(ts)==0){
          cat("No relationships exist in the required proportion of windows (",win.threshold,"); skipping simplified relationship summary.\n", sep="")
        } else {
          if (verbosity!=0 & output.files) cat('Drawing simplified summary diagram to file ', paste0(output.string,"_simplifiedRelationshipGraph.pdf"),'...\n', sep="")
          
          simplified.graph <- simplify.summary(ts, arrow.threshold, length(phyloscanner.trees), plot = T)
          
          if(output.files){
            simplified.graph$simp.diagram
            ggsave(file = file.path(output.dir, paste0(output.string,"_simplifiedRelationshipGraph.pdf")), width=simp.plot.dim, height=simp.plot.dim)
          }
        }
      } else {
        if (verbosity!=0 & output.files) cat('Drawing simplified summary diagram to file ', paste0(output.string,"_simplifiedRelationshipGraph.pdf"),'...\n', sep="")
        
        simplified.graph <- simplify.summary.multinomial(dc, win.threshold, arrow.threshold, contiguous = F, plot = T)
        
        simplified.graph$simp.diagram
        ggsave(file = file.path(output.dir, paste0(output.string,"_simplifiedRelationshipGraph.pdf")), width=simp.plot.dim, height=simp.plot.dim)

      }
    }
  }
}

if(output.blacklisting.report & output.files){
  if (verbosity!=0) cat('Saving blacklisting report to file', paste0(output.string,"_blacklistReport",csv.fe),'...\n', sep="")
  
  dfs <- lapply(phyloscanner.trees, function(x) {
    treebl.df <- x$bl.report
    if(nrow(treebl.df)>0){
      treebl.df$tree.id <- x$id
    } else {
      treebl.df$tree.id <- vector()
    }
    treebl.df <- treebl.df[,c(4,1,2,3)]
    treebl.df
  })
  
  output.bl.report <- do.call(rbind, dfs)
  
  write_csv(output.bl.report, path = file.path(output.dir, paste0(output.string,"_blacklistReport",csv.fe)))
}

if(output.rda){
  if (verbosity!=0) cat('Saving R workspace image to file ', paste0(output.string,"_workspace.rda"),'...\n', sep="")
  save.image(file=file.path(output.dir, paste0(output.string,"_workspace.rda")))
}

if(file.exists(file.path(output.dir, "Rplots.pdf"))) silent <- file.remove(file.path(output.dir, "Rplots.pdf"))

if(verbosity!=0) cat("Finished.\n")
