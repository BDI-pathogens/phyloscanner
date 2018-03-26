#!/usr/bin/env Rscript

options("warn"=1)

suppressMessages(require(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phyloscannerR, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(data.table, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ape, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phangorn, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ggtree, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phytools, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(network, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(scales, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gtable, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(grid, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gridExtra, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(kimisc, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(GGally, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(sna, quietly=TRUE, warn.conflicts=FALSE))

arg_parser		     <- ArgumentParser()

arg_parser$add_argument("tree", action="store", help="A string that begins the file names (including the path) of all input trees. e.g. path/to/RAxML_bestTree.InWindow_")
arg_parser$add_argument("alignment", action="store", help="A string that begins the file names (including the path) of the alignments. e.g. path/to/AlignedReadsInWindow_")
arg_parser$add_argument("outputString", action="store", help="A string that will be used to label all output files.")

arg_parser$add_argument("-og", "--outgroupName", action="store", help="The name of the tip in the phylogeny/phylogenies to be used as outgroup (if unspecified, trees will be assumed to be already rooted).")
arg_parser$add_argument("-m", "--multifurcationThreshold", help="If specified, short branches in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. If a number, this number will be used as the threshold, with all branches strictly smaller collapsed. If 'g', it will be guessed from the branch lengths and the width of the genomic window (if appropriate). It is recommended that trees are examined by eye to check that they do appear to have multifurcations if 'g' is used.")
arg_parser$add_argument("-b", "--blacklist", action="store", help="A path and string that begins all the file names for pre-existing blacklist files.")

# General, bland options

arg_parser$add_argument("-od", "--outputDir", action="store", help="All output will be written to this directory. If absent, current working directory.")
arg_parser$add_argument("-v", "--verbose", action="store", default=1, help="The level of verbosity in command line output. 0=none, 1=minimal, 2=detailed")
arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. This expects up to three capture groups, for host ID, read ID, and read count (in that order). If the latter two groups are missing then read information will not be used. If absent, the default is '^(.*)_read_([0-9]+)_count_([0-9]+)$', which matches input from the phyloscanner pipeline where the host ID is the BAM file name.")
arg_parser$add_argument("-y", "--fileNameRegex", action="store", default="^\\D*([0-9]+)_to_([0-9]+)\\D*$", help="Regular expression identifying window coordinates. Two capture groups: start and end; if the latter is missing then the first group is a single numerical identifier for the window. If absent, input will be assumed to be from the phyloscanner pipeline, and the host ID will be the BAM file name.")
arg_parser$add_argument("-tfe", "--treeFileExtension", action="store", default="tree", help="The file extension for tree files (default tree).")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default csv).")
arg_parser$add_argument("-sd", "--seed", action="store", help="Random number seed; used by the downsampling process, and also ties in some parsimony reconstructions can be broken randomly.")
arg_parser$add_argument("-ow", "--overwrite", action="store_true", help="Overwrite existing output files with the same names.")

# Normalisation options

arg_parser$add_argument("-nr", "--normRefFileName", action="store", help="Name of a file giving a normalisation constant for every genome position. Cannot be used simultaneously with -nc. If both -nr and -nc are absent then no normalisation will be performed.")
arg_parser$add_argument("-ns", "--normStandardiseGagPol", action="store_true", default=FALSE, help="An HIV-specific option: if true, the normalising constants are standardised so that the average on gag+pol equals 1. Otherwise they are standardised so the average on the whole genome equals 1.")
arg_parser$add_argument("-nc", "--normalisationConstants", action="store", help="Either a CSV file listing the file name for each tree (column 1) and the normalisation constant (column 2) or a single numerical normalisation constant to be applied to each window. If both -nr and -nc are absent then no normalisation will be performed.")

# Blacklisting

arg_parser$add_argument("-db", "--duplicateBlacklist", action="store", help="Perform blacklisting for likely contamination between samples in this data set, suggested by exact duplicate reads found in different samples. Use this option to specify a string that begins the file names (including their path) of the duplication data produced by phyloscanner_make_trees.py. For example, 'path/to/DuplicateReadCountsProcessed_InWindow_'. This option must be used in conjunction with the --rawBlacklistThreshold and/or --ratioBlacklistThreshold option.")
arg_parser$add_argument("-pbk", "--parsimonyBlacklistK", action="store", type="double", help="Perform parsimony-based blacklisting for likely contaminant sequences (including those originating from a sample outside the current data set). Use this option to specify the value of the within-host diversity penalty used (corresponding to the k parameter in the Sankoff option of the splitsRule argument). This option must be used in conjunction with the --rawBlacklistThreshold and/or --ratioBlacklistThreshold option.")
arg_parser$add_argument("-rwt", "--rawBlacklistThreshold", action="store", default=0, help="Used to specify a read count to be used as a raw threshold for blacklisting. --parsimonyBlacklistK and/or --duplicateBlacklist must also be used. If --parsimonyBlacklistK is used, we will blacklist any subgraph with a read count strictly less than this threshold. If --duplicateBlacklist is used, we will black list any duplicate read with a count strictly less than this threshold. The default value of 0 means nothing is blacklisted.")
arg_parser$add_argument("-rtt", "--ratioBlacklistThreshold", action="store", default=0, help="Used to specify a read count ratio (between 0 and 1) to be used as a threshold for blacklisting. --parsimonyBlacklistK and/or --duplicateBlacklist must also be used. If --parsimonyBlacklistK is used, we will blacklist a subgraph if the ratio of its read count to the total read count from the same host is strictly less than this threshold. If --duplicateBlacklist is used, we will black list a duplicate read if the ratio of its count to the count of the duplicate (from another host) is strictly less than this threshold.")
arg_parser$add_argument("-ub", "--dualBlacklist", action="store_true", default=F, help="Blacklist all reads from the minor subgraphs for all hosts established as dual by parsimony blacklisting.")
arg_parser$add_argument("-rcm", "--readCountsMatterOnZeroLengthBranches", default = FALSE, action="store_true", help="If present, read counts at tips will be taken into account in parsimony reconstructions at the parents of zero-length branches.")

# Blacklisting report

arg_parser$add_argument("-blr", "--blacklistReport", action="store_true", help="If present, output a CSV file of blacklisted tips from each tree.")

# Downsampling

arg_parser$add_argument("-dsl", "--maxReadsPerHost", action="store", type="integer", help="If given, blacklist to downsample read counts from each host to this number.")
arg_parser$add_argument("-dsb", "--blacklistUnderrepresented", action="store_true", help="If present and -dsl is given, blacklist hosts from trees where their total tip count does not reach the maximum.")

args                          <- arg_parser$parse_args()

verbosity                     <- args$verbose

overwrite                     <- args$overwrite

tree.input                    <- args$tree

tree.fe                       <- args$treeFileExtension
csv.fe                        <- args$csvFileExtension

tree.directory                <- dirname(tree.input)
tree.file.regex               <- paste0("^", basename(tree.input), "(.*)\\.",tree.fe,"$")

alignment.input               <- args$alignment

alignment.directory           <- dirname(alignment.input)
alignment.file.regex          <- paste0("^", basename(alignment.input), "(.*)\\.[A-Za-z]+$")

blacklist.input               <- args$blacklist

output.blacklisting.report    <- args$blacklistReport

if(!is.null(blacklist.input)){
  user.blacklist.directory    <- dirname(blacklist.input)
  user.blacklist.file.regex   <- paste0("^", basename(blacklist.input), "(.*)\\.",csv.fe,"$")
} else {
  user.blacklist.directory    <- NULL
  user.blacklist.file.regex   <- NULL
}


output.dir                    <- args$outputDir
if(is.null(output.dir)){
  output.dir                  <- getwd()
}
output.string                 <- args$outputString

if(!overwrite & length(list.files(path=output.dir, pattern=paste0("^", output.string, ".*")))>0 ){
  stop("Previous output with this output string (",output.string,") detected. Please re-run with --overwrite if you wish to overwrite this.")
}

outgroup.name                 <- args$outgroupName

if(is.null(outgroup.name)){
  warning("No outgroup name provided. Trees are assumed to be correctly rooted.")
}

use.m.thresh                  <- !is.null(args$multifurcationThreshold)

seed                          <- args$seed

if(is.null(seed)){
  seed                        <- sample.int(1000000000, 1)
} else {
  seed                        <- as.numeric(seed)
}
if(verbosity!=0) cat("Random number seed is",seed,"\n")
set.seed(seed)

if(use.m.thresh){
  if(args$multifurcationThreshold=="g"){
    m.thresh                  <- NA
  } else if(!is.na(as.numeric(args$multifurcationThreshold))){
    m.thresh                  <- as.numeric(args$multifurcationThreshold)
  } else {
    stop("Unknown argument for --multifurcationThreshold specified\n")
  }
} else {
  m.thresh                    <- -1       
}
norm.ref.file.name            <- args$normRefFileName
norm.standardise.gp           <- args$normStandardiseGagPol
norm.constants.input          <- args$normalisationConstants

tip.regex                     <- args$tipRegex
file.name.regex               <- args$fileNameRegex

if(!is.null(norm.ref.file.name) & !is.null(norm.constants.input)){
  warning("Normalisation reference file name and predetermined normalisation constants both specified. Only the latter will be used.")
}

do.dup.blacklisting           <- !is.null(args$duplicateBlacklist)
dup.input.file.name           <- args$duplicateBlacklist

if(!is.null(dup.input.file.name)){
  duplicate.file.directory    <- dirname(dup.input.file.name)
  duplicate.file.regex        <- paste0("^", basename(dup.input.file.name), "(.*)\\.[A-Za-z]+$")
} else {
  duplicate.file.directory    <- NULL
  duplicate.file.regex        <- NULL
}


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

downsample            <- !is.null(args$maxReadsPerHost)
downsampling.limit    <- args$maxReadsPerHost
if(is.null(downsampling.limit)){
  downsampling.limit <- Inf
}
blacklist.ur          <- args$blacklistUnderrepresented

read.counts.matter    <- args$readCountsMatterOnZeroLengthBranches

ptrees <- phyloscanner.generate.blacklist(
  tree.directory,
  tree.file.regex,
  outgroup.name,
  m.thresh,
  is.na(m.thresh),
  user.blacklist.directory,
  user.blacklist.file.regex,
  duplicate.file.directory,
  duplicate.file.regex,
  alignment.directory,
  alignment.file.regex,
  tip.regex,
  file.name.regex,
  seed,
  norm.ref.file.name,
  norm.standardise.gp,
  norm.constants.input,
  par.blacklisting.k,
  bl.raw.threshold,
  bl.ratio.threshold,
  do.dual.blacklisting,
  downsampling.limit,
  blacklist.ur,
  read.counts.matter,
  verbosity)


silent <- sapply(ptrees, function(ptree){
  
  if(!is.null(ptree$alignment.file.name)){
    
    seqs     <- read.dna(ptree$alignment.file.name, format="fasta")
    new.seqs <- seqs[which(!(labels(seqs) %in% ptree$original.tip.labels[ptree$blacklist])),]
    new.afn  <- paste0(output.string, "_", ptree$id, ".fasta")
    
    if(verbosity!=0) cat("Writing cleaned alignment to ",new.afn, "\n", sep="")
    write.dna(new.seqs, new.afn, format="fasta")

  } else {
    if(verbosity!=0) cat("No alignment file found for tree ID ",ptree$suffix, "\n", sep="")
  }
  
})

if(output.blacklisting.report){
  if (verbosity!=0) cat('Saving blacklisting report to file', paste0(output.string,"_blacklistReport.csv"),'\n')
  
  dfs <- lapply(ptrees, function(x) {
    window.df <- x$bl.report
    if(nrow(window.df)>0){
      window.df$id <- x$id
      window.df <- window.df[,c(3,1,2)]
    }
    window.df
  }
  )
  
  output.bl.report <- do.call(rbind, dfs)
  
  write.csv(output.bl.report, file = file.path(output.dir, paste0(output.string,"_blacklistReport.csv")), quote=F, row.names=F)
}

