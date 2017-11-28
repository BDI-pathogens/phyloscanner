#' Perform a phyloscanner analysis on a single tree
#'
#' This function performs a parsimony reconstruction and classification of pairwise host relationships.
#' @param tree.file.name The name of the tree file (Newick or NEXUS format).
#' @param splits.rule The rules by which the sets of hosts are split into groups in order to ensure that all groups can be members of connected subgraphs without causing conflicts. Options: s=Sankoff with optional within-host diversity penalty (slow, rigorous, recommended), r=Romero-Severson (quick, less rigorous with >2 hosts), f=Sankoff with continuation costs (experimental).
#' @param sankoff.k For \code{splits.rule = s} or \code{f} only. The \emph{k} parameter in the Sankoff reconstruction, representing the within-host diversity penalty.
#' @param sankoff.unassigned.switch.threshold For \code{splits.rule = s} only. Threshold at which a lineage reconstructed as infecting a host will transition to the unassigned state, if it would be equally parsimonious to remain in that host.
#' @param continuation.unassigned.proximity.cost For \code{splits.rule = f} only. The branch length at which an node is reconstructed as unassigned if all its neighbouring nodes are a greater distance away. The default is 1000, intended to be effectively infinite, such a node will never normally receive the unassigned state.
#' @param outgroup.name The name of the tip in the phylogeny/phylogenies to be used as outgroup (if unspecified, trees will be assumed to be already rooted). This should be sufficiently distant to any sequence obtained from a host that it can be assumed that the MRCA of the entire tree was not a lineage present in any sampled individual.
#' @param multifurcation.threshold If specified, branches shorter than this in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. 
#' @param guess.multifurcation.threshold Whether to guess the multifurcation threshold from the branch lengths of the trees and the width of the genomic window (if that information is available). It is recommended that trees are examined by eye to check that they do appear to have multifurcations if using this option.
#' @param user.blacklist.file.name The path of a text file containing the user-specified list of tips to be blacklisted
#' @param duplicate.file.name The path of a .csv file specifying which tree tips are from duplicate reads. Normally this is produced by \code{phyloscanner_make_trees.py}.
#' @param recombination.file.name The path for file containing the results of the \code{phyloscanner_make_trees.py} recombination metric analysis.
#' @param tip.regex Regular expression identifying tips from the dataset. This expects up to three capture groups, for host ID, read ID, and read count (in that order). If the latter two groups are missing then read information will not be used. The default matches input from the phyloscanner pipeline where the host ID is the BAM file name.
#' @param file.name.regex Regular expression identifying window coordinates. Two capture groups: start and end; if the latter is missing then the first group is a single numerical identifier for the window. The default matches input from the phyloscanner pipeline.
#' @param seed Random number seed; used by the downsampling process, and also ties in some parsimony reconstructions can be broken randomly.
#' @param norm.ref.file.name Name of a file giving a normalisation constant for every genome position. Cannot be used simultaneously with \code{norm.constants}. If neither is given then no normalisation will be performed.
#' @param norm.standardise.gag.pol Use only if \code{norm.ref.file.name} is given. An HIV-specific option: if true, the normalising constants are standardised so that the average on gag+pol equals 1. Otherwise they are standardised so the average on the whole genome equals 1.
#' @param norm.constants Either the path of a CSV file listing the file name for each tree (column 1) and the respective normalisation constant (column 2) or a single numerical normalisation constant to be applied to every tree. Cannot be used simultaneously with \code{norm.ref.file.name}. If neither is given then no normalisation will be performed.
#' @param parsimony.blacklist.k The \emph{k} parameter of the single-host Sankhoff parsimony reconstruction used to identify probable contaminants. A value of 0 is equivalent to not performing parsimony blacklisting. 
#' @param raw.blacklist.threshold Used to specify a read count to be used as a raw threshold for duplicate or parsimony blacklisting. Use with \code{parsimony.blacklist.k} or \code{duplicate.file.regex} or both. Parsimony blacklisting will blacklist any subgraph with a read count strictly less than this threshold. Duplicate blacklisting will black list any duplicate read with a count strictly less than this threshold. The default value of 0 means nothing is blacklisted.
#' @param ratio.blacklist.threshold Used to specify a read count ratio (between 0 and 1) to be used as a threshold for duplicate or parsimony blacklisting. Use with \code{parsimony.blacklist.k} or \code{duplicate.file.regex} or both. Parsimony blacklisting will blacklist a subgraph if the ratio of its read count to the total read count from the same host is strictly less than this threshold. Duplcate blacklisting will blacklist a duplicate read if the ratio of its count to the count of the duplicate (from another host) is strictly less than this threshold.
#' @param do.dual.blacklisting Blacklist all reads from the minor subgraphs for all hosts established as dual by parsimony blacklisting (which must have been done for this to do anything).
#' @param max.reads.per.host Used to turn on downsampling. If given, reads will be blacklisted such that read counts (or tip counts if no read counts are identified) from each host are equal (although see \code{blacklist.underrepresented}.
#' @param blacklist.underrepresented If TRUE and \code{max.reads.per.host} is given, blacklist hosts from trees where their total tip count does not reach the maximum.
#' @param use.ff Use the \code{ff} package to store parsimony reconstruction matrices. Use if you run out of memory.
#' @param prune.blacklist If TRUE, all blacklisted and reference tips (except the outgroup) are pruned away before starting parsimony-based reconstruction.
#' @param read.counts.matter.on.zero.length.tips If TRUE, read counts on tips will be taken into account in parsimony reconstructions at the parents of zero-length terminal branches. Not applicable for the Romero-Severson-like reconstruction method.
#' @param verbose Give verbose output.
#' @param no.progress.bars Hide the progress bars from verbose output.
#' @return A list of class \code{phyloscanner.trees} with a single item of class \code{phyloscanner.tree}.
#' @importFrom ape read.tree di2multi root node.depth.edgelength
#' @importFrom data.table data.table as.data.table set setnames
#' @importFrom ff ff
#' @importFrom phangorn Ancestors Descendants Children mrca.phylo getRoot
#' @export phyloscanner.analyse.tree

phyloscanner.analyse.tree <- function(
  tree.file.name,
  splits.rule = c("s", "r", "f"),
  sankoff.k = 0,
  sankoff.unassigned.switch.threshold = 0,
  continuation.unassigned.proximity.cost = 1000,
  outgroup.name = NULL,
  multifurcation.threshold = -1,
  guess.multifurcation.threshold = F,
  user.blacklist.file.name = NULL,
  duplicate.file.name = NULL,
  recombination.file.name = NULL,
  tip.regex = "^(.*)_read_([0-9]+)_count_([0-9]+)$",
  file.name.regex = "^\\D*([0-9]+)_to_([0-9]+)\\D*$",
  seed = sample(1:10000000,1),
  norm.ref.file.name = NULL,
  norm.standardise.gag.pol = F,
  norm.constants = NULL,
  parsimony.blacklist.k = 0,
  raw.blacklist.threshold = 0,
  ratio.blacklist.threshold = 0,
  do.dual.blacklisting = F,
  max.reads.per.host = Inf,
  blacklist.underrepresented = F,
  use.ff = F,
  prune.blacklist = F,
  read.counts.matter.on.zero.length.tips = T,
  verbose = F,
  no.progress.bars = F){
  
  tree.directory <- dirname(tree.file.name)
  user.blacklist.directory <- if(!is.null(user.blacklist.file.name)) dirname(user.blacklist.file.name) else NULL
  duplicate.file.directory <- if(!is.null(duplicate.file.name)) dirname(duplicate.file.name) else NULL
  recombination.file.directory <- if(!is.null(recombination.file.name)) dirname(recombination.file.name) else NULL
  
  phyloscanner.analyse.trees(tree.directory, basename(tree.file.name), splits.rule, sankoff.k, sankoff.unassigned.switch.threshold,
                             continuation.unassigned.proximity.cost, outgroup.name, multifurcation.threshold, guess.multifurcation.threshold, 
                             user.blacklist.directory, basename(user.blacklist.file.name), duplicate.file.directory, 
                             basename(duplicate.file.name), recombination.file.directory, basename(recombination.file.name), tip.regex, 
                             file.name.regex, seed, norm.ref.file.name, norm.standardise.gag.pol, norm.constants, parsimony.blacklist.k, 
                             raw.blacklist.threshold, ratio.blacklist.threshold, do.dual.blacklisting, max.reads.per.host, 
                             blacklist.underrepresented, use.ff, prune.blacklist, read.counts.matter.on.zero.length.tips, verbose, no.progress.bars)
}

#' Perform a phyloscanner analysis on a set of trees
#'
#' This function performs a parsimony reconstruction and classification of pairwise host relationships.
#' @param tree.directory The directory containing all input trees.
#' @param tree.file.regex A regular expression identifying every file in \code{tree.directory} that is to be included in the analysis. The first capture group, if present, gives a unique string identifying each tree. If this is NULL then \code{phyloscanner} will attempt to open every file in \code{tree.directory}.
#' @param splits.rule The rules by which the sets of hosts are split into groups in order to ensure that all groups can be members of connected subgraphs without causing conflicts. Options: s=Sankoff with optional within-host diversity penalty (slow, rigorous, recommended), r=Romero-Severson (quick, less rigorous with >2 hosts), f=Sankoff with continuation costs (experimental).
#' @param sankoff.k For \code{splits.rule = s} or \code{f} only. The \emph{k} parameter in the Sankoff reconstruction, representing the within-host diversity penalty.
#' @param sankoff.unassigned.switch.threshold For \code{splits.rule = s} only. Threshold at which a lineage reconstructed as infecting a host will transition to the unassigned state, if it would be equally parsimonious to remain in that host.
#' @param continuation.unassigned.proximity.cost For \code{splits.rule = f} only. The branch length at which an node is reconstructed as unassigned if all its neighbouring nodes are a greater distance away. The default is 1000, intended to be effectively infinite, such a node will never normally receive the unassigned state.
#' @param outgroup.name The name of the tip in the phylogeny/phylogenies to be used as outgroup (if unspecified, trees will be assumed to be already rooted). This should be sufficiently distant to any sequence obtained from a host that it can be assumed that the MRCA of the entire tree was not a lineage present in any sampled individual.
#' @param multifurcation.threshold If specified, branches shorter than this in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. 
#' @param guess.multifurcation.threshold Whether to guess the multifurcation threshold from the branch lengths of the trees and the width of the genomic window (if that information is available). It is recommended that trees are examined by eye to check that they do appear to have multifurcations if using this option.
#' @param user.blacklist.directory An optional path for a folder containing pre-existing blacklist files. These tips are specified by the user to be excluded from the analysis.
#' @param user.blacklist.file.regex A regular expression identifying every file in \code{user.blacklist.directory} that contains a blacklist. If a capture group is specified then its contents will uniquely identify the tree it belongs to, which must matches the IDs found by \code{tree.file.regex}. If these IDs cannot be identified then matching will be attempted using genome window coordinates.
#' @param duplicate.file.directory An optional path for a folder containing information on duplicate reads, to be used for duplicate blacklisting. Normally this is produced by \code{phyloscanner_make_trees.py}.
#' @param duplicate.file.regex A regular expression identifying every file in \code{duplicate.file.directory} that contains a duplicates file. If a capture group is specified then its contents will uniquely identify the tree it belongs to, which must matches the IDs found by \code{tree.file.regex}. If these IDs cannot be identified then matching will be attempted using genome window coordinates.
#' @param recombination.file.directory An optional path for a folder containing results of the \code{phyloscanner_make_trees.py} recombination metric analysis.
#' @param recombination.file.regex A regular expression identifying every file in \code{recombination.file.directory} that contains a recombination file. If a capture group is specified then its contents will uniquely identify the tree it belongs to, which must matches the IDs found by \code{tree.file.regex}. If these IDs cannot be identified then matching will be attempted using genome window coordinates.
#' @param tip.regex Regular expression identifying tips from the dataset. This expects up to three capture groups, for host ID, read ID, and read count (in that order). If the latter two groups are missing then read information will not be used. The default matches input from the phyloscanner pipeline where the host ID is the BAM file name.
#' @param file.name.regex Regular expression identifying window coordinates. Two capture groups: start and end; if the latter is missing then the first group is a single numerical identifier for the window. The default matches input from the phyloscanner pipeline.
#' @param seed Random number seed; used by the downsampling process, and also ties in some parsimony reconstructions can be broken randomly.
#' @param norm.ref.file.name Name of a file giving a normalisation constant for every genome position. Cannot be used simultaneously with \code{norm.constants}. If neither is given then no normalisation will be performed.
#' @param norm.standardise.gag.pol Use only if \code{norm.ref.file.name} is given. An HIV-specific option: if true, the normalising constants are standardised so that the average on gag+pol equals 1. Otherwise they are standardised so the average on the whole genome equals 1.
#' @param norm.constants Either the path of a CSV file listing the file name for each tree (column 1) and the respective normalisation constant (column 2) or a single numerical normalisation constant to be applied to every tree. Cannot be used simultaneously with \code{norm.ref.file.name}. If neither is given then no normalisation will be performed.
#' @param parsimony.blacklist.k The \emph{k} parameter of the single-host Sankhoff parsimony reconstruction used to identify probable contaminants. A value of 0 is equivalent to not performing parsimony blacklisting. 
#' @param raw.blacklist.threshold Used to specify a read count to be used as a raw threshold for duplicate or parsimony blacklisting. Use with \code{parsimony.blacklist.k} or \code{duplicate.file.regex} or both. Parsimony blacklisting will blacklist any subgraph with a read count strictly less than this threshold. Duplicate blacklisting will black list any duplicate read with a count strictly less than this threshold. The default value of 0 means nothing is blacklisted.
#' @param ratio.blacklist.threshold Used to specify a read count ratio (between 0 and 1) to be used as a threshold for duplicate or parsimony blacklisting. Use with \code{parsimony.blacklist.k} or \code{duplicate.file.regex} or both. Parsimony blacklisting will blacklist a subgraph if the ratio of its read count to the total read count from the same host is strictly less than this threshold. Duplcate blacklisting will blacklist a duplicate read if the ratio of its count to the count of the duplicate (from another host) is strictly less than this threshold.
#' @param do.dual.blacklisting Blacklist all reads from the minor subgraphs for all hosts established as dual by parsimony blacklisting (which must have been done for this to do anything).
#' @param max.reads.per.host Used to turn on downsampling. If given, reads will be blacklisted such that read counts (or tip counts if no read counts are identified) from each host are equal (although see \code{blacklist.underrepresented}.
#' @param blacklist.underrepresented If TRUE and \code{max.reads.per.host} is given, blacklist hosts from trees where their total tip count does not reach the maximum.
#' @param use.ff Use the \code{ff} package to store parsimony reconstruction matrices. Use if you run out of memory.
#' @param prune.blacklist If TRUE, all blacklisted and reference tips (except the outgroup) are pruned away before starting parsimony-based reconstruction.
#' @param read.counts.matter.on.zero.length.tips If TRUE, read counts on tips will be taken into account in parsimony reconstructions at the parents of zero-length terminal branches. Not applicable for the Romero-Severson-like reconstruction method.
#' @param verbose Give verbose output.
#' @param no.progress.bars Hide the progress bars from verbose output.
#' @return A list of class \code{phyloscanner.trees}.
#' @importFrom ape read.tree read.nexus di2multi root node.depth.edgelength
#' @importFrom data.table data.table as.data.table set setnames
#' @importFrom ff ff
#' @importFrom phangorn Ancestors Descendants Children mrca.phylo getRoot
#' @export phyloscanner.analyse.trees 

phyloscanner.analyse.trees <- function(
  tree.directory,
  tree.file.regex = "^RAxML_bestTree.InWindow_([0-9]+_to_[0-9]+)\\.tree$",
  splits.rule = c("s", "r", "f"),
  sankoff.k = 0,
  sankoff.unassigned.switch.threshold = 0,
  continuation.unassigned.proximity.cost = 1000,
  outgroup.name = NULL,
  multifurcation.threshold = -1,
  guess.multifurcation.threshold = F,
  user.blacklist.directory = NULL,
  user.blacklist.file.regex = NULL,
  duplicate.file.directory = NULL,
  duplicate.file.regex = "^DuplicateReadCountsProcessed_InWindow_([0-9]+_to_[0-9]+).csv$",
  recombination.file.directory = NULL,
  recombination.file.regex = "^RecombinantReads_InWindow_([0-9]+_to_[0-9]+).csv$",
  tip.regex = "^(.*)_read_([0-9]+)_count_([0-9]+)$",
  file.name.regex = "^\\D*([0-9]+)_to_([0-9]+)\\D*$",
  seed = sample(1:10000000,1),
  norm.ref.file.name = NULL,
  norm.standardise.gag.pol = F,
  norm.constants = NULL,
  parsimony.blacklist.k = 0,
  raw.blacklist.threshold = 0,
  ratio.blacklist.threshold = 0,
  do.dual.blacklisting = F,
  max.reads.per.host = Inf,
  blacklist.underrepresented = F,
  use.ff = F,
  prune.blacklist = F,
  read.counts.matter.on.zero.length.tips = T,
  verbose = F,
  no.progress.bars = F){
  
  splits.rule <- match.arg(splits.rule)
  
  use.m.thresh          <- multifurcation.threshold > 0 | guess.multifurcation.threshold
  
  if(guess.multifurcation.threshold & multifurcation.threshold >= 0){
    stop("Please either specify a multifurcation threshold, ask for a guess, or neither.")
  }
  
  if(!is.null(norm.ref.file.name) & !is.null(norm.constants)){
    stop("Please either ask for calculation of normalisation constants, provide your own, or neither.")
  }
  
  do.dup.blacklisting   <- !is.null(duplicate.file.directory)
  do.par.blacklisting   <- parsimony.blacklist.k > 0
  
  if(do.dual.blacklisting & !do.par.blacklisting){
    warning("Dual blacklisting requires parsimony blacklisting. Turning dual blacklisting off.")
    do.dual.blacklisting <- F
  }
  
  if(splits.rule == "r" & (sankoff.k > 0 | sankoff.unassigned.switch.threshold > 0 | continuation.unassigned.proximity.cost < Inf)){
    warning("Romero-Severson reconstuction has no parameters; specified values will be ignored.")
  }
  
  downsample            <- max.reads.per.host < Inf
  
  existing.bl           <- !is.null(user.blacklist.directory)
  do.recomb             <- !is.null(recombination.file.directory)
  
  # 1. Make the big list.
  
  all.tree.info <- list()
  
  full.tree.file.names <- list.files.mod(tree.directory, pattern=tree.file.regex, full.names=TRUE)
  tree.file.names <- list.files.mod(tree.directory,  pattern=tree.file.regex)
  
  if(length(tree.file.names)==0){
    stop("No tree files found.")
  }
  
  # There are two ways to match files - to the "suffixes" - which no longer literally have to be suffixes - and to window coordinates.
  
  match.mode <<- NA
  
  if(!is.null(tree.file.regex)){
    tree.identifiers <- sapply(tree.file.names, function(x) sub(tree.file.regex, "\\1", x))
    if(all(tree.identifiers!="")){
      match.mode <<- "suffix"
    } 
  } 
  
  if(is.na(match.mode)){
    match.mode <<- "coords"
    tree.identifiers <- tryCatch({
      sapply(tree.file.names, function(x) get.window.coords.string(x, file.name.regex))},
      error = function(e){
        match.mode <<- "none"
        tree.file.names
      })
  }
  
  if(length(unique(tree.identifiers)) != length(tree.identifiers)){
    stop("Some trees have duplicate IDs.")
  }
  
  for(id.no in 1:length(tree.identifiers)){
    suffix                        <- tree.identifiers[id.no] 
    
    tree.info                     <- list()
    class(tree.info)              <- append(class(tree.info), "phyloscanner.tree")
    
    tree.info$suffix              <- suffix
    tree.info$tree.file.name      <- full.tree.file.names[id.no]
    
    all.tree.info[[suffix]]       <- tree.info
  }
  
  # 2. Coordinates.
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) identify.window.coords(tree.info, file.name.regex), simplify = F, USE.NAMES = T)
  
  readable.coords <- all(sapply(all.tree.info, function(x) !is.null(x$window.coords)))
  
  # 3. Attach blacklist and recombination files
  
  if(existing.bl){
    user.blacklist.file.names <- list.files.mod(user.blacklist.directory, pattern=user.blacklist.file.regex)
    full.user.blacklist.file.names <- list.files.mod(user.blacklist.directory, pattern=user.blacklist.file.regex, full.names=TRUE)
    if(length(tree.file.names == 1) & length(user.blacklist.file.names)==1){
      user.blacklist.identifiers <- tree.identifiers
    } else {
      if(match.mode == "suffix"){
        user.blacklist.identifiers <- sapply(user.blacklist.file.names, function(x) sub(user.blacklist.file.regex, "\\1", x))
      } else {
        if(match.mode != "coords"){
          stop("Cannot match blacklist files with tree files using the information given.")
        } else {
          tryCatch({
            user.blacklist.identifiers <- sapply(user.blacklist.file.names, function(x) get.window.coords.string(x, file.name.regex))},
            error = function(e){
              stop("Cannot match blacklist files with tree files using the information given.")
            })
        }
      }
      if(!setequal(tree.identifiers, user.blacklist.identifiers)){
        warning("Tree files and blacklist files do not entirely match. Blacklist files with no tree files will be ignored; tree files with no blacklist file will have no tips blacklisted.")
      }
    }
    
    all.tree.info <- sapply(all.tree.info, function(tree.info) attach.file.names(tree.info, full.user.blacklist.file.names, user.blacklist.identifiers, "user.blacklist.file.name"), simplify = F, USE.NAMES = T)
  }
  
  if(do.recomb){
    recombination.file.names <- list.files.mod(recombination.file.directory, pattern=recombination.file.regex)
    full.recombination.file.names <- list.files.mod(recombination.file.directory, pattern=recombination.file.regex, full.names=TRUE)
    if(length(tree.file.names == 1) & length(recombination.file.names)==1){
      recomb.identifiers <- tree.identifiers
    } else {
      if(match.mode=="suffix"){
        recomb.identifiers <- sapply(recombination.file.names, function(x) sub(recombination.file.regex, "\\1", x))
      } else {
        if(match.mode != "coords"){
          stop("Cannot match recombination files with tree files using the information given.")
        } else {
          tryCatch({
            recomb.identifiers <- sapply(recombination.file.names, function(x) get.window.coords.string(x, file.name.regex))},
            error = function(e){
              stop("Cannot match recombination files with tree files using the information given.")
            })
        }
      }
      if(!setequal(tree.identifiers, recomb.identifiers)){
        stop("Tree files and recombination files do not match.")
      }
    }
    
    all.tree.info <- sapply(all.tree.info, function(tree.info) attach.file.names(tree.info, full.recombination.file.names, recomb.identifiers, "recombination.file.name"), simplify = F, USE.NAMES = T)
  }

  if(do.dup.blacklisting){
    duplicate.file.names <- list.files.mod(duplicate.file.directory, pattern=duplicate.file.regex)
    full.duplicate.file.names <- list.files.mod(duplicate.file.directory, pattern=duplicate.file.regex, full.names=TRUE)
    if(length(tree.file.names == 1) & length(duplicate.file.names)==1){
      duplicate.identifiers <- tree.identifiers
    } else {
      if(match.mode=="suffix"){
        duplicate.identifiers <- sapply(duplicate.file.names, function(x) sub(duplicate.file.regex, "\\1", x))
      } else {
        if(match.mode != "coords"){
          stop("Cannot match duplicate files with tree files using the information given.")
        } else {
          tryCatch({
            duplicate.identifiers <- sapply(duplicate.file.names, function(x) get.window.coords.string(x, file.name.regex))},
            error = function(e){
              stop("Cannot match duplicate files with tree files using the information given.")
            })
        }
      }
      if(!setequal(tree.identifiers, duplicate.identifiers)){
        stop("Tree files and duplicate files do not match.")
      }
    }
    
    all.tree.info <- sapply(all.tree.info, function(tree.info) attach.file.names(tree.info, full.duplicate.file.names, duplicate.identifiers, "duplicate.file.name"), simplify = F, USE.NAMES = T)
  }
  
  # 4. Read the trees
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) attach.tree(tree.info, verbose), simplify = F, USE.NAMES = T)
  
  # 5. Are there read counts at the tips?
  
  read.counts.check <- sapply(all.tree.info, function(tree.info) check.read.counts(tree.info, tip.regex))
  
  if(any(read.counts.check) & !all(read.counts.check)){
    warning("Read counts are not present in some trees; ignoring them throughout.\n")
  }
  has.read.counts <- all(read.counts.check)
  
  # 6. Root tree, collapse multifurcations, get host for each tip
  
  if(guess.multifurcation.threshold){
    warning("Attempting to guess a branch length threshold for multifurcations from the tree. Please visually examine the tree or trees for multifurcations before using the results of this analysis.")
  }
  

  all.tree.info <- sapply(all.tree.info, function(tree.info) prepare.tree(tree.info, outgroup.name, tip.regex, guess.multifurcation.threshold, multifurcation.threshold, verbose), 
                            simplify = F, USE.NAMES = T)
  
  all.tree.info[sapply(all.tree.info, is.null)] <- NULL
  
  if(length(all.tree.info)==0){
    stop("Cannot find any hosts on any tree that match this regular expression. Please check that it is correct.")
  }
  
  
  # 8. Read the blacklists
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) read.blacklist(tree.info, verbose), simplify = F, USE.NAMES = T)
  
  
  # 9. Rename tips from the prexisting blacklist
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) rename.user.blacklist.tips(tree.info), simplify = F, USE.NAMES = T)
  
  
  # 10. Get the normalisation constants
  
  if(!is.null(norm.constants)){
    if(!is.na(suppressWarnings(as.numeric(norm.constants)))){
      nc <- as.numeric(norm.constants)
      all.tree.info <- sapply(all.tree.info, function(tree.info) {
        tree.info$normalisation.constant  <- nc
        tree.info
      }, simplify = F, USE.NAMES = T)
    } else if(file.exists(norm.constants)){
      
      nc.df   <- read.csv(norm.constants, stringsAsFactors = F, header = F)
      all.tree.info <- sapply(all.tree.info, function(tree.info) {
        rows  <- which(nc.df[,1]==basename(tree.info$tree.file.name))
        
        if(length(rows)>0){
          # Take the first appearance of the file name
          
          row <- rows[1]
          nc  <- as.numeric(nc.df[row, 2])
          if(is.na(nc)){
            warning("Tree with suffix ",tree.info$suffix," given a normalisation constant in file ",norm.constants," which cannot be processed as a number; this tree will not have normalised branch lengths\n")
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
      
      if(grepl(paste0("csv", "$"), norm.ref.file.name)){
        norm.table	<- as.data.table(read.csv(norm.ref.file.name, stringsAsFactors=FALSE))
        
        if(ncol(norm.table)!=2){
          stop(paste0(norm.ref.file.name," is not formatted as expected for a normalisation lookup file; expecting two columns.\n"))
        } else {
          setnames(norm.table, 1, 'POSITION')
          setnames(norm.table, 2, 'NORM_CONST')
          
          if(norm.standardise.gag.pol){
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
  
  # 11. Apply the normalisation constants
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) apply.normalisation.constants(tree.info), simplify = F, USE.NAMES = T)
  
  
  # 12. Duplicate blacklisting
  
  if(do.dup.blacklisting){
    all.tree.info <- sapply(all.tree.info, find.duplicate.tips, simplify = F, USE.NAMES = T)
    
    all.tree.info <- sapply(all.tree.info, function(tree.info) blacklist.from.duplicates.vector(tree.info, verbose), simplify = F, USE.NAMES = T)
  }
  
  
  # 13. Parsimony blacklisting
  
  if(do.par.blacklisting){
    all.tree.info <- sapply(all.tree.info, function(tree.info) blacklist.using.parsimony(tree.info, tip.regex, outgroup.name, raw.blacklist.threshold, ratio.blacklist.threshold, parsimony.blacklist.k, has.read.counts, read.counts.matter.on.zero.length.tips, verbose), simplify = F, USE.NAMES = T)
  }
  
  # 14. Dual blacklisting
  
  if(do.dual.blacklisting){
    hosts.that.are.duals <- lapply(all.tree.info, function(tree.info){
      tree.info$duals.info$host
    })
    hosts.that.are.duals <- unique(unlist(hosts.that.are.duals))
    hosts.that.are.duals <- hosts.that.are.duals[order(hosts.that.are.duals)]
    
    dual.results <- blacklist.duals(all.tree.info, hosts.that.are.duals, summary.file = NULL, verbose)
    
    all.tree.info <- sapply(all.tree.info, function(tree.info) blacklist.from.duals.list(tree.info, dual.results, verbose), simplify = F, USE.NAMES = T)
  }
  
  
  # 15. Downsampling
  
  if(downsample){
    all.tree.info <- sapply(all.tree.info, function(tree.info) blacklist.from.random.downsample(tree.info, max.reads.per.host, blacklist.underrepresented, has.read.counts, tip.regex, seed, verbose), simplify = F, USE.NAMES = T)
  }
  
  # 16. All IDs can be safely gathered and mapped now. Windows with no patients should be removed.
  
  if(verbose) cat("Gathering host IDs...\n")
  
  hosts <- all.hosts.from.trees(all.tree.info)
  
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
  
  # 17. Prune away the blacklist if so requested
  
  if(prune.blacklist){
    all.tree.info <- sapply(all.tree.info, function(tree.info) {
      tree                <- tree.info$tree
      
      new.tree            <- process.tree(tree, outgroup.name, m.thresh = -1, tree.info$blacklist)
      
      tree.info$tree      <- new.tree
      tree.info$blacklist <- vector()
      
      tree.info
      
    }, simplify = F, USE.NAMES = T)
  }
  
  
  # 18. Parsimony reconstruction
  
  sankoff.p <- NA
  
  if(splits.rule == "s"){
    sankoff.p <- sankoff.unassigned.switch.threshold
  } else if(splits.rule == "f"){
    sankoff.p <- continuation.unassigned.proximity.cost
  } 
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) {
    # Do the reconstruction
    
    if(verbose) cat("Reconstructing internal node hosts on tree ID ",tree.info$suffix, "\n", sep="")
    
    tmp					     <- split.hosts.to.subgraphs(tree.info$tree, tree.info$blacklist, splits.rule, tip.regex, sankoff.k, sankoff.p, use.ff, read.counts.matter.on.zero.length.tips, multifurcation.threshold, hosts, verbose, no.progress.bars)
    tree					   <- tmp[['tree']]	
    
    # trees are annotated from now on
    
    rs.subgraphs		 <- tmp[['rs.subgraphs']]
    
    if(is.null(rs.subgraphs)){
      warning("No non-blacklisted tips for suffix ",tree.info$suffix,"; this window will not be analysed further.")
    }
    
    # Convert back to unnormalised branch lengths for output
    
    tree$edge.length <- tree$edge.length * tree.info$normalisation.constant
    
    if(!is.null(rs.subgraphs)){		
      
      if(has.read.counts){
        rs.subgraphs$reads <- sapply(rs.subgraphs$tip, function(x) as.numeric(read.count.from.label(x, tip.regex)))
      } else {
        rs.subgraphs$reads <- 1
      }
      # 
      tree.info$splits.table  <- rs.subgraphs
    }
    
    tree.info$tree <- tree
    
    tree.info
    
  }, simplify = F, USE.NAMES = T)
  
  
  # 19. For summary statistics
  
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) {
    clade.results                 <- resolveTreeIntoPatientClades(tree.info$tree, hosts, tip.regex, tree.info$blacklist, !has.read.counts)
    
    tree.info$clades.by.host      <- clade.results$clades.by.patient
    tree.info$clade.mrcas.by.host <- clade.results$clade.mrcas.by.patient
    
    tree.info
  }, simplify = F, USE.NAMES = T)
  
  # 20. Individual window classifications
  
  if(length(hosts)>1){
    all.tree.info <- sapply(all.tree.info, function(tree.info) {
      
      
      if(verbose) cat("Classifying pairwise host relationships for tree ID ",tree.info$suffix, ".\n", sep="")
      
      tree.info$classification.results <- classify(tree.info, verbose, no.progress.bars)
      
      tree.info
    }, simplify = F, USE.NAMES = T)
  }
  
  attr(all.tree.info, 'readable.coords') <- readable.coords
  attr(all.tree.info, 'match.mode')      <- match.mode
  attr(all.tree.info, 'has.read.counts') <- has.read.counts
  
  class(all.tree.info) <- append(class(all.tree.info), "phyloscanner.trees")
  
  all.tree.info
}

#' Perform blacklisting and clean the input alignments of blacklisted sequences without doing any further phyloscanner analysis
#'
#' This function performs all the blacklisting steps of a phyloscanner analysis and then produces new alignment files with blacklisted sequences removed
#' @param tree.directory The directory containing all input trees.
#' @param alignment.directory The directory containing the alignments. 
#' @param tree.file.regex A regular expression identifying every file in \code{tree.directory} that is to be included in the analysis. The first capture group, if present, gives a unique string identifying each tree. If this is NULL then \code{phyloscanner} will attempt to open every file in \code{tree.directory}.
#' @param alignment.file.regex A regular expression identifying every file in \code{alignment.directory} that is an alignment. If a capture group is specified then its contents will uniquely identify the tree it belongs to, which must matches the IDs found by \code{tree.file.regex}. If these IDs cannot be identified then matching will be attempted using genome window coordinates.
#' @param outgroup.name The name of the tip in the phylogeny/phylogenies to be used as outgroup (if unspecified, trees will be assumed to be already rooted). This should be sufficiently distant to any sequence obtained from a host that it can be assumed that the MRCA of the entire tree was not a lineage present in any sampled individual.
#' @param multifurcation.threshold If specified, branches shorter than this in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. 
#' @param guess.multifurcation.threshold Whether to guess the multifurcation threshold from the branch lengths of the trees and the width of the genomic window (if that information is available). It is recommended that trees are examined by eye to check that they do appear to have multifurcations if using this option.
#' @param user.blacklist.directory An optional path for a folder containing pre-existing blacklist files. These tips are specified by the user to be excluded from the analysis.
#' @param user.blacklist.file.regex A regular expression identifying every file in \code{user.blacklist.directory} that contains a blacklist. If a capture group is specified then its contents will uniquely identify the tree it belongs to, which must matches the IDs found by \code{tree.file.regex}. If these IDs cannot be identified then matching will be attempted using genome window coordinates.
#' @param duplicate.file.directory An optional path for a folder containing information on duplicate reads, to be used for duplicate blacklisting. Normally this is produced by \code{phyloscanner_make_trees.py}.
#' @param duplicate.file.regex A regular expression identifying every file in \code{duplicate.file.directory} that contains a duplicates file. If a capture group is specified then its contents will uniquely identify the tree it belongs to, which must matches the IDs found by \code{tree.file.regex}. If these IDs cannot be identified then matching will be attempted using genome window coordinates.
#' @param recombination.file.directory An optional path for a folder containing results of the \code{phyloscanner_make_trees.py} recombination metric analysis.
#' @param recombination.file.regex A regular expression identifying every file in \code{recombination.file.directory} that contains a recombination file. If a capture group is specified then its contents will uniquely identify the tree it belongs to, which must matches the IDs found by \code{tree.file.regex}. If these IDs cannot be identified then matching will be attempted using genome window coordinates.
#' @param tip.regex Regular expression identifying tips from the dataset. This expects up to three capture groups, for host ID, read ID, and read count (in that order). If the latter two groups are missing then read information will not be used. The default matches input from the phyloscanner pipeline where the host ID is the BAM file name.
#' @param file.name.regex Regular expression identifying window coordinates. Two capture groups: start and end; if the latter is missing then the first group is a single numerical identifier for the window. The default matches input from the phyloscanner pipeline.
#' @param seed Random number seed; used by the downsampling process, and also ties in some parsimony reconstructions can be broken randomly.
#' @param norm.ref.file.name Name of a file giving a normalisation constant for every genome position. Cannot be used simultaneously with \code{norm.constants}. If neither is given then no normalisation will be performed.
#' @param norm.standardise.gag.pol Use only if \code{norm.ref.file.name} is given. An HIV-specific option: if true, the normalising constants are standardised so that the average on gag+pol equals 1. Otherwise they are standardised so the average on the whole genome equals 1.
#' @param norm.constants Either the path of a CSV file listing the file name for each tree (column 1) and the respective normalisation constant (column 2) or a single numerical normalisation constant to be applied to every tree. Cannot be used simultaneously with \code{norm.ref.file.name}. If neither is given then no normalisation will be performed.
#' @param parsimony.blacklist.k The \emph{k} parameter of the single-host Sankhoff parsimony reconstruction used to identify probable contaminants. A value of 0 is equivalent to not performing parsimony blacklisting. 
#' @param raw.blacklist.threshold Used to specify a read count to be used as a raw threshold for duplicate or parsimony blacklisting. Use with \code{parsimony.blacklist.k} or \code{duplicate.file.regex} or both. Parsimony blacklisting will blacklist any subgraph with a read count strictly less than this threshold. Duplicate blacklisting will black list any duplicate read with a count strictly less than this threshold. The default value of 0 means nothing is blacklisted.
#' @param ratio.blacklist.threshold Used to specify a read count ratio (between 0 and 1) to be used as a threshold for duplicate or parsimony blacklisting. Use with \code{parsimony.blacklist.k} or \code{duplicate.file.regex} or both. Parsimony blacklisting will blacklist a subgraph if the ratio of its read count to the total read count from the same host is strictly less than this threshold. Duplcate blacklisting will blacklist a duplicate read if the ratio of its count to the count of the duplicate (from another host) is strictly less than this threshold.
#' @param do.dual.blacklisting Blacklist all reads from the minor subgraphs for all hosts established as dual by parsimony blacklisting (which must have been done for this to do anything).
#' @param max.reads.per.host Used to turn on downsampling. If given, reads will be blacklisted such that read counts (or tip counts if no read counts are identified) from each host are equal (although see \code{blacklist.underrepresented}.
#' @param blacklist.underrepresented If TRUE and \code{max.reads.per.host} is given, blacklist hosts from trees where their total tip count does not reach the maximum.
#' @param output.file.id A string identifying the cleaned alignments
#' @param verbose Give verbose output.
#' @param no.progress.bars Hide the progress bars from verbose output.
#' @importFrom ape read.tree read.nexus di2multi root node.depth.edgelength
#' @importFrom data.table data.table as.data.table set setnames
#' @importFrom phangorn Ancestors Descendants Children mrca.phylo getRoot
#' @export remove.blacklist.from.alignment 

remove.blacklist.from.alignment <- function(
  tree.directory,
  alignment.directory,
  tree.file.regex = "^RAxML_bestTree.InWindow_([0-9]+_to_[0-9]+)\\.tree$",
  alignment.file.regex ="^AlignedReadsInWindow_([0-9]+_to_[0-9]+)\\.fasta$",
  outgroup.name = NULL,
  multifurcation.threshold = -1,
  guess.multifurcation.threshold = F,
  user.blacklist.directory = NULL,
  user.blacklist.file.regex = NULL,
  duplicate.file.directory = NULL,
  duplicate.file.regex = "^DuplicateReadCountsProcessed_InWindow_([0-9]+_to_[0-9]+).csv$",
  tip.regex = "^(.*)_read_([0-9]+)_count_([0-9]+)$",
  file.name.regex = "^\\D*([0-9]+)_to_([0-9]+)\\D*$",
  seed = sample(1:10000000,1),
  norm.ref.file.name = NULL,
  norm.standardise.gag.pol = F,
  norm.constants = NULL,
  parsimony.blacklist.k = 0,
  raw.blacklist.threshold = 0,
  ratio.blacklist.threshold = 0,
  do.dual.blacklisting = F,
  max.reads.per.host = Inf,
  blacklist.underrepresented = F,
  output.file.id = "CleanedAlignment_InWindow_",
  verbose = F,
  no.progress.bars = T){
  
  use.m.thresh          <- multifurcation.threshold > 0 | guess.multifurcation.threshold
  
  if(guess.multifurcation.threshold & multifurcation.threshold >= 0){
    stop("Please either specify a multifurcation threshold, ask for a guess, or neither.")
  }
  
  if(!is.null(norm.ref.file.name) & !is.null(norm.constants)){
    stop("Please either ask for calculation of normalisation constants, provide your own, or neither.")
  }
  
  do.dup.blacklisting   <- !is.null(duplicate.file.directory)
  do.par.blacklisting   <- !is.null(parsimony.blacklist.k)
  if(do.par.blacklisting){
    do.par.blacklisting   <- parsimony.blacklist.k > 0
  }
  
  if(do.dual.blacklisting & !do.par.blacklisting){
    warning("Dual blacklisting requires parsimony blacklisting. Turning dual blacklisting off.")
    do.dual.blacklisting <- F
  }
  
  downsample            <- max.reads.per.host < Inf
  
  existing.bl           <- !is.null(user.blacklist.directory)
  
  # 1. Make the big list.
  
  all.tree.info <- list()
  
  full.tree.file.names <- list.files.mod(tree.directory, pattern=tree.file.regex, full.names=TRUE)
  tree.file.names <- list.files.mod(tree.directory,  pattern=tree.file.regex)
  
  if(length(tree.file.names)==0){
    stop("No tree files found.")
  }
  
  # There are two ways to match files - to the "suffixes" - which no longer literally have to be suffixes - and to window coordinates.
  
  match.mode <<- NA
  
  if(!is.null(tree.file.regex)){
    tree.identifiers <- sapply(tree.file.names, function(x) sub(tree.file.regex, "\\1", x))
    if(all(tree.identifiers!="")){
      match.mode <<- "suffix"
    } 
  } 
  
  if(is.na(match.mode)){
    match.mode <<- "coords"
    tree.identifiers <- tryCatch({
      sapply(tree.file.names, function(x) get.window.coords.string(x, file.name.regex))},
      error = function(e){
        match.mode <<- "none"
        tree.file.names
      })
  }
  
  if(length(unique(tree.identifiers)) != length(tree.identifiers)){
    stop("Some trees have duplicate IDs.")
  }
  
  for(id.no in 1:length(tree.identifiers)){
    suffix                        <- tree.identifiers[id.no] 
    
    tree.info                     <- list()
    class(tree.info)              <- append(class(tree.info), "phyloscanner.tree")
    
    tree.info$suffix              <- suffix
    tree.info$tree.file.name      <- full.tree.file.names[id.no]
    
    all.tree.info[[suffix]]       <- tree.info
  }
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) identify.window.coords(tree.info, file.name.regex), simplify = F, USE.NAMES = T)
  
  readable.coords <- all(sapply(all.tree.info, function(x) !is.null(x$window.coords)))
  
  # 2. Attach files
  
  if(existing.bl){
    user.blacklist.file.names <- list.files.mod(user.blacklist.directory, pattern=user.blacklist.file.regex)
    full.user.blacklist.file.names <- list.files.mod(user.blacklist.directory, pattern=user.blacklist.file.regex, full.names=TRUE)
    if(length(tree.file.names == 1) & length(user.blacklist.file.names)==1){
      user.blacklist.identifiers <- tree.identifiers
    } else {
      if(match.mode == "suffix"){
        user.blacklist.identifiers <- sapply(user.blacklist.file.names, function(x) sub(user.blacklist.file.regex, "\\1", x))
      } else {
        if(match.mode != "coords"){
          stop("Cannot match blacklist files with tree files using the information given.")
        } else {
          tryCatch({
            user.blacklist.identifiers <- sapply(user.blacklist.file.names, function(x) get.window.coords.string(x, file.name.regex))},
            error = function(e){
              stop("Cannot match blacklist files with tree files using the information given.")
            })
        }
      }
      if(!setequal(tree.identifiers, user.blacklist.identifiers)){
        warning("Tree files and blacklist files do not entirely match. Blacklist files with no tree files will be ignored; tree files with no blacklist file will have no tips blacklisted.")
      }
    }
    
    all.tree.info <- sapply(all.tree.info, function(tree.info) attach.file.names(tree.info, full.user.blacklist.file.names, user.blacklist.identifiers, "user.blacklist.file.name"), simplify = F, USE.NAMES = T)
  }
  
  if(do.dup.blacklisting){
    duplicate.file.names <- list.files.mod(duplicate.file.directory, pattern=duplicate.file.regex)
    full.duplicate.file.names <- list.files.mod(duplicate.file.directory, pattern=duplicate.file.regex, full.names=TRUE)
    if(length(tree.file.names == 1) & length(duplicate.file.names)==1){
      duplicate.identifiers <- tree.identifiers
    } else {
      if(match.mode=="suffix"){
        duplicate.identifiers <- sapply(duplicate.file.names, function(x) sub(duplicate.file.regex, "\\1", x))
      } else {
        if(match.mode != "coords"){
          stop("Cannot match duplicate files with tree files using the information given.")
        } else {
          tryCatch({
            duplicate.identifiers <- sapply(duplicate.file.names, function(x) get.window.coords.string(x, file.name.regex))},
            error = function(e){
              stop("Cannot match duplicate files with tree files using the information given.")
            })
        }
      }
      if(!setequal(tree.identifiers, duplicate.identifiers)){
        stop("Tree files and duplicate files do not match.")
      }
    }
    
    all.tree.info <- sapply(all.tree.info, function(tree.info) attach.file.names(tree.info, full.duplicate.file.names, duplicate.identifiers, "duplicate.file.name"), simplify = F, USE.NAMES = T)
  }
  

  alignment.file.names <- list.files.mod(alignment.directory, pattern=alignment.file.regex)
  full.alignment.file.names <- list.files.mod(alignment.directory, pattern=alignment.file.regex, full.names=TRUE)
  if(length(tree.file.names == 1) & length(alignment.file.names)==1){
    alignment.identifiers <- tree.identifiers
  } else {
    if(match.mode == "suffix"){
      alignment.identifiers <- sapply(alignment.file.names, function(x) sub(alignment.file.regex, "\\1", x))
    } else {
      if(match.mode != "coords"){
        stop("Cannot match alignment files with tree files using the information given.")
      } else {
        tryCatch({
          alignment.identifiers <- sapply(alignment.file.names, function(x) get.window.coords.string(x, file.name.regex))},
          error = function(e){
            stop("Cannot match blacklist files with tree files using the information given.")
          })
      }
    }
    if(!setequal(tree.identifiers, alignment.identifiers)){
      warning("Tree files and alignment files do not entirely match. Alignment files with no tree files will be ignored; if no alignment file exists then it will not be cleaned.")
    }
  }
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) attach.file.names(tree.info, full.alignment.file.names, alignment.identifiers, "alignment.file.name"), simplify = F, USE.NAMES = T)
  
  # 3. Read the trees
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) attach.tree(tree.info, verbose), simplify = F, USE.NAMES = T)
  
  # 4. Are there read counts at the tips?
  
  read.counts.check <- sapply(all.tree.info, function(tree.info) check.read.counts(tree.info, tip.regex))
  
  if(any(read.counts.check) & !all(read.counts.check)){
    warning("Read counts are not present in some trees; ignoring them throughout.\n")
  }
  has.read.counts <- all(read.counts.check)
  
  # 5. Root tree, collapse multifurcations, get host for each tip
  
  if(guess.multifurcation.threshold){
    warning("Attempting to guess a branch length threshold for multifurcations from the tree. Please visually examine the tree or trees for multifurcations before using the results of this analysis.")
  }

  all.tree.info <- sapply(all.tree.info, function(tree.info) prepare.tree(tree.info, outgroup.name, tip.regex, guess.multifurcation.threshold, multifurcation.threshold, verbose), 
                            simplify = F, USE.NAMES = T)
  
  all.tree.info[sapply(all.tree.info, is.null)] <- NULL
  
  if(length(all.tree.info)==0){
    stop("Cannot find any hosts on any tree that match this regular expression. Please check that it is correct.")
  }
  
  # 6. Read the blacklists
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) read.blacklist(tree.info, verbose), simplify = F, USE.NAMES = T)
  

  # 7. Get the normalisation constants
  
  if(!is.null(norm.constants)){
    if(!is.na(suppressWarnings(as.numeric(norm.constants)))){
      nc <- as.numeric(norm.constants)
      all.tree.info <- sapply(all.tree.info, function(tree.info) {
        tree.info$normalisation.constant  <- nc
        tree.info
      }, simplify = F, USE.NAMES = T)
    } else if(file.exists(norm.constants)){
      
      nc.df   <- read.csv(norm.constants, stringsAsFactors = F, header = F)
      all.tree.info <- sapply(all.tree.info, function(tree.info) {
        rows  <- which(nc.df[,1]==basename(tree.info$tree.file.name))
        
        if(length(rows)>0){
          # Take the first appearance of the file name
          
          row <- rows[1]
          nc  <- as.numeric(nc.df[row, 2])
          if(is.na(nc)){
            warning("Tree with suffix ",tree.info$suffix," given a normalisation constant in file ",norm.constants," which cannot be processed as a number; this tree will not have normalised branch lengths\n")
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
      
      if(grepl(paste0("csv", "$"), norm.ref.file.name)){
        norm.table	<- as.data.table(read.csv(norm.ref.file.name, stringsAsFactors=FALSE))
        
        if(ncol(norm.table)!=2){
          stop(paste0(norm.ref.file.name," is not formatted as expected for a normalisation lookup file; expecting two columns.\n"))
        } else {
          setnames(norm.table, 1, 'POSITION')
          setnames(norm.table, 2, 'NORM_CONST')
          
          if(norm.standardise.gag.pol){
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
  
  all.tree.info <- sapply(all.tree.info, function(tree.info) apply.normalisation.constants(tree.info), simplify = F, USE.NAMES = T)
  
  
  # 9. Duplicate blacklisting
  
  if(do.dup.blacklisting){
    all.tree.info <- sapply(all.tree.info, find.duplicate.tips, simplify = F, USE.NAMES = T)
    
    all.tree.info <- sapply(all.tree.info, function(tree.info) blacklist.from.duplicates.vector(tree.info, verbose), simplify = F, USE.NAMES = T)
  }
  
  
  # 10. Parsimony blacklisting
  
  if(do.par.blacklisting){
    all.tree.info <- sapply(all.tree.info, function(tree.info) blacklist.using.parsimony(tree.info, tip.regex, outgroup.name, raw.blacklist.threshold, ratio.blacklist.threshold, parsimony.blacklist.k, has.read.counts, read.counts.matter.on.zero.length.tips, verbose), simplify = F, USE.NAMES = T)
  }
  
  # 11. Dual blacklisting
  
  if(do.dual.blacklisting){
    hosts.that.are.duals <- lapply(all.tree.info, function(tree.info){
      tree.info$duals.info$host
    })
    hosts.that.are.duals <- unique(unlist(hosts.that.are.duals))
    hosts.that.are.duals <- hosts.that.are.duals[order(hosts.that.are.duals)]
    
    dual.results <- blacklist.duals(all.tree.info, hosts.that.are.duals, summary.file = NULL, verbose)
    
    all.tree.info <- sapply(all.tree.info, function(tree.info) blacklist.from.duals.list(tree.info, dual.results, verbose), simplify = F, USE.NAMES = T)
  }
  
  
  # 12. Downsampling
  
  if(downsample){
    all.tree.info <- sapply(all.tree.info, function(tree.info) blacklist.from.random.downsample(tree.info, max.reads.per.host, blacklist.underrepresented, has.read.counts, tip.regex, seed, verbose), simplify = F, USE.NAMES = T)
  }
  
  sapply(all.tree.info, function(tree.info){
 
    if(!is.null(tree.info$alignment.file.name)){
    
      seqs     <- read.dna(tree.info$alignment.file.name, format="fasta")
      new.seqs <- seqs[which(!(labels(seqs) %in% tree.info$original.tip.labels[tree.info$blacklist])),]
      new.afn  <- paste0(output.file.id, tree.info$suffix, ".fasta")
      
      if(verbose) cat("Writing cleaned alignment to ",new.afn, "\n", sep="")
      write.dna(new.seqs, new.afn, format="fasta")
    } else {
      if(verbose) cat("No alignment file found for tree ID ",tree.info$suffix, "\n", sep="")
    }
    
  })
  
  TRUE
}


#' Make a \code{data.table} of per-window host statistics
#'
#' This function collects per-window statistics on hosts
#' @param phyloscanner.trees A list of class \code{phyloscanner.trees}
#' @param hosts A list of hosts to record statistics for. If not specified, every identifiable host in \code{phyloscanner.trees}
#' @param tip.regex Regular expression identifying tips from the dataset. This expects up to three capture groups, for host ID, read ID, and read count (in that order). If the latter two groups are missing then read information will not be used. The default matches input from the phyloscanner pipeline where the host ID is the BAM file name.
#' @param verbose Produce verbose output
#' @return A \code{data.table}
#' @importFrom ape drop.tip unroot
#' @importFrom phytools nodeheight
#' @importFrom data.table rbindlist
#' @export gather.summary.statistics

gather.summary.statistics <- function(phyloscanner.trees, hosts = all.hosts.from.trees(phyloscanner.trees), tip.regex = "^(.*)_read_([0-9]+)_count_([0-9]+)$", verbose = F){
  
  has.read.counts <- attr(phyloscanner.trees, 'has.read.counts')
  
  pat.stats <- lapply(phyloscanner.trees, function(x) calc.all.stats.in.window(x, hosts, tip.regex, has.read.counts, verbose))
  pat.stats <- rbindlist(pat.stats)
  
  read.proportions <- lapply(phyloscanner.trees, function(y) sapply(hosts, function(x) get.read.proportions(x, y$suffix, y$splits.table), simplify = F, USE.NAMES = T))
  
  # Get the max split count over every window and host (the exact number of columns depends on this)
  
  max.splits <- max(sapply(read.proportions, function(x) max(sapply(x, function(y) length(y)))))
  
  read.prop.columns <- lapply(phyloscanner.trees, function(x){
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
  
  pat.stats
}



#' Graph summary statistics for a single host
#' @param phyloscanner.trees A list of class \code{phyloscanner.trees}
#' @param sum.stats The output of a call to \code{gather.summary.statistics}.
#' @param host The host to obtain graphs for.
#' @param verbose Verbose output
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import gtable
#' @import RColorBrewer
#' @importFrom data.table melt
#' @importFrom scales pretty_breaks
#' @export draw.summary.statistics

draw.summary.statistics <- function(phyloscanner.trees, sum.stats, host, verbose = F){
  
  readable.coords <- attr(phyloscanner.trees, 'readable.coords')
  
  if(length(unique(sum.stats$file.suffix))==1){
    stop("Only one tree, cannot draw summary statistics graph.")
  }
  
  if(!(host %in% sum.stats$id)){
    stop("Cannot find this host in this summary statistics table.")
  }
  
  coordinates <- lapply(phyloscanner.trees, "[[" , "window.coords")
  
  if(readable.coords){
    coordinates <- lapply(phyloscanner.trees, "[[" , "window.coords")
    starts <- sapply(coordinates, "[[", "start")
    ends <- sapply(coordinates, "[[", "end")
    ews <- min(starts)
    lwe <- max(ends)
    
  } else {
    coordinates <- sapply(phyloscanner.trees, "[[" , "xcoord")
    range <- max(coordinates) - min(coordinates)
    increment <- range/length(coordinates)
    
    ews <- min(coordinates) - 0.45*increment
    lwe <- max(coordinates) + 0.45*increment
  }
  
  x.limits <- c(ews, lwe)
  
  # Set up the boundaries of each window's region on the x-axis
  
  xcoords <- unique(sum.stats$xcoord)
  xcoords <- xcoords[order(xcoords)]
  
  missing.window.data <- find.gaps(xcoords)
  
  xcoords <- missing.window.data$x.coordinates
  regular.gaps <- missing.window.data$regular.gaps
  rectangles.for.missing.windows <- missing.window.data$rectangles.for.missing.windows
  bar.width <- missing.window.data$width
  
  produce.host.graphs(sum.stats, host, xcoords, x.limits, rectangles.for.missing.windows, bar.width, regular.gaps, readable.coords = readable.coords, verbose = verbose)
  
}

#' Draw summary statistics to file for many hosts as a multipage file
#' @param phyloscanner.trees A list of class \code{phyloscanner.trees}
#' @param sum.stats The output of a call to \code{gather.summary.statistics}.
#' @param hosts A vector of hosts to obtain graphs for. By default, all hosts detected in \code{phyloscanner.trees}.
#' @param file.name Output file name (should have a .pdf file extension)
#' @param height The height of each page of the output file in inches (defaults to A4 size)
#' @param width The width of each page of the output file in inches (defaults to A4 size)
#' @param verbose Verbose output
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import gtable
#' @import RColorBrewer
#' @importFrom data.table melt
#' @importFrom scales pretty_breaks
#' @export multipage.summary.statistics
multipage.summary.statistics <- function(phyloscanner.trees, sum.stats, hosts = all.hosts.from.trees(phyloscanner.trees), file.name, height=11.6929, width=8.26772, verbose = F){
  readable.coords <- attr(phyloscanner.trees, 'readable.coords')
  
  if(length(unique(sum.stats$file.suffix))==1){
    stop("Only one tree, cannot draw summary statistics graph.")
  }
  
  coordinates <- lapply(phyloscanner.trees, "[[" , "window.coords")
  
  if(readable.coords){
    coordinates <- lapply(phyloscanner.trees, "[[" , "window.coords")
    starts <- sapply(coordinates, "[[", "start")
    ends <- sapply(coordinates, "[[", "end")
    ews <- min(starts)
    lwe <- max(ends)
    
  } else {
    coordinates <- sapply(phyloscanner.trees, "[[" , "xcoord")
    range <- max(coordinates) - min(coordinates)
    increment <- range/length(coordinates)
    
    ews <- min(coordinates) - 0.45*increment
    lwe <- max(coordinates) + 0.45*increment
  }
  
  x.limits <- c(ews, lwe)
  
  # Set up the boundaries of each window's region on the x-axis
  
  xcoords <- unique(sum.stats$xcoord)
  xcoords <- xcoords[order(xcoords)]
  
  missing.window.data <- find.gaps(xcoords)
  
  xcoords <- missing.window.data$x.coordinates
  regular.gaps <- missing.window.data$regular.gaps
  rectangles.for.missing.windows <- missing.window.data$rectangles.for.missing.windows
  bar.width <- missing.window.data$width
  
  produce.pdf.graphs(file.name, sum.stats, hosts, xcoords, x.limits, rectangles.for.missing.windows, bar.width, regular.gaps = F, width, height, readable.coords = F, verbose = F)
}

#' Write the phylogeny with reconstructed host annotations to file
#' @param phyloscanner.tree A list of class \code{phyloscanner.tree} (usually an item in a list of class \code{phyloscanner.trees})
#' @param file.name The name of the output file
#' @param format The format - PDF or NEXUS - in which to write the output.
#' @param pdf.scale.bar.width The width, in substitutions per site, of the scale bar in PDF output
#' @param pdf.w The width of the output PDF file, in inches
#' @param pdf.hm The height, in inches per tip, of the output PDF file
#' @importFrom ggtree ggtree geom_point2 geom_tiplab geom_treescale
#' @import ggplot2
#' @import ggtree
#' @export write.annotated.tree
write.annotated.tree <- function(phyloscanner.tree, file.name, format = c("pdf", "nexus"), pdf.scale.bar.width = 0.01, pdf.w = 50, pdf.hm = 0.15){
  tree <- phyloscanner.tree$tree
  
  if(format == "pdf"){
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
    ggsave(file.name, device="pdf", height = pdf.hm*length(tree$tip.label), width = pdf.w, limitsize = F)
  } else if(format == "nexus"){
    write.ann.nexus(tree, file=file.name, annotations = c("INDIVIDUAL", "SPLIT"))
  }
}

#' Summarise the pairwise host relationships across all trees
#' @param phyloscanner.trees A list of class \code{phyloscanner.trees}
#' @param win.threshold The proportion of windows that a pair of hosts need to be related (adjacent and within \code{dist.threshold} of each other) in order for them to appear in the summary.
#' @param dist.threshold The patristic distance within which the subgraphs from two hosts need to be in order for them to be declared related (default is infinity, so adjacent hosts are always related).
#' @param allow.mt If FALSE, directionality is only inferred between pairs of hosts where a single clade from one host is nested in one from the other; this is more conservative.
#' @param close.sib.only If TRUE, then the distance threshold applies only to hosts on sibiling clades. Any ancestry is automatically a relationship.
#' @param verbose Give verbose output
#' @return A \code{data.table}, every line of which counts the number of pairwise relationships of a particular type between a pair of hosts
#' @importFrom data.table copy setkey setcolorder
#' @export transmission.summary

transmission.summary <- function(phyloscanner.trees, win.threshold=0, dist.threshold=Inf, allow.mt=T, close.sib.only = F, verbose = F){
  if(length(phyloscanner.trees)==1){
    stop("Can't summarise transmission information on a single tree. Use the collapsed tree instead?")
  }
  summarise.classifications(phyloscanner.trees, win.threshold*length(phyloscanner.trees), dist.threshold, allow.mt, close.sib.only, verbose, F)
}

#' Simplfy and visually display the pairwise host relationships across all trees
#' @param phyloscanner.trees A list of class \code{phyloscanner.trees}
#' @param trans.summary The output of \code{transmission.summary}; a \code{data.table}.
#' @param arrow.threshold The proportion of trees in which a pair of hosts need to show a direction of transmission for that direction to be indicated as an arrow. If both directions meet this threshold, the arrow is in the direction with the larger proportion of trees.
#' @param plot If TRUE, the returned list has an item called \code{simp.diagram}, a \code{ggplot} object plotting the simplified relationship diagram.
#' @importFrom GGally ggnet2
#' @export simplified.transmission.summary

simplified.transmission.summary <- function(phyloscanner.trees, transmission.summary, arrow.threshold, plot = F){
  total.trees <- length(phyloscanner.trees)
  
  transmission.summary$ancestry <- as.character(transmission.summary$ancestry)
  
  simplify.summary(transmission.summary, arrow.threshold, total.trees, plot)
}

#' @export
#' @keywords internal
identify.window.coords <- function(tree.info, file.name.regex) {
  tryCatch({
    coords                  <- get.window.coords(tree.info$suffix, file.name.regex)
    tree.info$window.coords <- coords
    tree.info$xcoord        <- (coords$end + coords$start)/2
    tree.info
  }, error = function(e){
    tree.info$xcoord        <- which(names(all.tree.info) == tree.info$suffix)
    tree.info
  })
}

#' @export
#' @keywords internal

attach.file.names <- function(tree.info, paths, identifiers, item.name){

  expected.file.name  <- paths[which(identifiers==tree.info$suffix)]
  if(length(expected.file.name)!=0){
    tree.info[[item.name]] <- expected.file.name
  }

  tree.info
}

# attach.file.names <- function(tree.info, 
#                               full.user.blacklist.file.names = NULL, 
#                               full.recombination.file.names = NULL, 
#                               full.duplicate.file.names = NULL, 
#                               full.alignment.file.names = NULL) {
#   if(!is.null(full.user.blacklist.file.names)){
#     expected.user.blacklist.file.name  <- full.user.blacklist.file.names[which(user.blacklist.identifiers==tree.info$suffix)]
#     if(length(expected.user.blacklist.file.name)!=0){
#       tree.info$user.blacklist.file.name <- expected.user.blacklist.file.name
#     }
#   }
#   
#   if(!is.null(full.recombination.file.names)){
#     expected.recomb.file.name  <- full.recombination.file.names[which(recomb.identifiers==tree.info$suffix)]
#     if(length(expected.recomb.file.name)!=0){
#       tree.info$recombination.file.name <- expected.recomb.file.name
#     }
#   }
#   
#   if(!is.null(full.duplicate.file.names)){
#     expected.duplicate.file.name  <- full.duplicate.file.names[which(duplicate.identifiers==tree.info$suffix)]
#     if(length(expected.duplicate.file.name)!=0){
#       tree.info$duplicate.file.name <- expected.duplicate.file.name
#     }
#   }
#   
#   if(!is.null(full.alignment.file.names)){
#     expected.alignment.file.name  <- full.alignment.file.names[which(alignment.identifiers==tree.info$suffix)]
#     if(length(expected.alignment.file.name)!=0){
#       tree.info$alignment.file.name <- expected.alignment.file.name
#     }
#   }
#   
#   tree.info
# }

#' @export
#' @keywords internal

attach.tree <- function(tree.info, verbose) {
  if(verbose){
    cat("Reading tree file",tree.info$tree.file.name,'\n')
  }
  
  first.line                     <- readLines(tree.info$tree.file.name, n=1)
  
  if(first.line == "#NEXUS"){
    tree                         <- read.nexus(tree.info$tree.file.name)
  } else {
    tree                         <- read.tree(tree.info$tree.file.name)
  }
  
  tree.info$tree                 <- tree
  
  tree.info
}

#' @export
#' @keywords internal

check.read.counts <- function(tree.info, tip.regex){
  tip.labels   <- tree.info$tree$tip.label
  read.counts  <- sapply(tip.labels, read.count.from.label, tip.regex)
  return(!all(is.na(read.counts)))
}

#' @export
#' @keywords internal

prepare.tree <- function(tree.info, outgroup.name, tip.regex, guess.multifurcation.threshold, multifurcation.threshold = -1, verbose = F){
  if(verbose){
    cat("Processing tree ID ",tree.info$suffix,"...\n", sep="")
  }
  
  tree <- tree.info$tree
  
  if(guess.multifurcation.threshold){
    minimum.bl                        <- min(tree$edge.length)
    if(readable.coords){
      window.width = tree.info$window.coords$end - tree.info$window.coords$start + 1
      one.snp <- 1/window.width
      
      if(minimum.bl > 0.25*one.snp){
        if(verbose){
          cat("In tree ID ",tree.info$suffix," the minimum branch length is ",minimum.bl,", which is equivalent to ",minimum.bl/one.snp," SNPs. Assuming this tree has no multifurcations.\n", sep="")
        }
        multifurcation.threshold                      <- -1
      } else {
        if(verbose) {
          cat("In tree ID ",tree.info$suffix," the minimum branch length is ",minimum.bl,", which is equivalent to ",minimum.bl/one.snp," SNPs. Using this branch length as a multifurcation threshold.\n", sep="")
        }
        if(minimum.bl==0){
          multifurcation.threshold                    <- 1E-9
        } else {
          multifurcation.threshold                    <- minimum.bl*1.0001
        }
      }
      
    } else {
      warning("Attempting to guess a branch length threshold for multifurcations from the tree. Please ensure that the tree has multifurcations before using the results of this analysis.")
      if(verbose){
        cat("In tree ID ",tree.info$suffix," the minimum branch length is ",minimum.bl,". Using this as a multifurcation threshold.\n", sep="")
      } 
      if(minimum.bl==0){
        multifurcation.threshold                    <- 1E-9
      } else {
        multifurcation.threshold                    <- minimum.bl*1.0001
      }
    }
  }
  
  new.tree <- process.tree(tree, outgroup.name, multifurcation.threshold)
  
  tree.info$tree                      <- new.tree
  tree.info$original.tip.labels       <- new.tree$tip.label
  
  hosts.for.tips                      <- sapply(tree$tip.label, function(x) host.from.label(x, tip.regex))
  tree.info$hosts.for.tips            <- hosts.for.tips
  
  tree.info
}

#' @export
#' @keywords internal

read.blacklist <- function(tree.info, verbose = F) {
  if(!is.null(tree.info$user.blacklist.file.name)){
    if(file.exists(tree.info$user.blacklist.file.name)){
      if (verbose) cat("Reading blacklist file ",tree.info$user.blacklist.file.name,'\n',sep="")
      blacklisted.tips                    <- read.table(tree.info$user.blacklist.file.name, sep=",", header=F, stringsAsFactors = F, col.names="read")
      blacklist                           <- vector()
      if(nrow(blacklisted.tips)>0){
        blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree.info$tree))
      }
      if(any(is.na(blacklist))){
        warning("Some tips listed in blacklist file ",tree.info$user.blacklist.file.name," are not tips of tree ",tree.info$tree.file.name, sep="")
      }
      blacklist <- blacklist[!is.na(blacklist)]
      
      tree.info$hosts.for.tips[blacklist] <- NA 
      
      if(verbose & length(blacklist)>0) {
        cat(length(blacklist), " tips pre-blacklisted for tree ID ",tree.info$suffix, ".\n", sep="")
      }
      
      tree.info$blacklist                 <- blacklist
    } else {
      cat(paste("WARNING: File ",tree.info$blacklist.input," does not exist; skipping.\n",sep=""))
    }
  } 
  
  tree.info
}

#' @export
#' @keywords internal

rename.user.blacklist.tips <- function(tree.info) {
  
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
}

#' @export
#' @keywords internal

apply.normalisation.constants <- function(tree.info) {
  tree              <- tree.info$tree
  tree$edge.length  <- tree$edge.length/tree.info$normalisation.constant
  tree.info$tree    <- tree
  tree.info
}

#' @export
#' @keywords internal

find.duplicate.tips <- function(tree.info) {
  if(file.exists(tree.info$duplicate.file.name)){
    tree.info$duplicate.tips <- strsplit(readLines(tree.info$duplicate.file.name, warn=F),",")
  } else {
    warning("No duplicates file found for tree suffix ",tree.info$suffix, "; skipping duplicate blacklisting.")
  }
  tree.info
}

#' @export
#' @keywords internal

blacklist.from.duplicates.vector <- function(tree.info, verbose) {
  tree <- tree.info$tree
  
  if(!is.null(tree.info$duplicate.tips)){
    
    duplicated                                   <- blacklist.exact.duplicates(tree.info, raw.blacklist.threshold, ratio.blacklist.threshold, tip.regex, verbose)
    
    duplicate.nos                                <- which(tree.info$original.tip.labels %in% duplicated)
    
    newly.blacklisted                            <- setdiff(duplicate.nos, tree.info$blacklist) 
    
    if(verbose & length(newly.blacklisted > 0)) cat(length(newly.blacklisted), " tips blacklisted as duplicates for tree ID ",tree.info$suffix, "\n", sep="")
    
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
}

#' @export
#' @keywords internal

blacklist.using.parsimony <- function(tree.info, tip.regex, outgroup.name, raw.blacklist.threshold, ratio.blacklist.threshold, parsimony.blacklist.k, has.read.counts, read.counts.matter, verbose){
  
  tree <- tree.info$tree
  tip.hosts <- sapply(tree$tip.label, function(x) host.from.label(x, tip.regex))
  tip.hosts[tree.info$blacklist] <- NA
  
  hosts <- unique(na.omit(tip.hosts))
  
  hosts <- hosts[order(hosts)]
  
  results <- sapply(hosts, function(x) get.splits.for.host(x, tip.hosts, tree, outgroup.name, raw.blacklist.threshold, ratio.blacklist.threshold, "s", parsimony.blacklist.k, 0, read.counts.matter, T, !has.read.counts, tip.regex, verbose), simplify = F, USE.NAMES = T)
  
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
}

#' @export
#' @keywords internal

blacklist.from.duals.list <- function(tree.info, dual.results, verbose) {
  tree <- tree.info$tree
  
  if(!is.null(dual.results[[tree.info$suffix]])){
    dual                                        <- dual.results[[tree.info$suffix]]
    dual.nos                                    <- which(tree.info$original.tip.labels %in% dual)
    
    newly.blacklisted                           <- setdiff(dual.nos, tree.info$blacklist) 
    
    if(verbose & length(newly.blacklisted)>0) cat(length(newly.blacklisted), " tips blacklisted for belonging to minor subgraphs in tree ID ",tree.info$suffix, "\n", sep="")
    
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
}

#' @export
#' @keywords internal

blacklist.from.random.downsample <- function(tree.info, max.reads.per.host, blacklist.underrepresented, has.read.counts, tip.regex, seed, verbose){
  tree.info <- downsample.tree(tree.info, NULL, max.reads.per.host, T, blacklist.underrepresented, !has.read.counts, tip.regex, seed, verbose)
  tree.info$hosts.for.tips[tree.info$blacklist] <- NA 
  
  tree.info
}