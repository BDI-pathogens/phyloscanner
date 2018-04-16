
#' @export
#' @keywords internal
initialise.phyloscanner <- function(
  tree.file.directory,
  tree.file.regex = "^RAxML_bestTree.InWindow_([0-9]+_to_[0-9]+)\\.tree$",
  outgroup.name = NULL,
  multifurcation.threshold = -1,
  guess.multifurcation.threshold = F,
  user.blacklist.directory = NULL,
  user.blacklist.file.regex = NULL,
  duplicate.file.directory = NULL,
  duplicate.file.regex = "^DuplicateReadCountsProcessed_InWindow_([0-9]+_to_[0-9]+).csv$",
  recombination.file.directory = NULL,
  recombination.file.regex = "^RecombinantReads_InWindow_([0-9]+_to_[0-9]+).csv$",
  alignment.file.directory = NULL,
  alignment.file.regex = NULL,
  tip.regex = "^(.*)_read_([0-9]+)_count_([0-9]+)$",
  file.name.regex = "^\\D*([0-9]+)_to_([0-9]+)\\D*$",
  norm.ref.file.name = NULL,
  norm.standardise.gag.pol = F,
  norm.constants = NULL,
  verbosity = 0){
  
  if(verbosity!=0){
    cat("Initialising...\n")
  }
  
  full.tree.file.names <- list.files.mod(tree.file.directory, pattern=tree.file.regex, full.names=TRUE)
  tree.file.names <- list.files.mod(tree.file.directory, pattern=tree.file.regex)
  
  if(length(tree.file.names)==0){
    stop("No tree files found.")
  }

  if(!is.null(norm.ref.file.name) & !is.null(norm.constants)){
    stop("Please either ask for calculation of normalisation constants, provide your own, or neither.")
  }
  
  do.dup.blacklisting   <- !is.null(duplicate.file.directory)
  include.alignment     <- !is.null(alignment.file.directory)
  existing.bl           <- !is.null(user.blacklist.directory)
  do.recomb             <- !is.null(recombination.file.directory)
  
  # Make the big list.
  
  ptrees <- list()
  
  # There are two ways to match files - to IDs and to window coordinates.
  
  match.mode <- NA
  
  if(!is.null(tree.file.regex)){
    tree.identifiers <- sapply(tree.file.names, function(x) sub(tree.file.regex, "\\1", x))
    if(all(tree.identifiers!="")){
      match.mode <- "ID"
    }
  }
  
  if(is.na(match.mode)){
    match.mode <- "coords"
    tree.identifiers <- tryCatch({
      sapply(tree.file.names, function(x) get.window.coords.string(x, file.name.regex))},
      error = function(e){
        match.mode <- "none"
        tree.file.names
      })
  }
  
  if(length(unique(tree.identifiers)) != length(tree.identifiers)){
    stop("Some trees have duplicate IDs.")
  }
  
  for(id.no in 1:length(tree.identifiers)){
    id                        <- tree.identifiers[id.no]
    
    ptree                     <- list()
    class(ptree)              <- append(class(ptree), "phyloscanner.tree")
    
    ptree$id                  <- id
    ptree$tree.file.name      <- full.tree.file.names[id.no]
    ptree$index               <- id.no
    
    ptree$bl.report           <- data.frame(tip = character(), reason = character(), row.names = NULL)
    
    ptrees[[id]]              <- ptree
  }
  
  # Coordinates.
  
  if(verbosity!=0){
    cat("Checking for window cooardinates in file names...")
  }
  
  ptrees <- sapply(ptrees, function(ptree) identify.window.coords(ptree, file.name.regex), simplify = F, USE.NAMES = T)
  
  readable.coords <- all(sapply(ptrees, function(x) !is.null(x$window.coords)))
  
  if(verbosity!=0){
    if(readable.coords){
      cat("found\n")
    } else {
      cat("not found\n")
    }
  }
  
  # Attach blacklist and recombination files
  
  if(existing.bl){
    
    user.blacklist.file.names <- list.files.mod(user.blacklist.directory, pattern=user.blacklist.file.regex)
    full.user.blacklist.file.names <- list.files.mod(user.blacklist.directory, pattern=user.blacklist.file.regex, full.names=TRUE)
    if(length(tree.file.names == 1) & length(user.blacklist.file.names)==1){
      user.blacklist.identifiers <- tree.identifiers
    } else {
      if(match.mode == "ID"){
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
        warning("Tree files and blacklist files do not entirely match. Blacklist files with no tree file will be ignored; tree files with no blacklist file will have no tips blacklisted.")
      }
    }
    
    ptrees <- sapply(ptrees, function(ptree) attach.file.names(ptree, full.user.blacklist.file.names, user.blacklist.identifiers, "user.blacklist.file.name"), simplify = F, USE.NAMES = T)
  }
  
  if(do.recomb){
    recombination.file.names <- list.files.mod(recombination.file.directory, pattern=recombination.file.regex)
    full.recombination.file.names <- list.files.mod(recombination.file.directory, pattern=recombination.file.regex, full.names=TRUE)
    if(length(tree.file.names == 1) & length(recombination.file.names)==1){
      recomb.identifiers <- tree.identifiers
    } else {
      if(match.mode=="ID"){
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
    
    ptrees <- sapply(ptrees, function(ptree) attach.file.names(ptree, full.recombination.file.names, recomb.identifiers, "recombination.file.name"), simplify = F, USE.NAMES = T)
  }
  
  if(do.dup.blacklisting){
    duplicate.file.names <- list.files.mod(duplicate.file.directory, pattern=duplicate.file.regex)
    full.duplicate.file.names <- list.files.mod(duplicate.file.directory, pattern=duplicate.file.regex, full.names=TRUE)
    if(length(tree.file.names == 1) & length(duplicate.file.names)==1){
      duplicate.identifiers <- tree.identifiers
    } else {
      if(match.mode=="ID"){
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
      if(length(intersect(tree.identifiers, duplicate.identifiers))==0){
        stop("Tree files identifiers and duplicate file identifiers do not match at all. Check file prefixes are correct.")
      }
      if(!setequal(tree.identifiers, duplicate.identifiers)){
        warning("Tree files identifiers and duplicate file identifiers do not entirely match.")
      }
    }
    
    ptrees <- sapply(ptrees, function(ptree) attach.file.names(ptree, full.duplicate.file.names, duplicate.identifiers, "duplicate.file.name"), simplify = F, USE.NAMES = T)
  }
  
  if(include.alignment){
    alignment.file.names <- list.files.mod(alignment.file.directory, pattern=alignment.file.regex)
    full.alignment.file.names <- list.files.mod(alignment.file.directory, pattern=alignment.file.regex, full.names=TRUE)
    if(length(tree.file.names == 1) & length(alignment.file.names)==1){
      alignment.identifiers <- tree.identifiers
    } else {
      if(match.mode=="ID"){
        alignment.identifiers <- sapply(alignment.file.names, function(x) sub(alignment.file.regex, "\\1", x))
      } else {
        if(match.mode != "coords"){
          stop("Cannot match alignment files with tree files using the information given.")
        } else {
          tryCatch({
            alignment.identifiers <- sapply(alignment.file.names, function(x) get.window.coords.string(x, file.name.regex))},
            error = function(e){
              stop("Cannot match alignment files with tree files using the information given.")
            })
        }
      }
      if(length(intersect(tree.identifiers, alignment.identifiers))==0){
        stop("Tree files identifiers and alignment file identifiers do not match at all. Check file prefixes are correct.")
      }
      if(!setequal(tree.identifiers, duplicate.identifiers)){
        warning("Tree files identifiers and alignment file identifiers do not entirely match.")
      }
    }
    
    ptrees <- sapply(ptrees, function(ptree) attach.file.names(ptree, full.alignment.file.names, alignment.identifiers, "alignment.file.name"), simplify = F, USE.NAMES = T)
  }
  
  # Read the trees
  
  if(verbosity!=0){
    cat("Reading trees...\n")
  }
  
  ptrees <- sapply(ptrees, function(ptree) attach.tree(ptree, verbosity==2), simplify = F, USE.NAMES = T)
  
  # Are there read counts at the tips?
  
  if(verbosity!=0){
    cat("Checking for read counts in tip names...")
  }  
  
  read.counts.check <- sapply(ptrees, function(ptree) check.read.counts(ptree, tip.regex))
  
  if(any(read.counts.check) & !all(read.counts.check)){
    warning("Read counts are not present in some trees; ignoring them throughout.\n")
  }
  has.read.counts <- all(read.counts.check)
  
  if(verbosity!=0){
    if(has.read.counts){
      cat("found\n")
    } else {
      cat("not found\n")
    }
  }  
  
  # Root tree, collapse multifurcations, get host for each tip
  
  if(guess.multifurcation.threshold){
    warning("Attempting to guess a branch length threshold for multifurcations from the tree. Please visually examine the tree or trees for multifurcations before using the results of this analysis.")
  }
  
  if(verbosity!=0){
    cat("Preparing trees...\n")
  }
  
  ptrees <- sapply(ptrees, function(ptree) prepare.tree(ptree, outgroup.name, tip.regex, guess.multifurcation.threshold, readable.coords, multifurcation.threshold, verbosity==2),
                   simplify = F, USE.NAMES = T)
  
  ptrees[sapply(ptrees, is.null)] <- NULL
  
  if(length(ptrees)==0){
    stop("Cannot find any hosts on any tree that match this regular expression. Please check that it is correct.")
  }
  
  # Get the normalisation constants
  
  if(!is.null(norm.constants)){
    
    if(!is.na(suppressWarnings(as.numeric(norm.constants)))){
      nc <- as.numeric(norm.constants)
      
      if (verbosity!=0) cat('Normalising branch lengths using specified constant (', nc, ")...\n", sep="")
      
      ptrees <- sapply(ptrees, function(ptree) {
        ptree$normalisation.constant  <- nc
        ptree
      }, simplify = F, USE.NAMES = T)
    } else if(file.exists(norm.constants)){
      
      if (verbosity!=0) cat('Normalising branch lengths using contents of file ', norm.constants, "...\n", sep="")
      
      
      nc.df   <- read.csv(norm.constants, stringsAsFactors = F, header = F)
      ptrees <- sapply(ptrees, function(ptree) {
        rows  <- which(nc.df[,1]==basename(ptree$tree.file.name))
        
        if(length(rows)>0){
          # Take the first appearance of the file name
          
          row <- rows[1]
          nc  <- as.numeric(nc.df[row, 2])
          if(is.na(nc)){
            warning("Tree with ID ",ptree$id," given a normalisation constant in file ",norm.constants," which cannot be processed as a number; this tree will not have normalised branch lengths\n")
            nc <- 1
          }
        } else {
          warning("Tree with ID ",ptree$id," has no normalisation constant in this lookup file and will not have normalised branch lengths\n")
          nc <- 1
        }
        
        ptree$normalisation.constant  <- nc
        ptree
      }, simplify = F, USE.NAMES = T)
      
    } else {
      ptrees <- sapply(ptrees, function(ptree) {
        ptree$normalisation.constant  <- 1
        ptree
      }, simplify = F, USE.NAMES = T)
      warning(paste0("Normalisation file ",norm.constants," not found. Tree branch lengths will not be normalised.\n"))
    }
  } else if(!is.null(norm.ref.file.name)){
    if(readable.coords){
      if (verbosity!=0) cat('Calculating normalisation constants from file ', norm.ref.file.name, "...\n", sep="")
      
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
            if (verbosity==2) cat('Standardising normalising constants to 1 on the gag+pol (prot + first part of RT in total 1300bp pol) region\n')
            
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
            if (verbosity==2) cat('Standardising normalising constants to 1 on the whole genome\n')
            
            tmp		<- norm.table[, mean(NORM_CONST)]
            if(!is.finite(tmp)){
              stop(paste0("Standardising constant is not finite"))
            }
            
            set(norm.table, NULL, 'NORM_CONST', norm.table[, NORM_CONST/tmp])
          }
          
          ptrees <- sapply(ptrees, function(ptree){
            ptree$normalisation.constant <- lookup.normalisation.for.tree(ptree, norm.table, lookup.column = "NORM_CONST")
            ptree
          }, simplify = F, USE.NAMES = T)
        }
      } else {
        warning(paste0("Unknown input file format for normalisation file; tree branch lengths will not be normalised.\n"))
        ptrees <- sapply(ptrees, function(ptree) {
          ptree$normalisation.constant  <- 1
          ptree
        }, simplify = F, USE.NAMES = T)
      }
    } else {
      warning(paste0("Cannot normalise branch lengths from file without window cooardinates in file suffixes; tree branch lengths will not be normalised.\n"))
      ptrees <- sapply(ptrees, function(ptree) {
        ptree$normalisation.constant  <- 1
        ptree
      }, simplify = F, USE.NAMES = T)
    }
  } else {
    ptrees <- sapply(ptrees, function(ptree) {
      ptree$normalisation.constant  <- 1
      ptree
    }, simplify = F, USE.NAMES = T)
  }
  
  list(ptrees = ptrees, has.read.counts = has.read.counts, readable.coords = readable.coords, match.mode = match.mode)
}


#' @export
#' @keywords internal
blacklist <- function(ptrees,
                      raw.blacklist.threshold,
                      ratio.blacklist.threshold,
                      parsimony.blacklist.k,
                      has.read.counts,
                      count.reads.in.parsimony,
                      tip.regex,
                      max.reads.per.host = 0,
                      blacklist.underrepresented = 0,
                      do.dup.blacklisting,
                      do.dual.blacklisting,
                      outgroup.name,
                      verbosity){
  
  
  # Read the user blacklists
  
  if(any(sapply(ptrees, function(x) !is.null(x$user.blacklist.file.name))) & verbosity!=0){
    cat("Reading user blacklists...\n")
  }
  
  ptrees <- sapply(ptrees, function(ptree) read.blacklist(ptree, verbosity==2), simplify = F, USE.NAMES = T)
  
  # Rename tips from the prexisting blacklist
  
  ptrees <- sapply(ptrees, function(ptree) rename.user.blacklist.tips(ptree), simplify = F, USE.NAMES = T)
  
  # Duplicate blacklisting
  
  if(do.dup.blacklisting){
    if (verbosity!=0) cat("Blacklisting duplicates from input files...\n", sep="")
    
    ptrees <- sapply(ptrees, find.duplicate.tips, simplify = F, USE.NAMES = T)
    
    ptrees <- sapply(ptrees, function(ptree) blacklist.from.duplicates.vector(ptree, raw.blacklist.threshold, ratio.blacklist.threshold, tip.regex, verbosity==2), simplify = F, USE.NAMES = T)
  }
  
  # Parsimony blacklisting
  
  if(parsimony.blacklist.k > 0){
    if (verbosity!=0) cat("Blacklisting probable contaminants using parsimony...\n", sep="")
    
    ptrees <- sapply(ptrees, function(ptree) blacklist.using.parsimony(ptree, tip.regex, outgroup.name, raw.blacklist.threshold, ratio.blacklist.threshold, parsimony.blacklist.k, has.read.counts, count.reads.in.parsimony, verbosity==2), simplify = F, USE.NAMES = T)
  }
  
  # Dual blacklisting
  
  if(do.dual.blacklisting){
    
    if(parsimony.blacklist.k == 0){
      stop("Dual blacklisting must be performed subsequent to parsimony blacklisting")
    }
    
    if (verbosity!=0) cat("Blacklisting minor subgraphs in probable dual infections...\n", sep="")
    
    hosts.that.are.duals <- lapply(ptrees, function(ptree){
      ptree$duals.info$host
    })
    hosts.that.are.duals <- unique(unlist(hosts.that.are.duals))
    
    if(!is.null(hosts.that.are.duals)){
      
      hosts.that.are.duals <- hosts.that.are.duals[order(hosts.that.are.duals)]
      
      dual.results <- blacklist.duals(ptrees, hosts.that.are.duals, summary.file = NULL, verbose = verbosity==2)
      
      ptrees <- sapply(ptrees, function(ptree) blacklist.from.duals.list(ptree, dual.results, verbose = verbosity==2), simplify = F, USE.NAMES = T)
    } else {
      if(verbosity!=0) cat("No dual infections identified in any tree.\n")
    }
  }
  
  # Downsampling
  
  if(max.reads.per.host < Inf){
    if (verbosity!=0) cat("Downsampling to equalise read counts across hosts...\n", sep="")
    
    ptrees <- sapply(ptrees, function(ptree) blacklist.from.random.downsample(ptree, max.reads.per.host, blacklist.underrepresented, has.read.counts, tip.regex, NA, verbose = verbosity==2), simplify = F, USE.NAMES = T)
  }
  
  ptrees
}


#' Perform a phyloscanner analysis on a tree or set of trees
#'
#' These functions perform a parsimony reconstruction and classification of pairwise host relationships.
#' @param tree.file.name The name of a single tree file (Newick or NEXUS format).
#' @param tree.file.directory The directory containing all input trees.
#' @param tree.file.regex A regular expression identifying every file in \code{tree.file.directory} that is to be included in the analysis. The first capture group, if present, gives a unique string identifying each tree. If this is NULL then \code{phyloscanner} will attempt to open every file in \code{tree.file.directory}.
#' @param splits.rule The rules by which the sets of hosts are split into groups in order to ensure that all groups can be members of connected subgraphs without causing conflicts. Options: s=Sankoff with optional within-host diversity penalty (slow, rigorous, recommended), r=Romero-Severson (quick, less rigorous with >2 hosts), f=Sankoff with continuation costs (experimental).
#' @param sankoff.k For \code{splits.rule = s} or \code{f} only. The \emph{k} parameter in the Sankoff reconstruction, representing the within-host diversity penalty.
#' @param sankoff.unassigned.switch.threshold For \code{splits.rule = s} only. Threshold at which a lineage reconstructed as infecting a host will transition to the unassigned state, if it would be equally parsimonious to remain in that host.
#' @param continuation.unassigned.proximity.cost For \code{splits.rule = f} only. The branch length at which an node is reconstructed as unassigned if all its neighbouring nodes are a greater distance away. The default is 1000, intended to be effectively infinite, such a node will never normally receive the unassigned state.
#' @param outgroup.name The name of the tip in the phylogeny/phylogenies to be used as outgroup (if unspecified, trees will be assumed to be already rooted). This should be sufficiently distant to any sequence obtained from a host that it can be assumed that the MRCA of the entire tree was not a lineage present in any sampled individual.
#' @param multifurcation.threshold If specified, branches shorter than this in the input tree will be collapsed to form multifurcating internal nodes. This is recommended; many phylogenetics packages output binary trees with short or zero-length branches indicating multifurcations. 
#' @param guess.multifurcation.threshold Whether to guess the multifurcation threshold from the branch lengths of the trees and the width of the genomic window (if that information is available). It is recommended that trees are examined by eye to check that they do appear to have multifurcations if using this option.
#' @param user.blacklist.file.name The path of a single text file containing the user-specified list of tips to be blacklisted
#' @param user.blacklist.directory An optional path for a folder containing pre-existing blacklist files. These tips are specified by the user to be excluded from the analysis.
#' @param user.blacklist.file.regex A regular expression identifying every file in \code{user.blacklist.directory} that contains a blacklist. If a capture group is specified then its contents will uniquely identify the tree it belongs to, which must matches the IDs found by \code{tree.file.regex}. If these IDs cannot be identified then matching will be attempted using genome window coordinates.
#' @param duplicate.file.name The path of a single \code{.csv} file specifying which tree tips are from duplicate reads. Normally this is produced by \code{phyloscanner_make_trees.py}.
#' @param duplicate.file.directory An optional path for a folder containing information on duplicate reads, to be used for duplicate blacklisting. Normally this is produced by \code{phyloscanner_make_trees.py}.
#' @param duplicate.file.regex A regular expression identifying every file in \code{duplicate.file.directory} that contains a duplicates file. If a capture group is specified then its contents will uniquely identify the tree it belongs to, which must matches the IDs found by \code{tree.file.regex}. If these IDs cannot be identified then matching will be attempted using genome window coordinates.
#' @param recombination.file.directory An optional path for a folder containing results of the \code{phyloscanner_make_trees.py} recombination metric analysis.
#' @param recombination.file.regex A regular expression identifying every file in \code{recombination.file.directory} that contains a recombination file. If a capture group is specified then its contents will uniquely identify the tree it belongs to, which must matches the IDs found by \code{tree.file.regex}. If these IDs cannot be identified then matching will be attempted using genome window coordinates.
#' @param recombination.file.name The path for a single file containing the results of the \code{phyloscanner_make_trees.py} recombination metric analysis.
#' @param alignment.directory The directory containing the alignments used to construct the phylogenies.
#' @param alignment.file.regex A regular expression identifying every file in \code{alignment.directory} that is an alignment. If a capture group is specified then its contents will uniquely identify the tree it belongs to, which must matches the IDs found by \code{tree.file.regex}. If these IDs cannot be identified then matching will be attempted using genome window coordinates. 
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
#' @param count.reads.in.parsimony If TRUE, read counts on tips will be taken into account in parsimony reconstructions at the parents of zero-length terminal branches. Not applicable for the Romero-Severson-like reconstruction method.
#' @param verbosity The type of verbose output. 0=none, 1=minimal, 2=complete
#' @param no.progress.bars Hide the progress bars from verbose output.
#' @details \code{phyloscanner.analyse.tree} is for a single phylogeny and \code{phyloscanner.analyse.trees} for a collection, while \code{phyloscanner.generate.blacklist} performs the blacklisting steps only.
#' @return A list of class \code{phyloscanner.trees}. Each element of this list is itself a list of class \code{phyloscanner.tree} and corresponds to a single tree, recording details of the \code{phyloscanner} reconstruction. The \code{names} of the \code{phyloscanner.trees} object are the tree IDs, usually derived from file suffixes. A list of class \code{phyloscanner.tree} may, depending on exact circumstances, have the following items:
#' \itemize{
#' \item{\code{id}}{ The tree ID.}
#' \item{\code{tree}}{ The tree as a \code{phylo} object. This will have been rooted and have multifurcations collapsed as requested, but branch lengths are original. It may have been pruned of blacklisted tips if \code{prune.blacklist} was specified.}
#' \item{\code{alignment}}{ The alignment as a \code{DNAbin} object.}
#' \item{\code{tree.file.name}}{ The file name from which the tree was loaded.}
#' \item{\code{alignment.file.name}}{ The file name for the alignment.}
#' \item{\code{user.blacklist.file.name}}{ The file name for the user-specified blacklist.}
#' \item{\code{duplicate.file.name}}{ The file name for the list of between-host duplicate tips.}
#' \item{\code{recombination.file.name}}{ The file name for the results of the \code{phyloscanner_make_trees.py} recombination metric analysis.}
#' \item{\code{index}}{ The index of this tree in the \code{phyloscanner.trees} list.}
#' \item{\code{bl.report}}{ A \code{data.frame} outlining the blacklisted tips in this tree and the reasons they were blacklisted.}
#' \item{\code{window.coords}}{ A vector giving the start and end of the genome cooardinates of the window from which the tree was built (if the windowed approach was used).}
#' \item{\code{xcoord}}{ A single genome position to locate this tree along the genome; generally the window midpoint in the windowed approach.}
#' \item{\code{duplicate.file.name}}{ The file name used to determine between-host duplicate tips}
#' \item{\code{original.tip.labels}}{ Blacklisting may lead to the pruinig of tips from the tree or their renaming. The original tip labels read from the tree file are recorded here.}
#' \item{\code{hosts.for.tips}}{ A vector mapping each tip onto its correspoinding hosts. Blacklisted tips are given NA.}
#' \item{\code{normalisation.constant}}{ The normalisation constant for this tree. This will be 1 if no normalisation was requested.}
#' \item{\code{duplicate.tips}}{ A list whose entries are vectors of tips whose sequences are exactly alike.}
#' \item{\code{blacklist}}{ A vector of numbers for all tips blacklisted for whatever reason. If the blacklist was pruned away, this will be empty.}
#' \item{\code{dual.detection.splits}}{ A \code{data.frame} determining the multiplicity of infection for each host as determined by parsimony blacklisting.}
#' \item{\code{duals.info}}{ A \code{data.frame} describing the subgraphs that each tip belong to in the dual infection detection, prior to parsimony and dual blacklisting.}
#' \item{\code{tips.for.hosts}}{ A list giving the tips numbers corresponding to each host}
#' \item{\code{read.counts}}{ A vector giving the read counts for each tip. Blacklisted tips and the outgroup have NAs. All non-NAs will be 1 if the data has no read count.}
#' \item{\code{splits.table}}{ A data frame giving the host and subgraph containing each tip, according to the parsimony reconstruction.}
#' \item{\code{clades.by.host}} {A list of lists of tips, each determining a monophyletic clade from one host.}
#' \item{\code{clade.mrcas.by.host}}{ A list of vectors containing the MRCA nodes of those clades.}
#' \item{\code{classification.results}}{ A \code{data.frame} desribing the pairwise topological classification of each pair of hosts in the tree.}
#' }
#' A \code{phyloscanner.trees} object has the following attributes:
#' \itemize{
#' \item{\code{readable.coords}}{ TRUE if genome window coordinates could be obtained from file names.}
#' \item{\code{match.mode}}{ Either "ID" (tree IDs were identified using \code{tree.file.regex}), "coords" (tree IDs were identified from what appear to be genome window coordinates in file names) or "none" (string IDs could not be determined).}
#' \item{\code{has.read.counts}}{ TRUE if \code{phyloscanner} detected read counts in tip labels.}
#' \item{\code{outgroup.name}}{ The tip label of the outgroup.}
#' }
#' @importFrom ape read.tree read.nexus di2multi root node.depth.edgelength
#' @importFrom data.table data.table as.data.table set setnames
#' @importFrom ff ff
#' @importFrom phangorn Ancestors Descendants Children mrca.phylo getRoot
#' @export phyloscanner.analyse.trees 

phyloscanner.analyse.trees <- function(
  tree.file.directory,
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
  alignment.file.directory = NULL,
  alignment.file.regex = NULL, 
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
  count.reads.in.parsimony = T,
  verbosity = 0,
  no.progress.bars = F){
  
  splits.rule <- match.arg(splits.rule)

  init <- initialise.phyloscanner(tree.file.directory,
                                    tree.file.regex,
                                    outgroup.name,
                                    multifurcation.threshold,
                                    guess.multifurcation.threshold,
                                    user.blacklist.directory,
                                    user.blacklist.file.regex,
                                    duplicate.file.directory,
                                    duplicate.file.regex,
                                    recombination.file.directory,
                                    recombination.file.regex,
                                    alignment.file.directory,
                                    alignment.file.regex,
                                    tip.regex,
                                    file.name.regex,
                                    norm.ref.file.name,
                                    norm.standardise.gag.pol,
                                    norm.constants,
                                    verbosity)
  
  ptrees <- init$ptrees
  readable.coords <- init$readable.coords
  has.read.counts <- init$has.read.counts
  match.mode <- init$match.mode
  
  do.par.blacklisting   <- parsimony.blacklist.k > 0
  
  if(do.dual.blacklisting & !do.par.blacklisting){
    warning("Dual blacklisting requires parsimony blacklisting. Turning dual blacklisting off.")
    do.dual.blacklisting <- F
  }
  
  if(splits.rule == "r"){
    if ((!is.na(sankoff.k) & sankoff.k > 0) | sankoff.unassigned.switch.threshold > 0 | continuation.unassigned.proximity.cost < 1000){
      warning("Romero-Severson reconstuction has no parameters; specified values will be ignored.")
    }
  }
  
  downsample            <- max.reads.per.host < Inf
  
  ptrees <- sapply(ptrees, function(ptree) apply.normalisation.constants(ptree), simplify = F, USE.NAMES = T)
  
  ptrees <- blacklist(ptrees,
                      raw.blacklist.threshold,
                      ratio.blacklist.threshold,
                      if(do.par.blacklisting) parsimony.blacklist.k else 0,
                      has.read.counts,
                      count.reads.in.parsimony,
                      tip.regex,
                      max.reads.per.host,
                      blacklist.underrepresented,
                      do.dup.blacklisting,
                      do.dual.blacklisting,
                      outgroup.name,
                      verbosity)
  
  # All IDs can be safely gathered and mapped now, and read counts calculated. Windows with no patients should be removed.
  
  if(verbosity!=0) cat("Gathering host IDs...\n")
  
  hosts <- all.hosts.from.trees(ptrees)
  
  if(length(hosts)==1){
    warning("Only one host detected in any tree, will skip classification of topological relationships between hosts.")
  }
  
  if(verbosity!=0){
    cat("***The full list of host IDs that will be used in this analysis follows. If this is not what you intended, check your --tipRegex***\n")
    for(host in hosts){
      cat(host, "\n")
    }
  }
  
  ptrees <- sapply(ptrees, function(ptree){
    if(all(is.na(ptree$hosts.for.tips))){
      warning("For tree ID ",ptree$id," no non-blacklisted tips remain; this window will be removed from the analysis.")
      NULL
    } else {
      tips.for.hosts <- sapply(hosts, function(x){
        
        which(ptree$hosts.for.tips == x)
        
      }, simplify = F, USE.NAMES = T)
      
      ptree$tips.for.hosts <- tips.for.hosts

      ptree
    }
  }, simplify = F, USE.NAMES = T)
  
  ptrees[sapply(ptrees, is.null)] <- NULL
  
  if(length(ptrees)==0){
    stop("All hosts have been blacklisted from all trees; nothing to do.")
  }
  
  # Prune away the blacklist if so requested
  
  if(prune.blacklist){
    if(verbosity!=0) cat("Pruning away all blacklisted tips (except the outgroup)...\n")
    
    ptrees <- sapply(ptrees, function(ptree) {
      tree                <- ptree$tree
      
      new.tree            <- process.tree(tree, outgroup.name, m.thresh = -1, ptree$blacklist)
      
      ptree$tree      <- new.tree
      ptree$blacklist <- vector()
      
      ptree
      
    }, simplify = F, USE.NAMES = T)
  }
  
  ptrees <- sapply(ptrees, function(ptree){
    if(has.read.counts){
      read.counts <- sapply(ptree$tree$tip.label, function(x) as.numeric(read.count.from.label(x, tip.regex)))
    } else {
      read.counts <- rep(1, length(ptree$tree$tip.label))
    }
    read.counts[ptree$blacklist] <- NA
    
    ptree$read.counts <- read.counts
    
    ptree
  }, simplify = F, USE.NAMES = T)
  
  # Parsimony reconstruction
  
  sankoff.p <- NA
  
  if(splits.rule == "s"){
    sankoff.p <- sankoff.unassigned.switch.threshold
  } else if(splits.rule == "f"){
    sankoff.p <- continuation.unassigned.proximity.cost
  }
  
  if(verbosity!=0) cat("Performing the parsimony reconstruction...\n", sep="")
  
  ptrees <- sapply(ptrees, function(ptree) {
    # Do the reconstruction

    if(verbosity==2) cat("Reconstructing internal node hosts on tree ID ",ptree$id, "\n", sep="")
    
    tmp					     <- split.hosts.to.subgraphs(ptree$tree, ptree$blacklist, splits.rule, tip.regex, sankoff.k, sankoff.p, use.ff, count.reads.in.parsimony, hosts, ptree$id, verbose = verbosity==2, no.progress.bars)
    tree					   <- tmp[['tree']]
    
    # trees are annotated from now on
    
    rs.subgraphs		 <- tmp[['rs.subgraphs']]
    
    if(is.null(rs.subgraphs)){
      warning("No non-blacklisted tips for ID ",ptree$id,"; this window will not be analysed further.")
    }
    
    # Convert back to unnormalised branch lengths for output
    
    tree$edge.length <- tree$edge.length * ptree$normalisation.constant
    
    if(!is.null(rs.subgraphs)){
      
      if(has.read.counts){
        rs.subgraphs$reads <- sapply(rs.subgraphs$tip, function(x) as.numeric(read.count.from.label(x, tip.regex)))
      } else {
        rs.subgraphs$reads <- 1
      }
      
      ptree$splits.table  <- rs.subgraphs
    }
    
    ptree$tree <- tree
    
    ptree
    
  }, simplify = F, USE.NAMES = T)
  
  
  # For summary statistics
  
  if(verbosity!=0) cat("Finding host MRCAs for summary statistics...\n", sep="")
  
  ptrees <- sapply(ptrees, function(ptree) {
    clade.results                 <- resolveTreeIntoPatientClades(ptree$tree, hosts, tip.regex, ptree$blacklist, !has.read.counts)
    
    ptree$clades.by.host      <- clade.results$clades.by.patient
    ptree$clade.mrcas.by.host <- clade.results$clade.mrcas.by.patient
    
    ptree
  }, simplify = F, USE.NAMES = T)
  
  # Individual window classifications
  
  if(length(hosts)>1){
    if(verbosity!=0) cat("Classifying pairwise host relationships...\n", sep="")
    
    ptrees <- sapply(ptrees, function(ptree) {
      
      if(verbosity==2) cat("Classifying host relationships for tree ID ",ptree$id, ".\n", sep="")
      
      ptree$classification.results <- classify(ptree, verbosity==2, no.progress.bars)
      
      ptree
    }, simplify = F, USE.NAMES = T)
  }
  
  class(ptrees) <- append(class(ptrees), "phyloscanner.trees")
  
  attr(ptrees, 'readable.coords') <- readable.coords
  attr(ptrees, 'match.mode')      <- match.mode
  attr(ptrees, 'has.read.counts') <- has.read.counts
  
  ptrees
}

#' @export
#' @rdname phyloscanner.analyse.trees

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
  alignment.file.name = NULL,
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
  count.reads.in.parsimony = T,
  verbosity = 0,
  no.progress.bars = F){
  
  tree.file.directory <- dirname(tree.file.name)
  user.blacklist.directory <- if(!is.null(user.blacklist.file.name)) dirname(user.blacklist.file.name) else NULL
  duplicate.file.directory <- if(!is.null(duplicate.file.name)) dirname(duplicate.file.name) else NULL
  recombination.file.directory <- if(!is.null(recombination.file.name)) dirname(recombination.file.name) else NULL
  alignment.file.directory <- if(!is.null(alignment.file.name)) dirname(alignment.file.name) else NULL
  
  phyloscanner.analyse.trees(tree.file.directory, basename(tree.file.name), splits.rule, sankoff.k, sankoff.unassigned.switch.threshold,
                             continuation.unassigned.proximity.cost, outgroup.name, multifurcation.threshold, guess.multifurcation.threshold,
                             user.blacklist.directory, basename(user.blacklist.file.name), duplicate.file.directory,
                             basename(duplicate.file.name), recombination.file.directory, basename(recombination.file.name),
                             alignment.file.directory, basename(alignment.file.name), tip.regex, file.name.regex, seed, norm.ref.file.name, 
                             norm.standardise.gag.pol, norm.constants, parsimony.blacklist.k, raw.blacklist.threshold, ratio.blacklist.threshold, 
                             do.dual.blacklisting, max.reads.per.host, blacklist.underrepresented, use.ff, prune.blacklist, count.reads.in.parsimony, 
                             verbosity, no.progress.bars)
}

#' @export
#' @rdname phyloscanner.analyse.trees

phyloscanner.generate.blacklist <- function(
  tree.file.directory,
  tree.file.regex = "^RAxML_bestTree.InWindow_([0-9]+_to_[0-9]+)\\.tree$",
  outgroup.name = NULL,
  multifurcation.threshold = -1,
  guess.multifurcation.threshold = F,
  user.blacklist.directory = NULL,
  user.blacklist.file.regex = NULL,
  duplicate.file.directory = NULL,
  duplicate.file.regex = "^DuplicateReadCountsProcessed_InWindow_([0-9]+_to_[0-9]+).csv$",
  alignment.file.directory = NULL,
  alignment.file.regex = NULL, 
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
  count.reads.in.parsimony = F,
  verbosity = 0){
 
  init <- initialise.phyloscanner(tree.file.directory,
                                    tree.file.regex,
                                    outgroup.name,
                                    multifurcation.threshold,
                                    guess.multifurcation.threshold,
                                    user.blacklist.directory,
                                    user.blacklist.file.regex,
                                    duplicate.file.directory,
                                    duplicate.file.regex,
                                    NULL,
                                    NULL,
                                    alignment.file.directory,
                                    alignment.file.regex,
                                    tip.regex,
                                    file.name.regex,
                                    norm.ref.file.name,
                                    norm.standardise.gag.pol,
                                    norm.constants,
                                    verbosity)
  
  ptrees <- init$ptrees
  readable.coords <- init$readable.coords
  has.read.counts <- init$has.read.counts
  match.mode <- init$match.mode
  
  ptrees <- sapply(ptrees, function(ptree) apply.normalisation.constants(ptree), simplify = F, USE.NAMES = T)
  
  ptrees <- blacklist(ptrees,
                      raw.blacklist.threshold,
                      ratio.blacklist.threshold,
                      if(do.par.blacklisting) parsimony.blacklist.k else 0,
                      has.read.counts,
                      count.reads.in.parsimony,
                      tip.regex,
                      max.reads.per.host,
                      blacklist.underrepresented,
                      do.dup.blacklisting,
                      do.dual.blacklisting,
                      outgroup.name,
                      verbosity)
  
  class(ptrees) <- append(class(ptrees), "phyloscanner.trees")
  
  attr(ptrees, 'readable.coords') <- readable.coords
  attr(ptrees, 'match.mode')      <- match.mode
  attr(ptrees, 'has.read.counts') <- has.read.counts
  
  ptrees
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
  
  read.proportions <- lapply(phyloscanner.trees, function(y) sapply(hosts, function(x) get.read.proportions(x, y$id, y$splits.table), simplify = F, USE.NAMES = T))
  
  # Get the max split count over every window and host (the exact number of columns depends on this)
  
  max.splits <- max(sapply(read.proportions, function(x) max(sapply(x, function(y) length(y)))))
  
  read.prop.columns <- lapply(phyloscanner.trees, function(x){
    window.props <- read.proportions[[x$id]]
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
  
  if(length(unique(sum.stats$file.id))==1){
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
  
  produce.pdf.graphs(file.name, sum.stats, hosts, xcoords, x.limits, rectangles.for.missing.windows, bar.width, regular.gaps, width, height, readable.coords, verbose)
}

#' Write the phylogeny with reconstructed host annotations to file
#' @param phyloscanner.tree A list of class \code{phyloscanner.tree} (usually an item in a list of class \code{phyloscanner.trees})
#' @param file.name The name of the output file
#' @param format The format - PDF or NEXUS - in which to write the output.
#' @param pdf.scale.bar.width The width, in substitutions per site, of the scale bar in PDF output
#' @param pdf.w The width of the output PDF file, in inches
#' @param pdf.hm The height, in inches per tip, of the output PDF file
#' @param verbose Verbose output
#' @importFrom ggtree ggtree geom_point2 geom_tiplab geom_treescale
#' @import ggplot2
#' @import ggtree
#' @export write.annotated.tree
write.annotated.tree <- function(phyloscanner.tree, file.name, format = c("pdf", "nex"), pdf.scale.bar.width = 0.01, pdf.w = 50, pdf.hm = 0.15, verbose = F){
  tree <- phyloscanner.tree$tree
  read.counts <- phyloscanner.tree$read.counts
  
  attr(tree, 'BLACKLISTED') <- c(is.na(read.counts), rep(FALSE, tree$Nnode) )
  
  read.counts <- c(read.counts, rep(1, tree$Nnode))
  read.counts[which(is.na(read.counts))] <- 1
  
  attr(tree, 'READ_COUNT') <- read.counts
  
  if(verbose) cat(paste0("Writing .",format," tree to file ",file.name,"\n"))
  
  if(format == "pdf"){

    tree.display <- ggtree(tree, aes(color=BRANCH_COLOURS)) +
      geom_point2(aes(subset=SUBGRAPH_MRCA, color=INDIVIDUAL), shape = 23, size = 3, fill="white") +
      geom_point2(aes(color=INDIVIDUAL, shape=BLACKLISTED), size=1) +
      scale_fill_hue(na.value = "black", drop=F) +
      scale_color_hue(na.value = "black", drop=F) +
      scale_shape_manual(values=c(16, 4)) +
      theme(legend.position="none") +
      geom_tiplab(aes(col=INDIVIDUAL)) +
      geom_treescale(width=pdf.scale.bar.width, y=-5, offset=1.5)
    
    if(any(read.counts!=1 & !is.na(read.counts))){
      tree.display <- tree.display + geom_point2(aes(color=INDIVIDUAL, size=READ_COUNT, shape=BLACKLISTED), alpha=0.5)
    }
    
    x.max <- ggplot_build(tree.display)$layout$panel_ranges[[1]]$x.range[2]
    
    # for ggplot2 version compatibility
    
    if(is.null(x.max)){
      x.max <- ggplot_build(tree.display)$layout$panel_params[[1]]$x.range[2]
    }
    
    tree.display <- tree.display + ggplot2::xlim(0, 1.1*x.max)
    tree.display
    ggsave(file.name, device="pdf", height = pdf.hm*length(tree$tip.label), width = pdf.w, limitsize = F)
    
  } else if(format == "nex"){
    write.ann.nexus(tree, file=file.name, annotations = c("INDIVIDUAL", "SPLIT", "READ_COUNT"))
  } else {
    stop("Unknown tree output format")
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
  out <- summarise.classifications(phyloscanner.trees, win.threshold*length(phyloscanner.trees), dist.threshold, allow.mt, close.sib.only, verbose, F)
  out$ancestry <- as.character(out$ancestry)
  out
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

attach.file.names <- function(ptree, paths, identifiers, item.name){
  
  expected.file.name  <- paths[which(identifiers==ptree$id)]
  if(length(expected.file.name)!=0){
    ptree[[item.name]] <- expected.file.name
  }
  
  ptree
}

# attach.file.names <- function(ptree,
#                               full.user.blacklist.file.names = NULL,
#                               full.recombination.file.names = NULL,
#                               full.duplicate.file.names = NULL,
#                               full.alignment.file.names = NULL) {
#   if(!is.null(full.user.blacklist.file.names)){
#     expected.user.blacklist.file.name  <- full.user.blacklist.file.names[which(user.blacklist.identifiers==ptree$suffix)]
#     if(length(expected.user.blacklist.file.name)!=0){
#       ptree$user.blacklist.file.name <- expected.user.blacklist.file.name
#     }
#   }
#
#   if(!is.null(full.recombination.file.names)){
#     expected.recomb.file.name  <- full.recombination.file.names[which(recomb.identifiers==ptree$suffix)]
#     if(length(expected.recomb.file.name)!=0){
#       ptree$recombination.file.name <- expected.recomb.file.name
#     }
#   }
#
#   if(!is.null(full.duplicate.file.names)){
#     expected.duplicate.file.name  <- full.duplicate.file.names[which(duplicate.identifiers==ptree$suffix)]
#     if(length(expected.duplicate.file.name)!=0){
#       ptree$duplicate.file.name <- expected.duplicate.file.name
#     }
#   }
#
#   if(!is.null(full.alignment.file.names)){
#     expected.alignment.file.name  <- full.alignment.file.names[which(alignment.identifiers==ptree$suffix)]
#     if(length(expected.alignment.file.name)!=0){
#       ptree$alignment.file.name <- expected.alignment.file.name
#     }
#   }
#
#   ptree
# }

#' @export
#' @keywords internal

attach.tree <- function(ptree, verbose) {
  if(verbose){
    cat("Reading tree file",ptree$tree.file.name,'\n')
  }
  
  first.line                     <- readLines(ptree$tree.file.name, n=1)
  
  if(first.line == "#NEXUS"){
    tree                         <- read.nexus(ptree$tree.file.name)
  } else {
    tree                         <- read.tree(ptree$tree.file.name)
  }
  
  ptree$tree                 <- tree
  
  ptree
}

#' @export
#' @keywords internal
#' @importFrom ape read.dna

attach.alignment <- function(ptree, verbose=F, ...) {
  if(verbose){
    cat("Reading alignment file",ptree$alignment.file.name,'\n')
  }
  
  alignment                 <- read.dna(ptree$alignment.file.name, ...)
  ptree$alignment           <- alignment
  
  ptree
}


#' @export
#' @keywords internal

check.read.counts <- function(ptree, tip.regex){
  tip.labels   <- ptree$tree$tip.label
  read.counts  <- sapply(tip.labels, read.count.from.label, tip.regex)
  return(!all(is.na(read.counts)))
}

#' @export
#' @keywords internal

prepare.tree <- function(ptree, outgroup.name, tip.regex, guess.multifurcation.threshold, readable.coords = F, multifurcation.threshold = -1, verbose = F){
  
  if(verbose){
    cat("Processing tree ID ",ptree$id,"...\n", sep="")
  }
  
  tree <- ptree$tree
  
  if(guess.multifurcation.threshold){
    minimum.bl                        <- min(tree$edge.length)
    if(readable.coords){
      window.width = ptree$window.coords$end - ptree$window.coords$start + 1
      one.snp <- 1/window.width
      
      if(minimum.bl > 0.25*one.snp){
        if(verbose){
          cat("In tree ID ",ptree$id," the minimum branch length is ",minimum.bl,", which is equivalent to ",minimum.bl/one.snp," SNPs. Assuming this tree has no multifurcations.\n", sep="")
        }
        multifurcation.threshold                      <- -1
      } else {
        if(verbose) {
          cat("In tree ID ",ptree$id," the minimum branch length is ",minimum.bl,", which is equivalent to ",minimum.bl/one.snp," SNPs. Using this branch length as a multifurcation threshold.\n", sep="")
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
        cat("In tree ID ",ptree$id," the minimum branch length is ",minimum.bl,". Using this as a multifurcation threshold.\n", sep="")
      }
      if(minimum.bl==0){
        multifurcation.threshold                    <- 1E-9
      } else {
        multifurcation.threshold                    <- minimum.bl*1.0001
      }
    }
  }
  
  new.tree <- process.tree(tree, outgroup.name, multifurcation.threshold)
  
  ptree$tree                      <- new.tree
  ptree$original.tip.labels       <- new.tree$tip.label
  ptree$m.thresh                  <- multifurcation.threshold
  
  hosts.for.tips                  <- sapply(tree$tip.label, function(x) host.from.label(x, tip.regex))
  ptree$hosts.for.tips            <- hosts.for.tips
  
  ptree
}

#' @export
#' @keywords internal

read.blacklist <- function(ptree, verbose = F) {
  if(!is.null(ptree$user.blacklist.file.name)){
    if(file.exists(ptree$user.blacklist.file.name)){
      if (verbose) cat("Reading blacklist file ",ptree$user.blacklist.file.name,'\n',sep="")
      blacklisted.tips                    <- read.table(ptree$user.blacklist.file.name, sep=",", header=F, stringsAsFactors = F, col.names="read")
      blacklist                           <- vector()
      if(nrow(blacklisted.tips)>0){
        blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=ptree$tree))
      }
      if(any(is.na(blacklist))){
        warning("Some tips listed in blacklist file ",ptree$user.blacklist.file.name," are not tips of tree ",ptree$tree.file.name, sep="")
      }
      blacklist <- blacklist[!is.na(blacklist)]
      
      ptree$hosts.for.tips[blacklist] <- NA
      
      if(verbose & length(blacklist)>0) {
        cat(length(blacklist), " tips pre-blacklisted for tree ID ",ptree$id, ".\n", sep="")
      }
      
      new.rows                        <- data.frame(tip = ptree$original.tip.labels[blacklist], reason="user_specified", row.names = NULL)
      ptree$bl.report                 <- rbind(ptree$bl.report, new.rows)
      
      ptree$blacklist                 <- blacklist
    } else {
      cat(paste("WARNING: File ",ptree$blacklist.input," does not exist; skipping.\n",sep=""))
    }
  }
  
  ptree
}

#' @export
#' @keywords internal

rename.user.blacklist.tips <- function(ptree) {
  tree <- ptree$tree
  
  if(is.null(ptree$tree)){
    stop("No tree for ID ",ptree$id,"\n")
  }
  
  if(!is.null(ptree$blacklist)){
    old.tip.labels                         <- tree$tip.label
    if(length(ptree$blacklist)>0){
      new.tip.labels                       <- old.tip.labels
      new.tip.labels[ptree$blacklist]  <- paste0(new.tip.labels[ptree$blacklist], "_X_USER")
      tree$tip.label                       <- new.tip.labels
      ptree$tree                       <- tree
      
    }
  }
  ptree
}

#' @export
#' @keywords internal

apply.normalisation.constants <- function(ptree) {
  if(!("phyloscanner.tree" %in% class(ptree))){
    stop("This is not a phyloscanner.tree object")
  }
  
  tree              <- ptree$tree
  tree$edge.length  <- tree$edge.length/ptree$normalisation.constant
  ptree$tree    <- tree
  ptree
}

#' @export
#' @keywords internal

find.duplicate.tips <- function(ptree) {
  if(!is.null(ptree$duplicate.file.name)){
    if(file.exists(ptree$duplicate.file.name)){
      ptree$duplicate.tips <- strsplit(readLines(ptree$duplicate.file.name, warn=F),",")
    } else {
      warning("No duplicates file found for tree ID ", ptree$id, "; skipping duplicate blacklisting.")
    }
  } else {
    warning("No duplicates file found for tree Id ", ptree$id, "; skipping duplicate blacklisting.")
  }
  ptree
}

#' @export
#' @keywords internal

blacklist.from.duplicates.vector <- function(ptree, raw.blacklist.threshold, ratio.blacklist.threshold, tip.regex, verbose = F) {
  tree <- ptree$tree
  
  if(!is.null(ptree$duplicate.tips)){
    
    duplicated                                   <- blacklist.exact.duplicates(ptree, raw.blacklist.threshold, ratio.blacklist.threshold, tip.regex, verbose)
    
    duplicate.nos                                <- which(ptree$original.tip.labels %in% duplicated)
    
    newly.blacklisted                            <- setdiff(duplicate.nos, ptree$blacklist)
    
    if(verbose & length(newly.blacklisted > 0)) cat(length(newly.blacklisted), " tips blacklisted as duplicates for tree ID ",ptree$id, "\n", sep="")
    
    ptree$hosts.for.tips[newly.blacklisted]  <- NA
    
    old.tip.labels    <- ptree$tree$tip.label
    
    if(length(newly.blacklisted)>0){
      new.tip.labels                             <- old.tip.labels
      new.tip.labels[newly.blacklisted]          <- paste0(new.tip.labels[newly.blacklisted], "_X_DUPLICATE")
      tree$tip.label                             <- new.tip.labels
      ptree$tree                             <- tree
    }
    
    if(length(duplicate.nos)>0){
      new.rows                                   <- data.frame(tip = ptree$original.tip.labels[duplicate.nos], reason="duplicate", row.names = NULL)
      ptree$bl.report                        <- rbind(ptree$bl.report, new.rows)
    }
    
    ptree$blacklist <- unique(c(ptree$blacklist, duplicate.nos))
    ptree$blacklist <- ptree$blacklist[order(ptree$blacklist)]
  }
  ptree
}

#' @export
#' @keywords internal

blacklist.using.parsimony <- function(ptree, tip.regex, outgroup.name, raw.blacklist.threshold, ratio.blacklist.threshold, parsimony.blacklist.k, has.read.counts, count.reads.in.parsimony, verbose){
  
  tree      <- ptree$tree
  tip.hosts <- sapply(tree$tip.label, function(x) host.from.label(x, tip.regex))
  tip.hosts[ptree$blacklist] <- NA
  
  hosts   <- unique(na.omit(tip.hosts))
  
  hosts   <- hosts[order(hosts)]
  
  results <- sapply(hosts, function(x) get.splits.for.host(x, tip.hosts, tree, outgroup.name, raw.blacklist.threshold, ratio.blacklist.threshold, "s", parsimony.blacklist.k, 0, count.reads.in.parsimony, T, !has.read.counts, tip.regex, verbose), simplify = F, USE.NAMES = T)
  
  contaminant                                 <- unlist(lapply(results, "[[", 2))
  contaminant.nos                             <- which(ptree$tree$tip.label %in% contaminant)
  
  newly.blacklisted                           <- setdiff(contaminant.nos, ptree$blacklist)
  
  if(verbose & length(newly.blacklisted)>0) cat(length(newly.blacklisted), " tips blacklisted as probable contaminants by parsimony reconstruction for tree ID ",ptree$id, "\n", sep="")
  
  ptree$hosts.for.tips[newly.blacklisted] <- NA
  
  old.tip.labels                              <- ptree$tree$tip.label
  
  if(length(newly.blacklisted)>0){
    new.tip.labels                            <- old.tip.labels
    new.tip.labels[newly.blacklisted]         <- paste0(new.tip.labels[newly.blacklisted], "_X_CONTAMINANT")
    tree$tip.label                            <- new.tip.labels
    ptree$tree                                <- tree
  }
  
  if(length(contaminant.nos)>0){
    new.rows                                  <- data.frame(tip = ptree$original.tip.labels[contaminant.nos], reason="parsimony_contaminant", row.names = NULL)
    ptree$bl.report                           <- rbind(ptree$bl.report, new.rows)
  }
  
  ptree$blacklist                             <- unique(c(ptree$blacklist, contaminant.nos))
  ptree$blacklist                             <- ptree$blacklist[order(ptree$blacklist)]
  
  which.are.duals                             <- which(unlist(lapply(results, "[[", 7)))
  mi.count                                    <- unlist(lapply(results, "[[", 8))
  
  multiplicity.table                          <- data.frame(host = hosts, count = mi.count, stringsAsFactors = F, row.names = NULL)
  
  ptree$dual.detection.splits                 <- multiplicity.table
  
  if(length(which.are.duals) > 0) {
    repeat.column <- as.vector(unlist(sapply(results[which.are.duals], function (x) rep(x$id, length(x$tip.names)) )))
    mi.df                                     <- data.frame(host = repeat.column,
                                                            tip.name = unlist(lapply(results[which.are.duals], "[[", 3)),
                                                            reads.in.subtree = unlist(lapply(results[which.are.duals], "[[", 4)),
                                                            split.ids = unlist(lapply(results[which.are.duals], "[[", 5)),
                                                            tips.in.subtree = unlist(lapply(results[which.are.duals], "[[", 6)),
                                                            stringsAsFactors = F, row.names = NULL
    )
    
    ptree$duals.info <- mi.df
  }
  
  ptree
}

#' @export
#' @keywords internal

blacklist.from.duals.list <- function(ptree, dual.results, verbose) {
  tree <- ptree$tree
  
  if(!is.null(dual.results[[ptree$id]])){
    dual                                        <- dual.results[[ptree$id]]
    dual.nos                                    <- which(ptree$original.tip.labels %in% dual)
    
    newly.blacklisted                           <- setdiff(dual.nos, ptree$blacklist)
    
    if(verbose & length(newly.blacklisted)>0) cat(length(newly.blacklisted), " tips blacklisted for belonging to minor subgraphs in tree ID ",ptree$id, "\n", sep="")
    
    ptree$hosts.for.tips[newly.blacklisted] <- NA
    
    old.tip.labels                              <- ptree$tree$tip.label
    
    if(length(newly.blacklisted)>0){
      new.tip.labels                            <- old.tip.labels
      new.tip.labels[newly.blacklisted]         <- paste0(new.tip.labels[newly.blacklisted], "_X_DUAL")
      tree$tip.label                            <- new.tip.labels
      ptree$tree                                <- tree
    }
    
    if(length(dual.nos)>0){
      new.rows                                  <- data.frame(tip = ptree$original.tip.labels[dual.nos], reason="dual_infection_minor_subgraph", row.names = NULL)
      ptree$bl.report                           <- rbind(ptree$bl.report, new.rows)
    }
    
    ptree$blacklist                             <- unique(c(ptree$blacklist, dual.nos))
    ptree$blacklist                             <- ptree$blacklist[order(ptree$blacklist)]
  }
  ptree
}

#' @export
#' @keywords internal

blacklist.from.random.downsample <- function(ptree, max.reads.per.host, blacklist.underrepresented, has.read.counts, tip.regex, seed, verbose){
  ptree <- downsample.tree(ptree, NULL, max.reads.per.host, T, blacklist.underrepresented, !has.read.counts, tip.regex, seed, verbose)
  ptree$hosts.for.tips[ptree$blacklist] <- NA
  
  ptree
}

#' @export
#' @keywords internal
identify.window.coords <- function(ptree, file.name.regex) {
  tryCatch({
    coords              <- get.window.coords(ptree$id, file.name.regex)
    ptree$window.coords <- coords
    ptree$xcoord        <- (coords$end + coords$start)/2
    ptree
  }, error = function(e){
    ptree$xcoord        <- ptree$index
    ptree
  })
}
