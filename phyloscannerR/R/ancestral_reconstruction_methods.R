#' Reconstruct the ancestral sequence at every node of the tree
#' @param phyloscanner.tree A list of class \code{phyloscanner.tree} (usually an item in a list of class \code{phyloscanner.trees})
#' @param verbose Verbose output
#' @param default If TRUE, the reconstruction is done according to the default model used in RAxML to build trees for phyloscanner. The ... below will be ignored.
#' @param ... Further arguments to be passed to pml and optim.pml
#' @importFrom phangorn optim.pml ancestral.pml
#' @export reconstruct.ancestral.sequences

reconstruct.ancestral.sequences <- function(phyloscanner.tree, verbose, default=F, ...){
  
  if(default){
    pml.opts <- list(k=4, model="GTR")
    opt.pml.opts <- list(model="GTR", optGamma=TRUE, optQ=T, optEdge=F, optBf=F)
  } else {

    allowed.pml <- names(formals(pml))
    allowed.optim.pml <- names(formals(optim.pml))
    dots <- list(...)
    pml.opts <- dots[names(dots) %in% allowed.pml]
    opt.pml.opts <- dots[names(dots) %in% allowed.optim.pml]
  
  }
  
  # No sequences = no reconstruction
  
  if(is.null(phyloscanner.tree$alignment)){
    stop(paste0("The phyloscanner.tree object with ID ",phyloscanner.tree$suffix," has no alignment item"))
  }
  
  phyl <- phyloscanner.tree$tree
  algn <- phyloscanner.tree$alignment
  
  # The tree may have had the blacklist pruned away, or it may have had tips renamed due to blacklisting
  
  phyl.tiplabels <- phyl$tip.label
  
  if(is.matrix(algn)){
    algn.seqlabels <- rownames(algn)
  } else {
    algn.seqlabels <- names(algn)
  }
  
  
  if(length(intersect(phyl.tiplabels, algn.seqlabels)) != length(phyl.tiplabels)){
    # this means pruning happened or the tree is the wrong tree
    if(length(phyl.tiplabels) < length(algn.seqlabels)){
      matches <- which(phyl.tiplabels %in% algn.seqlabels)
      
      if(length(matches)==0){
        stop("No tip labels in the tree appear in the sequence names of the alignment")
      }
      
      # there are a couple of formats here
      
      if(is.matrix(algn)){
        algn <- algn[matches,]
      } else {
        algn <- algn[matches]
      }
      
    } else if (length(phyl.tiplabels) == length(algn.seqlabels)) {
      
      # the original tip labels were stored
      
      old.tiplabels <- phyloscanner.tree$original.tip.labels
      
      if(is.null(old.tiplabels) | length(intersect(old.tiplabels, algn.seqlabels)) != length(old.tiplabels)){
        stop("Tip labels in the tree and names of the alignment do not match, and the discrepancy cannot be resolved")
      }
      
      phyl$tip.label <- old.tiplabels
      
    } else {
      stop("The tree has too many tips for this alignment")
    }
  }
    
  mfit <- do.call("pml", c(list(tree=phyl, data=as.phyDat(algn)), pml.opts))
  
  # I do not want this output. I also want this to work on Windows, however. This isn't exactly elegant,ne
  
  if(file.exists("/dev/null")){
    # *nix
    sink("/dev/null")
  } else {
    # Presumably Windows...
    sink("NUL")
  }
  
  optm.mfit <- do.call("optim.pml", c(list(object=mfit), opt.pml.opts))
  
  sink()

  reconst.aln <- phyDat2alignment(ancestral.pml(optim.mfit, return = "phyDat"))

  raw <- as.list(reconst.aln$seq)
  names(raw) <- reconst.aln$nam

  # this seems faintly crazy, but it's ape's problem
  
  alignment <- as.DNAbin(as.alignment(raw))
  
  # Need to remap the names to the unrooted tree.
  
  alignment
}