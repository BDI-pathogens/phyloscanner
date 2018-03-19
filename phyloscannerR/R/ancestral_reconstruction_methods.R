#' Reconstruct the ancestral sequence at every node of the tree
#' @param phyloscanner.tree A list of class \code{phyloscanner.tree} (usually an item in a list of class \code{phyloscanner.trees})
#' @param verbose Verbose output
#' @param ... Further arguments to be passed to optim.pml
#' @importFrom phangorn optim.pml ancestral.pml
#' @export reconstruct.ancestral.sequences

reconstruct.ancestral.sequences <- function(phyloscanner.tree, verbose, ...){
  
  # No sequences = no reconstruction
  
  if(is.null(phyloscanner.tree$alignment)){
    stop(paste0("The phyloscanner.tree object with ID ",phyloscanner.tree$suffix," has no alignment item"))
  }
  
  phyl <- phyloscanner.tree$tree
  algn <- phyloscanner.tree$alignment
  
  # The tree may have had the blacklist pruned away, or it may have had tips renamed due to blacklisting
  
  phyl.tiplabels <- phyl$tip.label
  algn.seqlabels <- names(algn)
  
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
    
    mfit <- pml(phyl, as.phyDat(algn))
    optim.mfit <- optim.pml(mfit, list(...))
    
    if(attr(phyl, "order")!="cladewise"){
      stop("Cannot map between node numbers when unrooting the tree")
    }
    
    ancestral.pml(optim.mfit)
    
  }
  
}