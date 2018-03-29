#' Reconstruct the ancestral sequence at every node of the tree
#' @param phyloscanner.tree A list of class \code{phyloscanner.tree} (usually an item in a list of class \code{phyloscanner.trees})
#' @param verbose Verbose output
#' @param default If TRUE, the reconstruction is done according to the default model used in RAxML to build trees for phyloscanner. The \code{...} below will be ignored.
#' @param ... Further arguments to be passed to \code{pml} and \code{optim.pml}
#' @return An alignment of the sequences at all nodes (in \code{DNAbin} format)
#' @importFrom phangorn pml optim.pml ancestral.pml as.phyDat phyDat2alignment
#' @importFrom ape as.DNAbin as.alignment
#' @export reconstruct.ancestral.sequences

reconstruct.ancestral.sequences <- function(phyloscanner.tree, verbose=F, default=F, ...){
  
  if(verbose){
    cat("Reconstructing ancestral sequences on tree ID ",phyloscanner.tree$id,"\n",sep="")
  }
  
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
    stop(paste0("The phyloscanner.tree object with ID ",phyloscanner.tree$id," has no alignment item"))
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
  
  # I do not want this output. I also want this to work on Windows, however. This isn't exactly elegant.
  
  if(file.exists("/dev/null")){
    # *nix
    sink("/dev/null")
  } else {
    # Presumably Windows...
    sink("NUL")
  }
  
  optim.mfit <- do.call("optim.pml", c(list(object=mfit), opt.pml.opts))
  
  sink()

  alignment <- as.DNAbin(ancestral.pml(optim.mfit, return = "phyDat"))

  alignment

}

#' Find the ancestral sequence at the MRCA of the tips from this host, or, if a dual infection was previously identified, of the MRCA of the tips making up each infection event
#' @param phyloscanner.tree A list of class \code{phyloscanner.tree} (usually an item in a list of class \code{phyloscanner.trees}). This must have an \code{ancestral.alignment} element (see \emph{reconstruct.ancestral.sequences})
#' @param host The host ID
#' @param individual.duals Whether to output multiple sequences for \code{host} based on the results of a previous dual infection analysis
#' @param verbose Verbose output
#' @importFrom phangorn mrca.phylo
#' @export reconstruct.host.ancestral.sequences

reconstruct.host.ancestral.sequences <- function(phyloscanner.tree, host, individual.duals = F, verbose = F){
  
  if(length(phyloscanner.tree$tips.for.hosts[[host]]) == 0 ){
    if(verbose) cat(paste0("A host with ID ",host," is not present in tree ID ",phyloscanner.tree$id,"\n"))
    return(NULL)
  }
  
  if(verbose){
    cat(paste0("Finding the MRCA sequence or sequences for host ",host," in tree ",phyloscanner.tree$id,"\n"))
  }
  
  # Currently multiplicity will evaluate to >1 even if only one subgraph remains after blacklisting.
  
  multiplicity <- 1
  
  if(individual.duals){
    if(!is.null(phyloscanner.tree$dual.detection.splits)){
      multiplicity <- phyloscanner.tree$dual.detection.splits$count[phyloscanner.tree$dual.detection.splits$host==host]
    } else {
      warning("Parsimony blacklisting was not performed on these trees. With no dual infections identified, only one reconstructed sequence will be returned")
    }
  }
  
  tree <- phyloscanner.tree$tree
  
  tips <- phyloscanner.tree$tips.for.hosts[[host]]
  
  aa <- phyloscanner.tree$ancestral.alignment
  
  if(multiplicity == 1){
    
    mrca <- mrca.phylo.or.unique.tip(tree, tips)
    
    if(length(tips)==1){
      if(is.matrix(aa)){
        sequence <- aa[which(rownames(aa)==tree$tip.label[tips]),]
        rownames(sequence) <- host 
      } else {
        sequence <- aa[which(names(aa)==tree$tip.label[tips])]
        names(sequence) <- host
      }
    } else {
      if(is.matrix(aa)){
        sequence <- aa[which(rownames(aa)==as.character(mrca)),]
        rownames(sequence) <- host 
      } else {
        sequence <- aa[which(names(aa)==as.character(mrca))]
        names(sequence) <- host
      }
    }
    
    setNames(list(sequence), host)
    
  } else {
    
    if(is.null(phyloscanner.tree$duals.info)){
      stop(paste0("Expecting a duals report for tree id ",phyloscanner.tree$id," but this is absent"))
    }
    
    # we only want the non-blacklisted tips
    
    di <- phyloscanner.tree$duals.info[which(phyloscanner.tree$duals.info$tip.name %in% tree$tip.label[tips]),]
    
    out <- lapply(unique(di$split.ids), function(split){
      tips <- di$tip.name[which(di$split.ids==split)]
      
      mrca <- mrca.phylo.or.unique.tip(tree, tips)
      
      if(length(tips)==1){
        if(is.matrix(aa)){
          sequence <- aa[which(rownames(aa)==tree$tip.label[tips]),]
          rownames(sequence) <- host 
        } else {
          sequence <- aa[which(names(aa)==tree$tip.label[tips])]
          names(sequence) <- host
        }
      } else {
        if(is.matrix(aa)){
          sequence <- aa[which(rownames(aa)==as.character(mrca)),]
          rownames(sequence) <- split
        } else {
          sequence <- aa[which(names(aa)==as.character(mrca))]
          names(sequence) <- split
        }
      }
      sequence
    })
    
    setNames(out, unique(di$split.ids))
  }

}
