library(phangorn)
library(optparse)
library(phytools)
library(ggplot2)

source("/Users/twoseventwo/Documents/phylotypes/TransmissionUtilityFunctions.R")

setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160307/")

file.name <- "RAxML_bestTree.InWindow_2175_to_2550.tree"
root.name <- "C.BW.00.00BW07621.AF443088"
zero.length.tips.count <- F

tree <-read.tree(file.name)
tree <- unroot(tree)
tree <- di2multi(tree, tol = 1E-5)
tree <- root(tree, outgroup = which(tree$tip.label==root.name), resolve.root = T)

tip.labels <- tree$tip.label

outgroup.no <- which(tip.labels == root.name)

split.tip.labels <- strsplit(tip.labels, label.separator)

patient.ids <- sapply(split.tip.labels, "[[", patient.id.position)

patient.ids[outgroup.no] <- "Ref"

patients <- unique(patient.ids)

patients <- patients[which(substr(patients, 1,3) == "BEE")]

patient.tips <-
  sapply(patients, get.tips.for.patient, tree = tree, patient.ids = patient.ids, blacklist=blacklist)

patient.mrcas <- lapply(patient.tips, function(node) mrca.phylo.or.unique.tip(tree, node, zero.length.tips.count))

node.assocs <- list()
node.assoc.counts <- rep(0, length(tree$tip.label) + tree$Nnode)

for(patient in patients){
  for(tip in patient.tips[[patient]]){
    keep.going <- T
    current.node <- tip
    while(keep.going){
      if(is.null(node.assocs[current.node][[1]])){
        node.assocs[current.node][1] <- patient
        node.assoc.counts[current.node] <- 1
      } else {
        if(patient %in% node.assocs[current.node][[1]]){
          break
        } else {
          node.assocs[current.node][[1]] <- c(node.assocs[current.node][[1]], patient)
          node.assoc.counts[current.node] <- node.assoc.counts[current.node] + 1
        }
      }
      if(current.node == patient.mrcas[[patient]]){
        keep.going <- F
      }
      current.node <- Ancestors(tree, current.node, type="parent")
    }
  }
}

root.no <- 3703


split.at <- function(tree, node, node.assocs, patient = NULL, careful = F){
  # NULL means split every single one
  if(!is.null(patient)){
    if(!(patient %in% node.assocs[[node]])){
      stop("Can't split based on this patient at this node")
    }
  }
  
  if(careful){
    # check you haven't failed to do this for any ancestor
    for(ancestor in Ancestors(tree, node, type="all")){
      if(!is.null(node.assocs[[node]]) & length(node.assocs[[node]])>1){
        stop("You haven't split an ancestor")
      }
    }
  }
  
  if(is.null(patient)){
    
  } else {
    
  }

}

split.patient.at <- function(tree, node, patient, patients, patient.tips, patient.mrcas, split.ref, careful = F){
  
  split.ref <<- split.ref
  
  current.split.count <- split.ref[[patient]]
  
  if(careful){
    if(patient.mrcas[[patient]] %in% Ancestors(tree, node, type="all")){
      stop("Can't split on a descendant of the MRCA")
    }
  }
  
  tip.nos <- patient.tips[[patient]]
  
  splits.to.make <- 0
  
  for(child in Children(tree, node)){
    desc.tips <- intersect(tip.nos, Descendants(tree, child, type="tips"))
    if(length(desc.tips)==length(tip.nos)){
      split.to.make <- 0
      break
    } else if(length(desc.tips>0)){
      splits.to.make <- splits.to.make + 1
    }
  }
  
  if(splits.to.make==0){
    stop("No splits here")
  }
  
  patients <- patients <- patients[patients != patient]
  split.ref[[patient]] <- split.ref[[patient]] + splits.to.make - 1
  patients <- c(patients,)
  
  for(child in Children(tree, node)){
    desc.tips <- intersect(tip.nos, Descendants(tree, child, type="tips"))
  
}