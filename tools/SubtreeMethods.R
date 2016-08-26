# CONSERVATIVE METHODS

# Work out which nodes must be assigned to each patient for full connectivity (ignoring conflicts)

annotate.internal <- function(tree, patients, patient.tips, patient.mrcas){
  # The patients associated with each node
  node.assocs <- list()
  # The number of patients associated with each node
  node.assoc.counts <- rep(0, length(tree$tip.label) + tree$Nnode)
  for(patient in patients){
    for(tip in patient.tips[[patient]]){
      keep.going <- T
      current.node <- tip
      while(keep.going){
        # if you're at this node, it should be given this association
        
        if(is.null(node.assocs[current.node][[1]])){
          # no existing associations
          node.assocs[current.node][1] <- patient
          node.assoc.counts[current.node] <- 1
        } else {
          if(patient %in% node.assocs[current.node][[1]]){
            # already associated, no need to go further up the tree
            break
          } else {
            # add it to the vector and increase the count
            node.assocs[current.node][[1]] <- c(node.assocs[current.node][[1]], patient)
            node.assoc.counts[current.node] <- node.assoc.counts[current.node] + 1
          }
        }
        # If the patient MRCA has been reached, can stop
        if(current.node == patient.mrcas[[patient]]){
          keep.going <- F
        }
        # Otherwise, go up
        current.node <- Ancestors(tree, current.node, type="parent")
      }
    }
  }
  return(list(details = node.assocs, counts = node.assoc.counts))
}

# Find all the splits for a given patient and update the patients vector and patient tips list

split.patient <- function(tree, patient, patients, patient.tips, associations, regexp){
#  cat("Splitting ",patient,"...\n",sep="")
  list.of.splits <- find.splits(tree, patient, patients, patient.tips, associations, regexp)
  
  new.patients <- vector()
  new.patient.tips <- list()
  
  if(length(list.of.splits) > 1){
    for(item.no in seq(1, length(list.of.splits))){
      new.name <- paste(patient,"-S",item.no, sep="")
      new.patients <- c(new.patients, new.name)
      new.patient.tips[[new.name]] <- which(tree$tip.label %in% list.of.splits[[item.no]])
    }
  }
  return(list(patients = new.patients, patient.tips = new.patient.tips))
}

find.splits <- function(tree, patient, patients, patient.tips, associations, regexp){
#  cat("Finding splits for ",patient,"...\n",sep="")
  
  tip.pool <- tree$tip.label[patient.tips[[patient]]]
  
  out <- list()
  
  last.ancestors <- sapply(tip.pool, function(x) last.nonconflicted.ancestor(tree, x, associations, regexp))
  
  for(unique.last.anc in unique(last.ancestors)){
    out[[paste("node_", unique.last.anc,sep="")]] <- tip.pool[which(last.ancestors==unique.last.anc)]
  }
  
  return(out)
  
}

# Tips in the same connected subtree will have the same earliest ancestor that is associated with
# that patient but is not conflicted

last.nonconflicted.ancestor <- function(tree, tip.label, associations, regexp){
#  cat("Finding last nonconflicted ancestor for ", tip.label, "...\n", sep="")
  tip <- which(tree$tip.label == tip.label)
  patient <- patient.from.label(tip.label, regexp)
  
  current.node <- tip
  
  repeat{
    next.node <- Ancestors(tree, current.node, type="parent")
    if(next.node == 0){
      return(current.node)
    }
    
    # associations will only contain keys up to the highest number node that has any associations
    
    if(next.node > length(associations)){
      return(current.node)
    }
    if((length(associations[[next.node]]) > 1) | !(patient %in% associations[[next.node]])){
      return(current.node)
    }
    current.node <- next.node
  }
}

# ROMERO-SEVERSON METHODS

# Annotate just the tips with their patients

annotate.tips <- function(tree, patients, patient.tips){
  node.assocs <- list()
  for(patient in patients){
    for(tip in patient.tips[[patient]]){
      node.assocs[[tip]] <- patient
    }
  }
  return(node.assocs)
}

# Used to find every connected subtree in the whole tree by finding the earliest node in them

count.splits <- function(tree, node, assocs, patients, counts.vec, first.nodes.list){
  print(node)
  
  for(child in Children(tree, node)){
    child.results <- count.splits(tree, child, assocs, patients, counts.vec, first.nodes.list) 
    
    counts.vec <- child.results$counts
    first.nodes.list <- child.results$first.nodes
  }
  
  parent <- Ancestors(tree, node, type="parent")
  change <- F
  if(parent == 0){
    change <- T
  } else if (assocs[[parent]]!=assocs[[node]]){
    change <- T
  }
  
  if(assocs[[node]]!="*"){
    if(change){
      patient <- assocs[[node]]
      pat.index <- which(patients == patient)
      counts.vec[pat.index] <- counts.vec[pat.index] + 1
      first.nodes.list[[patient]] <- c(first.nodes.list[[patient]], node)
    }
  }
  return(list(counts = counts.vec, first.nodes = first.nodes.list))
}

# Does the RS classification (todo move to TUF?)

classify.down <- function(node, tree, tip.assocs, patient.mrcas, blacklist, regex){
  
  if(is.tip(tree, node)){
    if(node %in% blacklist | is.na(patient.from.label(tree$tip.label[node], regex))){
      result <- "*"
    } else {
      result <- tip.assocs[[node]]
    }
  } else {
    
    child.assocs.vector <- sapply(Children(tree, node), classify.down, tree=tree, tip.assocs=tip.assocs, patient.mrcas = patient.mrcas, blacklist=blacklist, regex=regex)
    
    c.a.df <- as.data.frame(table(child.assocs.vector), stringsAsFactors=F )
    
    if(nrow(c.a.df)==1){
      result <- c.a.df[1,1]
    } else {
      c.a.df <- c.a.df[which(c.a.df[,1]!="*"),]
      
      winners <- c.a.df[which(c.a.df[,2]==max(c.a.df[,2])),]
      if(nrow(winners)==1){
        result <- winners[1,1]
      } else {
        result <- "*"
      }
    }
  }
  
  #check that this isn't above the MRCA for the patient, which leads to weird results
  
  if(result != "*"){
    if(node %in% Ancestors(tree, patient.mrcas[[result]], type="all")){
      result <- "*"
    }
  }
  
  new.assocs[[node]] <<- result
  return(result)
}

# "Star runs" are sequences of ancestors that are assigned "*" by the RS classification. Once the
# classification is done, they can be converted to having assigned nodes if the ancestor of the
# earliest node is a possible host for all nodes in the run

get.star.runs <- function(tree, assocs){
  out <- list()
  for(int.node in seq(length(tree$tip.label) + 1,  length(tree$tip.label) + tree$Nnode)){
    if(assocs[[int.node]]=="*"){
      #need to check that it got the star because of an ambiguity, not because there are only stars around there
      child.assocs <- unlist(assocs[Children(tree, int.node)])
      c.a.df <- as.data.frame(table(child.assocs), stringsAsFactors=F)
      c.a.df <- c.a.df[which(c.a.df[,1]!="*"),]
      if(nrow(c.a.df)>1){
        run <- int.node
        current.node <- Ancestors(tree, int.node, type="parent")
        if(current.node != 0){
          while(assocs[[current.node]]=="*"){
            run <- c(run, current.node)
            current.node <- Ancestors(tree, current.node, type="parent")
            if(current.node == 0){
              break
            }
          }
        }
        out[[int.node]] <- run
      }
    }
  }
  return(out)
}

output.trans.tree <- function(tree, assocs, file.name = NULL){
  # find the association of each node
  
  assocs.vec <- vector()
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    if(node.no > length(assocs)){
      # annoying and hacky
      assocs.vec[node.no] <- "none"
    } else if (is.null(assocs[[node.no]])){
      assocs.vec[node.no] <- "none"
    } else if (assocs[[node.no]]=="*") {
      assocs.vec[node.no] <- "none"
    } else {
      if(length(assocs[[node.no]])>1){
        stop("This method is intended for use on trees with conflicts resolved")
      }
      assocs.vec[node.no] <- assocs[[node.no]]
    }
  }
  
  first.of.split <- vector()
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    parent.no <- Ancestors(tree, node.no, type="parent")
    if(parent.no == 0){
      first.of.split[node.no] <- T
    } else {
      first.of.split[node.no] <- assocs.vec[parent.no] != assocs.vec[node.no]
    }
  }
  
  splits.vec <- assocs.vec
  
  splits.vec[which(splits.vec=="none" & first.of.split)] <- paste("none", which(splits.vec=="none" & first.of.split), sep="-S")
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    if(assocs.vec[node.no]=="none" & !first.of.split[node.no]){
      current.node.no <- node.no
      while(!first.of.split[current.node.no]){
        current.node.no <- Ancestors(tree, current.node.no, type="parent")
      }
      if(!startsWith(splits.vec[current.node.no], "none-")){
        print(node.no)
        stop("Uh-oh")
      }
      splits.vec[node.no] <- splits.vec[current.node.no]
    }
  }
  
  
  
  unique.splits <- unique(splits.vec)
  parent.splits <- vector()
  
  for(unique.split in unique.splits){
    root <- which(splits.vec == unique.split & first.of.split)
    parent.node <- Ancestors(tree, root, type="parent")
    if(parent.node != 0){
      parent.splits <- c(parent.splits, splits.vec[parent.node])
    } else {
      parent.splits <- c(parent.splits, "root")
    }
  }
  
  patients <-  unlist(lapply(strsplit(unique.splits, "-"), `[[`, 1)) 
  
  cytoscape.input <- data.frame(unique.splits, parent.splits, patients, stringsAsFactors = F)
  
  if(!is.null(file.name)){
    write.csv(cytoscape.input, file.name, row.names = F, quote=F)
  }
  
  return(cytoscape.input)
  
}

