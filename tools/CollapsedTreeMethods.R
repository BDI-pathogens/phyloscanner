# Collapsed tree navigation

# The parent of a node in the collapsed tree

get.tt.parent <- function(tt, label){
  if(length(which(tt$unique.splits==label))>0){
    return(tt$parent.splits[which(tt$unique.splits==label)])
  } else {
    return(NULL)
  }
}

# All ancestors of a node in the collapsed tree.

get.tt.ancestors <- function(tt, label){
  out <- vector()
  current.node <- get.tt.parent(tt, label)
  while(current.node != "root"){
    out <- c(out, current.node)
    current.node <- get.tt.parent(tt, current.node)
  }
  return(out)
}

# All children of a node in the collapsed tree

get.tt.children <- function(tt, label){
  return(tt$unique.splits[which(tt$parent.splits==label)])
}

# All adjacent nodes of a node in the collapsed tree (parent and children)

get.tt.adjacent <- function(tt, label){
  return(c(get.tt.children(tt, label), get.tt.parent(tt, label)))
}

# The grandparent of a node in the collapsed tree

get.tt.parent.patient <- function(tt, label){
  if(length(which(tt$unique.splits==label))>0){
    return(tt$patients[which(tt$unique.splits==label)])
  } else {
    return(NULL)
  }
}

# MRCA of a pair of nodes in the collapsed tree

get.tt.mrca <- function(tt, label1, label2){
  # sanity check
  ca.vec <- intersect(c(label1, get.tt.ancestors(tt, label1)), c(label2, get.tt.ancestors(tt, label2)))
  if(length(ca.vec)>1){
    for(i in seq(1, length(ca.vec)-1)){
      if(get.tt.parent(tt, ca.vec[i]) != ca.vec[i+1]){
        stop("get.tt.mrca not working as intended")
      }
    }
  }
  
  return(ca.vec[1])
}

# The patient corresponding to the MRCA of two nodes

get.tt.mrca.patient <- function(tt, label1, label2){
  node <- intersect(c(label1,get.tt.ancestors(tt, label1)), c(label2, get.tt.ancestors(tt, label2)))[1]
  
  return(tt$patients[which(tt$unique.splits==node)])
}

# The path from one node in the collapsed tree to another

get.tt.path <- function(tt, label1, label2){
  mrca <- get.tt.mrca(tt, label1, label2)
  
  if(mrca == label1){
    result <- vector()
    current.node <- label2
    while(current.node != label1){
      result <- c(result, current.node)
      current.node <- get.tt.parent(tt, current.node)
    }
    result <- c(result, label1)
    return(result)
  }
  
  if(mrca == label2){
    result <- vector()
    current.node <- label1
    while(current.node != label2){
      result <- c(result, current.node)
      current.node <- get.tt.parent(tt, current.node)
    }
    result <- c(result, label2)
    return(result)
  }
  
  first.half <- vector()
  current.node <- label1
  while(current.node != mrca){
    first.half <- c(first.half, current.node)
    current.node <- get.tt.parent(tt, current.node)
  }
  
  second.half <- vector()
  current.node <- label2
  while(current.node != mrca){
    second.half <- c(second.half, current.node)
    current.node <- get.tt.parent(tt, current.node)
  }
  
  return(c(first.half, mrca, rev(second.half)))
}

# Output the collapsed tree from a phylogeny and its node assocs (write to CSV file if file.name!=NULL)
# If prune.unsampled.tips = TRUE, the output collapsed tree will not have any "unsampled_region"s that have no children
# (which usually come from blacklisted tips, and looks rather odd if e.g. read downsampling has taken place)

output.trans.tree <- function(tree, assocs, file.name = NULL, prune.unsampled.tips = TRUE){
  # find the association of each node
  
  assocs.vec <- vector()
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    if(node.no > length(assocs)){
      # annoying and hacky
      assocs.vec[node.no] <- "none"
    } else if (is.null(assocs[[node.no]])){
      assocs.vec[node.no] <- "none"
    } else if (is.na(assocs[[node.no]])){
      assocs.vec[node.no] <- "none"
    } else if (assocs[[node.no]] %in% c("*", "unsampled")) {
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
  
  unsampled.roots <- which(splits.vec=="none" & first.of.split)
  
  splits.vec[unsampled.roots] <- 
    paste("unsampled_region-", 1:length(unsampled.roots), sep="")
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    if(assocs.vec[node.no]=="none" & !first.of.split[node.no]){
      current.node.no <- node.no
      while(!first.of.split[current.node.no]){
        current.node.no <- Ancestors(tree, current.node.no, type="parent")
      }
	  if(!grepl("^unsampled_region", splits.vec[current.node.no])){      
        print(node.no)
        stop("Uh-oh")
      }
      splits.vec[node.no] <- splits.vec[current.node.no]
    }
  }
  
  unique.splits <- unique(splits.vec)
  parent.splits <- vector()
  lengths <- vector()
  root.nos <- vector()
  
  for(unique.split in unique.splits){
    root <- which(splits.vec == unique.split & first.of.split)
    parent.node <- Ancestors(tree, root, type="parent")
    if(parent.node != 0){
      parent.splits <- c(parent.splits, splits.vec[parent.node])
      
      root.edge.no <- which(tree$edge[,2]==root)
      lengths <- c(lengths, tree$edge.length[root.edge.no])
    } else {
      parent.splits <- c(parent.splits, "root")
      lengths <- c(lengths, 0)
    }
    root.nos <- c(root.nos, root)
  }
  
  patients <-  unlist(lapply(strsplit(unique.splits, "-"), `[[`, 1)) 
  parent.patients <-  unlist(lapply(strsplit(parent.splits, "-"), `[[`, 1)) 
  
  tt.table <- data.frame(unique.splits, parent.splits, patients, parent.patients, lengths, root.nos, stringsAsFactors = F)
  
  for.output <- tt.table[,1:6]
  
  if(prune.unsampled.tips){
    unsampled.tips <- which(grepl("unsampled",for.output$unique.splits) &
                              !(for.output$unique.splits %in% for.output$parent.splits))
    
    if(length(unsampled.tips) > 0){
      for.output <- for.output[-unsampled.tips,]
    }
    #renumber
    unsampled.rows <- which(grepl("^unsampled_region",for.output$unique.splits))
    
    unsampled.labels <- for.output$unique.splits[which(grepl("^unsampled_region",for.output$unique.splits))]
    
    for(x in 1:length(unsampled.rows)) {
      old.label <- unsampled.labels[x]
      new.label <- paste("unsampled_region-",x,sep="")
      for.output$unique.splits[unsampled.rows[x]] <- new.label
      for.output$parent.splits[which(for.output$parent.splits==old.label)] <- new.label
    } 
  }
  
  if(!is.null(file.name)){
    write.csv(for.output[,1:6], file.name, row.names = F, quote=F)
  }
  return(for.output)
}

check.contiguous <- function(tt, patients, splits.for.patients, patients.for.splits){
  if(length(patients)!=2){
    stop("Not implemented")
  }
  pat.1.id <- patients[1]
  pat.2.id <- patients[2]
  
  OK <- TRUE
  all.nodes <-  c(splits.for.patients[[pat.1.id]], splits.for.patients[[pat.2.id]])
  
  for(node.1 in seq(1, length(all.nodes))){
    for(node.2 in seq(1, length(all.nodes))){
      if(node.1 < node.2){
        node.1.id <- all.nodes[node.1]
        node.2.id <- all.nodes[node.2]
        path <- get.tt.path(tt, node.1.id, node.2.id)
        for(node in path){
          if(!grepl("^unsampled_region",node)){
            if(!(patients.for.splits[[node]] %in% c(pat.1.id, pat.2.id))){
              OK <- FALSE
              break
            }
          }
        }
      }
      if(!OK){
        break
      }
    }
    if(!OK){
      break
    }
  }
  
  return(OK)
}

check.uninterrupted <- function(tt, patients, splits.for.patients, patients.for.splits){
  # this could certainly be faster
  if(length(patients)!=2){
    stop("Not implemented")
  }
  pat.1.id <- patients[1]
  pat.2.id <- patients[2]
  
  nodes.1 <- splits.for.patients[[pat.1.id]]
  nodes.2 <- splits.for.patients[[pat.2.id]]
  
  any.contiguity <- F
  any.interruption <- F
  
  for(node.1 in nodes.1){
    current.node <- get.tt.parent(tt, node.1)
    while(current.node!="root" & (grepl("^unsampled_region",current.node) | (if(is.null(patients.for.splits[[current.node]])) {T} else {patients.for.splits[[current.node]]==pat.1.id}))){

      current.node <- get.tt.parent(tt, current.node)
    }
    if(current.node != "root"){
      if(patients.for.splits[[current.node]]==pat.2.id){
        any.contiguity <- T
      } else {
        # this is a blocking node
        chain <- get.tt.ancestors(tt, current.node)
        if(length(intersect(c(pat.1.id, pat.2.id), unlist(patients.for.splits[chain])))>0){
          any.interruption <- T
          break
        }
      }
    }
  }
  
  if(!any.interruption){
    for(node.2 in nodes.2){
      current.node <- get.tt.parent(tt, node.2)
      while(current.node!="root" & (grepl("unsampled_region",current.node) | (if(is.null(patients.for.splits[[current.node]])) {T} else {patients.for.splits[[current.node]]==pat.2.id}))){
        current.node <- get.tt.parent(tt, current.node)
      }
      if(current.node != "root"){
        if(patients.for.splits[[current.node]]==pat.1.id){
          any.contiguity <- T
        } else {
          # this is a blocking node
          chain <- get.tt.ancestors(tt, current.node)
          if(length(intersect(c(pat.1.id, pat.2.id), unlist(patients.for.splits[chain])))>0){
            any.interruption <- T
            break
          }
        }
      }
      
    }
  }
  return(any.contiguity & !any.interruption)
}


extract.tt.subtree <- function(tt, patients, splits.for.patients, patients.for.splits){
  # for now, at least
  
  if(length(patients)!=2){
    stop("Not implemented")
  }
  if(!check.contiguous(tt, patients, splits.for.patients, patients.for.splits)){
    stop("Not contiguous")
  }
  
  pat.1.id <- patients[1]
  pat.2.id <- patients[2]
  
  pat.1.splts <- splits.for.patients[[pat.1.id]]
  pat.2.splts <- splits.for.patients[[pat.2.id]]
  
  sub.tt <- tt[which(tt$unique.splits %in% c(pat.1.splts, pat.2.splts)),]
  unsampled.below <- tt[which(tt$patients == "unsampled_region" & (tt$parent.splits %in% c(pat.1.splts, pat.2.splts))),]
  unsampled.above <- tt[which(tt$patients == "unsampled_region" & (tt$unique.splits %in% sub.tt$parent.splits)),]
  
  none.but.maybe.relevant <- c(unsampled.above$unique.splits, unsampled.below$unique.splits)
  
  adjacent.relevance.count <- sapply(none.but.maybe.relevant, function(x) length(intersect(c(pat.1.splts, pat.2.splts), get.tt.adjacent(tt, x) ) ))
  
  return(c(sub.tt$unique.splits, unique(none.but.maybe.relevant[which(adjacent.relevance.count > 1)])))
}

# every distance between subtrees

all.subtree.distances <- function(tree, tt, splits, assocs, slow=F){
  
  if(!slow){
    tree.dist <- dist.nodes(tree)
  } else {
    depths <- node.depth.edgelength(tree)
  }
  
  temp <- matrix(ncol = length(splits), nrow=length(splits))
  for(spt.1.no in 1:length(splits)){
    for(spt.2.no in 1:length(splits)){
      if(spt.1.no==spt.2.no){
        temp[spt.1.no, spt.2.no] <- 0
      } else if(spt.1.no<spt.2.no){
        spt.1 <- splits[spt.1.no]
        spt.2 <- splits[spt.2.no]
        
        chain.1 <- get.tt.ancestors(tt, spt.1)
        chain.2 <- get.tt.ancestors(tt, spt.2)
        
        if(spt.1 %in% chain.2){
          mrca.2 <- tt$root.nos[which(tt$unique.splits==spt.2)]
          current.node <- mrca.2
          length <- 0
          while(assocs[[current.node]]!=spt.1){
            length <- length + get.edge.length(tree, current.node)
            current.node <- Ancestors(tree, current.node, type="parent")
            if(is.root(tree, current.node)){
              stop("Big problem here (1)")
            }
          }
          temp[spt.1.no, spt.2.no] <- length
          
        } else if(spt.2 %in% chain.1){
          mrca.1 <- tt$root.nos[which(tt$unique.splits==spt.1)]
          current.node <- mrca.1
          length <- 0
          while(assocs[[current.node]]!=spt.2){
            length <- length + get.edge.length(tree, current.node)
            current.node <- Ancestors(tree, current.node, type="parent")
            if(is.root(tree, current.node)){
              stop("Big problem here (1)")
            }
          }
          temp[spt.1.no, spt.2.no] <- length
          
        } else {
          
          mrca.1 <- tt$root.nos[which(tt$unique.splits==spt.1)]
          mrca.2 <- tt$root.nos[which(tt$unique.splits==spt.2)]
          if(!slow){
            temp[spt.1.no, spt.2.no] <- tree.dist[mrca.1, mrca.2]
          } else {
            temp[spt.1.no, spt.2.no] <- pat.dist(tree, depths, mrca.1, mrca.2)
          }
        }
        temp[spt.2.no, spt.1.no] <- temp[spt.1.no, spt.2.no]
      }
    }
  }
  
  
  colnames(temp) <- splits
  rownames(temp) <- splits
  return(temp)
}

pat.dist <- function(tree, depths, node.1, node.2){
  node.1.chain <- Ancestors(tree, node.1, type="all")
  node.2.chain <- Ancestors(tree, node.2, type="all")
  if(node.1 %in% node.2.chain){
    mrca <- node.1
    dist <- depths[node.2] - depths[node.1]
  } else if(node.2 %in% node.1.chain){
    mrca <- node.2
    dist <- depths[node.1] - depths[node.2]
  } else {
    mrca <- node.1.chain[which(node.1.chain %in% node.2.chain)][1]
    dist <- (depths[node.1] - depths[mrca]) + (depths[node.2] - depths[mrca])
  }
  return(dist)
}

# are each pair of subtrees adjacent? If !none.matters then two nodes separated only by "none" are still adjacent

subtrees.adjacent <- function(tt, splits, none.matters = F){
  out <- matrix(ncol = length(splits), nrow=length(splits))
  for(spt.1.no in 1:length(splits)){
    for(spt.2.no in 1:length(splits)){
      if(spt.1.no==spt.2.no){
        out[spt.1.no, spt.2.no] <- NA
      } else if(spt.1.no<spt.2.no){
        spt.1 <- splits[spt.1.no]
        spt.2 <- splits[spt.2.no]
        
        if(spt.1 %in%  get.tt.adjacent(tt, spt.2)){
          out[spt.1.no, spt.2.no] <- T
          out[spt.2.no, spt.1.no] <- T
        } else if(!none.matters) {
          path <- get.tt.path(tt, spt.1, spt.2)
          internal.path <- path[2:(length(path)-1)]
          adj <- length(internal.path)==1 & grepl("unsampled_region",internal.path[1])
          out[spt.1.no, spt.2.no] <- adj
          out[spt.2.no, spt.1.no] <- adj
        } else {
          out[spt.1.no, spt.2.no] <- F
          out[spt.2.no, spt.1.no] <- F
        }
      }
    }
  }
  colnames(out) <- splits
  rownames(out) <- splits
  
  return(out)
}

# are pairs of subtrees from two patients not separated by any other subtrees from either of those patients?

subtrees.unblocked <- function(tt, splits){
  
  out <- matrix(ncol = length(splits), nrow=length(splits))
  for(spt.1.no in 1:length(splits)){
    for(spt.2.no in 1:length(splits)){
      if(spt.1.no==spt.2.no){
        out[spt.1.no, spt.2.no] <- NA
      } else if(spt.1.no<spt.2.no){
        
        spt.1 <- splits[spt.1.no]
        spt.2 <- splits[spt.2.no]
        
        pat.1 <- strsplit(spt.1, "-S")[[1]][1]
        pat.2 <- strsplit(spt.2, "-S")[[1]][1]
        
        if(spt.1 %in%  get.tt.adjacent(tt, spt.2)){
          out[spt.1.no, spt.2.no] <- T
          out[spt.2.no, spt.1.no] <- T
        } else {
          path <- get.tt.path(tt, spt.1, spt.2)
          internal.path <- path[2:(length(path)-1)]
          blockers <- which(grepl(paste0('^',pat.1),internal.path) | grepl(paste0('^',pat.2),internal.path))
          
          adj <- length(blockers) == 0
          out[spt.1.no, spt.2.no] <- adj
          out[spt.2.no, spt.1.no] <- adj
        } 
      }
    }
  }
  colnames(out) <- splits
  rownames(out) <- splits
  
  return(out)
}