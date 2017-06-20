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

# Output the collapsed tree from a phylogeny and its node associations

output.trans.tree <- function(tree, assocs){
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
    paste("unsampled_region-S", 1:length(unsampled.roots), sep="")
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    if(assocs.vec[node.no]=="none" & !first.of.split[node.no]){
      current.node.no <- node.no
      while(!first.of.split[current.node.no]){
        current.node.no <- Ancestors(tree, current.node.no, type="parent")
      }
      if(!grepl("^unsampled_region", splits.vec[current.node.no])){      
        stop("Parent of unsampled node is not one of the unsampled subgraph roots")
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
  
  patients <-  unlist(lapply(strsplit(unique.splits, "-S"), `[[`, 1)) 
  parent.patients <-  unlist(lapply(strsplit(parent.splits, "-S"), `[[`, 1)) 
  
  tt.table <- data.frame(unique.splits, parent.splits, patients, parent.patients, lengths, root.nos, stringsAsFactors = F)

  return(tt.table)
}

prune.unsampled.tips <- function(tt.table){
  
  for.output <- tt.table[,1:6]
  
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
    new.label <- paste("UnsampledRegion-S",x,sep="")
    for.output$unique.splits[unsampled.rows[x]] <- new.label
    for.output$parent.splits[which(for.output$parent.splits==old.label)] <- new.label
  } 
  
  for.output$patients[which(grepl("^unsampled_region",for.output$patients))] <- "UnsampledRegion"
  for.output$parent.patients[which(grepl("^unsampled_region",for.output$parent.patients))] <- "UnsampledRegion"

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
              stop("Reached the root?")
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
              stop("Reached the root?")
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

check.adjacency <- function(tt, patients, splits.for.patients){
  # this could certainly be faster
  if(length(patients)!=2){
    stop("Not implemented")
  }
  pat.1.id <- patients[1]
  pat.2.id <- patients[2]
  
  nodes.1 <- splits.for.patients[[pat.1.id]]
  nodes.2 <- splits.for.patients[[pat.2.id]]
  
  for(node.1 in nodes.1){
    for(node.2 in nodes.2){
      if(check.tt.node.adjacency(tt, node.1, node.2, T)){
        return(T)
      }
    }
  } 
  
  return(F)
}

check.tt.node.adjacency <- function(tt, label1, label2, allow.unsampled = F){
  path <- get.tt.path(tt, label1, label2)
  
  if(!allow.unsampled){
    return(length(path)==2)
  }
  
  if(length(path)==2){
    # they must be next to each other
    return(T)
  }
  
  if(length(path)>3){
    # they can't be - a path length greater than 2 can only go through an unsampled region, and adjacent unsampled regions are not allowed
    return(F)
  }
  
  return(substr(path[2], 1, 16) =="unsampled_region")
  
}

classify <- function(tree.info, verbose = F) {	
  
  if(is.null(tree.info$tree)){
    if (verbose) cat("Reading tree file ", tree.info$tree.file.name, "...\n", sep = "")
    
    pseudo.beast.import <- read.beast(tree.info$tree.file.name)
    tree <- attr(pseudo.beast.import, "phylo")
    
    if (verbose) cat("Reading annotations...\n")
    
    annotations <- attr(pseudo.beast.import, "stats")
    
    annotations$INDIVIDUAL <- as.character(annotations$INDIVIDUAL)
    annotations$SPLIT <- as.character(annotations$SPLIT)
  } else {
    tree <- tree.info$tree
    
    annotations <- data.frame(node = seq(1, length(tree$tip.label) + tree$Nnode), 
                              INDIVIDUAL = as.character(attr(tree, "INDIVIDUAL")), 
                              SPLIT = as.character(attr(tree, "SPLIT")), 
                              stringsAsFactors = F)
  }
  
  if(is.null(tree.info$splits.table)){
    if (verbose) cat("Reading splits file",splits.file.name,"...\n")
    
    splits <- read.csv(tree.info$splits.file.name, stringsAsFactors = F)
  } else {
    splits <- tree.info$splits.table
  }
  
  if (verbose) cat("Collecting tips for each patient...\n")
  
  patients <- unique(splits$patient)
  patients <- patients[patients!="unsampled"]
  
  all.splits <- unique(splits$subgraph)
  all.splits <- all.splits[all.splits!="unsampled"]
  
  in.order <- match(seq(1, length(tree$tip.label) + tree$Nnode), annotations$node)
  
  assocs <- annotations$SPLIT[in.order]
  assocs <- lapply(assocs, function(x) replace(x, is.na(x), "none"))
  assocs <- lapply(assocs, function(x) replace(x, x=="unsampled", "none"))
  
  splits.for.patients <- lapply(patients, function(x) unique(splits$subgraph[which(splits$patient==x)] ))
  names(splits.for.patients) <- patients
  
  patients.for.splits <- lapply(all.splits, function(x) unique(splits$patient[which(splits$subgraph==x)] ))
  names(patients.for.splits) <- all.splits
  
  patients.included <- patients
  
  total.pairs <- (length(patients.included) ^ 2 - length(patients.included))/2
  
  if (verbose) cat("Collapsing subtrees...\n")
  
  tt <- output.trans.tree(tree, assocs)
  
  if (verbose) cat("Identifying pairs of unblocked splits...\n")
  
  collapsed.adjacent <- subtrees.unblocked(tt, all.splits)
  
  if (verbose) cat("Calculating pairwise distances between splits...\n")
  
  split.distances <- tryCatch(
    all.subtree.distances(tree, tt, all.splits, assocs), warning=function(w){return(NULL)}, error=function(e){return(NULL)})
  
  if(is.null(split.distances)){
    split.distances <- all.subtree.distances(tree, tt, all.splits, assocs, TRUE)
  }
  
  if (verbose) cat("Testing pairs...\n")
  
  count <- 0
  adjacency.matrix <- matrix(NA, length(patients.included), length(patients.included))
  contiguity.matrix <- matrix(NA, length(patients.included), length(patients.included))
  path.matrix <- matrix(NA, length(patients.included), length(patients.included))
  nodes.1.matrix <- matrix(NA, length(patients.included), length(patients.included))
  nodes.2.matrix <- matrix(NA, length(patients.included), length(patients.included))
  dir.12.matrix <- matrix(NA, length(patients.included), length(patients.included))
  dir.21.matrix <- matrix(NA, length(patients.included), length(patients.included))
  min.distance.matrix <- matrix(NA, length(patients.included), length(patients.included))
  
  for(pat.1 in seq(1, length(patients.included))){
    for(pat.2 in  seq(1, length(patients.included))){
      if (pat.1 < pat.2) {
        
        count <- count + 1
        
        pat.1.id <- patients.included[pat.1]
        pat.2.id <- patients.included[pat.2]
        
        nodes.1 <- splits.for.patients[[pat.1.id]]
        nodes.2 <- splits.for.patients[[pat.2.id]]
        
        all.nodes <-  c(nodes.1, nodes.2)
        
        adjacency.matrix[pat.1, pat.2] <- check.adjacency(tt, c(pat.1.id, pat.2.id), splits.for.patients)		
        contiguity.matrix[pat.1, pat.2] <- check.contiguous(tt, c(pat.1.id, pat.2.id), splits.for.patients, patients.for.splits)		
        
        count.12 <- 0
        count.21 <- 0
        
        for(node.2 in nodes.2){
          ancestors <- get.tt.ancestors(tt, node.2)
          if(length(intersect(ancestors, nodes.1)) > 0){
            count.12 <- count.12 + 1
          }
        }
        
        for(node.1 in nodes.1){
          ancestors <- get.tt.ancestors(tt, node.1)
          if(length(intersect(ancestors, nodes.2)) > 0){
            count.21 <- count.21 + 1
          }
        }
        
        prop.12 <- count.12/length(nodes.2)
        prop.21 <- count.21/length(nodes.1)
        
        nodes.1.matrix[pat.1, pat.2] <- length(nodes.1)
        nodes.2.matrix[pat.1, pat.2] <- length(nodes.2)
        
        dir.12.matrix[pat.1, pat.2] <- count.12
        dir.21.matrix[pat.1, pat.2] <- count.21
        
        if(count.12 == 0 & count.21 == 0){
          path.matrix[pat.1, pat.2] <- "none"
        } else if(count.12 != 0 & count.21 == 0 & prop.12 == 1) {
          if(count.12 == 1){
            path.matrix[pat.1, pat.2] <- "anc"
          } else {
            path.matrix[pat.1, pat.2] <- "multiAnc"
          }
        } else if(count.21 != 0 & count.12 == 0 & prop.21 == 1) {
          if(count.21 == 1){
            path.matrix[pat.1, pat.2] <- "desc"
          } else {
            path.matrix[pat.1, pat.2] <- "multiDesc"
          }
        } else {
          path.matrix[pat.1, pat.2] <- "conflict"
        }
        
        pairwise.distances <- vector()
        
        for(node.1 in nodes.1){
          for(node.2 in nodes.2){
            if(collapsed.adjacent[node.1, node.2]){
              pairwise.distances <- c(pairwise.distances, split.distances[node.1, node.2])
            }
          }
        }
        
        min.distance.matrix[pat.1, pat.2] <- min(pairwise.distances)
        
        
        if (count %% 100 == 0) {
          if (verbose) cat(paste("Done ",format(count, scientific = F)," of ",format(total.pairs, scientific = F)," pairwise calculations\n", sep = ""))
        }
      }
    }
  }
  
  normalisation.constant <- tree.info$normalisation.constant
  
  normalised.distance.matrix<- min.distance.matrix/normalisation.constant
  
  adjacency.df <- convert.to.columns(adjacency.matrix, patients.included)
  contiguity.df <- convert.to.columns(contiguity.matrix, patients.included)
  dir.12.df <- convert.to.columns(dir.12.matrix, patients.included)
  dir.21.df <- convert.to.columns(dir.21.matrix, patients.included)
  nodes.1.df <- convert.to.columns(nodes.1.matrix, patients.included)
  nodes.2.df <- convert.to.columns(nodes.2.matrix, patients.included)
  path.df <- convert.to.columns(path.matrix, patients.included)
  min.distance.df <- convert.to.columns(min.distance.matrix, patients.included)
  if(normalisation.constant!=1){
    normalised.distance.df <- convert.to.columns(normalised.distance.matrix, patients.included)
  }
  
  classification <- cbind(adjacency.df, contiguity.df[,3], dir.12.df[,3], dir.21.df[,3], nodes.1.df[,3], nodes.2.df[,3], path.df[,3], min.distance.df[,3])
  
  column.names <- c("Patient_1", "Patient_2", "adjacent", "contiguous", "paths12", "paths21", "nodes1", "nodes2", "path.classification", "min.distance.between.subtrees")
  
  if(normalisation.constant!=1){
    classification <- cbind(classification, normalised.distance.df[,3])
    column.names <- c(column.names, "normalised.min.distance.between.subtrees")
  }
  
  colnames(classification) <- column.names
  return (list(classification = classification, collapsed=tt))
}

convert.to.columns <- function(matrix, names){
  a.table <- as.table(matrix)
  colnames(a.table) <- names
  rownames(a.table) <- names
  a.df <- as.data.frame(a.table, stringsAsFactors=F)
  keep <- complete.cases(a.df)
  a.df <- a.df[keep,]
  return(a.df)
}


summarise.classifications <- function(all.tree.info, hosts, min.threshold, dist.threshold, allow.mt = T, csv.fe = "csv", verbose = F){
  #
  # Read classification files
  #

  tt	<-	lapply(all.tree.info, function(tree.info){
    
    if(is.null(tree.info$classification.results$classification) & is.null(tree.info$classification.file.name)){
      NULL
    }
    
    if(!is.null(tree.info$classification.results$classification)){
      tt <- as.data.table(tree.info$classification.results$classification)
    } else {
      if (verbose) cat("Reading window input file ",tree.info$classification.file.name,"\n", sep="")
      tt <- as.data.table(read.table(tree.info$classification.file.name, sep=",", header=TRUE, stringsAsFactors=FALSE))
    }
    tt[, SUFFIX:=tree.info$suffix]			
    tt
  })	
  
  if(length(tt)==0){
    stop("No classification results present in any window; cannot continue.\n")
  }
  
  # need to worry about transmissions listed in the wrong direction in the input file
  # to avoid very large memory, consolidate by FILE
  
  if(verbose) cat("Rearranging patient pairs...\n")
  
  tt	<- lapply(tt, function(x){
    #x	<- tt[[1]]
    tmp	<- copy(x)
    setnames(tmp, c('Patient_1','Patient_2','paths12','paths21'), c('Patient_2','Patient_1','paths21','paths12'))
    set(tmp, tmp[, which(path.classification=="anc")], 'path.classification', 'TMP')
    set(tmp, tmp[, which(path.classification=="desc")], 'path.classification', 'anc')
    set(tmp, tmp[, which(path.classification=="TMP")], 'path.classification', 'desc')
    set(tmp, tmp[, which(path.classification=="multiAnc")], 'path.classification', 'TMP')
    set(tmp, tmp[, which(path.classification=="multiDesc")], 'path.classification', 'multiAnc')
    set(tmp, tmp[, which(path.classification=="TMP")], 'path.classification', 'multiDesc')
    x	<- rbind(x, tmp)
    setkey(x, Patient_1, Patient_2)
    subset(x, Patient_1 < Patient_2)
  })
  #
  # rbind consolidated files
  #
  
  if(verbose) cat("Consolidating file contents...\n")
  
  tt	<- do.call('rbind',tt)
  
  if(verbose) cat("Finding patristic distance columns...\n")
  
  # reset names depending on which Classify script was used
  if(any('normalised.min.distance.between.subtrees'==colnames(tt))){
    setnames(tt, 'normalised.min.distance.between.subtrees', 'PATRISTIC_DISTANCE')
  } else if(any('min.distance.between.subtrees'==colnames(tt))){
    setnames(tt, 'min.distance.between.subtrees', 'PATRISTIC_DISTANCE')
  }
  
  setnames(tt, c('Patient_1','Patient_2','path.classification','paths21','paths12','adjacent'), c('PAT.1','PAT.2','TYPE','PATHS.21','PATHS.12','ADJACENT'))
  
  # change type name depending on allow.mt
  if(!allow.mt){
    if(verbose) cat("Allowing only single lineage transmission...\n")
    set(tt, tt[, which(TYPE%in%c("multiAnc", "multiDesc"))], 'TYPE', 'conflict')
  }
  
  #	check we have patristic distances, paths
  
  stopifnot( !nrow(subset(tt, is.na(PATRISTIC_DISTANCE))) )
  stopifnot( !nrow(subset(tt, is.na(PATHS.12))) )
  stopifnot( !nrow(subset(tt, is.na(PATHS.21))) )
  
  #	set to numeric
  set(tt, NULL, 'PATRISTIC_DISTANCE', tt[, as.numeric(PATRISTIC_DISTANCE)])
  # 	add window coordinates
  
  #	rename TYPE
  set(tt, tt[,which(TYPE=='anc')], 'TYPE', 'anc_12')
  set(tt, tt[,which(TYPE=='desc')], 'TYPE', 'anc_21')
  set(tt, tt[,which(TYPE=='multiAnc')], 'TYPE', 'multi_anc_12')
  set(tt, tt[,which(TYPE=='multiDesc')], 'TYPE', 'multi_anc_21')
  
  if(verbose) cat("Reordering...\n")
  
  #	reorder
  setkey(tt, SUFFIX, PAT.1, PAT.2)
  
  if (verbose) cat("Making summary output table...\n")

  set(tt, NULL, c('PATHS.12','PATHS.21'),NULL)
  
  
  existence.counts <- tt[, list(both.exist=length(SUFFIX)), by=c('PAT.1','PAT.2')]

  tt <- merge(tt, existence.counts, by=c('PAT.1', 'PAT.2'))
  
  tt.close <- tt[which(tt$ADJACENT & tt$PATRISTIC_DISTANCE < dist.threshold ),]
   
  tt.close$NOT.SIBLINGS <- tt.close$ADJACENT & (tt.close$PATRISTIC_DISTANCE < dist.threshold) & tt.close$TYPE!="none"
  
  # How many windows have this relationship, ADJACENT and PATRISTIC_DISTANCE below the threshold?
  type.counts	<- tt.close[, list(windows=length(SUFFIX)), by=c('PAT.1','PAT.2','TYPE')]
  # How many windows have ADJACENT and PATRISTIC_DISTANCE below the threshold?
  any.counts <- tt.close[, list(all.windows=length(SUFFIX)), by=c('PAT.1','PAT.2')]
  # How many windows have a relationship other than "none", ADJACENT and PATRISTIC_DISTANCE below the threshold?
  ns.counts <- tt.close[, list(ns.windows=length(which(NOT.SIBLINGS))), by=c('PAT.1','PAT.2')]
  
  tt.close		<- merge(tt.close, type.counts, by=c('PAT.1','PAT.2','TYPE'))
  tt.close		<- merge(tt.close, any.counts, by=c('PAT.1','PAT.2'))
  tt.close		<- merge(tt.close, ns.counts, by=c('PAT.1','PAT.2'))
  
  tt.close[, fraction:=paste(windows,'/',both.exist,sep='')]
  #	convert "anc_12" and "ans_21" to "anc" depending on direction
  tt.close[, DUMMY:=NA_character_]
  tmp			<- tt.close[, which(TYPE=="anc_12")]
  set(tt.close, tmp, 'TYPE', "trans")
  tmp			<- tt.close[, which(TYPE=="anc_21")]
  set(tt.close, tmp, 'DUMMY', tt.close[tmp, PAT.1])
  set(tt.close, tmp, 'PAT.1', tt.close[tmp, PAT.2])
  set(tt.close, tmp, 'PAT.2', tt.close[tmp, DUMMY])
  set(tt.close, tmp, 'TYPE', "trans")
  
  tmp			<- tt.close[, which(TYPE=="multi_anc_12")]
  set(tt.close, tmp, 'TYPE', "multi_trans")
  tmp			<- tt.close[, which(TYPE=="multi_anc_21")]
  set(tt.close, tmp, 'DUMMY', tt.close[tmp, PAT.1])
  set(tt.close, tmp, 'PAT.1', tt.close[tmp, PAT.2])
  set(tt.close, tmp, 'PAT.2', tt.close[tmp, DUMMY])
  set(tt.close, tmp, 'TYPE', "multi_trans")
  
  tt.close[, DUMMY:=NULL]
  
  set(tt.close, NULL, c('SUFFIX', 'ADJACENT','PATRISTIC_DISTANCE'), NULL)
  tt.close <- tt.close[!duplicated(tt.close),]
  
  #	write to file
  setkey(tt.close, PAT.1, PAT.2, TYPE)
  #
  
  return(subset(tt.close, all.windows>=min.threshold))
}

