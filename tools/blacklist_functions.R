# Blacklister for exact duplicates. The entries argument is a list, the entries of each being groups of tips whose reads are identical.

blacklist.exact.duplicates <- function(tree.info, raw.threshold, ratio.threshold, tip.regex, verbose = F){
  
  if(is.null(tree.info$duplicate.tips)){
    stop("No duplicate tips for tree suffix ",tree.info$suffix)
  }
  
  entries <- tree.info$duplicate.tips
  
  pairs.table	<- data.frame(V1=character(0), V2=character(0))
  
  for(entry in entries){
    # get rid of anything that doesn't match the regexp (outgroup etc)
    tmp <-  entry[!is.na(sapply(entry, read.count.from.label, regexp = tip.regex))]
    if(length(tmp)>1){
      tmp <- t(combn(tmp,2))
      colnames(tmp) <- c('V1','V2')
      pairs.table <- rbind(pairs.table, tmp)  
    }  
  }
  
  pairs.table$V1 <- as.character(pairs.table$V1)
  pairs.table$V2 <- as.character(pairs.table$V2)
  
  tmp <- unlist(sapply(pairs.table$V1, function(x) read.count.from.label(x, tip.regex)))
  if(!is.null(tmp))
    pairs.table$reads.1 <- tmp
  if(is.null(tmp))
    pairs.table$reads.1 <- integer(0)
  
  tmp <- unlist(sapply(pairs.table$V2, function(x) read.count.from.label(x, tip.regex)))
  if(!is.null(tmp))
    pairs.table$reads.2 <- tmp
  if(is.null(tmp))
    pairs.table$reads.2 <- integer(0)
  
  tmp <- unlist(sapply(pairs.table$V1, function(x) host.from.label(x, tip.regex)))
  if(!is.null(tmp))
    pairs.table$host.1 <- tmp
  if(is.null(tmp))
    pairs.table$host.1 <- character(0)
  
  tmp <- unlist(sapply(pairs.table$V2, function(x) host.from.label(x, tip.regex)))
  if(!is.null(tmp))
    pairs.table$host.2 <- tmp
  if(is.null(tmp))
    pairs.table$host.2 <- character(0)
  
  # pairs from the same host aren't under consideration
  
  pairs.table <- pairs.table[which(pairs.table$host.1 != pairs.table$host.2),]
  
  # reverse the order so the read with the greater count is in the first column
  pairs.table[pairs.table$reads.2 > pairs.table$reads.1, c("V1", "V2", "reads.1", "reads.2")] <- pairs.table[pairs.table$reads.2 > pairs.table$reads.1, c("V2", "V1", "reads.2", "reads.1")] 
  
  # calculate the counts
  pairs.table$ratio <- pairs.table$reads.2 / pairs.table$reads.1
  if (verbose) cat("Tree suffix ",tree.info$suffix,": making blacklist with a ratio threshold of ",ratio.threshold," and a raw threshold of ",raw.threshold,"\n",sep="")
  blacklisted <- pairs.table[which(pairs.table$ratio < ratio.threshold | (pairs.table$reads.2<raw.threshold)),2]
  
  return(blacklisted)
}


# Strip the tree down to just reads from one host and an outgroup. Do a parsimony reconstruction on that tree and return the subgraphs.

get.splits.for.host <- function(host, tip.hosts, tree, root.name, raw.threshold, ratio.threshold, sankoff.k, check.duals, no.read.counts = T, verbose=F){
  if (verbose) cat("Identifying splits for host ", host, "\n", sep="")
  
  blacklist.items <- vector()
  
  dual <- F
  
  if(length(which(tip.hosts==host))>1){
    
    # Keep only the tips from this host and the outgroup
    if(!is.null(root.name)){
      outgroup.no <- which(tree$tip.label==root.name)
    } else {
      # Find a random outgroup if you can
      
      host.mrca <- mrca.phylo(tree, which(tip.hosts==host))
      if(host.mrca == getRoot(tree)){
        stop(paste0("No suitable outgroup found for host ",host,"; try adding one and specifying its tip name with --outgroupName"))
      }
      outgroup.no <- sample(1:length(tree$tip.label), 1)
      while(host.mrca %in% Ancestors(tree, outgroup.no)){
        outgroup.no <- sample(1:length(tree$tip.label), 1)
      }
      root.name <- tree$tip.label[outgroup.no]
    }
    
    tip.nos <- c(outgroup.no, which(tip.hosts==host))
    
    subtree <- drop.tip(tree, setdiff(1:length(tree$tip.label), tip.nos))
    
    # Re-find the root
    
    st.outgroup.no <- which(subtree$tip.label==root.name)
    
    # di2multi actually unroots the tree!
    
    subtree <- root(subtree, st.outgroup.no, resolve.root = T)
    
    # Set up split.and.annotate (which takes a list)
    
    host.tips <- list()
    host.tips[[host]] <- setdiff(1:length(subtree$tip.label), st.outgroup.no)
    
    tip.hosts <- rep(host, length(subtree$tip.label))
    tip.hosts[st.outgroup.no] <- "unsampled"
    
    # Perform split.and.annotate; get a list of splits
    
    split.results <- split.and.annotate(subtree, c(host, "unsampled"), tip.hosts, host.tips, NULL, vector(), tip.regex, "s", rep(1, length(subtree$tip.label)), sankoff.k, 0, useff = F, verbose)
    
    # vector of of split IDs
    
    host.split.ids <- split.results$split.hosts
    
    # list of tips that belong to each split ID
    
    tips.for.splits <- split.results$split.tips
    
    # The total number of reads is either got from the tips or assumed to be one per tip
    
    if(!no.read.counts){
      reads.per.tip <- sapply(subtree$tip.label, function(x) read.count.from.label(x, tip.regex))
    } else {
      reads.per.tip <- rep(1, length(subtree$tip.label))
      reads.per.tip[st.outgroup.no] <- NA
    }
    
    total.reads <- sum(reads.per.tip[!is.na(reads.per.tip)])
    
    if(length(host.split.ids)==1){
      
      # if the algorithm didn't split the tips
      if(verbose){
        cat("No splits for host ", host, "\n", sep="")
      }
      if(total.reads < raw.threshold){
        # blacklist the whole host if the read count from its only subgraph is below the raw threshold
        
        if(verbose){
          cat("Blacklisting ",host,"; not enough reads in total.\n", sep="")
        }
        blacklist.items <- subtree$tip.label[which(subtree$tip.label!=root.name)]
      }
      
    } else {
      # at least two subgraphs
      
      if (verbose) cat(length(host.split.ids), " splits for host ", host, "\n", sep="")
      
      props <- vector()
      too.small <- vector()
      
      tips.vector <- vector()
      read.count.vector <- vector()
      tip.count.vector <- vector()
      
      too.small <- sapply(host.split.ids, function(x) check.read.count.for.split(x, tips.for.splits, raw.threshold, ratio.threshold, reads.per.tip, total.reads))
      
      if(length(which(!too.small))>1 & check.duals){
        # there are at least two subgraphs which are too big to be contaminants
        
        if(verbose){
          cat(host, " looks like a dual infection\n", sep="")
        }
        
        dual <- T
      }
      
      # blacklist the tips from small subgraphs
      
      if(length(which(too.small))>0){
        blacklist.items <- subtree$tip.label[unlist(tips.for.splits[host.split.ids[which(too.small)]])]
      }
      
    }
    
    tips.vector <- subtree$tip.label[unlist(tips.for.splits[host.split.ids])]
    read.count.vector <- unlist(lapply(host.split.ids, function(x){
      rep(sum(read.count.from.label(subtree$tip.label[tips.for.splits[[x]]], tip.regex)), length(tips.for.splits[[x]]))
    }))
    tip.count.vector <- unlist(lapply(host.split.ids, function(x){
      rep(length(tips.for.splits[[x]]), length(tips.for.splits[[x]]))
    }))
    
  } else {
    
    # this covers the situation where there is only one host tip (because ape has trouble with two-tip trees)
    tip.no <- which(tip.hosts==host)
    
    if(!no.read.counts){
      reads <- read.count.from.label(tree$tip.label[tip.no], tip.regex)
    } else {
      reads <- 1
    }

    if(reads < raw.threshold){
      if(verbose){
        cat("Blacklisting ",host,"; not enough reads in total.\n", sep="")
      }
      blacklist.items <- tree$tip.label[tip.no]
    } else {
      if(verbose){
        cat("Keeping ",host,"; host has a single tip with sufficient reads.\n", sep="")
      }
    }
    
    tips.vector <- tree$tip.label[tip.no]
    read.count.vector <- reads
    tip.count.vector <- 1
  }
  
  
  list(id = host, blacklist.items = blacklist.items, tip.names = tips.vector, read.counts = read.count.vector, tip.counts = tip.count.vector, dual=dual)
}


# Return whether the number of reads in this split is below one of the thresholds

check.read.count.for.split <- function(split, tips.for.splits, raw.threshold, ratio.threshold, reads.per.tip, total.reads){
  # get the read counts
  
  count.in.split <- sum(reads.per.tip[tips.for.splits[[split]]])
  
  # find what proportion of all reads from this host are in this subgraph and check against the thresholds
  
  prop.in.split <- count.in.split/total.reads
  return (count.in.split < raw.threshold | prop.in.split < ratio.threshold)
}


# This 

blacklist.all.reads.for.dual.host <- function(host, tree.info){
  for(info in tree.info){
    labels <- info$tree$tip.label
    hosts <- sapply(labels, function(x) host.from.label(x, tip.regex))
    labels.to.go <- labels[hosts==host]
    info$blacklist <- c(info$blacklist, labels.to.go)
  }
}

blacklist.small.subgraphs.for.dual.host <- function(host, tree.info, max.reads){
  for(info in tree.info){
    labels <- info$tree$tip.label
    duals.table <- info$duals.info
    pat.rows <- duals.table[duals.table$host==host,]
    
    if(nrow(pat.rows)>0){
      max.reads <- max(pat.rows$reads.in.subtree)
      smaller.tips <- pat.rows$tip.name[pat.rows$reads.in.subtree!=max.reads]
      info$blacklist <- c(info$blacklist, smaller.tips)
    }
  }
}


# This takes a tree, a blacklist, and a string, appends that string to the blacklisted tips of the tree

blacklist.tip.rename <- function(tree, blacklist, reason){
  tree$tip.label[blacklist] <- paste0(tree$tip.label[blacklist], "_", reason)
  
  return(tree)
}

blacklist.duals <- function(all.tree.info, hosts, threshold = 1, summary.file=NULL, verbose=F){
  
  if(verbose) cat("Calculating dual window fractions...\n")
  
  #	Count number of potential dual windows by host
  
  fractions <- lapply(hosts, function(x){
    rowSums(sapply(all.tree.info, function(y){
      num <- x %in% y$duals.info$host
      denom <- x %in% y$hosts.for.tips
      c(num, denom)
    }))
  })
  
  names(fractions) <- hosts
  
  blacklists <- list()
  
  for(tree.info in all.tree.info){
    if(verbose) cat("Making new blacklists for tree suffix ",tree.info$suffix,"\n",sep="")
    
    new.bl <- vector()
    
    tip.names <- tree.info$tree$tip.label
    
    hosts.for.tips <- tree.info$hosts.for.tips
    
    for(a.host in hosts){
      if(fractions[[a.host]][1]/fractions[[a.host]][2] > threshold){
        # blacklist everything from this host
        if(verbose) cat("All tips from host ",a.host," blacklisted\n", sep="")
        
        new.bl <- unique(c(new.bl, na.omit(tip.names[hosts.for.tips==a.host])))
        
      } else {
        # blacklist smaller subgraphs
        if(verbose) cat("Tips from smaller subgraphs for host ",a.host," blacklisted\n", sep="")
        
        duals.table <- tree.info$duals.info
        if(!is.null(duals.table)){
          pat.rows <- duals.table[which(duals.table$host == a.host),]
          if(nrow(pat.rows)>0){
            max.reads <- max(pat.rows$reads.in.subtree)
            smaller.tips <- pat.rows$tip.name[pat.rows$reads.in.subtree!=max.reads]
            
            new.bl <- unique(c(new.bl, na.omit(smaller.tips)))
          }
        }
      }
    }
    blacklists[[tree.info$suffix]] <- new.bl
  }
  
  if(!is.null(summary.file)) {

    output <- lapply(hosts, function(x) c(fractions[[x]][1]/fractions[[a.host]][2], fractions[[x]][1]))
    
    proportions <- unlist(lapply(output, "[[", 1))
    counts <- unlist(lapply(output, "[[", 2))
    
    out.df <- data.frame(host = hosts, count = counts, proportion = proportions, stringsAsFactors = F)
    if(verbose) cat("Writing dual summary file to ", summary.file, "\n", sep="")
    write.csv(out.df, file=summary.file, row.names=FALSE, quote=FALSE)
  }
  
  return(blacklists)
}

