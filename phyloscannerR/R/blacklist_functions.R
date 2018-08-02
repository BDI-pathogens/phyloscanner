# Blacklister for exact duplicates. The entries argument is a list, the entries of each being groups of tips whose reads are identical.

#' @keywords internal
#' @export blacklist.exact.duplicates

blacklist.exact.duplicates <- function(ptree, raw.threshold, ratio.threshold, tip.regex, verbose = F){
  
  if(is.null(ptree$duplicate.tips)){
    stop("No duplicate tips for tree ID ",ptree$id)
  }
  
  entries <- ptree$duplicate.tips
  
  pairs.table <- entries %>% map(function(x){
    # get rid of anything that doesn't match the regexp (outgroup etc)
    tmp <-  x[!is.na(sapply(x, read.count.from.label, regexp = tip.regex))]
    if(length(tmp)>1){
      tmp <- as.tibble(t(combn(tmp,2)))
      names(tmp) <- c('tip.1','tip.2')
      
      tmp
    } else {
      NULL
    }
    
  })
  
  pairs.table	<- bind_rows(pairs.table)
  
  if(nrow(pairs.table) > 0){
    
    tmp <- unlist(sapply(pairs.table$tip.1, function(x) read.count.from.label(x, tip.regex)))
    if(!is.null(tmp))
      pairs.table$reads.1 <- tmp
    if(is.null(tmp))
      pairs.table$reads.1 <- integer(0)
    
    tmp <- unlist(sapply(pairs.table$tip.2, function(x) read.count.from.label(x, tip.regex)))
    if(!is.null(tmp))
      pairs.table$reads.2 <- tmp
    if(is.null(tmp))
      pairs.table$reads.2 <- integer(0)
    
    tmp <- unlist(sapply(pairs.table$tip.1, function(x) host.from.label(x, tip.regex)))
    if(!is.null(tmp))
      pairs.table$host.1 <- tmp
    if(is.null(tmp))
      pairs.table$host.1 <- character(0)
    
    tmp <- unlist(sapply(pairs.table$tip.2, function(x) host.from.label(x, tip.regex)))
    if(!is.null(tmp))
      pairs.table$host.2 <- tmp
    if(is.null(tmp))
      pairs.table$host.2 <- character(0)
    
    # pairs from the same host aren't under consideration
    
    pairs.table <- pairs.table %>% filter(host.1 != host.2)
    
    # reverse the order so the read with the greater count is in the first column
    
    pairs.table <- pairs.table %>% 
      mutate(new.tip.1 = ifelse(reads.2 > reads.1, tip.2, tip.1), 
             new.tip.2 = ifelse(reads.2 > reads.1, tip.1, tip.2),
             new.reads.1 = ifelse(reads.2 > reads.1, reads.2, reads.1),
             new.reads.2 = ifelse(reads.2 > reads.1, reads.1, reads.2),
             new.host.1 = ifelse(reads.2 > reads.1, host.2, host.1),
             new.host.2 = ifelse(reads.2 > reads.1, host.1, host.2),
             tip.1 = new.tip.1,
             tip.2 = new.tip.2,
             reads.1 = new.reads.1,
             reads.2 = new.reads.2,
             host.1 = new.host.1,
             host.2 = new.host.2) %>% 
      select(-new.tip.1, -new.tip.2, -new.reads.1, -new.reads.2, -new.host.1, -new.host.2)
    
    # calculate the counts
    pairs.table <- pairs.table %>% mutate(ratio = reads.2/reads.1)
    
    if (verbose) cat("Tree ID ",ptree$id,": making duplicate blacklist with a ratio threshold of ",ratio.threshold," and a raw threshold of ",raw.threshold,"\n",sep="")
    
    blacklisted <- pairs.table %>% filter((ratio < ratio.threshold) | (reads.2 < raw.threshold))
    
    blacklisted
  } else {

    pairs.table
  }
}


# Strip the tree down to just reads from one host and an outgroup. Do a parsimony reconstruction on that tree and return the subgraphs.

#' @keywords internal
#' @export get.splits.for.host
#' @importFrom ape root

get.splits.for.host <- function(host, 
                                tip.hosts, 
                                tree, 
                                outgroup.name, 
                                raw.threshold, 
                                ratio.threshold, 
                                sankoff.method = "s", 
                                sankoff.k, 
                                sankoff.p = 0, 
                                count.reads.in.parsimony, 
                                check.duals, 
                                no.read.counts = T, 
                                tip.regex, 
                                verbose = F, 
                                no.progress.bars = T, 
                                just.report.counts = F){
  
  if (verbose) cat("Identifying splits for host ", host, "\n", sep="")
  
  blacklist.items <- vector()
  
  dual <- F
  
  if(sum(tip.hosts==host, na.rm=T) > 1){
    
    # Keep only the tips from this host and the outgroup
    if(!is.null(outgroup.name)){
      outgroup.no <- which(tree$tip.label==outgroup.name)
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
      outgroup.name <- tree$tip.label[outgroup.no]
    }
    
    tip.nos <- c(outgroup.no, which(tip.hosts==host))
    
    subtree <- drop.tip(tree, setdiff(1:length(tree$tip.label), tip.nos))
    
    # Re-find the root
    
    st.outgroup.no <- which(subtree$tip.label==outgroup.name)
    
    # di2multi actually unroots the tree!
    
    subtree <- root(subtree, st.outgroup.no, resolve.root = T)
    
    # Set up split.and.annotate (which takes a list)
    
    host.tips <- list()
    host.tips[[host]] <- setdiff(1:length(subtree$tip.label), st.outgroup.no)
    
    tip.hosts <- rep(host, length(subtree$tip.label))
    tip.hosts[st.outgroup.no] <- "unassigned"
    
    if(!no.read.counts){
      reads.per.tip <- sapply(subtree$tip.label, function(x) read.count.from.label(x, tip.regex))
    } else {
      reads.per.tip <- rep(1, length(subtree$tip.label))
      reads.per.tip[st.outgroup.no] <- NA
    }
    
    multipliers <- rep(1, length(subtree$tip.label))
    
    if(count.reads.in.parsimony){
      multipliers <- reads.per.tip
      multipliers[which(is.na(multipliers))] <- 1
    }
    
    # Perform split.and.annotate; get a list of splits
    
    split.results <- split.and.annotate(subtree, c(host, "unassigned"), tip.hosts, host.tips, NULL, vector(), tip.regex, sankoff.method, multipliers, sankoff.k, sankoff.p, useff = F, verbose, no.progress.bars)
    
    # vector of of split IDs
    
    host.split.ids <- split.results$split.hosts
    
    # list of tips that belong to each split ID
    
    tips.for.splits <- split.results$split.tips
    
    # The total number of reads is either gotten from the tips or assumed to be one per tip
    
    total.reads <- sum(reads.per.tip[!is.na(reads.per.tip)])
    
    counts.per.split <- sapply(host.split.ids, function(x) sum(reads.per.tip[tips.for.splits[[x]]]), USE.NAMES=FALSE)
    
    if (just.report.counts) {return(counts.per.split)}
    
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
        blacklist.items <- subtree$tip.label[which(subtree$tip.label!=outgroup.name)]
        
        multiplicity <- 0
      } else {
        multiplicity <- 1
      }
      
    } else {
      # at least two subgraphs
      
      if (verbose) cat(length(host.split.ids), " splits for host ", host, "\n", sep="")
      
      too.small <- sapply(host.split.ids, function(x) check.read.count.for.split(x, tips.for.splits, raw.threshold, ratio.threshold, reads.per.tip, total.reads))
      
      if(length(which(!too.small))>1 & check.duals){
        # there are at least two subgraphs which are too big to be contaminants
        
        if(verbose){
          cat(host, " looks like a dual infection\n", sep="")
        }
        
        dual <- T
      }
      
      multiplicity <- length(which(!too.small))
      
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
    split.id.vector <- unlist(lapply(host.split.ids, function(x){
      rep(x, length(tips.for.splits[[x]]))
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
      multiplicity <- 0
    } else {
      if(verbose){
        cat("Keeping ",host,"; host has a single tip with sufficient reads.\n", sep="")
      }
      multiplicity <- 1
    }
    
    tips.vector <- tree$tip.label[tip.no]
    read.count.vector <- reads
    tip.count.vector <- 1
    split.id.vector <- host
  }
  
  if (just.report.counts) {return(c(reads))}
  
  list(id = host, 
       blacklist.items = blacklist.items, 
       tip.names = tips.vector, 
       read.counts = read.count.vector, 
       split.ids = split.id.vector,
       tip.counts = tip.count.vector, 
       dual = dual, 
       multiplicity = multiplicity)
}

# Return whether the number of reads in this split is below one of the thresholds

#' @keywords internal
#' @export check.read.count.for.split

check.read.count.for.split <- function(split, tips.for.splits, raw.threshold, ratio.threshold, reads.per.tip, total.reads){
  # get the read counts
  
  count.in.split <- sum(reads.per.tip[tips.for.splits[[split]]])
  
  # find what proportion of all reads from this host are in this subgraph and check against the thresholds
  
  prop.in.split <- count.in.split/total.reads
  
  count.in.split < raw.threshold | prop.in.split < ratio.threshold
}

# This takes a tree, a blacklist, and a string, appends that string to the blacklisted tips of the tree

#' @keywords internal
#' @export rename.blacklisted.tips

rename.blacklisted.tips <- function(tree, blacklist, reason){
  tree$tip.label[blacklist] <- paste0(tree$tip.label[blacklist], "_X_", reason)
  
  tree
}

#' @keywords internal
#' @export blacklist.duals

blacklist.duals <- function(ptrees, hosts, threshold = 1, summary.file=NULL, verbose=F){
  
  if(verbose) cat("Calculating dual window fractions...\n")
  
  #	Count number of potential dual windows by host
  
  fractions <- tibble(host = hosts)
  fractions <- fractions %>% mutate(numerator = map_int(host, function(x){
    sum(map_int(ptrees, function(y){
      x %in% y$duals.info$host
    }))
  }), 
  denominator = map_int(host, function(x){
    sum(map_int(ptrees, function(y){
      x %in% y$hosts.for.tips
    }))
  }), 
  proportion = numerator/denominator, 
  blacklist.entirely = proportion > threshold)
  
  blacklists <- list()
  
  for(ptree in ptrees){
    if(verbose) cat("Making new blacklists for tree ID ",ptree$id,"\n",sep="")
    
    tip.names <- ptree$tree$tip.label
    
    hosts.for.tips <- ptree$hosts.for.tips
    
    new.bl <- map(hosts, function(a.host){
      if(fractions$blacklist.entirely[fractions$host==a.host]){
        # blacklist everything from this host
        if(verbose) cat("All tips from host ",a.host," blacklisted\n", sep="")
        
        return(na.omit(tip.names[hosts.for.tips==a.host]))
        
      } else {
        # blacklist smaller subgraphs
        
        
        duals.table <- ptree$duals.info
        if(!is.null(duals.table)){
          
          pat.rows <- duals.table[which(duals.table$host == a.host),]
          if(nrow(pat.rows)>0){
            if(verbose) cat("Tips from smaller subgraphs for host ",a.host," blacklisted\n", sep="")
            max.reads <- max(pat.rows$reads.in.subtree)
            smaller.tips <- pat.rows$tip.name[pat.rows$reads.in.subtree!=max.reads]
            
            return(na.omit(smaller.tips))
          }
        }
        return(vector())
      }
    }
    )
    
    new.bl <- new.bl %>% unlist()
    
    blacklists[[ptree$id]] <- new.bl
  }
  
  if(!is.null(summary.file)) {
    
    out.tbl <- fractions %>% 
      select(-denominator, -blacklist.entirely) %>%
      rename(count = numerator)
    
    if(verbose) cat("Writing dual summary file to ", summary.file, "\n", sep="")
    write_csv(out.tbl, file=summary.file, row.names=FALSE, quote=FALSE)
  }
  
  blacklists
}

