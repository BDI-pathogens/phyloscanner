# Collects a variety of statistics about a single patient in a single tree

#' @keywords internal
#' @export calc.subtree.stats
#' @importFrom ape unroot

calc.subtree.stats <- function(host.id, tree.id, tree, tips.for.hosts, splits.table, tip.regex, no.read.counts, verbose = F){
  
  if(verbose) cat("Calculating statistics for host ",host.id,"\n", sep="")
  
  relevant.reads <- splits.table %>% filter(host == host.id)
  
  subgraphs <- length(unique(relevant.reads$subgraph))
  
  all.tips <- tips.for.hosts[[host.id]]
  
  subtree.all <- NULL

  if(length(all.tips)>1){
  
    subtree.all <- drop.tip(tree, tip=tree$tip.label[!(tree$tip.label %in% all.tips)])
    
    # just in case pruning changes the tip order
    
    if(!no.read.counts){
      reads.per.tip <- sapply(subtree.all$tip.label, function(x) as.numeric(read.count.from.label(x, tip.regex)))
    } else {
      reads.per.tip <- rep(1, length(subtree.all$tip.label))
    }
    
    names(reads.per.tip) <- subtree.all$tip.label
    
    overall.rtt <- calc.mean.root.to.tip(subtree.all, reads.per.tip)
    
    
    if(length(subtree.all$tip.label)>2){
      unrooted.subtree <- unroot(subtree.all)
      max.branch.length <- max(unrooted.subtree$edge.length)
    } else {
      # ape won't unroot two-tip trees
      max.branch.length <- sum(subtree.all$edge.length)
    }
    
    pat.distances <- cophenetic(subtree.all)
    max.pat.distance <- max(pat.distances)
    
    global.mean.pat.distance <- mean(pat.distances[upper.tri(pat.distances)])
    
    if(subgraphs==1){
      subgraph.mean.pat.distance <- global.mean.pat.distance
      largest.rtt <- overall.rtt
      
    } else {
      
      splits <- unique(relevant.reads$subgraph)
      
      reads.per.split <- sapply(splits, function(x) sum(splits.table$reads[which(splits.table$subgraph==x)] ) )
      
      winner <- splits[which(reads.per.split==max(reads.per.split))]
      
      if(length(winner)>1){
        #Rare but possible. Pick an arbitrary winner
        winner <- winner[1]
      }
      
      winner.tips <- relevant.reads$tip[which(relevant.reads$subgraph==winner)]
      if(length(winner.tips)==1){
        largest.rtt <- 0
        subgraph.mean.pat.distance <- NA
      } else {
        subtree <- drop.tip(tree, tip=tree$tip.label[!(tree$tip.label %in% winner.tips)])
        
        # just in case pruning changes the tip order
        
        if(!no.read.counts){
          reads.per.tip <- sapply(subtree$tip.label, function(x) as.numeric(read.count.from.label(x, tip.regex)))
        } else {
          reads.per.tip <- rep(1, length(subtree$tip.label))
        }
        
        names(reads.per.tip) <- subtree$tip.label
        
        largest.rtt <- calc.mean.root.to.tip(subtree, reads.per.tip)
        pat.distances <- cophenetic(subtree)
        subgraph.mean.pat.distance <- mean(pat.distances[upper.tri(pat.distances)])
      }
    }
    
  } else {
    
    # A patient with zero or one tips has no branches and patristic distances are meaningless.
    
    max.branch.length <- NA
    max.pat.distance <- NA
    global.mean.pat.distance <- NA
    subgraph.mean.pat.distance <- NA
    
    if(length(all.tips)==1){
      # A patient with one tip does have a root
      overall.rtt <- 0
      largest.rtt <- 0
    } else {
      # An absent patient has no root or tips
      overall.rtt <- NA
      largest.rtt <- NA
    }
  }
  
  return(list(overall.rtt = overall.rtt, largest.rtt = largest.rtt, max.branch.length = max.branch.length, max.pat.distance = max.pat.distance, 
              global.mean.pat.distance=global.mean.pat.distance, subgraph.mean.pat.distance = subgraph.mean.pat.distance))
}

# Calculates tip and read counts for all patients in a given window

#' @keywords internal
#' @export get.tip.and.read.counts

get.tip.and.read.counts <- function(ptree, hosts, tip.regex, has.read.counts, verbose = F){
  
  tree.id <- ptree$id
  
  if(verbose) cat("Calculating tip and read counts for tree ID ",ptree$id,"\n",sep="")
  
  tree <- ptree$tree
  blacklist <- ptree$blacklist
  
  # A vector of patients for each tip
  
  hosts.for.tips <- sapply(tree$tip.label, function(x) host.from.label(x, tip.regex))
  hosts.for.tips[blacklist] <- NA
  
  # If no patients from the input file are actually here, give a warning
  
  hosts.present <- intersect(hosts, unique(hosts.for.tips))
  if(length(hosts.present)==0){
    warning(paste("No listed hosts appear in tree ",tree.id,"\n",sep=""))
  }
  
  # A list of tips for each patient 
  
  tips.for.hosts <- map(setNames(hosts, hosts), function(x) tree$tip.label[which(hosts.for.tips==x)])
  
  # Make the tibble
  
  window.table <- tibble(host.id=hosts)
  window.table <- window.table %>% mutate(tree.id = tree.id) 
  window.table <- window.table %>% mutate(tips = map_int(host.id, function(x)  (length(tips.for.hosts[[x]]))))
  window.table <- window.table %>% mutate(reads = map_int(host.id, function(x){
    if(length(tips.for.hosts[[x]])==0){
      return(0L)
    } else {
      if(has.read.counts){
        return(sum(map_int(tips.for.hosts[[x]], function(y) read.count.from.label(y, tip.regex))))
      } else {
        return(as.integer(length(tips.for.hosts[[x]])))
      }
    }
  })) 
  
  window.table
}

#' @keywords internal
#' @export calc.alignment.stats.in.window
#' @importFrom pegas nuc.div

calc.alignment.stats.in.window <- function(ptree, hosts, tip.regex, has.read.counts, verbose = FALSE){
  if(verbose) cat("Calculating alignment host statistics for tree ID ",ptree$id,"\n",sep="")
  
  id <- ptree$id
  
  alignment <- ptree$alignment
  
  if(is.null(alignment)){
    stop("No alignment found.")
  }
  
  blacklist <- ptree$blacklist
  
  alignment <- alignment[which(!(rownames(alignment) %in% blacklist)),]
  
  hosts.for.tips <- sapply(rownames(alignment), function(x) host.from.label(x, tip.regex))
  
  hosts.present <- intersect(hosts, unique(hosts.for.tips))
  
  if(length(hosts.present)==0){
    warning(paste("No listed hosts appear in tree ",ptree$id,"\n",sep=""))
  }
  
  tips.for.hosts <- lapply(setNames(hosts, hosts), function(x) rownames(alignment)[which(hosts.for.tips==x)])
  
  if(has.read.counts){
    multiplicities <- map_dbl(rownames(alignment), function(x) as.numeric(read.count.from.label(x, tip.regex)))
  } else {
    multiplicities <- rep(1, nrow(alignment))
  }
    
  duplicated.alignment <- alignment
  
  for(i in 1:nrow(alignment)){
    if(!is.na(multiplicities[i]) & multiplicities[i] > 1){
      for(j in 2:multiplicities[i]){
        duplicated.alignment <- rbind(duplicated.alignment, alignment[i,])
      }
    }
  }
  
  window.table <- tibble(host.id = hosts, 
                         tree.id = id, 
                         nuc.div.raw = map_dbl(hosts, function(x){
                           sub.alignment <- alignment[which(rownames(alignment) %in% tips.for.hosts[[x]]),]
                           nuc.div(sub.alignment)
                         }),
                         nuc.div.weighted = map_dbl(hosts, function(x){
                           sub.alignment <- duplicated.alignment[which(rownames(duplicated.alignment) %in% tips.for.hosts[[x]]),]
                           nuc.div(sub.alignment)
                         }),
                         cumul.maf.raw = map_dbl(hosts, function(x){
                           sub.alignment <- alignment[which(rownames(alignment) %in% tips.for.hosts[[x]]),]
                           mean(map_dbl(1:ncol(sub.alignment), function(y) maf(sub.alignment[,y])) , na.rm = T)
                         }),
                         cumul.maf.weighted = map_dbl(hosts, function(x){
                           sub.alignment <- duplicated.alignment[which(rownames(duplicated.alignment) %in% tips.for.hosts[[x]]),]
                           mean(map_dbl(1:ncol(sub.alignment), function(y) maf(sub.alignment[,y])) , na.rm = T)
                         })
                         )
  window.table
}

# Minor allele frequency in one column

#' @keywords internal
#' @export maf

maf <- function(col.vec){
  col.vec <- as.character(col.vec)
  # just in case!
  
  col.vec[which(toupper(col.vec) == "U")] <- "T"
  
  
  nonmissing <- col.vec[which(col.vec != "-")]
  if(length(nonmissing)==0){
    return(NA)
  }
  
  denominator <- length(nonmissing)
  
  nucs <- c("A", "C", "G", "T")
  
  freqs <- map_dbl(nucs, function(x){
    sum(map_dbl(nonmissing, function(y){
      sum(nucleotide.probability(x)*nucleotide.probability(y))
    }))
  })
  
  if(sum(freqs) != denominator){
    stop("Your maths is wrong.")
  }
  return(1-(freqs[which(freqs == max(freqs))][1]/denominator))
}


# Takes a (potentially ambiguous) nucleotide code and returns a vector of probabilities of being A, C, G or T(/U) in that order

#' @keywords internal
#' @export nucleotide.probability
#' @importFrom glue glue

nucleotide.probability <- function(code){
  if(!(toupper(code) %in% c("A","C","G","T","U","R","Y","S","W","K","M","B","V","D","H","N","?"))){
    stop(glue("Unknown nucleotide code: {code}"))
  }
  
  switch(
    toupper(code),
    "A" = c(1,0,0,0),
    "C" = c(0,1,0,0),
    "G" = c(0,0,1,0),
    "T" = c(0,0,0,1),
    "U" = c(0,0,0,1),
    "R" = c(0.5,0,0.5,0),
    "Y" = c(0,0.5,0,0.5),
    "S" = c(0,0.5,0.5,0),
    "W" = c(0.5,0,0,0.5),
    "K" = c(0,0,0.5,0.5),
    "M" = c(0.5,0.5,0,0),
    "B" = c(0,1/3,1/3,1/3),
    "V" = c(1/3,1/3,1/3,0),
    "D" = c(1/3,0,1/3,1/3),
    "H" = c(1/3,1/3,0,1/3),
    "N" = c(0.25,0.25,0.25,0.25),
    "?" = c(0.25,0.25,0.25,0.25)
  )
}

#' @keywords internal
#' @export calc.phylo.stats.in.window

calc.phylo.stats.in.window <- function(ptree, hosts, tip.regex, has.read.counts, verbose = FALSE){

  if(verbose) cat("Calculating phylogenetic host statistics for tree ID ",ptree$id,"\n",sep="")
  
  id <- ptree$id
  
  # todo this needs to go 
  
  if(is.null(id)){
    id <- ptree$suffix
  }
  
  tree <- ptree$tree
  blacklist <- ptree$blacklist
  
  # Get the splits
  
  splits.table <- ptree$splits.table
  
  # Find the clades
  
  clade.mrcas.by.host <- ptree$clade.mrcas.by.host
  clades.by.host <- ptree$clades.by.host
  
  # A vector of patients for each tip
  
  hosts.for.tips <- sapply(tree$tip.label, function(x) host.from.label(x, tip.regex))
  hosts.for.tips[blacklist] <- NA
  
  # If no patients from the input file are actually here, give a warning
  
  hosts.present <- intersect(hosts, unique(hosts.for.tips))
  
  if(length(hosts.present)==0){
    warning(paste("No listed hosts appear in tree ",ptree$id,"\n",sep=""))
  }
  
  # A list of tips for each host 
  
  tips.for.hosts <- lapply(setNames(hosts, hosts), function(x) tree$tip.label[which(hosts.for.tips==x)])
  
  # Make the tibble
  
  window.table <- tibble(host.id = hosts, 
                         tree.id = id, 
                         xcoord = ptree$xcoord, 
                         tips = sapply(hosts, function(x) as.numeric(length(tips.for.hosts[[x]]))),
                         reads = sapply(hosts, function(x){
                           if(length(tips.for.hosts[[x]])==0){
                             return(0)
                           } else {
                             if(has.read.counts){
                               return(sum(sapply(tips.for.hosts[[x]], function(y) as.numeric(read.count.from.label(y, tip.regex)))))
                             } else {
                               return(as.numeric(length(tips.for.hosts[[x]])))
                             }
                           }
                         } 
                         ),
                         subgraphs =  sapply(hosts, function(x) length(unique(splits.table[which(splits.table$host==x),]$subgraph))),
                         clades = sapply(hosts, function(x) length(clades.by.host[[x]]))
  )
  
  new.cols <- map_dfr(hosts, function(x) calc.subtree.stats(x, id, tree, tips.for.hosts, splits.table, tip.regex, !has.read.counts, verbose))
  
  normalised.new.cols <- new.cols %>% transmute_all(funs(./ptree$normalisation.constant))  
  
  colnames(normalised.new.cols) <- paste0("normalised.", colnames(new.cols))
  
  window.table <- bind_cols(window.table, new.cols) 
  window.table <- bind_cols(window.table, normalised.new.cols) 
  
  recomb.file.name <- ptree$recombination.file.name
  
  
  
  if (!is.null(recomb.file.name)) {
    
    recomb.df <- read_csv(recomb.file.name)
    col.names <- colnames(recomb.df)
    for (expected.col.name in c("Bam file","Recombination metric")) {
      if (! expected.col.name %in% col.names) {
        stop("Expected column name ", expected.col.name, " missing from ",
             recomb.file.name, ". Quitting.")
      }
    }
    missing.IDs <- setdiff(hosts, recomb.df$`Bam file`)
    extra.IDs <- setdiff(recomb.df$`Bam file`, hosts)
    if (length(missing.IDs) > 0) {
      stop(paste0("Error: the following bam file IDs are missing from ",
                  recomb.file.name, ": ", paste(missing.IDs, collapse = " "), "\nQuitting."))
    }
    if (length(extra.IDs) > 0) {
      stop(paste0("Error: the following bam file IDs were unexpected in ",
                  recomb.file.name, ": ", paste(extra.IDs, collapse = " "), "\nQuitting."))
    }
    recomb.df <- recomb.df %>% 
      select(`Bam file`, `Recombination metric`) %>% 
      rename(host.id = `Bam file`, recombination.metric = `Recombination metric`)

    window.table <- window.table %>% inner_join(recomb.df, by = "host.id")
  }
  
  # If you did dual detection, add that in
  
  if(!is.null(ptree$dual.detection.splits)){
    # Q: is there a tidyverse match?
    
    window.table$solo.dual.count <- ptree$dual.detection.splits$count[match(window.table$host.id, ptree$dual.detection.splits$host)]
    
  }
  

  window.table
}

# Get the proportion of reads from a given patient in each of its subgraphs in a single tree

#' @keywords internal
#' @export get.read.proportions

get.read.proportions <- function(host.id, tree.id, splits.table){
  
  this.pat.splits <- splits.table[which(splits.table$host==host.id),]
  
  if(nrow(this.pat.splits)>0){
    this.pat.reads.by.split <- aggregate(this.pat.splits$reads, by=list(Category=this.pat.splits$subgraph), sum)
    
    this.pat.reads.by.split <- this.pat.reads.by.split[order(this.pat.reads.by.split$x, decreasing = T),]
    return(this.pat.reads.by.split$x/sum(this.pat.reads.by.split$x))
    
  } else {
    return(NA)
  }
}

