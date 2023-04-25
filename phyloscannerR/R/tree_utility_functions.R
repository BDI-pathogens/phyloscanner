
# Get the tip number of this label

#' @keywords internal
#' @export get.tip.no

get.tip.no <- function(tree, name) {
  return(match(name, tree$tip.label))
}

# Get all tips associated with a host

#' @keywords internal
#' @export get.tips.for.host

get.tips.for.host <- function(tree, host.string, host.ids, blacklist) {
  node.numbers <- which(grepl(paste0('^',host.string),host.ids))
  
  node.numbers <-
    node.numbers[which(!(node.numbers %in% blacklist))]
  return(node.numbers)
}

#' @keywords internal
#' @export get.tips.for.sample

get.tips.for.sample <- function(tree, sample.string){
  return(which(grepl(paste0('^',sample.string),tree$tip.label)))
}

# Edge length by node number (probably in some library somewhere)

#' @keywords internal
#' @export get.edge.length

get.edge.length <- function(tree, node) {
  index <- which(tree$edge[,2] == node)
  return(tree$edge.length[index])
}

# Node is a tip (likewise)

#' @keywords internal
#' @export is.tip

is.tip <- function(tree, node) {
  return(node <= length(tree$tip.label))
}

# Node is root (likewise)

#' @keywords internal
#' @export is.root

is.root <- function(tree, node) {
  return(Ancestors(tree, node, type = "parent") == 0)
}

# Get the sequence of ancestors from one node to the other

#' @keywords internal
#' @export get.ancestral.sequence

get.ancestral.sequence <- function(tree, desc, anc) {
  if (!(anc %in% Ancestors(tree, desc, type = "all"))) {
    stop("anc is not an ancestor of desc at all")
  }
  current <- desc
  out <- desc
  
  while (current != anc) {
    current <- Ancestors(tree, current, type = "parent")
    
    out <- c(out, current)
    
  }
  return(out)
}

# Output the host that this tip belongs to. The vector host.ids 
# records hosts in the same order as the tips in the tree

#' @keywords internal
#' @export get.host.from.tip

get.host.from.tip <- function(tree, node, host.ids) {
  # is it a tip?
  if (!(is.tip(tree, node))) {
    stop("Not a tip")
  }
  host.id <- host.ids[node]
  return(host.id)
  
}

# Output the set of hosts whose tips are descended from a node 
# (or the host of the node itself if a tip)

#' @keywords internal
#' @export get.hosts.from.this.clade

get.hosts.from.this.clade <- function(tree, node, host.ids) {
  if (is.tip(tree,node)) {
    return(get.host.from.tip(tree, node, host.ids))
  }
  
  out <- vector()
  tips <- Descendants(tree, node, type = "tips")[[1]]
  for (tip in tips) {
    if (!(get.host.from.tip(tree,tip) %in% out)) {
      out <- c(out, get.host.from.tip(tree, tip, host.ids))
    }
  }
  return(out)
}

# Is desc _unambiguously_ a descendant of anc? I.e. is the MRCA node of desc 
# a descendant of the MRCA node of anc, but no tips of anc are descended from
# the MRCA node of desc?
# mrca.list is a list of MRCA nodes indexed by host
# tip.list is a list of vectors of tips indexed by host

#' @keywords internal
#' @export is.descendant.of

is.descendant.of <- function(tree, desc, anc, mrca.list, tip.list) {
  mrca.desc <- mrca.list[[desc]]
  mrca.anc <- mrca.list[[anc]]
  
  if (!(mrca.anc %in% Ancestors(tree, mrca.desc, type = "all"))) {
    return(FALSE)
  } else {
    for (tip in tip.list[[anc]]) {
      if (mrca.desc %in% Ancestors(tree, tip, type = "all")) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

# Is the intersection between vectors x and y empty?

#' @keywords internal
#' @export nonempty.intersection

nonempty.intersection <- function(x,y) {
  return(length(intersect(x,y)) > 0)
}

# Returns whether the tips corresponding to hosts pat.1 and pat.2 are intermingled
# in tree, i.e. if the MRCA of one set is both descended from (or equal to) the MRCA of 
# the other and an ancestor of at least one tip from the other
# mrca.list is a list of MRCA nodes indexed by host
# tip.list is a list of vectors of tips indexed by host

#' @keywords internal
#' @export are.intermingled

are.intermingled <-
  function(tree, pat.1, pat.2, mrca.list, tip.list) {
    mrca.1 <- mrca.list[[pat.1]]
    mrca.2 <- mrca.list[[pat.2]]
    
    if (mrca.1 == mrca.2) {
      return(TRUE)
    }
    
    if ((mrca.1 %in% Ancestors(tree, mrca.2, type = "all"))) {
      for (tip in tip.list[[pat.1]]) {
        if (mrca.2 %in% Ancestors(tree, tip, type = "all")) {
          return(TRUE)
        }
      }
    }
    
    if ((mrca.2 %in% Ancestors(tree, mrca.1, type = "all"))) {
      for (tip in tip.list[[pat.2]]) {
        if (mrca.1 %in% Ancestors(tree, tip, type = "all")) {
          return(TRUE)
        }
      }
    }
    
    return(FALSE)
  }

# Output the set of hosts whose mrca is this node

#' @keywords internal
#' @export get.hosts.with.these.mrcas

get.hosts.with.these.mrcas <-
  function(tree, node, host.mrcas) {
    out <- vector()
    mrca.vec <- sapply(host.mrcas, "[[", 1)
    
    numbers <- which(mrca.vec == node)
    
    return(names(numbers))
  }

# Get the distance between the MRCA nodes of these two hosts

#' @keywords internal
#' @export get.mrca.distance

get.mrca.distance <- function(tree, pat.1, pat.2, mrca.list) {
  mrca.1 <- mrca.list[[pat.1]]
  mrca.2 <- mrca.list[[pat.2]]
  
  return(total.length.between(tree, mrca.1, mrca.2))
  
}

# Get the distance between these two nodes

#' @keywords internal
#' @export total.length.between

total.length.between <- function(tree, node.1, node.2) {
  if (node.1 == node.2) {
    return(0)
  }
  total.length <- 0
  if (node.1 %in% Descendants(tree, node.2, type = "all")) {
    current.node <- node.1
    while (current.node != node.2) {
      total.length <- total.length + get.edge.length(tree, current.node)
      current.node <- Ancestors(tree, current.node, type = "parent")
    }
    return(total.length)
    
  } else if (node.2 %in% Descendants(tree, node.1, type = "all")) {
    current.node <- node.2
    while (current.node != node.1) {
      total.length <- total.length + get.edge.length(tree, current.node)
      current.node <- Ancestors(tree, current.node, type = "parent")
    }
    return(total.length)
    
  } else {
    joint.mrca <- mrca.phylo(tree, c(node.1, node.2))
    return(
      total.length.between(tree, node.1, joint.mrca) + total.length.between(tree, node.2, joint.mrca)
    )
  }
  
}

# mrca.phylo applied to a tip only will not return that tip. This function will.

#' @keywords internal
#' @export mrca.phylo.or.unique.tip

mrca.phylo.or.unique.tip <-
  function(tree, node, zero.length.tips.count = FALSE) {
    if (length(node) == 1) {
      if (!zero.length.tips.count | !is.tip(tree,node)) {
        return(node)
      } else {
        length <- get.edge.length(tree, node)
        while (length < 1E-5) {
          node <- Ancestors(tree, node, type = "parent")
          if (is.root(tree, node)) {
            break
          }
          length <- get.edge.length(tree, node)
          
        }
        return(node)
      }
    } else {
      mrca <- mrca.phylo(tree, node)
      
      return(mrca)
      
    }
  }

# This function returns TRUE wherever elements are the same, including NA's, and false 
# everywhere else.

#' @keywords internal
#' @export compareNA

compareNA <- function(v1,v2) {
  same <- (v1 == v2)  |  (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

# Last element in a vector

#' @keywords internal
#' @export get.last

get.last <- function(vec)
  return(vec[length(vec)])

# Get host id from tip label

#' @keywords internal
#' @export host.from.label

host.from.label <- function(label, regexp){
  if(length(grep(regexp, label)>0)) {
    return(sub(regexp, "\\1", label))
  } else {
    return(NA)
  }
}

# Get read count from label

#' @keywords internal
#' @export read.count.from.label

read.count.from.label <- function(label, regexp){
  if(length(grep(regexp, label)>0)) {
    return(as.integer(sub(regexp, "\\3", label)))
  } else {
    return(NA)
  }
}

# Starting at node 1, determine whether the path to node 2 is blocked by node 3
# Should work on both rooted and unrooted trees I _hope_.

#' @keywords internal
#' @export path.exists

path.exists <- function(tree, node.1, node.2, node.3, last.node = -1){
  #  cat(node.1," ")
  if(node.1 == node.2){
    return(TRUE)
  }
  if(node.1 == node.3){
    return(FALSE)
  }
  neighbours = vector()
  if(!is.tip(tree, node.1)){
    neighbours <- Children(tree, node.1)
  }
  if(length(Ancestors(tree, node.1, type="parent"))!=0){
    if(Ancestors(tree, node.1, type="parent")!=0){
      neighbours <- c(neighbours, Ancestors(tree, node.1, type="parent"))
    }
  }
  for(neighbour in neighbours){
    if(neighbour != last.node){
      if(path.exists(tree, neighbour, node.2, node.3, node.1)){
        return(TRUE)
        break
      }
    }
  }
  #you've been everywhere and you can't find a way through
  return(FALSE)
}

#' @keywords internal
#' @export tips.reachable

tips.reachable <- function(tree, node, blocking.edges, last.node = -1){
  out <- vector()
  if(is.tip(tree, node)){
    out <- c(out, node)
  }
  neighbours = neighbouring.nodes(tree, node)
  for(neighbour in neighbours){
    if(neighbour != last.node){
      edge.vec <- c(node, neighbour)
      if(is.na(row.match(edge.vec, t(as.matrix(blocking.edges)))) &  is.na(row.match(rev(edge.vec), t(as.matrix(blocking.edges))))){
        out <- c(out, tips.reachable(tree, neighbour, blocking.edges, node))
      }
    }
  }
  return(out)
}

#' @keywords internal
#' @export neighbouring.nodes

neighbouring.nodes <- function(tree, node){
  part.1 <- tree$edge[which(tree$edge[,2]==node), 1]
  part.2 <- tree$edge[which(tree$edge[,1]==node), 2]
  return(c(part.1, part.2))
}


# Drop a set of tips and return a vector which maps nodes from the full tree to the subtree

#' @keywords internal
#' @export drop.tip.get.map

drop.tip.get.map <- function(phy, tip){
  if(length(unique(tree$tip.label))!=length(tree$tip.label)){
    stop("This won't work if there are duplicate tip names")
  }
  reference <- vector()
  
  phy.2 <- drop.tip(phy, tip)
  
  for(new.label in seq(1, length(phy.2$tip.label))){
    reference[new.label] = which(phy$tip.label==phy.2$tip.label[new.label])
  }
  
  mrcas.1 <- mrca(phy)
  mrcas.2 <- mrca(phy.2)
  
  for(tip.1 in seq(1, length(phy.2$tip.label))){
    for(tip.2 in seq(1, length(phy.2$tip.label))){
      if(tip.1 < tip.2){
        new.mrca <- mrcas.2[tip.1,tip.2]
        old.mrca <- mrcas.1[reference[tip.1], reference[tip.2]]
        reference[new.mrca] = old.mrca
      }
    }
  }
  return(list(tree = phy.2, reference=reference))
}

#' @keywords internal
#' @export extract.subtrees.for.hosts

extract.subtrees.for.hosts <- function(tree, hosts, splits){
  labels.to.keep <- splits$tip.names[which(splits$orig.hosts %in% hosts)]
  
  pruned.tree <- drop.tip(tree, which(!(tree$tip.label %in% labels.to.keep)))
  
  return(pruned.tree)
}

#' @keywords internal
#' @export prop.internal.longer.than.root

prop.internal.longer.than.root <- function(tree, split, splits){
  tips <- splits$tip.names[which(splits$host.splits==split)]
  
  if(length(tips)==1){
    return(0)
  }
  
  split.root <- mrca.phylo.or.unique.tip(tree, which(tree$tip.label %in% tips))
  dist.to.root <- 0 
  
  current.node <- split.root
  while(!is.root(tree, current.node)){
    dist.to.root <- dist.to.root + get.edge.length(tree, split.root)
    current.node <- Ancestors(tree, current.node, type="parent")
  }
  
  just.clade <- drop.tip(tree, which(!(tree$tip.label %in% tips)))
  
  edges <- just.clade$edge.length
  
  return(sum(edges > dist.to.root)/length(edges))
}

#' @keywords internal
#' @export process.tree
#' @importFrom ape di2multi root

process.tree <- function(tree, root.name=NULL, m.thresh=-1, blacklist.for.pruning = vector(), normalisation.constant = 1) {
  
  if(m.thresh != -1){
    tree <- di2multi(tree, tol = m.thresh)
  }
  
  if(!is.null(root.name)){
    tree <- root(tree, outgroup = root.name, resolve.root = T)
  } else {
    # There are problems with a non-binary root.  Resolve it arbitrarily; most times a root should be specified.
    if(length(Children(tree, getRoot(tree)))>2){
      first.child <- Children(tree, getRoot(tree))[1]
      if(first.child <= length(tree$tip.label)){
        tree <- root(tree, outgroup=first.child, resolve.root = T)
      } else {
        tree <- root(tree, node=first.child, resolve.root = T)
      }
    }
  }
  
  if(length(blacklist.for.pruning) > 0){
    tree <- drop.tip(tree, blacklist.for.pruning)
  }
  
  tree$edge.length <- tree$edge.length/normalisation.constant
  
  return(tree)
}

#' @title Select for further analysis relationship classifications by read and tip counts 
#' @export select.windows.by.read.and.tip.count
#' @author Oliver Ratmann, Matthew Hall
#' @param ptrees A list of class \code{phyloscanner.trees} produced by \code{phyloscanner.analyse.trees}.
#' @param dwin A data frame produced by \code{classify.pairwise.relationships}.
#' @param tip.regex The regular expression used to identify host IDs in tip names
#' @param min.reads The minimum number of reads from a host in a window needed in order for that window to count in determining relationships involving that patient
#' @param min.tips The minimum number of tips from a host in a window needed in order for that window to count in determining relationships involving that patient
#' @param verbose Verbose output
#' @return A data frame with viral phylogenetic classifications of pairwise host relationships in each deep sequence phylogeny
#' @examples
#' \dontrun{
#' require(phyloscannerR)
#' #
#' #	continue Rakai example,
#' #	load phyloscanner output from 'phyloscanner.analyse.trees'
#' #	
#' file	<- system.file(file.path('extdata','ptyr192_phsc_analyse_trees_output.R'),package='phyloscannerR')
#' load(file)	#loads 'phsc', output from 'phyloscanner.analyse.trees'
#' #	use distance thresholds found in analysis of Rakai couples
#' close.threshold <- 0.025
#' distant.threshold <- 0.05	
#' #	use relationship types based on adjacency
#' #	this also considers linkage etc between individuals who have dual infections, recombinants etc 
#' #	..and thus may not have *all* their subgraphs adjacent to each other
#' relationship.types <- c('proximity.3.way',					
#' 		'close.and.adjacent',					
#' 		'close.and.adjacent.and.directed',					
#' 		'close.and.adjacent.and.ancestry.cat')	
#' dwin <- classify.pairwise.relationships(phsc, allow.mt=TRUE, close.threshold=close.threshold, distant.threshold=distant.threshold,relationship.types=relationship.types, verbose=TRUE)
#' tip.regex <- "^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$"
#' min.reads <- 30
#' min.tips <- 1
#' dwin <- select.windows.by.read.and.tip.count(phsc, dwin, tip.regex, min.reads, min.tips)
#' #
#' # 	end of Rakai example
#' #
#' }	 	
select.windows.by.read.and.tip.count <- function(ptrees, dwin, tip.regex, min.reads, min.tips, verbose=F)			
{
  host.tips.and.reads <- map(ptrees, function(x) phyloscannerR:::get.tip.and.read.counts(x, all.hosts.from.trees(ptrees), tip.regex, attr(ptrees, 'has.read.counts'), verbose = F))
  host.tips.and.reads <- bind_rows(host.tips.and.reads)
  
  if(verbose) cat('Merging tip and read counts...\n')
  
  dwin <- dwin %>% 
    inner_join(host.tips.and.reads, by=c("host.1"="host.id", "tree.id")) %>% 
    rename(tips.1 = tips, reads.1=reads)
  dwin <- dwin %>% 
    inner_join(host.tips.and.reads, by=c("host.2"="host.id", "tree.id")) %>% 
    rename(tips.2 = tips, reads.2=reads)	
  
  if(verbose) cat('Reducing transmission window stats to windows with at least',min.reads,'reads and at least',min.tips,'tips...\n')
  
  dwin	<- dwin %>% 
    filter(reads.1 >= min.reads & reads.2 >= min.reads & tips.1 >= min.tips & tips.2 >= min.tips)
  
  if(verbose) cat('Total number of windows with transmission assignments is ',nrow(dwin),'.\n', sep="")		
  
  dwin
}
