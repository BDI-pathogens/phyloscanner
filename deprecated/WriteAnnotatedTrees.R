# This is borrowed from ape with minor modifications

library(ape)

write.ann.tree <-
  function(phy, file = "", annotations = NULL, append = FALSE, digits = 10, tree.names = FALSE)
  {
    if (!(inherits(phy, c("phylo", "multiPhylo"))))
      stop("object \"phy\" has no trees")
    
    if (inherits(phy, "phylo")) phy <- c(phy)
    N <- length(phy)
    res <- character(N)
    
    if (is.logical(tree.names)) {
      if (tree.names) {
        tree.names <-
          if (is.null(names(phy))) character(N)
        else names(phy)
      } else tree.names <- character(N)
    }
    
    for (i in 1:N)
      res[i] <- .write.ann.tree2(phy[[i]], annotations, digits = digits,
                             tree.prefix = tree.names[i])
    
    if (file == "") return(res)
    else cat(res, file = file, append = append, sep = "\n")
  }

.write.ann.tree2 <- function(phy, annotations = NULL, digits = 10, tree.prefix = "")
{
  brl <- !is.null(phy$edge.length)
  nodelab <- !is.null(phy$node.label)
  phy$tip.label <- checkLabel(phy$tip.label)
  if (nodelab) phy$node.label <- checkLabel(phy$node.label)
  f.d <- paste("%.", digits, "g", sep = "")
  cp <- function(x){
    STRING[k] <<- x
    k <<- k + 1
  }
  add.internal <- function(i) {
    cp("(")
    desc <- kids[[i]]
    for (j in desc) {
      if (j > n) add.internal(j)
      else add.terminal(ind[j])
      if (j != desc[length(desc)]) cp(",")
    }
    cp(")")
    if(!is.null(annotations)){
      ann.string <- process.annotations(phy, i, annotations)
      if(!is.na(ann.string)){
        cp(ann.string)
      }
    }
    if (nodelab && i > n) cp(phy$node.label[i - n]) # fixed by Naim Matasci (2010-12-07)
    if (brl) {
      cp(":")
      cp(sprintf(f.d, phy$edge.length[ind[i]]))
    }
  }
  add.terminal <- function(i) {
    cp(phy$tip.label[phy$edge[i, 2]])
    if(!is.null(annotations)){
      ann.string <- process.annotations(phy, which(ind==i), annotations)
      if(!is.na(ann.string)){
        cp(ann.string)
      }
    }
    if (brl) {
      cp(":")
      cp(sprintf(f.d, phy$edge.length[i]))
    }
  }
  
  n <- length(phy$tip.label)
  
  ## borrowed from phangorn:
  parent <- phy$edge[, 1]
  children <- phy$edge[, 2]
  kids <- vector("list", n + phy$Nnode)
  for (i in 1:length(parent))
    kids[[parent[i]]] <- c(kids[[parent[i]]], children[i])
  
  ind <- match(1:max(phy$edge), phy$edge[, 2])
  
  LS <- 4*n + 5
  if (brl) LS <- LS + 4*n
  if (nodelab)  LS <- LS + n
  STRING <- character(LS)
  k <- 1
  cp(tree.prefix)
  cp("(")
  getRoot <- function(phy)
    phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]
  root <- getRoot(phy) # replaced n+1 with root - root has not be n+1
  desc <- kids[[root]]
  for (j in desc) {
    if (j > n) add.internal(j)
    else add.terminal(ind[j])
    if (j != desc[length(desc)]) cp(",")
  }
  
  if (is.null(phy$root.edge)) {
    cp(")")
    if (nodelab) cp(phy$node.label[1])
    cp(";")
  }
  else {
    cp(")")
    if (nodelab) cp(phy$node.label[1])
    cp(":")
    cp(sprintf(f.d, phy$root.edge))
    cp(";")
  }
  paste(STRING, collapse = "")
}

process.annotations <- function(tree, index, ann.names){
  any.not.na <- F
  value.list <- list()
  for(name in ann.names){
    tree.att <- attr(tree, name)
    value.list[[name]] <- tree.att[index]
    if(!is.na(value.list[[name]])){
      any.not.na <- T
    }
  }
 
  if(any.not.na){
    string <- "[&"
    first <- T
    for(item in ann.names) {
      if(first){
        first <- F
      } else {
        string <- paste(string, ",", sep="")
      }
      string <- paste(string, item, "=", sep="")
      if(length(value.list[[item]])==1){
        string <- paste(string, quote.if.char(value.list[[item]]), sep="")
      } else {
        string <- paste(string, "{", paste(quote.if.char(value.list[[item]]), collapse=","), "}", sep="")
      }
    }
    
    string <- paste(string, "]", sep="")
    return(string)
  }
  return(NA)
}

quote.if.char <- function(x){
  if(is.numeric(x)){
    return(x)
  } else {
    return(paste("\"", x, "\"", sep=""))
  }
}

write.ann.nexus <- function(..., annotations = NULL, file = "", translate = TRUE)
{
  obj <- ape:::.getTreesFromDotdotdot(...)
  ntree <- length(obj)
  cat("#NEXUS\n", file = file)
  cat(paste("[R-package APE, ", date(), "]\n\n", sep = ""),
      file = file, append = TRUE)
  
  N <- length(obj[[1]]$tip.label)
  
  cat("Begin TAXA;\n", file = file, append = TRUE)
  cat(paste("\tDIMENSIONS NTAX = ", N, ";\n", sep = ""),
      file = file, append = TRUE)
  cat("\tTAXLABELS\n", file = file, append = TRUE)
  cat(paste("\t\t'", obj[[1]]$tip.label, "'", sep = ""),
      sep = "\n", file = file, append = TRUE)
  cat("\t;\n", file = file, append = TRUE)
  cat("End;\n", file = file, append = TRUE)
  
  cat("Begin trees;\n", file = file, append = TRUE)
  if (translate) {
    cat("\tTRANSLATE\n", file = file, append = TRUE)
    obj <- .compressTipLabel(obj)
    X <- paste("\t\t", 1:N, " '", attr(obj, "TipLabel"), "',", sep = "")
    ## We remove the last comma:
    X[length(X)] <- gsub(",", "", X[length(X)])
    cat(X, file = file, append = TRUE, sep = "\n")
    cat("\t;\n", file = file, append = TRUE)
    class(obj) <- NULL
    for (i in 1:ntree)
      obj[[i]]$tip.label <- as.character(1:N)
  } else {
    if (is.null(attr(obj, "TipLabel"))) {
      for (i in 1:ntree)
        obj[[i]]$tip.label <- checkLabel(obj[[i]]$tip.label)
    } else {
      attr(obj, "TipLabel") <- checkLabel(attr(obj, "TipLabel"))
      obj <- .uncompressTipLabel(obj)
    }
  }
  
  title <- names(obj)
  if (is.null(title))
    title <- rep("UNTITLED", ntree)
  else {
    if (any(s <- title == "")) title[s] <- "UNTITLED"
  }
  
  for (i in 1:ntree) {
    if (class(obj[[i]]) != "phylo") next
    root.tag <- if (is.rooted(obj[[i]])) "= [&R] " else "= [&U] "
    cat("\tTree *", title[i], root.tag, file = file, append = TRUE)
    cat(write.ann.tree(obj[[i]], annotations = annotations, file = ""),
        "\n", sep = "", file = file, append = TRUE)
  }
  cat("End;\n", file = file, append = TRUE)
}
