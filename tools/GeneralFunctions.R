# list.files seems to behave differently on different systems

list.files.mod <- function(path = ".", pattern = NULL, all.files = FALSE,
                           full.names = FALSE, recursive = FALSE,
                           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE){
  
  original <- list.files(path, pattern, all.files, full.names, recursive, ignore.case, include.dirs, no..)
  return(sapply(original, function(x) if(substr(x, 1, 2)=="./") {substr(x, 3, nchar(x))} else {x}))
}


get.suffix <- function(file.name, prefix, extension){
  
  file.name <- basename(file.name)
  prefix <- basename(prefix)

  substr(file.name, nchar(prefix)+1, nchar(file.name)-nchar(extension)-1)
}

prepare.tree <- function(file.name, root.name = NULL, multi.threshold = NULL, normalisation.constant = 1){
  tree <- read.tree(file.name)
  tree <- unroot(tree)
  if(!is.null(multi.threshold)){
    if(is.na(multi.threshold)){
      multi.threshold <- 1.0001*min(tree$edge.length)
    }
    tree <- di2multi(tree, tol = multi.threshold)
  }
  if(!is.null(root.name)){
    tree <- root(tree, outgroup = which(tree$tip.label==root.name), resolve.root = T)
  }
  tree$edge.length <- tree$edge.length/normalisation.constant
  
  return(tree)
  
}


get.window.coords <- function(string){

  regex <- "^\\D*([0-9]+)_to_([0-9]+).*$"
  
  start <- if(length(grep(regex, string))>0) as.numeric(sub(regex, "\\1", string)) else NA
  end <- if(length(grep(regex, string))>0) as.numeric(sub(regex, "\\2", string)) else NA
  
  if(any(is.na(start)) | any(is.na(end))) {
    stop(paste0("ERROR: cannot determine window coordinates"))
  }
  
  return(list(start=start, end = end))
}
