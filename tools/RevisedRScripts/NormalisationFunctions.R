# This function returns a normalisation constant for a genome window bounded by start and end.

lookup.normalisation.for.tree <- function(tree.info, lookup.df, lookup.column = "NORM_CONST"){
  if(is.null(tree.info$window.coords)){
    stop(paste0("Tree ",tree.info$name, "has no window coordinates" ))
  }
  if(length(tree.info$window.coords)!=2){
    stop(paste0("Tree ",tree.info$name, "has malformed window coordinates" ))
  }
  return(.lookup.normalisation.for.coords(tree.info$window.coords$start, tree.info$window.coords$end, lookup.df, lookup.column))
}

.lookup.normalisation.for.coords <- function(start, end, lookup.df, lookup.column = "NORM_CONST"){
  if(!(lookup.column %in% colnames(lookup.df))){
    stop(paste0("No column called ",lookup.column," in lookup data frame."))
  }
  middle <- (end + start)/2
  
  window.rows <- lookup.df[which(lookup.df$W_FROM<=middle & lookup.df$W_TO>=middle),]

  return(colMeans(window.rows)[lookup.column]) 
  
}
