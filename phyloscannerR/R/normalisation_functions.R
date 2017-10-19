# This function returns a normalisation constant for a genome window bounded by start and end.

#' @keywords internal
#' @export lookup.normalisation.for.tree

lookup.normalisation.for.tree <- function(tree.info, lookup.df, lookup.column = "NORM_CONST"){

  if(is.null(tree.info$window.coords)){
    stop(paste0("Tree ",tree.info$name, "has no window coordinates" ))
  }
  if(length(tree.info$window.coords)!=2){
    stop(paste0("Tree ",tree.info$name, "has malformed window coordinates" ))
  }
  return(.lookup.normalisation.for.coords(tree.info$window.coords$start, tree.info$window.coords$end, lookup.df, lookup.column))
}

#' @keywords internal
#' @export .lookup.normalisation.for.coords

.lookup.normalisation.for.coords <- function(start, end, lookup.df, lookup.column = "NORM_CONST"){
  if(!(lookup.column %in% colnames(lookup.df))){
    stop(paste0("No column called ",lookup.column," in lookup data frame."))
  }

  window.rows <- lookup.df[which(lookup.df$POSITION<=end & lookup.df$POSITION>=start),]
  
  return(colMeans(window.rows)[lookup.column]) 
}
