# list.files seems to behave differently on different systems

#' @keywords internal
#' @export list.files.mod

list.files.mod <- function(path = ".", pattern = NULL, all.files = FALSE,
                           full.names = FALSE, recursive = FALSE,
                           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE){
  
  original <- list.files(path, pattern, all.files, full.names, recursive, ignore.case, include.dirs, no..)
  return(sapply(original, function(x) if(substr(x, 1, 2)=="./") {substr(x, 3, nchar(x))} else {x}))
}

#' @keywords internal
#' @export get.suffix

get.suffix <- function(file.name, prefix, extension){
  
  file.name <- basename(file.name)
  prefix <- basename(prefix)

  substr(file.name, nchar(prefix)+1, nchar(file.name)-nchar(extension)-1)
}

#' @keywords internal
#' @export get.window.coords

get.window.coords <- function(string, regex = "^\\D*([0-9]+)_to_([0-9]+).*$"){
  start <- if(length(grep(regex, string))>0) as.numeric(sub(regex, "\\1", string)) else NA
  end   <- if(length(grep(regex, string))>0) as.numeric(sub(regex, "\\2", string)) else NA

  if(any(is.na(start)) | any(is.na(end))) {
    stop(paste0("ERROR: cannot determine window coordinates"))
  }
  
  return(list(start=start, end = end))
}

#' @keywords internal
#' @export get.window.coords.string

get.window.coords.string <- function(string, regex = "^\\D*([0-9]+)_to_([0-9]+).*$"){
  wc <- get.window.coords(string, regex)
  return(paste0(wc$start, "_to_", wc$end))
}

#' @keywords internal
#' @export all.hosts.from.trees

all.hosts.from.trees <- function(ptrees){
  hosts <- ptrees %>% map("hosts.for.tips")
  hosts <- unique(unlist(hosts))
  hosts <- hosts[!is.na(hosts)]
  hosts <- hosts[order(hosts)]
  hosts
}
