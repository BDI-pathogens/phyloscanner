args = commandArgs(trailingOnly=TRUE)

tree.file.root <- args[1]
norm.file.name <- args[2]
output.file.name <- args[3]

list.files.mod <- function(path = ".", pattern = NULL, all.files = FALSE,
                           full.names = FALSE, recursive = FALSE,
                           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE){
  
  original <- list.files(path, pattern, all.files, full.names, recursive, ignore.case, include.dirs, no..)
  return(sapply(original, function(x) if(substr(x, 1, 2)=="./") {substr(x, 3, nchar(x))} else {x}))
}


norm.table <- read.csv(norm.file.name, stringsAsFactors = F)

all.tree.files <- sort(list.files.mod(dirname(tree.file.root), pattern=paste(basename(tree.file.root),'.*\\.tree$',sep=''), full.names=TRUE))

suffixes <- substr(all.tree.files, nchar(tree.file.root) + 1, nchar(all.tree.files))
suffixes <- gsub('\\.tree','.csv',suffixes)

regex <- "^\\D*([0-9]+)_to_([0-9]+).*$"
starts <- sapply(suffixes, function(x) if(length(grep(regex, x))>0) as.numeric(sub(regex, "\\1", x)) else NA)
ends <- sapply(suffixes, function(x) if(length(grep(regex, x))>0) as.numeric(sub(regex, "\\2", x)) else NA)

window.df <- data.frame(file.name = basename(all.tree.files), start=starts, end=ends)

window.df$norm.constant <- sapply(1:nrow(window.df), function(x){
  rows <- which(norm.table$position >= window.df$start[x] & norm.table$position <= window.df$end[x])
  mean(norm.table$MEDIAN_PWD[rows])
} )

write.table(window.df[,c(1,4)], output.file.name, quote=F, col.names = F, row.names = F, sep=",")