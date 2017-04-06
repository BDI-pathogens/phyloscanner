args = commandArgs(trailingOnly=TRUE)

#	tree.file.root <- '/Users/Oliver/duke/tmp/pty_17-04-06-11-39-28/ptyr19_InWindow_'
#	norm.file.name <- '/Users/Oliver/Library/R/3.3/library/phyloscan/data/hiv.hxb2.norm.constants.rda'
tree.file.root <- args[1]
norm.file.name <- args[2]
output.file.name <- args[3]
norm.var <- args[4]				#column name to be used in normalisation reference file

require(data.table)
w.regex	<- "^\\D*([0-9]+)_to_([0-9]+).*$"

list.files.mod <- function(path = ".", pattern = NULL, all.files = FALSE,
                           full.names = FALSE, recursive = FALSE,
                           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE){
  
  original <- list.files(path, pattern, all.files, full.names, recursive, ignore.case, include.dirs, no..)
  return(sapply(original, function(x) if(substr(x, 1, 2)=="./") {substr(x, 3, nchar(x))} else {x}))
}

cat('\nload normalising constants reference file', norm.file.name)
if(grepl('csv$',norm.file.name))
	norm.table	<- as.data.table(read.csv(norm.file.name, stringsAsFactors = F))
if(grepl('rda$',norm.file.name))
{
	tmp			<- load(norm.file.name)
	if(length(tmp)!=1)	stop("Expected one R data.table in file",norm.file.name)
	eval(parse(text=paste("norm.table<- ",tmp,sep='')))
	stopifnot( c('W_FROM','W_TO',norm.var)%in%colnames(norm.table) )
} 

setnames(norm.table, norm.var, 'NORM_CONST')
norm.table[, W_MID:= (W_FROM+W_TO)/2]
norm.table	<- subset(norm.table, select=c(W_MID, NORM_CONST))

wdf		<- data.table(F=sort(list.files.mod(dirname(tree.file.root), pattern=paste('^',basename(tree.file.root),'.*\\.tree$',sep=''), full.names=TRUE)))
wdf[, SUFFIX:= gsub('\\.tree','',substring(F, nchar(tree.file.root) + 1L))]
wdf[, W_FROM:= as.integer(gsub(w.regex, '\\1', SUFFIX))]
wdf[, W_TO:= as.integer(gsub(w.regex, '\\2', SUFFIX))]
setkey(wdf, W_FROM)
tmp		<- wdf[, {
			list( NORM_CONST=subset(norm.table, W_MID>=W_FROM & W_MID<=W_TO)[, mean(NORM_CONST)]	)
		}, by=c('W_FROM','W_TO')]
wdf		<- merge(wdf, tmp, by=c('W_FROM','W_TO'))
wdf		<- subset(wdf, select=c('F','SUFFIX','NORM_CONST'))

cat('\nwrite normalising constants to file', output.file.name)
write.table(window.df[,c(1,4)], output.file.name, quote=F, col.names = F, row.names = F, sep=",")
