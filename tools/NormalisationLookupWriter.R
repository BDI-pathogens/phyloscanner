suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(data.table, quietly=TRUE, warn.conflicts=FALSE))
# 	Define arguments
arg_parser		<- ArgumentParser(description="Calculate normalising constants for each tree file from reference table.")
arg_parser$add_argument("tree.file.root", action="store", type="character", help="Start of tree file names.")
arg_parser$add_argument("norm.file.name", action="store", type="character", help="File name of reference table.")
arg_parser$add_argument("output.file.name", action="store", type="character", help="Output file name of csv file that has two columns, the base tree file name and the corresponding normalising constant.")
arg_parser$add_argument("norm.var", action="store", type="character", help="Column name in reference table that is to be used as normalising constant.")
arg_parser$add_argument("--standardize", action="store_true", default=FALSE, help="If true, the normalising constants are standardized so that the average on gag+pol equals 1. This way the normalised branch lengths are interpretable as typical distances on gag+pol")
arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")

# 	Parse arguments
args 			<- arg_parser$parse_args()
script.dir <- args$scriptdir

source(file.path(script.dir, "GeneralFunctions.R"))

tree.file.root 	<- args$tree.file.root
norm.file.name 	<- args$norm.file.name
output.file.name<- args$output.file.name
norm.var 		<- args$norm.var
norm.standardize<- args$standardize
verbose <- args$verbose
#	tree.file.root <- '/Users/Oliver/duke/tmp/pty_17-04-06-21-47-20/ptyr22_InWindow_'
#	norm.file.name <- '/Users/Oliver/git/phyloscan/data/hiv.hxb2.norm.constants.rda'
#	norm.var		<- "MEDIAN_PWD"

#	start script
w.regex	<- "^\\D*([0-9]+)_to_([0-9]+).*$"
#	load reference table
#	and define normalising constant
#	then define midpoints of windows
if (verbose) cat('Loading normalising constants reference file ', norm.file.name, "\n")
if(grepl('csv$',norm.file.name))
	norm.table	<- as.data.table(read.csv(norm.file.name, stringsAsFactors=FALSE))
if(grepl('rda$',norm.file.name)){
	tmp			<- load(norm.file.name)
	if(length(tmp)!=1)	stop("Expected one R data.table in file",norm.file.name)
	eval(parse(text=paste("norm.table<- ",tmp,sep='')))
	stopifnot( c('W_FROM','W_TO',norm.var)%in%colnames(norm.table) )
} 
setnames(norm.table, norm.var, 'NORM_CONST')
norm.table[, W_MID:= (W_FROM+W_TO)/2]
#	standardize to 1 on gag+pol ( prot + first part of RT in total 1300bp )
if(norm.standardize){
	#790 - 3385
	if (verbose) cat('Standardising normalising constants to 1 on the gag+pol (prot + first part of RT in total 1300bp pol) region\n')
	tmp		<- subset(norm.table, W_MID>=790L & W_MID<=3385L)
	stopifnot( nrow(tmp)>0 )	# norm.table must contain gag+pol region	
	tmp		<- tmp[, mean(NORM_CONST)]
	stopifnot( is.finite(tmp) )
	set(norm.table, NULL, 'NORM_CONST', norm.table[, NORM_CONST/tmp])
}
norm.table	<- subset(norm.table, select=c(W_MID, NORM_CONST))
#	guess window coordinates of tree files from file name
if (verbose) cat('Reading tree files starting with ', tree.file.root, "\n")
wdf		<- list.files.mod(dirname(tree.file.root), pattern=paste('^',basename(tree.file.root),'.*\\.tree$',sep=''), full.names=FALSE)
stopifnot( length(wdf)>0 )	# expect some tree files
wdf		<- data.table(F=sort(wdf))
wdf[, SUFFIX:= gsub('\\.tree','',substring(F, nchar(basename(tree.file.root)) + 1L))]
wdf[, W_FROM:= as.integer(gsub(w.regex, '\\1', SUFFIX))]
wdf[, W_TO:= as.integer(gsub(w.regex, '\\2', SUFFIX))]
setkey(wdf, W_FROM)
#	assign to each tree file 
#	the average normalising constant of those windows 
#	whose midpoint falls into W_FROM to W_TO
tmp		<- wdf[, {
			list( NORM_CONST=subset(norm.table, W_MID>=W_FROM & W_MID<=W_TO)[, mean(NORM_CONST)]	)
		}, by=c('W_FROM','W_TO')]
wdf		<- merge(wdf, tmp, by=c('W_FROM','W_TO'))
wdf		<- subset(wdf, select=c('F','NORM_CONST'))
#	write to file
cat('Writing normalising constants to file', output.file.name,'\n')
write.table(wdf, output.file.name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
