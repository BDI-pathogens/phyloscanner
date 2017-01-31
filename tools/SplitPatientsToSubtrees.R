command.line <- T
list.of.packages <- c("phangorn", "argparse", "phytools", "ggplot2", "gdata", "mvtnorm", "expm")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

if(!("ggtree" %in% installed.packages()[,"Package"])){
  source("https://bioconductor.org/biocLite.R")
  biocLite("ggtree")
}

if(command.line){
  suppressMessages(require(argparse, quietly=TRUE, warn.conflicts=FALSE))
  
  arg_parser = ArgumentParser(description="Split the tips from each patient up, according to one of several schemes, such that the tree can subsequently be partitioned into connected subtrees, one per split.")  
  arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
  arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
  arg_parser$add_argument("-z", "--zeroLengthTipsCount", action="store_true", default=FALSE, help="If present, a zero length terminal branch associates the patient at the tip with its parent node, interrupting any inferred transmission pathway between another pair of hosts that goes through the node.")
  arg_parser$add_argument("-s", "--splitsRule", action="store", default="r", help="The rules by which the sets of patients are split into groups in order to ensure that all groups can be members of connected subtrees without causing conflicts. Currently available: c=conservative, r=Romero-Severson (default).")
  arg_parser$add_argument("-b", "--blacklist", action="store", help="A .csv file listing tips to ignore. Alternatively, a base name that identifies blacklist files. Blacklist files are assumed to end in .csv.")
  arg_parser$add_argument("-k", "--kParam", action="store", default=35, help="The k parameter in the cost matrix for Sankhoff reconstruction (see documentation)")
  arg_parser$add_argument("-i", "--inputFile", metavar="inputTreeFileName", help="Tree file name. Alternatively, a base name that identifies a group of tree file names can be specified. Tree files are assumed to end in .tree.")  
  arg_parser$add_argument("-t", "--breakTiesUnsampled", action="store_true", default=FALSE)
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
  arg_parser$add_argument("-OD", "--outputdir", action="store", help="Full path of the directory for output; if absent, current working directory")
  arg_parser$add_argument("-OF", "--outputfileid", action="store", help="A string identifying output files.")
  arg_parser$add_argument("-pw", "--pdfwidth", action="store", default=100, help="Width of tree pdf in inches.")
  arg_parser$add_argument("-ph", "--pdfrelheight", action="store", default=0.15, help="Relative height of tree pdf.")
  arg_parser$add_argument("-db", "--debug", action="store_true", default=FALSE, help="Debugging mode.")
  
  
  
  args <- arg_parser$parse_args()

  script.dir <- args$scriptdir
  zero.length.tips.count <- args$zeroLengthTipsCount
  break.ties.unsampled <- args$breakTiesUnsampled
  tree.file.names <- args$inputFile
  output.dir <- args$outputdir
  if(is.null(output.dir)){
    output.dir <- getwd()
  }
  sankhoff.k <- as.numeric(args$kParam)
  out.identifier <- args$outputfileid
  blacklist.files <- args$blacklist
  root.name <- args$outgroupName
  tip.regex <- args$tipRegex
  mode <- args$splitsRule
  pdf.hm <- as.numeric(args$pdfrelheight)
  pdf.w <- as.numeric(args$pdfwidth)
  debug <- args$debug
  
  if(!(mode %in% c("c", "r", "f", "s"))){
    stop(paste("Unknown split classifier: ", mode, "\n", sep=""))
  }
#  if(debug)
#	  cat(paste(script.dir, zero.length.tips.count, tree.file.names, out.identifier, blacklist.files, root.name, tip.regex, mode, pdf.hm, pdf.w, sep='\n'))
  
  
} else {
  # script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools/"
  # mode <- "r"
  # zero.length.tips.count <- F
  # sankhoff.k <- 0
  # pdf.w = 10
  # pdf.hm = 0.2
  # 
  # # BEEHIVE example
  # 
  # setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517_clean/")
  # output.dir <- getwd()
  # tree.file.names <- "RAxML_bestTree.InWindow_1550_to_1900.tree"
  # blacklist.files <- "PatientBlacklist_InWindow_1550_to_1900.csv"
  # out.identifier <- "test_s"
  # mode <- "s"
  # root.name <- "C.BW.00.00BW07621.AF443088"
  # tip.regex <- "^(.*)-[0-9].*_read_([0-9]+)_count_([0-9]+)$"
  # sankhoff.k <- 25
  # break.ties.unsampled <- TRUE
  # 
  # # Rakai example 
  # 
  # setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/RakaiSeqs/")
  # output.dir <- getwd()
  # tree.file.names <- "RAxML_bestTree.InWindow_550_to_900.tree"
  # blacklist.files <- "FullBlacklist_InWindow_550_to_900.csv"
  # out.identifier <- "test_pytr1"
  # root.name <- "B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
  # tip.regex <- "^(.*)-[0-9].*_read_([0-9]+)_count_([0-9]+)$"
  # mode <- "s"
  # zero.length.tips.count <- F
  # sankhoff.k <- 10
  # break.ties.unsampled <- F
  #  
  # MRSA example
  
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/Thai MRSA 6/Matthew/")
  output.dir <- "/Users/twoseventwo/Dropbox (Infectious Disease)/Thai MRSA 6/Matthew/refinement"
  tree.file.names <- "RAxML_bipartitions.ST239_no_bootstraps_T056corr.tree"
  blacklist.files <- NULL
  out.identifier <- "mrsa_k10_bp"
  root.name <- "TW20"
  tip.regex <- "^([ST][0-9][0-9][0-9]_[A-Z])[A-Z0-9]*_[A-Z][0-9][0-9]$"
  mode <- "s"
  zero.length.tips.count <- F
  sankhoff.k <- 10
  break.ties.unsampled <- F
  
  if(0)
  {
	  script.dir		<- '/Users/Oliver/git/phylotypes/tools'
	  tree.file.names 	<- '/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/pty_16-09-15-09-43-21/ptyr10_'
	  blacklist.files	<- '/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/pty_16-09-15-09-43-21/ptyr10_blacklist_'	  
	  out.identifier 	<- "ptyr10_"
	  output.dir		<- "/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/pty_16-09-15-09-43-21"
	  root.name 		<- "REF_CPX_AF460972"
	  tip.regex 		<- "^(.*)_read_([0-9]+)_count_([0-9]+)$"	  
	  mode 				<- "r"
	  pdf.hm 			<- 0.15
	  pdf.w 			<- 30	  
	  zero.length.tips.count	<- FALSE
  }
}

suppressMessages(require(phangorn, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phytools, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ggplot2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(ggtree, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(data.table, quietly=TRUE, warn.conflicts=FALSE))

source(file.path(script.dir, "TransmissionUtilityFunctions.R"))
source(file.path(script.dir, "SubtreeMethods.R"))
source(file.path(script.dir, "WriteAnnotatedTrees.R"))

#
#	define internal functions
#
split.patients.to.subtrees<- function(tree.file.name, mode, blacklist.file, root.name, tip.regex, zero.length.tips.count, sankhoff.k, break.ties.unsampled)
{
	cat("SplitPatientsToSubtrees.R run on: ", tree.file.name,", rules = ",mode,"\n", sep="")
	
	tree <- read.tree(tree.file.name)
	tree <- unroot(tree)
	tree <- di2multi(tree, tol = 1E-5)
	if(!is.null(root.name)){
		tree <- root(tree, outgroup = which(tree$tip.label==root.name), resolve.root = T)
	}
	
	tip.labels <- tree$tip.label
	
	blacklist <- vector()
	
	if(!is.null(blacklist.file)){
		if(file.exists(blacklist.file)){
			cat("Reading blacklist file",blacklist.file,'\n')
			blacklisted.tips <- read.table(blacklist.file, sep=",", header=F, stringsAsFactors = F, col.names="read")
			if(nrow(blacklisted.tips)>0){
				blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
			}
		} else {
			warning(paste("File ",blacklist.file," does not exist; skipping.",paste=""))
		}
	} 
	
	
	patient.ids <- sapply(tip.labels, function(x) patient.from.label(x, tip.regex))
	
	if(length(blacklist)>0){
		patients <- unique(na.omit(patient.ids[-blacklist]))
	} else {
		patients <- unique(na.omit(patient.ids))
	}
	
	if(length(patients)==0){
	  stop("No patient IDs detected")
	} 

	patient.tips <-
		 lapply(patients, function(x)  setdiff(which(patient.ids==x), blacklist))
	names(patient.tips) <- patients
	
	patient.mrcas <- lapply(patient.tips, function(node) mrca.phylo.or.unique.tip(tree, node, zero.length.tips.count))
	
	results <- split.and.annotate(tree, patients, patient.tips, patient.mrcas, blacklist, tip.regex, mode, sankhoff.k, break.ties.unsampled)
	
	rm(patient.tips)
	rm(patient.mrcas)
	
	node.shapes <- rep(FALSE, length(tree$tip.label) + tree$Nnode)
	for(mrca in results$first.nodes){
		node.shapes[mrca] <- TRUE
	}
			
	temp.ca <- rep(NA, length(tree$tip.label) + tree$Nnode)
	
	for(item in seq(1, length(results$assocs))){
		if(!is.null(results$assocs[[item]])){
			if(!(results$assocs[[item]] %in% c("*", "unsampled"))) {
				temp.ca[item] <- results$assocs[[item]]
			}
		}
	}
	
	temp.ca.pat <- sapply(temp.ca, function(x) unlist(strsplit(x, "-S"))[1] )
	names(temp.ca.pat) <- NULL
	temp.ca.pat <- factor(temp.ca.pat, levels = sample(levels(as.factor(temp.ca.pat))))
	attr(tree, 'SPLIT') <- factor(temp.ca)
	attr(tree, 'INDIVIDUAL') <- temp.ca.pat
	attr(tree, 'NODE_SHAPES') <- node.shapes
			
	orig.patients <- vector()
	patient.splits <- vector()
	tip.names <- vector()	
	for(patient.plus in results$split.patients){
		tips <- results$split.tips[[patient.plus]]
		orig.patients <- c(orig.patients, rep(unlist(strsplit(patient.plus, "-S"))[1], length(tips)))
		patient.splits <- c(patient.splits, rep(patient.plus, length(tips)))
		tip.names <- c(tip.names, tree$tip.label[tips])
	}	
	
	rm(results)
	
	rs.subtrees <- data.frame(orig.patients, patient.splits, tip.names)
			
	list(tree=tree, rs.subtrees=rs.subtrees)
}
#
#	start script
#

#	check if 'tree.file.names' is tree
single.file	<- file.exists(tree.file.names)
#
#	if is tree, single file mode:
#
if(single.file) {
	#	if 'tree.file.names' is tree, process just one tree
	tree.file.name		<- tree.file.names[1]	
	blacklist.file		<- blacklist.files[1]
	tmp					<- split.patients.to.subtrees(tree.file.name, mode, blacklist.file, root.name, tip.regex, zero.length.tips.count, sankhoff.k, break.ties.unsampled)
	tree				<- tmp[['tree']]	
	rs.subtrees			<- tmp[['rs.subtrees']]
	#
	#	write rda file (potentially before plotting fails so we can recover)
	#
	tmp 				<- file.path(output.dir,paste('Subtrees_',mode,'_',out.identifier,'.rda',sep=''))
	cat("Writing output to file",tmp,"...\n")
	save(rs.subtrees, tree, file=tmp)		
	#
	#	plot tree
	#
	cat("Drawing tree...\n")
	tree.display 		<- ggtree(tree, aes(color=INDIVIDUAL)) +
									geom_point2(shape = 16, size=3, aes(subset=NODE_SHAPES)) +
									scale_fill_hue(na.value = "black") +
									scale_color_hue(na.value = "black") +
									theme(legend.position="none") +
									geom_tiplab(aes(col=INDIVIDUAL)) + 
	                				geom_treescale(width=0.01, y=-5, offset=1.5)
	
	x.max <- ggplot_build(tree.display)$layout$panel_ranges[[1]]$x.range[2]
	
	tree.display <- tree.display + ggplot2::xlim(0, 1.1*x.max)
	tree.display	
	tmp	<- file.path(output.dir,paste('Tree_',mode,'_',out.identifier,'.pdf',sep=''))
	cat("Plot to file",tmp,"...\n")
	ggsave(tmp, device="pdf", height = pdf.hm*length(tree$tip.label), width = pdf.w, limitsize = F)

	tmp 				<- file.path(output.dir, paste('ProcessedTree_',mode,'_',out.identifier,'.tree',sep=''))
	cat("Writing rerooted, multifurcating, annotated tree to file",tmp,"...\n")	
	write.ann.nexus(tree, file=tmp, annotations = c("INDIVIDUAL", "SPLIT"))

	#
	#	write csv file
	#
	
	tmp 				<- file.path(output.dir, paste('Subtrees_',mode,'_',out.identifier,'.csv',sep=''))
	cat("Writing output to file",tmp,"...\n")	
	write.csv(rs.subtrees, file=tmp, row.names = F, quote=F)	
} else {	
	prefix.wfrom 		<- 'Window_'
	prefix.wto 			<- 'Window_[0-9]+_to_'
	tree.file.names		<- sort(list.files(dirname(tree.file.names), pattern=paste(basename(tree.file.names),'.*\\.tree$',sep=''), full.names=TRUE))
	#	code below makes the script more robust to previous errors in gen blacklist files / no trees
	#	read tree file names and define windows
	df	<- data.table(	TF=tree.file.names, 
						W_FROM= as.integer(gsub(prefix.wfrom,'',regmatches(tree.file.names, regexpr(paste(prefix.wfrom,'[0-9]+',sep=''),tree.file.names)))),
						W_TO= as.integer(gsub(prefix.wto,'',regmatches(tree.file.names, regexpr(paste(prefix.wto,'[0-9]+',sep=''),tree.file.names))))
						)					
	if(is.null(blacklist.files))
		df[, BF:=NA_character_]
	if(!is.null(blacklist.files))
	{
		#	read blacklist file names and define windows
		tmp <- sort(list.files(dirname(blacklist.files), pattern=paste(basename(blacklist.files),'.*\\.csv$',sep=''), full.names=TRUE))	
		tmp	<- data.table(	BF=tmp, 
				W_FROM= as.integer(gsub(prefix.wfrom,'',regmatches(tmp, regexpr(paste(prefix.wfrom,'[0-9]+',sep=''),tmp)))),
				W_TO= as.integer(gsub(prefix.wto,'',regmatches(tmp, regexpr(paste(prefix.wto,'[0-9]+',sep=''),tmp))))
				)
		#	merge so that each black list file corresponds to a tree file (by window coordinates)
		df	<- merge(df,tmp,by=c('W_FROM','W_TO'),all=1)	
		#	handle missing tree and blacklist files
		tmp	<- subset(df, is.na(TF))[, paste(W_FROM, collapse=',')]
		if(nchar(tmp))	cat('\nNo tree files for windows',tmp,'Ignoring blacklist files.')
		df	<- subset(df, !is.na(TF))
		tmp	<- subset(df, is.na(BF))[, paste(W_FROM, collapse=',')]
		if(nchar(tmp))	cat('\nNo blacklist files for windows',tmp)			
	}
	for(tree.i in seq_len(nrow(df)))
	{
		tree.file.name		<- df[tree.i, TF]
		blacklist.file		<- df[tree.i, BF]
		if(is.na(blacklist.file))	
			blacklist.file	<- NULL
		out.identifier		<- gsub('\\.tree$','',basename(tree.file.name))
		tmp					<- split.patients.to.subtrees(tree.file.name, mode, blacklist.file, root.name, tip.regex, zero.length.tips.count, sankhoff.k, break.ties.unsampled)
		tree				<- tmp[['tree']]	
		rs.subtrees			<- tmp[['rs.subtrees']]
		#
		#	write rda file (potentially before plotting fails so we can recover)
		#
		tmp 				<- file.path(output.dir,paste('Subtrees_',mode,'_',out.identifier,'.rda',sep=''))
		cat("Writing output to file",tmp,"...\n")
		save(rs.subtrees, tree, file=tmp)	
	
		#
		#	plot tree
		#
		cat("Drawing tree...\n")
		tree.display 		<- ggtree(tree, aes(color=INDIVIDUAL)) +
				geom_point2(shape = 16, size=3, aes(subset=NODE_SHAPES)) +
				scale_fill_hue(na.value = "black") +
				scale_color_hue(na.value = "black") +
				theme(legend.position="none") +
				geom_tiplab(aes(col=INDIVIDUAL)) + 
		    geom_treescale(width = 0.01, y=-5, offset=1.5)
		
		x.max <- ggplot_build(tree.display)$layout$panel_ranges[[1]]$x.range[2]
		
		tree.display <- tree.display + ggplot2::xlim(0, 1.1*x.max)
		tree.display	
		tmp					<- file.path(output.dir,paste('Tree_',mode,'_',out.identifier,'.pdf',sep=''))
		cat("Plot to file",tmp,"...\n")
		ggsave(tmp, device="pdf", height = pdf.hm*length(tree$tip.label), width = pdf.w, limitsize = F)
		#
		#	write csv file
		#
		tmp 				<- file.path(output.dir, paste('Subtrees_',mode,'_',out.identifier,'.csv',sep=''))
		cat("Writing output to file",tmp,"...\n")	
		write.csv(rs.subtrees, file=tmp, row.names = F, quote=F)
		
		tmp 				<- file.path(output.dir, paste('ProcessedTree_',mode,'_',out.identifier,'.tree',sep=''))
		cat("Writing rerooted, multifurcating, annotated tree to file",tmp,"...\n")	
		write.ann.nexus(tree, file=tmp, annotations = c("INDIVIDUAL", "SPLIT"))		
	}
}
  
