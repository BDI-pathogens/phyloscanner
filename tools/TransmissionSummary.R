#!/usr/bin/env Rscript

list.of.packages <- c("prodlim","reshape2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos="http://cran.ma.imperial.ac.uk/")

tree.fe <- ".tree"
csv.fe <- ".csv"

command.line <- T
if(command.line){
  suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
  #	OR line breaks result in error "rjson::fromJSON(output) : unexpected character 'F'"
  tmp	<- "Summarise topological relationships suggesting direction of transmission across windows. Outputs a .csv file of relationships between patient IDs. Relationship types: 'anc' indicates a topology that suggests the patient in column 1 infected the patient in column 2, 'cher' where there is only one subtree per patient and the MRCAs of both have the same parent in the phylogeny, 'unint' where cher is not present reads from the patients form sibling clades from an unsampled common ancestor patient (and that neither has any sampled ancestors occuring later in the transmission chain than that common ancestor), 'int' one where reads from both patients are intermingled and the direction of transmission cannot be established. (More documentation to follow.)"
  arg_parser = ArgumentParser(description=tmp)
  arg_parser$add_argument("-s", "--summaryFile", action="store", help="The full output file from SummaryStatistics.R; necessary only to identify windows in which no reads are present from each patient. If absent, window counts will be given without denominators.")
  arg_parser$add_argument("-m", "--minThreshold", action="store", default=1, type="integer", help="Relationships between two patients will only appear in output if a transmission chain between them appears in at least these many windows (default 1). High numbers are useful for drawing figures in e.g. Cytoscape with few enough arrows to be comprehensible. The script is extremely slow if this is set to 0.")
  arg_parser$add_argument("-c", "--distanceThreshold", action="store", default=-1, help="Maximum distance threshold on a window for a relationship to be reconstructed between two patients on that window.")
  arg_parser$add_argument("-p", "--allowMultiTrans", action="store_true", default=FALSE, help="If absent, directionality is only inferred between pairs of patients where a single clade from one patient is nested in one from the other; this is more conservative")
  arg_parser$add_argument("-d", "--detailedOutput", action="store", help="If present, a file describing the relationships between each pair of patients on each window will be written to the specified path in .rda format")
  arg_parser$add_argument("idFile", action="store", help="A file containing a list of the IDs of all the patients to calculate and display statistics for.")
  arg_parser$add_argument("inputFiles", action="store", help="Either (if -l is present) a list of all input files (output from ClassifyRelationships.R), separated by colons, or (if not) a single string that begins every input file name.")
  arg_parser$add_argument("outputFile", action="store", help="A .csv file to write the output to.")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
  arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")
  args <- arg_parser$parse_args()
  summary.file <- args$summaryFile
  id.file <- args$idFile
  output.file <- args$outputFile
  script.dir <- args$scriptdir  
  verbose <- args$verbose
  min.threshold <- as.numeric(args$minThreshold)
  dist.threshold <- as.numeric(args$distanceThreshold)
  if(dist.threshold==-1){
    dist.threshold <- Inf
  }
  detailed.output <- args$detailed
  if(is.null(min.threshold)){
    split.threshold <- 1L
  }
  allow.mt <- args$allowMultiTrans
  input.file.name <- args$inputFiles
  
  source(file.path(script.dir, "TreeUtilityFunctions.R"))
  source(file.path(script.dir, "GeneralFunctions.R"))
  
  input.files <- list.files.mod(dirname(input.file.name), pattern=paste(basename(input.file.name)), full.names=TRUE)

  # We only want the *classification*.csv files produced by
  # ClassifyRelationships.R, not the collapsedTree files if present.
  # NB this regex needs to match the hard-coded file naming for collapsed trees.
  input.files <- input.files[! grepl('collapsedTree', input.files)]
  
} else {
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20161013/")
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools"
  summary.file	<- NULL
  id.file 					<- "patientIDList.txt"
  input.file.name			<- "/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20161013/Classifications_s/Classification_s_run20161013D_inWindow_"
  min.threshold				<- 16.5
  dist.threshold <- 0.05183
  allow.mt 				<- TRUE
  output.file 				<- "run20161013D.csv"
  detailed.output				<- NULL
  input.files 				<- sort(list.files.mod(dirname(input.file.name), pattern=paste(basename(input.file.name)), full.names=TRUE))
  
  setwd("/Users/twoseventwo/Documents/Croucher alignments")
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools"
  summary.file <- "stats_patStatsFull.csv"
  id.file <- "patients.txt"
  input.file.name <- "pneumo_classification_"
  dist.threshold <- 0.005
  min.threshold <- 50
  allow.mt <- T
  output.file <- "testn_0.005.csv"
  detailed.output <- "test2.csv"
  
  
  setwd("/Users/twoseventwo/Downloads/PossibleDuals/")
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools"
  summary.file <- "summary_patStatsFull.csv"
  id.file <- "PossibleDualsForPaper.txt"
  input.file.name <- "Classification_s_classification_InWindow_"
  dist.threshold <- 1
  min.threshold <- 0
  allow.mt <- T
  output.file <- "test1.csv"
  detailed.output <- "test2.csv"
  
  
  input.files <- list.files.mod(dirname(input.file.name), pattern=paste(basename(input.file.name)), full.names=TRUE)
  
  input.files <- input.files[! grepl('collapsedTree', input.files)]
  input.files <- input.files[! grepl('_classification\\.', input.files)]
  
  if(0)
  {
    script.dir					<- "/Users/Oliver/git/phylotypes/tools"
    summary.file				<- "/Users/Oliver/duke/tmp/pty_17-04-04-16-26-50/ptyr22_patStatsFull.csv"
    id.file 					<- "/Users/Oliver/duke/tmp/pty_17-04-04-16-26-50/ptyr22_patient.txt"
	input.file.name				<- "/Users/Oliver/duke/tmp/pty_17-04-04-16-26-50/ptyr22_classification_InWindow_"
    min.threshold				<- 1
    allow.mt 				<- TRUE
    output.file 				<- "/Users/Oliver/duke/tmp/pty_17-04-04-16-26-50/ptyr22_trmStats.csv"
	detailed.output				<- "/Users/Oliver/duke/tmp/pty_17-04-04-16-26-50/ptyr22_patStatsPerWindow.csv"
    input.files 				<- sort(list.files(dirname(input.file.name), pattern=basename(input.file.name), full.names=TRUE))
  }
}
#
suppressMessages(library(prodlim, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(reshape2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(gdata, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ggplot2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(data.table, quietly=TRUE, warn.conflicts=FALSE))

#
# Get the denominators
#
reads.table	<- NULL		# can replace prev give.denom by is.null(reads.table)
if(!is.null(summary.file)){
 	if (verbose) cat("Getting window counts per patient from ",summary.file,"...\n", sep="")
	reads.table 	<- read.csv(summary.file, stringsAsFactors = F)
  required.columns <- c("file.suffix", "id", "reads", "tips")
  missing.columns <-
  required.columns[! required.columns %in% colnames(reads.table)]
  if (length(missing.columns) > 0) {
    stop(paste0("The following required columns are missing from the summary",
    " csv file: ", paste(missing.columns, collapse=', '), ". Quitting."))
  }
	reads.table 	<- reads.table[,c("file.suffix", "id", "reads", "tips")]
	reads.table 	<- as.data.table(reads.table)
	setnames(reads.table, c("file.suffix"), c("SUFFIX"))
	reads.table$SUFFIX <- as.character(reads.table$SUFFIX)
	reads.table[, present:= reads>0]  
} 
#
# Read in the IDs. Remove duplicates. Shuffle their order if desired.
#
patient.ids	<- unique(scan(id.file, what="", sep="\n", quiet=TRUE))
if(!length(patient.ids))
  stop(paste("No IDs found in ", id.file, ". Quitting.\n", sep=""))  
#
# Read classification files
#
tt	<-	lapply(input.files, function(x){
      if (verbose) cat("Reading window input file ",x,"\n", sep="")
			tt <- as.data.table(read.table(x, sep=",", header=TRUE, stringsAsFactors=FALSE))
			tt[, SUFFIX:=get.suffix(x, input.file.name, csv.fe)]			
			tt
		})	
#
# need to worry about transmissions listed in the wrong direction in the input file
# to avoid very large memory, consolidate by FILE
#
tt	<- lapply(tt, function(x){
			#x	<- tt[[1]]
			tmp	<- copy(x)
			setnames(tmp, c('Patient_1','Patient_2','paths12','paths21'), c('Patient_2','Patient_1','paths21','paths12'))
			set(tmp, tmp[, which(path.classification=="anc")], 'path.classification', 'TMP')
			set(tmp, tmp[, which(path.classification=="desc")], 'path.classification', 'anc')
			set(tmp, tmp[, which(path.classification=="TMP")], 'path.classification', 'desc')
			set(tmp, tmp[, which(path.classification=="multiAnc")], 'path.classification', 'TMP')
			set(tmp, tmp[, which(path.classification=="multiDesc")], 'path.classification', 'multiAnc')
			set(tmp, tmp[, which(path.classification=="TMP")], 'path.classification', 'multiDesc')
			x	<- rbind(x, tmp)
			setkey(x, Patient_1, Patient_2)
			subset(x, Patient_1 < Patient_2)
		})
#
# rbind consolidated files
#
tt	<- do.call('rbind',tt)
# reset names depending on which Classify script was used
if(any('normalised.min.distance.between.subtrees'==colnames(tt))){
  setnames(tt, 'normalised.min.distance.between.subtrees', 'PATRISTIC_DISTANCE')
} else if(any('min.distance.between.subtrees'==colnames(tt)))
	{setnames(tt, 'min.distance.between.subtrees', 'PATRISTIC_DISTANCE')}
setnames(tt, c('Patient_1','Patient_2','path.classification','paths21','paths12','adjacent','contiguous'), c('pat.1','pat.2','TYPE','PATHS.21','PATHS.12','ADJACENT','CONTIGUOUS'))
# change type name depending on allow.mt
if(!allow.mt){
	set(tt, tt[, which(TYPE%in%c("multiAnc", "multiDesc"))], 'TYPE', 'conflict')
}

#	check we have patristic distances, paths
stopifnot( !nrow(subset(tt, is.na(PATRISTIC_DISTANCE))) )
stopifnot( !nrow(subset(tt, is.na(PATHS.12))) )
stopifnot( !nrow(subset(tt, is.na(PATHS.21))) )
#	set to numeric
set(tt, NULL, 'PATRISTIC_DISTANCE', tt[, as.numeric(PATRISTIC_DISTANCE)])
# 	add window coordinates

#	calculate #unique reads and #reads for each window
#
dp	<-  reads.table[,{
			#	combine by window only for present patients to avoid large mem
			dp			<- data.table(pat.1=NA_character_, pat.2=NA_character_)
			if(length(which(present))>1){
				dp			<- t( combn( id[present], 2 ) )
				colnames(dp)<- c('pat.1','pat.2')
				dp			<- as.data.table(dp)	
				# make sure combinations are in the direction as in 'tt'
				tmp			<- copy(dp)
				setnames(tmp, c('pat.1','pat.2'), c('pat.2','pat.1'))
				dp			<- rbind(dp, tmp)				
				dp			<- subset(dp, pat.1<pat.2)
			}
			dp						
		}, by=c('SUFFIX')]
dp <- subset(dp, !is.na(pat.1) & !is.na(pat.2))
tmp	<- subset(reads.table, present, c(SUFFIX, id, reads, tips))
setnames(tmp, c("id","reads","tips"), c("pat.1","pat.1_reads","pat.1_tips"))
dp <- merge(dp, tmp, by=c('SUFFIX','pat.1'))
setnames(tmp, c("pat.1","pat.1_tips","pat.1_reads"), c("pat.2","pat.2_tips","pat.2_reads"))
dp <- merge(dp, tmp, by=c('SUFFIX','pat.2'))
#
#	merge reads/leaves with tt
#
tt			<- merge(tt, dp, by=c('pat.1','pat.2','SUFFIX'))

if(nrow(tt)==0){
  cat("Failed to merge tables; e.g. file suffix",tt$SUFFIX[1]," not found in ",summary.file,"\n",sep="")
  quit(save="\n")
}

#	rename TYPE
set(tt, tt[,which(TYPE=='anc')], 'TYPE','anc_12')
set(tt, tt[,which(TYPE=='desc')], 'TYPE','anc_21')
set(tt, tt[,which(TYPE=='multiAnc')], 'TYPE','multi_anc_12')
set(tt, tt[,which(TYPE=='multiDesc')], 'TYPE','multi_anc_21')

#	reorder
setkey(tt, SUFFIX, pat.1, pat.2)
tt			<- subset(tt, select=c('SUFFIX','pat.1','pat.2','TYPE','PATRISTIC_DISTANCE','ADJACENT','CONTIGUOUS','PATHS.12','PATHS.21','pat.1_tips','pat.1_reads','pat.2_tips','pat.2_reads'))
#tt			<- subset(tt, select=c('SUFFIX','pat.1','pat.2','TYPE','PATRISTIC_DISTANCE','ADJACENT','CONTIGUOUS','PATHS.12','PATHS.21'))
setnames(tt, colnames(tt),toupper(colnames(tt)))
#	write to file
if(!is.null(detailed.output))
{	
 	if (verbose) cat("Save detailed per window table:",detailed.output)
 	save(tt, file=gsub('\\.csv','\\.rda',detailed.output))
}

if (verbose) cat("Making summary output table...\n")
set(tt, NULL, c('PAT.1_TIPS','PAT.1_READS','PAT.2_TIPS','PAT.2_READS','PATHS.12','PATHS.21'),NULL)
#set(tt, NULL, c('PATHS.12','PATHS.21'),NULL)

existence.counts <- tt[, list(both.exist=length(SUFFIX)), by=c('PAT.1','PAT.2')]
tt <- merge(tt, existence.counts, by=c('PAT.1', 'PAT.2'))

tt.close <- tt[which(tt$ADJACENT & tt$PATRISTIC_DISTANCE < dist.threshold ),]
tt.close$NOT.SIBLINGS <- tt.close$ADJACENT & (tt.close$PATRISTIC_DISTANCE < dist.threshold) & tt.close$TYPE!="none"

# How many windows have this relationship, ADJACENT and PATRISTIC_DISTANCE below the threshold?
type.counts	<- tt.close[, list(windows=length(SUFFIX)), by=c('PAT.1','PAT.2','TYPE')]
# How many windows have ADJACENT and PATRISTIC_DISTANCE below the threshold?
any.counts <- tt.close[, list(all.windows=length(SUFFIX)), by=c('PAT.1','PAT.2')]
# How many windows have a relationship other than "none", ADJACENT and PATRISTIC_DISTANCE below the threshold?
ns.counts <- tt.close[, list(ns.windows=length(which(NOT.SIBLINGS))), by=c('PAT.1','PAT.2')]

tt.close		<- merge(tt.close, type.counts, by=c('PAT.1','PAT.2','TYPE'))
tt.close		<- merge(tt.close, any.counts, by=c('PAT.1','PAT.2'))
tt.close		<- merge(tt.close, ns.counts, by=c('PAT.1','PAT.2'))

tt.close[, fraction:=paste(windows,'/',both.exist,sep='')]
#	convert "anc_12" and "ans_21" to "anc" depending on direction
tt.close[, DUMMY:=NA_character_]
tmp			<- tt.close[, which(TYPE=="anc_12")]
set(tt.close, tmp, 'TYPE', "trans")
tmp			<- tt.close[, which(TYPE=="anc_21")]
set(tt.close, tmp, 'DUMMY', tt.close[tmp, PAT.1])
set(tt.close, tmp, 'PAT.1', tt.close[tmp, PAT.2])
set(tt.close, tmp, 'PAT.2', tt.close[tmp, DUMMY])
set(tt.close, tmp, 'TYPE', "trans")

tmp			<- tt.close[, which(TYPE=="multi_anc_12")]
set(tt.close, tmp, 'TYPE', "multi_trans")
tmp			<- tt.close[, which(TYPE=="multi_anc_21")]
set(tt.close, tmp, 'DUMMY', tt.close[tmp, PAT.1])
set(tt.close, tmp, 'PAT.1', tt.close[tmp, PAT.2])
set(tt.close, tmp, 'PAT.2', tt.close[tmp, DUMMY])
set(tt.close, tmp, 'TYPE', "multi_trans")

tt.close[, DUMMY:=NULL]

set(tt.close, NULL, c('SUFFIX', 'ADJACENT','CONTIGUOUS','PATRISTIC_DISTANCE'), NULL)
tt.close <- tt.close[!duplicated(tt.close),]

#	write to file
setkey(tt.close, PAT.1, PAT.2, TYPE)
#
if (verbose) cat('Writing summary to file',output.file,'\n')
write.csv(subset(tt.close, all.windows>=min.threshold), file=output.file, row.names=FALSE, quote=FALSE)

