#!/usr/bin/env Rscript

list.of.packages <- c("prodlim","reshape2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos="http://cran.ma.imperial.ac.uk/")
#
#	constants
#
prefix.wfrom 		<- 'Window_'
prefix.wto 			<- 'Window_[0-9]+_to_'
#
#
#
command.line <- TRUE
if(command.line){
  suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
  #	OR line breaks result in error "rjson::fromJSON(output) : unexpected character 'F'"
  tmp	<- "Summarise topological relationships suggesting direction of transmission across windows. Outputs a .csv file of relationships between patient IDs. Relationship types: 'anc' indicates a topology that suggests the patient in column 1 infected the patient in column 2, 'cher' where there is only one subtree per patient and the MRCAs of both have the same parent in the phylogeny, 'unint' where cher is not present reads from the patients form sibling clades from an unsampled common ancestor patient (and that neither has any sampled ancestors occuring later in the transmission chain than that common ancestor), 'int' one where reads from both patients are intermingled and the direction of transmission cannot be established. (More documentation to follow.)"
  arg_parser = ArgumentParser(description=tmp)
  arg_parser$add_argument("-l", "--filesAreLists", action="store_true", default=FALSE, help="If present, arguments specifying input files will be parsed as lists of files separated by colons. If absent, they will be parsed as a string which each file of that type begins with.")
  arg_parser$add_argument("-s", "--summaryFile", action="store", help="The full output file from SummaryStatistics.R; necessary only to identify windows in which no reads are present from each patient. If absent, window counts will be given without denominators.")
  arg_parser$add_argument("-m", "--minThreshold", action="store", default=1, type="integer", help="Relationships between two patients will only appear in output if a transmission chain between them appears in at least these many windows (default 1). High numbers are useful for drawing figures in e.g. Cytoscape with few enough arrows to be comprehensible. The script is extremely slow if this is set to 0.")
  arg_parser$add_argument("-p", "--allowSplits", action="store_true", default=FALSE, help="If absent, directionality is only inferred between pairs of patients whose reads are not split; this is more conservative.")
  arg_parser$add_argument("-d", "--detailedOutput", action="store", help="If present, a file describing the relationships between each pair of patients on each window will be written to the specified path in .rda format")
  arg_parser$add_argument("idFile", action="store", help="A file containing a list of the IDs of all the patients to calculate and display statistics for.")
  arg_parser$add_argument("inputFiles", action="store", help="Either (if -l is present) a list of all input files (output from LikelyTransmissions.R), separated by colons, or (if not) a single string that begins every input file name.")
  arg_parser$add_argument("outputFile", action="store", help="A .csv file to write the output to.")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
  args <- arg_parser$parse_args()
  files.are.lists <- args$filesAreLists
  summary.file <- args$summaryFile
  id.file <- args$idFile
  output.file <- args$outputFile
  script.dir <- args$scriptdir  
  min.threshold <- args$minThreshold
  detailed.output <- args$detailed
  if(is.null(min.threshold)){
    split.threshold <- 1L
  }
  allow.splits <- args$allowSplits
  input.files.name <- args$inputFiles
  
  	
  if(!files.are.lists){    
    input.files <- sort(list.files(dirname(input.files.name), pattern=paste(basename(input.files.name),".*LikelyTransmissions.*csv$",sep=""), full.names=TRUE))
  } else {
    input.files <- unlist(strsplit(input.files.name, ":"))
  }
  
} else {
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161007_couples_w270_rerun/")
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools"
  summary.file <- "ptyr4_patStatsFull.csv"
  id.file <- "ptyr4_patients.txt" 
  detailed.output <- "quickout.csv"
  min.threshold <- 0
  allow.splits <- TRUE
  output.file <- "test.csv"  
  input.files <- sort(list.files("ptyr4_trmStats_LikelyTransmissions", pattern="LikelyTransmissions.csv"))
  if(0)
  {
    script.dir					<- "/Users/Oliver/git/phylotypes/tools"
    summary.file				<- "~/duke/tmp/pty_16-12-16-20-41-31/ptyr22_patStatsFull.csv"
    id.file 					<- "~/duke/tmp/pty_16-12-16-20-41-31/ptyr22_patients.txt"
	input.files.name			<- "~/duke/tmp/pty_16-12-16-20-41-31/ProcessedTree_r_ptyr22_"
    min.threshold				<- 1
    allow.splits 				<- TRUE
    output.file 				<- "~/duke/tmp/pty_16-12-16-20-41-31/ptyr22_trmStats.csv"
	detailed.output				<- "~/duke/tmp/pty_16-12-16-20-41-31/ptyr22_patStatsPerWindow.csv"
    input.files 				<- sort(list.files(dirname(input.files.name), pattern=paste(basename(input.files.name),".*LikelyTransmissions.csv$",sep=''), full.names=TRUE))
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
if(!is.null(summary.file))
{
 	cat("Getting window counts per patient...\n")
	reads.table 	<- as.data.table(read.csv(summary.file, stringsAsFactors = F)[,c("window.start", "window.end", "patient", "reads", "leaves")])
	setnames(reads.table, c("window.start","window.end"), c("W_FROM","W_TO"))
	reads.table[, present:= reads>0]  
} 
#
# Read in the IDs. Remove duplicates. Shuffle their order if desired.
#
patient.ids	<- unique(scan(id.file, what="", sep="\n", quiet=TRUE))
if(!length(patient.ids))
  stop(paste("No IDs found in ", id.file, ". Quitting.\n", sep=""))  
#
# Read likely transmissions files
#
tt	<-	lapply(input.files, function(x){
			tt <- as.data.table(read.table(x, sep=",", header=TRUE, stringsAsFactors=FALSE))
			tt[, FILE:=basename(x)]
			setnames(tt, c('paths21'),c('paths.21'))
			tt
		})	
#
# need to worry about transmissions listed in the wrong direction in the input file
# to avoid very large memory, consolidate by FILE
#
tt	<- lapply(tt, function(x){
			#x	<- tt[[1]]
			tmp	<- copy(x)
			setnames(tmp, c('Patient_1','Patient_2','paths.12','paths.21'), c('Patient_2','Patient_1','paths.21','paths.12'))
			set(tmp, tmp[, which(path.classification=="anc")], 'path.classification', 'desc')
			set(tmp, tmp[, which(path.classification=="desc")], 'path.classification', 'anc')
			set(tmp, tmp[, which(path.classification=="multiAnc")], 'path.classification', 'multiDesc')
			set(tmp, tmp[, which(path.classification=="multiDesc")], 'path.classification', 'multiAnc')
			x	<- rbind(x, tmp)
			setkey(x, Patient_1, Patient_2)
			subset(x, Patient_1<Patient_2)
		})
#
# rbind consolidated files
#
tt	<- do.call('rbind',tt)
setnames(tt, c('Patient_1','Patient_2','path.classification','contiguous','paths.21','paths.12','mean.distance.between.subtrees'), c('pat.1','pat.2','TYPE','CONTIGUOUS','PATHS.21','PATHS.12','PATRISTIC_DISTANCE'))
# change type name depending on allow.splits
if(!allow.splits)
{
	set(tt, tt[, which(TYPE%in%c("multiAnc", "multiDesc") & CONTIGUOUS)], 'TYPE', 'conflict')
	set(tt, tt[, which(TYPE%in%c("multiAnc", "multiDesc") & !CONTIGUOUS)], 'TYPE', 'none')
}
#	check we have patristic distances, paths
stopifnot( !nrow(subset(tt, is.na(PATRISTIC_DISTANCE))) )
stopifnot( !nrow(subset(tt, is.na(PATHS.12))) )
stopifnot( !nrow(subset(tt, is.na(PATHS.21))) )
#	set to numeric
set(tt, NULL, 'PATRISTIC_DISTANCE', tt[, as.numeric(PATRISTIC_DISTANCE)])
# 	add window coordinates
tt[, W_FROM:= tt[,as.integer(gsub(prefix.wfrom,'',regmatches(FILE, regexpr(paste(prefix.wfrom,'[0-9]+',sep=''),FILE))))]]
tt[, W_TO:= tt[, as.integer(gsub(prefix.wto,'',regmatches(FILE, regexpr(paste(prefix.wto,'[0-9]+',sep=''),FILE))))]]
#
#	calculate #unique reads and #reads for each window
#
dp			<-  reads.table[,{
			#	combine by window only for present patients to avoid large mem
			dp			<- data.table(pat.1=NA_character_, pat.2=NA_character_)
			if(length(which(present))>1)
			{
				dp			<- t( combn( patient[present], 2 ) )
				colnames(dp)<- c('pat.1','pat.2')
				dp			<- as.data.table(dp)	
				# make sure combinations are in the direction as in 'tt'
				tmp			<- copy(dp)
				setnames(tmp, c('pat.1','pat.2'), c('pat.2','pat.1'))
				dp			<- rbind(dp, tmp)				
				dp			<- subset(dp, pat.1<pat.2)
			}
			dp						
		}, by=c('W_FROM','W_TO')]
dp			<- subset(dp, !is.na(pat.1) & !is.na(pat.2))
tmp			<- subset(reads.table, present, c(W_FROM, patient, reads, leaves))
setnames(tmp, c("patient","reads","leaves"), c("pat.1","pat.1_reads","pat.1_leaves"))
dp			<- merge(dp, tmp, by=c('W_FROM','pat.1'))
setnames(tmp, c("pat.1","pat.1_leaves","pat.1_reads"), c("pat.2","pat.2_leaves","pat.2_reads"))
dp			<- merge(dp, tmp, by=c('W_FROM','pat.2'))
#
#	merge reads/leaves with tt
#
stopifnot(nrow(dp)==nrow(tt))
tt			<- merge(tt, dp, by=c('pat.1','pat.2','W_FROM','W_TO'))
stopifnot(nrow(dp)==nrow(tt))
#	rename TYPE
set(tt, tt[,which(TYPE=='anc')], 'TYPE','anc_12')
set(tt, tt[,which(TYPE=='desc')], 'TYPE','anc_21')
set(tt, tt[,which(TYPE=='multiAnc')], 'TYPE','multi_anc_12')
set(tt, tt[,which(TYPE=='multiDesc')], 'TYPE','multi_anc_21')

#	reorder
setkey(tt, W_FROM, W_TO, pat.1, pat.2)
tt			<- subset(tt, select=c('W_FROM','W_TO','pat.1','pat.2','TYPE','PATRISTIC_DISTANCE','CONTIGUOUS','PATHS.12','PATHS.21','pat.1_leaves','pat.1_reads','pat.2_leaves','pat.2_reads'))
setnames(tt, colnames(tt),toupper(colnames(tt)))
#	write to file
if(!is.null(detailed.output))
{	
	cat("Save detailed per window table:",detailed.output)
	save(tt, file=gsub('\\.csv','\\.rda',detailed.output))
}
# 	Not currently in use
#	make summary of transmission assignments across all columns
#
#cat("Making summary output table...\n")
#set(tt, NULL, c('pat.1_leaves','pat.1_reads','pat.2_leaves','pat.2_reads'),NULL)
#tmp		<- tt[, list(windows=length(W_FROM)), by=c('pat.1','pat.2','TYPE')]
#tt		<- merge(tt, tmp, by=c('pat.1','pat.2','TYPE'))
#	evaluate total.trans, denominator, and fraction
#tmp		<- tt[, list(denominator=length(W_FROM), total.notdisc=length(which(TYPE!='disconnected')), total.trans=length(which(grepl('anc',TYPE)))), by=c('pat.1','pat.2')]
#tmp		<- subset(tmp, total.notdisc>0)
#tt		<- merge(tt, tmp, by=c('pat.1','pat.2'))
#tt[, fraction:=paste(windows,'/',denominator,sep='')]
#	convert "anc_12" and "ans_21" to "anc" depending on direction
#tt[, DUMMY:=NA_character_]
#tmp			<- tt[, which(TYPE=="anc_12")]
#set(tt, tmp, 'TYPE', "anc")
#tmp			<- tt[, which(TYPE=="anc_21")]
#set(tt, tmp, 'DUMMY', tt[tmp, pat.1])
#set(tt, tmp, 'pat.1', tt[tmp, pat.2])
#set(tt, tmp, 'pat.2', tt[tmp, DUMMY])
#set(tt, tmp, 'TYPE', "anc")
#tt[, DUMMY:=NULL]
#	write to file
#setkey(tt, pat.1, pat.2, TYPE)
#tt			<- subset(tt, TYPE!='disconnected')
#write.csv(subset(tt, total.trans>=min.threshold), file=output.file, row.names=FALSE, quote=FALSE)

# Not currently in use - for dating

# patient.data <- read.table("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/Collective_notebook/BEEHIVE_summary29.02.2016.csv", 
#                            sep=",", stringsAsFactors = F, header=T)
# 
# date.info <- patient.data[,c("PATIENT", "Date.first.pos", "Date.sampled")]
# date.info$Date.first.pos <- parse_date_time(date.info$Date.first.pos, orders="ymd")
# date.info$Date.first.pos <- decimal_date(date.info$Date.first.pos)
# date.info$Date.sampled <- parse_date_time(date.info$Date.sampled, orders="ymd")
# date.info$Date.sampled <- decimal_date(date.info$Date.sampled)
# date.info$PATIENT <- paste(date.info$PATIENT,"-1",sep="")
# date.info$PATIENT[7] <- "BEE0006-2"
# date.info <- date.info[which(date.info$PATIENT %in% patient.ids),]
# 
# 
# first <- T
# 
# earlier.patient <- vector()
# later.patient <- vector()
# time.difference <- vector()
# net.transmission.windows <- vector()
# transmissions.inferred <- vector()
# 
# for(early in patient.ids){
#   for(late in patient.ids){
#     if(early!=late){
#       early.date <- date.info$Date.sampled[which(date.info$PATIENT==early)]
#       late.date <- date.info$Date.sampled[which(date.info$PATIENT==late)]
#       if(length(early.date)>0 & length(late.date)>0){
#         if(!is.na(early.date) & !is.na(late.date)){
#           if(early.date<=late.date){
#             temporal <- new.out.2[which(new.out.2$pat.1==early & new.out.2$pat.2==late & new.out.2$TYPE=="trans"),]
#             contratemporal <- new.out.2[which(new.out.2$pat.1==late & new.out.2$pat.2==early & new.out.2$TYPE=="trans"),]
#             if(nrow(temporal) > 0){
#               forwards.windows <- temporal$windows
#               total.trans <- temporal$total.trans
#             } else {
#               forwards.windows <- 0
#             }
#             if(nrow(contratemporal) > 0){
#               backwards.windows <- contratemporal$windows
#               total.trans <- contratemporal$total.trans
#             } else {
#               backwards.windows <- 0
#             }
#             
#             if(forwards.windows > 0 | backwards.windows > 0){
#               cat("Calculating for ",early," and ",late,"\n",sep="")
#               earlier.patient <- c(earlier.patient, early)
#               later.patient <- c(later.patient, late)
#               time.difference <- c(time.difference, late.date - early.date)
#               net.transmission.windows <- c(net.transmission.windows, forwards.windows - backwards.windows)
#               transmissions.inferred <- c(transmissions.inferred, total.trans)
#             }
#           }
#         }
#       }
#     }
#   }
# }
# 
# scatter.data <- data.frame(earlier.patient, later.patient, time.difference, net.transmission.windows, transmissions.inferred)
# 
# 
# graph.diffs <- ggplot(scatter.data, aes(time.difference, net.transmission.windows))
# 
# graph.diffs + geom_point(aes(size=transmissions.inferred), alpha=0.33) +
#   theme_bw() +
#   labs(x="Years between sampling dates", y="Net windows suggesting transmission from early\n sample to late sample (R-S classification)") +
#   scale_size("Windows with\n inferred\n transmissions", range = c(0, 5), breaks=c(1,5,10,15,20,25,30),)
# 
# ggsave("DiffSamp_RS.pdf", width=10, height=6)
# 
# 
# 
# first <- T
# 
# earlier.patient <- vector()
# later.patient <- vector()
# time.difference <- vector()
# net.transmission.windows <- vector()
# transmissions.inferred <- vector()
# 
# for(early in patient.ids){
#   for(late in patient.ids){
#     if(early!=late){
#       early.date <- date.info$Date.sampled[which(date.info$PATIENT==early)]
#       late.date <- date.info$Date.sampled[which(date.info$PATIENT==late)]
#       if(length(early.date)>0 & length(late.date)>0){
#         if(!is.na(early.date) & !is.na(late.date)){
#           if(early.date<=late.date){
#             temporal <- new.out[which(new.out$pat.1==early & new.out$pat.2==late & new.out$TYPE=="MP"),]
#             contratemporal <- new.out[which(new.out$pat.1==late & new.out$pat.2==early & new.out$TYPE=="MP"),]
#             if(nrow(temporal) > 0){
#               forwards.windows <- temporal$windows
#               total.trans <- temporal$total.trans
#             } else {
#               forwards.windows <- 0
#             }
#             if(nrow(contratemporal) > 0){
#               backwards.windows <- contratemporal$windows
#               total.trans <- contratemporal$total.trans
#             } else {
#               backwards.windows <- 0
#             }
#             
#             if(forwards.windows > 0 | backwards.windows > 0){
#               cat("Calculating for ",early," and ",late,"\n",sep="")
#               earlier.patient <- c(earlier.patient, early)
#               later.patient <- c(later.patient, late)
#               time.difference <- c(time.difference, late.date - early.date)
#               net.transmission.windows <- c(net.transmission.windows, forwards.windows - backwards.windows)
#               transmissions.inferred <- c(transmissions.inferred, total.trans)
#             }
#           }
#         }
#       }
#     }
#   }
# }
# 
# scatter.data <- data.frame(earlier.patient, later.patient, time.difference, net.transmission.windows, transmissions.inferred)
# 
# graph.diffs <- ggplot(scatter.data, aes(time.difference, net.transmission.windows))
# 
# graph.diffs + geom_point(aes(size=transmissions.inferred), alpha=0.33) +
#   theme_bw() +
#   labs(x="Years between sampling dates", y="Net windows suggesting transmission from early\n sample to late sample (old classification)") +
#   scale_size("Windows with\n inferred\n transmissions", range = c(0, 5), breaks=c(0,5,10,15,20,25,30),)
# 
# ggsave("DiffSamp_old.pdf", width=10, height=6)
# 
# 
# first <- T
# 
# earlier.patient <- vector()
# later.patient <- vector()
# time.difference <- vector()
# net.transmission.windows <- vector()
# transmissions.inferred <- vector()
# 
# for(early in patient.ids){
#   for(late in patient.ids){
#     if(early!=late){
#       early.date <- date.info$Date.first.pos[which(date.info$PATIENT==early)]
#       late.date <- date.info$Date.first.pos[which(date.info$PATIENT==late)]
#       if(length(early.date)>0 & length(late.date)>0){
#         if(!is.na(early.date) & !is.na(late.date)){
#           if(early.date<=late.date){
#             temporal <- new.out.2[which(new.out.2$pat.1==early & new.out.2$pat.2==late & new.out.2$TYPE=="trans"),]
#             contratemporal <- new.out.2[which(new.out.2$pat.1==late & new.out.2$pat.2==early & new.out.2$TYPE=="trans"),]
#             if(nrow(temporal) > 0){
#               forwards.windows <- temporal$windows
#               total.trans <- temporal$total.trans
#             } else {
#               forwards.windows <- 0
#             }
#             if(nrow(contratemporal) > 0){
#               backwards.windows <- contratemporal$windows
#               total.trans <- contratemporal$total.trans
#             } else {
#               backwards.windows <- 0
#             }
#             
#             if(forwards.windows > 0 | backwards.windows > 0){
#               cat("Calculating for ",early," and ",late,"\n",sep="")
#               earlier.patient <- c(earlier.patient, early)
#               later.patient <- c(later.patient, late)
#               time.difference <- c(time.difference, late.date - early.date)
#               net.transmission.windows <- c(net.transmission.windows, forwards.windows - backwards.windows)
#               transmissions.inferred <- c(transmissions.inferred, total.trans)
#             }
#           }
#         }
#       }
#     }
#   }
# }
# 
# scatter.data <- data.frame(earlier.patient, later.patient, time.difference, net.transmission.windows, transmissions.inferred)
# 
# 
# graph.diffs <- ggplot(scatter.data, aes(time.difference, net.transmission.windows))
# 
# graph.diffs + geom_point(aes(size=transmissions.inferred), alpha=0.33) +
#   theme_bw() +
#   labs(x="Years between dates of first positive test", y="Net windows suggesting transmission from early\n sample to late sample (R-S classification)") +
#   scale_size("Windows with\n inferred\n transmissions", range = c(0, 5), breaks=c(1,5,10,15,20,25,30),)
# 
# ggsave("DiffFP_RS.pdf", width=10, height=6)
# 
# first <- T
# 
# earlier.patient <- vector()
# later.patient <- vector()
# time.difference <- vector()
# net.transmission.windows <- vector()
# transmissions.inferred <- vector()
# 
# for(early in patient.ids){
#   for(late in patient.ids){
#     if(early!=late){
#       early.date <- date.info$Date.first.pos[which(date.info$PATIENT==early)]
#       late.date <- date.info$Date.first.pos[which(date.info$PATIENT==late)]
#       if(length(early.date)>0 & length(late.date)>0){
#         if(!is.na(early.date) & !is.na(late.date)){
#           if(early.date<=late.date){
#             temporal <- new.out[which(new.out$pat.1==early & new.out$pat.2==late & new.out$TYPE=="MP"),]
#             contratemporal <- new.out[which(new.out$pat.1==late & new.out$pat.2==early & new.out$type=="MP"),]
#             if(nrow(temporal) > 0){
#               forwards.windows <- temporal$windows
#               total.trans <- temporal$total.trans
#             } else {
#               forwards.windows <- 0
#             }
#             if(nrow(contratemporal) > 0){
#               backwards.windows <- contratemporal$windows
#               total.trans <- contratemporal$total.trans
#             } else {
#               backwards.windows <- 0
#             }
#             
#             if(forwards.windows > 0 | backwards.windows > 0){
#               cat("Calculating for ",early," and ",late,"\n",sep="")
#               earlier.patient <- c(earlier.patient, early)
#               later.patient <- c(later.patient, late)
#               time.difference <- c(time.difference, late.date - early.date)
#               net.transmission.windows <- c(net.transmission.windows, forwards.windows - backwards.windows)
#               transmissions.inferred <- c(transmissions.inferred, total.trans)
#             }
#           }
#         }
#       }
#     }
#   }
# }
# 
# scatter.data <- data.frame(earlier.patient, later.patient, time.difference, net.transmission.windows, transmissions.inferred)
# 
# graph.diffs <- ggplot(scatter.data, aes(time.difference, net.transmission.windows))
# 
# graph.diffs + geom_point(aes(size=transmissions.inferred), alpha=0.33) +
#   theme_bw() +
#   labs(x="Years between sampling dates", y="Net windows suggesting transmission from early\n sample to late sample (old classification)") +
#   scale_size("Windows with\n inferred\n transmissions", range = c(0, 5), breaks=c(0,5,10,15,20,25,30),)
# 
# ggsave("DiffFP_old.pdf", width=10, height=6)
