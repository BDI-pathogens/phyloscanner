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
    summary.file				<- "~/duke/tmp/pty_16-11-17-15-04-28/ptyr22_patStatsFull.csv"
    id.file 					<- "~/duke/tmp/pty_16-11-17-15-04-28/ptyr22_patients.txt"
	input.files.name			<- "~/duke/tmp/pty_16-11-17-15-04-28/ptyr22_"
    min.threshold				<- 1
    allow.splits 				<- TRUE
    output.file 				<- "~/duke/tmp/pty_16-11-17-15-04-28/ptyr22_trmStats.csv"
	detailed.output				<- "~/duke/tmp/pty_16-11-17-15-04-28/ptyr22_patStatsPerWindow.csv"
    input.files 				<- sort(list.files(dirname(input.files.name), pattern=paste(basename(input.files.name),".*LikelyTransmissions.csv$",sep=''), full.names=TRUE))
  }
}

num.windows <- length(input.files)

suppressMessages(library(prodlim, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(reshape2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(gdata, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ggplot2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(data.table, quietly=TRUE, warn.conflicts=FALSE))

# Get the denominators

if(!is.null(summary.file)){
  cat("Getting window counts per patient...\n")
  full.summary <- read.csv(summary.file, stringsAsFactors = F)
  reads.table <- full.summary[,c("window.start", "window.end", "patient", "reads", "leaves")]
  reads.table$present <- reads.table$reads>0
  give.denoms <- T
} else {
  give.denoms <- F
}

# Read in the IDs. Remove duplicates. Shuffle their order if desired.
ids <- scan(id.file, what="", sep="\n", quiet=TRUE)

num.ids <- length(ids)
if (num.ids == 0) {
  cat(paste("No IDs found in ", id.file, ". Quitting.\n", sep=""))
  quit("no", 1)
}
patient.ids <- unique(ids)

window.count <- 0

relationships.ever.suggested <- vector()

cat("Determining pairs for which relationships are ever suggested...\n") 

for(window in seq(1, num.windows)){
  
  window.count <- window.count + 1
  cat("Reading file", input.files[window],"...\n")
  transmissions.table <- read.table(input.files[window], sep=",", header=TRUE, stringsAsFactors=FALSE)
  
  if(min.threshold == 0){
    relationships.ever.suggested <- rbind(relationships.ever.suggested, transmissions.table[,1:2])
  } else {
    relationships.ever.suggested <- rbind(relationships.ever.suggested, transmissions.table[which(transmissions.table$Relationship %in% c("anc", "desc", "intAnc", "intDesc")),1:2])
  }
  
}
relationships.ever.suggested <- as.data.frame(unique(relationships.ever.suggested), stringsAsFactors=FALSE)

# need to worry about transmissions listed in the wrong direction in the input file

reversed.t.e.s <- as.data.frame(cbind(relationships.ever.suggested[,2],relationships.ever.suggested[,1]), stringsAsFactors = FALSE)
matches <- row.match(relationships.ever.suggested, reversed.t.e.s)

to.go <- rep(FALSE, nrow(relationships.ever.suggested))

for(comp in seq(1, length(matches))){
  if(!(is.na(matches[comp])) & (comp < matches[comp])){
    to.go[comp] <- TRUE
  }
}

relationships.ever.suggested <- relationships.ever.suggested[which(!to.go), ]

window.count <- 0

first <- TRUE
#first.sib <- TRUE

cat("Identifying relationships in each window...\n")

for(window in seq(1, num.windows)){
  
  window.vector <- rep(NA, nrow(relationships.ever.suggested))
  dist.vector <- rep(NA, nrow(relationships.ever.suggested))
  p1.vector <- rep(NA, nrow(relationships.ever.suggested))
  p2.vector <- rep(NA, nrow(relationships.ever.suggested))
  #  sib.window.vector <- rep(NA, nrow(transmissions.ever.suggested))
  
  window.count <- window.count + 1
  
  cat("Opening file: ",input.files[window],"\n", sep="")
  
  transmissions.table <- read.table(input.files[window], 
                                    sep=",", 
                                    header=TRUE,
                                    stringsAsFactors=FALSE)
	# transform 'details' column into usable format: max p for Patient_1 and the same for Patient_2
	if(any(!is.na(transmissions.table$details)))
	{
		tmp	<- strsplit(transmissions.table$details,';')
		tmp	<- do.call('rbind',lapply(seq_along(tmp), function(i){
							z	<- data.table(ROW=i, ID=NA_character_, P=NA_character_)
							if(!any(is.na(tmp[[i]])))
								z<- data.table(ROW=i, ID=gsub('(.*):[0-9]*\\.?[0-9]*','\\1',tmp[[i]]), P=gsub('.*:([0-9]*\\.?[0-9]*)','\\1',tmp[[i]]))
							z
						} ))
		set(tmp, NULL, 'P', tmp[, as.numeric(P)])
		tmp	<- subset(tmp, !is.na(ID))[, list(P=max(P)), by=c('ROW','ID')]
		tmp	<- tmp[, list(Patient_1=ID[1], Patient_2=ID[2], Patient_1_P=P[1], Patient_2_P=P[2]), by='ROW']
		transmissions.table$ROW	<- seq_len(nrow(transmissions.table))
		transmissions.table		<- merge(transmissions.table, subset(tmp, select=c(ROW, Patient_1_P, Patient_2_P)), by="ROW", all=1)
		transmissions.table$ROW	<- NULL
	}
	else
	{
		transmissions.table$Patient_1_P <- NA_real_
		transmissions.table$Patient_2_P <- NA_real_
	}
	
  # Not used right now
  
  # sib.dist.table <- read.table(paste("SibDists_",start,"_to_",end,".csv",sep=""), 
  #                              sep=",", 
  #                              header=TRUE,
  #                              stringsAsFactors=FALSE)
  
  # if(nrow(sib.dist.table)>0){
  #   for(row in seq(1,nrow(sib.dist.table))){
  #     if(!is.na(row.match(sib.dist.table[row,1:2], transmissions.ever.suggested))){
  #       sib.window.vector[row.match(sib.dist.table[row,1:2], transmissions.ever.suggested)] <- sib.dist.table[row,3]
  #     } else if(!is.na(row.match(rev(sib.dist.table[row,1:2]), transmissions.ever.suggested))){
  #       sib.window.vector[row.match(rev(sib.dist.table[row,1:2]), transmissions.ever.suggested)] <- sib.dist.table[row,3]
  #     }
  #   }
  # }
  # 
	
  #	OR: replacing FOR loop with a faster merge							
  #transmissions.table <- as.data.table(transmissions.table)
  	
  
  if(nrow(transmissions.table)>0){
    for(row in seq(1,nrow(transmissions.table))){
      # if the row exists
      if(!is.na(row.match(transmissions.table[row,1:2], relationships.ever.suggested))){
        #whatever the third column says, things are the right way round        
        window.vector[row.match(transmissions.table[row,1:2], relationships.ever.suggested)] <- transmissions.table[row,3] 		
		dist.vector[row.match(transmissions.table[row,1:2], relationships.ever.suggested)] <- transmissions.table[row,'Distance']
		p1.vector[row.match(transmissions.table[row,1:2], relationships.ever.suggested)] <- transmissions.table[row,'Patient_1_P']
		p2.vector[row.match(transmissions.table[row,1:2], relationships.ever.suggested)] <- transmissions.table[row,'Patient_2_P']
		
      } else if(!is.na(row.match(rev(transmissions.table[row,1:2]), relationships.ever.suggested))){
		
		dist.vector[row.match(rev(transmissions.table[row,1:2]), relationships.ever.suggested)] <- transmissions.table[row,'Distance']
		#p's will be the wrong way round
		p1.vector[row.match(rev(transmissions.table[row,1:2]), relationships.ever.suggested)] <- transmissions.table[row,'Patient_2_P']
		p2.vector[row.match(rev(transmissions.table[row,1:2]), relationships.ever.suggested)] <- transmissions.table[row,'Patient_1_P']
		#transmissions will be the wrong way round        
        if(transmissions.table[row,3] == "anc"){
          window.vector[row.match(rev(transmissions.table[row,1:2]), relationships.ever.suggested)] <- "desc" 
        } else if (transmissions.table[row,3] == "desc"){
          window.vector[row.match(rev(transmissions.table[row,1:2]), relationships.ever.suggested)] <- "anc" 
        } else if (transmissions.table[row,3] == "intAnc"){
          window.vector[row.match(rev(transmissions.table[row,1:2]), relationships.ever.suggested)] <- "intDesc" 
        } else if (transmissions.table[row,3] == "intDesc"){
          window.vector[row.match(rev(transmissions.table[row,1:2]), relationships.ever.suggested)] <- "intAnc" 
        } else {
          window.vector[row.match(rev(transmissions.table[row,1:2]), relationships.ever.suggested)] <- transmissions.table[row,3]
        }
      }
      if(row %% 100 == 0){
        cat("Processed ",row," out of ", nrow(transmissions.table), " rows\n", sep="")
      }
    }
  }
  if(first){
    first <- FALSE
    window.table <- window.vector
	dist.table <- dist.vector
	p1.table <- p1.vector
	p2.table <- p2.vector
  } else {
    window.table <- cbind(window.table, window.vector)
	dist.table <- cbind(dist.table, dist.vector)
	p1.table <- cbind(p1.table, p1.vector)
	p2.table <- cbind(p2.table, p2.vector)
  }
  # if(first.sib){
  #   first.sib <- FALSE
  #   out.sib <- sib.window.vector
  # } else {
  #   out.sib <- cbind(out.sib, sib.window.vector) 
  # }
}

window.table <- cbind(relationships.ever.suggested, window.table)
dist.table <- cbind(relationships.ever.suggested, dist.table)
p1.table <- cbind(relationships.ever.suggested, p1.table)
p2.table <- cbind(relationships.ever.suggested, p2.table)
#out.sib <- cbind(transmissions.ever.suggested, out.sib)

colnames(window.table) <- c("pat.1", "pat.2", basename(input.files))
colnames(dist.table) <- c("pat.1", "pat.2", basename(input.files))
colnames(p1.table) <- c("pat.1", "pat.2", basename(input.files))
colnames(p2.table) <- c("pat.1", "pat.2", basename(input.files))
#colnames(out.sib) <- c("pat.1", "pat.2", paste("window.start.",seq(800,9050,by=250), sep=""))
#write.table(out.sib, file="quickout_sib.csv", sep=",", col.names=NA)
# out <- read.table("quickout.csv", sep=",", stringsAsFactors = F, header = T)
# out.sib <- read.table("quickout_sib.csv", sep=",", stringsAsFactors = F, header = T)

#first <- TRUE
#for(row in seq(1, nrow(window.table))){
  
#  row.list <- list()
  
#  row.list["anc"] <- 0
#  row.list["desc"] <- 0
#  row.list["int"] <- 0
#  row.list["cher"] <- 0
#  row.list["unint"] <- 0
  
#  if(give.denoms){
#    first.patient <- window.table$pat.1[row]
#    first.present <- reads.table[which(reads.table$patient==first.patient & reads.table$present),]
#    second.patient <- window.table$pat.2[row]
#    second.present <- reads.table[which(reads.table$patient==second.patient & reads.table$present),]
#    
#    denominator <- length(intersect(first.present$window.start, second.present$window.start))
#  }
  
#  for(col in seq(3, ncol(window.table))){
#    if(!is.na(window.table[row, col])){
#      if(allow.splits){
#        if(as.character(window.table[row,col])=="trueInt"){
#          row.list[["int"]] <- row.list[["int"]] + 1
#        } else if(as.character(window.table[row,col]) %in% c("anc", "intAnc")) {
#          row.list[["anc"]] <- row.list[["anc"]] + 1
#        } else if(as.character(window.table[row,col]) %in% c("desc", "intDesc")) {
#          row.list[["desc"]] <- row.list[["desc"]] + 1
#        } else if(as.character(window.table[row,col])=="cher") {
#          row.list[["cher"]] <- row.list[["cher"]] + 1
#        } else if(as.character(window.table[row,col])=="unint") {
#          row.list[["unint"]] <- row.list[["unint"]] + 1
#        }
#      } else {
#        if(as.character(window.table[row,col]) %in% c("trueInt", "intAnc", "intDesc")){
#          row.list[["int"]] <- row.list[["int"]] + 1
#        } else if(as.character(window.table[row,col])=="anc") {
#          row.list[["anc"]] <- row.list[["anc"]] + 1
#        } else if(as.character(window.table[row,col])=="desc") {
#          row.list[["desc"]] <- row.list[["desc"]] + 1
#        } else if(as.character(window.table[row,col])=="cher") {
#          row.list[["cher"]] <- row.list[["cher"]] + 1
#        } else if(as.character(window.table[row,col])=="unint") {
#          row.list[["unint"]] <- row.list[["unint"]] + 1
#        }
#      }
#    }
#  }
  
#  cher.row <- c(window.table[row,1], window.table[row,2], row.list[["cher"]], "cher", row.list[["anc"]]+row.list[["desc"]])
#  unint.row <- c(window.table[row,1], window.table[row,2], row.list[["unint"]], "unint", row.list[["anc"]]+row.list[["desc"]])
#  anc.row <- c(window.table[row,1], window.table[row,2], row.list[["anc"]], "anc", row.list[["anc"]]+row.list[["desc"]])
#  desc.row <- c(window.table[row,2], window.table[row,1], row.list[["desc"]], "anc", row.list[["anc"]]+row.list[["desc"]])
#  int.row <- c(window.table[row,1], window.table[row,2], row.list[["int"]], "int", row.list[["anc"]]+row.list[["desc"]])
  
#  if(give.denoms){
#    cher.row <- c(cher.row, paste(row.list[["cher"]], "/", denominator, sep=""))
#    unint.row <- c(unint.row, paste(row.list[["unint"]], "/", denominator, sep=""))
#    anc.row <- c(anc.row, paste(row.list[["anc"]], "/", denominator, sep=""))
#    desc.row <- c(desc.row, paste(row.list[["desc"]], "/", denominator, sep=""))
#    int.row <- c(int.row, paste(row.list[["int"]], "/", denominator, sep=""))
# }
  
  # intAnc.row <- c(out[row,1], out[row,2], row.list[["intAnc"]], "PPU", row.list[["anc"]]+row.list[["desc"]]+row.list[["intAnc"]]+row.list[["intDesc"]]+row.list[["trueInt"]], denominator, paste(row.list[["intAnc"]], "/", denominator, sep=""), NA)
  # intDesc.row <- c(out[row,2], out[row,1], row.list[["intDesc"]], "PPU", row.list[["anc"]]+row.list[["desc"]]+row.list[["intAnc"]]+row.list[["intDesc"]]+row.list[["trueInt"]], denominator, paste(row.list[["intDesc"]], "/", denominator, sep=""), NA)
  # trueInt.row <- c(out[row,1], out[row,2], row.list[["trueInt"]], "PPE", row.list[["anc"]]+row.list[["desc"]]+row.list[["intAnc"]]+row.list[["intDesc"]]+row.list[["trueInt"]], denominator, paste(row.list[["trueInt"]], "/", denominator, sep=""), NA)
  
#  new.rows <- data.frame(rbind(cher.row, unint.row, anc.row, desc.row, int.row), stringsAsFactors = F)
  
#  colnames(new.rows) <- c("pat.1", "pat.2", "windows", "type", "total.trans")
#  if(give.denoms){
#    colnames(new.rows)[6] <- "fraction"
#  }
  
  # sib.row.2 <- c(out[row,1], out[row,2], row.list.2[["sib"]], "sib", row.list.2[["anc"]]+row.list.2[["desc"]]+row.list.2[["int"]], denominator, paste(row.list.2[["sib"]], "/", denominator, sep=""), mean.sib)
  # anc.row.2 <- c(out[row,1], out[row,2], row.list.2[["anc"]], "trans", row.list.2[["anc"]]+row.list.2[["desc"]]+row.list.2[["int"]], denominator, paste(row.list.2[["anc"]], "/", denominator, sep=""), NA)
  # desc.row.2 <- c(out[row,2], out[row,1], row.list.2[["desc"]], "trans", row.list.2[["anc"]]+row.list.2[["desc"]]+row.list.2[["int"]], denominator, paste(row.list.2[["desc"]], "/", denominator, sep=""), NA)
  # int.row.2 <- c(out[row,1], out[row,2], row.list.2[["int"]], "int", row.list.2[["anc"]]+row.list.2[["desc"]]+row.list.2[["int"]], denominator, paste(row.list.2[["int"]], "/", denominator, sep=""), NA)
  # 
  # new.rows.2 <- data.frame(rbind(sib.row.2, anc.row.2, desc.row.2, int.row.2), stringsAsFactors = F)
  
  #  colnames(new.rows.2) <- c("pat.1", "pat.2", "windows", "type", "total.trans", "present", "fraction", "mean.sib")
  
#  if(first){
#    first <- FALSE
#    new.out <- new.rows[which(new.rows$windows>0),]
    #    new.out.2 <- new.rows.2[which(new.rows.2$windows>0),]
#  } else {
#    new.out <- rbind(new.out, new.rows[which(new.rows$windows>0),])
    #    new.out.2 <- rbind(new.out.2, new.rows.2[which(new.rows.2$windows>0),])
#  }
#  
#  if(row %% 100 == 0){
#    cat("Written ",row," out of ",nrow(window.table)," rows\n",sep="")
#  }
#}
#new.out$windows <- as.numeric(new.out$windows)
#new.out$total.trans <- as.numeric(new.out$total.trans)
#write.table(new.out[which(new.out$total.trans>=min.threshold),], file=output.file, row.names = F, sep=",", quote=FALSE)
#
#	OR: 	I vectorized the code above this line-- hopefulling speeding things up a bit
#			suggest to remove and replace with the lines below
#

cat("Making per-window output table...\n")
#	get tables into long format
window.table[] 	<- lapply(window.table, as.character)
window.table	<- melt(window.table, id.vars=c('pat.1','pat.2'), value.name="TYPE")
dist.table[] 	<- lapply(dist.table, as.character)
dist.table		<- melt(dist.table, id.vars=c('pat.1','pat.2'), value.name="PATRISTIC_DISTANCE")
p1.table[] 		<- lapply(p1.table, as.character)
p1.table		<- melt(p1.table, id.vars=c('pat.1','pat.2'), value.name="ID1_P")
p2.table[] 		<- lapply(p2.table, as.character)
p2.table		<- melt(p2.table, id.vars=c('pat.1','pat.2'), value.name="ID2_P")

window.table <- merge(window.table, dist.table, by=c('pat.1','pat.2','variable'))
window.table <- merge(window.table, p1.table, by=c('pat.1','pat.2','variable'))
window.table <- merge(window.table, p2.table, by=c('pat.1','pat.2','variable'))
dist.table <- NULL
p1.table <- NULL
p2.table <- NULL
# convert "anc" and "desc" to "trans" depending on direction

window.table[which(window.table$TYPE=="desc"), c("pat.1", "pat.2")] <- window.table[which(window.table$TYPE=="desc"), c("pat.2", "pat.1")]
window.table[which(window.table$TYPE=="anc"), "TYPE"] <- "trans"
window.table[which(window.table$TYPE=="desc"), "TYPE"] <- "trans"

window.table[which(window.table$TYPE=="intDesc"), c("pat.1", "pat.2")] <- window.table[which(window.table$TYPE=="intDesc"), c("pat.2", "pat.1")]
window.table[which(window.table$TYPE=="intAnc"), "TYPE"] <- "intTrans"
window.table[which(window.table$TYPE=="intDesc"), "TYPE"] <- "intTrans"

window.table	<- as.data.table(window.table)
setnames(window.table, c('variable'), c('SOURCE_FILE'))
#	add window coordinates
window.table[, W_FROM:= window.table[,as.integer(gsub(prefix.wfrom,'',regmatches(SOURCE_FILE, regexpr(paste(prefix.wfrom,'[0-9]+',sep=''),SOURCE_FILE))))]]
window.table[, W_TO:= window.table[, as.integer(gsub(prefix.wto,'',regmatches(SOURCE_FILE, regexpr(paste(prefix.wto,'[0-9]+',sep=''),SOURCE_FILE))))]]	
#	change type name
if(allow.splits)
{
	set(window.table, window.table[, which(TYPE=='trueInt')], 'TYPE', 'int')
	set(window.table, window.table[, which(TYPE %in% c("intTrans"))], 'TYPE', 'trans')	
}
if(!allow.splits)
{
	set(window.table, window.table[, which(TYPE %in% c("trueInt", "intTrans"))], 'TYPE', 'int')	
}
#	set to numeric
set(window.table, NULL, 'PATRISTIC_DISTANCE', window.table[, as.numeric(PATRISTIC_DISTANCE)])
set(window.table, NULL, 'ID1_P', window.table[, as.numeric(ID1_P)])
set(window.table, NULL, 'ID2_P', window.table[, as.numeric(ID2_P)])

#	OR: calculate #unique reads and #reads for each window 
reads.table	<- as.data.table(reads.table)
setnames(reads.table, c("window.start","window.end"), c("W_FROM","W_TO"))
dp			<-  reads.table[,{
			#	combine by window only for present patients to avoid large mem
			dp			<- data.table(pat.1=NA_character_, pat.2=NA_character_)
			if(length(which(present))>1)
			{
				dp			<- t( combn( patient[present], 2 ) )
				colnames(dp)<- c('pat.1','pat.2')
				dp			<- as.data.table(dp)	
			}
			dp						
		}, by=c('W_FROM','W_TO')]
dp			<- subset(dp, !is.na(pat.1) & !is.na(pat.2))
tmp			<- subset(reads.table, present, c(W_FROM, patient, reads, leaves))
setnames(tmp, c("patient","reads","leaves"), c("pat.1","pat.1_reads","pat.1_leaves"))
dp			<- merge(dp, tmp, by=c('W_FROM','pat.1'))
setnames(tmp, c("pat.1","pat.1_leaves","pat.1_reads"), c("pat.2","pat.2_leaves","pat.2_reads"))
dp			<- merge(dp, tmp, by=c('W_FROM','pat.2'))
#	OR: add window.table to this table of pairs. keep track of direction
tmp			<- subset(window.table, !is.na(TYPE), select=c('pat.1','pat.2','W_FROM','TYPE','PATRISTIC_DISTANCE','ID1_P','ID2_P'))
set(tmp, tmp[, which(TYPE=='trans')], 'TYPE', 'trans_12')
#	merge (1,2) with dp
dp			<- merge(dp, tmp, by=c('W_FROM','pat.1','pat.2'),all.x=1)
set(tmp, tmp[, which(TYPE=='trans_12')], 'TYPE', 'trans_21')
setnames(tmp, c('pat.1','pat.2','TYPE','PATRISTIC_DISTANCE','ID1_P','ID2_P'), c('pat.2','pat.1','TYPE2','PATRISTIC_DISTANCE2','ID1_P2','ID2_P2'))
#	merge (2,1) with dp -- we now have all trm assignments to the unordered pair in dp either in col TYPE or TYPE2
dp			<- merge(dp, tmp, by=c('W_FROM','pat.1','pat.2'),all.x=1)
#	consolidate assignments
stopifnot(dp[, length(which(!is.na(TYPE) & !is.na(TYPE2)))==0L])
stopifnot(dp[, length(which(!is.na(PATRISTIC_DISTANCE) & !is.na(PATRISTIC_DISTANCE2)))==0L])
tmp			<- dp[, which(!is.na(TYPE2))]
set(dp, tmp, 'TYPE', dp[tmp,TYPE2])
dp[, TYPE2:=NULL]
tmp			<- dp[, which(!is.na(PATRISTIC_DISTANCE2))]
set(dp, tmp, 'PATRISTIC_DISTANCE', dp[tmp,PATRISTIC_DISTANCE2])
dp[, PATRISTIC_DISTANCE2:=NULL]
tmp			<- dp[, which(!is.na(ID1_P2))]
set(dp, tmp, 'ID1_P', dp[tmp,ID1_P2])
dp[, ID1_P2:=NULL]
tmp			<- dp[, which(!is.na(ID2_P2))]
set(dp, tmp, 'ID2_P', dp[tmp,ID2_P2])
dp[, ID2_P2:=NULL]
#	check we have patristic distances for all cherries
stopifnot(dp[, length(which(TYPE=='cher' & is.na(PATRISTIC_DISTANCE)))==0L])
stopifnot(dp[, length(which(TYPE%in%c('cher','unint') & is.na(ID1_P)))==0L])
stopifnot(dp[, length(which(TYPE%in%c('cher','unint') & is.na(ID2_P)))==0L])
#	set disconnected only when both patients present
set(dp, dp[,which(is.na(TYPE))], 'TYPE','disconnected')
setkey(dp, W_FROM, W_TO, pat.1, pat.2)
dp			<- subset(dp, select=c('W_FROM','W_TO','pat.1','pat.2','TYPE','PATRISTIC_DISTANCE','ID1_P','ID2_P','pat.1_leaves','pat.1_reads','pat.2_leaves','pat.2_reads'))
#	write to file
if(!is.null(detailed.output))
{	
	save(dp, file=gsub('\\.csv','\\.rda',detailed.output))
}
#
#	make summary of transmission assignments across all columns
#
window.table	<- subset(window.table, select=c('pat.1','pat.2','W_FROM','W_TO','TYPE','SOURCE_FILE'))
setkey(window.table, pat.1, pat.2, W_FROM)
#
cat("Making summary output table...\n")
ans		<- window.table[, list(windows=length(W_FROM)), by=c('pat.1','pat.2','TYPE')]
#	evaluate total.trans, denominator, and fraction
tmp		<- dp[, list(denominator=length(W_FROM), total.trans=length(W_FROM[TYPE=='trans_12'|TYPE=='trans_21'])), by=c('pat.1','pat.2')]
tmp2	<- copy(tmp)
setnames(tmp2, c('pat.1','pat.2'), c('pat.2','pat.1'))
tmp		<- rbind(tmp, tmp2)
ans		<- merge(subset(ans, !is.na(TYPE)), tmp, by=c('pat.1','pat.2'))
ans[, fraction:=paste(windows,'/',denominator,sep='')]

#if(give.denoms){
#  determine.denominator <- function(a, b){
#    first.present <- reads.table[which(reads.table$patient==a & reads.table$present),]
#    second.present <- reads.table[which(reads.table$patient==b & reads.table$present),]
#    
#    return (length(intersect(first.present$window.start, second.present$window.start)))
#  }
#  
#  ans$denominator <- mapply(determine.denominator, ans$pat.1, ans$pat.2)
#}
#tmp	<- ans[, list(TYPE=TYPE, fraction=paste(windows,'/',denominator,sep='')), by=c('pat.1','pat.2')]
#ans	<- merge(ans, tmp, by=c('pat.1','pat.2','TYPE'))
#	evaluate total.trans

#get.total.trans <- function(a,b){
#  if(nrow(ans[which(ans$pat.1==a & ans$pat.2==b & ans$TYPE=="trans"),])==0 &
#     nrow(ans[which(ans$pat.1==b & ans$pat.2==a & ans$TYPE=="trans"),])==0){
#    return(0)
#  } else if (nrow(ans[which(ans$pat.1==a & ans$pat.2==b & ans$TYPE=="trans"),])==1 &
#             nrow(ans[which(ans$pat.1==b & ans$pat.2==a & ans$TYPE=="trans"),])==1){
#    return(ans[which(ans$pat.1==a & ans$pat.2==b & ans$TYPE=="trans"), windows] + ans[which(ans$pat.1==b & ans$pat.2==a & ans$TYPE=="trans"), windows])
#  } else if (nrow(ans[which(ans$pat.1==a & ans$pat.2==b & ans$TYPE=="trans"),])==1)  {
#    return (ans[which(ans$pat.1==a & ans$pat.2==b & ans$TYPE=="trans"), windows])
#  } else if (nrow(ans[which(ans$pat.1==b & ans$pat.2==a & ans$TYPE=="trans"),])==1){
#    return (ans[which(ans$pat.1==b & ans$pat.2==a & ans$TYPE=="trans"), windows])
#  } else {
#    stop("Too many transmission rows")
#  }
#}

#ans$total.trans <- mapply(get.total.trans, ans$pat.1, ans$pat.2)

setkey(ans, pat.1, pat.2, TYPE)
#	write to file
write.csv(subset(ans, total.trans>=min.threshold), file=output.file, row.names=FALSE, quote=FALSE)


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
