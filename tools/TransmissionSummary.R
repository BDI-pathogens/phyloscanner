#!/usr/bin/env Rscript

list.of.packages <- c("prodlim")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

command.line <- T
if(command.line){
  library(argparse)
  
  arg_parser = ArgumentParser(description="Summarise topological relationships suggesting direction of transmission across windows. Outputs a .csv file of relationships between patient IDs. Relationship types: 'anc' indicates a topology that suggests the patient in column 1 infected the patient in column 2, 'sib' one where reads from the patients form sibling clades from an unsampled common ancestor patient (and that neither has any sampled ancestors occuring later in the transmission chain than that common ancestor), 'int' one where reads from both patients are intermingled and the direction of transmission cannot be established. (More documentation to follow.)")
  
  arg_parser$add_argument("-l", "--filesAreLists", action="store_true", default=FALSE, help="If present, arguments specifying input files will be parsed as lists of files separated by colons. If absent, they will be parsed as a string which each file of that type begins with.")
  arg_parser$add_argument("-s", "--summaryFile", action="store", help="The full output file from SummaryStatistics.R; necessary only to identify windows in which no reads are present from each patient. If absent, window counts will be given without denominators.")
  arg_parser$add_argument("-m", "--minThreshold", action="store", type="integer", help="Relationships between two patients will only appear in output if a transmission chain between them appears in at least these many windows (default 1). High numbers are useful for drawing figures in e.g. Cytoscape with few enough arrows to be comprehensible. The script is extremely slow if this is set to 0.")
  arg_parser$add_argument("-p", "--allowSplits", action="store_true", default=FALSE, help="If absent, directionality is only inferred between pairs of patients whose reads are not split; this is more conservative.")
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
  if(is.null(min.threshold)){
    split.threshold <- 1L
  }
  allow.splits <- args$allowSplits
  input.files.name <- args$inputFiles
  
  	
  if(!files.are.lists){
    input.files <- sort(list.files(dirname(input.files.name), pattern=paste(basename(input.files.name),".*LikelyTransmissions.csv$",sep=""), full.names=TRUE))
  } else {
    input.files <- unlist(strsplit(input.files.name, ":"))
  }
  
} else {
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517_clean/")
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools"
  summary.file <- "ss_r_patStatsFull.csv"
  id.file <- "patientIDList.txt"  
  min.threshold <- 0
  allow.splits <- TRUE
  output.file <- "test.csv"  
  input.files <- sort(list.files(getwd(), pattern="LikelyTransmissions_r.*\\.csv"))
  if(0)
  {
    script.dir					<- "/Users/Oliver/git/phylotypes/tools"
    summary.file				<- "/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptoutput/ptyr115_patStatsFull.csv"
    id.file 					<- "/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptoutput/ptyr115_patients.txt"
	input.files.name			<- "/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptoutput/ptyr115_"
    min.threshold				<- 1
    allow.splits 				<- TRUE
    output.file 				<- "/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptoutput/ptyr115_trmStats.csv"  
    input.files 				<- sort(list.files(dirname(input.files.name), pattern=paste(basename(input.files.name),".*LikelyTransmissions.csv$",sep=''), full.names=TRUE))
  }
}

num.windows <- length(input.files)

library(prodlim)
library(gdata)
library(ggplot2)

# Get the denominators

if(!is.null(summary.file)){
  cat("Getting window counts per patient...\n")
  full.summary <- read.csv(summary.file, stringsAsFactors = F)
  reads.table <- full.summary[,c("window.start", "patient", "reads")]
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
  #  sib.window.vector <- rep(NA, nrow(transmissions.ever.suggested))
  
  window.count <- window.count + 1
  
  cat("Opening file: ",input.files[window],"\n", sep="")
  
  transmissions.table <- read.table(input.files[window], 
                                    sep=",", 
                                    header=TRUE,
                                    stringsAsFactors=FALSE)
  
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
  if(nrow(transmissions.table)>0){
    for(row in seq(1,nrow(transmissions.table))){
      # if the row exists
      if(!is.na(row.match(transmissions.table[row,1:2], relationships.ever.suggested))){
        #whatever the third column says, things are the right way round
        
        window.vector[row.match(transmissions.table[row,1:2], relationships.ever.suggested)] <- transmissions.table[row,3] 
      } else if(!is.na(row.match(rev(transmissions.table[row,1:2]), relationships.ever.suggested))){
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
  } else {
    window.table <- cbind(window.table, window.vector) 
  }
  # if(first.sib){
  #   first.sib <- FALSE
  #   out.sib <- sib.window.vector
  # } else {
  #   out.sib <- cbind(out.sib, sib.window.vector) 
  # }
}

window.table <- cbind(relationships.ever.suggested, window.table)
#out.sib <- cbind(transmissions.ever.suggested, out.sib)

colnames(window.table) <- c("pat.1", "pat.2", input.files)
#colnames(out.sib) <- c("pat.1", "pat.2", paste("window.start.",seq(800,9050,by=250), sep=""))

#write.table(out, file="quickout.csv", sep=",", col.names=NA)
#write.table(out.sib, file="quickout_sib.csv", sep=",", col.names=NA)

# out <- read.table("quickout.csv", sep=",", stringsAsFactors = F, header = T)
# out.sib <- read.table("quickout_sib.csv", sep=",", stringsAsFactors = F, header = T)

cat("Making output table...\n")

first <- TRUE

for(row in seq(1, nrow(window.table))){
  
  row.list <- list()
  
  row.list["anc"] <- 0
  row.list["desc"] <- 0
  row.list["int"] <- 0
  row.list["sib"] <- 0
  
  if(give.denoms){
    first.patient <- window.table$pat.1[row]
    first.present <- reads.table[which(reads.table$patient==first.patient & reads.table$present),]
    second.patient <- window.table$pat.2[row]
    second.present <- reads.table[which(reads.table$patient==second.patient & reads.table$present),]
    
    denominator <- length(intersect(first.present$window.start, second.present$window.start))
  }
  
  for(col in seq(3, ncol(window.table))){
    if(!is.na(window.table[row, col])){
      if(allow.splits){
        if(as.character(window.table[row,col])=="trueInt"){
          row.list[["int"]] <- row.list[["int"]] + 1
        } else if(as.character(window.table[row,col]) %in% c("anc", "intAnc")) {
          row.list[["anc"]] <- row.list[["anc"]] + 1
        } else if(as.character(window.table[row,col]) %in% c("desc", "intDesc")) {
          row.list[["desc"]] <- row.list[["desc"]] + 1
        } else if(as.character(window.table[row,col])=="sib"){
          row.list[["sib"]] <- row.list[["sib"]] + 1
        }
      } else {
        if(as.character(window.table[row,col]) %in% c("trueInt", "intAnc", "intDesc")){
          row.list[["int"]] <- row.list[["int"]] + 1
        } else if(as.character(window.table[row,col])=="anc") {
          row.list[["anc"]] <- row.list[["anc"]] + 1
        } else if(as.character(window.table[row,col])=="desc") {
          row.list[["desc"]] <- row.list[["desc"]] + 1
        } else if(as.character(window.table[row,col])=="sib"){
          row.list[["sib"]] <- row.list[["sib"]] + 1
        }
      }
    }
  }
  
  sib.row <- c(window.table[row,1], window.table[row,2], row.list[["sib"]], "sib", row.list[["anc"]]+row.list[["desc"]])
  anc.row <- c(window.table[row,1], window.table[row,2], row.list[["anc"]], "anc", row.list[["anc"]]+row.list[["desc"]])
  desc.row <- c(window.table[row,2], window.table[row,1], row.list[["desc"]], "anc", row.list[["anc"]]+row.list[["desc"]])
  int.row <- c(window.table[row,1], window.table[row,2], row.list[["int"]], "int", row.list[["anc"]]+row.list[["desc"]])
  
  if(give.denoms){
    sib.row <- c(sib.row, paste(row.list[["sib"]], "/", denominator, sep=""))
    anc.row <- c(anc.row, paste(row.list[["anc"]], "/", denominator, sep=""))
    desc.row <- c(desc.row, paste(row.list[["desc"]], "/", denominator, sep=""))
    int.row <- c(int.row, paste(row.list[["int"]], "/", denominator, sep=""))
    
  }
  
  # intAnc.row <- c(out[row,1], out[row,2], row.list[["intAnc"]], "PPU", row.list[["anc"]]+row.list[["desc"]]+row.list[["intAnc"]]+row.list[["intDesc"]]+row.list[["trueInt"]], denominator, paste(row.list[["intAnc"]], "/", denominator, sep=""), NA)
  # intDesc.row <- c(out[row,2], out[row,1], row.list[["intDesc"]], "PPU", row.list[["anc"]]+row.list[["desc"]]+row.list[["intAnc"]]+row.list[["intDesc"]]+row.list[["trueInt"]], denominator, paste(row.list[["intDesc"]], "/", denominator, sep=""), NA)
  # trueInt.row <- c(out[row,1], out[row,2], row.list[["trueInt"]], "PPE", row.list[["anc"]]+row.list[["desc"]]+row.list[["intAnc"]]+row.list[["intDesc"]]+row.list[["trueInt"]], denominator, paste(row.list[["trueInt"]], "/", denominator, sep=""), NA)
  
  new.rows <- data.frame(rbind(sib.row, anc.row, desc.row, int.row), stringsAsFactors = F)
  
  colnames(new.rows) <- c("pat.1", "pat.2", "windows", "type", "total.trans")
  if(give.denoms){
    colnames(new.rows)[6] <- "fraction"
  }
  
  # sib.row.2 <- c(out[row,1], out[row,2], row.list.2[["sib"]], "sib", row.list.2[["anc"]]+row.list.2[["desc"]]+row.list.2[["int"]], denominator, paste(row.list.2[["sib"]], "/", denominator, sep=""), mean.sib)
  # anc.row.2 <- c(out[row,1], out[row,2], row.list.2[["anc"]], "trans", row.list.2[["anc"]]+row.list.2[["desc"]]+row.list.2[["int"]], denominator, paste(row.list.2[["anc"]], "/", denominator, sep=""), NA)
  # desc.row.2 <- c(out[row,2], out[row,1], row.list.2[["desc"]], "trans", row.list.2[["anc"]]+row.list.2[["desc"]]+row.list.2[["int"]], denominator, paste(row.list.2[["desc"]], "/", denominator, sep=""), NA)
  # int.row.2 <- c(out[row,1], out[row,2], row.list.2[["int"]], "int", row.list.2[["anc"]]+row.list.2[["desc"]]+row.list.2[["int"]], denominator, paste(row.list.2[["int"]], "/", denominator, sep=""), NA)
  # 
  # new.rows.2 <- data.frame(rbind(sib.row.2, anc.row.2, desc.row.2, int.row.2), stringsAsFactors = F)
  
  #  colnames(new.rows.2) <- c("pat.1", "pat.2", "windows", "type", "total.trans", "present", "fraction", "mean.sib")
  
  if(first){
    first <- FALSE
    new.out <- new.rows[which(new.rows$windows>0),]
    #    new.out.2 <- new.rows.2[which(new.rows.2$windows>0),]
  } else {
    new.out <- rbind(new.out, new.rows[which(new.rows$windows>0),])
    #    new.out.2 <- rbind(new.out.2, new.rows.2[which(new.rows.2$windows>0),])
  }
  
  if(row %% 100 == 0){
    cat("Written ",row," out of ",nrow(window.table)," rows\n",sep="")
  }
}

new.out$windows <- as.numeric(new.out$windows)
new.out$total.trans <- as.numeric(new.out$total.trans)

write.table(new.out[which(new.out$total.trans>=min.threshold),], file=output.file, row.names = F, sep=",", quote=FALSE)

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
#             temporal <- new.out.2[which(new.out.2$pat.1==early & new.out.2$pat.2==late & new.out.2$type=="trans"),]
#             contratemporal <- new.out.2[which(new.out.2$pat.1==late & new.out.2$pat.2==early & new.out.2$type=="trans"),]
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
#             temporal <- new.out[which(new.out$pat.1==early & new.out$pat.2==late & new.out$type=="MP"),]
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
#             temporal <- new.out.2[which(new.out.2$pat.1==early & new.out.2$pat.2==late & new.out.2$type=="trans"),]
#             contratemporal <- new.out.2[which(new.out.2$pat.1==late & new.out.2$pat.2==early & new.out.2$type=="trans"),]
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
#             temporal <- new.out[which(new.out$pat.1==early & new.out$pat.2==late & new.out$type=="MP"),]
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
