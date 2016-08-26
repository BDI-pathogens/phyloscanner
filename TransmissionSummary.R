# summarises transmission information

command.line <- F
if(command.line){
  library(argparse)
  
  arg_parser = ArgumentParser(description="Summarise topological relationships suggesting direction of transmission across windows.")
  
  arg_parser$add_argument("-l", "--filesAreLists", action="store_true", default=FALSE, help="If present, arguments specifying input files will be parsed as lists of files separated by colons. If absent, they will be parsed as a string which each file of that type begins with.")
  arg_parser$add_argument("-s", "--summaryFile", action="store", help="The full output file from SummaryStatistics.R; necessary only to identify windows in which no reads are present from each patient. If absent, window counts will be given without denominators.")
  arg_parser$add_argument("-w", "--windows", action="store", help="The window in the genome which each tree file, in the same order as the tree files (user-specified if -l is present, in the order provided by the file system if now), is constructed from. In each window this is given as n-m where n is the start and m the end; coordinates for each window are separated by a colon. If not given, the script will attempt to obtain these from the file names.")
  arg_parser$add_argument("-m", "--minThreshold", action="store", default=1, type="integer", help="Integer; relationships will only appear if directionality between these patients appears on at least these many windows (default=1)")
  arg_parser$add_argument("-p", "--allowSplits", action="store_true", default=FALSE,help="If absent, directionality is only inferred between pairs of patients whose reads are not split; this is more conservative.")
  arg_parser$add_argument("idFile", action="store", help="A file containing a list of the IDs of all the patients to calculate and display statistics for.")
  arg_parser$add_argument("inputFiles", action="store", help="Either (if -l is present) a list of all input files, separated by colons, or (if not) a single string that begins every input file name.")
  
  
  args <- arg_parser$parse_args()
  
  files.are.lists <- args$filesAreLists
  summary.file <- args$summaryFile
  id.file <- args$idFile
  
  min.threshold <- args$minThreshold
  allow.splits <- args$allowSplits
  
  if(!files.are.lists){
    input.files <- sort(list.files(getwd(), pattern=paste(tree.file.name,".*\\.tree",sep="")))
  } else {
    input.files <- unlist(strsplit(tree.file.name, ":"))
  }
  
  if(!is.null(args$windows)){
    windows <- unlist(strsplit(args$windows, ":"))
    if(length(input.files)!=length(windows)){
      stop("Number of input files and number of windows differ")
    }
  } else {
    windows <- NULL
  }
} else {
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517_clean/")
  
  summary.file <- "ss_c_patStatsFull.csv"
  id.file <- "patientIDList.txt"
  
  min.threshold <- 6
  allow.splits <- TRUE
  
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517/")
  input.files <- sort(list.files(getwd(), pattern="LikelyTransmissions.*\\.csv"))
  
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517_clean/")
}

num.windows <- length(input.files)

library(prodlim)
#library(lubridate)
library(gdata)
library(ggplot2)

# Find the window numbers

window.starts <- vector()
window.middles <- vector()
window.ends <- vector()

setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517/")

for(window.no in seq(1, length(tree.files))){
  input.file.name <- input.files[window.no]
  
  if(!is.null(windows)){
    window <- windows[window.no]
    start <- as.numeric(unlist(strsplit(window, "-"))[1])
    end <- as.numeric(unlist(strsplit(window, "-"))[2])
    middle <- (start+end)/2
  } else {
    split.name <- strsplit(input.file.name, "[_\\.]")[[1]]
    start <- as.numeric(split.name[length(split.name)-3])
    end <- as.numeric(split.name[length(split.name)-1])
    middle <- (start+end)/2
  }
  window.starts <- c(window.starts, start)
  window.middles <- c(window.middles, middle)
  window.ends <- c(window.ends, end)
}

setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517_clean/")

# Get the denominators

full.summary <- read.csv(summary.file, stringsAsFactors = F)
reads.table <- full.summary[,c("window.start", "patient", "reads")]
reads.table$present <- reads.table$reads>0

# Read in the IDs. Remove duplicates. Shuffle their order if desired.
ids <- scan(id.file, what="", sep="\n", quiet=TRUE, skip = 1)

num.ids <- length(ids)
if (num.ids == 0) {
  cat(paste("No IDs found in ", id.file, ". Quitting.\n", sep=""))
  quit("no", 1)
}
patient.ids <- unique(ids)

parent.list <- vector("list", length(patient.ids))

names(parent.list) <- patient.ids

for(patient.id in patient.ids){
  parent.list[[patient.id]] <- rep(NA, num.windows)
}

window.count <- 0

transmissions.ever.suggested <- vector()

setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517/")

for(window in seq(1, num.windows)){
  window.count <- window.count + 1
  
  transmissions.table <- read.table(input.files[window], sep=",", header=TRUE, stringsAsFactors=FALSE)
  
  if(nrow(transmissions.table)>0){
    for(row in seq(1,nrow(transmissions.table))){
      if(transmissions.table$Relationship[row] %in% c("anc", "desc", "intAnc", "intDesc", "trueInt")){
        
        pat.1 <- transmissions.table[row,1]
        pat.2 <- transmissions.table[row,2]
        
        #if not ignoring split patients
        
        #           pat.1 <- substr(pat.1, 1, 9)
        #           pat.2 <- substr(pat.2, 1, 9)
        
        transmissions.ever.suggested <- rbind(transmissions.ever.suggested, c(pat.1, pat.2))
      } 
    }
  }
}

transmissions.ever.suggested <- as.data.frame(unique(transmissions.ever.suggested), stringsAsFactors=FALSE)

reversed.t.e.s <- as.data.frame(cbind(transmissions.ever.suggested[,2],transmissions.ever.suggested[,1]), stringsAsFactors = FALSE)

matches <- row.match(transmissions.ever.suggested, reversed.t.e.s)

to.go <- rep(FALSE, nrow(transmissions.ever.suggested))

for(comp in seq(1, length(matches))){
  if(!(is.na(matches[comp])) & (comp < matches[comp])){
    to.go[comp] <- TRUE
  }
}

transmissions.ever.suggested <- transmissions.ever.suggested[which(!to.go), ]

window.count <- 0

first <-TRUE
#first.sib <- TRUE

for(window in seq(1, num.windows)){
  window.vector <- rep(NA, nrow(transmissions.ever.suggested))
#  sib.window.vector <- rep(NA, nrow(transmissions.ever.suggested))
  
  window.count <- window.count + 1

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
      if(!is.na(row.match(transmissions.table[row,1:2], transmissions.ever.suggested))){
        #whatever the third column says, things are the right way round
        
        window.vector[row.match(transmissions.table[row,1:2], transmissions.ever.suggested)] <- transmissions.table[row,3] 
      } else if(!is.na(row.match(rev(transmissions.table[row,1:2]), transmissions.ever.suggested))){
        #transmissions will be the wrong way round
        
        if(transmissions.table[row,3] == "anc"){
          window.vector[row.match(rev(transmissions.table[row,1:2]), transmissions.ever.suggested)] <- "desc" 
        } else if (transmissions.table[row,3] == "desc"){
          window.vector[row.match(rev(transmissions.table[row,1:2]), transmissions.ever.suggested)] <- "anc" 
        } else if (transmissions.table[row,3] == "intAnc"){
          window.vector[row.match(rev(transmissions.table[row,1:2]), transmissions.ever.suggested)] <- "intDesc" 
        } else if (transmissions.table[row,3] == "intDesc"){
          window.vector[row.match(rev(transmissions.table[row,1:2]), transmissions.ever.suggested)] <- "intAnc" 
        } else {
          window.vector[row.match(rev(transmissions.table[row,1:2]), transmissions.ever.suggested)] <- transmissions.table[row,3]
        }
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

window.table <- cbind(transmissions.ever.suggested, window.table)
#out.sib <- cbind(transmissions.ever.suggested, out.sib)

colnames(window.table) <- c("pat.1", "pat.2", input.files)
#colnames(out.sib) <- c("pat.1", "pat.2", paste("window.start.",seq(800,9050,by=250), sep=""))

#write.table(out, file="quickout.csv", sep=",", col.names=NA)
#write.table(out.sib, file="quickout_sib.csv", sep=",", col.names=NA)

# out <- read.table("quickout.csv", sep=",", stringsAsFactors = F, header = T)
# out.sib <- read.table("quickout_sib.csv", sep=",", stringsAsFactors = F, header = T)

first <- TRUE

for(row in seq(1, nrow(window.table))){
  row.list <- list()
  
  row.list["anc"] <- 0
  row.list["desc"] <- 0
  row.list["int"] <- 0
  row.list["sib"] <- 0
  
  for(col in seq(3, ncol(window.table))){
    if(allow.splits){
      if(as.character(window.table[row,col])=="trueInt"){
        row.list[["int"]] <- row.list[["int"]] + 1
      } else if(as.character(window.table[row,col]) %in% c("anc", "intAnc")) {
        row.list[["anc"]] <- row.list[["anc"]] + 1
      } else if(as.character(window.table[row,col]) %in% c("desc", "intDesc")) {
        row.list[["desc"]] <- row.list["desc"]] + 1
      } else if(as.character(window.table[row,col])=="sib"){
        row.list[["sib"]] <- row.list.["sib"]] + 1
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


  denominator <- sum(window.rows[1,2:35] & window.rows[2,2:35])
  
  
  sib.row <- c(out[row,1], out[row,2], row.list[["sib"]], "MM", row.list[["anc"]]+row.list[["desc"]]+row.list[["intAnc"]]+row.list[["intDesc"]]+row.list[["trueInt"]], denominator, paste(row.list[["sib"]], "/", denominator, sep=""), mean.sib)
  anc.row <- c(out[row,1], out[row,2], row.list[["anc"]], "MP", row.list[["anc"]]+row.list[["desc"]]+row.list[["intAnc"]]+row.list[["intDesc"]]+row.list[["trueInt"]], denominator, paste(row.list[["anc"]], "/", denominator, sep=""), NA)
  desc.row <- c(out[row,2], out[row,1], row.list[["desc"]], "MP", row.list[["anc"]]+row.list[["desc"]]+row.list[["intAnc"]]+row.list[["intDesc"]]+row.list[["trueInt"]], denominator, paste(row.list[["desc"]], "/", denominator, sep=""), NA)
  intAnc.row <- c(out[row,1], out[row,2], row.list[["intAnc"]], "PPU", row.list[["anc"]]+row.list[["desc"]]+row.list[["intAnc"]]+row.list[["intDesc"]]+row.list[["trueInt"]], denominator, paste(row.list[["intAnc"]], "/", denominator, sep=""), NA)
  intDesc.row <- c(out[row,2], out[row,1], row.list[["intDesc"]], "PPU", row.list[["anc"]]+row.list[["desc"]]+row.list[["intAnc"]]+row.list[["intDesc"]]+row.list[["trueInt"]], denominator, paste(row.list[["intDesc"]], "/", denominator, sep=""), NA)
  trueInt.row <- c(out[row,1], out[row,2], row.list[["trueInt"]], "PPE", row.list[["anc"]]+row.list[["desc"]]+row.list[["intAnc"]]+row.list[["intDesc"]]+row.list[["trueInt"]], denominator, paste(row.list[["trueInt"]], "/", denominator, sep=""), NA)
  
  new.rows <- data.frame(rbind(sib.row, anc.row, desc.row, intAnc.row, intDesc.row, trueInt.row), stringsAsFactors = F)
  
  colnames(new.rows) <- c("pat.1", "pat.2", "windows", "type", "total.trans", "present", "fraction", "mean.sib")
  
  sib.row.2 <- c(out[row,1], out[row,2], row.list.2[["sib"]], "sib", row.list.2[["anc"]]+row.list.2[["desc"]]+row.list.2[["int"]], denominator, paste(row.list.2[["sib"]], "/", denominator, sep=""), mean.sib)
  anc.row.2 <- c(out[row,1], out[row,2], row.list.2[["anc"]], "trans", row.list.2[["anc"]]+row.list.2[["desc"]]+row.list.2[["int"]], denominator, paste(row.list.2[["anc"]], "/", denominator, sep=""), NA)
  desc.row.2 <- c(out[row,2], out[row,1], row.list.2[["desc"]], "trans", row.list.2[["anc"]]+row.list.2[["desc"]]+row.list.2[["int"]], denominator, paste(row.list.2[["desc"]], "/", denominator, sep=""), NA)
  int.row.2 <- c(out[row,1], out[row,2], row.list.2[["int"]], "int", row.list.2[["anc"]]+row.list.2[["desc"]]+row.list.2[["int"]], denominator, paste(row.list.2[["int"]], "/", denominator, sep=""), NA)
  
  new.rows.2 <- data.frame(rbind(sib.row.2, anc.row.2, desc.row.2, int.row.2), stringsAsFactors = F)
  
  colnames(new.rows.2) <- c("pat.1", "pat.2", "windows", "type", "total.trans", "present", "fraction", "mean.sib")
  
  if(first){
    first <- FALSE
    new.out <- new.rows[which(new.rows$windows>0),]
    new.out.2 <- new.rows.2[which(new.rows.2$windows>0),]
  } else {
    new.out <- rbind(new.out, new.rows[which(new.rows$windows>0),])
    new.out.2 <- rbind(new.out.2, new.rows.2[which(new.rows.2$windows>0),])
  }
  
}

new.out$windows <- as.numeric(new.out$windows)
new.out$total.trans <- as.numeric(new.out$total.trans)
new.out$mean.sib <- as.numeric(new.out$mean.sib)

new.out$rel.width = new.out$windows/34
new.out$rel.width[which(new.out$type=="MM")] <- 1-((new.out$mean.sib[which(new.out$type=="MM")]/max(new.out$mean.sib[!is.na(new.out$mean.sib)])))

new.out$arrow.label <- new.out$fraction
new.out$arrow.label[which(new.out$type=="MM")] <- format(new.out$mean.sib[which(new.out$type=="MM")], scientific = T, digits = 3)

new.out.2$windows <- as.numeric(new.out.2$windows)
new.out.2$total.trans <- as.numeric(new.out.2$total.trans)
new.out.2$mean.sib <- as.numeric(new.out.2$mean.sib)

new.out.2$rel.width = new.out.2$windows/34
new.out.2$rel.width[which(new.out$type=="sib")] <- 1-((new.out.2$mean.sib[which(new.out$type=="sib")]/max(new.out.2$mean.sib[!is.na(new.out.2$mean.sib)])))
new.out.2$arrow.label <- new.out.2$fraction
new.out.2$arrow.label[which(new.out.2$type=="sib")] <- format(new.out.2$mean.sib[which(new.out.2$type=="sib")], scientific = T, digits = 3))


write.table(new.out, file="transSummaryRS.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=2),], file="transSummaryRS_2.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=3),], file="transSummaryRS_3.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=4),], file="transSummaryRS_4.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=5),], file="transSummaryRS_5.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=6),], file="transSummaryRS_6.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=7),], file="transSummaryRS_7.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=10),], file="transSummaryRS_10.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=16),], file="transSummaryRS_16.csv", row.names = F, sep=",", quote=FALSE)

write.table(new.out.2, file="transSummaryRS_merged.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out.2[which(new.out.2$total.trans>=2),], file="transSummaryRS_merged_2.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out.2[which(new.out.2$total.trans>=3),], file="transSummaryRS_merged_3.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out.2[which(new.out.2$total.trans>=4),], file="transSummaryRS_merged_4.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out.2[which(new.out.2$total.trans>=5),], file="transSummaryRS_merged_5.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out.2[which(new.out.2$total.trans>=6),], file="transSummaryRS_merged_6.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out.2[which(new.out.2$total.trans>=9),], file="transSummaryRS_merged_9.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out.2[which(new.out.2$total.trans>=7),], file="transSummaryRS_merged_7.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out.2[which(new.out.2$total.trans>=10),], file="transSummaryRS_merged_10.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out.2[which(new.out.2$total.trans>=16),], file="transSummaryRS_merged_16.csv", row.names = F, sep=",", quote=FALSE)

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
