#!/usr/bin/env Rscript

list.of.packages <- c("reshape2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Please run package_install.R to continue\n")
  quit(save="no", status=1)
}

suppressMessages(library(reshape2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))

tmp	<- "Simplifies host relationships for graph visualisation"
arg_parser = ArgumentParser(description=tmp)
arg_parser$add_argument("-dt", "--directionThreshold", action="store", default=0.33, type="double", help="Directed edges appear between hosts in the output if they are related and a direction of transmission is indicated in at least this proportion of windows. If both directions meet the threshold, then the direction is the one with the greater proportion. This should be smaller than the window threshold used to determine whether links appear at all (-swt in phyloscanner_analyse_trees.R, -m in transmission_summary.R)")
arg_parser$add_argument("inputFile", action="store", help="The transmission summary output from a multiple-tree phyloscanner_analyse_trees.R run or transmission_summary.R")
arg_parser$add_argument("outputFile", action="store", help="A .csv file to write the output to.")
arg_parser$add_argument("totalTrees", action="store", type="double", help="The total number of trees in the phyloscanner analysis.")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")

args <- arg_parser$parse_args()

input.file.name       <- args$inputFile
direction.threshold   <- args$directionThreshold
total.trees           <- args$totalTrees
output.file.name      <- args$outputFile

input <- read.csv(input.file.name, stringsAsFactors = F)

done <- rep(FALSE, nrow(input))

for(line in 1:nrow(input)){
  if(!done[line]){
    forwards.rows <- which(input$host.1 == input$host.1[line] & input$host.2 == input$host.2[line])
    backwards.rows <- which(input$host.1 == input$host.2[line] & input$host.2 == input$host.1[line])
    
    done[forwards.rows] <- T
    done[backwards.rows] <- T
    
    input$ancestry[intersect(forwards.rows, which(input$ancestry=="trans"))] <- "trans12"
    input$ancestry[intersect(forwards.rows, which(input$ancestry=="multi_trans"))] <- "multi_trans12"
    
    input$ancestry[intersect(backwards.rows, which(input$ancestry=="trans"))] <- "trans21"
    input$ancestry[intersect(backwards.rows, which(input$ancestry=="multi_trans"))] <- "multi_trans21"
    
    input[backwards.rows,c(1,2)] <- input[backwards.rows,c(2,1)]
  }
}

input.wide <- reshape(input, direction="w", idvar=c("host.1", "host.2"), timevar = "ancestry", v.names = "ancestry.tree.count",  drop=c("fraction"))

input.wide[is.na(input.wide)] <- 0

input.wide$total.equiv <- input.wide$ancestry.tree.count.none + input.wide$ancestry.tree.count.complex
input.wide$total.12 <- input.wide$ancestry.tree.count.trans12
if(!is.null(input.wide$ancestry.tree.count.multi_trans12)){
  input.wide$total.12 <- input.wide$total.12 + input.wide$ancestry.tree.count.multi_trans12
}
input.wide$total.21 <- input.wide$ancestry.tree.count.trans21
if(!is.null(input.wide$ancestry.tree.count.multi_trans21)){
  input.wide$total.21 <- input.wide$total.21+ input.wide$ancestry.tree.count.multi_trans21
}
input.wide$total <- input.wide$total.21 + input.wide$total.12 + input.wide$total.equiv

dir <- input.wide$total.12 >= direction.threshold*total.trees | input.wide$total.21 >= direction.threshold*total.trees

input.wide$arrow[!dir] <- "none"

if(length(which(dir)>0)){
  input.wide$arrow[dir] <- sapply(which(dir), function(x)  if(input.wide$total.12[x]>input.wide$total.21[x]) "forwards" else "backwards")
  input.wide$label[dir] <- paste(round(pmax(input.wide[dir, "total.12"],input.wide[dir, "total.21"])/total.trees, 2),"/", round(input.wide[dir, "total"]/total.trees, 2), sep="")
}

input.wide$label[!dir] <- as.character(round(input.wide[!dir, "total"]/total.trees, 2))

write.csv(input.wide[,c(1,2,13,14)], output.file.name, quote=F, row.names=F)
