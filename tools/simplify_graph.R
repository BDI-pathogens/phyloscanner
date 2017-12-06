#!/usr/bin/env Rscript

suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(phyloscannerR, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(reshape2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(GGally, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(kimisc, quietly=TRUE, warn.conflicts=FALSE))

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

results <- simplify.summary(input, direction.threshold, total.trees, plot = F)

write.csv(results$simp.table, output.file.name, quote=F, row.names=F)
