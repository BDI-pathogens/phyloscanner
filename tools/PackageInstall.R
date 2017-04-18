every.package <- c("argparse", 
                   "ape", 
                   "phangorn", 
                   "data.table", 
                   "ggplot2", 
                   "ff",
                   "phytools", 
                   "dplyr", 
                   "reshape", 
                   "dtplyr", 
                   "gtable", 
                   "grid", 
                   "gridExtra", 
                   "RColorBrewer", 
                   "scales", 
                   "pegas",
                   "prodlim",
                   "reshape2")

new.packages <- every.package[!(every.package %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

if(!("ggtree" %in% installed.packages()[,"Package"])){
  source("https://bioconductor.org/biocLite.R")
  biocLite("ggtree")
}