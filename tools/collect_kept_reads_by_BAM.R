args = commandArgs(trailingOnly=TRUE)

library(tidyverse)

regex.finder <- "(.*)_keptReads_[0-9]*_to_[0-9]*_(.*)\\.txt"

setwd(args[1])

bam.reference <- read_csv(args[2], col_names = F)

final.dir = args[3]

run.name = args[4]

all.in <- list.files()

bams.tbl <- all.in %>% str_match(regex.finder) %>% as_tibble() %>% filter(!is.na(V1))

unique.ids <- bams.tbl %>% pull("V3") %>% unique()

for(uid in unique.ids){
  relevant.names <- bams.tbl %>% filter(V3 == uid) %>% pull("V1")
  
  to.keep <- relevant.names %>% map(function(a.file){
    hi <- read_csv(a.file, col_names = F, col_types = cols())
    if(nrow(hi) > 0){
      hi
    } else {
      NULL
    }
    
  }) %>% bind_rows %>% unique
  
  keptReads.fn <- paste0(unique(bams.tbl$V2), "_keptReads_allWindows_", uid, ".txt")
  
  original.bam.name <- bam.reference$X1[which(bam.reference$X3 == uid)]
  new.bam.name <- paste0(final.dir, "/", map_chr(original.bam.name, function(x){
    just.fn = basename(x)
    paste0(substr(just.fn, 1, nchar(just.fn) - 4), "_", run.name, "_cleaned.bam")
  }))
  
  write(paste0(original.bam.name, " ", new.bam.name, " -F ", args[1], "/", keptReads.fn), file = "ENRFB_input.txt", append = T)
  
  write_csv(to.keep, keptReads.fn, col_names = F)
  
}

