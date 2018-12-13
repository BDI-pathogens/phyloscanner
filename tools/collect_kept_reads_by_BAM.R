args = commandArgs(trailingOnly=TRUE)

library(tidyverse)

regex.finder <- "(.*)_keptReads_[0-9]*_to_[0-9]*_(.*)\\.txt"

setwd(args[1])

all.in <- list.files()

bams.tbl <- all.in %>% str_match(regex.finder) %>% as.tibble() %>% filter(!is.na(V1))

unique.ids <- bams.tbl %>% pull("V3") %>% unique()

for(uid in unique.ids){
  relevant.bams <- bams.tbl %>% filter(V3 == uid) %>% pull("V1")
  
  to.keep <- relevant.bams %>% map(function(a.file){
    hi <- read_csv(a.file, col_names = F, col_types = cols())
    if(nrow(hi) > 0){
      hi
    } else {
      NULL
    }
    
  }) %>% bind_rows %>% unique()
  
  write_csv(to.keep, paste0(unique(bams.tbl$V2), "_keptReads_allWindows_", uid, ".txt"), col_names =F)
  
}

