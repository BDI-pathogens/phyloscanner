keep.likely.transmission.pairs <- function(dchain, threshold){
  
  dchain <- dchain[SCORE_LINKED>threshold]
  dchain[SCORE_DIR_12 <= threshold & SCORE_DIR_21 <= threshold, EST_DIR:='unclear']
  dchain[SCORE_DIR_12 > threshold, EST_DIR:='12']
  dchain[SCORE_DIR_21 > threshold, EST_DIR:='21']
  
  # find source recipient
  dchain <- dchain[EST_DIR != 'unclear']
  dchain[, `:=` (SOURCE=H1, RECIPIENT=H2)]
  dchain[EST_DIR == '21', `:=` (SOURCE=H2, RECIPIENT=H1) ]
  dchain[, `:=` (H1=NULL, H2=NULL)]
}


get.meta.data <- function(meta.rccs.1, meta.rccs.2, meta.mrc, time.first.positive, anonymisation.keys, community.keys){
  
  colnames(meta.rccs.1) <- tolower(colnames(meta.rccs.1))
  colnames(meta.rccs.2) <- tolower(colnames(meta.rccs.2))
  colnames(meta.mrc) <- tolower(colnames(meta.mrc))
  
  #
  # process RCCS meta-data
  
  # process first meta data for RCCS
  meta.rccs.1[, pt_id := paste0('RK-', rccs_studyid)]
  tmp <- meta.rccs.1[, list(date = min(date), comm_num = min(comm_num)), by = 'pt_id'] # keep comm of first visit
  meta.rccs.1 <- unique( merge(meta.rccs.1[, .(pt_id, date, comm_num)], tmp, by = c('pt_id', 'date', 'comm_num')) )
  meta.rccs.1 <- meta.rccs.1[, .(pt_id, comm_num)] 
  meta.rccs.1[, cohort_round := 'R15-R16']
  
  # process second meta data for RCCS
  meta.rccs.2 <- unique( meta.rccs.2[, .(pt_id, sex)] ) 
  
  # merge two data sources for RCCS
  meta.rccs <- merge(meta.rccs.2, meta.rccs.1, by = 'pt_id', all.x = T)
  meta.rccs[is.na(cohort_round), cohort_round := 'R17-R18']
  
  # add time first positive
  time.first.positive <- time.first.positive[pt_id %in% unique(meta.rccs$pt_id)]
  time.first.positive[, visit_dt := as.Date(visit_dt, '%d/%m/%Y') ]
  time.first.positive[, date_birth := visit_dt - age_at_visit * 365]
  time.first.positive[, age_first_positive := as.numeric(as.Date(first_pos_dt, '%d/%m/%Y') - date_birth) / 365]
  # keep only first visit
  tmp <- time.first.positive[, list(visit_dt = min(visit_dt)), by = 'pt_id']
  time.first.positive <- unique( merge(time.first.positive, tmp, by = c('pt_id', 'visit_dt')) )
  # merge
  meta.rccs <- merge(meta.rccs, time.first.positive[, .(pt_id, age_first_positive)], by = 'pt_id')
  
  stopifnot(length(unique(meta.rccs$pt_id)) == nrow(meta.rccs))
  cat('There is ', nrow(meta.rccs), ' individuals included in the RCCS meta-data\n')
  
  
  #
  # process MRC meta-data
  meta.mrc <- unique( meta.mrc[, .(pt_id, sex, age_enrol)] ) 
  
  stopifnot(length(unique(meta.mrc$pt_id)) == nrow(meta.mrc))
  cat('There is ', nrow(meta.mrc), ' individuals included in the MRC meta-data\n\n')
  
  
  #
  # merge RCCS and MRC meta data
  meta.rccs[, cohort := 'RCCS']
  meta.mrc[, cohort := 'MRC']
  meta <- rbind(meta.rccs, meta.mrc, fill = TRUE)
  cat('There is ', nrow(meta), ' individuals included in the meta-data\n')
  
  
  # transform variables
  meta[, age_infection := as.numeric(ifelse(is.na(age_first_positive), age_enrol, age_first_positive)) - 1]
  cat('There is ', nrow(meta[!is.na(age_infection)]), ' individuals included in the meta-data with proxy for the age at infection\n')
  
  
  #
  # last changes
  # anonymisation keys
  colnames(anonymisation.keys) <- tolower(colnames(anonymisation.keys))
  meta <- merge(meta, anonymisation.keys, by = 'pt_id')
  
  # community key
  colnames(community.keys) <- tolower(colnames(community.keys))
  community.keys[, comm := ifelse(strsplit(comm_num_a, '')[[1]][1] == 'f', 'fishing', 'island'), by = 'comm_num_a']
  meta <- merge(meta, community.keys, by.x = 'comm_num', by.y = 'comm_num_raw', all.x = T)
  
  # keep only variable of interest
  meta <- meta[, .(pt_id, aid, sex, age_infection, age_first_positive, age_enrol, cohort_round, cohort, comm)]
  
  # remove very young patients (bug?)
  cat('\nExcluding very young patients"')
  print_table(meta[age_infection < 11, .(aid, age_infection)])
  meta <- meta[is.na(age_infection) | age_infection >= 11]

  # rm duplicate
  meta <- unique(meta)
  
  return(meta)
}

pairs.get.meta.data <- function(dchain, meta){
  
  # stopifnot(unique(dchain$SOURCE) %in% unique(meta_data$aid))
  # stopifnot(unique(dchain$RECIPIENT) %in% unique(meta_data$aid))
  
  # merge by source
  tmp <- copy(meta)
  names(tmp) = paste0(names(tmp), '.SOURCE')
  tmp1 <- merge(dchain, tmp, by.x = 'SOURCE', by.y = 'aid.SOURCE')
  
  # merge by reicpient
  tmp <- copy(meta)
  names(tmp) = paste0(names(tmp), '.RECIPIENT')
  tmp1 <- merge(tmp1, tmp, by.x = 'RECIPIENT', by.y = 'aid.RECIPIENT')
  
  # individuals without meta data
  missing_indiv = unique(c(dchain$SOURCE, dchain$RECIPIENT)[!(c(dchain$SOURCE, dchain$RECIPIENT) %in% meta$aid)])
  cat('There are ', length(missing_indiv), 'indivs without meta data:\n' )
  cat(missing_indiv)
  
  return(tmp1)
}

print.statements.about.pairs <- function(pairs, outdir){
  
  cat('\nThere is ', nrow(pairs), ' source-recipient pairs\n\n')
  
  cat(nrow(pairs[!is.na(age_infection.SOURCE) & !is.na(age_infection.RECIPIENT)]), ' pairs have a proxy for the time of infection of the source and recipient\n')
  cat(nrow(pairs[((sex.SOURCE == 'F' & sex.RECIPIENT == 'M') | (sex.SOURCE == 'M' & sex.RECIPIENT == 'F')) & (!is.na(age_infection.SOURCE) & !is.na(age_infection.RECIPIENT))]), ' pairs are heteroxuals have a proxy for the time of infection of the source and recipient\n\n')                
  
  cat('\nPairs by cohort')
  tab <- pairs[, list(count = .N), by = c('cohort.SOURCE', 'cohort.RECIPIENT')]
  print_table(tab)
  
  cat('\nPairs enrolled in RCCS by cohort round')
  tab <- pairs[cohort.RECIPIENT == 'RCCS' & cohort.SOURCE == 'RCCS', list(count = .N), by = c('cohort_round.SOURCE', 'cohort_round.RECIPIENT')]
  print_table(tab)
  
  cat('\nPairs by sex')
  tab <- pairs[, list(count = .N), by = c('sex.SOURCE', 'sex.RECIPIENT')]
  print_table(tab)
  
  cat('\nPairs by community')
  tab <- pairs[, list(count = .N), by = c('comm.SOURCE', 'comm.RECIPIENT')]
  print_table(tab)
  
  cat('\nPairs by age of infection\n')
  plot_hist_age_infection(pairs, outdir)

}

print_table <- function(table) print(knitr::kable(table))

get.age.map <- function(pairs){
  ages <- pairs[, {
    min_age = floor(min(c(age_infection.SOURCE, age_infection.RECIPIENT)))
    max_age = floor(max(c(age_infection.SOURCE, age_infection.RECIPIENT)))
    list(age = min_age:max_age)}]
  
  
  age_map <- data.table(expand.grid(age_infection.SOURCE = ages$age, age_infection.RECIPIENT = ages$age))
  age_map <- age_map[order(age_infection.SOURCE, age_infection.RECIPIENT)]
}
