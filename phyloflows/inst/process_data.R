library(data.table)
library(phyloscannerR)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(knitr)

# change as appropriate
indir.repository <- '~/git/phyloscanner/phyloflows/inst'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
outdir <- '~/Box\ Sync/2021/phyloflows/'

# indicators 
include.mrc <- F
include.only.heterosexual.pairs <- T
threshold.likely.connected.pairs <- 0.6
lab <- paste0('MRC_', include.mrc, '_OnlyHTX_', include.only.heterosexual.pairs, '_threshold_', threshold.likely.connected.pairs)

# file paths
file.path.chains.data <- file.path(indir.deepsequence_analyses,'210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks.rda')
file.path.meta.data.rccs.1 <- file.path(indir.deepsequence_analyses, 'RakaiPangeaMetaData_v2.rda')
file.path.meta.data.rccs.2 <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '200316_pangea_db_sharing_extract_rakai.csv')
file.path.meta.data.mrc <- file.path(indir.deepsequencedata, 'PANGEA2_MRC','200319_pangea_db_sharing_extract_mrc.csv')
file.anonymisation.keys <- file.path(indir.deepsequence_analyses,'important_anonymisation_keys_210119.csv')
file.community.keys <- file.path(indir.deepsequence_analyses,'community_names.csv')
file.time.first.positive <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '211111_pangea_db_sharing_extract_rakai_age_firstpos_lastneg.csv')
outdir.lab <- file.path(outdir, lab); dir.create(outdir.lab)
  
# load functions
source(file.path(indir.repository, 'functions', 'summary_functions.R'))
source(file.path(indir.repository, 'functions', 'plotting_functions.R'))
source(file.path(indir.repository, 'functions', 'stan_utils.R'))

# load files
load(file.path.chains.data)
load(file.path.meta.data.rccs.1)
meta.rccs.1 <- as.data.table(rccsData)
meta.rccs.2 <- as.data.table( read.csv(file.path.meta.data.rccs.2))
meta.mrc <- as.data.table( read.csv(file.path.meta.data.mrc))

# load keys
anonymisation.keys <- read.csv(file.anonymisation.keys)
community.keys <- as.data.table( read.csv(file.community.keys) )
time.first.positive <- as.data.table( read.csv(file.time.first.positive) )

# get meta data
meta_data <- get.meta.data(meta.rccs.1, meta.rccs.2, meta.mrc, time.first.positive, anonymisation.keys, community.keys)

# get likely transmission pairs
chain <- keep.likely.transmission.pairs(as.data.table(dchain), threshold.likely.connected.pairs)

# merge meta data to source and recipient
pairs <- pairs.get.meta.data(chain, meta_data)

if(!include.mrc){
  cat('Keep only pairs in RCCS')
  pairs <- pairs[cohort.RECIPIENT == 'RCCS' & cohort.SOURCE == 'RCCS']
}
if(include.only.heterosexual.pairs){
  cat('Keep only heterosexual pairs')
  pairs <- pairs[(sex.RECIPIENT == 'M' & sex.SOURCE == 'F') | (sex.RECIPIENT == 'F' & sex.SOURCE == 'M')]
}

print.statements.about.pairs(copy(pairs), outdir.lab)

# keep only pairs with source-recipient with proxy for the time of infection
pairs <- pairs[!is.na(age_infection.SOURCE) & !is.na(age_infection.RECIPIENT)]
plot_age_source_recipient(pairs[sex.SOURCE == 'M' & sex.RECIPIENT == 'F'], 'Male -> Female', 'MF', outdir.lab)
plot_age_source_recipient(pairs[sex.SOURCE == 'F' & sex.RECIPIENT == 'M'], 'Female -> Male', 'FM', outdir.lab)

# prepare age map
df_age <- get.age.map(pairs)

# prepare stan data
stan_data <- prepare_stan_data(pairs, df_age)
stan_data <- add_2D_splines_stan_data(stan_data, spline_degree = 3, n_knots_rows = 15, n_knots_columns = 15, unique(df_age$age_infection.SOURCE))

## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=file.path(outdir.lab, paste0("stanin_",lab,".RData")) )




