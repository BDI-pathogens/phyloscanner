library(rstan)
library(data.table)	
library(ggplot2)

.indir <- '/rds/general/user/mm3218/home/git/phyloscanner/phyloflows/inst'
.outdir <- '/rds/general/user/mm3218/home/projects/2021/phyloflows'

if(0){
  .indir <- '~/git/phyloscanner/phyloflows/inst'
  .outdir <- '~/Box\ Sync/2021/RCCS/outputs'
}

lab <- "MRC_FALSE_OnlyHTX_TRUE_threshold_0.6"
stan_model <- 'gp_211117'

path.to.stan.data <- file.path(.outdir, lab, paste0("stanin_",lab,".RData"))
path.to.stan.output <- file.path(.outdir, lab, paste0(stan_model, '_', lab, '.rds'))

load(path.to.stan.data)
indir <- .indir; outdir <- .outdir
outdir.lab <- file.path(outdir, lab)

source(file.path(indir, 'functions', 'postprocessing_functions.R'))

# table
df_direction <- data.table(index_direction = 1:2, is_mf = stan_data$is_mf)
df_direction[, label_direction := ifelse(is_mf == 1, 'Male -> Female', 'Female -> Male')]
df_age$index_age <- 1:nrow(df_age)

# samples 
fit <- readRDS(path.to.stan.output)
samples <- rstan::extract(fit)

# convergence diagnostics 
make_convergence_diagnostics_stats(fit, outdir.lab)

# intensity of the poisson process
intensity_PP <- summarise_var_by_age_direction(samples, 'log_lambda', df_direction, df_age, outdir.lab, transform = 'exp')
count_data <- prepare_count_data(stan_data)
plot_intensity_PP(intensity_PP, count_data, outdir.lab)
  

