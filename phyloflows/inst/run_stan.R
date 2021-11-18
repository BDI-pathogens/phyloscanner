library(rstan)
library(data.table)	

.indir <- '/rds/general/user/mm3218/home/git/phyloscanner/phyloflows/inst'
.outdir <- '/rds/general/user/mm3218/home/projects/2021/phyloflows'

if(0){
  .indir <- '~/git/phyloscanner/phyloflows/inst'
  .outdir <- '~/Box\ Sync/2021/RCCS/outputs'
}

lab <- "MRC_FALSE_OnlyHTX_TRUE_threshold_0.6"
stan_model <- 'gp_211117'
DEBUG <- F

path.to.stan.data <- file.path(.outdir, paste0("stanin_",lab,".RData"))
path.to.stan.model <- file.path(indir, 'stan_models', paste0(stan_model, '.stan'))

load(path.to.stan.data)
indir <- .indir; outdir <- .outdir; outdir.lab <- file.path(outdir, lab)

# run stan model
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
model = rstan::stan_model(path.to.stan.model)

if(DEBUG){
  fit <- sampling(model, data = stan_data, iter = 10, warmup = 5, chains=1, thin=1)
}else{
  fit <- sampling(model, data = stan_data,
                  iter = 3000, warmup = 500, chains=4, thin=1, seed = 5,
                  verbose = FALSE, control = list(adapt_delta = 0.999,max_treedepth=15))
}

sum = summary(fit)
sum$summary[which(sum$summary[,9] < 100),]

file = file.path(outdir.lab, paste0(stan_model, '_', lab, '.rds'))
saveRDS(fit,file = file)
     