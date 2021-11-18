library(rstan)
library(data.table)	

indir <- '/rds/general/user/mm3218/home/git/phyloscanner/phyloflows/inst'
outdir <- '/rds/general/user/mm3218/home/projects/2021/phyloflows'

if(0){
  indir <- '~/git/phyloscanner/phyloflows/inst'
  outdir <- '~/Box\ Sync/2021/phyloflows'
}

DEBUG <- F
stan_model <- '211115'
lab <- stan_model

path.to.stan.data <- file.path(outdir, 'input.rda')
path.to.stan.model <- file.path(indir, 'stan_models', paste0('gp_', stan_model, '.stan'))
outdir.lab <- file.path(outdir, lab); dir.create(outdir.lab)

load(path.to.stan.data)
source(file.path(indir, 'functions', 'stan_utils.R'))

# stan data
stan_data <- ageanalysis(outdir = outdir)
stan_data <- add_2D_splines_stan_data(stan_data, spline_degree = 3, n_knots_rows = 15, n_knots_columns = 15, AGES = floor(sort(unique(stan_data$x[,1]))))

# save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=file.path(outdir.lab, paste0("stanin_",lab,".RData")) )

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

save(fit,file = file.path(outdir.lab, 'gp_211115.rda'))



