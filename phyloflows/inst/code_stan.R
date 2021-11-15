library(rstan)
library(data.table)	

indir <- '~/git/phyloscanner/phyloflows/inst'
outdir <- '~/Box\ Sync/2021/phyloflows'

DEBUG <- F
path.to.stan.data <- file.path(outdir, 'input.rda')
path.to.stan.model <- file.path(indir, 'stan_models', 'gp_211115.stan')

source(file.path(indir, 'input.R'))

# stan data
stan_data <- ageanalysis(outdir = outdir)
stan_data <- add_2D_splines_stan_data(stan_data, spline_degree = 3, n_knots_rows = 12, n_knots_columns = 10)

# run stan model
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
model = rstan::stan_model(path.to.stan.model)

if(DEBUG){
  fit <- sampling(model, data = data.fit, iter = 10, warmup = 5, chains=1, thin=1)
}else{
  fit <- sampling(model, data = data.fit,
                  iter = 3000, warmup = 500, chains=4, thin=1, seed = 5,
                  verbose = FALSE, control = list(adapt_delta = 0.999,max_treedepth=15))
}

save(fit,file = file.path(outdir, 'gp_211115.rda'))


