library(rstan)
setwd('~/ageanalysis/')
load('input.rda')

# I usually run for 4 times to get 1e4 samples
fit <- stan(file = file.path('gpxx.stan'),
            data = data.fit,
            iter = 3000,  warmup = 500, chains=1, thin=1, seed = 5,
            algorithm = "NUTS", verbose = FALSE,
            control = list(adapt_delta = 0.999,max_treedepth=15))
save(fit,file = 'gpxx.rda')


