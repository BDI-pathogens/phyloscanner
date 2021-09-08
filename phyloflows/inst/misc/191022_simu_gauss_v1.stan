data {
  int<lower=1> N;
  int<lower=1> D;
  // the row vector code is taken from here:
  // https://github.com/stan-dev/example-models/blob/master/misc/gaussian-process/gp-sim-multi.stan
  vector[D] x[N];  
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;    
}

transformed data {
  // code for multi-dimensional kernel by hand is here: 
  // https://aheblog.com/2016/12/07/geostatistical-modelling-with-r-and-stan/   
  matrix[N, N] cov =   cov_exp_quad(x, alpha, rho)
                     + diag_matrix(rep_vector(1e-10, N));
  matrix[N, N] L_cov = cholesky_decompose(cov);
}

parameters {}
model {}

generated quantities {
  vector[D] mu[N]; 
  vector[D] f[N]; 
  vector[D] y[N];
   
  for (d in 1:D)
    mu[d] = rep_vector(0, N); 
  f = multi_normal_cholesky_rng(mu, L_cov);
  for (d in 1:D)
  	for (n in 1:N)
    	y[d][n] = normal_rng(f[d][n], sigma);
}
