functions {
  matrix L_cov_exp_quad_ARD(vector[] x,
                            real alpha,
                            vector rho) {
    int N = size(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + 1e-10;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha
                      * exp(-0.5 * dot_self((x[i] - x[j]) ./ rho));
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + 1e-10;
    return cholesky_decompose(K);
  }
}
data {
  int<lower=1> N_group;
  int<lower=1> N_per_group;
  int<lower=1> D;
  vector[D] x[N_per_group];
  int y[N_per_group, N_group];
	
  int<lower=1> N_xi;	
  matrix[N_xi,2] shape;
  int xi_id_src[N_per_group,N_group];
  int xi_id_rec[N_per_group,N_group];
  int id_mf[4];
  int id_fm[4];
	
}

transformed data {
}

parameters {
  vector<lower=0>[D] rho[2];
  vector<lower=0>[2] alpha;
  vector[N_group] mu;
  vector[N_per_group] beta[N_group];
  vector<lower=0,upper=1>[N_xi] xi;
}

transformed parameters {
  vector[N_per_group] lxi_pair[N_group];  
  vector[N_per_group] llbd[N_group];  
  vector<lower=0>[N_group] alphav;
  vector<lower=0>[D] rhov[N_group];
  matrix[N_per_group, N_per_group] L_cov[N_group];
  vector[N_per_group] f[N_group];

  alphav[id_mf] = rep_vector(alpha[1],4);
  alphav[id_fm] = rep_vector(alpha[2],4);

  for (i in 1:4){
    rhov[id_mf[i]] = rho[1];
    rhov[id_fm[i]] = rho[2];
  }

  for (i in 1:N_group){
    L_cov[i] = L_cov_exp_quad_ARD(x,alphav[i],rhov[i]);
    f[i] = L_cov[i] * beta[i]; 
    lxi_pair[i] = log(xi[xi_id_src[,i]]) + log(xi[xi_id_rec[,i]]);
    llbd[i] = mu[i] + f[i] + lxi_pair[i];
  }

}


model {
  for (i in 1:N_xi){
    xi[i] ~ beta(shape[i,1],shape[i,2]);
  }
  for (i in 1:N_group){
    beta[i] ~ normal(0,1);
  }
  rho[1][1] ~ inv_gamma(8.91924,17.2902);
  rho[1][2] ~ inv_gamma(4.62909,11.0337); 
  rho[2][1] ~ inv_gamma(4.62909,11.0337);
  rho[2][2] ~ inv_gamma(8.91924,17.2902);
  
  alpha ~ normal(0, 10);
  mu ~ normal(0, 10);

  for (i in 1:N_group){
    y[,i] ~ poisson_log(llbd[i]);
  }
}

