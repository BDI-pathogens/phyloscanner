data {
  int<lower=1> N;  // number of observations
  int Y[N];  // response variable
  // sampling fractions
  int<lower=1> N_xi; 
  matrix[N_xi,2] shape;
  int xi_id[N,2];
  real alpha;
}
transformed data {
}
parameters {
  vector<lower=0>[N] lambda;
  vector<lower=0,upper=1>[N_xi] xi;
}
transformed parameters {
  vector[N] lxi_pair; 
  vector[N] betav;
  real beta;
  for (i in 1:N){
   lxi_pair[i] = log(xi[xi_id[i,1]]) + log(xi[xi_id[i,2]]);
   if (Y[i]==0){
     betav[i] = 1/(xi[xi_id[i,1]]*xi[xi_id[i,2]])-1;}
   else{
     betav[i] = Y[i]/(xi[xi_id[i,1]]*xi[xi_id[i,2]]);}
  }
  beta = sum(betav);
}
model {
  for (i in 1:N_xi){
    target += beta_lpdf(xi[i]|shape[i,1],shape[i,2]);
  }
  for (i in 1:N){
    target += gamma_lpdf(lambda[i]|alpha,beta); 
  }
  // likelihood including all constants
  target += poisson_log_lpmf(Y | lxi_pair + log(lambda));
}


