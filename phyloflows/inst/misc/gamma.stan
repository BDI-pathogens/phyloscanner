data {
  int<lower=1> L;  // number of pairs
  int n[L];  // transmission counts
  int<lower=1> N_xi; // number of sampling strata 
  matrix[N_xi,2] shape; // parameters of sampling fraction distributions
  int xi_id_src[L]; //indices of source strata
  int xi_id_rec[L]; //indices of recipient strata
}
transformed data {
  real alpha = 0.8/L; // gamma prior parameter
}
parameters {
  vector<lower=0>[L] lambda; // poisson rate
  vector<lower=0,upper=1>[N_xi] xi;   // sampling fractions 
}
transformed parameters {
  vector[L] xi_pair; // sampling fractions of pairs
  vector[L] betav; // average actual transmission counts for pairs
  real beta;
  for (i in 1:L){
   xi_pair[i] = xi[xi_id_src[i]] * xi[xi_id_rec[i]];
   if (Y[i]==0){
     betav[i] = 1/xi_pair[i]-1;}
   else{
     betav[i] = Y[i]/xi_pair[i];}
  }
  beta = 0.8/sum(betav); //  gamma prior parameter
}
model {
  for (i in 1:N_xi){
    xi[i]~ beta(shape[i,1],shape[i,2]); // sampling fraction priors
  }
  for (i in 1:L){
    lambda[i] ~ gamma(alpha,beta); // poisson rate priors
    n[i] ~  poisson(lambda[i] * xi_pair[i]); // likelihood
  }
}

