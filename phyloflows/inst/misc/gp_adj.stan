functions {
  vector lambda_nD(real[] L, int[] m, int D) {
    vector[D] lam;
    for(i in 1:D){
      lam[i] = ((m[i]*pi())/(2*L[i]))^2; }
    return lam;
  }

  real spd_nD(real alpha, row_vector rho, vector w, int D) {
    real S;
    S = alpha^2 * sqrt(2*pi())^D * prod(rho) * exp(-0.5*((rho .* rho) * (w .* w)));
    return S;
  }

  vector phi_nD(real[] L, int[] m, matrix x) {
    int c = cols(x);
    int r = rows(x);
    matrix[r,c] fi;
    vector[r] fi1;
    for (i in 1:c){
      fi[,i] = 1/sqrt(L[i])*sin(m[i]*pi()*(x[,i]+L[i])/(2*L[i]));
    }
    fi1 = fi[,1];
    for (i in 2:c){
      fi1 = fi1 .* fi[,i];
    }
    return fi1;
  }
}
data {
  int<lower=1> N;
  int<lower=1> D;
  matrix[N,D] x;
  int<lower=0> y[N];

  real L[D];
  int<lower=1> M;	
  int<lower=1> M_nD;			
  int indices[M_nD,D];	
}

transformed data {
  matrix[N,M_nD] PHI;
	
  for (m in 1:M_nD){ 
    PHI[,m] = phi_nD(L, indices[m,], x); 
  }
}

parameters {
  row_vector<lower=0>[D] rho;
  real<lower=0> alpha;
  real mu;
  vector[M_nD] beta;
}

transformed parameters {
  vector[N] f;
  vector[M_nD] diagSPD;
  vector[M_nD] SPD_beta;

  for(m in 1:M_nD){ 
    diagSPD[m] =  sqrt(spd_nD(alpha, rho, sqrt(lambda_nD(L, indices[m,], D)), D)); 
  }
  SPD_beta = diagSPD .* beta;
  f= PHI * SPD_beta;
}


model {
  beta ~ normal(0,1);
  rho[2] ~ inv_gamma(2.60877,7.73394);
  rho[1] ~ inv_gamma(3.48681,9.21604);
  alpha ~ normal(0, 10);
  mu ~ normal(0, 10);
  y ~ poisson_log(f + mu);
}


