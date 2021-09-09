\\spectral density
    functions {
  vector spd_cov_exp_quad(vector[] x, real alpha, vector rho) {
    int L_per_group = dims(x)[1];
    int D = dims(x)[2];
    int Drho = rows(rho);
    vector[L_per_group] out;
    if (Drho == 1) {
      // one dimensional or isotropic GP
      real constant = square(alpha) * (sqrt(2 * pi()) * rho[1])^D;
      real neg_half_rho2 = -0.5 * square(rho[1]);
      for (m in 1:L_per_group) {
        out[m] = constant * exp(neg_half_rho2 * dot_self(x[m]));
      }
    } else {
      // multi-dimensional non-isotropic GP
      real constant = square(alpha) * sqrt(2 * pi())^D * prod(rho);
      vector[Drho] neg_half_rho2 = -0.5 * square(rho);
      for (m in 1:L_per_group) {
        out[m] = constant * exp(dot_product(neg_half_rho2, square(x[m])));
      }
    }
    return out;
  }

  vector gpa_coefficient(int M_nD, real alpha, vector rho,
  vector beta, vector[] slambda) { 
    vector[M_nD] diag_spd = sqrt(spd_cov_exp_quad(slambda, 
    alpha, rho));
    return diag_spd .* beta;
  }
}
data {
  int<lower=1> L; // number of pairs
  int<lower=1> L_group; // male->female and female->male
  int<lower=1> L_per_group; // number of age group pairs
  int<lower=1> D; // 2
  matrix[L_per_group,D] x; // age of source and recipients
  int n[L_per_group, L_group]; // transmission counts
	
  int<lower=1> M_nD; // total number of basis functions
  matrix[L_per_group, M_nD] Xgp; // eigenfunctions
  vector[D] slambda[M_nD]; // eigenvectors	

  int<lower=1> N_xi;	// number of sampling strata
  matrix[N_xi,2] shape; // parameters of sampling fraction distributions
  int xi_id_src[L_per_group,L_group]; // indices of groups for sources
  int xi_id_rec[L_per_group,L_group]; // indices of groups for recipients
  int<lower=0> id_mf; // male -> female pair indices
  int<lower=0> id_fm; // female -> male pair indices
  
  real d1_invg_par1_mf;  
  // parameters of inverse gamma prior (male -> female, source)
  real d1_invg_par2_mf;  
  // parameters of inverse gamma prior (male -> female, source)
  real d2_invg_par1_mf;  
  // parameters of inverse gamma prior (male -> female, recipient)
  real d2_invg_par2_mf;  
  // parameters of inverse gamma prior (male -> female, recipient)
  real d1_invg_par1_fm;  
  // parameters of inverse gamma prior (female -> male, source)
  real d1_invg_par2_fm;  
  // parameters of inverse gamma prior (female -> male, source)
  real d2_invg_par1_fm;  
  // parameters of inverse gamma prior (female -> male, recipient)
  real d2_invg_par2_fm;  
  // parameters of inverse gamma prior (female -> male, recipient)
}

parameters {
  vector<lower=0>[D] rho[L_group];// lengthscale  
  vector<lower=0>[L_group] alpha; // marginal standard deviation
  real mu_mf; // differences of male -> female flows
  real mu; // overall intercept
  vector[M_nD] beta[2];
  vector<lower=0,upper=1>[N_xi] xi; // sampling proportions
}

transformed parameters {
  vector[L_per_group] lxi_pair[L_group];
  vector[M_nD] coefficient[2];  
  vector[L_per_group] llbd[L_group];  // log poisson rate of actual flows
  vector[L_group] muv;
  vector[L_per_group] f_mf; // gaussian process male -> female
  vector[L_per_group] f_fm; // gaussian process female -> male
  vector[L_per_group] f[L_group];

  muv = rep_vector(0,L_group);
  muv[id_mf] = muv[id_mf] + mu_mf;


  coefficient[1] = gpa_coefficient(M_nD, alpha[1], rho[1],
  beta[1], slambda);
  f_mf = Xgp * coefficient[1];
  coefficient[2] = gpa_coefficient(M_nD, alpha[2], rho[2],
  beta[2], slambda);
  f_fm = Xgp * coefficient[2];
  f[id_mf]=f_mf;
  f[id_fm]=f_fm;

  for (i in 1:L_group){
    lxi_pair[i] = log(xi[xi_id_src[,i]]) + log(xi[xi_id_rec[,i]]);
    llbd[i] =mu +  muv[i] + f[i] + lxi_pair[i];
  }

}

model {
  for (i in 1:N_xi){
    xi[i] ~ beta(shape[i,1],shape[i,2]);
  }
  for (i in 1:2){
    beta[i] ~ normal(0,1);
  }
  rho[1][1] ~ inv_gamma(d1_invg_par1_mf, d1_invg_par2_mf); 
  rho[1][2] ~ inv_gamma(d2_invg_par1_mf, d2_invg_par2_mf); 
  rho[2][1] ~ inv_gamma(d1_invg_par1_fm, d1_invg_par2_fm); 
  rho[2][2] ~ inv_gamma(d2_invg_par1_fm, d2_invg_par2_fm);
  
  alpha ~ normal(0, 10);
  mu_mf ~ normal(0, 10);
  mu ~ normal(0, 10);

  for (i in 1:L_group){
    n[,i] ~ poisson_log(llbd[i]);
  }
}
